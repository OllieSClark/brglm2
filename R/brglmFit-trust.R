# Copyright (C) 2016- Ioannis Kosmidis

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# Trust-region iteration for brglmFit

#' Solve the adjusted score equations using a trust-region iteration
#'
#' Each iteration approximately minimises a quadratic model of the
#' squared norm of the adjusted score for the regression parameters
#' within a trust region, accepts or rejects the step according to the
#' ratio of the actual to the predicted reduction, and updates the radius
#' (Nocedal & Wright, 2006, Algorithm 4.1). The dispersion is then
#' updated with a Newton step, halved until the dispersion is positive.
#'
#' @param betas Starting values for the regression parameters
#' @param dispersion Starting value for the dispersion
#' @param y Response vector
#' @param x Design matrix (dense, or a CsparseMatrix)
#' @param weights Observation weights
#' @param offset Offset vector
#' @param family GLM family object
#' @param fixed_totals NULL or grouping vector for fixed totals
#' @param row_totals NULL or vector of row totals
#' @param no_dispersion Logical indicating if dispersion is fixed
#' @param nobs Number of observations
#' @param nvars Number of variables
#' @param keep Logical vector indicating which observations to keep
#' @param df_residual Residual degrees of freedom
#' @param adjustment_function Score adjustment function
#' @param d1_transformed_dispersion First derivative of the dispersion transformation
#' @param d2_transformed_dispersion Second derivative of the dispersion transformation
#' @param control Control parameters, as returned by brglmControl
#'
#' @return List with the estimates, the adjusted gradients, the step
#'     components at the estimates, the number of iterations and
#'     whether the iteration failed, as for fit_fisher_scoring
#'
#' @references
#'
#' Nocedal J, Wright S J (2006). *Numerical Optimization*, 2nd edition.
#' Springer.
#'
#' Eisenstat S C, Walker H F (1996). Choosing the forcing terms in an
#' inexact Newton method. *SIAM Journal on Scientific Computing*, **17**,
#' 16-32.
fit_trust_region <- function(betas, dispersion, y, x, weights, offset, family,
                             fixed_totals = NULL, row_totals = NULL,
                             no_dispersion = FALSE, nobs, nvars, keep, df_residual,
                             adjustment_function,
                             d1_transformed_dispersion, d2_transformed_dispersion,
                             control) {
    is_sparse <- inherits(x, "sparseMatrix")

    ## Evaluate at the starting values
    theta <- c(betas, dispersion)
    transformed_dispersion <- eval(control$Trans)

    ## The inverse of the information about the regression parameters is
    ## only needed by the median adjustment and when the dispersion is updated
    need_inverse <- !no_dispersion || control$type == "AS_median"

    fit <- compute_fit(pars = theta,
                       y = y,
                       x = x,
                       weights = weights,
                       offset = offset,
                       family = family,
                       fixed_totals = fixed_totals,
                       row_totals = row_totals,
                       no_dispersion = no_dispersion,
                       nobs = nobs,
                       nvars = nvars,
                       keep = keep,
                       need_qr = TRUE,
                       need_hatvalues = TRUE,
                       need_inverse = need_inverse)

    step_components_beta <- compute_step_components(theta, fit = fit, level = 0,
                                                    adjustment_function = adjustment_function,
                                                    x = x, nobs = nobs, nvars = nvars,
                                                    weights = weights, df_residual = df_residual,
                                                    d1_transformed_dispersion = d1_transformed_dispersion,
                                                    d2_transformed_dispersion = d2_transformed_dispersion)
    step_components_zeta <- compute_step_components(theta, fit = fit, level = 1,
                                                    adjustment_function = adjustment_function,
                                                    x = x, nobs = nobs, nvars = nvars,
                                                    weights = weights, df_residual = df_residual,
                                                    d1_transformed_dispersion = d1_transformed_dispersion,
                                                    d2_transformed_dispersion = d2_transformed_dispersion)
    if (step_components_beta$failed_inversion) {
        warning("failed to invert the information matrix")
    }
    if (step_components_beta$failed_adjustment) {
        warning("failed to calculate score adjustment")
    }
    adjusted_grad_beta <- with(step_components_beta, grad + adjustment)
    adjusted_grad_zeta <- if (no_dispersion) NA_real_ else with(step_components_zeta, grad + adjustment)

    ## Trust-region constants
    Delta <- 1
    Delta_max <- 1e10
    radius_shrink <- 0.25
    radius_expand <- 0.75
    shrink_factor <- 0.25
    expand_factor <- 2
    eta <- 0.1
    cg_maxiter <- min(nvars, 50)

    ## Steps are assessed by the reduction in a merit function. For
    ## maximum likelihood and maximum penalized likelihood, the adjusted
    ## score is the gradient of the (penalized) log-likelihood, so its
    ## negative is used; the quadratic model minimised for the step is then
    ## a model of the merit itself, and every step is a descent direction.
    ## Otherwise, the squared norm of the adjusted score is used
    use_objective <- control$type %in% c("ML", "MPL_Jeffreys") && is.null(fixed_totals)
    objective <- function(fit) {
        penalty <- if (control$type == "MPL_Jeffreys") control$a * 2 * sum(log(abs(diag(fit$R_matrix)))) else 0
        sum(fit$deviance_residuals) / (2 * fit$dispersion) - penalty
    }

    ## Candidates with invalid linear predictors or means are rejected, as
    ## in glm.fit()
    valid_eta <- unless_null(family$valideta, function(eta) TRUE)
    valid_mu <- unless_null(family$validmu, function(mu) TRUE)

    failed <- FALSE
    if (control$maxit == 0) {
        iter <- 0L
    } else {
        for (iter in seq.int(control$maxit)) {
            r_adj_sq <- sum(adjusted_grad_beta^2)
            if (use_objective) {
                merit <- objective(fit)
            }

            ## Tolerance for the CG iteration (Eisenstat & Walker, 1996, choice 2)
            cg_tol <- if (iter == 1) 0.5 else min(0.5, 0.9 * (r_adj_sq / r_adj_sq_previous)^2)
            r_adj_sq_previous <- r_adj_sq

            ## For sparse x, products with the information are formed
            ## without the information matrix
            if (is_sparse) {
                precision <- fit$precision
                working_weights <- fit$working_weights
                Fmatvec <- function(d) {
                    precision * .crossprod_vec(working_weights * drop(x %*% d), x, nobs, nvars)
                }
                diag_F <- precision * .crossprod_vec(working_weights, x^2, nobs, nvars)
            } else {
                Fmatvec <- diag_F <- NULL
            }

            step <- cg_steihaug_pcg(r_adj = adjusted_grad_beta,
                                    F_info = if (is_sparse) NULL else fit$info_beta,
                                    Fmatvec = Fmatvec,
                                    diag_F = diag_F,
                                    Delta = Delta,
                                    tol = cg_tol,
                                    maxiter = cg_maxiter)

            ## Predicted reduction in the merit
            Fp <- if (is_sparse) Fmatvec(step$p) else drop(fit$info_beta %*% step$p)
            pred_reduction <- if (use_objective) {
                                  sum(adjusted_grad_beta * step$p) - 0.5 * sum(step$p * Fp)
                              } else {
                                  r_adj_sq - sum((adjusted_grad_beta - Fp)^2)
                              }

            betas_candidate <- betas + step$p
            theta_candidate <- c(betas_candidate, dispersion)
            fit_candidate <- try(compute_fit(pars = theta_candidate,
                                             y = y,
                                             x = x,
                                             weights = weights,
                                             offset = offset,
                                             family = family,
                                             fixed_totals = fixed_totals,
                                             row_totals = row_totals,
                                             no_dispersion = no_dispersion,
                                             nobs = nobs,
                                             nvars = nvars,
                                             keep = keep,
                                             need_qr = TRUE,
                                             need_hatvalues = FALSE,
                                             need_inverse = need_inverse), silent = TRUE)

            if (inherits(fit_candidate, "try-error") ||
                !valid_eta(fit_candidate$etas) || !valid_mu(fit_candidate$mus)) {
                Delta <- shrink_factor * Delta
                if (Delta < 1e-16) {
                    warning("trust region radius became too small")
                    failed <- TRUE
                    break
                }
                next
            }

            ## Actual reduction in the merit. For the squared norm of the
            ## adjusted score, as in quasi Fisher scoring, the hat values in
            ## the adjustment are held at the current point; with the hat
            ## values at the candidate, the model, which ignores the
            ## derivatives of the adjustment, can fail to predict any
            ## reduction for small steps
            if (use_objective) {
                actual_reduction <- merit - objective(fit_candidate)
            } else {
                fit_candidate$hatvalues <- fit$hatvalues
                adjustment_candidate <- adjustment_function(theta_candidate, fit = fit_candidate, level = 0,
                                                            x, nobs, nvars, weights)
                r_adj_candidate <- fit_candidate$grad_beta + adjustment_candidate
                actual_reduction <- r_adj_sq - sum(r_adj_candidate^2)
            }

            ## A non-finite reduction rejects the step. Close to the
            ## solution, the reduction in the objective is lost in rounding
            ## error, and steps are accepted unless the objective clearly
            ## increases
            noise <- if (use_objective) 100 * .Machine$double.eps * max(1, abs(merit)) else 1e-16
            rho <- if (!is.finite(actual_reduction)) {
                -Inf
            } else if (!is.finite(pred_reduction) || abs(pred_reduction) <= noise) {
                if (actual_reduction > -noise) 1 else 0
            } else {
                actual_reduction / pred_reduction
            }

            ## Update the radius
            if (rho < radius_shrink) {
                Delta <- shrink_factor * Delta
            } else if (rho > radius_expand) {
                Delta <- min(expand_factor * Delta, Delta_max)
            }

            ## Accept or reject the step
            if (rho > eta) {
                betas <- betas_candidate
                theta <- theta_candidate
                fit <- fit_candidate
                fit$hatvalues <- compute_hatvalues(fit, x, nobs, nvars)
                step_components_beta <- compute_step_components(theta, fit = fit, level = 0,
                                                                adjustment_function = adjustment_function,
                                                                x = x, nobs = nobs, nvars = nvars,
                                                                weights = weights, df_residual = df_residual,
                                                                d1_transformed_dispersion = d1_transformed_dispersion,
                                                                d2_transformed_dispersion = d2_transformed_dispersion)
                adjusted_grad_beta <- with(step_components_beta, grad + adjustment)
                accept_status <- "ACCEPT"
            } else {
                accept_status <- "REJECT"
            }

            if (Delta < 1e-16) {
                warning("trust region radius became too small")
                failed <- TRUE
                break
            }

            step_components_zeta <- compute_step_components(theta, fit = fit, level = 1,
                                                            adjustment_function = adjustment_function,
                                                            x = x, nobs = nobs, nvars = nvars,
                                                            weights = weights, df_residual = df_residual,
                                                            d1_transformed_dispersion = d1_transformed_dispersion,
                                                            d2_transformed_dispersion = d2_transformed_dispersion)
            adjusted_grad_zeta <- if (no_dispersion) NA_real_ else with(step_components_zeta, grad + adjustment)

            step_zeta <- if (no_dispersion || df_residual < 1 ||
                             step_components_zeta$failed_inversion ||
                             step_components_zeta$failed_adjustment) {
                             NA_real_
                         } else {
                             as.vector(adjusted_grad_zeta * step_components_zeta$inverse_info)
                         }

            ## Newton step for the dispersion, halved until the dispersion is positive
            if (!is.na(step_zeta)) {
                transformed_dispersion_previous <- transformed_dispersion
                step_factor <- 0L
                valid_dispersion <- FALSE
                while (step_factor <= control$max_step_factor) {
                    transformed_dispersion <- transformed_dispersion_previous + 2^(-step_factor) * step_zeta
                    dispersion_new <- eval(control$inverseTrans)
                    if (valid_dispersion <- is.finite(dispersion_new) && dispersion_new > 0) break
                    step_factor <- step_factor + 1L
                }
                if (!valid_dispersion) {
                    transformed_dispersion <- transformed_dispersion_previous
                    warning("failed to update the dispersion parameter")
                    failed <- TRUE
                    break
                }

                fit_new <- compute_fit(pars = c(betas, dispersion_new),
                                       y = y,
                                       x = x,
                                       weights = weights,
                                       offset = offset,
                                       family = family,
                                       fixed_totals = fixed_totals,
                                       row_totals = row_totals,
                                       no_dispersion = no_dispersion,
                                       nobs = nobs,
                                       nvars = nvars,
                                       keep = keep,
                                       need_qr = FALSE,
                                       need_hatvalues = FALSE,
                                       need_inverse = FALSE)

                ## The information about the regression parameters scales
                ## with the dispersion, and the hat values do not depend on it
                scale <- dispersion_new / dispersion
                fit_new$info_beta <- fit$info_beta / scale
                fit_new$inverse_info_beta <- fit$inverse_info_beta * scale
                fit_new$hatvalues <- fit$hatvalues
                fit_new$qr_decomposition <- fit$qr_decomposition
                fit_new$R_matrix <- fit$R_matrix
                fit <- fit_new
                dispersion <- dispersion_new
                theta <- c(betas, dispersion)
                step_components_beta <- compute_step_components(theta, fit = fit, level = 0,
                                                                adjustment_function = adjustment_function,
                                                                x = x, nobs = nobs, nvars = nvars,
                                                                weights = weights, df_residual = df_residual,
                                                                d1_transformed_dispersion = d1_transformed_dispersion,
                                                                d2_transformed_dispersion = d2_transformed_dispersion)
                adjusted_grad_beta <- with(step_components_beta, grad + adjustment)
                step_components_zeta <- compute_step_components(theta, fit = fit, level = 1,
                                                                adjustment_function = adjustment_function,
                                                                x = x, nobs = nobs, nvars = nvars,
                                                                weights = weights, df_residual = df_residual,
                                                                d1_transformed_dispersion = d1_transformed_dispersion,
                                                                d2_transformed_dispersion = d2_transformed_dispersion)
                adjusted_grad_zeta <- with(step_components_zeta, grad + adjustment)
            }

            ## The Fisher scoring step at the current point, used for the
            ## convergence criterion as in fit_fisher_scoring(); the trust-region
            ## step is bounded by the radius and cannot be used for that
            step_beta <- fit$dispersion * backsolve(fit$R_matrix,
                                                    backsolve(fit$R_matrix, adjusted_grad_beta, transpose = TRUE))

            if (control$trace) {
                cat(sprintf("Iter %3d: ||r_adj|| = %.4e  ||p|| = %.4e  |dispersion step| = %.3e  Delta = %.3e  rho = %6.3f  %s%s\n",
                            iter, sqrt(r_adj_sq), sqrt(sum(step$p^2)),
                            if (is.na(step_zeta)) 0 else abs(step_zeta),
                            Delta, rho, accept_status,
                            if (step$on_boundary) " [boundary]" else ""))
            }

            ## Convergence criteria
            failed <- step_components_beta$failed_inversion || step_components_beta$failed_adjustment
            beta_converged <- isTRUE(max(abs(step_beta)) < control$epsilon)
            zeta_converged <- is.na(step_zeta) || abs(step_zeta) < control$epsilon
            if (failed || (beta_converged && zeta_converged)) {
                break
            }
        }
    }

    list(betas = betas,
         dispersion = dispersion,
         transformed_dispersion = transformed_dispersion,
         adjusted_grad_beta = adjusted_grad_beta,
         adjusted_grad_zeta = adjusted_grad_zeta,
         step_components_beta = step_components_beta,
         step_components_zeta = step_components_zeta,
         iter = iter,
         failed = failed)
}

#' Compute the hat values from the decomposition in a fit
#'
#' @param fit Object of class "brglmFit_quantities", computed with need_qr = TRUE
#' @param x Design matrix (dense, or a CsparseMatrix)
#' @param nobs Number of observations
#' @param nvars Number of variables
#'
#' @return The hat values, as computed by compute_fit
compute_hatvalues <- function(fit, x, nobs, nvars) {
    if (inherits(x, "sparseMatrix")) {
        Z <- backsolve(fit$R_matrix, t(as.matrix(x)), transpose = TRUE)
        fit$working_weights * .colSums(Z * Z, nvars, nobs, TRUE)
    } else {
        Q_matrix <- qr.Q(fit$qr_decomposition)
        .rowSums(Q_matrix * Q_matrix, nobs, nvars, TRUE)
    }
}

#' Solve the trust-region subproblem using preconditioned CG-Steihaug
#'
#' Approximately minimises -r_adj' p + 0.5 p' F p subject to ||p|| <= Delta,
#' using conjugate gradients with a Jacobi (diagonal) preconditioner. The
#' iteration stops at the boundary if a step leaves the trust region or a
#' direction of non-positive curvature is found, and otherwise when the
#' preconditioned residual has been reduced by the factor tol.
#'
#' @param r_adj Adjusted score vector
#' @param F_info Expected information matrix; can be NULL if Fmatvec is given
#' @param Fmatvec NULL or a function returning F %*% d for a vector d
#' @param diag_F NULL or the diagonal of F; computed from F_info if NULL
#' @param Delta Trust-region radius
#' @param tol Relative tolerance for the preconditioned residual
#' @param maxiter Maximum number of CG iterations
#'
#' @return List with the step p, whether it is on the boundary of the
#'     trust region (on_boundary), and the number of CG iterations (cg_iter)
#'
#' @references
#'
#' Steihaug T (1983). The conjugate gradient method and trust regions in
#' large scale optimization. *SIAM Journal on Numerical Analysis*, **20**,
#' 626-637.
#'
#' Nocedal J, Wright S J (2006). *Numerical Optimization*, 2nd edition.
#' Springer. Algorithms 5.3 and 7.2.
cg_steihaug_pcg <- function(r_adj, F_info, Fmatvec = NULL, diag_F = NULL, Delta, tol = 0.1, maxiter = 50) {
    # Jacobi preconditioner M = diag(F), guarding against near-zero entries
    d_F <- if (!is.null(diag_F)) diag_F else diag(F_info)
    d_F <- pmax(d_F, .Machine$double.eps * max(d_F))
    Minv <- function(v) v / d_F

    z <- numeric(length(r_adj))
    r <- r_adj
    y <- Minv(r)
    d <- y
    ry <- sum(r * y)
    ry0 <- ry

    for (j in seq_len(maxiter)) {
        Fd <- if (!is.null(Fmatvec)) Fmatvec(d) else drop(F_info %*% d)
        dFd <- sum(d * Fd)

        # Non-positive curvature: go to the boundary along d
        if (dFd <= 0) {
            tau <- find_boundary_step(z, d, Delta)
            return(list(p = z + tau * d, on_boundary = TRUE, cg_iter = j))
        }

        alpha <- ry / dFd
        z_new <- z + alpha * d

        # The step leaves the trust region: stop at the boundary
        if (sqrt(sum(z_new^2)) >= Delta) {
            tau <- find_boundary_step(z, d, Delta)
            return(list(p = z + tau * d, on_boundary = TRUE, cg_iter = j))
        }

        z <- z_new
        r <- r - alpha * Fd
        y <- Minv(r)
        ry_new <- sum(r * y)

        # Convergence of the preconditioned residual, relative to the start
        if (sqrt(abs(ry_new)) < tol * sqrt(abs(ry0)))
            return(list(p = z, on_boundary = FALSE, cg_iter = j))

        d <- y + (ry_new / ry) * d
        ry <- ry_new
    }
    list(p = z, on_boundary = FALSE, cg_iter = maxiter)
}

#' Find the step length to the trust-region boundary
#'
#' @param z Current point, with ||z|| <= Delta
#' @param d Search direction
#' @param Delta Trust-region radius
#'
#' @return The smallest tau > 0 with ||z + tau d|| = Delta, or 0 if there is none
find_boundary_step <- function(z, d, Delta) {
    a <- sum(d^2)
    b <- 2 * sum(z * d)
    cc <- sum(z^2) - Delta^2
    disc <- b^2 - 4 * a * cc
    if (disc < 0) return(0)
    t1 <- (-b + sqrt(disc)) / (2 * a)
    t2 <- (-b - sqrt(disc)) / (2 * a)
    pos <- c(t1, t2)[c(t1, t2) > 0]
    if (length(pos) == 0) 0 else min(pos)
}
