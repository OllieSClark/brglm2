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
#' Each iteration takes the Fisher scoring step for the regression
#' parameters if it is inside the trust region, and otherwise a step to
#' its boundary, accepts or rejects the step according to the ratio of
#' the actual to the predicted reduction in a merit function, and updates
#' the radius (Nocedal & Wright, 2006, Algorithm 4.1). After an accepted
#' step, the dispersion is updated with a Newton step, halved until the
#' dispersion is positive.
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
#' Powell M J D (1970). A hybrid method for nonlinear equations. In
#' Rabinowitz P (ed), *Numerical Methods for Nonlinear Algebraic
#' Equations*, 87-114. Gordon and Breach.
#'
#' More J J (1978). The Levenberg-Marquardt algorithm: implementation and
#' theory. In Watson G A (ed), *Numerical Analysis*, Lecture Notes in
#' Mathematics 630, 105-116. Springer.
#'
#' Deuflhard P (2011). *Newton Methods for Nonlinear Problems*. Springer.
#'
#' Firth D (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika*, **80**, 27-38.
fit_trust_region <- function(betas, dispersion, y, x, weights, offset, family,
                             fixed_totals = NULL, row_totals = NULL,
                             no_dispersion = FALSE, nobs, nvars, keep, df_residual,
                             adjustment_function,
                             d1_transformed_dispersion, d2_transformed_dispersion,
                             control) {
    ## Evaluate at the starting values
    theta <- c(betas, dispersion)
    transformed_dispersion <- eval(control$Trans)

    ## Gradient, inverse information and adjustment for the regression parameters
    beta_components <- function(theta, fit) {
        compute_step_components(theta, fit = fit, level = 0,
                                adjustment_function = adjustment_function,
                                x = x, nobs = nobs, nvars = nvars,
                                weights = weights, df_residual = df_residual,
                                d1_transformed_dispersion = d1_transformed_dispersion,
                                d2_transformed_dispersion = d2_transformed_dispersion)
    }

    ## The inverse of the information about the regression parameters is
    ## only needed by the median adjustment
    need_inverse <- control$type == "AS_median"

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

    step_components_beta <- beta_components(theta, fit)
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

    ## Newton step for the (transformed) dispersion; it is updated only
    ## after accepted steps for the regression parameters
    dispersion_step <- function(components, adjusted_grad) {
        if (no_dispersion || df_residual < 1 ||
            components$failed_inversion || components$failed_adjustment) {
            NA_real_
        } else {
            as.vector(adjusted_grad * components$inverse_info)
        }
    }
    step_zeta <- dispersion_step(step_components_zeta, adjusted_grad_zeta)

    ## The Fisher scoring step F^{-1} r for the regression parameters, from
    ## the decomposition of the information F = R'R / dispersion
    newton_step <- function(fit, r) {
        fit$dispersion * backsolve(fit$R_matrix, backsolve(fit$R_matrix, r, transpose = TRUE))
    }
    step_beta <- newton_step(fit, adjusted_grad_beta)

    ## Steps are assessed by the reduction in a merit function. For
    ## maximum likelihood and maximum penalized likelihood, the adjusted
    ## score is the gradient of the (penalized) log-likelihood, so its
    ## negative is used. With the canonical link, the mean bias-reducing
    ## adjustment for the regression parameters is that of the
    ## Jeffreys-prior penalty with a = 1/2 (Firth, 1993), so mean (and
    ## mixed) bias reduction is assessed in the same way. Otherwise, the
    ## natural level function r' F^{-1} r is used, with the information F
    ## held at the current point (Deuflhard, 2011); this is the squared
    ## length of the Fisher scoring step, which the iteration reduces even
    ## when the norm of the adjusted score does not decrease
    canonical_links <- c(binomial = "logit", poisson = "log", Gamma = "inverse",
                         gaussian = "identity", inverse.gaussian = "1/mu^2")
    is_canonical <- identical(unname(canonical_links[family$family]), family$link)
    penalty_power <- switch(control$type,
                            "MPL_Jeffreys" = control$a,
                            "AS_mean" = ,
                            "AS_mixed" = if (is_canonical) 1/2 else NA_real_,
                            "ML" = 0,
                            NA_real_)
    use_objective <- !is.na(penalty_power) && is.null(fixed_totals)
    objective <- function(fit) {
        penalty <- if (penalty_power != 0) penalty_power * 2 * sum(log(abs(diag(fit$R_matrix)))) else 0
        sum(fit$deviance_residuals) / (2 * fit$dispersion) - penalty
    }

    ## The trust-region step:
    ## * with the objective, the dogleg step (Powell, 1970) for the model
    ##   -r'p + p'F p / 2, with the radius bounding the scaled length ||D p||,
    ##   D the square roots of the diagonal of the information, kept
    ##   non-decreasing (More, 1978): the Fisher scoring step if it is inside
    ##   the trust region, and otherwise the point on the boundary along the
    ##   path from the minimizer of the model along the scaled score to the
    ##   Fisher scoring step. Every direction along the path reduces the
    ##   objective for small steps
    ## * with the natural level function, which is reduced by short enough
    ##   steps along the Fisher scoring direction but not necessarily along
    ##   other directions, the Fisher scoring step shortened to the
    ##   boundary, with the radius bounding the length in the norm defined by
    ##   the information, ||p||_F = sqrt(p' F p)
    ## In both cases the step does not depend on the scale of the
    ## covariates, and the initial radius is the length of the first Fisher
    ## scoring step
    info_times <- function(fit, v) drop(crossprod(fit$R_matrix, fit$R_matrix %*% v)) / fit$dispersion
    scaling <- sqrt(colSums(fit$R_matrix^2) / fit$dispersion)
    newton_length <- if (use_objective) {
                         sqrt(sum((scaling * step_beta)^2))
                     } else {
                         sqrt(sum(adjusted_grad_beta * step_beta))
                     }
    Delta <- if (isTRUE(newton_length > 0)) newton_length else 1
    Delta_max <- 1e10
    radius_shrink <- 0.25
    radius_expand <- 0.75
    shrink_factor <- 0.25
    expand_factor <- 2
    eta <- 0.1

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
            ## The natural level function r' F^{-1} r at the current point
            level <- sum(adjusted_grad_beta * step_beta)
            if (use_objective) {
                merit <- objective(fit)
                ## Dogleg step in the scaled variables q = D p
                scaling <- pmax(scaling, sqrt(colSums(fit$R_matrix^2) / fit$dispersion))
                q_newton <- scaling * step_beta
                newton_length <- sqrt(sum(q_newton^2))
            } else {
                newton_length <- sqrt(level)
            }
            if (!is.finite(newton_length)) {
                warning("failed to invert the information matrix")
                failed <- TRUE
                break
            }
            if (newton_length <= Delta) {
                step <- step_beta
                step_length <- newton_length
            } else if (use_objective) {
                g <- adjusted_grad_beta / scaling
                gHg <- sum((g / scaling) * info_times(fit, g / scaling))
                q_cauchy <- (sum(g^2) / gHg) * g
                cauchy_length <- sqrt(sum(q_cauchy^2))
                q <- if (!is.finite(cauchy_length) || cauchy_length >= Delta) {
                         (Delta / sqrt(sum(g^2))) * g
                     } else {
                         q_cauchy + find_boundary_step(q_cauchy, q_newton - q_cauchy, Delta) * (q_newton - q_cauchy)
                     }
                step <- q / scaling
                step_length <- Delta
            } else {
                step <- (Delta / newton_length) * step_beta
                step_length <- Delta
            }
            on_boundary <- step_length >= 0.99 * Delta

            ## Predicted reduction in the merit, from the model -r'p + p'F p / 2
            ## for the objective; the model for the natural level function,
            ## from the linear model r - F p, is twice that plus a constant.
            ## Close to the solution, the reduction in the objective is lost in
            ## rounding error, and the natural level function is used instead
            model_reduction <- sum(adjusted_grad_beta * step) - 0.5 * sum(step * info_times(fit, step))
            use_level <- !use_objective || !isTRUE(abs(model_reduction) > 1e-10 * max(1, abs(merit)))
            pred_reduction <- if (use_level) 2 * model_reduction else model_reduction

            betas_candidate <- betas + step
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
                Delta <- shrink_factor * min(Delta, step_length)
                if (Delta < 1e-16) {
                    warning("trust region radius became too small")
                    failed <- TRUE
                    break
                }
                next
            }

            ## Actual reduction in the merit. For the natural level function,
            ## as in quasi Fisher scoring, the hat values in the adjustment
            ## are held at the current point
            if (!use_level) {
                actual_reduction <- merit - objective(fit_candidate)
            } else {
                fit_candidate$hatvalues <- fit$hatvalues
                adjustment_candidate <- adjustment_function(theta_candidate, fit = fit_candidate, level = 0,
                                                            x, nobs, nvars, weights)
                r_adj_candidate <- fit_candidate$grad_beta + adjustment_candidate
                level_candidate <- fit$dispersion * sum(backsolve(fit$R_matrix, r_adj_candidate, transpose = TRUE)^2)
                actual_reduction <- level - level_candidate
            }

            ## A non-finite reduction rejects the step; reductions within
            ## rounding error accept it unless the merit clearly increases
            noise <- 100 * .Machine$double.eps * (if (use_level) level else max(1, abs(merit)))
            rho <- if (isTRUE(max(abs(step_beta)) < control$epsilon)) {
                ## The Fisher scoring step is already within the tolerance,
                ## and the reductions are rounding error
                1
            } else if (!is.finite(actual_reduction)) {
                -Inf
            } else if (!is.finite(pred_reduction) || abs(pred_reduction) <= noise) {
                if (actual_reduction > -noise) 1 else 0
            } else {
                actual_reduction / pred_reduction
            }

            ## Update the radius: shrink relative to the length of the step,
            ## which can be much shorter than the radius, and expand only if
            ## the step reached the boundary (Lin & Moré, 1999; Conn et al.,
            ## 2000, Section 10.5.2)
            if (rho < radius_shrink) {
                Delta <- shrink_factor * min(Delta, step_length)
            } else if (rho > radius_expand && step_length >= 0.99 * Delta) {
                Delta <- min(expand_factor * Delta, Delta_max)
            }

            ## Accept or reject the step
            accepted <- rho > eta
            if (accepted) {
                betas <- betas_candidate
                theta <- theta_candidate
                fit <- fit_candidate
                fit$hatvalues <- compute_hatvalues(fit, x, nobs, nvars)
                step_components_beta <- beta_components(theta, fit)
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

            ## Newton step for the dispersion after an accepted step, halved
            ## until the dispersion is positive
            if (accepted) {
                step_components_zeta <- compute_step_components(theta, fit = fit, level = 1,
                                                                adjustment_function = adjustment_function,
                                                                x = x, nobs = nobs, nvars = nvars,
                                                                weights = weights, df_residual = df_residual,
                                                                d1_transformed_dispersion = d1_transformed_dispersion,
                                                                d2_transformed_dispersion = d2_transformed_dispersion)
                adjusted_grad_zeta <- if (no_dispersion) NA_real_ else with(step_components_zeta, grad + adjustment)
                step_zeta <- dispersion_step(step_components_zeta, adjusted_grad_zeta)
            }
            if (accepted && !is.na(step_zeta)) {
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
                if (!is.null(fit$inverse_info_beta)) {
                    fit_new$inverse_info_beta <- fit$inverse_info_beta * scale
                }
                fit_new$hatvalues <- fit$hatvalues
                fit_new$qr_decomposition <- fit$qr_decomposition
                fit_new$R_matrix <- fit$R_matrix
                fit <- fit_new
                dispersion <- dispersion_new
                theta <- c(betas, dispersion)
                ## The adjustments for the regression parameters do not
                ## depend on the dispersion; only the score does
                step_components_beta$grad <- fit$grad_beta
                adjusted_grad_beta <- with(step_components_beta, grad + adjustment)
            }

            ## The Fisher scoring step at the current point, for the next
            ## iteration and for the convergence criterion, as in
            ## fit_fisher_scoring(); the trust-region step is bounded by the
            ## radius and cannot be used for that
            step_beta <- newton_step(fit, adjusted_grad_beta)

            if (control$trace) {
                cat(sprintf("Iter %3d: ||r_adj|| = %.4e  ||p||_F = %.4e  |dispersion step| = %.3e  Delta = %.3e  rho = %6.3f  %s%s\n",
                            iter, sqrt(r_adj_sq), step_length,
                            if (is.na(step_zeta)) 0 else abs(step_zeta),
                            Delta, rho, accept_status,
                            if (on_boundary) " [boundary]" else ""))
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

    ## Dispersion components at the final values
    step_components_zeta <- compute_step_components(theta, fit = fit, level = 1,
                                                    adjustment_function = adjustment_function,
                                                    x = x, nobs = nobs, nvars = nvars,
                                                    weights = weights, df_residual = df_residual,
                                                    d1_transformed_dispersion = d1_transformed_dispersion,
                                                    d2_transformed_dispersion = d2_transformed_dispersion)
    adjusted_grad_zeta <- if (no_dispersion) NA_real_ else with(step_components_zeta, grad + adjustment)

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
    pivot <- if (is.null(fit$qr_decomposition)) seq_len(nvars) else fit$qr_decomposition$pivot
    .hatvalues(x, fit$working_weights, fit$R_matrix, pivot, nobs, nvars)
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
