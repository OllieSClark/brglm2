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

# Quasi Fisher scoring iteration for brglmFit

#' Solve the adjusted score equations using quasi Fisher scoring with step halving
#'
#' @param betas Starting values for the regression parameters
#' @param dispersion Starting value for the dispersion
#' @param y Response vector
#' @param x Design matrix
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
#'     whether the iteration failed
fit_fisher_scoring <- function(betas, dispersion, y, x, weights, offset, family,
                               fixed_totals = NULL, row_totals = NULL,
                               no_dispersion = FALSE, nobs, nvars, keep, df_residual,
                               adjustment_function,
                               d1_transformed_dispersion, d2_transformed_dispersion,
                               control) {
    trace_iteration <- function() {
        if (iter %% control$trace == 0) {
            st <-  max(abs(step_beta), na.rm = TRUE)
            gr <- max(abs(adjusted_grad_beta), na.rm = TRUE)
            cat("Coefficients update:\t")
            cat("Outer/Inner iteration:\t", sprintf("%03d", iter), "/", sprintf("%03d", step_factor), "\n", sep = "")
            if (!no_dispersion) {
                st <- abs(step_zeta)
                gr <- abs(adjusted_grad_zeta)
                cat("Dispersion update:\t")
                cat("Outer iteration:\t", sprintf("%03d", iter), "\n")
            }
            cat("max |step|:", format(round(st, 6), nsmall = 6, scientific = FALSE), "\t",
                "max |gradient|:", format(round(gr, 6), nsmall = 6, scientific = FALSE), "\n")
        }
    }

    ## Evaluate at the starting values
    theta <- c(betas, dispersion)
    transformed_dispersion <- eval(control$Trans)

    ## Compute fit
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
                       need_hatvalues = TRUE)

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
    adjusted_grad_beta <- with(step_components_beta, {
        grad + adjustment
    })
    step_beta <- drop(step_components_beta$inverse_info %*% adjusted_grad_beta)
    ## Dispersion quantities
    if (no_dispersion) {
        adjusted_grad_zeta <- step_zeta <- NA_real_
    } else {
        if (step_components_zeta$failed_inversion) {
            warning("failed to invert the information matrix")
        }
        if (step_components_zeta$failed_adjustment) {
            warning("failed to calculate score adjustment")
        }
        adjusted_grad_zeta <- with(step_components_zeta, {
            grad + adjustment
        })
        step_zeta <- as.vector(adjusted_grad_zeta * step_components_zeta$inverse_info)
    }

    ## Main iterations
    slowit <- control$slowit
    if (control$maxit == 0) {
        iter <- 0
        failed <- FALSE
    } else {
        ## Outer iteration
        for (iter in seq.int(control$maxit)) {
            step_factor <- 0
            testhalf <- TRUE
            ## Inner iteration
            while (testhalf & step_factor < control$max_step_factor) {
                ## store previous values
                betas0 <- betas
                dispersion0 <- dispersion
                step_beta_previous <- step_beta
                step_zeta_previous <- step_zeta

                ## Update betas
                betas <- betas + slowit * 2^(-step_factor) * step_beta

                ## Update zetas
                if (!no_dispersion & df_residual > 0) {
                    transformed_dispersion <- eval(control$Trans)
                    transformed_dispersion <- transformed_dispersion + 2^(-step_factor) * step_zeta
                    dispersion <- eval(control$inverseTrans)
                }

                ## Compute key quantities
                theta <- c(betas, dispersion)
                transformed_dispersion <- eval(control$Trans)

                ## Mean quantities
                fit <- try(compute_fit(pars = theta,
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
                                       need_hatvalues = TRUE), silent = TRUE)

                ## This is to capture qr failing and revering to previous estimates
                if (failed_adjustment_beta <- inherits(fit, "try-error")) {
                    betas <- betas0
                    dispersion <- dispersion0
                    warning("failed to calculate score adjustment")
                    break
                }
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
                if (failed_inversion_beta <- step_components_beta$failed_inversion) {
                    betas <- betas0
                    dispersion <- dispersion0
                    warning("failed to invert the information matrix")
                    break
                }
                if (failed_adjustment_beta <- step_components_beta$failed_adjustment) {
                    betas <- betas0
                    dispersion <- dispersion0
                    warning("failed to calculate score adjustment")
                    break
                }
                adjusted_grad_beta <- with(step_components_beta, grad + adjustment)
                step_beta <- drop(step_components_beta$inverse_info %*% adjusted_grad_beta)

                ## Dispersion quantities
                if (no_dispersion) {
                    adjusted_grad_zeta <- step_zeta <- NA_real_
                    failed_inversion_zeta <- failed_adjustment_zeta <- FALSE
                } else {
                    if (failed_inversion_zeta <- step_components_zeta$failed_inversion) {
                        warning("failed to invert the information matrix")
                        break
                    }
                    if (failed_adjustment_zeta <- step_components_zeta$failed_adjustment) {
                        warning("failed to calculate score adjustment")
                        break
                    }
                    adjusted_grad_zeta <- with(step_components_zeta, grad + adjustment)
                    step_zeta <- as.vector(adjusted_grad_zeta * step_components_zeta$inverse_info)
                }

                ## Convergence criteria
                linf_current <- max(abs(c(step_beta, step_zeta)), na.rm = TRUE)
                linf_previous <- max(abs(c(step_beta_previous, step_zeta_previous)), na.rm = TRUE)
                testhalf <- linf_current > linf_previous

                ## Continue inner loop
                ## if (step_factor == 0 & iter == 1)  {
                ##     testhalf <- TRUE
                ## }
                step_factor <- step_factor + 1

                ##  Trace here
                if (control$trace) {
                    trace_iteration()
                }
            }
            failed <- failed_adjustment_beta || failed_inversion_beta || failed_adjustment_zeta || failed_inversion_zeta
            if (failed || linf_current < control$epsilon) {
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
