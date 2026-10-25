data("endometrial", package = "brglm2")
endo_fit <- function(link, type, ...) {
    glm(HG ~ NV + PI + EH, family = binomial(link), data = endometrial,
        method = "brglmFit", type = type, ...)
}

## with the cauchit link, the iteration used to diverge to estimates of
## the order of 1e295 and report convergence; it now reaches the root of
## the adjusted score equations
fit <- endo_fit("cauchit", "AS_mean")
expect_true(fit$converged)
expect_equal(unname(coef(fit)), c(6.638223, 3.421147, -0.09168754, -4.151625), tolerance = 1e-05)
expect_true(max(abs(fit$grad), na.rm = TRUE) < 1e-04)

fit <- endo_fit("cauchit", "MPL_Jeffreys")
expect_true(fit$converged)
expect_equal(unname(coef(fit)), c(6.115486, 2.604345, -0.08664978, -3.799252), tolerance = 1e-05)
expect_true(max(abs(fit$grad), na.rm = TRUE) < 1e-04)

## convergence is not reported when the estimates are infinite
expect_warning(fit <- endo_fit("logit", "ML"), pattern = "algorithm did not converge")
expect_false(fit$converged)

## from poor starting values, step halving from the current estimates
## reaches the same estimates as from the default starting values
for (s in 1:3) {
    set.seed(s)
    X <- cbind(1, matrix(rnorm(800), 200))
    y <- rbinom(200, 1, plogis(X %*% c(0.3, 0.5, -0.5, 0.5, -0.5)))
    set.seed(10000 + s)
    start <- rnorm(5, 0, 3)
    fit_default <- brglmFit(x = X, y = y, family = binomial())
    fit_start <- suppressWarnings(brglmFit(x = X, y = y, family = binomial(), start = start,
                                           control = list(maxit = 1000)))
    expect_true(fit_start$converged)
    expect_equal(fit_start$coefficients, fit_default$coefficients, tolerance = 1e-05)
}

## if no halved step is shorter than the current one, the full step is
## taken, which is what the Poisson log-linear model needs here
set.seed(4)
X <- cbind(1, matrix(rnorm(800), 200))
y <- rpois(200, exp(X %*% c(0.5, 0.3, -0.3, 0.3, -0.3)))
set.seed(20004)
start <- rnorm(5, 0, 1)
fit_default <- brglmFit(x = X, y = y, family = poisson())
fit_start <- brglmFit(x = X, y = y, family = poisson(), start = start, control = list(maxit = 1000))
expect_true(fit_start$converged)
expect_equal(fit_start$coefficients, fit_default$coefficients, tolerance = 1e-05)
