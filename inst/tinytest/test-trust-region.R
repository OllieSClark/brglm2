find_boundary_step <- brglm2:::find_boundary_step

## the boundary step lands on the boundary
z <- c(0.1, -0.2, 0.3)
d <- c(1, 2, -1)
tau <- find_boundary_step(z, d, 2)
expect_true(tau > 0)
expect_equal(sqrt(sum((z + tau * d)^2)), 2)
## from the centre, the step is Delta / ||d||
expect_equal(find_boundary_step(numeric(3), d, 2), 2 / sqrt(sum(d^2)))

## The trust-region algorithm in brglmFit
data("endometrial", package = "brglm2")
endo_fit <- function(link, type, ...) {
    glm(HG ~ NV + PI + EH, family = binomial(link), data = endometrial,
        method = "brglmFit", type = type, ...)
}

## with the cauchit link, Fisher scoring diverges from the default
## starting values, and the trust region finds the root of the mean
## bias-reducing adjusted score equations
fit_tr <- endo_fit("cauchit", "AS_mean", algorithm = "trust_region")
expect_true(fit_tr$converged)
expect_equal(unname(coef(fit_tr)), c(6.638223, 3.421147, -0.09168754, -4.151625), tolerance = 1e-05)
## the adjusted score at the estimates, as computed for Fisher scoring, is zero
expect_warning(
    fit_check <- endo_fit("cauchit", "AS_mean", algorithm = "fisher_scoring",
                          start = coef(fit_tr), maxit = 0)
)
expect_true(max(abs(fit_check$grad), na.rm = TRUE) < 1e-04)

## for maximum penalized likelihood, steps are assessed with the
## penalized log-likelihood, and the trust region finds the maximum, where
## Fisher scoring also diverges
fit_tr <- endo_fit("cauchit", "MPL_Jeffreys", algorithm = "trust_region")
expect_true(fit_tr$converged)
expect_equal(unname(coef(fit_tr)), c(6.115486, 2.604345, -0.08664978, -3.799252), tolerance = 1e-05)
expect_warning(
    fit_check <- endo_fit("cauchit", "MPL_Jeffreys", algorithm = "fisher_scoring",
                          start = coef(fit_tr), maxit = 0)
)
expect_true(max(abs(fit_check$grad), na.rm = TRUE) < 1e-04)

## convergence is not reported for maximum likelihood with infinite
## estimates (NV separates the response)
expect_warning(fit_tr <- endo_fit("logit", "ML", algorithm = "trust_region"),
               pattern = "algorithm did not converge")
expect_false(fit_tr$converged)

## from poor starting values, the trust region reaches the same estimates
## as Fisher scoring from the default ones
fit_fs <- endo_fit("logit", "AS_mean", algorithm = "fisher_scoring")
fit_tr <- endo_fit("logit", "AS_mean", algorithm = "trust_region", start = c(5, -5, 5, -5))
expect_true(fit_tr$converged)
expect_equal(coef(fit_tr), coef(fit_fs), tolerance = 1e-05)

## the bias correction is always a single Fisher scoring step
fit_fs <- endo_fit("logit", "correction", algorithm = "fisher_scoring")
fit_tr <- endo_fit("logit", "correction", algorithm = "trust_region")
expect_identical(coef(fit_tr), coef(fit_fs))
expect_identical(fit_tr$control$algorithm, "fisher_scoring")

## Fisher scoring and the trust region give the same estimates
tol <- 1e-05
both <- function(fit_fun) {
    lapply(c(fs = "fisher_scoring", tr = "trust_region"), fit_fun)
}
expect_agree <- function(fits, model = "mean") {
    expect_true(fits$fs$converged)
    expect_true(fits$tr$converged)
    expect_equal(coef(fits$tr, model = model), coef(fits$fs, model = model), tolerance = tol)
}

## binomial, all links and types
data("lizards", package = "brglm2")
for (link in c("logit", "probit", "cloglog", "cauchit")) {
    for (type in c("AS_mixed", "AS_mean", "AS_median", "MPL_Jeffreys", "ML")) {
        fits <- both(function(a) glm(cbind(grahami, opalinus) ~ height + diameter + light + time,
                                     family = binomial(link), data = lizards,
                                     method = "brglmFit", type = type, algorithm = a))
        expect_agree(fits)
    }
}

## Poisson
for (type in c("AS_mixed", "AS_median", "ML")) {
    fits <- both(function(a) glm(breaks ~ wool + tension, family = poisson(), data = warpbreaks,
                                 method = "brglmFit", type = type, algorithm = a))
    expect_agree(fits)
}

## Gamma, with transformations of the dispersion
data("coalition", package = "brglm2")
for (transformation in c("identity", "log", "inverse")) {
    for (type in c("AS_mixed", "AS_mean", "AS_median", "MPL_Jeffreys", "ML")) {
        fits <- both(function(a) glm(duration ~ fract + numst2, family = Gamma(), data = coalition,
                                     method = "brglmFit", type = type,
                                     transformation = transformation, algorithm = a))
        expect_agree(fits, model = "full")
    }
}

## inverse Gaussian and Gaussian
fits <- both(function(a) glm(duration ~ fract + numst2, family = inverse.gaussian("log"),
                             data = coalition, method = "brglmFit", algorithm = a))
expect_agree(fits, model = "full")
data("anorexia", package = "MASS")
fits <- both(function(a) glm(Postwt ~ Prewt + Treat + offset(Prewt), family = gaussian(),
                             data = anorexia, method = "brglmFit", algorithm = a))
expect_agree(fits, model = "full")

## aliased columns; the coefficients are small, so both iterations are
## run to a tighter tolerance than the default
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
fits <- both(function(a) glm(conc ~ lot*log(u) + I(2*log(u)), data = clotting, family = Gamma(),
                             method = "brglmFit", epsilon = 1e-10, algorithm = a))
expect_agree(fits, model = "full")
expect_true(is.na(coef(fits$tr)["I(2 * log(u))"]))

## multinomial logistic regression, through the Poisson log-linear model
## with fixed totals
data("hepatitis", package = "brglm2")
hepat <- hepatitis
hepat$type <- with(hepat, factor(1 - HCV * nonABC + HCV + 2 * nonABC))
hepat$type <- factor(hepat$type, labels = c("noDisease", "C", "nonABC"))
fits <- both(function(a) brmultinom(type ~ group * time, data = hepat, weights = counts,
                                    type = "AS_mean", algorithm = a))
expect_true(fits$tr$converged)
expect_equal(coef(fits$tr), coef(fits$fs), tolerance = tol)

## sparse model matrices
set.seed(123)
n <- 600
dat <- data.frame(f1 = factor(sample(15, n, replace = TRUE)),
                  f2 = factor(sample(12, n, replace = TRUE)),
                  z = rnorm(n))
X <- model.matrix(~ f1 + f2 + z, data = dat)
dat$y <- rbinom(n, 1, plogis(drop(X %*% rnorm(ncol(X), 0, 0.3))))
for (type in c("AS_mixed", "AS_median")) {
    fits <- both(function(a) glm(y ~ f1 + f2 + z, family = binomial(), data = dat,
                                 method = "brglmFit", type = type, sparse = "always", algorithm = a))
    expect_agree(fits)
}
