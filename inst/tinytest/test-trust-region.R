cg_steihaug_pcg <- brglm2:::cg_steihaug_pcg
find_boundary_step <- brglm2:::find_boundary_step

## A small positive definite information matrix and score
set.seed(321)
A <- matrix(rnorm(25), 5, 5)
F_info <- crossprod(A) + diag(5)
r_adj <- rnorm(5)
p_newton <- solve(F_info, r_adj)
model <- function(p) -sum(r_adj * p) + 0.5 * sum(p * (F_info %*% p))

## the boundary step lands on the boundary
z <- c(0.1, -0.2, 0.3)
d <- c(1, 2, -1)
tau <- find_boundary_step(z, d, 2)
expect_true(tau > 0)
expect_equal(sqrt(sum((z + tau * d)^2)), 2)
## from the centre, the step is Delta / ||d||
expect_equal(find_boundary_step(numeric(3), d, 2), 2 / sqrt(sum(d^2)))

## with a large radius and a tight tolerance, the solution is the Newton step
step <- cg_steihaug_pcg(r_adj, F_info, Delta = 100, tol = 1e-12, maxiter = 50)
expect_false(step$on_boundary)
expect_equal(step$p, p_newton, tolerance = 1e-08)

## with a small radius, the step is on the boundary and reduces the model
Delta <- 0.5 * sqrt(sum(p_newton^2))
step <- cg_steihaug_pcg(r_adj, F_info, Delta = Delta, tol = 1e-12, maxiter = 50)
expect_true(step$on_boundary)
expect_equal(sqrt(sum(step$p^2)), Delta)
expect_true(model(step$p) < 0)

## with non-positive curvature, the step goes to the boundary
F_indefinite <- diag(c(2, -1, 3))
step <- cg_steihaug_pcg(c(1, 1, 1), F_indefinite, Delta = 1.5, tol = 1e-12, maxiter = 50)
expect_true(step$on_boundary)
expect_equal(sqrt(sum(step$p^2)), 1.5)

## for a diagonal information matrix, the Jacobi preconditioner gives the
## solution in one iteration, however badly scaled the matrix is
F_diagonal <- diag(c(1, 1e2, 1e4, 1e6))
step <- cg_steihaug_pcg(c(1, 1, 1, 1), F_diagonal, Delta = 100, tol = 1e-12, maxiter = 50)
expect_equal(step$cg_iter, 1)
expect_equal(step$p, c(1, 1e-2, 1e-4, 1e-6))

## a function for F %*% d and the diagonal of F give the same step as F
step_matrix <- cg_steihaug_pcg(r_adj, F_info, Delta = Delta, tol = 1e-12, maxiter = 50)
step_function <- cg_steihaug_pcg(r_adj, NULL, Fmatvec = function(d) drop(F_info %*% d),
                                 diag_F = diag(F_info), Delta = Delta, tol = 1e-12, maxiter = 50)
expect_equal(step_function, step_matrix)

## the number of CG iterations is capped by maxiter
step <- cg_steihaug_pcg(r_adj, F_info, Delta = 100, tol = 1e-12, maxiter = 2)
expect_equal(step$cg_iter, 2)
expect_false(step$on_boundary)

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

## convergence is not reported when the adjusted score is not zero
expect_warning(fit_tr <- endo_fit("cauchit", "MPL_Jeffreys", algorithm = "trust_region"),
               pattern = "algorithm did not converge")
expect_false(fit_tr$converged)

## nor for maximum likelihood with infinite estimates (NV separates the response)
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
