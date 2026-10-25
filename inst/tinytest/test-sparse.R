## Fits using a sparse representation of the model matrix agree with
## the dense ones
set.seed(123)
n <- 600
dat <- data.frame(f1 = factor(sample(15, n, replace = TRUE)),
                  f2 = factor(sample(12, n, replace = TRUE)),
                  z = rnorm(n))
X <- model.matrix(~ f1 + f2 + z, data = dat)
eta <- drop(X %*% rnorm(ncol(X), 0, 0.3))
dat$y_bin <- rbinom(n, 1, plogis(eta))
dat$y_gam <- rgamma(n, shape = 3, rate = 3 / exp(eta))

## 27 columns, about 15% non-zero entries, so sparse = "auto" uses the
## sparse representation
for (type in c("AS_mixed", "AS_mean", "AS_median", "MPL_Jeffreys", "ML", "correction")) {
    fit_dense <- glm(y_bin ~ f1 + f2 + z, family = binomial(), data = dat,
                     method = "brglmFit", type = type, sparse = "never")
    fit_sparse <- update(fit_dense, sparse = "always")
    fit_auto <- update(fit_dense, sparse = "auto")
    expect_equal(coef(fit_sparse), coef(fit_dense), tolerance = 1e-06)
    expect_equal(vcov(fit_sparse), vcov(fit_dense), tolerance = 1e-06)
    expect_identical(coef(fit_auto), coef(fit_sparse))
}

for (type in c("AS_mixed", "AS_mean", "AS_median", "ML")) {
    fit_dense <- glm(y_gam ~ f1 + f2 + z, family = Gamma("log"), data = dat,
                     method = "brglmFit", type = type, sparse = "never")
    fit_sparse <- update(fit_dense, sparse = "always")
    expect_equal(coef(fit_sparse, model = "full"), coef(fit_dense, model = "full"),
                 tolerance = 1e-06)
}

## A sparse model matrix supplied to brglmFit is used as such with
## sparse = "auto", and densified with sparse = "never"
X_sparse <- Matrix::Matrix(X, sparse = TRUE)
fit_sparse <- brglmFit(x = X, y = dat$y_bin, family = binomial(),
                       control = list(sparse = "always"))
fit_input <- brglmFit(x = X_sparse, y = dat$y_bin, family = binomial())
fit_never <- brglmFit(x = X_sparse, y = dat$y_bin, family = binomial(),
                      control = list(sparse = "never"))
fit_dense <- brglmFit(x = X, y = dat$y_bin, family = binomial(),
                      control = list(sparse = "never"))
expect_identical(fit_input$coefficients, fit_sparse$coefficients)
expect_identical(fit_never$coefficients, fit_dense$coefficients)

## Aliased columns are handled as for dense model matrices
X_aliased <- cbind(X, z2 = 2 * X[, "z"])
fit_dense <- brglmFit(x = X_aliased, y = dat$y_bin, family = binomial(),
                      control = list(sparse = "never"))
fit_sparse <- brglmFit(x = X_aliased, y = dat$y_bin, family = binomial(),
                       control = list(sparse = "always"))
expect_true(is.na(fit_sparse$coefficients["z2"]))
expect_equal(fit_sparse$coefficients, fit_dense$coefficients, tolerance = 1e-06)

## Models with few columns are not affected by sparse = "auto"
data("lizards", package = "brglm2")
liz_never <- glm(cbind(grahami, opalinus) ~ height + diameter + light + time,
                 family = binomial(), data = lizards, method = "brglmFit", sparse = "never")
liz_auto <- update(liz_never, sparse = "auto")
expect_identical(coef(liz_auto), coef(liz_never))
