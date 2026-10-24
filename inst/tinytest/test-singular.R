## A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
    u = c(5,10,15,20,30,40,60,80,100, 5,10,15,20,30,40,60,80,100),
    conc = c(118,58,42,35,27,25,21,19,18,69,35,26,21,18,16,13,12,12),
    lot = factor(c(rep(1, 9), rep(2, 9))))
mod <- glm(conc ~ lot*log(u) + I(2*log(u)), data = clotting, family = Gamma)

X <- model.matrix(mod)
Y <- mod$y
## brglmFit returns an error if singular.ok = TRUE
expect_error(brglm_fit(X, Y, family = Gamma(), singular.ok = FALSE),
             pattern = "singular fit encountered")
expect_true(is.na(coef(brglm_fit(X, Y, family = Gamma(), singular.ok = TRUE))["I(2 * log(u))"]))

## aliased columns of zeros (from an unused factor level) are handled
## like any other aliased column
clotting$lot3 <- factor(clotting$lot, levels = c("1", "2", "3"))
X3 <- model.matrix(~ lot3*log(u), data = clotting)
fit3 <- brglm_fit(X3, Y, family = Gamma())
fit2 <- brglm_fit(model.matrix(~ lot*log(u), data = clotting), Y, family = Gamma())
expect_true(all(is.na(coef(fit3)[c("lot33", "lot33:log(u)")])))
expect_equal(unname(coef(fit3)[!is.na(coef(fit3))]), unname(coef(fit2)), tolerance = 1e-08)
