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
