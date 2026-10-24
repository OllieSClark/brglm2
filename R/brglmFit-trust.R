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
