############################################################
## Simultaneous Functional Quantile Regression (SFQR)
## Refactored end-to-end script (NO copy-paste blocks)
## - Consistent parameterization:
##   * beta = (slope in active subspace Q, intercept spline in u)
##   * optional add_x enters via theta (nb x q) and bb0(u)
## - Handles add_x = NULL safely
## - Centralizes construction of X(u) blocks
## - Centralizes CV loops and fitting
############################################################

suppressPackageStartupMessages({
  library(BPST)
  library(Matrix)
  library(MASS)
  library(tensr)
  library(CVXR)
  library(fda)
})

## -----------------------------
## 0) User settings / constants
## -----------------------------
tmpid <- 1
set.seed(tmpid)

tau <- seq(0.05, 0.95, length.out = 19)
lt  <- length(tau)

nobs <- 200
nk   <- 10
ncv  <- 5

cand1 <- c(1e-2, 1e-3, 1e-4)           # group-lasso scale candidates
cand2 <- c(1e-10, 1e-11, 1e-12)        # roughness candidates (initial + sparse)
cand3 <- cand2                         # roughness candidates (refit on active)

d <- 2
r <- 1

quant_loss <- function(u, t) 0.5 * abs(u) + (t - 0.5) * u

## ------------------------------------------
## 1) Helpers: CV folds, Simpson, safe checks
## ------------------------------------------
make_folds <- function(n, K, seed = 1) {
  set.seed(seed)
  idx <- sample.int(n)
  split(idx, rep(1:K, length.out = n))
}

simpson_weights <- function(Mobs) {
  if (Mobs %% 2 == 0) stop("Simpson rule requires odd Mobs.")
  c(1, rep(c(4, 2), (Mobs - 3) / 2), 4, 1)
}

## ---------------------------------------------
## 2) Precompute triangulation / spline machinery
## ---------------------------------------------
# Expect DataGen_I.R creates: daty, mtx, t_grid, newgrids, add_x (possibly NULL), etc.
source("DataGen_I.R")

V.est  <- readRDS("nodes.RDS")
Tr.est <- readRDS("tri.RDS")

# Bernstein basis over triangulation for (t,u)
z_full <- cbind(rep(newgrids, each = lt), rep(tau, times = length(newgrids)))
Bfull <- suppressWarnings(basis(V.est, Tr.est, d, r, z_full))
B.est <- Bfull$B
q2    <- Bfull$Q2
Kmat  <- as.matrix(Bfull$K)

# Smoothness constraints (used later to build Q22)
H  <- suppressWarnings(as.matrix(smoothness(V.est, Tr.est, d, r)))
H1 <- suppressWarnings(as.matrix(smoothness(V.est, Tr.est, d, r - 1)))

# Keep only inside-domain (t,u) points
ind_inside <- Bfull$Ind.inside
z <- z_full[ind_inside, , drop = FALSE]

# Roughness in reduced coordinates
rp_full <- t(q2) %*% Kmat %*% q2

# B-spline basis in u for intercept and varying add_x effect
knots  <- quantile(tau, seq(0, 1, length.out = nk))
norder <- 2
nb     <- length(knots) + norder - 2
basisobj <- create.bspline.basis(
  rangeval = c(min(tau), max(tau)),
  nbasis   = nb,
  norder   = norder,
  breaks   = knots
)
bb0 <- getbasismatrix(tau, basisobj)  # lt x nb

# Simpson setup for integral in t
Mobs  <- ncol(mtx)
tmin  <- min(t_grid); tmax <- max(t_grid)
hDM   <- (tmax - tmin) / Mobs
wSim  <- simpson_weights(Mobs)

# Convenience: Bernstein slices per tau
B_by_u <- lapply(seq_len(lt), function(i) B.est[z[,2] == tau[i], , drop = FALSE])

## ---------------------------------------------
## 3) Build X(u) given a basis Q (q2 or Q22)
## ---------------------------------------------
build_X_u <- function(i, idx, Q, mtx, add_x, bb0, B_by_u, hDM, wSim) {
  bb <- B_by_u[[i]]
  p  <- (hDM / 3) * mtx[idx, , drop = FALSE] %*% diag(wSim) %*% bb
  
  X <- cbind(
    p %*% Q,
    matrix(bb0[i, ], nrow = length(idx), ncol = ncol(bb0), byrow = TRUE)
  )
  
  bbtheta <- matrix(bb0[i, ], nrow = 1, ncol = ncol(bb0), byrow = TRUE)
  
  list(X = X, bbtheta = bbtheta)
}

## ----------------------------------------------------
## 4) Generic fitter: given Q and roughness matrix rp_Q
## ----------------------------------------------------
fit_sfqr <- function(y, idx, tau, Q, rp_Q, mtx, add_x, bb0, B_by_u, hDM, wSim, gamma) {
  
  y_sub <- y[idx]
  n_sub <- length(idx)
  
  has_addx <- !is.null(add_x)
  q_addx   <- if (has_addx) ncol(add_x) else 0
  
  beta  <- CVXR::Variable(ncol(Q) + ncol(bb0))
  if (has_addx) theta <- CVXR::Variable(rows = ncol(bb0), cols = q_addx)
  
  obj <- 0
  for (i in seq_along(tau)) {
    tmp <- build_X_u(i, idx, Q, mtx, add_x, bb0, B_by_u, hDM, wSim)
    X   <- tmp$X
    bbtheta <- tmp$bbtheta
    
    if (!has_addx) {
      obj <- obj + sum(quant_loss(y_sub - X %*% beta, tau[i]))
    } else {
      ax <- add_x[idx, , drop = FALSE]
      obj <- obj + sum(quant_loss(y_sub - X %*% beta - ax %*% t(bbtheta %*% theta), tau[i]))
    }
  }
  
  roughness <- quad_form(beta[1:ncol(Q)], rp_Q) * n_sub
  prob <- Problem(Minimize(obj + gamma * roughness))
  
  # Use your psolve if you have it; otherwise CVXR::solve
  fit <- psolve(prob, solver='SCS')
  out <- list(beta = fit$getValue(beta), theta = NULL)
  if (has_addx) out$theta <- fit$getValue(theta)
  out
}

## ----------------------------------------------------
## 5) CV for gamma (roughness) given Q and rp_Q
## ----------------------------------------------------
cv_sfqr_gamma <- function(y, tau, Q, rp_Q, mtx, add_x, bb0, B_by_u, hDM, wSim, cand_gamma, ncv = 5, seed = 1) {
  
  n <- length(y)
  folds <- make_folds(n, ncv, seed = seed)
  
  cvloss <- numeric(length(cand_gamma))
  
  for (lg in seq_along(cand_gamma)) {
    gamma <- cand_gamma[lg]
    loss_sum <- 0
    
    for (k in seq_along(folds)) {
      te <- folds[[k]]
      tr <- setdiff(seq_len(n), te)
      
      fit <- fit_sfqr(
        y = y, idx = tr, tau = tau,
        Q = Q, rp_Q = rp_Q,
        mtx = mtx, add_x = add_x, bb0 = bb0, B_by_u = B_by_u,
        hDM = hDM, wSim = wSim,
        gamma = gamma
      )
      
      # evaluate loss on test
      has_addx <- !is.null(add_x)
      y_te <- y[te]
      
      for (i in seq_along(tau)) {
        tmp <- build_X_u(i, te, Q, mtx, add_x, bb0, B_by_u, hDM, wSim)
        X   <- tmp$X
        bbtheta <- tmp$bbtheta
        
        if (!has_addx) {
          resid <- y_te - X %*% fit$beta
        } else {
          ax <- add_x[te, , drop = FALSE]
          resid <- y_te - X %*% fit$beta - ax %*% t(bbtheta %*% fit$theta)
        }
        loss_sum <- loss_sum + sum(quant_loss(resid, tau[i]))
      }
    }
    
    cvloss[lg] <- loss_sum
    message("gamma idx=", lg, " gamma=", signif(gamma, 3), " cvloss=", signif(loss_sum, 6))
  }
  
  cvloss
}

## ----------------------------------------------------
## 6) Initial (non-sparse) fit: choose gamma via CV
## ----------------------------------------------------
# Load precomputed CV if you already have it; otherwise run it:
tmpname <- paste0(tmpid, "_gammacv.RDS")
gammacv <- readRDS(tmpname)
optid <- which.min(gammacv)
gamma0 <- cand2[optid]
message("Initial-fit optimal gamma = ", gamma0)

fit_init <- fit_sfqr(
  y = daty, idx = seq_along(daty), tau = tau,
  Q = q2, rp_Q = rp_full,
  mtx = mtx, add_x = add_x, bb0 = bb0, B_by_u = B_by_u,
  hDM = hDM, wSim = wSim,
  gamma = gamma0
)
g1 <- fit_init$beta
saveRDS(g1, file = paste0("estg1_", tmpid, ".RDS"))

## ----------------------------------------------------
## 6) Sparse fit (adaptive group lasso on triangles)
##     - we keep your exact logic, but remove copy-paste
## ----------------------------------------------------
nbern <- choose(d + 2, 2)
ngrp  <- nrow(q2) / nbern

# Choose your integration matrix option
int_mtx <- diag(1, nbern)
# int_mtx <- matrix(c(...), nrow=6, ncol=6, byrow=TRUE) # your alternative

# Compute adaptive weights lambda_k based on initial slope coefficients
lambda <- numeric(ngrp)
for (k in 1:ngrp) {
  rows_k <- (1:nbern) + (k - 1) * nbern
  # this matches your expression
  lambda[k] <- sqrt(
    g1[1:ncol(q2)] %*% t(q2[rows_k, , drop = FALSE]) %*% int_mtx %*% q2[rows_k, , drop = FALSE] %*% g1[1:ncol(q2)]
  )
}
alpha <- -1

# Pick tuning indices (you can loop over cand1/cand2 later if needed)
l1 <- 1
l2 <- 2
scale <- cand1[l1]
gamma_sparse <- cand2[l2]

# Fit with group-lasso + roughness:
# We implement sparse penalty as a CVXR expression on beta[1:ncol(q2)] (slope part only)
beta  <- CVXR::Variable(ncol(q2) + ncol(bb0))
has_addx <- !is.null(add_x)
if (has_addx) theta <- CVXR::Variable(rows = ncol(bb0), cols = ncol(add_x))

obj <- 0
idx_all <- seq_along(daty)

for (i in seq_along(tau)) {
  tmp <- build_X_u(i, idx_all, q2, mtx, add_x, bb0, B_by_u, hDM, wSim)
  X <- tmp$X
  bbtheta <- tmp$bbtheta
  if (!has_addx) {
    obj <- obj + sum(quant_loss(daty - X %*% beta, tau[i]))
  } else {
    obj <- obj + sum(quant_loss(daty - X %*% beta - add_x %*% t(bbtheta %*% theta), tau[i]))
  }
}
# adaptive group lasso on slope coefficients only (beta[1:ncol(q2)])
grp_lasso <- 0
for (k in 1:ngrp) {
  rows_k <- (1:nbern) + (k - 1) * nbern
  grp_lasso <- grp_lasso + (lambda[k]^alpha) *
    norm(int_mtx %*% q2[rows_k, , drop = FALSE] %*% beta[1:ncol(q2)], type = "2")
}

roughness <- quad_form(beta[1:ncol(q2)], rp_full) * length(daty)

prob <- Problem(Minimize(obj + scale * nobs * grp_lasso + gamma_sparse * roughness))

fit <- psolve(prob,solver="SCS")

out <- list(beta = fit$getValue(beta), theta = NULL)
if (has_addx) out$theta <- fit$getValue(theta)

tmpg <- out$beta
saveRDS(tmpg, file = paste0("estg_", tmpid, ".RDS"))

## ----------------------------------------------------
## 8) Identify active triangles & construct Q22, rp1
## ----------------------------------------------------
b1 <- q2 %*% tmpg[1:ncol(q2)]
norm_tri <- sapply(1:ngrp, function(k) {
  rows_k <- (1:nbern) + (k - 1) * nbern
  sqrt(sum(b1[rows_k]^2))
})
message("# of zero triangles (<=1e-3): ", sum(norm_tri <= 1e-3))

inactive1 <- rep(norm_tri, each = nbern)

Q22 <- NULL
rp1 <- NULL
z1  <- NULL

if (sum(norm_tri >= 1e-3) > 0) {
  if (sum(inactive1 <= 1e-3) > 0) {
    indicator1 <- apply(abs(H1[, inactive1 <= 1e-3, drop = FALSE]), 1, sum)
    indicator  <- apply(abs(H[,  inactive1 <= 1e-3, drop = FALSE]), 1, sum)
    
    tmpH  <- H[indicator == 0, , drop = FALSE]
    tmpHH <- rbind(tmpH, H1[indicator1 != 0, , drop = FALSE])
    
    # additionally constrain inactive Bernstein coefficients to 0
    for (j in which(inactive1 <= 1e-3)) {
      e <- rep(0, length(inactive1)); e[j] <- 1
      tmpHH <- rbind(tmpHH, e)
    }
    
    Q22 <- qrH(tmpHH)
    rp1 <- t(Q22) %*% Kmat %*% Q22
    z1  <- sum(norm_tri >= 1e-3)
  } else {
    Q22 <- q2
    rp1 <- t(Q22) %*% Kmat %*% Q22
    z1  <- 0
  }
  saveRDS(Q22, file = paste0("Q22_", tmpid, ".RDS"))
}

## ----------------------------------------------------
## 9) CV for refit roughness gamma on active region (Q22)
## ----------------------------------------------------
if (!is.null(z1)) {
  out_cand3 <- cv_sfqr_gamma(
    y = daty, tau = tau, Q = Q22, rp_Q = rp1,
    mtx = mtx, add_x = add_x, bb0 = bb0, B_by_u = B_by_u,
    hDM = hDM, wSim = wSim,
    cand_gamma = cand3, ncv = ncv, seed = tmpid
  )
  saveRDS(out_cand3, file = paste0(tmpid, "_refit_gammacv.RDS"))
  
  l3 <- which.min(out_cand3)
  gamma_refit <- cand3[l3]
  message("Refit optimal gamma = ", gamma_refit)
  
  fit_refit <- fit_sfqr(
    y = daty, idx = seq_along(daty), tau = tau,
    Q = Q22, rp_Q = rp1,
    mtx = mtx, add_x = add_x, bb0 = bb0, B_by_u = B_by_u,
    hDM = hDM, wSim = wSim,
    gamma = gamma_refit
  )
  
  g_final <- fit_refit$beta
  saveRDS(g_final, file = paste0("fg_", tmpid, ".RDS"))
  
  if (!is.null(add_x)) {
    gg_final <- fit_refit$theta
    saveRDS(gg_final, file = paste0("fgg_", tmpid, ".RDS"))
  }
}

## ----------------------------------------------------
## 10) Plot example (same as yours)
## ----------------------------------------------------
# plot true slope (assumes rhos(t_grid) exists)
plot(t_grid, rhos(t_grid), col = "red", type = "l")

if (!is.null(z1)) {
  # example: plot estimated slope at tau[9]
  quantile_index <- 9
  lines(
    t_grid,
    B_by_u[[quantile_index]] %*% Q22 %*% g_final[1:ncol(Q22)],
    type = "l"
  )
} else {
  # if no active region refit, plot sparse estimate using q2
  quantile_index <- 9
  lines(
    t_grid,
    B_by_u[[quantile_index]] %*% q2 %*% tmpg[1:ncol(q2)],
    type = "l"
  )
}
