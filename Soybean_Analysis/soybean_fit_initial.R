## ============================================================
## Functional quantile regression w/ roughness penalty (CV tuning)
## Cleaned + modularized version of your script
## ============================================================

suppressPackageStartupMessages({
  library(stats)
  library(BPST)
  library(Matrix)
  library(MASS)
  library(tensr)
  library(CVXR)
  library(fda)
})

## ----------------------------
## 0) Global settings / helpers
## ----------------------------
tau <- seq(0.05, 0.95, length.out = 19)
lt  <- length(tau)
quant_loss <- function(u, t) {0.5 * abs(u) + (t - 0.5) * u}


# Build K-fold splits (balanced, deterministic given set.seed)
make_folds <- function(n, K = 10, seed = 1) {
  set.seed(seed)
  id <- sample.int(n)
  fold_id <- rep(seq_len(K), length.out = n)
  split(id, fold_id)
}

## ----------------------------
## 1) Load data + preprocess
## ----------------------------
tmpx    <- readRDS("soy_avg_x.RDS")
daty    <- readRDS("soy_avg_y.RDS")
add_var <- readRDS("soy_avg_add.RDS")

add_var[, 3] <- add_var[, 3] / add_var[, 2]
add_x <- add_var[, -2, drop = FALSE]
colnames(add_x) <- c("prcp", "irr_ratio")

mtx2 <- tmpx[,  30:334, drop = FALSE]  # min temp
mtx1 <- tmpx[, 396:700, drop = FALSE]  # max temp
mtx  <- (mtx1 + mtx2) / 2

n     <- nrow(mtx)
Mobs  <- ncol(mtx)

t_grid <- seq(0, 1, length.out = Mobs)

stopifnot(length(daty) == nrow(mtx), nrow(add_x) == n)
Y <- daty

## ----------------------------
## 2) Triangulation basis pieces
## ----------------------------
V.est  <- readRDS("nodes.RDS")
Tr.est <- readRDS("tri.RDS")

ntr   <- nrow(Tr.est)
d     <- 2
r     <- 1
nbern <- choose(d + 2, 2)

# (Optional) triangle "area" computation kept but not used downstream
area <- numeric(ntr)
for (i in seq_len(ntr)) {
  v.index <- Tr.est[i, ]
  a <- V.est[v.index[1], ]
  b <- V.est[v.index[2], ]
  c <- V.est[v.index[3], ]
  area[i] <- abs(a[1] * (b[2] - c[2]) +
                   b[1] * (c[2] - a[2]) +
                   c[1] * (a[2] - b[2])) / 2
}
area <- sqrt(area)

# Build (t, tau) grid for basis()
newgrids <- t_grid
z <- cbind(rep(newgrids, each = lt), rep(tau, times = length(newgrids)))

Bfull.est <- basis(V.est, Tr.est, d, r, z)
B.est     <- Bfull.est$B
q2        <- Bfull.est$Q2
Kmat      <- as.matrix(Bfull.est$K)
ind       <- Bfull.est$Ind.inside

# keep only inside points
z     <- z[ind, , drop = FALSE]
B.est <- B.est[ind, , drop = FALSE]

# sanity: z now aligned with rows of B.est
stopifnot(nrow(z) == nrow(B.est))

# group count (if you need it later)
ngrp <- nrow(q2) / nbern

## ----------------------------
## 3) Simpson rule on t-grid
## ----------------------------
# Simpson requires odd number of points (even number of sub-intervals)
if ((Mobs - 1) %% 2 != 0) {
  stop("Simpson rule needs (Mobs-1) even (i.e., Mobs odd). Current Mobs = ", Mobs)
}

hDM   <- (max(t_grid) - min(t_grid)) / (Mobs - 1)
cefDM <- c(1, rep(c(4, 2), (Mobs - 3) / 2), 4, 1)  # length Mobs
stopifnot(length(cefDM) == Mobs)

Wsimpson <- diag(cefDM, nrow = Mobs, ncol = Mobs)

## ----------------------------
## 4) B-spline basis in tau
## ----------------------------
knots  <- quantile(tau, seq(0, 1, length.out = 31))
norder <- 2
nb     <- length(knots) + norder - 2

basisobj <- create.bspline.basis(
  rangeval = c(min(tau), max(tau)),
  nbasis   = nb,
  norder   = norder,
  breaks   = knots
)

bb0 <- getbasismatrix(tau, basisobj)  # lt x nb

## ----------------------------
## 5) Precompute B(t, tau_i) slices
## ----------------------------
# Each element is a (length(t_grid) x pB) matrix for fixed tau_i
B_by_tau <- vector("list", lt)
for (i in seq_len(lt)) {
  B_by_tau[[i]] <- B.est[z[, 2] == tau[i], , drop = FALSE]
}

## ----------------------------
## 6) Roughness penalty matrix
## ----------------------------
rp <- t(q2) %*% Kmat %*% q2  # p_q2 x p_q2

## ----------------------------
## 7) CV tuning loop
## ----------------------------
ncv  <- 10
cand <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)

folds <- make_folds(n, K = ncv, seed = 1)

# helper to build design at a fixed tau index i for a given row-index set (idx)
build_design <- function(idx, i) {
  bb_tau <- bb0[i, , drop = FALSE]           # 1 x nb
  Btau   <- B_by_tau[[i]]                    # Mobs x pB
  p      <- (hDM / 3) * mtx[idx, , drop = FALSE] %*% Wsimpson %*% Btau  # |idx| x pB
  X      <- cbind(p %*% q2, matrix(bb_tau, nrow = nrow(p), ncol = ncol(bb0), byrow = TRUE))
  list(X = X, bb_tau = bb_tau)
}

gammacv <- matrix(NA_real_, nrow = length(cand), ncol = 2)
colnames(gammacv) <- c("cand_index", "cv_loss")

for (l in seq_along(cand)) {
  cvloss <- 0
  
  for (cv in seq_len(ncv)) {
    cat("cand=", l, " fold=", cv, "\n")
    
    tind  <- folds[[cv]]
    trind <- setdiff(seq_len(n), tind)
    
    # Variables
    beta  <- Variable(ncol(q2) + ncol(bb0))
    theta <- Variable(rows = ncol(bb0), cols = ncol(add_x))  # nb x p_add
    
    roughness <- quad_form(beta[1:ncol(q2)], rp)
    
    # objective across all taus
    obj <- 0
    for (i in seq_len(lt)) {
      des  <- build_design(trind, i)
      X    <- des$X
      bb_i <- des$bb_tau
      # add_x %*% t(bb_i %*% theta) gives n x 1 because bb_i is 1 x nb
      mu_add <- add_x[trind, , drop = FALSE] %*% t(bb_i %*% theta)
      obj <- obj + sum(quant_loss(Y[trind] - X %*% beta - mu_add, tau[i]))
    }
    
    prob <- Problem(Minimize(obj + cand[l] * n * roughness))
    
    fit <- psolve(prob, solver = "ECOS")
    
    g  <- as.matrix(fit$getValue(beta))   # (p_q2 + nb) x 1
    gg <- as.matrix(fit$getValue(theta))  # nb x p_add
    
    # compute CV loss on test fold
    fold_loss <- 0
    for (i in seq_len(lt)) {
      des  <- build_design(tind, i)
      X    <- des$X
      bb_i <- des$bb_tau
      mu_add <- add_x[tind, , drop = FALSE] %*% t(bb_i %*% gg)
      res <- Y[tind] - X %*% g - mu_add
      fold_loss <- fold_loss + sum(quant_loss(res, tau[i]))
    }
    
    cvloss <- cvloss + fold_loss
  }
  
  gammacv[l, ] <- c(l, cvloss)
  saveRDS(gammacv, file = "gammacv.RDS")
}

print(gammacv)
best_l <- gammacv[which.min(gammacv[, 2]), 1]
cat("Best cand index =", best_l, "  cand =", cand[best_l], "\n")
