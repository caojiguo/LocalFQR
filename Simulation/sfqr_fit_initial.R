## =========================
## 0) Packages
## =========================
ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  invisible(TRUE)
}

ensure_github_pkg <- function(repo, pkg = basename(repo)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    ensure_pkg("devtools")
    devtools::install_github(repo)
  }
  invisible(TRUE)
}

load_pkgs <- function() {
  ensure_github_pkg("FIRST-Data-Lab/BPST", "BPST")
  ensure_pkg("Matrix")
  ensure_pkg("fda")
  ensure_pkg("tensr")
  ensure_pkg("CVXR")
  ensure_pkg("MASS") 
  ensure_pkg("e1071")# you use MASS but didn't ensure it above
  # fdaPDE is not used in this snippet; load only if you need it elsewhere:
  # ensure_pkg("fdaPDE")
  
  suppressPackageStartupMessages({
    library(BPST)
    library(Matrix)
    library(fda)
    library(tensr)
    library(CVXR)
    library(MASS)
  })
  invisible(TRUE)
}

load_pkgs()

## =========================
## 1) Small utilities
## =========================
quant_loss <- function(u, tau) 0.5 * abs(u) + (tau - 0.5) * u

# Split indices into K folds (deterministic if set.seed called outside)
make_folds <- function(n, K) {
  idx <- sample.int(n)
  split(idx, rep(1:K, length.out = n))
}

# Simpson weights (Mobs must be odd; checks help avoid silent bugs)
simpson_weights <- function(Mobs) {
  if (Mobs %% 2 == 0) stop("Simpson rule requires an odd number of grid points (Mobs odd).")
  c(1, rep(c(4, 2), (Mobs - 3) / 2), 4, 1)
}

## =========================
## 2) Precompute design blocks that depend only on (tau, triangulation, X grid)
## =========================
precompute_sfqr_blocks <- function(tau, V.est, Tr.est, d, r,
                                   newgrids, t_grid, mtx) {
  lt <- length(tau)
  
  # Build z grid and Bernstein basis over triangulation
  z_full <- cbind(rep(newgrids, each = lt), rep(tau, times = length(newgrids)))
  Bfull <- suppressWarnings(basis(V.est, Tr.est, d, r, z_full))
  B.est <- Bfull$B
  q2    <- Bfull$Q2
  Kmat  <- as.matrix(Bfull$K)
  ind   <- Bfull$Ind.inside
  z     <- z_full[ind, , drop = FALSE]
  
  # Roughness penalty in reduced coordinates
  nbern <- choose(d + 2, 2)
  rp <- t(q2) %*% Kmat %*% q2
  
  # B-spline basis in u
  nk <- 10
  knots  <- quantile(tau, seq(0, 1, length.out = nk))
  norder <- 2
  nb     <- length(knots) + norder - 2
  basisobj <- create.bspline.basis(rangeval = c(min(tau), max(tau)),
                                   nbasis = nb, norder = norder, breaks = knots)
  bb0 <- getbasismatrix(tau, basisobj)
  
  # Simpson setup
  Mobs <- ncol(mtx)
  hDM  <- (max(t_grid) - min(t_grid)) / Mobs
  wSim <- simpson_weights(Mobs)
  
  # For each quantile level, store the Bernstein basis slice bb_u (rows correspond to grid t)
  # Note: your code uses B.est[z[,2]==tau[r],]; keep the same logic.
  B_by_u <- lapply(seq_len(lt), function(rq) B.est[z[,2] == tau[rq], , drop = FALSE])
  
  list(
    tau = tau, lt = lt,
    z = z,
    B_by_u = B_by_u,
    q2 = q2,
    bb0 = bb0,
    rp = rp,
    hDM = hDM,
    wSim = wSim
  )
}

## =========================
## 3) Build design matrices for a set of indices (train or test)
## =========================
build_design_list <- function(idx, mtx, add_x, blocks) {
  lt  <- blocks$lt
  q2  <- blocks$q2
  bb0 <- blocks$bb0
  hDM <- blocks$hDM
  wSim <- blocks$wSim
  B_by_u <- blocks$B_by_u
  
  # Each element: list(X = ..., bb1 = ...)
  designs <- vector("list", lt)
  
  for (rq in seq_len(lt)) {
    bb <- B_by_u[[rq]]
    # Simpson integral approx: h/3 * X * diag(w) * bb
    p <- (hDM / 3) * mtx[idx, , drop = FALSE] %*% diag(wSim) %*% bb
    X <- cbind(p %*% q2, matrix(bb0[rq, ], nrow = length(idx), ncol = ncol(bb0), byrow = TRUE))
    bb1 <- matrix(bb0[rq, ], nrow = 1, ncol = ncol(bb0), byrow = TRUE)
    designs[[rq]] <- list(X = X, bb1 = bb1)
  }
  
  designs
}

## =========================
## 4) Fit model for ONE training split and ONE gamma (cand2 value)
## =========================
fit_sfqr_cvxr <- function(y_tr, designs_tr, tau, rp, gamma,
                          add_x_tr = NULL) {
  lt <- length(tau)
  p_beta <- ncol(designs_tr[[1]]$X)
  beta <- CVXR::Variable(p_beta)
  
  if (is.null(add_x_tr)) {
    obj <- 0
    for (rq in seq_len(lt)) {
      X <- designs_tr[[rq]]$X
      obj <- obj + sum(quant_loss(y_tr - X %*% beta, tau[rq]))
    }
    # roughness only on the q2 part (first ncol(rp) entries)
    roughness <- quad_form(beta[1:ncol(rp)], rp) * length(y_tr)
    prob <- Problem(Minimize(obj + gamma * roughness))
    fit  <- psolve(prob)
    list(beta = fit$getValue(beta), theta = NULL)
  } else {
    # If you truly need theta: keep it, but make it explicit
    # theta dims: (ncol(bb1) x ncol(add_x)) per your original code
    # NOTE: Your original model uses add_x %*% t(bb1 %*% theta)
    bb1 <- designs_tr[[1]]$bb1
    theta <- CVXR::Variable(rows = ncol(bb1), cols = ncol(add_x_tr))
    
    obj <- 0
    for (rq in seq_len(lt)) {
      X <- designs_tr[[rq]]$X
      bb1 <- designs_tr[[rq]]$bb1
      obj <- obj + sum(quant_loss(y_tr - X %*% beta - add_x_tr %*% t(bb1 %*% theta), tau[rq]))
    }
    roughness <- quad_form(beta[1:ncol(rp)], rp) * length(y_tr)
    prob <- Problem(Minimize(obj + gamma * roughness))
    fit  <- psolve(prob)
    list(beta = fit$getValue(beta), theta = fit$getValue(theta))
  }
}

## =========================
## 5) Compute CV loss on a held-out split
## =========================
cv_loss_split <- function(y_te, designs_te, tau, beta_hat,
                          add_x_te = NULL, theta_hat = NULL) {
  lt <- length(tau)
  loss <- 0
  for (rq in seq_len(lt)) {
    X <- designs_te[[rq]]$X
    bb1 <- designs_te[[rq]]$bb1
    if (is.null(add_x_te)) {
      resid <- y_te - X %*% beta_hat
    } else {
      resid <- y_te - X %*% beta_hat - add_x_te %*% t(bb1 %*% theta_hat)
    }
    loss <- loss + sum(quant_loss(resid, tau[rq]))
  }
  loss
}

## =========================
## 6) Cross-validation over gamma candidates
## =========================
cv_sfqr_gamma <- function(y, mtx, add_x, blocks, cand2, ncv = 10, seed = 1) {
  set.seed(seed)
  folds <- make_folds(length(y), ncv)
  
  out <- numeric(length(cand2))
  
  for (lg in seq_along(cand2)) {
    gamma <- cand2[lg]
    cvloss <- 0
    
    for (cv in seq_len(ncv)) {
      tind  <- folds[[cv]]
      trind <- setdiff(seq_along(y), tind)
      
      designs_tr <- build_design_list(trind, mtx, add_x, blocks)
      designs_te <- build_design_list(tind,  mtx, add_x, blocks)
      
      fit <- fit_sfqr_cvxr(
        y_tr = y[trind],
        designs_tr = designs_tr,
        tau = blocks$tau,
        rp = blocks$rp,
        gamma = gamma,
        add_x_tr = if (is.null(add_x)) NULL else add_x[trind, , drop = FALSE]
      )
      
      cvloss <- cvloss + cv_loss_split(
        y_te = y[tind],
        designs_te = designs_te,
        tau = blocks$tau,
        beta_hat = fit$beta,
        add_x_te = if (is.null(add_x)) NULL else add_x[tind, , drop = FALSE],
        theta_hat = fit$theta
      )
    }
    
    out[lg] <- cvloss
    message("gamma idx = ", lg, ", cvloss = ", cvloss)
  }
  
  out
}

## =========================
## 7) Usage (matches your current workflow)
## =========================
tmpid <- 1
tau <- seq(0.05, 0.95, length.out = 19)
cand2 <- c(1e-8, 1e-9, 1e-10, 1e-11, 1e-12)

set.seed(tmpid)
nobs <- 200
source("DataGen_I.R")
V.est  <- readRDS("nodes.RDS")
Tr.est <- readRDS("tri.RDS")

blocks <- precompute_sfqr_blocks(
  tau = tau,
  V.est = V.est, Tr.est = Tr.est,
  d = 2, r = 1,
  newgrids = newgrids,
  t_grid = t_grid,
  mtx = mtx
)

gammacv <- cv_sfqr_gamma(
  y = daty,
  mtx = mtx,
  add_x = add_x,      # can be NULL
  blocks = blocks,
  cand2 = cand2,
  ncv = 5,
  seed = tmpid
)

saveRDS(gammacv, file = paste0(tmpid, "_gammacv.RDS"))

