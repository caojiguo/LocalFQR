## ============================================================
## Soybean data analysis (minimal revision; keep original logic)
## - Adds comments
## - Wraps repeated code in small helpers
## - Fixes: B.est/z alignment
## - Keeps: your quant_loss and CVXR::norm(type="2") usage
## ============================================================

suppressPackageStartupMessages({
  library(stats)
  library(BPST)
  library(Matrix)
  library(MASS)
  library(tensr)
  library(CVXR)
  library(fda)
  library(plotly)
})

## ----------------------------
## 0) Basic setup + data
## ----------------------------
tau <- seq(0.05, 0.95, length.out = 19)
lt  <- length(tau)

tmpx    <- readRDS("soy_avg_x.RDS")
daty    <- readRDS("soy_avg_y.RDS")
add_var <- readRDS("soy_avg_add.RDS")

add_var[, 3] <- add_var[, 3] / add_var[, 2]
add_x <- add_var[, -2]
colnames(add_x) <- c("prcp", "irr_ratio")

mtx2 <- tmpx[,  30:334]  # min temp
mtx1 <- tmpx[, 396:700]  # max temp
mtx  <- (mtx1 + mtx2) / 2

n    <- nrow(mtx)
Mobs <- ncol(mtx)

t_grid <- newgrids <- seq(0, 1, length.out = Mobs)
Y <- daty

## ----------------------------
## 1) Triangulation basis
## ----------------------------
V.est  <- readRDS("nodes.RDS")
Tr.est <- readRDS("tri.RDS")

ntr   <- nrow(Tr.est)
d     <- 2
r     <- 1
nbern <- choose(d + 2, 2)

## triangle areas (kept for completeness; not used later)
area <- numeric(ntr)
for (i in 1:ntr) {
  v.index <- Tr.est[i, ]
  a <- V.est[v.index[1], ]
  b <- V.est[v.index[2], ]
  c <- V.est[v.index[3], ]
  area[i] <- abs(a[1]*(b[2]-c[2]) + b[1]*(c[2]-a[2]) + c[1]*(a[2]-b[2]))/2
}
area <- sqrt(area)

## (t,u) grid for basis evaluation
z <- cbind(rep(newgrids, each = lt), rep(tau, times = length(newgrids)))

Bfull.est <- basis(V.est, Tr.est, d, r, z)
B.est <- Bfull.est$B
q2    <- Bfull.est$Q2
K     <- as.matrix(Bfull.est$K)

ind <- Bfull.est$Ind.inside

## IMPORTANT FIX: subset BOTH z and B.est
z     <- z[ind, , drop = FALSE]
B.est <- B.est[ind, , drop = FALSE]

ngrp <- nrow(q2) / nbern

## ----------------------------
## 2) Simpson rule weights (keep your original definition)
## NOTE: you used hDM=(tmax-tmin)/Mobs originally; keep it to match old output.
## ----------------------------
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM  <- (tmax - tmin) / Mobs

cefDM <- c(1, rep(c(4,2), (Mobs-3)/2), 4, 1)

## ----------------------------
## 3) B-spline basis in tau
## ----------------------------
knots  <- quantile(tau, seq(0, 1, length = 10))
norder <- 2
nb     <- length(knots) + norder - 2

basisobj <- create.bspline.basis(
  rangeval = c(min(tau), max(tau)),
  nbasis   = nb,
  norder   = norder,
  breaks   = knots
)
bb0 <- getbasismatrix(tau, basisobj)

## ----------------------------
## 4) Loss + roughness matrix
## ----------------------------
quant_loss <- function(u, t) {
  0.5 * abs(u) + (t - 0.5) * u
}
rp <- t(q2) %*% K %*% q2

## ------------------------------------------------------------
## Helper: build X, bbtheta at quantile index i (1..lt)
## ------------------------------------------------------------
build_X <- function(i, idx = NULL) {
  # idx: row indices of mtx/add_x/daty; if NULL uses all rows
  if (is.null(idx)) idx <- seq_len(n)
  
  bb <- B.est[z[,2] == tau[i], , drop = FALSE]  # Mobs x pB
  p  <- (hDM/3) * mtx[idx, , drop = FALSE] %*% diag(cefDM) %*% bb
  X  <- cbind(p %*% q2,
              matrix(bb0[i, ], nrow = nrow(p), ncol = ncol(bb0), byrow = TRUE))
  bbtheta <- matrix(bb0[i, ], nrow = 1, ncol = ncol(bb0), byrow = TRUE)
  
  list(X = X, bbtheta = bbtheta)
}

## ------------------------------------------------------------
## Helper: objective across all taus for given beta/theta
## ------------------------------------------------------------
build_obj_all_tau <- function(beta, theta, idx = NULL, yvec = daty) {
  if (is.null(idx)) idx <- seq_len(n)
  
  out <- 0
  for (i in 1:lt) {
    tmp <- build_X(i, idx)
    X  <- tmp$X
    bbtheta <- tmp$bbtheta
    
    out <- out + sum(
      quant_loss(
        yvec[idx] - X %*% beta - add_x[idx, , drop = FALSE] %*% t(bbtheta %*% theta),
        tau[i]
      )
    )
  }
  out
}

## ============================================================
## Part A: Roughness-only fit (to compute lambda weights)
## ============================================================
cand <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
gammacv <- readRDS("gammacv.RDS")
optid <- which.min(gammacv[,2])
print(gammacv[optid, ])

beta  <- Variable(ncol(q2) + ncol(bb0))
theta <- Variable(rows = ncol(bb0), cols = ncol(add_x))
roughness <- quad_form(beta[1:ncol(q2)], rp)

obj1 <- build_obj_all_tau(beta, theta, idx = NULL, yvec = daty)

gamma0 <- cand[gammacv[optid, 1]] * n
prob1  <- Problem(Minimize(obj1 + gamma0 * roughness))

## keep your original solver behavior: solve(prob)
fit1 <- solve(prob1, solver='SCS')
g1   <- fit1$getValue(beta)

## lambda weights (keep your formula)
int_mtx <- diag(1, nbern)
lambda <- numeric(ngrp)

for (k in 1:ngrp) {
  rows_k <- (1:nbern) + (k - 1) * nbern
  lambda[k] <- as.numeric(
    (g1[1:ncol(q2), ] %*% t(q2[rows_k, , drop = FALSE]) %*% int_mtx %*%
       q2[rows_k, , drop = FALSE] %*% g1[1:ncol(q2), ])^0.5
  )
}

## ============================================================
## Part B: Group-lasso + roughness fit
## ============================================================
alpha <- -1

cand1 <- c(1e-4, 1e-3, 1e-2, 1e-1, 1)
cand2 <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
cand3 <- cand2

l1 <- 4
l2 <- 5
l3 <- 2

scale <- cand1[l1] * n
gamma <- cand2[l2] * n

beta  <- Variable(ncol(q2) + ncol(bb0))
theta <- Variable(rows = ncol(bb0), cols = ncol(add_x))
roughness <- quad_form(beta[1:ncol(q2)], rp)

## IMPORTANT: keep EXACTLY the norm syntax that worked for you:
## use norm(..., type="2") and optionally force CVXR::norm
grp_lasso <- (lambda[1]^alpha) * CVXR::norm(int_mtx %*% q2[1:nbern, , drop = FALSE] %*% beta[1:ncol(q2)],
                                            type = "2")
if (ngrp >= 2) {
  for (k in 2:ngrp) {
    rows_k <- (1:nbern) + (k - 1) * nbern
    grp_lasso <- grp_lasso +
      (lambda[k]^alpha) * CVXR::norm(int_mtx %*% q2[rows_k, , drop = FALSE] %*% beta[1:ncol(q2)],
                                     type = "2")
  }
}

obj2 <- build_obj_all_tau(beta, theta, idx = NULL, yvec = daty)

prob2 <- Problem(Minimize(obj2 + scale * grp_lasso + gamma * roughness))

## keep your choice (SCS) to match old behavior
fit2 <- psolve(prob2, solver = "SCS")

g  <- fit2$getValue(beta)
gg <- fit2$getValue(theta)

## ============================================================
## Part C: Identify active regions + refit
## ============================================================
b1 <- q2 %*% g[1:ncol(q2)]
norm1 <- numeric(ngrp)

for (k in 1:ngrp) {
  rows_k <- (1:nbern) + (k - 1) * nbern
  norm1[k] <- sqrt(sum(b1[rows_k, ]^2))
}

print(sum(norm1 <= 1e-3))
inactive1 <- rep(norm1, each = nbern)

if (sum(norm1 >= 1e-3) > 0) {
  if (sum(inactive1 <= 1e-3) > 0) {
    qs  <- q2[inactive1 <= 1e-3, , drop = FALSE]
    q21 <- qrH(qs)
    q21 <- as.matrix(q21)
  } else {
    q21 <- diag(1, nrow = ncol(q2))
  }
  
  rp2 <- t(q21) %*% t(q2) %*% K %*% q2 %*% q21
  z1  <- sum(norm1 >= 1e-3)
} else {
  z1 <- NULL
}

if (!is.null(z1)) {
  beta <- Variable(ncol(q21) + ncol(bb0))
  roughness3 <- quad_form(beta[1:ncol(q21)], rp2)
  
  theta <- Variable(rows = ncol(bb0), cols = ncol(add_x))
  
  ## objective with screened design
  obj3 <- 0
  for (i in 1:lt) {
    bb <- B.est[z[,2] == tau[i], , drop = FALSE]
    p  <- (hDM/3) * mtx %*% diag(cefDM) %*% bb
    X  <- cbind(p %*% q2 %*% q21,
                matrix(bb0[i, ], nrow = nrow(p), ncol = ncol(bb0), byrow = TRUE))
    bbtheta <- matrix(bb0[i, ], nrow = 1, ncol = ncol(bb0), byrow = TRUE)
    
    obj3 <- obj3 + sum(
      quant_loss(daty - X %*% beta - add_x %*% t(bbtheta %*% theta), tau[i])
    )
  }
  
  prob3 <- Problem(Minimize(obj3 + cand3[l3] * roughness3))
  fit3  <- psolve(prob3)   # keep your default choice
  
  g  <- fit3$getValue(beta)
  gg <- fit3$getValue(theta)
}

## ============================================================
## Part D: Plots (kept basically same)
## ============================================================
betas <- NULL
for (u in tau) {
  betas <- rbind(betas, as.numeric(B.est[z[,2] == u, ] %*% q2 %*% q21 %*% g[1:ncol(q21), ]))
}

fig <- plot_ly(x = 30:335, y = tau, z = betas, colors = "YlOrRd") %>% add_surface()

axx <- list(
  title = "t (Month)",
  ticktext = list("Feb","Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec"),
  tickvals = c(30,61,91,122,152,183,213,243,273,304,335)
)
axy <- list(
  title = "u (Quantile)",
  tickvals = seq(0.05, 0.95, by = 0.15),
  ticktext = as.character(round(seq(0.05, 0.95, by = 0.15), 2))
)
axz <- list(title = "beta(t,u)")

fig <- fig %>% layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))
fig

par(mar=c(5.1, 4.5, 4.1, 2.1))
plot(tau,
     B.est[z[,1] == (unique(z[,1])[291]), ] %*% q2 %*% q21 %*% g[1:ncol(q21), ],
     type = "l", xlab = "Quantile", ylab = "beta(t,u)",
     cex.lab = 1.5, cex.axis = 1.5, xaxt = "n", lwd = 2, col = "red")
axis(1, at = c(tau[1], tau[5], tau[9], tau[13], tau[17]),
     labels = c("5%", "25%", "50%", "75%", "95%"),
     cex.axis = 1.5)
