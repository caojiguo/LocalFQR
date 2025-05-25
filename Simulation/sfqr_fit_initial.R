if (!requireNamespace("BPST", quietly = TRUE)) {
  devtools::install_github("FIRST-Data-Lab/BPST")
}
if (!requireNamespace("stats", quietly = TRUE)) {
  install.packages("stats")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix")
}
if (!requireNamespace("fda", quietly = TRUE)) {
  install.packages("fda")
}
if (!requireNamespace("tensr", quietly = TRUE)) {
  install.packages("tensr")
}
if (!requireNamespace("CVXR", quietly = TRUE)) {
  install.packages("CVXR")
}
if (!requireNamespace("fdaPDE", quietly = TRUE)) {
  install.packages("fdaPDE")
}
library(stats)
library(BPST)
library(Matrix)
library(fda)
library(MASS)
library(tensr)
library(CVXR)
# index of simulations
tmpid <- 1
# quantiles of interest
tau <- seq(0.05,0.95,length.out = 19)
lt <- length(tau)
# sample size
nobs <- 500
# number of knots for B-splines used to approximate the coefficients for non-functional covariate
nk <- 10
# number of folds for cross validation
ncv <- 10
# candidate values for roughness penalty
cand2 <- c(1e-8,1e-9,1e-10, 1e-11, 1e-12)
# degree of Bernstein polynomials
d <- 2
# smoothness conditions for Bernstein polynomial approximation over the triangulation
r <- 1
set.seed(tmpid)
source('DataGen_I.R')
V.est <- readRDS(file = 'nodes.RDS')
Tr.est <- readRDS(file='tri.RDS')
# number of triangles
ntr <- nrow(Tr.est)
# number of basis functions on each triangle
nbern <- choose(d+2,2)
# area of each triangle
area <- c()
for (i in 1:ntr)
{
  # index of vertex
  v.index <- Tr.est[i,]
  a <- V.est[v.index[1],]
  b <- V.est[v.index[2],]
  c <- V.est[v.index[3],]
  tmp <- abs(a[1]*(b[2]-c[2]) + b[1]*(c[2]-a[2]) + c[1]*(a[2]-b[2]))/2
  area <- c(area,tmp)
}
area <- sqrt(area)


z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))
# compute the values of Bernstein polynomials based on the given triangulation
Bfull.est <- basis(V.est,Tr.est,d,r,z)
B.est <- Bfull.est$B
# matrix derived from QR decomposition
q2 <- Bfull.est$Q2
saveRDS(object=B.est,file='B.est.RDS')
saveRDS(object=q2,file='q2.RDS')
# symmetric and positive semi-definite matrix corresponding to the roughness penalty
K <- as.matrix(Bfull.est$K)


# Define the B-splines
knots <- quantile(tau,seq(0,1,length=nk))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
basisobj=create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)

ind <- Bfull.est$Ind.inside
z <- z[ind,]

tt <- z[z[,2]==tau[1],1]
l <- length(tt)
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs
cefDM   = c(1, rep(c(4,2), (Mobs-3)/2), 4,  1)
n <- nrow(mtx)
y <- daty

ngrp <- nrow(q2)/nbern
quant_loss <- function(u, t) {0.5 * abs(u) + (t - 0.5) * u}
rp <-t(q2)%*%K%*%q2

# use cross validation to tune roughness parameter
gammacv <- c()
for (lg in 1:length(cand2)){
  cvloss <- 0
  # cross validation starts
  for (cv in 1:ncv){
    lc <- length(daty)/ncv
    tind <- (cv-1)*lc + c(1:lc)
    trind <- c(1:length(daty))[-tind]
    bb <- B.est[z[,2]==tau[1],]
    bb1 <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined over the triangulation
    p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
    X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(mtx[trind,]),ncol=ncol(bb0),byrow=T))
    beta <- Variable(ncol(q2) + ncol(bb0))
    if (is.null(add_x)){
      obj <- sum(quant_loss(daty[trind] - X%*% beta,tau[1]))
    }else{
      theta <- Variable(rows = ncol(bb1),cols = ncol(add_x))
      obj <- sum(quant_loss(daty[trind] - X%*% beta-add_x[trind,]%*%t(bb1%*%theta),tau[1]))
    }
    roughness <- quad_form(beta[1:ncol(q2)],rp)*length(trind)
    for (i in 1:(lt-1))
    {
      tt <- z[z[,2]==tau[i+1],1]
      l <- length(tt)
      bb <- B.est[z[,2]==tau[i+1],]
      bb1 <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
      # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined over the triangulation
      p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
      X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(mtx[trind,]),ncol=ncol(bb0),byrow=T))
      if (is.null(add_x)){
        obj <- obj + sum(quant_loss(daty[trind] - X %*% beta,tau[i+1]))
      }else{
        obj <- obj + sum(quant_loss(daty[trind] - X %*% beta-add_x[trind,]%*%t(bb1%*%theta),tau[i+1]))
      }
    }
    print(lg)
    prob <- Problem(Minimize(obj+cand2[lg]*roughness))
    fit <- psolve(prob)
    g <- fit$getValue(beta)
    if (is.null(add_x)==0){
      gg <- fit$getValue(theta)
    }

    bb <- B.est[z[,2]==tau[1],]
    # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined over the triangulation
    p <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
    X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(mtx[tind,]),ncol=ncol(bb0),byrow=T))
    bb1 <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    if (is.null(add_x)){
      cvres <- daty[tind]- X%*%g
    }else{
      cvres <- daty[tind]- X%*%g -add_x[tind,]%*%t(bb1%*%gg)
    }
    # compute the quantile loss on the testing set
    cvloss <- cvloss + sum(quant_loss(cvres,tau[1]))
    for (i in 1:(lt-1))
    {
      tt <- z[z[,2]==tau[i+1],1]
      l <- length(tt)
      bb <- B.est[z[,2]==tau[i+1],]
      # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined over the triangulation
      p <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
      X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(mtx[tind,]),ncol=ncol(bb0),byrow=T))
      bb1 <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
      if (is.null(add_x)){
        cvres <- daty[tind]- X%*%g
      }else{
        cvres <- daty[tind]- X%*%g - add_x[tind,]%*%t(bb1%*%gg)
      }
      cvloss <- cvloss + sum(quant_loss(cvres,tau[i+1]))
    }
  }
  gammacv <- rbind(gammacv, cvloss)
  tmpname <- paste0(tmpid,'gammacv.RDS')
  saveRDS(object = gammacv, file=tmpname)
}
