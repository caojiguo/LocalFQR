library(stats)
library(BPST)
library(Matrix)
library(fda)
library(MASS)
library(tensr)
library(CVXR)
library(Rmosek)
# index of simulations
tmpid <- 1
# quantiles of interest
tau <- seq(0.25,0.75,length.out = 17)
lt <- length(tau)
# sample size
nobs <- 300
# number of knots for B-splines
nk <- 10
# number of folds for cross validation
ncv <- 10
# candidate values for roughness penalty
cand2 <- c(5e-4,1e-4,5e-5,1e-5,5e-6)
# degree of Bernstein polynomials
d <- 2
# smoothness conditions for Bernstein polynomial approximation over the triangulation
r <- 1

# load the data generating model
source('simple_DataGen_addx.R')
# input the triangulation
V.est <- readRDS(file = 'nodes.RDS')
Tr.est <- readRDS(file='tri.RDS')
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

n <- nrow(mtx)
y <- daty
ind <- Bfull.est$Ind.inside
z <- z[ind,]
tt <- z[z[,2]==tau[1],1]
l <- length(tt)
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs
nbern <- choose(d+2,2)
ngrp <- nrow(q2)/nbern

# define quantile loss 
quant_loss <- function(u, t) {0.5 * abs(u) + (t - 0.5) * u}
rp <-t(q2)%*%K%*%q2

# use cross validation to decide the tuning parameter corresponding to the roughness penalty
gammacv <- c()
for (lg in 1:length(cand2)){
  # lg is the index of tuning parameter within the candidate set
  cvloss <- 0
  for (cv in 1:ncv){
    # cv is the index of validation set within the cross validation
    lc <- length(daty)/ncv
    tind <- (cv-1)*lc + c(1:lc)
    trind <- c(1:length(daty))[-tind]
    bb <- B.est[z[,2]==tau[1],]
    bb1 <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    p <- hDM*mtx[trind,]%*%bb
    X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(mtx[trind,]),ncol=ncol(bb0),byrow=T))
    # define the unknown approximation coefficients
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
      p <- hDM*mtx[trind,]%*%bb
      X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(mtx[trind,]),ncol=ncol(bb0),byrow=T))
      if (is.null(add_x)){
        obj <- obj + sum(quant_loss(daty[trind] - X %*% beta,tau[i+1]))
      }else{
        obj <- obj + sum(quant_loss(daty[trind] - X %*% beta-add_x[trind,]%*%t(bb1%*%theta),tau[i+1]))
      }
    }
    print(lg)
    prob <- Problem(Minimize(obj+cand2[lg]*roughness))
    fit <- psolve(prob,solver="MOSEK")
    g <- fit$getValue(beta)
    if (is.null(add_x)==0){
      gg <- fit$getValue(theta)
    }

    bb <- B.est[z[,2]==tau[1],]
    p <- hDM*mtx[tind,]%*%bb
    X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(mtx[tind,]),ncol=ncol(bb0),byrow=T))
    bb1 <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    if (is.null(add_x)){
      cvres <- daty[tind]- X%*%g
    }else{
      cvres <- daty[tind]- X%*%g -add_x[tind,]%*%t(bb1%*%gg)
    }
    cvloss <- cvloss + sum(quant_loss(cvres,tau[1]))
    for (i in 1:(lt-1))
    {
      tt <- z[z[,2]==tau[i+1],1]
      l <- length(tt)
      bb <- B.est[z[,2]==tau[i+1],]
      p <- hDM*mtx[tind,]%*%bb
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

