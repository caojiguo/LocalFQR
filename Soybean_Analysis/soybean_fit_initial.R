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
if (!requireNamespace("plotly", quietly = TRUE)) {
  install.packages("plotly")
}
if (!requireNamespace("fdaPDE", quietly = TRUE)) {
  install.packages("fdaPDE")
}
library(stats)
library(BPST)
library(Matrix)
library(MASS)
library(tensr)
library(CVXR)
library(fda)
tau <- seq(0.05,0.95,length.out = 19)

tmpx <- readRDS('soy_avg_x.RDS')
daty <- readRDS('soy_avg_y.RDS')
add_var <- readRDS('soy_avg_add.RDS')
add_var[,3] <- add_var[,3]/add_var[,2]
add_x <- add_var[,-2]
names(add_x) <- c('prcp','irr_ratio')
mtx2 <- tmpx[,30:334] # minimum temperature
mtx1 <- tmpx[,396:700] # maximum temperature
mtx <- (mtx1 + mtx2)/2

t_grid <- newgrids <- seq(0,1,length.out = ncol(mtx1))
Y <- daty
V.est <- readRDS(file = 'nodes.RDS')
Tr.est <- readRDS(file='tri.RDS')

# number of triangles
ntr <- nrow(Tr.est)
# number of basis functions on each triangle
d <- 2
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

lt <- length(tau)
z <- cbind(rep(newgrids,each=length(tau)),rep(tau,length(newgrids)))
r <- 1
Bfull.est <- basis(V.est,Tr.est,d,r,z)
B.est <- Bfull.est$B
q2 <- Bfull.est$Q2
K <- as.matrix(Bfull.est$K)
dim(B.est)
ind <- Bfull.est$Ind.inside
z <- z[ind,]
Mobs <- ncol(mtx1)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM  <- (tmax-tmin)/Mobs
cefDM <- c(1, rep(c(4,2), (Mobs-3)/2), 4,  1)

knots <- quantile(tau,seq(0,1,length=31))
norder <- 2
nknots <- length(knots)
nb <- nknots + norder - 2
basisobj <- create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)

n <- nrow(mtx)
y <- daty

ngrp <- nrow(q2)/nbern
quant_loss <- function(u, t) {0.5 * abs(u) + (t - 0.5) * u}
rp <-t(q2)%*%K%*%q2

ncv <- 10
# use cross validation to tune roughness parameter
cand <- c(1e-6,1e-5,1e-4,1e-3,1e-2)
gammacv <- c()
for (l in 1:length(cand)){
  cvloss <- 0
  for (cv in 1:ncv){
    print(c(l,cv))
    lc <- length(daty)/ncv
    lc <- round(lc)
    tind <- c(1:length(daty))[c(1:length(daty))%%ncv==1]
    tind <- tind[-length(tind)]
    tind <- tind + cv-1
    trind <- c(1:length(daty))[-tind]
    bb <- B.est[z[,2]==tau[1],]
    p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
    X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
    bb1 <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    beta <- Variable(ncol(q2) + ncol(bb0))
    theta <- Variable(rows = ncol(bb1),cols = ncol(add_x))
    roughness <- quad_form(beta[1:ncol(q2)],rp)
    obj <- sum(quant_loss(daty[trind] - X%*% beta -  add_x[trind,]%*%t(bb1%*%theta),tau[1]))
    for (i in 1:(lt-1))
    {
      # Simpson's rule
      bb <- B.est[z[,2]==tau[i+1],]
      p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
      X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
      bb1 <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
      obj <- obj + sum(quant_loss(daty[trind] - X %*% beta - add_x[trind,]%*%t(bb1%*%theta),tau[i+1]))
    }
    prob <- Problem(Minimize(obj+cand[l]*n*roughness))
    fit <- psolve(prob,solver='ECOS')
    g <- fit$getValue(beta)
    g <- as.matrix(g)
    print(g[1:5])
    gg <- fit$getValue(theta)
    print(gg[1:5])
    gg <- as.matrix(gg)
    bb <- B.est[z[,2]==tau[1],]
    # calculate loss in testing dataset
    p <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
    X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
    bb1 <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    cvres <- daty[tind]- X%*%g - add_x[tind,]%*%t(bb1%*%gg)
    cvloss <- cvloss + sum(quant_loss(cvres,tau[1]))
    for (i in 1:(lt-1))
    {
      # Simpson's rule
      bb <- B.est[z[,2]==tau[i+1],]
      p <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
      X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
      bb1 <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
      cvres <-  daty[tind]- X%*%g - add_x[tind,]%*%t(bb1%*%gg)
      cvloss <- cvloss + sum(quant_loss(cvres,tau[i+1]))
    }
    # print(cvloss)
  }
  gammacv <- rbind(gammacv, c(l,cvloss))
  saveRDS(object = gammacv, file='gammacv.RDS')
} 
