library(BPST)
library(Matrix)
library(MASS)
library(tensr)
library(CVXR)
library(fda)
library(Rmosek)
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
# candidate values for group lasso type penalty
cand1 <- c(1e-2,1e-3,1e-4)
# candidate values for roughness penalty
cand3 <- cand2 <- c(1e-10, 1e-11, 1e-12)
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
# find the linear constraints corresponding to the desired smoothness of the estimator
H <- smoothness(V.est, Tr.est, d, r)
H <- as.matrix(H)
H1 <- smoothness(V.est, Tr.est, d, r-1)
H1 <- as.matrix(H1)
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

# load the cross validation results for the initial estimate of bivariate slope function and decide the optimal tuning parameter value for the roughness penalty
tmpname <- paste0(tmpid,'gammacv.RDS')
gammacv <- readRDS(tmpname)
optid <- which(gammacv==min(gammacv))[1]
bb <- B.est[z[,2]==tau[1],]

# Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
p <- hDM/3*mtx%*%diag(cefDM)%*%bb
X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
beta <- Variable(ncol(q2) + ncol(bb0)*(1+ncol(add_x)))
for (jj in 1:ncol(add_x)){
  X <- cbind(X,add_x[,jj]%*%bbtheta)
}
obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
for (i in 1:(lt-1))
{
  bb <- B.est[z[,2]==tau[i+1],]
  # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
  p <- hDM/3*mtx%*%diag(cefDM)%*%bb
  X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
  bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
  if (is.null(add_x)){
    obj <- obj + sum(quant_loss(daty - X%*% beta,tau[i+1]))
  }else{
    for (jj in 1:ncol(add_x)){
      X <- cbind(X,add_x[,jj]%*%bbtheta)
    }
    obj <- obj + sum(quant_loss(daty - X %*% beta, tau[i+1]))
  }
}

gamma <- cand2[optid]
roughness <- quad_form(beta[1:ncol(q2)],rp)*length(daty)
print(c('optimal gamma:', gamma))
prob <- Problem(Minimize(obj+gamma*roughness))
fit1 <- psolve(prob)
# obtain the initial estimate for the bivariate slope function
g1 <- fit1$getValue(beta)
saveRDS(object = g1, file=paste0('estg1_',tmpid,'.RDS'))
rm(prob)
rm(beta)
rm(obj)
rm(roughness)
rm(fit1)

detach(package:CVXR,unload = TRUE)
library(CVXR)

# 1st option of sparse penalty
int_mtx <- diag(1,nbern)

# 2nd option of sparse penalty
int_mtx <- matrix(c(6,3,3,1,1,1,3,4,2,3,2,1,3,2,4,1,2,3,1,3,1,6,3,1,1,2,2,3,4,3,1,1,3,1,3,6),nrow=6,ncol=6,byrow=T)

# calculate the weights of group lasso type penalty by treating each triangle of the triangulation as a group
lambda <- c()
for (k in 1:ngrp){
  lambda <- c(lambda, (g1[1:ncol(q2),]%*%t(q2[1:nbern + (k-1)*nbern,])%*%int_mtx%*%q2[1:nbern + (k-1)*nbern,]%*%g1[1:ncol(q2),])^0.5)
}
alpha <- -1

# Fit the model based on given tuning parameter values, f
l1 <- 2
l2 <- 2
scale <- cand1[l1]
gamma <- cand2[l2]
bb <- B.est[z[,2]==tau[1],]
p <- hDM/3*mtx%*%diag(cefDM)%*%bb
X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
beta <- Variable(ncol(q2) + ncol(bb0)*(1+ncol(add_x)))
for (jj in 1:ncol(add_x)){
  X <- cbind(X,add_x[,jj]%*%bbtheta)
}
obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
for (i in 1:(lt-1))
{
  bb <- B.est[z[,2]==tau[i+1],]
  # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
  p <- hDM/3*mtx%*%diag(cefDM)%*%bb
  X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
  bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
  for (jj in 1:ncol(add_x)){
    X <- cbind(X,add_x[,jj]%*%bbtheta)
  }
  obj <- obj + sum(quant_loss(daty - X %*% beta,tau[i+1]))
}

# adaptive group lasso penalties
grp_lasso <- lambda[1]^alpha*norm(int_mtx%*%q2[1:nbern,]%*%beta[1:ncol(q2)],type="2")
for (k in 2:ngrp){
  grp_lasso <- grp_lasso + lambda[k]^alpha*norm(int_mtx%*%q2[1:nbern + (k-1)*nbern,]%*%beta[1:ncol(q2)],type="2")
}

# roughness penalties
roughness <- quad_form(beta[1:ncol(q2)],rp)*length(daty)
prob <- Problem(Minimize(obj+ scale*nobs*grp_lasso + gamma*roughness))
fit2 <- psolve(prob)
tmpg <- fit2$getValue(beta)
saveRDS(object = tmpg, file=paste0('estg_',tmpid,'.RDS'))
rm(prob)
rm(beta)
rm(obj)
rm(grp_lasso)
rm(roughness)
rm(fit2)
detach(package:CVXR,unload = TRUE)
library(CVXR)

# identify the active and inactive regions
b1 <- q2%*%tmpg[1:ncol(q2)]
norm1 <- sum(b1[1:nbern,]^2)^0.5
for (k in 2:ngrp){
  norm1 <- c(norm1, sum(b1[1:nbern + (k-1)*nbern,]^2)^0.5)
}
print(c('# of zero triangles:',sum(norm1<=1e-3)))
inactive1 <- rep(norm1, each = nbern)
if (sum(norm1>=1e-3)>0){
  if (sum(inactive1<=1e-3)>0){
    indicator1 <-  apply(abs(H1[,inactive1<=1e-3]),1,sum)
    indicator <-  apply(abs(H[,inactive1<=1e-3]),1,sum)
    tmpH <- H[indicator==0,]
    tmpHH <- rbind(tmpH,H1[indicator1!=0,])
    for (j in 1:length(inactive1)){
      if (inactive1[j]<=1e-3){
        tmpinactive <- rep(0,length(inactive1))
        tmpinactive[j] <- 1
        tmpHH <- rbind(tmpHH,tmpinactive)
      }
    }
    Q22 <- qrH(tmpHH)
    rp1 <- t(Q22)%*%K%*%Q22
    z1 <- sum(norm1>=1e-3)
  }else{
    Q22 <- q2
    rp1 <- t(Q22)%*%K%*%Q22
    z1 <- 0
  }
}else{z1 <- NULL}
# Use cross validation to decide the optimal tuning parameter value for roughness penalty when refitting the model only using data witin signal regions
output <- c()
for (l3 in 1:length(cand3)){
  cvloss <- 0
  lc <- nobs/ncv
  for (cv in 1:ncv){
    tind <- (cv-1)*lc + c(1:lc)
    trind <- c(1:length(daty))[-tind]
    bb <- B.est[z[,2]==tau[1],]
    # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
    p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
    bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    
    X <- cbind(p%*%Q22,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
    beta <- Variable(ncol(Q22) + ncol(bb0))
    
    if (is.null(add_x)){
      print('no add_x')
      obj <- sum(quant_loss(daty[trind] - X%*% beta,tau[1])) 
    }else{
      beta <- Variable(ncol(Q22) + ncol(bb0)*(1+ncol(add_x)))
      for (jj in 1:ncol(add_x)){
        X <- cbind(X,add_x[trind,jj]%*%bbtheta)
      }
      obj <-sum(quant_loss(daty[trind]-X%*%beta,tau[1]))
    }
    
    for (i in 1:(length(tau)-1))
    {
      bb <- B.est[z[,2]==tau[i+1],]
      # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
      p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
      bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
      X <- cbind(p%*%Q22,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
      if (is.null(add_x)){
        obj <- obj + sum(quant_loss(daty[trind] - X %*% beta,tau[i+1])) 
      }else{
        for (jj in 1:ncol(add_x)){
          X <- cbind(X,add_x[trind,jj]%*%bbtheta)
        }
        obj <- obj + sum(quant_loss(daty[trind] - X %*% beta,tau[i+1]))
      }
    }
    
    roughness <- quad_form(beta[1:ncol(Q22)],rp1)*length(trind)
    prob <- Problem(Minimize(obj+cand3[l3]*roughness))
    fit <- psolve(prob)
    print(c('l1,l2,l3,cv_id:',l1,l2,l3,cv))
    g <- fit$getValue(beta)
    rm(prob)
    rm(beta)
    rm(obj)
    rm(roughness)
    rm(fit)
    
    # calculate loss in testing data
    bb <- B.est[z[,2]==tau[1],]
    # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
    lp <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
    lX <- cbind(lp%*%Q22,matrix(bb0[1,],nrow=nrow(lp),ncol=ncol(bb0),byrow=T))
    bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    if (is.null(add_x)){
      cvres <- daty[tind]- lX%*%g
    }else{
      for (jj in 1:ncol(add_x)){
        lX <- cbind(lX,add_x[tind,jj]%*%bbtheta)
      }
      cvres <- daty[tind]- lX%*%g
    }
    cvloss <- cvloss + sum(quant_loss(cvres,tau[1]))
    for (i in 1:(length(tau)-1))
    {
      bb <- B.est[z[,2]==tau[i+1],]
      # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
      lp <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
      lX <- cbind(lp%*%Q22,matrix(bb0[i+1,],nrow=nrow(lp),ncol=ncol(bb0),byrow=T))
      bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
      if (is.null(add_x)){
        cvres <-  daty[tind]- lX%*%g
      }else{
        for (jj in 1:ncol(add_x)){
          lX <- cbind(lX,add_x[tind,jj]%*%bbtheta)
        }
        cvres <-  daty[tind]- lX%*%g
      }
      cvloss <- cvloss + sum(quant_loss(cvres,tau[i+1]))
    }
  }
  output <- rbind(output,c(l3,cvloss,sum(norm1<=1e-3)))
}

# Refit the model by only using the information within the estimated active region
optid <- which(output[,2]==min(output[,2]))[1]
l3 <- output[optid,1]
g <- readRDS(paste0('estg_',tmpid,'.RDS'))

b1 <- q2%*%g[1:ncol(q2)]
norm1 <- sum(b1[1:nbern,]^2)^0.5
for (k in 2:ngrp){
  norm1 <- c(norm1, sum(b1[1:nbern + (k-1)*nbern,]^2)^0.5)
}
print(sum(norm1<=1e-3))
inactive1 <- rep(norm1, each = nbern)
if (sum(norm1>=1e-3)>0){
  if (sum(inactive1<=1e-3)>0){
    indicator1 <-  apply(abs(H1[,inactive1<=1e-3]),1,sum)
    indicator <-  apply(abs(H[,inactive1<=1e-3]),1,sum)
    # tmpH <- H[indicator==0,]
    # tmpHH <- rbind(tmpH,H1[indicator1!=0,])
    # for (j in 1:length(inactive1)){
    #   if (inactive1[j]<=1e-3){
    #     tmpinactive <- rep(0,length(inactive1))
    #     tmpinactive[j] <- 1
    #     tmpHH <- rbind(tmpHH,tmpinactive)
    #   }
    # }
    Q22 <- qrH(tmpHH)
    rp1 <- t(Q22)%*%K%*%Q22
    z1 <- sum(norm1>=1e-3)
  }else{
    Q22 <- q2
    rp1 <- t(Q22)%*%K%*%Q22
  }
  saveRDS(object = Q22, file=paste0('Q22_',tmpid,'.RDS'))
}else{z1 <- NULL}

if (is.null(z1)==0){
  bb <- B.est[z[,2]==tau[1],]
  # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
  p <- hDM/3*mtx%*%diag(cefDM)%*%bb
  X <- cbind(p%*%Q22,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
  bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
  beta <- Variable(ncol(Q22) + ncol(bb0))
  roughness <- quad_form(beta[1:ncol(Q22)],rp1)*length(daty)
  if (is.null(add_x)){
    obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
  }else{
    theta <- Variable(rows = ncol(bb0),cols = ncol(add_x))
    obj <- sum(quant_loss(daty - X%*% beta -  add_x%*%t(bbtheta%*%theta),tau[1]))
  }
  
  for (i in 1:(length(tau)-1))
  {
    bb <- B.est[z[,2]==tau[i+1],]
    # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
    p <- hDM/3*mtx%*%diag(cefDM)%*%bb
    X <- cbind(p%*%Q22,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
    bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
    if (is.null(add_x)){
      obj <- obj + sum(quant_loss(daty - X %*% beta,tau[i+1])) 
    }else{
      obj <- obj + sum(quant_loss(daty - X %*% beta - add_x%*%t(bbtheta%*%theta),tau[i+1]))
    }
  }
  prob <- Problem(Minimize(obj+cand3[l3]*roughness))
  fit <- psolve(prob)
  g <- fit$getValue(beta)
  if (is.null(add_x)==0){
    gg <- fit$getValue(theta)
    saveRDS(object = gg, file=paste0('fgg_',tmpid,'.RDS'))
  }
  saveRDS(object = g, file=paste0('fg_',tmpid,'.RDS'))
}

# plot the true slope
plot(t_grid,rhos(t_grid),col='red',type='l')
# plot the estimated slope function across quantiles
lines(t_grid,B.est[z[,2]==tau[9],]%*%Q22%*%g[1:ncol(Q22),],type='l')

