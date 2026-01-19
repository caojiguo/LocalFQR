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
t_grid <- newgrids <- seq(0,1,length.out = ncol(mtx))
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
ind <- Bfull.est$Ind.inside
z <- z[ind,]
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs
cefDM <- c(1, rep(c(4,2), (Mobs-3)/2), 4,  1)

knots <- quantile(tau,seq(0,1,length=10))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
basisobj=create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)

n <- nrow(mtx)
y <- daty

ngrp <- nrow(q2)/nbern
quant_loss <- function(u, t) {0.5 * abs(u) + (t - 0.5) * u}
rp <-t(q2)%*%K%*%q2

cand <- c(1e-6,1e-5,1e-4,1e-3,1e-2)
gammacv <- readRDS('gammacv.RDS')
optid <- which(gammacv[,2]==min(gammacv[,2]))[1]
print(gammacv[optid,])

beta <- Variable(ncol(q2) + ncol(bb0))
theta <- Variable(rows = ncol(bb0),cols = ncol(add_x))
roughness <- quad_form(beta[1:ncol(q2)],rp)
bb <- B.est[z[,2]==tau[1],]
p <- hDM/3*mtx%*%diag(cefDM)%*%bb
X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
obj <- sum(quant_loss(y - X%*% beta - add_x%*%t(bbtheta%*%theta),tau[1]))
for (i in 1:(lt-1))
{
  bb <- B.est[z[,2]==tau[i+1],]
  # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
  p <- hDM/3*mtx%*%diag(cefDM)%*%bb
  X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
  bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
  obj <- obj + sum(quant_loss(y - X %*% beta - add_x%*%t(bbtheta%*%theta),tau[i+1]))
}

gamma <- cand[gammacv[optid,1]]*n
prob <- Problem(Minimize(obj+gamma*roughness))
fit1 <- solve(prob)
g1 <- fit1$getValue(beta)

int_mtx <- diag(1,nbern)

lambda <- c()
for (k in 1:ngrp){
  lambda <- c(lambda, (g1[1:ncol(q2),]%*%t(q2[1:nbern + (k-1)*nbern,])%*%int_mtx%*%q2[1:nbern + (k-1)*nbern,]%*%g1[1:ncol(q2),])^0.5)
}

alpha <- -1

cand1 <- c(1e-4,1e-3,1e-2,1e-1,1)
cand3 <- cand2 <- c(1e-6,1e-5,1e-4,1e-3,1e-2)

for (l1 in 1:length(cand1)){
  for (l2 in 1:length(cand2)){
    scale <- cand1[l1]*n
    gamma <- cand2[l2]*n
    bb <- B.est[z[,2]==tau[1],]
    p <- hDM/3*mtx%*%diag(cefDM)%*%bb
    X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
    bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
    beta <- Variable(ncol(q2) + ncol(bb0))
    theta <- Variable(rows = ncol(bb0),cols = ncol(add_x))
    
    # roughness penalties
    roughness <- quad_form(beta[1:ncol(q2)],rp)
    # weights for group lasso type penalty
    grp_lasso <- lambda[1]^alpha*norm(int_mtx%*%q2[1:nbern,]%*%beta[1:ncol(q2)],type="2")
    for (k in 2:ngrp){
      grp_lasso <- grp_lasso + lambda[k]^alpha*norm(int_mtx%*%q2[1:nbern + (k-1)*nbern,]%*%beta[1:ncol(q2)],type="2")
    }
    
    obj <- sum(quant_loss(daty - X%*% beta -  add_x%*%t(bbtheta%*%theta),tau[1]))
    for (i in 1:(lt-1))
    {
      bb <- B.est[z[,2]==tau[i+1],]
      # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
      p <- hDM/3*mtx%*%diag(cefDM)%*%bb
      X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
      bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
      obj <- obj + sum(quant_loss(daty - X %*% beta - add_x%*%t(bbtheta%*%theta),tau[i+1]))
    }
    prob <- Problem(Minimize(obj+ scale*grp_lasso + gamma*roughness))
    fit2 <- psolve(prob)
    g <- fit2$getValue(beta)
    gg <- fit2$getValue(theta)
    
    
    # identify the active regions
    b1 <- q2%*%g[1:ncol(q2)]
    norm1 <- sum(b1[1:nbern,]^2)^0.5
    for (k in 2:ngrp){
      norm1 <- c(norm1, sum(b1[1:nbern + (k-1)*nbern,]^2)^0.5)
    }
    
    print(sum(norm1<=1e-3))
    inactive1 <- rep(norm1, each = nbern)
    
    if (sum(norm1>=1e-3)>0){
      if (sum(inactive1<=1e-3)>0){
        qs <- q2[inactive1<=1e-3,]
        q21 <- qrH(qs)
        q21 <- as.matrix(q21)
      }else{
        q21 <- diag(1,nrow=ncol(q2))
      }
      rp <- t(q21)%*%t(q2)%*%K%*%q2%*%q21
      z1 <- sum(norm1>=1e-3)
    }else{z1 <- NULL}
    
    # cross-validation for roughness penalties based on the selected regions
    output <- c()
    print('cross-validation starts')
    for (l3 in 1:length(cand3)){
      cvloss <- 0
      ncv <- 10
      for (cv in 1:ncv){
        lc <- length(daty)/ncv
        lc <- floor(lc)
        tind <- (cv-1)*lc + c(1:lc)
        tind <- tind[tind<=nrow(mtx)]
        trind <- c(1:length(daty))[-tind]
        if (is.null(z1)==0){
          bb <- B.est[z[,2]==tau[1],]
          # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
          p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
          X <- cbind(p%*%q2%*%q21,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
          bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
          beta <- Variable(ncol(q21)+ ncol(bb0))
          theta <- Variable(rows = ncol(bb0),cols = ncol(add_x))
          roughness <- quad_form(beta[1:ncol(q21)],rp)
          obj <- sum(quant_loss(daty[trind] - X%*% beta -  add_x[trind,]%*%t(bbtheta%*%theta),tau[1]))
        }
        if (is.null(z1)){
          X <- matrix(bb0[1,],nrow=length(trind),ncol=ncol(bb0),byrow=T)
          beta <- Variable(ncol(bb0))
          bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
          obj <- sum(quant_loss(daty[trind] -  X%*% beta - add_x[trind,]%*%t(bbtheta%*%theta),tau[1]))
        }
        for (i in 1:(length(tau)-1))
        {
          if (is.null(z1)==0){
            bb <- B.est[z[,2]==tau[i+1],]
            # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
            p <- hDM/3*mtx[trind,]%*%diag(cefDM)%*%bb
            X <- cbind(p%*%q2%*%q21,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
            bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
            obj <- obj + sum(quant_loss(daty[trind] - X %*% beta - add_x[trind,]%*%t(bbtheta%*%theta),tau[i+1]))
          }
          if (is.null(z1)){
            X <- matrix(bb0[i+1,],nrow=length(trind),ncol=ncol(bb0),byrow=T)
            beta <- Variable(ncol(bb0))
            bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
            obj <- obj + sum(quant_loss(daty[trind] -X%*% beta- add_x[trind,]%*%t(bbtheta%*%theta),tau[i+1]))
          }
        }
        if (is.null(z1)==0){
          prob <- Problem(Minimize(obj+cand3[l3]*n*roughness))
          fit <- psolve(prob,solver='ECOS')
          g <- fit$getValue(beta)
          gg <- fit$getValue(theta)
        }
        if (is.null(z1)){
          prob <- Problem(Minimize(obj))
          fit <- psolve(prob,solver='ECOS')
          g <- fit$getValue(beta)
          gg <- fit$getValue(theta)
        }
        
        # calculate loss in testing dataset
        if (is.null(z1)==0){
          bb <- B.est[z[,2]==tau[1],]
          # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
          lp <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
          lX <- cbind(lp%*%q2%*%q21,matrix(bb0[1,],nrow=nrow(lp),ncol=ncol(bb0),byrow=T))
          bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
          cvres <- daty[tind]- lX%*%g - add_x[tind,]%*%t(bbtheta%*%gg)
        }
        if (is.null(z1)){
          lX <- matrix(bb0[1,],nrow=length(tind),ncol=ncol(bb0),byrow=T)
          bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
          cvres <- daty[tind]- lX%*%g - add_x[tind,]%*%t(bbtheta%*%gg)
        }
        cvloss <- cvloss + sum(quant_loss(cvres,tau[1]))
        
        for (i in 1:(length(tau)-1))
        {
          if (is.null(z1)==0){
            bb <- B.est[z[,2]==tau[i+1],]
            # Using Simpson rule to approximate the integral of x(t) * b_k(u,t) over the domain of t for given u, where b_k(u,t) represents Bernstein polynomials defined ove the triangulation
            lp <- hDM/3*mtx[tind,]%*%diag(cefDM)%*%bb
            lX <- cbind(lp%*%q2%*%q21,matrix(bb0[i+1,],nrow=nrow(lp),ncol=ncol(bb0),byrow=T))
            bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
            cvres <-  daty[tind]- lX%*%g - add_x[tind,]%*%t(bbtheta%*%gg)
          }
          
          if (is.null(z1)){
            lX <- matrix(bb0[i+1,],nrow=length(tind),ncol=ncol(bb0),byrow=T)
            bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
            cvres <-  daty[tind]- lX%*%g - add_x[tind,]%*%t(bbtheta%*%gg)
          }
          cvloss <- cvloss + sum(quant_loss(cvres,tau[i+1]))
        }
        print(cvloss)
      }
      output <- rbind(output,c(l1,l2,l3,cvloss,sum(norm1<=1e-3)))
    }
  }
}

saveRDS(object = output, file=paste0('output.RDS'))
