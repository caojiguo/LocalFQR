library(stats)
library(BPST)
library(Matrix)
library(MASS)
library(tensr)
library(CVXR)
library(fda)
library(Rmosek)
library(plotly)

# load the cross validation output
output <- readRDS('cv_output.RDS')
# find the optimal tuning parameter values
bestid <- which(output[,4]==min(output[,4]))

# quantiles of interest
tau <- seq(0.25,0.75,length.out = 17)
lt <- length(tau)
# number of knots for B-splines
nk <- 10
# number of folds for cross validation
ncv <- 10
# degree of Bernstein polynomials
d <- 2
# smoothness conditions for Bernstein polynomial approximation over the triangulation
r <- 1

# candiate values for tuning parameters
cand1 <- c(exp(6),exp(5),exp(4))
cand3 <- cand2 <- c(exp(0),exp(-2),exp(-4),exp(-6),exp(-8),exp(-10))

# load the soybean yield data
tmpx <- readRDS('soy_avg_x.RDS')
daty <- readRDS('soy_avg_y.RDS')
add_var <- readRDS('soy_add_avg.RDS')
add_var[,3] <- add_var[,3]/add_var[,2]
add_x <- add_var[,-2]
names(add_x) <- c('prcp','irr_ratio')
# minimum temperature
mtx2 <- tmpx[,30:335]
# maximum temperature
mtx1 <- tmpx[,395:700]  
# average temperature
mtx <- (mtx1 + mtx2)/2
t_grid <- seq(0,1,length.out = ncol(mtx))
Y = daty
# input the triangulation
V.est <- readRDS(file = 'nodes.RDS')
Tr.est <- readRDS(file='tri.RDS')
# number of triangles
ntr <- nrow(Tr.est)
# number of basis functions on each triangle
nbern <- choose(d+2,2)
z <- cbind(rep(t_grid,each=length(tau)),rep(tau,length(t_grid)))
# compute the values of Bernstein polynomials based on the given triangulation
Bfull.est <- basis(V.est,Tr.est,d,r,z)
H1 <- smoothness(V.est, Tr.est, d, r-1)
H1 <- as.matrix(H1)
H <- smoothness(V.est, Tr.est, d, r)
H <- as.matrix(H)
B.est <- Bfull.est$B
# matrix derived from QR decomposition
q2 <- Bfull.est$Q2
# symmetric and positive semi-definite matrix corresponding to the roughness penalty
K <- as.matrix(Bfull.est$K)
ind <- Bfull.est$Ind.inside
z <- z[ind,]
tt <- z[z[,2]==tau[1],1]
l <- length(tt)
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs

# Define the B-splines
knots <- quantile(tau,seq(0,1,length=nk))
norder= 2
nknots = length(knots)
nb = nknots + norder - 2
basisobj=create.bspline.basis(rangeval = c(min(tau), max(tau)), nbasis = nb, norder = norder,breaks = knots)
bb0 <- getbasismatrix(tau, basisobj)
n <- nrow(mtx)
y <- daty
ngrp <- nrow(q2)/nbern

# define quantile loss 
quant_loss <- function(u, t) {0.5 * abs(u) + (t - 0.5) * u}
rp <-t(q2)%*%K%*%q2

# load the cross validation output regarding the tuning parameter corresponding to the roughness penalty when we seek for the initial estimate of the slope function
gammacv <- readRDS('gammacv.RDS')
optid <- which(gammacv[,2]==min(gammacv[,2]))[1]

beta <- Variable(ncol(q2) + ncol(bb0))
bb <- B.est[z[,2]==tau[1],]
p <- hDM*mtx%*%bb
X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
if (is.null(add_x)){
  obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
}else{
  beta <- Variable(ncol(q2) + ncol(bb0)*(1+ncol(add_x)))
  if (ncol(add_x)==1){
    tmpb <- bbtheta
  }else{
    tmpb <- bbtheta
    for (jj in 2:ncol(add_x)){
      tmpb <- bdiag(tmpb,bbtheta)
    }
  }
  X <- cbind(X,add_x%*%tmpb)
  obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
}
for (i in 1:(lt-1))
{
  bb <- B.est[z[,2]==tau[i+1],]
  p <- hDM*mtx%*%bb
  X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
  bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
  if (is.null(add_x)){
    obj <- obj + sum(quant_loss(daty - X%*% beta,tau[i+1]))
  }else{
    tmpb <- bbtheta
    for (jj in 2:ncol(add_x)){
      tmpb <- bdiag(tmpb,bbtheta)
    }
    X <- cbind(X,add_x%*%tmpb)
    obj <- obj + sum(quant_loss(daty - X %*% beta, tau[i+1]))
  }
}
gamma <- cand2[optid]
roughness <- quad_form(beta[1:ncol(q2)],rp)
prob <- Problem(Minimize(obj+gamma*roughness))

# find the initial estimate for the slope function
fit1 <- psolve(prob,solver="MOSEK")
g1 <- fit1$getValue(beta)

# set up the sparse penalty based on the initial estimates
int_mtx <- matrix(c(6,3,3,1,1,1,3,4,2,3,2,1,3,2,4,1,2,3,1,3,1,6,3,1,1,2,2,3,4,3,1,1,3,1,3,6),nrow=6,ncol=6,byrow=T)

lambda <- c()
for (k in 1:ngrp){
  lambda <- c(lambda, (g1[1:ncol(q2),]%*%t(q2[1:nbern + (k-1)*nbern,])%*%int_mtx%*%q2[1:nbern + (k-1)*nbern,]%*%g1[1:ncol(q2),])^0.5)
}
alpha <- -1

# weights for group lasso penalties
grp_lasso <- lambda[1]^alpha*norm(int_mtx%*%q2[1:nbern,]%*%beta[1:ncol(q2)],type="2")
for (k in 2:ngrp){
  grp_lasso <- grp_lasso + lambda[k]^alpha*norm(int_mtx%*%q2[1:nbern + (k-1)*nbern,]%*%beta[1:ncol(q2)],type="2")
}

# fit the model based on the initial estimate derived above and the optimal tuning parameter values
l1 <- output[bestid,1]
l2 <- output[bestid,2]
scale <- cand1[l1]
gamma <- cand2[l2]
bb <- B.est[z[,2]==tau[1],]
p <- hDM*mtx%*%bb
X <- cbind(p%*%q2,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
beta <- Variable(ncol(q2) + ncol(bb0))

if (is.null(add_x)){
  obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
}else{
  beta <- Variable(ncol(q2) + ncol(bb0)*(1+ncol(add_x)))
  for (jj in 1:ncol(add_x)){
    X <- cbind(X,add_x[,jj]%*%bbtheta)
  }
  obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
}

for (i in 1:(lt-1))
{
  bb <- B.est[z[,2]==tau[i+1],]
  p <- hDM*mtx%*%bb
  X <- cbind(p%*%q2,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
  bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
  if (is.null(add_x)){
    obj <- obj + sum(quant_loss(daty - X%*% beta,tau[i+1]))
  }else{
    for (jj in 1:ncol(add_x)){
      X <- cbind(X,add_x[,jj]%*%bbtheta)
    }
    obj <- obj + sum(quant_loss(daty - X %*% beta,tau[i+1]))
  }
}

roughness <- quad_form(beta[1:ncol(q2)],rp)
prob <- Problem(Minimize(obj+ scale*grp_lasso + gamma*roughness))
fit2 <- psolve(prob,solver="MOSEK")
g <- fit2$getValue(beta)

b1 <- q2%*%g[1:ncol(q2)]
norm1 <- sum(b1[1:nbern,]^2)^0.5
for (k in 2:ngrp){
  norm1 <- c(norm1, sum(b1[1:nbern + (k-1)*nbern,]^2)^0.5)
}


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


# refit the model only using the information within the selected active regions
l3 <- output[bestid,3]
if (is.null(z1)==0){
  bb <- B.est[z[,2]==tau[1],]
  p <- hDM*mtx%*%bb
  X <- cbind(p%*%Q22,matrix(bb0[1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
  bbtheta <- matrix(bb0[1,],nrow=1,ncol=ncol(bb0),byrow=T)
  beta <- Variable(ncol(Q22) + ncol(bb0))
  roughness <- quad_form(beta[1:ncol(Q22)],rp1)
  if (is.null(add_x)){
    obj <- sum(quant_loss(daty - X%*% beta,tau[1]))
  }else{
    theta <- Variable(rows = ncol(bb0),cols = ncol(add_x))
    obj <- sum(quant_loss(daty - X%*% beta -  add_x%*%t(bbtheta%*%theta),tau[1]))
  }
  
  for (i in 1:(length(tau)-1))
  {
    bb <- B.est[z[,2]==tau[i+1],]
    p <- hDM*mtx%*%bb
    X <- cbind(p%*%Q22,matrix(bb0[i+1,],nrow=nrow(p),ncol=ncol(bb0),byrow=T))
    bbtheta <- matrix(bb0[i+1,],nrow=1,ncol=ncol(bb0),byrow=T)
    if (is.null(add_x)){
      obj <- obj + sum(quant_loss(daty - X %*% beta,tau[i+1])) 
    }else{
      obj <- obj + sum(quant_loss(daty - X %*% beta - add_x%*%t(bbtheta%*%theta),tau[i+1]))
    }
  }
  prob <- Problem(Minimize(obj+cand3[l3]*roughness))
  fit <- psolve(prob,solver="MOSEK")
  g <- fit$getValue(beta)
  if (is.null(add_x)==0){
    gg <- fit$getValue(theta)
  }
}

# the following codes are used to generate the plots presented in the paper
betas <- c()
for (i in tau){
  betas <- rbind(betas,as.numeric(B.est[z[,2]==i,]%*%Q22%*%g[1:ncol(Q22),]))
}
fig <- plot_ly(x = c(30:335), y = tau, z = betas, colors='YlOrRd') %>% add_surface()

axx <- list(title = "Month",
            ticktext=list("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),
            tickvals = c(30,61,91,122,152,183,213,243,273,304,335))
axy <- list(title = "Quantile")
axz <- list(title = "Value")
fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig

k=1
plot(c(30:335), B.est[z[,2]==tau[k],]%*%Q22%*%g[1:ncol(Q22),],type='l',xlab='Month',ylab='',cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=2,col='red',ylim=c(-10,8))
axis(1, at=c(30,61,91,122,152,183,213,243,273,304,335), lab=c("Feb","Mar","Apr","May","June", "July", "Aug","Sep", "Oct", "Nov", "Dec"),cex.axis=1.5)
k=9
lines(c(30:335),B.est[z[,2]==tau[k],]%*%Q22%*%g[1:ncol(Q22),],cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=2,col='blue',lty=2)
k=17
lines(c(30:335),B.est[z[,2]==tau[k],]%*%Q22%*%g[1:ncol(Q22),],cex.lab=1.5,cex.axis=1.5,xaxt='n',lwd=2,lty=3)

legend(50, -1,legend=c(tau[1],tau[9],tau[17]),col=c('red','blue','black'), lwd=2 ,lty=c(1:3),cex=1.5,bty='n',x.intersp = 0.5, pt.cex=0.5)


