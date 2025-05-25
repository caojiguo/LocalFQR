library(e1071)

# time grid
t_grid <- newgrids <- seq(0,1,length.out = 400)


rho1 <- function(t) 20*(t-0.5)^2*(t<=0.5)
rho2 <- function(t) (t-0.5)^2*(t<=0.5)
rhos <- function(t) 20*(t-0.5)^2*(t<=0.5) +(t-0.5)^2*(t<=0.5)*qnorm(u)

mtx <- matrix(0,nobs,length(t_grid))
add_x <- matrix(rnorm(nobs*2,0,0.05),nrow=nobs)
daty <- c()
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs
cefDM   = c(1, rep(c(4,2), (Mobs-2)/2), 1)

for (j in 1:nobs)
{
  X <- rwiener(end=1,frequency = length(t_grid))
  mtx[j,] <- c(X)
}
err <- rnorm(nobs,0,1)
# Use Simpson rule to approximate the integral of x(t) * rho1(t)
tmpy <- hDM/3*mtx%*%diag(cefDM)%*%rho1(t_grid)
# Use Simpson rule to approximate the integral of (x(t)+4) * rho2(t)
sigmax <- hDM/3*(mtx+4)%*%diag(cefDM)%*%rho2(t_grid)
daty <- tmpy + sigmax*err + add_x%*%c(1,1)
saveRDS(object = mtx, file='mtx.RDS')
saveRDS(object = add_x, file='add_x.RDS')
saveRDS(object = daty, file='daty.RDS')
