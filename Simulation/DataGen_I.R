library(e1071)

# time grid
t_grid <- newgrids <- seq(0,1,length.out = 399)

rhos <- function(t) 20*(sin(t*4*pi)*(t-0.5)^2)*(t<=0.5)

mtx <- matrix(0,nobs,length(t_grid))
add_x <- matrix(rnorm(nobs*2,0,0.05),nrow=nobs)
daty <- c()
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs
cefDM   = c(1, rep(c(4,2), (Mobs-3)/2), 4, 1)
#set.seed(id)
for (j in 1:nobs)
{
  X <- rwiener(end=1,frequency = length(t_grid))
  mtx[j,] <- c(X)
}
err <- rnorm(nobs,0,0.05)
# Use Simpson rule to approximate the integral of x(t) * rhos(t)
daty <- hDM/3*mtx%*%diag(cefDM)%*%rhos(t_grid)
daty <- daty + add_x%*%c(1,1) + err
