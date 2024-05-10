library(e1071)

#tmpid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# set.seed(id)
# grids for T
t_grid <- newgrids <- seq(0,1,length.out = 400)

tmp1 <- '20*(t-0.5)^2*(t<=0.5)'
tmp2 <- '(t-0.5)^2*(t<=0.5)'
text <- paste0('function(t)',tmp1)
rho1 <- eval(parse(text=text))
text <- paste0('function(t)',tmp2)
rho2 <- eval(parse(text=text))
text <- paste0('function(t)','20*(t-0.5)^2*(t<=0.5) +(t-0.5)^2*(t<=0.5)*qnorm(u)')
rhos <- eval(parse(text=text))

mtx <- matrix(0,nobs,length(t_grid))
add_x <- matrix(rnorm(nobs*2,0,0.05),nrow=nobs)
daty <- c()
Mobs    = ncol(mtx)
tmin <- range(t_grid)[1]
tmax <- range(t_grid)[2]
hDM     = (tmax-tmin)/Mobs
cefDM   = c(1, rep(c(4,2), (Mobs-2)/2), 1)
#set.seed(id)
for (j in 1:nobs)
{
  X <- rwiener(end=1,frequency = length(t_grid))
  mtx[j,] <- c(X)
}
err <- rnorm(nobs,0,1)
tmpy <- hDM/3*mtx%*%diag(cefDM)%*%rho1(t_grid)
sigmax <- hDM/3*(mtx+4)%*%diag(cefDM)%*%rho2(t_grid)
daty <- tmpy + sigmax*err + add_x%*%c(1,1)
