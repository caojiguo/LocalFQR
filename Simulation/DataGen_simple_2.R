library(e1071)

#tmpid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# set.seed(id)
# grids for T
t_grid <- newgrids <- seq(0,1,length.out = 400)
tmp2 <- '5*sin(2*pi*(t-0.25))*(t<=0.75)*(t>=0.25)'
text <- paste0('function(t)',tmp2)
# rho(t)
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
err <- rnorm(nobs,0,0.05)
daty <- hDM/3*mtx%*%diag(cefDM)%*%rhos(t_grid)
daty <- daty + add_x%*%c(1,1) + err
