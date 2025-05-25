library(fdaPDE)
dtau <- 0.05
border1 <- c(0,0)
n1 <- 4
n2 <- 10
for (i in 1:n1)
{
  border1 <- rbind(border1,c(0,i/n1))
}

border2 <- c()
for (i in 1:n2)
{
  border2 <- rbind(border2,c(0+i/n2,1))
}

border3 <- c()
for (i in 1:n1)
{
  border3 <- rbind(border3,c(1,1-i/n1))
}

border4 <- c()
for (i in 1:n2)
{
  border4 <- rbind(border4,c(1-i/n2,0))
}


border <- rbind(border1,border2,border3,border4[-nrow(border4),])
#plot(border)
l_b <- nrow(border)
vborder <- cbind(c(1:l_b),c(2:l_b,1))

nodes <- c()
for (j in 1:(n1-1))
{
  for (k in 1:(n2-1))
  {
    nodes <- rbind(nodes,c(1-k/n2,j/n1))
  }
}
l_in <- nrow(nodes)
nodes <- rbind(border,nodes)
nodes[,1] <- nodes[,1]
dat <- cbind(c(1:(l_in + l_b)),nodes)
dat[,3] <- dat[,3]*(max(tau)-min(tau) + dtau*2)+ min(tau)-dtau
newnodes <- dat[,c(2,3)]
unique(newnodes[,2])
mesh <- create.mesh.2D(nodes = newnodes, segments = vborder, order = 1)
saveRDS(object = newnodes,file = 'nodes.RDS')
saveRDS(object = mesh$triangles,file = 'tri.RDS')
