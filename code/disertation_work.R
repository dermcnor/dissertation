require(lattice)
require(RColorBrewer)
source("dmvpe.R")
x <- seq(-10,10,.5)
a <- expand.grid(x=x,y=x)
a$z <- dmvpe(expand.grid(x,x),beta=2,sigma=matrix(c(10, 0, 0, 10),ncol=2))
z <- matrix(a$z,ncol=length(x))
nrz <- nrow(z)
ncz <- ncol(z)
nbcol <- 20
color <- terrain.colors(nbcol)
color2 <- colorRampPalette(rev(brewer.pal(5, "YlOrBr")), space = "Lab")(nbcol)
color3 <- colorRampPalette(rev(brewer.pal(5, "Blues")), space = "Lab")(nbcol)

zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
persp(x=x, y=x, z=z, col=color[facetcol], phi=30, theta=-30, box = FALSE)

filled.contour(z,color.palette = colorRampPalette(rev(brewer.pal(5, "YlOrBr")), space = "Lab"))
filled.contour(z,color.palette = colorRampPalette(rev(brewer.pal(5, "Blues")), space = "Lab"))



ll <- function(theta, x){
		mu <- theta[1:2]
		covmat <- matrix(c(theta[3:4],0,theta[5,ncol=2)
		beta1 <- theta[6]
        -sum(dmvpe(x, mean=mu, sigma=covmat, beta=beta1, log=TRUE))
}



x <- rnorm(100,mean=3,sd=10)
x <- cbind(x, 2*x+rnorm(100))
nlm(ll,c(3,7,95,189,380,1),x=x)
a <- chol(matrix(c(95,189,189,380),ncol=2))
t(a)%*%a

optim(c(c(apply(x,2,mean)),cov(x)[upper.tri(diag(2),diag=TRUE)],1),ll,x=x,method = "L-BFGS-B", lower=c(rep(-Inf,2),0,-Inf,0,0))