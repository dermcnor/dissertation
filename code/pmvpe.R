pmvpe <- function(x,mu=rep(0,p),sigma=diag(p),b=1){
# MV Power Exponential probability density function (pdf).
# Y = pmvpe(X,MU,SIGMA,B) Returns the multivariate pe pdf with mean, MU, 
# covariance matrix, SIGMA, and kurtosis parameter, B,
# at the values in X. 
# Default values for MU, SIGMA, and B are [0 0]', [1 0;0 1], and 1 respectively.

n,p] = size(x);

y = zeros(n,p);


xn = -.5*abs((x-repmat(mu,n,1))/sigma)^(2*b);
yn = 1/((sigma)*gamma(1+.5/b)*2^(1+.5/b));
y  = exp( xn ) .* yn;
}