# The logistic model to be used to cumulative cases
model_logistic <- expression(
    K * N / (1 + (K*N/x0-1) * exp(-lambda*times))
    )

model_logistic_baseline <- expression(
    B + K * N / (1 + (K*N/x0-1) * exp(-lambda*times))
    )

# The Richards model to be used to cumulative cases
model_richards <- expression(
  K * N / (1 + ((K*N/x0)^alpha-1) * 
             exp(-lambda / (1 - (x0 / K / N) ^ alpha) *alpha*times)) ^(1/alpha) 
)

model_richards_baseline <- expression(
  B + K * N / (1 + ((K*N/x0)^alpha-1) * 
             exp(-lambda / (1 - (x0 / K / N) ^ alpha) *alpha*times)) ^(1/alpha) 
)

# The exponential growth model
model_exp <- expression(x0 * exp(lambda * times))
 
convolute.exp <- function(X, mu) {
    lX = length(X)
    lk = max(min(lX, round(10/mu)),3)
    kern = diff(pexp(0:lk, rate=mu))
    return(filter(c(rep(0, lk-1), X), kern, sides=1)[(1:lX)+lk-1])
}
 
model_logistic_delay <- expression(
    convolute.exp(eval(model_logistic), mu)
)
 
derivatives <- function(Rcum, pnames) {
    deriv(Rcum, pnames, function.arg=c(pnames,"N","times"))
}
 
deriv.logistic <- derivatives(model_logistic, c('x0', 'lambda', 'K'))
deriv.logistic.baseline <- derivatives(model_logistic_baseline, 
                                      c('x0', 'lambda', 'K', 'B'))

deriv.richards <- derivatives(model_richards, 
                              c('x0', 'lambda', 'K', 'alpha'))

deriv.richards.baseline <- derivatives(model_richards_baseline, 
                              c('x0', 'lambda', 'K', 'alpha', 'B'))

derivatives.logistic.delay <- function(x0, lambda, K, mu, N, times) {
    d <- attr(deriv.logistic(x0, lambda, K, N, times), 'gradient')
    X = convolute.exp(eval(model_logistic_delay), mu)
    grad = array(0, c(length(X), 4), list(NULL, c('x0', 'lambda', 'K', 'mu')))
    grad[, 1:3] = apply(d, 2, convolute.exp, mu=mu)
    grad[, 'mu'] = X*mu
    attr(X, 'gradient') <- grad
    return(X)
}
 
## deriv of P(x,lambda) = -sum(dpois(x,lambda,log=TRUE)) 
##     wrt lambda == sum(1-lambda/x) = N - lambda*(sum(1/x))
## deriv of P(x,lambda) wrt p = dP/d(lambda) * d(lambda)/dp
## compute gradient vector
gradlikfun <- function(p, derivatives, dat,times,N,incid=TRUE) {
    gcall <- do.call(derivatives,c(as.list(p),list(times=times,N=N))) ## values + gradient matrix
    lambda <- gcall
    attr(lambda,"gradient") <- NULL
    if (incid) lambda <- diff(lambda)
    gmat <- attr(gcall,"gradient") ## extract gradient
    if (incid) gmat <- apply(gmat,2,diff)  ## differences
    ## apply chain rule (multiply columns of gmat by dP/dlambda)
    totderiv <- sweep(gmat,MARGIN=1,(1-dat/lambda),"*") 
    ## deriv of summed likelihood = sum of derivs of likelihood
    csums <- colSums(totderiv)
    csums[!is.finite(csums)] <- -1 ## ????
    csums
}
