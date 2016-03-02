library('bbmle');
source("models.R")

# function exp_fit
# parameters:
#   - times, X: the time series
#   - theta0: the initial guess of the parameters
#   - start: the time as the start of the fitting window 
#   - is.cumulative: whether the time series is cumulaative incidence
exp_fit = function(times, X, theta0, start, is.cumulative, N=NULL) {
	start.time = times[start]
	# convert cumulative cases to new cases
	if (is.cumulative) X = pmax(diff(c(0,X)), 0)
	# fit from start to the peak
	m = which.max(X[start:length(X)])
	t = times[start-1+1:m] - start.time
	d = X[start-1+1:m]
	use.richards = exists("alpha", envir=as.environment(as.list(theta0)))
	if (use.richards) {
		model = model_richards
		derivatives = deriv.richards
	} else {
		model = model_logistic
		derivatives = deriv.logistic
	}
	
	# N is the total number of cases, used by some models 
	# to scale some parameters
	if(is.null(N)) N <- sum(d)
	
	#whether a fit is valid
	is.valid <- function(fit) {
		return (min(eigen(fit@details$hessian, only.values = TRUE)$values) > 0)
	}
	
	## equivalent (using model_*): 
	calc_mean <- function(p, times, N) {
		pp <- c(as.list(p),list(times=times,N=N))
		vals <- eval(model,envir=pp)
		if (is.cumulative) diff(c(0,vals)) else vals
	}
	
	likfun <- function(p, derivatives, dat,times, N) {
		lambda=calc_mean(p, times, N)
		lik = -sum(dpois(dat,lambda,log=TRUE))
		if (is.infinite(lik) || is.na(lik)) return(1e10)
		return(lik)
	}
	parnames(likfun) <- names(theta0)
	
	# sometimes confint may find a better solution. 
	# In this case we need to refit the model with the better
	# solution as a starting paoint, because the parameter
	# names in the returned better fit can be wrong
	while (TRUE) {
		if (is.function(derivatives)) {
			gr <- gradlikfun
		} else gr <- NULL
		# fit using the better answer between nlminb and BFGS
		fit.nlminb <- mle2(likfun, start=theta0,gr=gr,
						   data=list(derivatives = derivatives, times=t, dat=d, N=N),
						   optimizer='nlminb',
						   upper = Inf,
						   lower = 0,
						   vecpar=TRUE)
		fit.BFGS <- mle2(likfun, start=theta0,gr=gr,
						 data=list(derivatives = derivatives, times=t, dat=d, N=N),
						 method='L-BFGS-B',
						 upper = Inf,
						 lower = 0,
						 vecpar=TRUE)
		if (logLik(fit.nlminb) < logLik(fit.BFGS)) {
			fit = fit.BFGS
		} else fit = fit.nlminb
		
		# likelihood profile
		p = coef(fit)
		if (is.valid(fit)) {
			prof <- profile(fit, which='lambda')
		} else prof <- profile(fit, which='lambda', 
							   std.err=p['lambda']/100,
							   maxsteps=5000)
		# check if the profile returns a better fit
		if (class(prof) == 'profile.mle2') {
			conf = confint(prof)
			break;
		}
		# if we reach here, confint must have found a 
		# better solution. we need to repeat
		theta0 = coef(prof)
		names(theta0) = names(p)
	}
	
	k = length(theta0)
	n = length(d)
	logL = logLik(fit)
	AICc = 2 * k - 2 * logL + 2 * k * (k + 1) / (n - k - 1)
	result = c(p['lambda'], lower=conf[[1]], upper=conf[[2]], 
			   AICc=AICc, logL=logL)
	
	is.cumulative = FALSE
	mean = calc_mean(p, t, N)
	
	# Calculate cases for confidence interval:
	p.upper <- p.lower <- p
	p.lower['lambda'] <- result['lower']
	p.upper['lambda'] <- result['upper']
	mean.lower <- calc_mean(p.lower, t, N)
	mean.upper <- calc_mean(p.upper, t, N)
	
	list(result = result,
		 param = list(p=p,p.lower=p.lower,p.upper=p.upper),
		 fit = data.frame(times=t+start.time,
		 				 cases=mean,
		 				 cases.lower = mean.lower,
		 				 cases.upper = mean.upper)
	)
}

