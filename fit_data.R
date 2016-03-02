library(ggplot2)
source("exp_fit.R")

set.seed(1234)

richards.inc <- function(N,prm,times){

	K <- prm["K"]
	x0 <- prm["x0"]
	lambda <- prm["lambda"]
	alpha <- prm['alpha']

	cuminc <- 	K * N / (1 + ((K*N/x0)^alpha-1) * 
							exp(-lambda / (1 - (x0 / K / N) ^ alpha) *alpha*times)) ^(1/alpha) 

	inc <- c(x0, diff(cuminc))	
	
	return(list(inc=inc, cuminc=cuminc))
}

fit_data <- function(data, starttime, theta0, N, is.cumulative){
	
	fit <- data.frame()
	res <- exp_fit(times = data$times, 
				   X = data$cases,
				   theta0 = theta0,
				   start = starttime,
				   is.cumulative = is.cumulative,
				   N = N)
	
	lambdas <- c(r=res$result["lambda"],
						  lower=res$result["lower"],
						  upper=res$result["upper"])
	
	fit = rbind(fit, cbind(data,
						   cases.lower=data$cases,
						   cases.upper=data$cases, 
						   legend="cases"))
	fit = rbind(fit, cbind(res$fit, legend="fit"))
	
	return(list(lambdas=lambdas, fit=fit, param=res$param))
}


tt <- 1:30
data <- data.frame(times=tt, cases=(rpois(n = length(tt),lambda = exp(0.18*tt))))
# plot(data)
# text(tt,data$cases,labels = data$cases,pos = 3)
N <- max(data$cases)

starttime <- 1
is.cumulative <- T

theta0 <- c(x0=50, lambda=0.05, K=12, alpha=0.1)
myfit <- fit_data(data, starttime, theta0, is.cumulative, N = N)
df <- myfit[["fit"]]
prm <- myfit[["param"]]
prm$p

calc_mean <- function(p, times, N,is.cumulative) {
	pp <- c(as.list(p),list(times=times,N=N))
	vals <- eval(model_richards,envir=pp)
	if (!is.cumulative) diff(c(0,vals)) else vals
}

# rich.inc <- richards.inc(N = N, prm=prm$p, times=tt)
# rich.inc2 <- calc_mean(p = prm$p, times = tt, N = N,is.cumulative=F)
# plot(data$times, data$cases)
# lines(rich.inc[["inc"]])
# lines(rich.inc2,col="red")

g <- ggplot(df)+geom_pointrange(aes(x=times,y=cases,ymin=cases.lower,ymax=cases.upper,
									colour=legend, shape = legend))
plot(g)


