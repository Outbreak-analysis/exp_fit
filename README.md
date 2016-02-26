# exp_fit
Fitting the exponential growth rate

## The files:
* exp_fit.R: defines the exp_fit function that fits the exponential growth rate
* models.R: defines some commonly used phenomenological models that are used by exp_fit
  * namely: the logistic model and the Richards model.
* test.R: an example of using exp_fit
* Ebola.csv: some Ebola data used by test.R

## The exp_fit function
`exp_fit <- function(times, X, theta0, start, is.cumulative)`

This function fits the time series (times, X), i.e., the epidemic curve, to either the Richards model or the logistic model using the mle2 function from the bbmle package. The data used for fitting is from the data indexed by start to the peak of the epidemic curve.
* See [this paper](http://link.springer.com/article/10.1007/s11538-013-9918-2)

### parameters:
* times, X: the time series
* theta0: the initial guess of the parameters
  * If it has an alpha element, then the Richards model will be used. Otherwise, the logistic model will be used. 
* start: the time as the start of the fitting window 
* is.cumulative: whether the time series is cumulaative incidence

### return values
A list with elements:
* $result: a vector containing the following named elements
  * lambda: the point estimate of the exponential growth rate
  * lower: the lower boundary of 95% confidence interval
  * upper: the upper boundary of the 95% confidence interval
  * AICc: AICc
  * logL: the log-likelihood at the point estimate
* $fit a data frame with the best fit time series
