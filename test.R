source("exp_fit.R")

capitalize <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), substring(s, 2),
          sep="", collapse=" ")
}

pick.country <- function(ebola, country) {
    country = paste(unlist(strsplit(country, " ")), collapse="")
    name = capitalize(paste("Cases", country, sep="_"))
    if (!(name %in% colnames(ebola))) {
        NA
    } else {
        data <- data.frame(times = ebola$Day, cases = ebola[,name])
        ind = rev(which(!is.na(data$cases)))
        data[ind,]
    }
}

fit.all <- function(countries) {
  ebola = read.csv("ebola.csv")
  lambdas = data.frame()
  fit = data.frame()
  country.names <- names(countries)
  for (name in country.names) {
    Name = capitalize(name)
    data = pick.country(ebola, Name)
    res=exp_fit(data$times, data$cases,
                theta0=countries[[name]]$theta0,
                countries[[name]]$start,
                is.cumulative=TRUE 
    )
    lambdas = rbind(lambdas, data.frame(r=res$result["lambda"],
                                        lower=res$result["lower"],
                                        upper=res$result["upper"]))
    fit = rbind(fit, cbind(data, legend="cases", country=Name))
    fit = rbind(fit, cbind(res$fit, legend="fit", country=Name))
  }
  list(lambdas=lambdas, fit=fit)
}

add.country <- function(countries=list(), name, start=1, theta0=c(x0=1, lambda=0.1, K=10, alpha=0.5)) {
    countries[[name]] = list(start=start, theta0=theta0)
    countries
}

countries <- add.country(name="Guinea", theta0 = c(x0=50, lambda=0.05, K=12, alpha=0.1))
countries <- add.country(countries, "Sierra Leone")
countries <- add.country(countries, "Liberia")

fit <- fit.all(countries)

print(fit$lambdas)

require(ggplot2)

# The palette with black:
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fig <- ggplot(data=fit$fit, aes(x=times, y=cases, shape=legend, col=legend)) +
    geom_point() + scale_y_log10() + facet_grid(country ~ .) +theme_bw() +
    scale_fill_manual(values=cbPalette) +
    scale_colour_manual(values=cbPalette)


print(fig)

