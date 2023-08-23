library(R6)

GammaComp = R6Class("GammaComp", list(
  alpha = 3,
  beta = 3,
  total_time = 1,
  ai = c(1:10)*0.1,
  d = 0.1,
  
  initialize = function(alpha, beta)
  {
    self$alpha = alpha
    self$beta = beta
  },
  
  generateSample = function(n = 1)
  {
    g = rgamma(n, shape = self$alpha, rate = self$beta)
    return(g)
  },
  
  density = function(x, log=TRUE)
  {
    alpha = self$alpha
    beta = self$beta
    if(log)
    {
      d = alpha*log(beta)-log(gamma(alpha)) + (alpha-1)*log(x) - beta*x
    } else
    {
      d = beta^alpha/gamma(alpha)*x^(alpha-1)*exp(-beta*x)
    }
    return(d)
  },
  
  phi = function(x)
  {
    alpha = self$alpha
    beta = self$beta
    return(0.5* (beta^2 - 2*(alpha-1)*beta/x + (alpha-1)*(alpha-2)/x^2))
  },
  
  lowerBound = function(interval)
  {
    alpha = self$alpha
    beta = self$beta
    if (alpha == 1)
    {
      return(0.5*beta^2)
    } else
    {
      return(min(self$phi(interval)))
    }
  },
  
  upperBound = function(interval)
  {
    return(max(self$phi(interval)))
  }
)
)