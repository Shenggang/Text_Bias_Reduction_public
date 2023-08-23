pestimate = function(l,t,c,x,y,comp)
{
  # decide number of time points
  k = rpois(1, l*t)
  
  # decide where
  times = runif(k, min = 0, max = t)
  
  # sample a BB
  bb = brownian_bridge(t, x, y, times)
  
  # calculate the estimate
  phi = c - comp$phi(bb[2,2:(k+1)])
  return(exp((l-c)*t)*l^(-k)*prod(phi))
}

plestimate = function(l,t,c,x,y,comp)
{
  # decide number of time points
  k = rpois(1, l*t)
  
  # decide where
  times = runif(k, min = 0, max = t)
  
  # sample a BB
  bb = brownian_bridge(t, x, y, times)
  
  # calculate the estimate
  phi = c - comp$phi(bb[2,2:(k+1)])
  return((l-c)*t-k*log(l)+sum(log(phi)))
}

plestimate_full = function(l,t,ub,bb,tp = c(),comp)
{
  if (length(tp) == 0)
  {
    tp = bb[1,]
  }
  tp = sort(tp)
  # calculate the estimate
  k = length(tp)
  phi= -comp$phi(bb[-1, 2:(length(tp)-1)])+ub
  i = 1; j = 1
  s = sum(log(phi))
  return((l-ub)*t-k*log(l)+s)
}

p_parameter = function(t,x,y,comp)
{
  # choose l,c
  k = rpois(1, 50)
  times = runif(k, min = 0, max = t)
  bb = brownian_bridge(t, x, y, times)
  phi = comp$phi(bb[2,])
  c = abs(max(phi))
  l = c - comp$lowerBound()
  
  # decide number of time points
  k = rpois(1, l*t)
  
  # decide where
  times = runif(k, min = 0, max = t)
  
  # sample a BB
  bb = brownian_bridge(t, x, y, times)
  
  # calculate the estimate
  phi = comp$phi(bb[2,2:(k+1)])
  c = abs(max(phi))*1.1
  return(list(bb = bb, c = c, l = l))
}
