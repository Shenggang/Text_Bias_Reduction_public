compute_layer_acc_weight_log = function(Xs, times, bound1, bound2)
{
  diffs = diff(times)
  prob1 = 0; prob2 = 1;
  for (i in c(2:length(times)))
  {
    prob1 = prob1 + log(compute_delta(diffs[i-1], Xs[i-1], Xs[i], bound2))
    prob2 = prob2 * compute_delta(diffs[i-1], Xs[i-1], Xs[i], bound1, bound2)
  }
  return(prob1 + log(prob2*0.5+0.5))
}

compute_layer_weight = function(total_time, start, end, K, L)
{
  x = (start-end)*0.5
  y = abs(x)
  return(compute_gamma(total_time, x, -x, y+L) - compute_gamma(total_time, x, -x, y+K))
}

compute_bb_weight = function(total_time,x,y,comp,margin = 10)
{
  # for each dimension, simulate layer
  # combine and choose ub, generate ppp and then bb
  # simulate layer
  ai = comp$ai; d = comp$d
  interval = matrix(nrow=length(x), ncol=2)
  K = c(); L = c()
  # simualte layers for each bridge
  for (j in 1:length(x))
  {
    I = layer(total_time, x[j], y[j], ai, d)
    K[j] = fetch_bound(I-1, ai, d); L[j] = fetch_bound(I, ai, d)
    interval[j, ] = c(min(x[j],y[j])-L[j], max(x[j],y[j])+L[j])
  }
  # choose c
  lb = comp$lowerBound(interval)
  ub = comp$upperBound(interval) + margin
  lambda = ub - lb
  if (is.nan(lambda) || is.na(lambda) || is.infinite(lambda))
    return(-Inf)
  # simulate layered BB
  k = rpois(1, lambda*total_time)
  if (is.nan(k) || is.na(k) || is.infinite(k))
    return(-Inf)
  if (k > 3000)
    return(-Inf) # avoid processing long samples
  times = sort(runif(k, min = 0, max = total_time))
  # simulate bb individually
  bb = matrix(nrow=length(x)+1, ncol=k+2)
  bb[1,] = c(0, times, total_time)
  for (j in 1:length(x)) {
    #simulate BB
    temp = layered_BB(total_time, x[j], y[j], times, K[j], L[j])
    bb[j+1,] = extract_essential(temp, c(0,times, total_time))[2,]
  }
  # compute Poisson estimator 
  est = plestimate_full(lambda, total_time, ub, bb, times, comp)
  return(est)
}

compute_bb_weight_approx = function(total_time,x,y,comp,margin = 10, k = 50)
{
  # for each dimension, simulate layer
  # combine and choose ub, generate ppp and then bb
  # simulate layer
  ai = comp$ai; d = comp$d
  interval = matrix(nrow=length(x), ncol=2)
  K = c(); L = c()
  # simualte layers for each bridge
  for (j in 1:length(x))
  {
    I = layer(total_time, x[j], y[j], ai, d)
    K[j] = fetch_bound(I-1, ai, d); L[j] = fetch_bound(I, ai, d)
    interval[j, ] = c(min(x[j],y[j])-L[j], max(x[j],y[j])+L[j])
  }
  # choose c
  lb = comp$lowerBound(interval)
  ub = comp$upperBound(interval) + margin
  lambda = ub - lb
  # simulate layered BB
  times = (1:k)/(k+1)*total_time
  # simulate bb individually
  bb = matrix(nrow=length(x)+1, ncol=k+2)
  bb[1,] = c(0, times, total_time)
  for (j in 1:length(x)) {
    #simulate BB
    temp = layered_BB(total_time, x[j], y[j], times, K[j], L[j])
    bb[j+1,] = extract_essential(temp, c(0,times, total_time))[2,]
  }
  # compute Poisson estimator 
  est = plestimate_full(lambda, total_time, ub, bb, times, comp)
  return(est)
}

compute_gamma = function(total_time, start, end, bound)
{
  precision = 10
  total = max(ceiling((sqrt(total_time*(precision*log(10)+(start^2+end^2)/4)/2)+(bound+(start+end)/2))/2/bound),1)
  prob = 1
  for(j in 1:total)
  {
    prob = prob - (sigma_term(total_time, start, end, bound, j) - tau_term(total_time, start, end, bound, j))
  }
  return(prob)
}

compute_delta = function(total_time, start, end, bound1, bound2 = FALSE)
{
  if (max(abs(start),abs(end))> bound1)
  {
    return(FALSE)
  }
  if (end == 0)
  {
    end = start; start = 0
  }
  if (start == 0)
  {
    if (bound2 == FALSE)
    {
      return(compute_delta_zero_1(total_time, end, bound1))
    } else
    {
      return(compute_delta_zero_2(total_time, end, bound1, bound2))
    }
  } else
  {
    if (bound2 == FALSE)
    {
      return(compute_delta_1(total_time, start, end, bound1))
    } else
    {
      return(compute_delta_2(total_time, start, end, bound1, bound2))
    }
  }
}

compute_delta_zero_1 = function(total_time, end, bound1)
{
  precision = 10
  total = max(ceiling((sqrt(total_time*precision*log(10)+end^2/4)+end/2)/bound1),1)
  prob = 1
  for (count in 1:total)
  {
    prob = prob - ksi_term(total_time, -end, bound1, count)/end + ksi_term(total_time, end, bound1, count)/end
  }
  return(prob)
}


compute_delta_zero_2 = function(total_time, end, bound1, bound2)
{
  precision = 10
  total = max(ceiling((sqrt(total_time*precision*log(10)+end^2/4)+end/2)/bound1),1)
  nom = end
  den = end
  for (count in 1:total)
  {
    nom = nom - ksi_term(total_time, -end, bound1, count) + ksi_term(total_time, end, bound1, count)
    den = den - ksi_term(total_time, -end, bound2, count) + ksi_term(total_time, end, bound2, count)
  }
  return(nom/den)
}

compute_delta_1 = function(total_time, start, end, bound1)
{
  precision = 10
  Z = 1/(1-exp(-2*start*end/total_time))
  limit = bound1*0.5
  start = start - limit; end = end - limit
  total = max(ceiling((sqrt(total_time*(precision*log(10)+(start^2+end^2)/4)/2)+(limit+(start+end)/2))/2/limit),1)
  prob = Z
  for(j in 1:total)
  {
    prob = prob - (sigma_term(total_time, start, end, limit,j) - tau_term(total_time, start, end, limit,j))*Z
  }
  return(prob)
}

compute_delta_2 = function(total_time, start, end, bound1, bound2)
{
  precision = 10
  limit1 = bound1*0.5
  limit2 = bound2*0.5
  u1_K = start - limit1; u2_K = end - limit1
  u1_L = start - limit2; u2_L = end - limit2
  total = max(ceiling((sqrt(total_time*(precision*log(10)+(u1_K^2+u2_K^2)/4)/2)+(limit1+(u1_K+u2_K)/2))/2/limit1),1)
  nom = 1
  den = 1
  for (j in 1:total)
  {
    nom = nom - sigma_term(total_time, u1_K, u2_K, limit1, j) + tau_term(total_time, u1_K, u2_K, limit1, j)
    den = den - sigma_term(total_time, u1_L, u2_L, limit2, j) + tau_term(total_time, u1_L, u2_L, limit2, j)
  }
  return(nom/den)
}
  