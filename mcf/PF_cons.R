MHFC_gen_sample_layered = function(comps, A, S, total_time, positive = FALSE, margin = 10)
{
  C = length(comps);
  B = A%*%(sqrt(total_time)*diag(C))
  x = c()
  # generate initial points
  sampled = FALSE
  while (!sampled)
  {
    for (i in 1:C)
    {
      x[i] = comps[[i]]$generateSample()
    }
    # generate end points
    rhs = S-A%*%x
    y = (sqrt(total_time)*diag(C))%*%t(rlcnorm(1, B, rhs)) + x
    sampled = !positive
    if (positive)
    {
      sampled = all(y>0)
    }
    ap = norm_lcnst(as.matrix(x), diag(C)*total_time, A, S)
    if (is.na(ap) || is.nan(ap))
    {
      ap = -Inf
    }
    if (sampled && log(runif(1)) >= ap)
    {
      sampled = FALSE
    }
  }
  Z = 1
  for (j in c(1:C))
  {
    if (is.nan(Z)) {break;}
    Z = Z + compute_bb_weight(total_time, x[j], y[j], comps[[j]], margin)
    #print(Z)
  }
  
  return(list(x = x, y = y, est = Z))
}

pf_cons = function(total_samples, comps, A, S, total_time, margin = 10, positive=FALSE, max_trial=100)
{
  C = length(comps)
  result = matrix(nrow=total_samples, ncol=C)
  weights = matrix(nrow=total_samples, ncol=1)
  ptm = proc.time()
  print(paste("Starting at ", Sys.time(), sep = "-----"))
  # produce samples
  for (i in 1:total_samples)
  {
    #print(i)
    if (i %% 2000 == 0)
      print(paste(Sys.time(),i,sep = "-----"))
    
    trials = 0
    rst = MHFC_gen_sample_layered(comps, A, S, total_time, positive = positive, margin = 10)
    while ((is.infinite(rst$est) || is.nan(rst$est)) && (trials < max_trial))
    {
      rst = MHFC_gen_sample_layered(comps, A, S, total_time, positive = positive, margin = 10)
      trials = trials+1
    }
    result[i,] = rst$y
    weights[i] = rst$est
    #print(weights[i])
  }
  # normalize weights
  norm_weights = weights - max(weights)
  nW = exp(norm_weights)/sum(exp(norm_weights))
  # Compute ESS
  ESS = sum(nW)^2/sum(nW^2)
  
  # compute mean and variance
  mean_est = t(result)%*%nW
  var_est = t(result^2)%*%nW - mean_est
  results = list(result=result, weights=nW, time=proc.time()-ptm, mean=mean_est, var=var_est, ESS=ESS)
  return(results)
}