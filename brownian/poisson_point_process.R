poisson_point_process = function(total_time, M)
{
  N = rpois(1, total_time*M);
  if (N > 3000 || is.na(N))
  {
    ppp = NaN
  } else{
    ppp = c()
    if (N > 0)
    {
      ppp = rbind(runif(N,0,total_time), runif(N,0,M));
    }
    if (length(ppp)>2)
    {
      ppp = ppp[, order(ppp[1,])]
    }
  }
  return(ppp);
}