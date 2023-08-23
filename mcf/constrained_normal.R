rlcnorm = function(n, A, b)
{
  cols = ncol(A)
  rows = nrow(A)
  udv = svd(A, nu = rows, nv = cols)
  b = t(udv[['u']])%*%b
  W = udv[['d']]
  Y = matrix(nrow = n, ncol = cols)
  for (i in 1:rows)
  {
    Y[,i] = b[i]/W[i]
  }
  Y[,(rows+1):cols] = rnorm(n*(cols-rows))
  return(t(udv[['v']]%*%t(Y)))
}

norm_cnst = function(mu, Sigma, A, b)
{
  d = length(b)
  S = A%*%Sigma%*%t(A)
  dif = b - A%*%mu
  if (d == 1)
  {
    return(exp(-0.5*dif^2/S)) 
  } else {
    return(exp(-0.5*t(dif)%*%inv(S)%*%dif))
  }
}

norm_lcnst = function(mu, Sigma, A, b)
{
  d = length(b)
  S = A%*%Sigma%*%t(A)
  dif = b - A%*%mu
  if (d == 1)
  {
    return(-0.5*dif^2/S) 
  } else {
    return(-0.5*t(dif)%*%inv(S)%*%dif)
  }
}