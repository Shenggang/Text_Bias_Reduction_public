layer = function(total_time, start, end, ai, d)
{
  x = (start - end)*0.5
  y = abs(x)
  u = runif(1)
  tf = check_gamma(u,total_time, x, -x, y+ai[1]);
  if (tf)
  {
    return(1);
  }
  else
  {
    min = 1; max = 2; l = length(ai)
    K = fetch_bound(max, ai, d);
    while (!check_gamma(u,total_time, x, -x, y+K))
    {
      min = max;
      max = max*2;
      K = fetch_bound(max, ai, d)
    }
    while (min < max -1)
    {
      mid = floor((min+max)/2);
      K = fetch_bound(mid, ai, d)
      if (check_gamma(u,total_time, x, -x, y+K))
      {
        max = mid;
      } else
      {
        min = mid;
      }
    }
    return(max)
  }
}

layered_BB = function(total_time, start, end, time_points, K, L)
{
  if (length(time_points) > 0)
  {
    time_points = sort(time_points)
  }
  K_old = K; L_old = L
  while(TRUE)
  {
    K = K_old; L = L_old
    rejected = FALSE
    xb = min(start, end); yb = max(start,end)
    if (runif(1) <= 0.5)
    {
      interval = c(xb - L, xb - K)
      mid = r_min_brownian_bridge(total_time, start, end, interval)
      bb = brownian_bridge_with_min(total_time, start, end, mid[1], mid[2], time_points)
      # set the starts, ends and intervals
      Xs = bb[2,] - mid[1]
      L = yb + L - mid[1]; K = yb + K - mid[1]
    } else
    {
      interval = c(yb + K, yb + L)
      mid = r_max_brownian_bridge(total_time, start, end, interval)
      bb = brownian_bridge_with_max(total_time, start, end, mid[1], mid[2], time_points)
      # set the starts, ends and intervals
      Xs = mid[1] - bb[2,]
      L = mid[1] - xb + L; K = mid[1] - xb + K
    }
    #check if accepted
    if ((max(bb[2,]) > yb+L_old)||(min(bb[2,]) < xb-L_old))
    {
      rejected=TRUE
    }
    th = c(sum(time_points < mid[2]), sum(time_points <= mid[2])+1)
    tps = c(0, time_points[consc(1,th[1])], mid[2], time_points[consc(th[2],length(time_points))], total_time)
    diffs = diff(tps)
    for (i in c(2:length(tps)))
    {
      if (!check_delta(diffs[i-1], Xs[i-1], Xs[i], L))
      {
        # rejected
        rejected = TRUE
        break
      }
    }
    if (!rejected)
    {
      for (i in c(2:length(tps)))
      {
        if (!check_delta(diffs[i-1], Xs[i-1], Xs[i], K, L))
        {
          if (runif(1) <= 0.5 && !rejected)
          {
            return(bb)
          } else
          {
            #rejected
            rejected = TRUE
            break
          }
        }
      }
      if (!rejected)
        return(bb)
    }
  }
}

generate_LBB = function(total_time, start, end, time_points, K, L)
{
  xb = min(start, end); yb = max(start,end)
  if (runif(1) <= 0.5)
  {
    interval = c(xb - L, xb - K)
    mid = r_min_brownian_bridge(total_time, start, end, interval)
    bb = brownian_bridge_with_min(total_time, start, end, mid[1], mid[2], time_points)
    # set the starts, ends and intervals
    Xs = bb[2,] - mid[1]
    L = yb + L - mid[1]; K = yb + K - mid[1]
  } else
  {
    interval = c(yb + K, yb + L)
    mid = r_max_brownian_bridge(total_time, start, end, interval)
    bb = brownian_bridge_with_max(total_time, start, end, mid[1], mid[2], time_points)
    # set the starts, ends and intervals
    Xs = mid[1] - bb[2,]
    L = mid[1] - xb + L; K = mid[1] - xb + K
  }
  return(list(Xs=Xs, mid=mid, bb=bb, K=K, L=L))
}

layered_multi_BB = function(total_time, start, end, time_points, K, L)
{
  if (length(time_points) > 0)
  {
    time_points = sort(time_points)
  }
  C = length(start)
  while(TRUE)
  {
    rejected = c(FALSE, FALSE)
    bb = matrix(nrow=length(start)+1, ncol=length(time_points)+2)
    bb[1,] = c(0, time_points, total_time)
    for (i in 1:C)
    {
      res = generate_LBB(total_time, start[i], end[i], time_points, K[i], L[i])
      mid = res$mid; Xs = res$Xs; bL = res$L; bK = res$K; temp = res$bb
      bb[i+1,] = extract_essential(temp, bb[1,])[2,]
      #check if accepted
      th = c(sum(time_points < mid[2]), sum(time_points <= mid[2])+1)
      tps = c(0, time_points[consc(1,th[1])], mid[2], time_points[consc(th[2],length(time_points))], total_time)
      diffs = diff(tps)
      for (j in c(2:length(tps)))
      {
        if (!check_delta(diffs[j-1], Xs[j-1], Xs[j], bL))
        {
          # rejected
          rejected[i] = TRUE
          break
        }
      }
      if (!rejected[i])
      {
        for (j in c(2:length(tps)))
        {
          if (!check_delta(diffs[j-1], Xs[j-1], Xs[j], bK, bL))
          {
            if (runif(1) <= 0.5 && !rejected[i])
            {
              break
            } else
            {
              #rejected
              rejected[i] = TRUE
              break
            }
          }
        }
      }
    }
    if (sum(rejected)==0)
      return(bb)
  }
}

extract_essential = function(bb, time_points)
{
  result = matrix(nrow=2, ncol=length(time_points))
  result[1,] = time_points
  i = 1
  for (j in 1:length(bb[2,]))
  {
    if (bb[1,j] == time_points[i])
    {
      result[2,i] = bb[2,j]
      i = i+1
    }
    if (i == length(time_points)+1)
      break
  }
  return(result)
}





