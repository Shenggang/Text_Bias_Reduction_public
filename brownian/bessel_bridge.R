bessel_bridge = function(total_time, end, time_points)
{
  if (length(time_points) == 0)
  {
    return(rbind(c(0,total_time),c(0, end)))
  }
  time_points = sort(time_points)
  if (tail(time_points, n=1) >= total_time)
  {
    print("Required time point exceeds total time. Unable to simulate Bessel bridge")
    return()
  }
  bb1 = brownian_bridge(total_time, 0,0, time_points)
  bb2 = brownian_bridge(total_time, 0,0, time_points)
  bb3 = brownian_bridge(total_time, 0,0, time_points)
  t = c(0, time_points, total_time)/total_time
  beb = sqrt((end*t+bb1[2,])^2+bb2[2,]^2+bb3[2,]^2)
  return(rbind(c(0,time_points, total_time),beb))
}