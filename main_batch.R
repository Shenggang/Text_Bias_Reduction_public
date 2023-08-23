source('./fit_data.R')


# constrained sampling
total_samples = 2000
total_entries = 5
total_time = 0.1

A = matrix(beta, nrow = 1)
S = -intercept

results = list()
result.mean = matrix(0, nrow=total_entries, ncol=4)
result.ESS = c()

cores = 5
cl = makeCluster(cores, outfile="log")
registerDoParallel(cl)

log_file = 'log'

cat("\n\n\n===============================================\n", file=log_file, append=TRUE)
start_time = Sys.time()
summary <- foreach (i = 1001:(1000+total_entries), .packages = c("R6","statmod","mvtnorm")) %dopar%
{
  #cat("======== Computing entries ", i," / ", total_entries," ============\n")
  cat("======== Computing entries ", i," / ", total_entries," ============\n", file=log_file, append=TRUE)
  cat(paste("Starting at ", Sys.time(), sep = "-----"), '\n', file=log_file, append=TRUE)
  comps = list()
  word_counts = as.double(dataframe[i, 2:5])
  word_counts = word_counts + (word_counts==0)
  shapes = word_counts*rates
  for (j in 1:4)
  {
    comps[[j]] = GammaComp$new(shapes[j], rates[j])
  }
  results = pf_cons(total_samples, comps, A, S, total_time, positive = TRUE)
  before = beta%*%t(dataframe[i, 2:5])+intercept
  after = beta%*%round(results$mean)+intercept
  # cat("Bias orig = ", dataframe$score[i], "   Bias before = ", before, 
  #     ",   Bias after = ", after, ",   Improv = ", abs(before) - abs(after),";\n", collapse='')
  cat("Bias orig = ", dataframe$score[i], "   Bias before = ", before, 
      ",   Bias after = ", after, ",   Improv = ", abs(before) - abs(after),";\n", collapse='', file=log_file, append=TRUE)
  return(list(r.mean = results$mean, r.ESS = results$ESS))
}
stopCluster(cl)
cat("Run complete, time elapsed = ", round(Sys.time()-start_time), " seconds.\n\n", file=log_file, append=TRUE)
save(summary, file='results_10000.rdata') 
