library("ggplot2")
setwd("~/Dropbox/Research/Dissertation Project/analysis/seriationtrials")

# Input values
assemblages <- 64
biggest_solution <- 64
df <- data.frame(rbind(c(0, 0, 0, 0, 0, 0)))
colnames(df) <- c("SolutionSize", "Permutations", "Seconds", "Days", "Months", "Years")

# functions

permutations_to_test <- function(biggest_solution, assemblages) {
  perm <- numeric(biggest_solution)
  for(n in 3:biggest_solution) {
    perm[n]<- choose(assemblages, n) * ( factorial(n) / 2 )
  }
  perm
}




times_for_solution_size <- function(biggest_solution, assemblages, df) {
  secs_per_trial <- 0.01
  servers <- 16
  cores <- 8
  comb <- permutations_to_test(biggest_solution, assemblages)
  
  
  
  
  parallelism <- servers * cores
  
  
  total_trials <- sum(comb)
  trial_batches <- total_trials / parallelism
  
  days_in_avg_month <- 30.4368
  time_batches_sec <- trial_batches * secs_per_trial
  time_batches_days <- time_batches_sec / 86400
  time_batches_months <- time_batches_days / days_in_avg_month
  time_batches_years <- time_batches_days / 365.25
  
  times <- c(biggest_solution, total_trials, time_batches_sec, time_batches_days, time_batches_months, time_batches_years)
  df <- rbind(df, times)
  df
}


permutations <- permutations_to_test(32, 64)



for(size in 3:biggest_solution) {
  df <- times_for_solution_size(size, assemblages, df)
}

unscaled_title <- "Years of Processing Time for 64 Assemblages on 64x8 Cluster"
u <- ggplot(df, aes(x = SolutionSize, y = Years)) 
u + geom_line() + scale_x_continuous(name = "Solution Size") + scale_y_log10(name = "Years on Parallel 64x8 Cluster") + opts(title=unscaled_title) 
ggsave(filename = "graphics/proc-time-seriation-solutions-brute-force-64-years.pdf")

unscaled_title <- "Months of Processing Time for 64 Assemblages on 64x8 Cluster"
u <- ggplot(df, aes(x = SolutionSize, y = Months)) 
u + geom_line() + scale_x_continuous(name = "Solution Size") + scale_y_log10(name = "Months on Parallel 64x8 Cluster") + opts(title=unscaled_title) 
ggsave(filename = "graphics/proc-time-seriation-solutions-brute-force-64-months.pdf")


