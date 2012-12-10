library(ggplot2)
library(extrafont)

### Already have read in fonts (see previous answer on how to do this)
loadfonts()
setwd("~/Dropbox/Research/Dissertation Project/src/Seriation/R")

# Input values
assemblages <- 63
biggest_solution <- 63
df <- data.frame(rbind(c(0, 0, 0, 0, 0, 0, 0, 0)))
colnames(df) <- c("SolutionSize", "Combinations", "Permutations", "NonemptyPartitions", "Seconds", "Days", "Months", "Years")


partitions_63_nonempty <- c(1., 4.61169*10^18, 1.9076*10^29, 3.54461*10^36, 9.03498*10^41, 
                            1.46612*10^46, 3.4559*10^49, 1.94257*10^52, 3.59062*10^54, 
                            2.71973*10^56, 9.87897*10^57, 1.93258*10^59, 2.22363*10^60, 
                            1.61137*10^61, 7.76382*10^61, 2.59817*10^62, 6.25869*10^62, 
                            1.11783*10^63, 1.51746*10^63, 1.59892*10^63, 1.33135*10^63, 
                            8.89653*10^62, 4.83527*10^62, 2.16251*10^62, 8.04026*10^61, 
                            2.50763*10^61, 6.61279*10^60, 1.48488*10^60, 2.85687*10^59, 
                            4.73578*10^58, 6.79727*10^57, 8.48437*10^56, 9.2455*10^55, 
                            8.82582*10^54, 7.40286*10^53, 5.4702*10^52, 3.569*10^51, 
                            2.05998*10^50, 1.05351*10^49, 4.77992*10^47, 1.92583*10^46, 
                            6.89448*10^44, 2.19378*10^43, 6.20403*10^41, 1.5587*10^40, 
                            3.47628*10^38, 6.87377*10^36, 1.20301*10^35, 1.85934*10^33, 
                            2.53061*10^31, 3.02205*10^29, 3.15235*10^27, 2.85621*10^25, 
                            2.2322*10^23, 1.49155*10^21, 8.42613*10^18, 3.96605*10^16, 
                            1.52526*10^14, 4.66453*10^11, 1.09007*10^9, 1.82671*10^6, 1953., 1.)





# functions

permutations_to_test <- function(biggest_solution, assemblages) {
  perm <- numeric(biggest_solution)
  for(n in 3:biggest_solution) {
    perm[n] <- choose(assemblages, n) * ( factorial(n) / 2 )
  }
  perm
}




times_for_solution_size <- function(biggest_solution, assemblages, df) {
  secs_per_trial <- 0.01
  servers <- 16
  cores <- 8
  perm <- permutations_to_test(biggest_solution, assemblages)

  
  parallelism <- servers * cores
  
  total_trials <- sum(perm)
  trial_batches <- total_trials / parallelism
  
  days_in_avg_month <- 30.4368
  time_batches_sec <- trial_batches * secs_per_trial
  time_batches_days <- time_batches_sec / 86400
  time_batches_months <- time_batches_days / days_in_avg_month
  time_batches_years <- time_batches_days / 365.25
  
  comb <- choose(assemblages, biggest_solution)
  part <- partitions_63_nonempty[biggest_solution]
  
  times <- c(biggest_solution, comb, total_trials, part, time_batches_sec, time_batches_days, time_batches_months, time_batches_years)
  df <- rbind(df, times)
  
  df
}


theme_xkcd <- theme(
  panel.background = element_rect(fill="white"),
  axis.ticks = element_line(colour=NA),
  panel.grid = element_line(colour="white"),
  axis.text.y = element_text(colour=NA),
  axis.text.x = element_text(colour="black"),
  text = element_text(size=16, family="Humor Sans")
)



## actual analysis

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

unscaled_title <- "Number of Combinations versus Solution Size (63 assemblages)"
u <- ggplot(df, aes(x = SolutionSize, y = Combinations)) 
u + geom_line() + scale_x_continuous(name = "Solution Size") + scale_y_log10(name = "Combinations") + opts(title=unscaled_title) 
ggsave(filename = "graphics/combinations.pdf")


num <- numeric(29)
for(k in 3:32) {
  num[k] <- format(choose(64, k), scientific = TRUE)
}

