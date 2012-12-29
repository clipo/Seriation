library(ggplot2)
library(randtoolbox)
library(xtable)
setwd("~/Dropbox/Research/Dissertation Project/src/Seriation/analysis/R")

# Input values
assemblages <- 64
biggest_solution <- 20
servers <- 8
cores <- 8
parallelism <- servers * cores
secs_per_trial <- 0.005   # 5 ms per test
partitions_63_nonempty <- stirling(assemblages)





# functions

permutations_to_test <- function(biggest_solution, assemblages) {
  perm <- numeric(biggest_solution)
  for(n in 3:biggest_solution) {
    perm[n] <- choose(assemblages, n) * ( factorial(n) / 2 )
  }
  perm
}

times_for_permutation_test <- function(assemblages) {
  
  
  total_trials <- (factorial(assemblages) / 2)
  
  trial_batches <- total_trials / parallelism
  
  days_in_avg_month <- 30.4368
  time_batches_sec <- trial_batches * secs_per_trial
  time_batches_days <- time_batches_sec / 86400
  time_batches_months <- time_batches_days / days_in_avg_month
  time_batches_years <- time_batches_days / 365.25
  
  times <- c(k, total_trials, time_batches_sec, time_batches_days, time_batches_months, time_batches_years)
}


times_for_solution_size <- function(k, assemblages) {


  total_trials <- choose(assemblages, k) * (factorial(k) / 2)
  

  
  trial_batches <- total_trials / parallelism
  
  days_in_avg_month <- 30.4368
  time_batches_sec <- trial_batches * secs_per_trial
  time_batches_days <- time_batches_sec / 86400
  time_batches_months <- time_batches_days / days_in_avg_month
  time_batches_years <- time_batches_days / 365.25
  
  comb <- choose(assemblages, k)
  part <- partitions_63_nonempty[k]
  
  times <- c(k, comb, total_trials, part, time_batches_sec, time_batches_days, time_batches_months, time_batches_years)
}





## actual analysis
df <- data.frame(matrix(ncol=8,nrow=0))

for(k in 3:biggest_solution) {
  row <- times_for_solution_size(k, assemblages)
  df <- rbind(df, row)
}
colnames(df) <- c("SolutionSize", "Combinations", "PossibleSolutions", "NonemptyPartitions", "Seconds", "Days", "Months", "Years")

# unscaled_title <- "Years of Processing Time for 64 Assemblages on 64x8 Cluster"
# u <- ggplot(df, aes(x = SolutionSize, y = Years)) 
# u + geom_line() + scale_x_continuous(name = "Solution Size") + scale_y_log10(name = "Years on Parallel 64x8 Cluster") + opts(title=unscaled_title) 
# ggsave(filename = "graphics/proc-time-seriation-solutions-brute-force-64-years.pdf")

y_caption <- paste("Months To Test K-Solutions on", servers, "servers with", cores, "cores", sep=" " )
unscaled_title <- paste("Processing Time and Solution Scaling for 64 Assemblages")
u <- ggplot(df, aes(x = SolutionSize, y = Months)) 
u + geom_line() + scale_x_continuous(name = "Solution Size") + scale_y_log10(name = y_caption) + opts(title=unscaled_title) 
ggsave(filename = "graphics/proc-time-seriation-solutions-brute-force-64-months.pdf")

# unscaled_title <- "Number of Combinations versus Solution Size (63 assemblages)"
# u <- ggplot(df, aes(x = SolutionSize, y = Combinations)) 
# u + geom_line() + scale_x_continuous(name = "Solution Size") + scale_y_log10(name = "Combinations") + opts(title=unscaled_title) 
# ggsave(filename = "graphics/combinations.pdf")



seriation_combinatorics_for_k <- function(assemblages,k) {
  comb <- choose(assemblages,k)
  poss_solns <- choose(assemblages, k) * ( factorial(k) / 2 )
  c(k, comb, poss_solns)
}

assem <- 64
soln_comb <- data.frame(matrix(ncol=3,nrow=0))
for(k in 3:20) {
  soln_comb <- rbind(soln_comb, seriation_combinatorics_for_k(assem, k))
}
colnames(soln_comb) <- c("K", "Combinations", "Solutions of Size K")

xt <- xtable(soln_comb, align="|c|c|r|r|", display=c("d","d","e","e"), caption="Combinations of size K and unique possible seriations of size K for N = 64 assemblages")
print(xt, include.rownames=FALSE)


########## Single Seriation Solution Permutations 
single_seriation_stats <- function(n) {
  total_trials <- (factorial(n) / 2)
  trial_batches <- total_trials / parallelism
  
  days_in_avg_month <- 30.4368
  time_batches_sec <- trial_batches * secs_per_trial
  time_batches_days <- time_batches_sec / 86400
  time_batches_months <- time_batches_days / days_in_avg_month
  time_batches_years <- time_batches_days / 365.25
  times <- c(n, total_trials, time_batches_sec,time_batches_years)
}

sscomb <- data.frame(matrix(ncol=3,nrow=0))

nseq <- c(seq(from=4, to=10, by=2), 12, 13, 14, 15, seq(from=20, to=100, by=10))

for(n in nseq) {
  sscomb <- rbind(sscomb, single_seriation_stats(n))
}

colnames(sscomb) <- c("N", "Seriation Solutions", "Seconds","Years")
capt <- paste("Number of unique seriation solutions and processing time for sets of assemblages 4 < n < 100, testing solutions across",parallelism,"cores",sep=" ")
xt2 <- xtable(sscomb, align="|c|c|r|r|r|", display=c("d","d","g","g","g"),caption=capt)
print(xt2, include.rownames=FALSE)


%%%%%%%%% analysis of subsetting %%%%%%%%%%

ex_num_sol_init <- c(3,4,6,8)
assem <- c(20, 40, 60)
max_assem <- max(assem)
subset_incr <- 5

max_subset_labels <- c(ex_num_sol_init, seq(from=10, to=max_assem/2, by=subset_incr))
max_rows <- length(max_subset_labels)

cnames <-c("Subset Size (m)")
subsets <- data.frame(matrix(ncol=0,nrow=max_rows))
# The subset sizes are the first column
subsets <- cbind(subsets, max_subset_labels)

# Create columns for each assemblage size, picking out the {n,k} for the desired subset sizes
for(n in assem) {
  # add the assemblage size as a colname
  cnames <- c(cnames, as.character(n))
  
  # calculate additional subset sizes based on the size of the assemblage, up to 1/2*N since {n,k} peaks there
  more_subsets <- seq(from=10, to=n/2, by=subset_incr)
  subsets_selected <- c(ex_num_sol_init, more_subsets)
  
  # calculate the stirling numbers of the second kind
  # the output of stirling() is a zero-based array for substantive reasons, 
  # so to pull out the appropriate values in a 1-based array index language, add one to 
  # all indices...
  stirl <- stirling(n)
  col <- stirl[subsets_selected + 1]
  
  # since this is a data frame, pad the rest of the column with NA
  num_NA_needed <- max_rows - length(col)
  col <- c(col, rep(NA_integer_, num_NA_needed))
  
  subsets <- cbind(subsets, col)
}

colnames(subsets) <- cnames

capt2 <- paste("Number of ways to form m subsets (seriation solutions) from 20, 40, and 60 assemblages")
xt3 <- xtable(subsets, align="|c|c|r|r|r|", display=c("d","d","g","g","g"),caption=capt2)
print(xt3, include.rownames=FALSE)


# graph of stirling2 for 40 assemblages to show asymmetry

assem <- 60
number_subsets <- stirling(assem)
subset_size <- seq(from=0, to=assem, by=1)
subsetgraph <- data.frame(matrix(ncol=0,nrow=(assem+1)))
subsetgraph <- cbind(subsetgraph, subset_size)
subsetgraph <- cbind(subsetgraph, number_subsets)

total_subsets <- sum(number_subsets)
print(total_subsets)
print(factorial(assem))

y_caption <- paste("Number of Unique Subsets")
u <- ggplot(subsetgraph, aes(x = subset_size, y = number_subsets)) 
u + geom_bar(stat="identity") + scale_x_continuous(name = "Number of Subsets", limits=c(0,assem+1)) + scale_y_log10(name = y_caption)  

###################################################
# analysis of total permutations for solution subsets given Stirling sums (Equation 3)

assem <- 40
numsubset <- stirling(assem)
subset_size <- seq(from=0, to=assem, by=1)

perm_per_subset_size <- factorial(assem-subset_size-1)
# remove the final NaN for subset {N,N}
perm_per_subset_size[is.nan(perm_per_subset_size)] <- 0

total_soln_per_subset_size <- numsubset * perm_per_subset_size
total_solutions <- sum(total_soln_per_subset_size)
print(total_solutions)





### total solutions for multiple seriation groups, with timing for enumerative tests

total_solutions_mult_groups <- function(n) {
  num_subsets <- stirling(n)
  subset_size <- seq(from=0, to=n, by=1)
  
  # permutations
  perm_per_subset_size <- factorial(n-subset_size-1)
  # remove the final NaN for subset {N,N}
  perm_per_subset_size[is.nan(perm_per_subset_size)] <- 0
  total_soln_per_subset_size <- num_subsets * perm_per_subset_size
  total_solutions <- sum(total_soln_per_subset_size)
  total_solutions
}

timing_mult_groups <- function(sol) {
  trial_batches <- sol / parallelism
  
  days_in_avg_month <- 30.4368
  time_batches_sec <- trial_batches * secs_per_trial
  time_batches_days <- time_batches_sec / 86400
  time_batches_months <- time_batches_days / days_in_avg_month
  time_batches_years <- time_batches_days / 365.25
  times <- c(time_batches_sec,time_batches_years)
}

multgroupsol_timing <- data.frame(matrix(ncol=3,nrow=0))

nseq <- c(seq(from=4, to=10, by=2), 12, 13, 14, 15, 16, seq(from=20, to=100, by=20))

for(n in nseq) {
  sol <- total_solutions_mult_groups(n)
  timing <- timing_mult_groups(sol)
  row <- c(n, sol, timing)
  multgroupsol_timing <- rbind(multgroupsol_timing, row)
}

colnames(multgroupsol_timing) <- c("N", "Total Solutions", "Seconds","Years")
capt <- paste("Number of total solutions with multiple seriation groups and processing time for sets of assemblages 4 < n < 100, testing solutions across",parallelism,"cores",sep=" ")
xt3 <- xtable(multgroupsol_timing, align="|c|c|r|r|r|", display=c("d","d","g","g","g"),caption=capt)
print(xt3, include.rownames=FALSE)




