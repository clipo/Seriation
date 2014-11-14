library("devtools")
install_github("gianmarcoalberti/CAseriation")
library(CAseriation)

# read data from choosen table
mydata <- read.table(file="./pfg-cpl-headers.txt", row.names=1, header=T)
pfg <- data.frame(mydata)
check.ca.plot(pfg,1,2)
sort.table(pfg,1)
plot.clusters.rows(pfg,1,2) 
plot.clusters.cols(pfg,1,2) 