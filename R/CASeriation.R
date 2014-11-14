library("devtools")
install_github("gianmarcoalberti/CAseriation")
library(CAseriation)

setwd("/Volumes/Macintosh HD/Users/clipo/PycharmProjects/Seriation/R")

plot.clusters.rows <- function (data, x,y){
  res.ca <- CA(data, axes=c(x, y), graph=FALSE)
  resclust.rows<-HCPC(res.ca, nb.clust=-1, metric="euclidean", method="ward.D2", order=TRUE, graph.scale="inertia", graph=FALSE, cluster.CA="rows")
  plot(resclust.rows, axes=c(x,y), choice="map", draw.tree=FALSE, ind.names=TRUE, new.plot=TRUE)
  write.xlsx(resclust.rows$data.clust, "CAseriation_rows_clust.xlsx", sheetName="rowClusters", col.names=TRUE, row.names=TRUE, append=FALSE)
}

sort.table <- function (data, dim){
  # get some details about the input table in order to calculate the table's dimensionality
  nrows <- nrow(data)
  ncols <- ncol(data)
  numb.dim.cols<-ncol(data)-1
  numb.dim.rows<-nrow(data)-1
  a <- min(numb.dim.cols, numb.dim.rows) #dimensionality of the table
  #get the CA dataframe after the 'ca' package, selecting a number of dimensions equal to the table's dimensionality
  res.ca<-ca(data, nd=a)
  #get the coordinates on the selected CA axis
  row.c<-res.ca$rowcoord[,dim]
  col.c<-res.ca$colcoord[,dim]
  #seriation
  #sort the table according to the coord of the selected CA dimension and plot seriation chart (Bertin plot)
  print(sorted.table<-data[order(row.c), order(col.c)]) #sort the table
  print(sorted.table.PA<-apply(sorted.table, 2, function(x) ifelse(x>0, 1, 0))) #transform the sorted table into incidence table (presence/absence) to be used for Bertin plot
  bertinplot(as.matrix(sorted.table.PA), options = list(panel=panel.squares, spacing = 0)) #plot the sorted incidence table
  bertinplot(as.matrix(sorted.table.PA), options = list(panel=panel.squares, spacing = 0, reverse=TRUE)) #plot the sorted incidence table reversed
  #plot the battleship chart for the seriated table
  battleship.plot(sorted.table)
  #plot the battleship chart for the trasposed seriated table
  battleship.plot(t(sorted.table))
  #export relevant data to Excel
  write.xlsx(data, "CAseriation.xlsx", sheetName="originalDATA", col.names=TRUE, row.names=TRUE, append=FALSE)
  write.xlsx(sorted.table, "CAseriation.xlsx", sheetName="sorted_table", col.names=TRUE, row.names=TRUE, append=TRUE)
  write.xlsx(sorted.table.PA, "CAseriation.xlsx", sheetName="sorted_table_PA", col.names=TRUE, row.names=TRUE, append=TRUE)
}

# read data from choosen table
mydata <- read.table(file="./pfg-cpl-headers.txt", row.names=1, header=T)
pfg <- data.frame(mydata)
check.ca.plot(pfg,1,2)
sort.table(pfg,1)
plot.clusters.rows(pfg,1,2) 
plot.clusters.cols(pfg,1,2) 