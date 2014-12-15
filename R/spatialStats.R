library(plyr)
setwd("/Volumes/Macintosh HD/Users/clipo/PycharmProjects/Seriation/R")
data.path <- "../testdata/pfg.txt" 
xydata.path <- "../testdata/xy.txt" 
# We read the whole file
PFG <- read.table(data.path, header = FALSE, sep="\t")
xy <- read.table(xydata.path, header = TRUE, sep="\t")

#xy <- as.numeric(xy[-1])
#pfgpercent <- list()

pfgpercent <- PFG[-1] / rowSums(PFG[-1], na.rm = T) * 100
coords <- coordinates(xy[-1])
Sy8_nb <- knn2nb(knearneigh(coords, k = 1), row.names = xy[,1])
Sy9_nb <- knn2nb(knearneigh(coords, k = 2), row.names = xy[,1])
Sy10_nb <- knn2nb(knearneigh(coords, k = 4), row.names = xy[,1])
nb_l <- list(k1 = Sy8_nb, k2 = Sy9_nb, k4 = Sy10_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose = FALSE, force = TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)

dsts <- unlist(nbdists(Sy8_nb, coords))
summary(dsts)
max_1nn <- max(dsts)
max_1nn

Sy11_nb <- dnearneigh(coords, d1 = 0, d2 = 0.75 * max_1nn, row.names = xy[,1])
Sy12_nb <- dnearneigh(coords, d1 = 0, d2 = 1 * max_1nn, row.names = xy[,1])
Sy13_nb <- dnearneigh(coords, d1 = 0, d2 = 1.5 * max_1nn, row.names = xy[,1])
nb_l <- list(d1 = Sy11_nb, d2 = Sy12_nb, d3 = Sy13_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose = FALSE, force = TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)