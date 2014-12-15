library(sgeostat)
setwd("/Volumes/Macintosh HD/Users/clipo/PycharmProjects/Seriation/R")
data.path <- "../testdata/pfg-cpl.txt" 
xydata.path <- "../testdata/pfgXY.txt" 
# We read the whole file
PFG <- read.table(data.path, header = TRUE, sep="\t")
xy <- read.table(xydata.path, header = TRUE, sep="\t")

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}
GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      return(distM <- sqrt(   (g1$Northing-g2$Northing)^2 + (g1$Easting-g2$Easting)^2  )/1000)
      # if lat and long
      #  require("Imap")
      #  return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$Northing, lon.1=g1$Easting, lat.2=g2$Northing, lon.2=g2$Easting, units="m")))
    }
    return(mapply(DistM, g1, g2))

  }
  n.geopoints <- nrow(df.geopoints)
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("Assemblage", "Northing", "Easting")], 1:n.geopoints, function(x){return(list(x))})
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$Assemblage
  colnames(mat.distances) <- df.geopoints$Assemblage
  return(mat.distances)
}


pfgpercent <- PFG[-1] / rowSums(PFG[-1], na.rm = T) * 100
rownames(pfgpercent)<-PFG[,1]
colnames(pfgpercent)[1]="Assemblage"
xyPFG <- merge(x=xy, y=PFG, by="Assemblage",all=FALSE)

## distance based on percentages
distanceMatrix <- dist(pfgpercent, method = "euclidean") # distance matrix
typeFit <- hclust(distanceMatrix, method="ward.D2") 
plot(typeFit) # display dendogram

## distances based on XY location
#xyDistanceMatrix <-GeoDistanceInMetresMatrix(xyPFG)
xyTemp<-merge(y=pfgpercent,x=xy, by="row.names", all=FALSE)
xyAssemblage<-xyTemp[,2:4]
rownames(xyAssemblage)<-xyAssemblage[,1]
xyAssemblage<- subset(xyAssemblage, select = c("Easting","Northing"))
xyDistanceMatrix <- dist(xyAssemblage, method = "euclidean")
distFit <- hclust(xyDistanceMatrix, method="ward.D2") 
plot(distFit) # display dendogram

rownames(xy)<-xy[,1]
xyPFGPercent <- merge(x=xy, y=pfgpercent, by="row.names",all=FALSE)
xyPFGPercent$Easting<-xyPFGPercent$Easting/1000
xyPFGPercent$Northing<-xyPFGPercent$Northing/1000
location<-point(xyPFGPercent,x="Easting",y="Northing")
pfgPair <-pair(location, num.lags=9, maxdist=120)
pfg.v<-est.variogram(location,pfgPair,'pfg') 
library(spatstat)
#xyPFGPercent$Easting<-xyPFGPercent$Easting/10000
#xyPFGPercent$Northing<-xyPFGPercent$Northing/10000
minX<-min(xyPFGPercent$Easting)
maxX<-max(xyPFGPercent$Easting)
minY<-min(xyPFGPercent$Northing)
maxY<-max(xyPFGPercent$Northing)
pfgPattern<-ppp(xyPFGPercent$Easting,xyPFGPercent$Northing,c(minX,maxX),c(minY,maxY))
summary(pfgPattern)
plot(density(pfgPattern))

dists<-dist(xyPFGPercent[,3:4])
summary(dists)
breaks=seq(0,75,l=11)
v1<-variog(coords=xyPFGPercent[,3:4],data=xyPFGPercent[,5],breaks=breaks)
v1.summary<-cbind(c(1:10),v1$v,v1$n)
colnames(v1.summary) <- c("lag", "semi-variance", "# of pairs")
plot(v1,type="b",main="Variogram: Parkin_Punctate")

library(aqp)
library(cluster)
library(ape)
library(igraph)
library(vegan)
library(circular)
library(sharpshootR)
tFit <- as.phylo(as.hclust(typeFit))
dFit <- as.phylo(as.hclust(distFit))
dueling.dendrograms(tFit, dFit, lab.1='Types', lab.2='Geographic Distance')


typeFit <- hclust(dist(pfgpercent[order(row.names(pfgpercent)),]), "ward.D2")
distFit <- hclust(dist(xyAssemblage[order(row.names(xyAssemblage)),]), "ward.D2")
# And the plot:
layout(matrix(1:5,nrow=1),width=c(5,.75,3,.75,5))

# The first dendrogram:
l <- length(typeFit$order)
# The matrix to draw the arrows:
cbind((1:l)[order(typeFit$order)],(1:l)[order(distFit$order)]) -> ord_arrow
# The two vectors of ordered leave labels:
typeFit$labels[typeFit$order]->leaves1
distFit$labels[distFit$order]->leaves2

# The first dendrogram:
par(mar=c(3,3,3,0))
plot(as.dendrogram(typeFit),horiz=TRUE,leaflab="none", ylim=c(0,l), main="Type % Distance")

# The first serie of labels (i draw them separately because, for the second serie, I didn't find a simple way to draw them nicely on the cluster):
par(mar=c(3,0,3,0))
plot(NA, bty="n",axes=FALSE,xlim=c(0,1), ylim=c(0,l),ylab="",xlab="")
sapply(1:l,function(x)text(x=0,y=x,labels=leaves1[x], pos=4, cex=0.8))

# The arrows:
par(mar=c(3,0,3,0))
plot(NA, bty="n",axes=FALSE,xlim=c(0,1), ylim=c(0,l),ylab="",xlab="")
apply(ord_arrow,1,function(x){arrows(0,x[1],1,x[2],code=3, length=0.05, col="blue")})

# The second serie of labels:
par(mar=c(3,0,3,0))
plot(NA, bty="n",axes=FALSE, xlim=c(0,1), ylim=c(0,l), ylab="",xlab="")
sapply(1:l,function(x)text(x=1,y=x,labels=leaves2[x], pos=2, cex=0.8))

# And the second dendrogram (to reverse it I reversed the xlim vector:
par(mar=c(3,0,3,3))
plot(as.dendrogram(distFit),horiz=TRUE, xlim=c(0,max(dist(xyAssemblage))), leaflab="none", ylim=c(0,l), main="Geographic Distance")
