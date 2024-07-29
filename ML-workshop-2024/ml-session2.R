########################################################################
# Machine Learning Methods in R: Boot Camp
# July 30, 2024
########################################################################

# Do this only once:
# install.packages(c("caTools", "dendextend", "ggbeeswarm", "ggplot2", "pheatmap", "pROC", "randomForest", "readr", "rpart", "rpart.plot", "Rtsne"))

library(caTools)
library(dendextend)
library(ggbeeswarm)
library(ggplot2)
library(pheatmap)
library(pROC)
library(randomForest)
library(readr)
library(rpart)
library(rpart.plot)
library(Rtsne)


########################################################################

# Set your own path - different from mine

setwd('~/Desktop/bioinfo-tutorials/full-course-2024-ML-workshop/')

library(readr)
clu <- read_table2("r-sample-data.txt")
mm = data.matrix(clu[,-(1:2)])
rownames(mm) = clu$Gene

###########################################################################

tt = t(mm)

kmeans(tt, centers = 2)
km = kmeans(tt, centers = 2)
str(km)

plot(tt, col = km$cluster, pch = km$cluster)
plot(tt[,2:3], col = km$cluster, pch = km$cluster)

pairs(tt[,1:4], col = km$cluster, pch = km$cluster)


km = kmeans(tt, centers = 2)
pairs(tt[,1:4], col = km$cluster, pch = km$cluster)

km = kmeans(tt, centers = 3)
pairs(tt[,1:4], col = km$cluster, pch = km$cluster)

# k-means algorithm with 3 centers, run 10 times
km = kmeans(tt, centers = 3, nstart = 10)
pairs(tt[,1:4], col = km$cluster, pch = km$cluster)

#####

km = kmeans(tt, centers = 2, nstart = 1)
km$tot.withinss
sum(km$withinss)
km$betweenss


km = kmeans(tt, centers = 1)
km$tot.withinss

km = kmeans(tt, centers = 2)
km$tot.withinss

km = kmeans(tt, centers = 3)
km$tot.withinss

km = kmeans(tt, centers = 4)
km$tot.withinss


plot( sapply(1:10, function(k) kmeans(tt, centers = k, nstart = 10)$tot.withinss))
barplot( sapply(1:10, function(k) kmeans(tt, centers = k, nstart = 10)$tot.withinss))


#####

tt0 = tt
tt0[1:20,1:20] = 0

pheatmap(tt0)

# EXERCISE: plot tot.withinss for k in 1:10, use apply

# scree plot

plot( sapply(1:10, function(k) kmeans(tt0, centers = k, nstart = 10)$tot.withinss))
barplot( sapply(1:10, function(k) kmeans(tt0, centers = k, nstart = 10)$tot.withinss))

########################################################################
# Principal components 

km = kmeans(tt, centers = 2)

pca = prcomp(tt)
plot(pca)
summary(pca)
plot(pca$x, col = km$cluster, pch = km$cluster)

pca = prcomp(tt0)
plot(pca)
summary(pca)
plot(pca$x, col = km$cluster, pch = km$cluster)


########################################################################

plot( hclust( dist(tt) ) )
plot( as.dendrogram( hclust( dist(tt) ) ) )


#install.packages('dendextend')
library(dendextend)

km = kmeans(tt, centers = 2)

hc = hclust(dist(tt, method = 'euclidean'), method = 'complete')
dend = as.dendrogram(hc)
labels_colors(dend) =  km$cluster [ order.dendrogram(dend) ]
plot(dend)



########################################################################
# Different dendrograms

dist_method  = "manhattan" # "euclidean", "maximum", "manhattan"
clust_method = "ward.D2" # "complete", "average", "single", "ward.D", "ward.D2"

# Compute Euclidean distance between samples
dist=dist(tt , diag=TRUE, method = dist_method)

# Perform clustering with hclust

hc <- hclust(dist, method = clust_method, )

summary(hc)
plot(hc)

cutree(hc, h = 60)
cutree(hc, h = 3)

table(km$cluster, cutree(hc, h = 15))


########################################################################
# Pretty dendrogram: modified from https://r-graph-gallery.com/31-custom-colors-in-dendrogram.html

# To color each leaf of the Tree, change the attribute of each leaf. 
# This can be done using the dendrapply function. 

dhc <- as.dendrogram(hc)
plot(dhc)

colLab = function(n){
  if(is.leaf(n)){
    #take the current attributes
    a=attributes(n)
    
    #extract the first letter of the column label 
    specie=substr(attributes(n)$label, 1, 1);
    if(specie=="d"){col_specie="red"};if(specie=="c"){col_specie="green3"};if(specie=="u"){col_specie="blue"}
    if(specie=="d"){pchs=17};if(specie=="c"){pchs=20};if(specie=="u"){pchs=4}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar, list(lab.cex=.7, pch=pchs, col=col_specie, lab.col=col_specie))
  }
  return(n)
}

# apply this to  dendrogram
dL <- dendrapply(dhc, colLab)


plot(dL , main= paste("Clustering with unknown samples: dist =", dist_method, ", clu =", clust_method ))

legend("topright", 
       legend = c("Disease", "Control", "Unknown"), 
       col = c("red", "green3", "blue"), 
       pch = c(17, 20, 4), 
       bty = "n")



########################################################################
# Clustering with heatmap

library(pheatmap)
pheatmap(mm)
pheatmap(tt , show_rownames = F, show_colnames = F, border_color = NA)
?pheatmap

pheatmap(tt , border_color = NA, cex = .8)
pheatmap(tt , border_color = NA, cex = .8, cutree_rows = 6)

pheatmap(tt , border_color = NA, cex = .8, clustering_distance_cols = "euclidean")
pheatmap(tt , border_color = NA, cex = .8, clustering_distance_cols = "manhattan")
pheatmap(tt , border_color = NA, cex = .8, clustering_distance_cols = "manhattan", cutree_cols = 4)


########################################################################
# tSNE


library(Rtsne)
rtsne_out <- Rtsne(tt , perplexity = .4)
plot(rtsne_out$Y, pch=substr(rownames(tt),1,1), col = km$cluster)


########################################################################
# END OF SCRIPT

