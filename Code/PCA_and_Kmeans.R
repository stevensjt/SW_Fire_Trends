#Install packages
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(corrplot)
library(cluster)
library(missMDA)
library(Hmisc)
library(tidyverse)
library(cluster)


options(scipen=999)

#Attach data
setwd("D:/Patch_Metrics/R/PCA/")
all_patch_metrics.df <- read.csv(file="all_metrics_1984_2015_my_apriori_wfu_corrected.csv", header=TRUE, sep=",")
#all_patch_metrics.df <- read.csv(file="all_metrics_1984_2015_my_apriori_email.csv", header=TRUE, sep=",")
names(all_patch_metrics.df)
attach(all_patch_metrics.df)

#identify which columns to use
all_metrics_active <- all_patch_metrics.df[,5:48]
#all_metrics_active <- all_patch_metrics.df[,5:43]
names(all_metrics_active)
#remove core area metrics 
all_metrics_active <- all_metrics_active[-c(26:38)]
#remove perimeter-area metrics
all_metrics_active <- all_metrics_active[-c(20:23)]
names(all_metrics_active)

#commpute PCA with factoextra
patch.pca <- prcomp(all_metrics_active, scale = TRUE)
summary(patch.pca)

#get info on PCA 
get_pca(patch.pca)
#loadings - columns are eigenvectors, linear combination ::prcomp
patch.pca$rotation
# the coordinates of the individuals on the principal coordinates
patch.pca$x[,1:7]

#Get eigenvalues and visulize eigenvalues with scree plot
get_eig(patch.pca)
fviz_eig(patch.pca, addlabels = TRUE, ylim = c(0,50))

#keep components 1-7 bc they have eigenvalue > 1. In this case 7 components explain 86.14% of variability, 
patch.var <- get_pca_var(patch.pca)
patch.var$contrib[,1:7] #the larger the value of the variable contribution, the more it explains its variability
patch.var$cor[,1:7] #correlation between variables and dimensions


#biplot of variables
fviz_pca_var(patch.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


# Visualize Contributions of variables to PC1
fviz_contrib(patch.pca, choice = "var", axes = 1, top = 44)
# Visualize Contributions of variables to PC2
fviz_contrib(patch.pca, choice = "var", axes = 2, top = 44)
#Visualize Total Contributions of variables to PC1 - PC7
fviz_contrib(patch.pca, choice = "var", axes = 1:7, top = 44)


#cluster analysis - scale data to get variables in same units
kmeans_patch <- scale(all_metrics_active)
kmeans_patch
#cluster analysis - get optimal number of clusters. returns 10 clusters
set.seed(123)
fviz_nbclust(kmeans_patch, kmeans, method = "gap_stat", k.max = 20, iter.max=30)


########cluster analysis - kmeans cluster 
res.pca <- PCA(all_metrics_active, scale.unit = TRUE, ncp = 10, graph = FALSE)
var <- get_pca_var(res.pca)
var
#coordinates of variables to create a scatter plot
head(var$coord)

#####kmeans cluster visualization#####
# Create a grouping variable using kmeans
# Create 10 groups (centers = 10)
set.seed(123)
res.km <- kmeans(var$coord, centers = 10, nstart = 25)
grp <- as.factor(res.km$cluster)
# Color variables by groups
fviz_pca_var(res.pca, col.var = grp, 
             palette = c("black", "red", "blue", "darkorange", "coral4", 
                         "darkgreen", "purple", "darkslategray", "gold"),
             repel = TRUE,
             legend.title = "Cluster")

fviz_pca_var(res.pca, col.var = grp, 
             palette = c("black", "red", "blue", "darkorange", "coral4", 
                         "darkgreen", "purple", "darkslategray", "gold",
                         "brown","deeppink3", "navy"),
             repel = TRUE,
             legend.title = "Cluster")

###########################################################
