library(factoextra)
library(FactoMineR)

metadata.t<- data.frame(sample_data(ps.pruned)) 
meta.df <- data.frame(sample_data(ps.pruned))
metadata.t[,1:5]<-NULL # we can get rid of the categorical variables with this
metadata.t$Cu.1 <- NULL
metadata.t$Hg <- NULL
metadata.t$Dist <- NULL

res.pca <- PCA(metadata.t, scale = TRUE)
quartz()
plot(res.pca) #provides a pca plot of the metals data 


fviz_screeplot (res.pca, ncp=10) # % variance explained by each PC
fviz_pca_var(res.pca) # Variables factor map
var.PCA1.plot <-fviz_contrib(res.pca, choice = "var", axes = 1)  # set axes=1 to consider only PC1
var.PCA2.plot <- fviz_contrib(res.pca, choice = "var", axes = 2)  # set axes=1 to consider only PC1
pca.plot <- fviz_pca_biplot(res.pca)
