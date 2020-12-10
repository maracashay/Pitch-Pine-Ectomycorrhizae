library(vegan)

rem.ps <- subset_samples(ps.rare, Treatment == "NR")
rem.ps.df<- data.frame(otu_table(rem.ps))
env.data<- data.frame(sample_data(rem.ps))
env.data[,1:5]<-NULL
cca.all<- cca(esv.rare ~ ., data = env.data) # build a CCA model with all variables

mod0<- cca(esv.full ~ 1, env.data) # create a null model
mod <- step(mod0, scope = formula(cca.all), test = "perm", perm.max = 100)

rem.ps.df <- subset_samples(ps, Treatment == "NR")
rem.df <- data.frame(otu_table(rem.ps.df))
rem.meta.1 <- data.frame(sample_data(rem.ps.df))
pca.p <- rda(rem.df, scale=TRUE)
ef <- envfit(pca.p, rem.meta, permu=99)
ef

dfjd <- vegdist(rem.df, method = "bray", scale = TRUE)
dfpca <- cmdscale(dfjd, eig = TRUE, k = 2)
explainvar1 <- round(dfpca$eig[1] / sum(dfpca$eig), 2) * 100
explainvar2 <- round(dfpca$eig[2] / sum(dfpca$eig), 2) * 100
sum.eig <- sum(explainvar1, explainvar2)


colvecC <-c("black", "gray60", "red4", "dodgerblue", "cyan", "magenta", "chartreuse", "goldenrod", "purple", "springgreen4")

plot(dfpca$points[ ,1], dfpca$points[ ,2],
     xlab = paste("PCoA 1 (", explainvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainvar2, "%)", sep = ""),
     pch = 16, cex = 4, type = "n", cex.lab = 1.5, cex.axis = 1.2, axes = FALSE)
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)
points(dfpca$points[ ,1], dfpca$points[ ,2],pch = 19, cex = 2, col = colvecCv[rem.meta.1$Site])
legend("topleft", inset=c(0,0), legend=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), col=c("black", "gray60", "red4", "dodgerblue", "cyan", "magenta", "chartreuse", "goldenrod", "purple", "springgreen4"), pch=19, cex=2, ncol=2, title="Site")
