library(vegan)

rem.ps <- subset_samples(ps.rare, Treatment == "NR")
rem.ps.df<- data.frame(otu_table(rem.ps))
env.data<- data.frame(sample_data(rem.ps))
env.data1 <- data.frame(env.data$Zn)
env.data1$Cd <- env.data$Cd

env.data1 <- log(env.data)


m <- rda(decostand(rem.ps.df, "hell") ~., env.data1)

safe_colorblind <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                     "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


plot(m2, type="n", display="sites",xlab = paste("PCA1"),
     ylab = paste("PCA2"), cex.lab = 1.5, cex.axis = 1.2) 
points(m2, display="bp", lwd=2, col="black")
text(m2, display="bp", col="black", font=2, cex = 1)
points(m2, display="sites", col=safe_colorblind[env.data$Site],cex=3, font=1, pch = 16)
text(m2, dislpay = "sites", col = "red", pos = 1)
legend("topleft", legend = levels(env.data$Site), col = safe_colorblind[env.data$Site], pch = 16, cex = 1, title = "Site")
