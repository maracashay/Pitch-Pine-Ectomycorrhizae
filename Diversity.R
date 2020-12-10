# Alpha-diversity
set.seed(500)
ps.rare <- rarefy_even_depth(ps.phyloseq, 1413)

sample_sums(ps.rare)

library(microbiome)

alpha.div <- estimate_richness(ps.rare, measures = c("Shannon", "Observed"))

meta.rare <- data.frame(sample_data(ps.rare))

meta.rare$Shannon <- alpha.div$Shannon
meta.rare$Richness <- alpha.div$Richness

shannon.model <- kruskal.test(meta.rare$Richness ~ meta.rare$Group)
dunnTest(meta.rare$Richness ~ meta.rare$Group, method = "bh")

# Beta-diversity heatmap

# use taxa_level function from MicrobiomeSeq package

ps.gen <- taxa_level(ps.rare, "Genus")

plot_heatmap(ps.gen, "PCoA", "bray", "Site.Treatment", na.value = "grey", low = "#66CCFF", high = "#000033")

# Beta-diversity PERMANOVA
library(vegan)
bc <- phyloseq::distance(otu_table(ps.rare), "bray")
adonis(bc ~ Treatment, data = meta.rare)

library(RVAideMemoire)
pairwise.perm.anova(bc, meta.rare$Treatmet, nperm = 999)

