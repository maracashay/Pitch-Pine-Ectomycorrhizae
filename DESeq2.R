gen.fe <- phyloseq_to_deseq2(full.phyloseq, ~ Treatment) # converts the phyloseq object into a Deseq object
gen.fe.counts <- counts(gen.fe) # contains sequence counts for each genera
geoMeans = apply(gen.fe.counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.gen.fe <-estimateSizeFactors(gen.fe, geoMeans=geoMeans) 
dds.gen.fe <- DESeq(dds.gen.fe, test="Wald", fitType="parametric")

res.gen.fe = results(dds.gen.fe, cooksCutoff = FALSE) # extracts results table from DESeq analysis
alpha=0.05 # set alpha value 
sigtab.gen.fe = res.gen.fe[which(res.gen.fe$padj < alpha), ] # deseq output that excludes taxa that have p-values greater than alpha
sigtab.gen.fe #this file includes only the genera that were significantly 
