# These are example graphics for a timeline (Dr. Moisés Selman paper)

# Figure 1: Heatmap

# scale data to mean=0, sd=1 and convert to matrix
mtscaled <- as.matrix(scale(mtcars))

jpeg("Heatmap.jpeg", width=10, height=10, units="cm", res=300)
  palette <- colorRampPalette(c('#19FF62','#000000', '#F50106'))(256)
  heatmap(mtscaled, scale='none', col=palette, labRow=NA, labCol=NA)
dev.off()

#palette <- colorRampPalette(c('#00A835','white', '#F50106'))(256)
#heatmap(mtscaled, scale='none', col=palette, labRow=NA, labCol=NA)

#palette <- colorRampPalette(c('#F50106', 'white', '#00A835' ))(256)
#heatmap(mtscaled, scale='none', col=palette, labRow=NA, labCol=NA)

#library(gplots)

#palette <- colorRampPalette(c('#19FF62','#000000', '#F50106'))(256)
#heatmap.2(mtscaled, scale='none', col=palette, labRow=NA, labCol=NA, 
#          Rowv=TRUE, Colv=TRUE,tracecol='none')

#palette <- colorRampPalette(c('#00A835','white', '#F50106'))(256)
#heatmap.2(mtscaled, scale='none', col=palette, labRow=NA, labCol=NA, 
#          Rowv=TRUE, Colv=TRUE,tracecol='none')

#palette <- colorRampPalette(c('#F50106', 'white', '#00A835' ))(256)
#heatmap.2(mtscaled, scale='none', col=palette, labRow=NA, labCol=NA, 
#          Rowv=TRUE, Colv=TRUE,tracecol='none')



# Figure 2: Manhattan plot
library(qqman)

# Make the Manhattan plot on the gwas Results dataset
jpeg("manhattan.jpeg", width=32, height=10, units="cm", res=300)
  manhattan(gwasResults, chr="CHR", bp="BP", snp="SNP", p="P", col=c(
  "dodgerblue2", "#E31A1C", # red
               "green4",
               "#6A3D9A", # purple
               "#FF7F00", # orange
               "black", "gold1",
               "skyblue2", "#FB9A99", # lt pink
               "palegreen2",
               "#CAB2D6", # lt purple
               "#FDBF6F", # lt orange
               "gray70", "khaki2",
               "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
               "darkturquoise", "green1", "yellow4", "yellow3",
               "darkorange4", "brown"))
dev.off()


# Figure 3: TSNE plot (Single cell)

# We used the tsne of the OSCA book

library(DropletTestFiles)
raw.path <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/raw.tar.gz")
out.path <- file.path(tempdir(), "pbmc4k")
untar(raw.path, exdir=out.path)

library(DropletUtils)
fname <- file.path(out.path, "raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)

library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
  rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol)

library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce.pbmc)$ID, 
                   column="SEQNAME", keytype="GENEID")

set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]
unfiltered <- sce.pbmc

stats <- perCellQCMetrics(sce.pbmc, subsets=list(Mito=which(location=="MT")))
high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
sce.pbmc <- sce.pbmc[,!high.mito]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- high.mito

library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster=clusters)
sce.pbmc <- logNormCounts(sce.pbmc)

set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)

set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc, subset.row=top.pbmc, technical=dec.pbmc)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred="PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred="PCA")

g <- buildSNNGraph(sce.pbmc, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce.pbmc) <- factor(clust)

#plotTSNE(sce.pbmc, colour_by="label")
#plotTSNE(sce.pbmc, colour_by="label", text_by="label")
#plotTSNE(sce.pbmc, colour_by="label", add_legend=F)

jpeg("TSNE.jpeg", width=10, height=8, units="cm", res=300)
  plotTSNE(sce.pbmc, colour_by="label", text_by="label", add_legend=F) +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()







