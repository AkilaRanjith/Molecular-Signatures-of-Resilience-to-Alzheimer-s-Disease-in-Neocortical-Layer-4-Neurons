---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library("SingleCellExperiment"))
library(ComplexHeatmap)

#### read the RDS file low, intermediate, high disease.Group file seperately and run the below code
###reading low
sce<-readRDS("/scratch/groups/icobos/Akila/recompute/cellchat/BA9/low_BA9.rds")
sf <- 2^rnorm(ncol(sce))
sf <- sf/mean(sf)
normcounts(sce) <- t(t(counts(sce))/sf)
dim(normcounts(sce))
# computing log-counts
logcounts(sce) <- log2(normcounts(sce)+1)
dim(normcounts(sce))
meta<-as.data.frame(colData(sce))
cellchat <- createCellChat(object = sce, meta = meta, group.by = "Author_Annotation")

### read the receptor-ligand database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

### compute communication probability 
cellchat <- computeCommunProb(cellchat,"triMean",population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,"/scratch/groups/icobos/Akila/recompute/cellchat/BA9/low/cellchat_network_low1.csv")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
mat <- cellchat@net$weight
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,15)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
saveRDS(cellchat,"/BA9/low/cellchat_low.rds")
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 10, units = 'in', dpi = 300)
}
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.



```{r}
#### Differential interaction Analysis


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library("SingleCellExperiment"))
library(ComplexHeatmap)


###read the cellchat object from low,int,high disease group based on brain region
cellchat_low<-readRDS("BA9/low/cellchat_low.rds")
cellchat_int<-readRDS("BA9/int/cellchat_int.rds")
cellchat_high<-readRDS("BA9/high/cellchat_high.rds")

### merge low,int cellchat object to find the differential interactions

### do the same thing for comparing the int and high pathology change the object.list accordingly

object.list <- list(low = cellchat_low, int = cellchat_int)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


pdf(file = "low_int/numberofinteraction_low_int.pdf")
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()

#### differentail interaction between neuronal popultions
pdf(file = "BA9/low_int/Neuron_heatmap.pdf")
gg1 <- netVisual_heatmap(cellchat, font.size = 6, row.show = c(5:41), col.show = c(5:41),width=15,height=15)
gg1
dev.off()


#### differential interaction between major_celltypes

group.cellType1 <- factor(cellchat@meta$major_celltypes)
group.cellType<-unique(group.cellType1)
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf(file = "/major_celltype_differential_interaction1.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = F, vertex.weight = 10,vertex.label.cex = 0.7,margin = 0.3,arrow.width = 0.7,edge.width.max = 5)

#### plot the enriched pathway

pdf(file = "/BA9/Enriched_pathway.pdf")
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

#### comparison between low/int
pos.dataset = "int"
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up1 <- subsetCommunication(cellchat, net = net, datasets = "int",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down1 <- subsetCommunication(cellchat, net = net, datasets = "low",ligand.logFC = -0.1, receptor.logFC = -0.1)

#### comparison between int/high
pos.dataset = "high"
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up1 <- subsetCommunication(cellchat, net = net, datasets = "high",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down1 <- subsetCommunication(cellchat, net = net, datasets = "int",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up1, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down1, cellchat)

#### save the upregulated and down regulated L-R
write.csv(gene.down,"int_high_genelist_down.csv")
write.csv(gene.up,"int_high_genelist_up.csv")

write.csv(net.down1,"low_int_downregulated_pairs.csv")
write.csv(net.up1,"low_int_upregultaed_pairs.csv")

#### get the number of interactions based on low,int,high disease group for BA9
low<-cellchat_low@net$count
int1<-cellchat_int@net$count
high<-cellchat_high@net$count
write.csv(low,"low_interaction_counts.csv")
write.csv(int1,"int_interaction_counts.csv")
write.csv(high,"high_interaction_counts.csv")

```

