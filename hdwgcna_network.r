---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
library('Seurat')
library('scater')
library('patchwork')
library('hdWGCNA')
library('SingleCellExperiment')
library('zellkonverter',lib="/home/users/akila20/R/x86_64-pc-linux-gnu-library/4.2/")
library("WGCNA")
library('cowplot')
library('harmony')
library('igraph')
library('ggrepel')
library('Matrix')
library('viridis')
library('RColorBrewer')
library('dplyr')
library('reshape2')
library('igraph')
fig_dir<-"/scratch/groups/icobos/Akila/recompute/hdwgcna/celltype_inter/figure/"
set.seed(12345)

###convert h5ad file to seurat object using reticulate

###read the converted seurat object
seurat_obj1<-readRDS("/hdwgcna/seurat_obj.rds")

#### find HVG
seurat_obj1 <- FindVariableFeatures(object = seurat_obj1)
groups <- as.character(unique(seurat_obj1$Author_Annotation))

rename1<-c(	"Astro1","Astro2","Astro3","Astro4","OPC","Micro1","Micro2","Micro3","Micro4","oligo1","oligo2","VLMC","endothelial","pericytes","Ex8",	"Ex9",	"Ex16",	"Ex13",	"Ex10",	"Ex14",	"Ex3",	"Ex15",	"Ex1",	"Ex2",	"Ex18",	"Ex17",	"Ex7",	"Ex4",	"Ex5",	"Ex12",	"Ex6",	"Ex11",	"In2",	"In17",	"In7",	"In1",	"In3",	"In6",	"In15",	"In16",	"In14",	"In19",	"In18",	"In11",	"In4",	"In9",	"In5",	"In8",	"In12",	"In13",	"In10")

for(i in 1:51)
{
  print(i)
  g <-as.character(groups[i])
  print(g)
  cur_group <- as.character(groups[i])
  rename <-as.character(rename1[i])
  print(rename)
  seurat_obj <- subset(x = seurat_obj1, subset = (author_annotation == as.character(groups[i])))
    
  # only run metacell aggregation one time
    print('Setting up and running metacells')
    # setup training data for scWGCNA:
    seurat_obj <- SetupForWGCNA(
      seurat_obj,
      gene_select = "fraction",
      fraction = 0.05,
      wgcna_name = rename
    )
    saveRDS(seurat_obj, paste0("/hdwgcna/celltype_inter/",cur_group,"_setup.rds"))
    print('Done set up')
      # construct metacells:
    seurat_obj <- MetacellsByGroups(
      seurat_obj = seurat_obj,
      group.by = c("Author_Annotation", "subject"),
      k = 20,
      max_shared = 15,
      min_cells = 50,
      ident.group = 'Author_Annotation',
      reduction="X_pca_harmony"
     
    )
    print('done metacell')
    seurat_obj <- NormalizeMetacells(seurat_obj)
  #setup expression matrix
   seurat_obj <- SetDatExpr(seurat_obj,assay = "originalexp",group.by = "Author_Annotation",group_name=cur_group )
   print('expression set up')
    
  # construct network
  seurat_obj <- TestSoftPowers(seurat_obj)
  seurat_obj <- ConstructNetwork(seurat_obj, tom_name = rename)
  seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
  seurat_obj <- ModuleEigengenes(seurat_obj)
  seurat_obj <- ModuleConnectivity(seurat_obj, group.by = 'Author_Annotation',group_name = cur_group, wgcna_name = rename)
  print("done connnectivity")
  # rename the modules
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste0(rename, '-M'),
    wgcna_name=rename
  )
   saveRDS(seurat_obj, paste0("/hdwgcna/celltype_inter/",rename,".rds"))
### plot figures
   pdf(paste0(fig_dir,rename," dendrogram.pdf"),height=6,width=10)
   p<-PlotDendrogram(seurat_obj, main=paste0(rename,' interaction with other celltypes'))
   p
   dev.off
   
   ### plot KMe
   pdf(paste0(fig_dir,rename," plot_KME.pdf"),height=15,width=12)
   p<-PlotKMEs(seurat_obj, ncol=3,text_size=3)
   p
   dev.off()
   
   ### feature plot
   seurat_obj <- ModuleExprScore(seurat_obj,n_genes = 25,method='Seurat')
   mat <- seurat_obj[['X_pca_harmony']]@cell.embeddings *(-1)
   seurat_obj[['UMAP.rev']] <- CreateDimReducObject(embeddings = mat, key = 'TsneRe_', assay = 'originalexp')
   DimPlot(seurat_obj)
   
   ### plot hME (hetereogenous expression pattern with in module)
   plot_list <- ModuleFeaturePlot(seurat_obj,features='hMEs',order=TRUE ,reduction = "X_pca_harmony",raster = FALSE)
   wrap_plots(plot_list, ncol=2,width=8,height=10)
   pdf(paste0(fig_dir,rename," hME_module_feature_plots.pdf"),height=12,width=10)
   p <- wrap_plots(plot_list, ncol=2)
   p
   dev.off()
    
   # Plot hME (hetereogenous expression pattern within module)
   plot_list <- ModuleFeaturePlot(seurat_obj, features = 'hMEs', order = TRUE, reduction = "X_pca_harmony", raster = FALSE)
   p <- wrap_plots(plot_list, ncol = 2, width = 8, height = 10)

    # Save the plot as a PNG file
   png(paste0(fig_dir, rename, "_hME_module_feature_plots.png"), width = 800, height = 1000)
   print(p)
   dev.off()

   
   # get hMEs from seurat object
   hMEs <- GetMEs(seurat_obj, harmonized=TRUE)
   mods <- colnames(hMEs); mods <- mods[mods != 'grey']

   # add hMEs to Seurat meta-data:
   seurat_obj@meta.data <- cbind(seurat_obj@meta.data, hMEs)
   # get MEs from seurat object
   MEs <- GetMEs(seurat_obj, harmonized=TRUE)
   mods <- colnames(MEs); mods <- mods[mods != 'grey']

    # add MEs to Seurat meta-data:
   seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
    
    # plot with Seurat's DotPlot function
    pdf(paste0(fig_dir,rename," module_dotplots.pdf"),height=15,width=15)
    p <- DotPlot(seurat_obj, features=mods, group.by = 'Author_Annotation')
    
    # flip the x/y axes, rotate the axis labels, and change color scheme:
    p <- p +
    coord_flip() +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')

    # plot output
    p
    dev.off()
    ##save png file
    p <- DotPlot(seurat_obj, features = mods, group.by = 'Author_Annotation')
    print(p)
    png(paste0(fig_dir, rename, "_module_dotplots.png"), height = 1500, width = 1500)
    p
    dev.off()
    
    ### Module networks
    ModuleNetworkPlot(seurat_obj)
    old_folder_name <- "ModuleNetworks"
    new_folder_name <- paste0("ModuleNetworks",rename)
    file.rename(old_folder_name, new_folder_name)
    
    
    ### HUb genes plot
    # hubgene network
    HubGeneNetworkPlot(seurat_obj,n_hubs = 3, n_other=5,edge_prop = 0.75,mods = 'all',raster = FALSE)
    g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
    seurat_obj <- RunModuleUMAP(seurat_obj,n_hubs = 10, n_neighbors=15, min_dist=0.5,spread=1.4)
    
    
    # get the hub gene UMAP table from the seurat object
    umap_df <- GetModuleUMAP(seurat_obj)

    # plot with ggplot 
    ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) + geom_point(color=umap_df$color,size=umap_df$kME*2) + umap_theme()
    
    
    ### Hub gene in UMAP
    pdf(paste0(fig_dir,rename," hub_UMAP.pdf"),height=15,width=18)
    U<-ModuleUMAPPlot(seurat_obj,edge.alpha=0.25,sample_edges=TRUE,edge_prop=0.1, label_hubs=3,keep_grey_edges=FALSE,     vertex.label.cex = 0.4)
    U
    dev.off()
    
    ### DME analysis
    group1 <- seurat_obj@meta.data %>% subset(author_annotation == cur_group & Disease.Group == "low") %>% rownames
    group2 <- seurat_obj@meta.data %>% subset(author_annotation == cur_group & Disease.Group == "int") %>% rownames
    group3 <- seurat_obj@meta.data %>% subset(author_annotation == cur_group & Disease.Group == "high") %>% rownames
    
    ### early DME low/int
    
    early <- FindDMEs(seurat_obj, barcodes1 = group2, barcodes2 = group1,test.use='wilcox',wgcna_name=rename,min.pct=0.2)

    #### Late DME int/high
   
    late <- FindDMEs(seurat_obj,barcodes1 = group3, barcodes2 = group2, test.use='wilcox', wgcna_name=rename,min.pct=0.2)
    
    ### write csv files for DME and hub genes list
    
    write.csv(early,paste0(fig_dir,rename,"early_DME.csv"))
    write.csv(late,paste0(fig_dir,rename,"late_DME.csv"))
    
    hub_df <- GetHubGenes(seurat_obj, n_hubs = 50)
    write.csv(hub_df,paste0(fig_dir,rename,"_hubgenes.csv"))
    saveRDS(seurat_obj, paste0("/hdwgcna/celltype_inter/",rename,"processed.rds"))
}
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.
