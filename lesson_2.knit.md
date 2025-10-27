---
title: "Lesson 2"
subtitle: ""
author:
  - "Tutorial written by Kevin Burfeind"
  - "Tutorial revised by Jonathan Nelson"
  - "Tutorial adapted by Sienna Blanche"
date: "27 October 2025"
output:
  html_document: 
    toc: yes
    toc_depth: 5
    toc_float: true
    number_sections: no
    df_print: paged
    code_folding: hide
    highlight: pygments
---



At this point everyone should have RStudio installed along with Seurat. Please refer to the [Downloading R ](https://usckrc.github.io/Education_Website/downloading_R.html) section of the website if you need to do that again.

# Getting Started

Create the "snRNA_seq_learning" folder on your desktop if you have not already done so. Properly set your working directory to the **snRNA_seq_learning folder** by clicking on **Session>Set Working Directory>Choose Directory**. This tells R studio where our downloaded files are.

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_1/ROC_1_Screenshot 0.png)

When the file browser opens, navigate to and  **select the scRNAseq_learning folder**, then select open. This sets the working directory so RStudio knows where to look for your files.

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_1/ROC_1_Screenshot 1.1.png)

## Necessary Materials

You can download the necessary file from the [GEO entry](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111107) and click the sample titled **"Glom_rep1"(GSM3022239)** then download it by clicking on the "ftp" link next to **GSM3022239_dge_glom_rep1.txt.gz** at the bottom of the page.

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 1.png)

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 2.png)

Drag and drop the file into the "snRNA_seq_learning" folder. Double click on the file to unzip it. Do the same for the sample titled **"Glom_rep2" (GSM3022240)**.

Alternatively, you can download it from [Google Drive] NEED (https://drive.google.com/file/d/16JzkZN2qmUzK8bg1_D2x3vRstWrqvHUE/view?usp=sharing) (**recommended**). To follow along with this tutorial you can download the [R Script] NEED (https://drive.google.com/file/d/1RODN4wE31qgnaK6_WjnC1BihtVgX7fMu/view?usp=sharing) and run it or create your own and copy the code from the tutorial into the script.

# Step 1 

## Activate the R packages needed to run the analysis.

Copy the following lines of code into your file and run them:

**library(Seurat)**

**library(dplyr)**

**library(patchwork)**

# Step 2

## Upload the data 

We will be uploading the data slightly differently this time, by using the read.table function. Copy and paste the following code which converts the .txt file into a data file.

**Glom_rep1.data <- read.table(file = "snRNAseq_learning/GSM3022239_deg_glom_rep1.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric")))**

Glom_rep1.data is a data frame that can be utilized in several different ways. We will next create a seurat object out of the data frame using the following code.

**Glom_rep1 <- CreateSeuratObject(counts = Glom_rep1.data, project = "Glom_merged", min.cells = 3, min.features = 200)**

**Glom_rep1**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 3.png)

Now Glom_rep1 is an object of class Seurat with 17511 features (genes) across 4556 samples (cells). Since we are planning to merge this dataset with another, we will need to check to see if its identity is coded in the metadata.

**head(Glom_rep1@meta.data)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 4.png)

As you can see, the column orig.ident contains "Glom_merged", which is is the project name. Once we add another dataset, that will be overwritten, and can't be used to pull out data. We will add a column to the metadata that denotes the replicate.

**Glom_rep1 <- AddMetaData(object = Glom_rep1, metadata = "rep1", col.name = "replicate")**

**head(Glom_rep1@meta.data)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 5.png)

Next, make the second second replicate a data frame.

**Glom_rep2.data <- read.table(file = "Glom/GSM3022240_dge_glom_rep2.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric")))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 6.png)

We will now create a seurat object out of Glom_rep2 as we did for Glom_rep1.

**Glom_rep2 <- CreateSeuratObject(counts = Glom_rep2.data, project = "Glom_merged", min.cells = 3, min.features = 200)**

**Glom_rep2**

Now Glom_rep2 is an object of class Seurat with 18309 features (genes) across 4366 samples (cells). We will next add a column to the metadata indicate the replicate.

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 7.png)

**Glom_rep2 <- AddMetaData(object = Glom_rep2, metadata = "rep2", col.name = "replicate")**

**head(Glom_rep2@meta.data)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 8.png)

# Step 3:

## Merge the two Seurat objects. 

To merge the two objects, we will use the "merge" function, but there are several different ways to merge Seurat objects. Note: the "merge" function does not work when datasets are very different from one another - i.e., stim cell vs. unstim cell. However, It does seem to work when comparing disease vs. no disease

**Glom_merged <- merge(Glom_rep1, y = Glom_rep2, add.cell.ids = c("Glom_rep1", "Glom_rep2"), project = "Glom_merged")**

**Glom_merged**

Glom_merged is now an object with 19183 features (genes) across 8922 samples (cells). Note that the sample number is the sum of the two reps: 4556 + 4366 = 8922.

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 9.png)

Also notice the cell names now have an added identifier

**head(colnames(Glom_merged))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 10.png)

We can also double check the number of cells in each group using the "table" function, pulling out "replicate".

**table(Glom_merged$replicate)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 11.png)

Now, begin the QC that we did last lesson, starting with filtering out samples with a high % of mito genes. This adds a column in the metadata called "percent.mito".

**Glom_merged[["percent.mt"]] <- PercentageFeatureSet(Glom_merged, pattern = "^mt-")**

**head(Glom_merged@meta.data)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 12.png)

As with last lesson, visualize the dataset for quality control using 3 different metrics: 

nFeature_RNA = Number of differing genes detected in each cell
nCount_RNA = Number of mRNA molecules detected in each cell
percent.mt = Percent of counts that come from mitochondrial genes

**VlnPlot(Glom_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 13.png)

We will interpret this information to filter out doublets and damaged cells. High nFeature_RNA or nCount_RNA could equate to doublets or mixed cell debris. We will filter out cells with more than 2000 genes and 4000 mRNA molecules. High percent.mt could equate to damaged cells. We will filter out cells with more than 20% mitochondrial content.

**Glom_merged <- subset(Glom_merged, subset = nFeature_RNA < 2000 & nFeature_RNA < 4000 & percent.mt < 20)**

Now we will replot the graphs to visualize the plots again with the newly filtered data

**VlnPlot(Glom_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)**

**Glom_merged**

Glom_merged is now an object with 19183 features (genes) across 8361 samples (cells). This filtered out 561 cells and 0 genes!

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 14.png)

# Step 4

## Data Normalization and Identification of Genes Driving Variability (cluster formation)

Normalize the matrix and create a plot of the genes that are the most variable in the dataset.

**Glom_merged <- NormalizeData(Glom_merged, normalization.method = "LogNormalize", scale.factor = 10000)**

**Glom_merged <- FindVariableFeatures(Glom_merged, selection.method = "vst", nfeatures = 2000)**

This line identifies the 10 most highly variable genes and creates a new object in the enviornment called "top10".

**top10 <- head(VariableFeatures(Glom_merged), 10)**

These lines creates a plot of the variable features within the dataset and labels the top 10 genes to see which are driving the clustering.

**plot1 <- VariableFeaturePlot(Glom_merged)**

**plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)**

**plot2**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 15.png)

# Step 5

## Creation and Visualization of Principal Components

Calculate the principal components and then print the 5 genes that drive the first 5 principal components. Can you get a sense for how the cells might be clustering based off of these genes? Note: the "ScaleData" and "RunPCA" steps can take very long with large datasets

**all.genes <- rownames(Glom_merged)**

**Glom_merged <- ScaleData(Glom_merged, features = all.genes)**

**Glom_merged <- RunPCA(Glom_merged, features = VariableFeatures(object = Glom_merged))**

**print(Glom_merged[["pca"]], dims = 1:5, nfeatures = 5)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 16.png)

Now we will visualize the genes driving the PCA a few ways. This line visualizes the genes driving the first 2 PCA dimensions.

**VizDimLoadings(Glom_merged, dims = 1:2, reduction = "pca")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 17.png)

This line visualizes the scRNAseq dataset as a PCA (hopefully you can appreciate the poor clustering).

**DimPlot(Glom_merged, reduction = "pca")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 18.png)

This line creates a heatmap of the top 500 genes that drive the first 2 PCA dimensions.

**DimHeatmap(Glom_merged, dims = 1:2, cells = 500, balanced = TRUE)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 19.png)

This plot shows the variability and helps to decide how many PCA dimensions to use when measuring the clusters. When the standard deviation in the PC dimensions becomes small it no longer has a big effect. In this case we will use 10 PC's.

**ElbowPlot(Glom_merged)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 20.png)

# Step 6

## Calculation and Creation of a UMAP visualization

Determine the shape of the clusters and the number of populations in our scRNASeq visualization. These lines control the resolution of the clusters and the number of separate populations.

"FindNeighbors" changes the separation between clusters
"FindClusters" changes the number of differntial populations

**Glom_merged <- FindNeighbors(Glom_merged, dims = 1:10)**

**Glom_merged <- FindClusters(Glom_merged, resolution = 0.15)**

**Glom_merged <- RunUMAP(Glom_merged, dims = 1:10)**

This line generates the UMAP visualization. Take a moment to reflect on how do you think the clusters look.

**DimPlot(Glom_merged, reduction = "umap")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 21.png)

The first sanity/QC check for merged datasets is to see if cells are clustering based on replicate. Since we have replicate info in the metadata, we can check this using the "group.by" function.

**DimPlot(Glom_merged, reduction = "umap", group.by = "replicate")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 22.png)

The replicates are integrated nicely. It does not appear there is any batch effect. If there is a batch effect, consider using the ["FindIntegrationAnchors"](https://satijalab.org/seurat/archive/v3.0/immune_alignment.html) function, which we will not address in today's lesson You can also decrease the number of variable features in the "FindVariableFeatures" function. Right now it is set to 2000, decreasing it will help deal with batch effect, but will decrease the resolution of clustering. 

We will next go through the cell identification features we covered in the last lesson. First generation Heatmaps of differential genes between clusters to identify cell types. This can take a VERY long time with large datasets

**Glom_merged[["RNA"]] <- JoinLayers(Glom_merged[["RNA"]])**

**Glom.markers <- FindAllMarkers(Glom_merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)**

**Glom.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)**

**top5 <- Glom.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)**

**DoHeatmap(Glom_merged, features = top5$gene)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 23.png)

This looks very similar to the heatmap we generated from the dataset last lesson, let's just double check that the identities are the same.

Nphs2 = Podocyte
Pecam1 = Endothelial
Slc12a3 = Tubule
Tagln = Mural

**VlnPlot(Glom_merged, features = c("Nphs2", "Pecam1", "Slc12a3", "Tagln"))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 24.png)

It looks like there are two different podocyte populations. We can check out there differences between clusters 0 and 1 using the "FindMarkers" function.

**cluster0.markers <- FindMarkers(Glom_merged, ident.1 = 0, ident.2 = 1, min.pct = 0.25)**

**head(cluster0.markers, n = 25)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 25.png)

Using this information, we can visualize the differences using the "Violin Plot" function.

**VlnPlot(Glom_merged, features = c("Ctgf"))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 26.png)

We won't spend a lot of time trying to identify different podocyte populations, but you could look at more distinct populations, such as endothelial (population 2) and epithelial (population 3) cells.

**cluster2.markers <- FindMarkers(Glom_merged, ident.1 = 2, ident.2 = 3, min.pct = 0.25)**

**head(cluster2.markers, n = 10)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 27.png)

# Step 7

## Re-name the clusters by cell-type, and add that to the metadata

**new.cluster.ids <- c("Podo1", "Podo2", "Endo", "Tubule", "Mural")**

**names(new.cluster.ids) <- levels(Glom_merged)**

**Glom_merged <- RenameIdents(Glom_merged, new.cluster.ids)**

**DimPlot(Glom_merged, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE) + NoLegend()**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 28.png)

Note that Podo1 and Podo2 are right next to each other, yet the Endo population is split apart. I personally think this means the two podocyte populations are not significant. The other Endo "cluster" may represent doublets, but we will not spend much time addressing this. If you want, you can spend time with the "FeaturePlot" function to visualize what genes are expressed in the other Endo "Cluster", then adjust your "FindNeighbors" and "FindClusters" values a bit.

**FeaturePlot(Glom_merged, features = c("Nphs2"))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 29.png)

# Step 8

## Subclustering

Let's do subclustering of the Tubule cells to demonstrate how one can identify sub populations by isolating a certain cluster then reclustering those cells. First, create a subset out of the tubule cells using the "subset" function. We will make a new Seurat object called "Glom_tubule".

**Glom_tubule <-subset(Glom_merged, idents = c("Tubule"))**

**Glom_tubule**

Glom_tuble is an object with 19183 features (genes) across 556 samples (cells).

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 30.png)

Next we will re-cluster the cells in Glom_tubule, starting with data normalization.

**Glom_tubule <- NormalizeData(Glom_tubule, normalization.method = "LogNormalize", scale.factor = 10000)**

**Glom_tubule <- FindVariableFeatures(Glom_tubule, selection.method = "vst", nfeatures = 2000)**

Similarly, we will re-scale the data. Note that the row names are now based on the row names of the new Seurat object.

**tubule.genes <- rownames(Glom_tubule)**

**Glom_tubule <- ScaleData(Glom_tubule, features = tubule.genes)**

**Glom_tubule <- RunPCA(Glom_tubule, features = VariableFeatures(object = Glom_tubule))**

**print(Glom_tubule[["pca"]], dims = 1:5, nfeatures = 5)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 31.png)

Visualize the genes driving the first two PCAs with the following code.

**VizDimLoadings(Glom_tubule, dims = 1:2, reduction = "pca")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 32.png)

Also visualize the PCA.

**DimPlot(Glom_tubule, reduction = "pca")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 33.png)

Create a heatmap of the top 500 genes that drive the first 2 PCA dimensions. 

**DimHeatmap(Glom_tubule, dims = 1:2, cells = 500, balanced = TRUE)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 34.png)

Visualize the variance induced by each PC using the "Elbow Plot" function.

**ElbowPlot(Glom_tubule)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 35.png)

Now we can recluster the cells. Based on the elbow plot, we will set the dims to 10 again.

**Glom_tubule <- FindNeighbors(Glom_tubule, dims = 1:10)**

**Glom_tubule <- FindClusters(Glom_tubule, resolution = 0.1)**

**Glom_tubule <- RunUMAP(Glom_tubule, dims = 1:10)**

This line generates the UMAP visualization. How many clusters do you see? Do they seem like they are real clusters?

**DimPlot(Glom_tubule, reduction = "umap")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 36.png)

Check to make sure the replicates aren't driving the clustering.

**DimPlot(Glom_tubule, reduction = "umap", group.by = "replicate")**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 44.png)

Let's do a sanity check just to make sure we aren't doing the same cluster all over again. Do a violin plot of the same genes that defined the original clusters.

**VlnPlot(Glom_tubule, features = c("Nphs2", "Pecam1", "Slc12a3", "Tagln"))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 37.png)

It looks like there may be one population that contains podocyte-tubule cell doublets, both otherwise looks pretty good. Now let's generate another heatmap of differential genes between clusters to identify cell types.

**Tubule.markers <- FindAllMarkers(Glom_tubule, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)**

**Tubule.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)**

**top5_tubule <- Tubule.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)**

**DoHeatmap(Glom_tubule, features = top5_tubule$gene)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 38.png)

Pvalb, Egf, and Umod are all expressed by DCT.

Aqp2 is a marker of principal cells (PC) of the Collecting Duct.

Wt1 is expressed by podocytes.

Slc4a1 is expressed by intercalated cells (IC) of the collecting duct.

Let's check with a violin plot.

**VlnPlot(Glom_tubule, features = c("Aqp2", "Egf", "Nphs2", "Atp6v1g3"))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 39.png)

Note the Egf is expressed in population 0 and 2, which suggests that population 2 containes DCT/Podocyte doublets. Next, re-name the clusters by cell-type, and add that to the metadata.

**new.cluster.ids.tubule <- c("DCT", "PC", "Doublet", "IC")**

**names(new.cluster.ids.tubule) <- levels(Glom_tubule)**

**Glom_tubule <- RenameIdents(Glom_tubule, new.cluster.ids.tubule)**

**DimPlot(Glom_tubule, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 7, repel = FALSE) + NoLegend()**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 40.png)

Now let's get to know the clusters a little better. To determine how many cells are in each cluster we can run the following code.

**table(Idents(Glom_tubule))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 41.png)

To determine the proportion of cells are in each cluster, run the following code.

**prop.table(table(Idents(Glom_tubule)))**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 42.png)

To determine the how does cluster membership vary by replicate, run the following code.

**table(Idents(Glom_merged), Glom_merged$replicate)**

**table(Idents(Glom_tubule), Glom_tubule$replicate)**

![](/Users/sblanche/Desktop/GitHub/Education_Website/images/lesson_2/roc_2_screenshot 43.png)

We can also save the "Tubule.markers" as a .csv so we can look over it in excel, or use it for downstream analysis.

**write.csv(Tubule.markers, "Glom/Tubule.markers.csv")**

Lastly, we will save the "Glom_tubule" object so it can be loaded later. 

# Step 9

## Save Seurat object. 

Will save into the folder that the Rscript file is located in.

**saveRDS(Glom_tubule, file = "Glom/Glom_tutorial.rds")**

To load a SeuratObject that was saved as an RDS, run the following code.

**Glom_tubule <- readRDS("Glom/Glom_tutorial.rds")**

# Extra code for uploading different data files to Seuat

All three of these examples start with turning the data into a data frame which Seurat can read, then creating a Seurat object out of the data frame. Barcodes, features, and matrix files - common 10X output. Make sure the three files (usually zipped) are in one folder. They must be titled "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz".

First upload the data using the "Read10X" function.

**Object.data <- Read10X(data.dir = "Path to file")**

## H5 files

To load H5 files, you will need to install the "hd5r" package.

**install.packages("hdf5r")**

**library(hdf5r)**

to upload the data, it must end in ".h5".

**Object.data <- Read10X_h5("Path to .h5 file")**

## DGE table

To read DGE tables, use the the "read.table" function.The file must be in ".txt" format. It will likely be in ".dge.txt" format.

**Object.data <- read.table(file = "Path to .dge.txt", header = TRUE, row.names = 1, colClasses =c("character", rep("numeric")))**

For all of these, make a Seurat object out of the data frame.

**Object_name <- CreateSeuratObject(counts = Object.data, min.cells = 3, min.genes = 200, project = "Project_name")**

# Session Info


``` r
sessionInfo()
```

```
## R version 4.4.2 (2024-10-31)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sequoia 15.6.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Los_Angeles
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] kableExtra_1.4.0   openxlsx_4.2.8     viridis_0.6.5      viridisLite_0.4.2 
##  [5] gplots_3.2.0       reshape2_1.4.4     pheatmap_1.0.13    ggvenn_0.1.10     
##  [9] lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1      purrr_1.1.0       
## [13] readr_2.1.5        tidyr_1.3.1        tibble_3.3.0       tidyverse_2.0.0   
## [17] devtools_2.4.5     usethis_3.1.0      here_1.0.1         ggpmisc_0.6.2     
## [21] ggpp_0.5.9         ggplot2_3.5.2      knitr_1.50         patchwork_1.3.1   
## [25] Seurat_5.3.0       SeuratObject_5.1.0 sp_2.2-0           dplyr_1.1.4       
## 
## loaded via a namespace (and not attached):
##   [1] RcppAnnoy_0.0.22       splines_4.4.2          later_1.4.2           
##   [4] bitops_1.0-9           polyclip_1.10-7        fastDummies_1.7.5     
##   [7] lifecycle_1.0.4        rprojroot_2.1.0        globals_0.18.0        
##  [10] lattice_0.22-7         MASS_7.3-65            magrittr_2.0.3        
##  [13] plotly_4.11.0          sass_0.4.10            rmarkdown_2.29        
##  [16] jquerylib_0.1.4        yaml_2.3.10            remotes_2.5.0         
##  [19] httpuv_1.6.16          sctransform_0.4.2      spam_2.11-1           
##  [22] zip_2.3.3              sessioninfo_1.2.3      pkgbuild_1.4.8        
##  [25] spatstat.sparse_3.1-0  reticulate_1.43.0      cowplot_1.2.0         
##  [28] pbapply_1.7-4          RColorBrewer_1.1-3     abind_1.4-8           
##  [31] pkgload_1.4.0          Rtsne_0.17             ggrepel_0.9.6         
##  [34] irlba_2.3.5.1          listenv_0.9.1          spatstat.utils_3.1-5  
##  [37] MatrixModels_0.5-4     goftest_1.2-3          RSpectra_0.16-2       
##  [40] spatstat.random_3.4-1  fitdistrplus_1.2-4     parallelly_1.45.1     
##  [43] svglite_2.2.1          codetools_0.2-20       xml2_1.3.8            
##  [46] tidyselect_1.2.1       farver_2.1.2           matrixStats_1.5.0     
##  [49] spatstat.explore_3.5-2 jsonlite_2.0.0         ellipsis_0.3.2        
##  [52] progressr_0.15.1       ggridges_0.5.6         survival_3.8-3        
##  [55] systemfonts_1.2.3      tools_4.4.2            ica_1.0-3             
##  [58] Rcpp_1.1.0             glue_1.8.0             gridExtra_2.3         
##  [61] xfun_0.52              withr_3.0.2            fastmap_1.2.0         
##  [64] SparseM_1.84-2         caTools_1.18.3         digest_0.6.37         
##  [67] timechange_0.3.0       R6_2.6.1               mime_0.13             
##  [70] textshaping_1.0.1      colorspace_2.1-1       scattermore_1.2       
##  [73] gtools_3.9.5           tensor_1.5.1           spatstat.data_3.1-6   
##  [76] generics_0.1.4         data.table_1.17.8      httr_1.4.7            
##  [79] htmlwidgets_1.6.4      uwot_0.2.3             pkgconfig_2.0.3       
##  [82] gtable_0.3.6           lmtest_0.9-40          htmltools_0.5.8.1     
##  [85] profvis_0.4.0          dotCall64_1.2          scales_1.4.0          
##  [88] png_0.1-8              spatstat.univar_3.1-4  rstudioapi_0.17.1     
##  [91] tzdb_0.5.0             nlme_3.1-168           cachem_1.1.0          
##  [94] zoo_1.8-14             KernSmooth_2.23-26     parallel_4.4.2        
##  [97] miniUI_0.1.2           pillar_1.11.0          vctrs_0.6.5           
## [100] RANN_2.6.2             urlchecker_1.0.1       promises_1.3.3        
## [103] xtable_1.8-4           cluster_2.1.8.1        evaluate_1.0.4        
## [106] cli_3.6.5              compiler_4.4.2         rlang_1.1.6           
## [109] future.apply_1.20.0    plyr_1.8.9             fs_1.6.6              
## [112] stringi_1.8.7          deldir_2.0-4           lazyeval_0.2.2        
## [115] spatstat.geom_3.5-0    quantreg_6.1           Matrix_1.7-3          
## [118] RcppHNSW_0.6.0         hms_1.1.3              future_1.67.0         
## [121] shiny_1.11.1           ROCR_1.0-11            igraph_2.1.4          
## [124] memoise_2.0.1          bslib_0.9.0            polynom_1.4-1
```
