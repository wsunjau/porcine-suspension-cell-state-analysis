Analysis in Seurat 

##01_packages_and_settings

library(dplyr)
library(Seurat)
library(patchwork)
library(future)
library(ggplot2)

library(ComplexHeatmap)
library(scCustomize)
library(cowplot)

set.seed(123)


##02_load_data 

P5.data <-Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
P5 <- CreateSeuratObject(counts = P5.data, project = "P5", min.cells = 50, min.features = 500)
P5@meta.data[,"Sample"]<-"P5"
P5

P7.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
P7 <- CreateSeuratObject(counts = P7.data, project = "P7", min.cells = 50, min.features = 500)
P7@meta.data[,"Sample"]<-"P7"
P7
#13749 features across 11235 samples within 1 assay

MC.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
MC <- CreateSeuratObject(counts = MC.data, project = "MC", min.cells = 50, min.features = 500)
MC@meta.data[,"Sample"]<-"MC"
MC
#13394 features across 13503 samples within 1 assay 

SD7_3S6.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
SD7_3S6 <- CreateSeuratObject(counts = SD7_3S6.data, project = "SD7_3S6", min.cells = 50, min.features = 500)
SD7_3S6@meta.data[,"Sample"]<-"SD7_3S6"
SD7_3S6
#13289 features across 16347 samples within 1 assay

SD17_3S6.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
SD17_3S6 <- CreateSeuratObject(counts = SD17_3S6.data, project = "SD17_3S6", min.cells = 50, min.features = 500)
SD17_3S6@meta.data[,"Sample"]<-"SD17_3S6"
SD17_3S6
#13064 features across 11769 samples within 1 assay 

SD38_3S6.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
SD38_3S6 <- CreateSeuratObject(counts = SD38_3S6.data, project = "SD38_3S6", min.cells = 50, min.features = 500)
SD38_3S6@meta.data[,"Sample"]<-"SD38_3S6"
SD38_3S6
#13552 features across 11216 samples within 1 assay 

SD147_CS11.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
SD147_CS11 <- CreateSeuratObject(counts = SD147_CS11.data, project = "SD147_CS11", min.cells = 50, min.features = 500)
SD147_CS11@meta.data[,"Sample"]<-"SD147_CS11"
SD147_CS11
#13140 features across 12136 samples within 1 assay  

SD277_CS11.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
SD277_CS11 <- CreateSeuratObject(counts = SD277_CS11.data, project = "SD277_CS11", min.cells = 50, min.features = 500)
SD277_CS11@meta.data[,"Sample"]<-"SD277_CS11"
SD277_CS11
#13646 features across 11205 samples within 1 assay

SD300_SWCB.data <- Read10X(data.dir = "path/to/data./filtered_feature_bc_matrix")
SD300_SWCB <- CreateSeuratObject(counts = SD300_SWCB.data, project = "SD300_SWCB", min.cells = 50, min.features = 500)
SD300_SWCB@meta.data[,"Sample"]<-"SD300_SWCB"
SD300_SWCB
#13628 features across 7441 samples within 1 assay

Combined <- merge(P5, c(P7,MC,SD7_3S6,SD17_3S6,SD38_3S6,SD147_CS11,SD277_CS11,SD300_SWCB), add.cell.ids = c("P5","P7","MC", "SD7_3S6", "SD17_3S6","SD38_3S6","SD147_CS11","SD277_CS11","SD300_SWCB"), project = "Ding2")

Combined <-JoinLayers(Combined)
Combined   #17093 features across 109589 samples within 1 assay

saveRDS(Combined, file = "Ding2.suspention-1.merge9.rds")


## 03_quality_control

Combined[["percent.mt"]] <- PercentageFeatureSet(Combined, pattern = "^MT-")
head(Combined@meta.data, 5)
VlnPlot(Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
Combined.filter <- subset(Combined, subset = nFeature_RNA > 500 & nCount_RNA < 80000 & percent.mt < 25)  # percent.mt none  #  select 700 generates too low median features per cell
Combined.filter 
#17093 features across 108960 samples within 1 assay


##04_normalization_and_scaling
 
Combined.filter <- SCTransform(Combined.filter, vars.to.regress = "percent.mt", verbose = FALSE)
Combined.filter <- RunPCA(Combined.filter, features = VariableFeatures(object = Combined.filter), npcs = 100)
DimPlot(Combined.filter, reduction = "pca")
ElbowPlot(Combined.filter, ndims = 100, reduction = "pca")  # select 30

saveRDS(Combined.filter, file = "Ding2.suspention-3.scTransform.merge9.rds")

new.cluster.ids <- c("P5", "P7", "P9", "SD7", "SD17", "SD38", "SD147", "SD277", "SD300") #transitional
names(new.cluster.ids) <- levels(Combined.filter.Harm)

Combined.filter.Harm <- RenameIdents(Combined.filter.Harm, new.cluster.ids)
levels(Combined.filter.Harm)


Combined.filter.Harm[["Final.id"]] <- Idents(object = Combined.filter.Harm)

library(harmony)
Combined.filter.Harm <- RunHarmony(object = Combined.filter,group.by.vars="new.cluster.ids")
Combined.filter.Harm <- FindNeighbors(Combined.filter.Harm, reduction = "harmony",dims = 1:50)  
Combined.filter.Harm <- FindClusters(Combined.filter.Harm, resolution = 0.04)
Combined.filter.Harm <- RunUMAP(Combined.filter.Harm, reduction = "harmony", reduction.name="umap.harm", dims = 1:50) 

## 05_dimensionality_reduction_and_clustering

DimPlot(Combined.filter.Harm, reduction = "umap.harm", label = TRUE, pt.size = 0.1, label.size = 3.5, raster=FALSE) + NoLegend()+FontSize(x.title = 8, y.title = 8)+theme(axis.text = element_text(size = 8)) + theme(strip.text = element_text(size=8, face="bold"))
#reduction = Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca

DimPlot(Combined.filter.Harm, reduction = "umap.harm",split.by = "Sample", label = TRUE, pt.size = 0.1, label.size = 3.5, raster=FALSE) + NoLegend()+FontSize(x.title = 8, y.title = 8)+theme(axis.text = element_text(size = 8)) + theme(strip.text = element_text(size=8, face="bold"))
#reduction = Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca



## 06_marker_analysis
Combined.filter.Harm.clusterMarkers.POS <- FindAllMarkers(Combined.filter.Harm,only.pos = TRUE,min.pct = 0.1)
write.csv(Combined.filter.Harm.clusterMarkers.POS,"Combined.filter.Harm.clusterMarkers.POS.csv")

## 07 cell_state_annotation
levels(Combined.filter.Harm)
new.cluster.ids <- c("Adip","Mes","ECM-R","Myo") 
names(new.cluster.ids) <- levels(Combined.filter.Harm)
Combined.filter.Harm <- RenameIdents(Combined.filter.Harm, new.cluster.ids)

levels(Combined.filter.Harm)

Combined.filter.Harm[["NF.id"]] <- Idents(object = Combined.filter.Harm)

DimPlot(Combined.filter.Harm, reduction = "umap.harm", label = TRUE, pt.size = 0.1, label.size = 3.5,raster=FALSE) + NoLegend()+FontSize(x.title = 8, y.title = 8)+theme(axis.text = element_text(size = 8)) + theme(strip.text = element_text(size=8, face="bold"))

ggsave("nf.Combined.filter.Harm.UMAP.seuratclusters.tiff", dpi=300, width=2.5, height=2.2)


levels(Combined.filter.Harm)
Idents(Combined.filter.Harm) <-"Cluster.id"
levels(Combined.filter.Harm)
features <- c("DES","MYL1","FLNC","THY1","ENG","DPT","MMP1","MMP9","LUM","ACSS3","PPARG","CD36")
#(ENG/CD105, thy1/cd86)
levels(Combined.filter.Harm) <- c("Myo","Mes","ECM-R","Adip")
DotPlot(Combined.filter.Harm, features = features, dot.scale = 3,cols = c("grey", "red3"))+ RotatedAxis()+theme(axis.text.x = element_text(size = 8))+theme(axis.text.y = element_text(size = 8, face="italic"))+theme(legend.text = element_text(size = 8))+theme(legend.title = element_text(size = 8))+theme(axis.title.x = element_text(size = 10))+theme(axis.title.y = element_text(size = 10))+coord_flip()+theme(legend.key.size = unit(0.2,"cm"))

DotPlot(Combined.filter.Harm, features = features, dot.scale = 3,cols = c("grey", "red3"))+ RotatedAxis()+theme(axis.text.x = element_text(size = 8, face="italic"))+theme(axis.text.y = element_text(size = 8))+theme(legend.text = element_text(size = 8))+theme(legend.title = element_text(size = 8))+theme(axis.title.x = element_text(size = 10))+theme(axis.title.y = element_text(size = 10))+theme(legend.key.size = unit(0.2,"cm"))

ggsave("Combined.filter.Harm.AllCluster.dotplot.tiff", dpi=300, width=4.6, height=2.1)

## 08_pathway_analysis

library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(enrichplot)
library(org.Hs.eg.db)  

library(stringr)
set.seed(1234)
# Adip
AdiDEs <- read.csv("...path_to_file/Adip.csv")

#geneID.TOP.up <- subset(AllDEs, p_val_adj < 0.05 & avg_log2FC > 1.0)
geneID.TOP.up <- AdiDEs$gene[1:150]
geneID.up <- as.character(geneID.TOP.up)

GOPsDEs.BP <- enrichGO(geneID.up,   
                    OrgDb         = 'org.Hs.eg.db',
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    keyType = 'SYMBOL',
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,              
                    readable      = TRUE)

dim(GOPsDEs.BP)
head(GOPsDEs.BP)
#16  12


write.csv(GOPsDEs.BP, "Adip.marker150.GOPsDEs.BP.csv")

a <- c("fatty acid metabolic process","fatty acid derivative biosynthetic process","regulation of lipid localization","lipid import into cell","regulation of fatty acid oxidation")

p1 <- dotplot(GOPsDEs.BP, showCategory = a, font.size=8, x = "GeneRatio")+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color="black", size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(legend.text = element_text(size = 6))+theme(legend.title = element_text(size = 6))+theme(plot.title = element_text(size=8,face = "bold"))+theme(axis.title = element_text(size = 8))+theme(axis.text.x = element_text(angle = 40, size=8))+theme(plot.title = element_text(size=8,face = "bold"))+scale_size_area(max_size = 4)+theme(legend.key.width = unit(0.1,"cm"))+theme(legend.key.height = unit(0.1,"cm"))+scale_y_discrete(labels=function(y) str_wrap(y, width=42))

p1

ggsave('Adip.marker150.GOPsDEs.BP.tiff', dpi=300, width=3.5, height=2.0)





## 09_cell cycle_analysis

cellcycle <- read.csv("..path to file/CellCycle.genes.x.csv")
s.genes <- cellcycle$S
g2m.genes <- cellcycle$G2
Combined.filter.Harm <- CellCycleScoring(Combined.filter.Harm, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Combined.filter.Harm <- RunPCA(Combined.filter.Harm, features = c(s.genes, g2m.genes))


DimPlot(Combined.filter.Harm)+theme(axis.text = element_text(size = 8)) + theme(strip.text = element_text(size=8, face="bold"))+FontSize(x.title = 8, y.title = 8)+theme(axis.text = element_text(size = 8))+theme(legend.text = element_text(size=8))

ggsave("Combined.filter.Harm.cellcycle.all.tiff", dpi=300, width=3.0, height=2.2)

DimPlot(Combined.filter.Harm, split.by = "Sample")+theme(axis.text = element_text(size = 8)) + theme(strip.text = element_text(size=8, face="bold"))+FontSize(x.title = 8, y.title = 8)+theme(axis.text = element_text(size = 8))+theme(legend.text = element_text(size=8))
ggsave("Combined.filter.Harm.cellcycle.Combined.id.tiff", dpi=300, width=10, height=2.5)

cell.ratio <-table(Combined.filter.Harm@meta.data$Sample, Combined.filter.Harm$Phase)
print(cell.ratio)
write.csv(cell.ratio,"cellCycle.cell.ratio.csv")

cell.ratio1 <-table(Combined.filter.Harm@meta.data$Cluster.id, Combined.filter.Harm$Phase)
print(cell.ratio1)
write.csv(cell.ratio1,"cellCycle.celltype.ratio.csv")



## 10_pseudotime_analysis

library(monocle3)
library(SeuratWrappers)
set.seed(1234)

Idents(Combined.filter.Harm) <-"orig.ident"
levels(Combined.filter.Harm)

D7 <- subset(Combined.filter.Harm, idents="SD7_3S6")
levels(D7)

Idents(D7) <-"Cluster.id"
DimPlot(D7, reduction = "umap.harm", label = TRUE, pt.size = 0.05, label.size = 2)+ NoLegend()+FontSize(x.title = 8, y.title = 8, main = 8)+theme(axis.text = element_text(size = 8)) + theme(strip.text = element_text(size=8, face="bold"))
ggsave("SD7_3S6.subPopulation.tiff",dpi=300, width = 2.0, height = 1.95)

D7_Diet <- DietSeurat(D7, counts=TRUE,data=TRUE, dimreducs = "umap.harm")

D7_Diet[["UMAP"]] <- D7_Diet[["umap.harm"]]  ## If THE UMAP from seurat needs to be kept, it only accepts "UMAP' or "PCA" and others. As we want to use the WNN umap, we need to name it as "UMAP" to be processed in Monocle3.

D7_Diet.cds <- as.cell_data_set(D7_Diet)
D7_Diet.cds <- cluster_cells(cd=D7_Diet.cds, reduction_method = "UMAP")
D7_Diet.cds <- learn_graph(D7_Diet.cds,use_partition = TRUE)

D7_Diet.cds <- order_cells(D7_Diet.cds, reduction_method = "UMAP")


plot_cells(cds = D7_Diet.cds, label_cell_groups = FALSE, color_cells_by = "ident", show_trajectory_graph = TRUE, label_branch_points = FALSE, label_leaves = FALSE, label_roots = TRUE, cell_size = 0.5)+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.55))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.55))+theme(axis.ticks.y = element_line(size = 0.55))+theme(legend.title = element_text(size = 8))+theme(legend.text = element_text(size = 8))

ggsave("D7_Diet.PSEDO.tiff", dpi=300, width=3.0, height=1.8)


plot_cells(cds = D7_Diet.cds, label_cell_groups = FALSE, color_cells_by = "pseudotime", show_trajectory_graph = TRUE, label_branch_points = FALSE, label_leaves = FALSE, label_roots = TRUE, cell_size = 0.5)+theme(axis.title = element_text(size=8))+theme(axis.text = element_text(color='black', size=8))+theme(legend.title = element_text(size=8)+theme(legend.text = element_text(size = 8)))+theme(axis.line.x = element_line(size=0.55))+theme(axis.line.y = element_line(size=0.55))+theme(axis.ticks.x = element_line(size=0.55))+theme(axis.ticks.y = element_line(size = 0.55))+theme(legend.title = element_text(size = 8))+theme(legend.text = element_text(size = 8))
ggsave("D7_Diet.PSEDO2.tiff", dpi=300, width=2.5, height=1.8)


