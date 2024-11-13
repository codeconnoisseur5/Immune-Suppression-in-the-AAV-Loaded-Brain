# Immune-Suppression-in-the-AAV-Loaded-Brain
library(Rcpp)
library(harmony)

logFCfilter=1           
adjPvalFilter=0.05      

workDir="。。。。。"
setwd(workDir)

dirs=list.dirs(workDir)
dirs_sample=dirs[-1]
names(dirs_sample)=gsub(".+\\/(.+)", "\\1", dirs_sample)
counts <- Read10X(data.dir = dirs_sample)
pbmc = CreateSeuratObject(counts, min.cells=3, min.features=100)

 c[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#??? (file="01.featureViolin.pdf", width=10, height=6.5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 15)    #??? ?? (file="01.featureCor.pdf", width=13, height=7)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#??? c <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=1500)
#??? 10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf", width=10, height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


pbmc=RunHarmony(pbmc, "orig.ident")

#????ÿ??P "02.pcaGene.pdf", width=10, height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures=20)
dev.off()

#????PCA? "02.PCA.pdf", width=7.5, height=5)
DimPlot(object=pbmc, reduction="pca")
dev.off()

#PCA????? "02.pcaHeatmap.pdf", width=10, height=8)
DimHeatmap(object=pbmc, dims=1:4, cells=500, balanced=TRUE, nfeatures=30, ncol=2)
dev.off()

#?õ?ÿ??P ackStraw(object=pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object=pbmc, dims=1:20)
dev.off()

pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)     #?????ڽ??? indClusters(pbmc, resolution=seq(0.5, 1.2, by=0.1))
pbmc <- FindClusters(object = pbmc, resolution=0.6)
#????????? "03.cluster.pdf", width=7, height=6)
#pbmc <-RunUMAP(object = pbmc, dims = 1:pcSelect)        #UMAP????
pbmc, reduction = "umap", pt.size = 2, label = TRUE)   #UMAP????
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             object = pbmc, pt.size = 2, label = TRUE)     #TSNE???ӻ
write.table(pbmc$seurat_clusters,file="03.Cluster.txt",quote=F,sep="\t",col.names=F)

##????ÿ? ers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#????marke "03.clusterHeatmap.pdf",width=15, height=15)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()



######### SingleR <- GetAssayData(pbmc, layer="data")
clusters<-pbmc@meta.data$seurat_clusters
ref1=get(load("ref_Human_all.RData"))
ref2=get(load("ref_Hematopoietic.RData"))
ref3=get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4=get(load("ImmGenData.Rdata"))
singler=SingleR(test=pbmc_for_SingleR, ref =list(ref1, ref2, ref3, ref4),
                labels=list(ref1$label.main,ref2$label.main,ref3$label.main,ref4$label.main), clusters = clusters)

clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
clusterAnn$labels=gsub("_", " ", clusterAnn$labels)
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
#????ϸ??? ()
for(i in 1:length(pbmc$seurat_clusters)){
	index=pbmc$seurat_clusters[i]
	cellAnn=c(cellAnn, clusterAnn[index,2])
}
cellAnnOut=cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut)=c("id", "labels")
write.table(cellAnnOut, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#ϸ??ע?ͺ =gsub("_", " ", singler$labels)
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.cellAnn.pdf", width=7.5, height=6)
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)     #UMAP???object = pbmc, pt.size = 2, label = TRUE)                #TSNE????
#???????ӻ ("(.*?)\\..*", "\\1", colnames(pbmc))
names(Type)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=Type, col.name="Type")
pdf(file="04.group.cellAnn.pdf", width=11, height=6)
#DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE, split.by="Type")       #UMAP???ӻ object = pbmc, pt.size = 1, label = TRUE, split.by="Type")     #TSNE???ӻ

#ϸ??????? ub("(.*?)\\..*", "\\1", colnames(pbmc))
groups=paste0(groups, "_", cellAnn)
names(groups)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
for(cellName in unique(cellAnn)){
	conName=paste0("Control_", cellName)
	treatName=paste0("Treat_", cellName)
	pbmc.markers=FindMarkers(pbmc, ident.1=treatName, ident.2=conName, group.by='group', logfc.threshold=0.1)
	sig.markersGroup=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
	sig.markersGroup=cbind(Gene=row.names(sig.markersGroup), sig.markersGroup)
	write.table(sig.markersGroup,file=paste0("05.", cellName, ".diffGene.txt"),sep="\t",row.names=F,quote=F)
}

#???浥ϸ? , cellAnn, file="Seurat.Rdata")


######### atrix=as.matrix(pbmc@assays$RNA$data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name=row.names(monocle.matrix), row.names=row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

#??Seurat? s(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#????ϸ??? n=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

cds <- estgc()
imateDispersions(cds)
cds <- setgc()
OrderingFilter(cds, as.vector(sig.markers$gene))
cds <- redgc()
uceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- ordgc()
erCells(cds)
#????ϸ??? "06.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()
#?????ֻ?? "06.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#????ϸ??? "06.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#????????? "06.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()

library(usethis)
library(devtools)
library(limma)
library(registry)
library(rngtools)
library(cluster)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(dplyr)
library(igraph)
library(CellChat)
setwd("/")

 
load("Seurat.Rdata")
expMatrix=as.matrix(pbmc@assays$RNA$data)
meta=as.data.frame(cellAnn)
colnames(meta)[1]="labels"
row.names(meta)=names(pbmc$seurat_clusters)

 
cellchat <- createCellChat(object = expMatrix, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use="labels")
groupSize <- as.numeric(table(cellchat@idents))      

 
CellChatDB <- CellChatDB.human       
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

 
pdf(file="COMM01.DatabaseCategory.pdf", width=7, height=5)
showDatabaseCategory(CellChatDB)
dev.off()

 
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)       
cellchat <- identifyOverExpressedInteractions(cellchat)      
cellchat <- projectData(cellchat, PPI.human)  

 
cellchat <- computeCommunProb(cellchat)
 
cellchat <- filterCommunication(cellchat, min.cells = 10)
 
df.net=subsetCommunication(cellchat)
write.table(file="COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)

 
cellchat <- computeCommunProbPathway(cellchat)

 
cellchat <- aggregateNet(cellchat)
 
pdf(file="COMM03.cellNetworkCount.pdf", width=7, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
 
pdf(file="COMM04.cellNetworkWeight.pdf", width=7, height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()

 (file="COMM05.singleCell.pdf", width=9, height=7.5)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
}
dev.off()

#??? l in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  pdf(file=paste0("COMM05.", cel, ".pdf"), width=6.5, height=5.5)
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
  dev.off()
}

#????ϸ? e="COMM06.bubble.pdf", width=9.5, height=6)
netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x=45)
dev.off()


