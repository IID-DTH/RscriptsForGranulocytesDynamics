library(Seurat)
library(magrittr)
library(dplyr)
library(monocle2)

nlcl_group=data.frame(row.names = colnames(nlcl))
nlcl_group$group="Naive"
nlcl_group[grep("-2",rownames(nlcl_group)),"group"]="LPS"
nlcl_group[grep("-3",rownames(nlcl_group)),"group"]="CLP"
nlcl_group[grep("-4",rownames(nlcl_group)),"group"]="CLP-LPS"

nlcl_filtered <- Read10X(data.dir = "CLP_LPS_aggr/outs/filtered_gene_bc_matrices_mex/mm10/")
nlcl <- CreateSeuratObject(counts = nlcl_filtered, project = "nlcl", min.cells = 3, min.features = 200, meta.data = nlcl_group)
nlcl[["percent.mt"]] <- PercentageFeatureSet(nlcl, pattern = "^mt-")
nlcl=subset(nlcl,subset=nFeature_RNA >200 & percent.mt < 5)

nlcl.list <- SplitObject(nlcl, split.by = "group")
for (i in 1:length(nlcl.list)) {
    nlcl.list[[i]] <- NormalizeData(nlcl.list[[i]], verbose = FALSE)
    nlcl.list[[i]] <- FindVariableFeatures(nlcl.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}
nlcl.anchors <- FindIntegrationAnchors(object.list = nlcl.list, dims = 1:30)
nlcl.integrated <- IntegrateData(anchorset = nlcl.anchors, dims = 1:30)
DefaultAssay(nlcl.integrated) <- "integrated"
nlcl.integrated <- ScaleData(nlcl.integrated, verbose = FALSE)
nlcl.integrated <- RunPCA(nlcl.integrated, npcs = 30, verbose = FALSE)w
nlcl.integrated <- RunTSNE(nlcl.integrated, reduction = "pca", dims = 1:30)
#ElbowPlot(nlcl.integrated)

nlcl.integrated <- FindNeighbors(nlcl.integrated, dims = 1:15)
nlcl.integrated <- FindClusters(nlcl.integrated)
nlcl.integrated$group=factor(nlcl.integrated$group,levels=c("Naive","LPS","CLP","CLP-LPS"))
nlcl_clust=data.frame(cluster=Idents(nlcl.integrated))
nlcl_clust$group=nlcl.integrated$group
nlcl_clust$celltype="Granulocyte"
nlcl_clust[nlcl_clust$cluster %in% c(8,10,17),]$celltype = "Monocyte"
nlcl_clust[nlcl_clust$cluster %in% c(14,15),"celltype"] = "Erthrocyte"
nlcl_clust[nlcl_clust$cluster ==16 ,"celltype"] = "T cell"
nlcl.integrated$celltype=nlcl_clust$celltype

nlcl_granu=subset(nlcl.integrated,celltype=="Granulocyte")
nlcl_granu <- FindClusters(nlcl_granu,resolution = 0.2)
nlcl.markers_0.2 <- FindAllMarkers(nlcl_granu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nlcl_granu$cluster_0.2=factor(nlcl_granu$integrated_snn_res.0.2,levels = c(4,3,1,0,2),labels = c("G1","G2","G3","G4","G5"))
nlcl.markers_0.2=nlcl.markers_0.2[order(nlcl.markers_0.2$cluster_0.2),]
top10_0.2 <- nlcl.markers_0.2 %>% group_by(clusterG) %>% top_n(n = 5, wt = avg_logFC)

#Figure 3B
 DoHeatmap(nlcl_granu, features = top10_0.2$gene,group.by="cluster_0.2") + NoLegend()
 FeaturePlot(nlcl_granu, features = c("Elane", "Rps2", "Ube2c","Top2a","Cebpe","Chil3","Mmp8","Retnlg","Wfdc17","Ifitm1"),min.cutoff = 0) 

select_gene=c("Elane", "Top2a","Cebpe","Mmp8","Ccl6","Ly6g","Itgam")

gene_nlcl<-load_cellranger_matrix("./CLP_LPS_aggr")
fd_nlcl <- fData(gene_nlcl)
colnames(fd_nlcl)[2] <- "gene_short_name"
cds_nlcl <- newCellDataSet(exprs(gene_nlcl),
                           phenoData = new("AnnotatedDataFrame", data = pData(gene_nlcl)),
                           featureData = new("AnnotatedDataFrame", data = fd_nlcl),
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
cds_nlcl <- estimateSizeFactors(cds_nlcl)
cds_nlcl <- estimateDispersions(cds_nlcl)

cds_nlcl_granu <- cds_nlcl[,pData(cds_nlcl)$barcode %in% granu_cluster$id ]

pData(cds_nlcl_granu)=cbind(pData(cds_nlcl_granu),granu_cluster[rownames(pData(cds_nlcl_granu)),c("cluster","group")])
colnames(pData(cds_nlcl_granu))[3]="cluster_seurat"
cds_nlcl_granu <- detectGenes(cds_nlcl_granu, min_expr = 0.1)
fData(cds_nlcl_granu)$use_for_ordering <- fData(cds_nlcl_granu)$num_cells_expressed > 0.05 * ncol(cds_nlcl_granu)

granu_diff<-read.delim("granu_diff.txt",stringsAsFactors = F)
granu_diff<- merge(granu_diff,fData(cds_nlcl_granu)[,c("id","gene_short_name","num_cells_expressed")],by.x="gene",by.y="gene_short_name",all.x=T)
cds_nlcl_granu_ps=cds_nlcl_granu
fData(cds_nlcl_granu_ps)$use_for_ordering <- fData(cds_nlcl_granu_ps)$gene_short_name %in% subset(granu_diff,num_cells_expressed>10)$gene
cds_nlcl_granu_ps <- reduceDimension(cds_nlcl_granu_ps,method = 'DDRTree',cores=8)
cds_nlcl_granu_ps <- orderCells(cds_nlcl_granu_ps,reverse = T)
plot_cell_trajectory(cds_nlcl_granu_ps,color_by = "cluster_seurat",show_branch_points=F)+facet_wrap(~group)

#Figure 3C
ggplot(pData(cds_nlcl_granu_ps),aes(x=Pseudotime,fill=cluster_seurat))+geom_histogram(bins = 100)+theme_classic()+facet_grid(group~.)

granu_nl_diff=data.frame()
for (i in 1:nlevels(granu_cluster_0.2$cluster)){
    tmp_nl = differentialGeneTest(cds_nlcl_granu[fData(cds_nlcl_granu)$num_cells_expressed >= 10, subset(granu_cluster_0.2,cluster == levels(granu_cluster_0.2$cluster)[i] & group %in% c("Naive", "LPS"))$id],fullModelFormulaStr = '~group', cores = 8)
    tmp_nl$cluster = levels(granu_cluster_0.2$cluster)[i]
    #tmp_nl=subset(tmp_nl,qval<0.01)
    granu_nl_diff=rbind(granu_nl_diff,tmp_nl)
}

granu_cl_diff=data.frame()
for (i in 1:nlevels(granu_cluster_0.2$cluster)){
    tmp_nl = differentialGeneTest(cds_nlcl_granu[fData(cds_nlcl_granu)$num_cells_expressed >= 10, subset(granu_cluster_0.2,cluster == levels(granu_cluster_0.2$cluster)[i] & group %in% c("CLP", "CLP-LPS"))$id],fullModelFormulaStr = '~group', cores = 8)
    tmp_nl$cluster = levels(granu_cluster_0.2$cluster)[i]
#    tmp_nl=subset(tmp_nl,qval<0.001)
    granu_cl_diff=rbind(granu_cl_diff,tmp_nl)
}

nl_group_exp_logFC=data.frame(row.names = rownames(nlcl_cluster_group_exp2),
    G1=log2(nlcl_cluster_group_exp2$G1_LPS+0.1)-log2(nlcl_cluster_group_exp2$G1_Naive+0.1),
    G2=log2(nlcl_cluster_group_exp2$G2_LPS+0.1)-log2(nlcl_cluster_group_exp2$G2_Naive+0.1),
    G3=log2(nlcl_cluster_group_exp2$G3_LPS+0.1)-log2(nlcl_cluster_group_exp2$G3_Naive+0.1),
    G4=log2(nlcl_cluster_group_exp2$G4_LPS+0.1)-log2(nlcl_cluster_group_exp2$G4_Naive+0.1),
    G5=log2(nlcl_cluster_group_exp2$G5_LPS+0.1)-log2(nlcl_cluster_group_exp2$G5_Naive+0.1))
cl_group_exp_logFC=data.frame(row.names = rownames(nlcl_cluster_group_exp2),
    G1=log2(nlcl_cluster_group_exp2$G1_CLP.LPS+0.1)-log2(nlcl_cluster_group_exp2$G1_CLP+0.1),
    G2=log2(nlcl_cluster_group_exp2$G2_CLP.LPS+0.1)-log2(nlcl_cluster_group_exp2$G2_CLP+0.1),
    G3=log2(nlcl_cluster_group_exp2$G3_CLP.LPS+0.1)-log2(nlcl_cluster_group_exp2$G3_CLP+0.1),
    G4=log2(nlcl_cluster_group_exp2$G4_CLP.LPS+0.1)-log2(nlcl_cluster_group_exp2$G4_CLP+0.1),
    G5=log2(nlcl_cluster_group_exp2$G5_CLP.LPS+0.1)-log2(nlcl_cluster_group_exp2$G5_CLP+0.1))
mt_nl_logFC<-melt(as.matrix(nl_group_exp_logFC))
colnames(mt_nl_logFC) =c("id","cluster","logFC")
mt_nl_logFC=merge(granu_nl_diff,mt_nl_logFC)
mt_nl_logFC$group="non"
mt_nl_logFC[mt_nl_logFC$pval<0.05&mt_nl_logFC$logFC>1,"group"]="Up"
mt_nl_logFC[mt_nl_logFC$pval<0.05&mt_nl_logFC$logFC< -1,"group"]="Down"

mt_cl_logFC<-melt(as.matrix(cl_group_exp_logFC))
colnames(mt_cl_logFC) =c("id","cluster","logFC")
mt_cl_logFC=merge(granu_cl_diff,mt_cl_logFC)
mt_cl_logFC$group="non"
mt_cl_logFC[mt_cl_logFC$pval<0.05&mt_cl_logFC$logFC>1,"group"]="Up"
mt_cl_logFC[mt_cl_logFC$pval<0.05&mt_cl_logFC$logFC< -1,"group"]="Down"

nl_sig=unique(subset(mt_nl_logFC,num_cells_expressed>100&group=="Up")$id)
cl_sig=unique(subset(mt_cl_logFC,num_cells_expressed>100&group=="Up")$id)
cds_nlcl_granu_lps<-setOrderingFilter(cds_nlcl_granu,union(nl_sig,cl_sig))
cds_nlcl_granu_lps<-setOrderingFilter(cds_nlcl_granu,rownames(list_lps_response))
cds_nlcl_granu_lps <- reduceDimension(cds_nlcl_granu_lps,  reduction_method = 'DDRTree',cores=8)
cds_nlcl_granu_lps <-orderCells(cds_nlcl_granu_lps,reverse = T)
plot_cell_trajectory(cds_nlcl_granu_lps,color_by = "cluster_0.2",show_branch_points = F)+facet_wrap(~group)

pseudotime_lps=(pData(cds_nlcl_granu_lps))[,c(1,4,13,14,15)]
colnames(pseudotime_lps)[5]="pseudotime_lps"
pseudotime_diff=pData(cds_nlcl_granu_ps)[,c(1,4,12,14)]
colnames(pseudotime_diff)[3]="pseudotime_diff"
 pseudotime_lps=merge(pseudotime_lps,pseudotime_diff)
 ggplot(pseudotime_lps,aes(x=pseudotime_diff,y=pseudotime_lps,col=cluster_0.2))+geom_point()+facet_wrap(~group)+theme_classic()+geom_hline(yintercept = c(7.0,10.5))
ggplot(pseudotime_lps,aes(x=pseudotime_lps,fill=cluster_0.2))+geom_histogram(bins = 100)+facet_grid(group~.)+theme_classic()+geom_vline(xintercept = c(7.0,10.5))
