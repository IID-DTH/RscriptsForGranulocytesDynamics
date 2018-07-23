
library(cellrangerRkit)
library(ggplot2)
library(reshape2)
library(pheatmap)

gene_nl<-load_cellranger_matrix("./naive_lps")
gene_nl_res<-load_cellranger_analysis_results("./naive_lps")
gene_naive<-load_cellranger_matrix("./naive/")
gene_naive_res<-load_cellranger_analysis_results("./naive/")
gene_lps<-load_cellranger_matrix("./lps/")
gene_lps_res<-load_cellranger_analysis_results("./lps/")

gene_exp_nl<-exprs(gene_nl)
gene_exp_naive<-exprs(gene_naive)
gene_exp_lps<-exprs(gene_lps)

##Data processing in Naive cells
raw_naive_exp<-read.delim("naive/outs/raw_gene_bc_matrices/mm10/matrix.mtx",sep=" ",skip = 1,header=F)
raw_naive_exp<-raw_naive_exp[3:nrow(raw_naive_exp),]
 colnames(raw_naive_exp)<-c("geneID","BarcodeID","UMICount")
 raw_naive_exp$geneID<-factor(raw_naive_exp$geneID)
 raw_naive_exp$BarcodeID<-factor(raw_naive_exp$BarcodeID)

gene_ids<-read.delim("naive/outs/raw_gene_bc_matrices/mm10/genes.tsv",header=F)
colnames(gene_ids) <- c("ENSG_id","gene_symbol")
  
naive_UMI_cell<-data.frame(UMI=table(tapply(raw_naive_exp$UMICount,raw_naive_exp$BarcodeID,sum)))
colnames(naive_UMI_cell)<-c("UMI_per_cell","CellCount")
naive_UMI_cell$UMI_per_cell<-as.numeric(as.character(naive_UMI_cell$UMI_per_cell))
  
naive_UMI_cell$group<-"Backgroud"
naive_UMI_cell[naive_UMI_cell$UMI_per_cell>=naive_UMI_cell[30,1]*0.1,"group"]<-"Cells"
naive_UMI_cell[naive_UMI_cell$UMI_per_cell>30000,"group"]<-"Burn-in"

naive_gene_count<-data.frame(geneCount=apply(gene_exp_naive,2,function(a){length(a[a>0])}),UMICount=apply(gene_exp_naive,2,sum))
naive_gene_count<-naive_gene_count[order(naive_gene_count$UMICount,decreasing=T),]

use_genes <- get_nonzero_genes(gene_naive) 
gene_naive_bcnorm <- normalize_barcode_sums_to_median(gene_naive[use_genes,])
gene_naive_log <- log_gene_bc_matrix(gene_naive_bcnorm,base=10)

##Cell clustering in Naive cells

gene_naive_tsne<-data.frame(gene_naive_res$tsne,row.names = "Barcode")
gene_naive_tsne_dist<-as.matrix(dist(gene_naive_tsne[rownames(gene_naive_res)[6:nrow(gene_naive_res)]]))
gene_naive_tsne_dist<-max(gene_naive_tsne_dist)-gene_naive_tsne_dist
naive_hc<-hclust(dist(gene_naive_tsne_dist))
gsHc <- clusGap(gene_naive_tsne, FUN = hclusCut, K.max = 20, B = 60)
plot(gsHc)
# Result shows that the best cluster number is 9
naive_cluster_tsne<-data.frame(clust=cutree(naive_hc,k=9))
naive_cluster_tsne$clust<-factor(naive_cluster_tsne$clust)
naive_cluster_tsne$clust<-factor(naive_cluster_tsne$clust,levels=c(1,6,5,8,3,2,7,4,9),labels=c("NM","NG1","NG2","NG3","NG4","NG5","NG6","NG7","NE"))

naive_cluster_tsne_type<-naive_cluster_tsne
naive_cluster_tsne_type$CellType<-"Granulocyte"
naive_cluster_tsne_type[naive_cluster_tsne_type$clust=="N1","CellType"] <- "Monocyte"
naive_cluster_tsne_type[naive_cluster_tsne_type$clust=="N9","CellType"] <- "Megakaryocyte"
naive_cluster_tsne_type<-cbind(naive_cluster_tsne_type,gene_naive_tsne[rownames(naive_cluster_tsne_type),])

 
##Differential expression in Naive clusters
naive_gene_count$clust_hclust<-factor(naive_cluster_tsne[naive_gene_count$Barcode,"clust"])
select_naive_cell<-subset(naive_gene_count,clust_hclust!="NM"&clust_hclust!="NE")
select_naive_exp<-gene_exp_naive[,select_naive_cell]
select_naive_exp<-select_naive_exp[apply(select_naive_exp,1,sum)>0,]

select_naive_exp_tpm<-t(t(select_naive_exp)/apply(select_naive_exp,2,sum)*1000000)
select_naive_exp_tpm_log<-log2(select_naive_exp_tpm+1)
select_naive_gene_tpm_sig<-data.frame(p.value=apply(select_naive_exp_tpm_log,1,function(a){kruskal.test(a~select_naive_cell$clust_hclust)$p.value}))
select_naive_gene_tpm_fdr<-apply(select_naive_gene_tpm_sig,2,function(a){p.adjust(a,method = "fdr")})
select_naive_gene_tpm_fdr<-data.frame(select_naive_gene_tpm_fdr)
select_naive_exp_tpm_mean<-t(apply(select_naive_exp_tpm_log,1,function(a){tapply(a,select_naive_cell$clust_hclust,mean)}))
scale_select<-select_naive_exp_tpm_mean[rownames(subset(select_naive_gene_tpm_fdr,p.value<0.0001)),]

scale_select<-t(scale(t(scale_select)))
naive_gene_clust<-data.frame(clust=(cutree(hclust(dist(scale_select)),h=4.8)))
naive_gene_clust$clust<-factor(naive_gene_clust$clust,labels = c("G1 - 5460 genes","G2 - 571 genes","G3 - 164 genes","G4 - 222 genes"))

select_naive_exp_tpm_pos<-t(apply(select_naive_exp_tpm_log,1,function(b){tapply(b,select_naive_cell$clust_hclust,function(a){length(a[a>0])})}))
select_naive_exp_tpm_pos<-t(t(select_naive_exp_tpm_pos)/c(table(select_naive_cell$clust_hclust)))

##Data Processing and cell clustering in LPS cells
raw_lps_exp<-read.delim("lps/outs/raw_gene_bc_matrices/mm10/matrix.mtx",sep=" ",skip = 1,header=F)
raw_lps_exp<-raw_lps_exp[3:nrow(raw_lps_exp),]
colnames(raw_lps_exp)<-c("geneID","BarcodeID","UMICount")
raw_lps_exp$geneID<-factor(raw_lps_exp$geneID)
raw_lps_exp$BarcodeID<-factor(raw_lps_exp$BarcodeID)
lps_UMI_cell<-data.frame(UMI=table(tapply(raw_lps_exp$UMICount,raw_lps_exp$BarcodeID,sum)))
colnames(lps_UMI_cell)<-c("UMI_per_cell","CellCount")
lps_UMI_cell$UMI_per_cell<-as.numeric(as.character(lps_UMI_cell$UMI_per_cell))
lps_UMI_cell$group<-"Backgroud"
lps_UMI_cell[lps_UMI_cell$UMI_per_cell>=lps_UMI_cell[30,1]*0.1,"group"]<-"Cells"
lps_UMI_cell[lps_UMI_cell$UMI_per_cell>30000,"group"]<-"Burn-in"
lps_gene_count<-data.frame(geneCount=apply(gene_exp_lps,2,function(a){length(a[a>0])}),UMICount=apply(gene_exp_lps,2,sum))
lps_gene_count<-lps_gene_count[order(lps_gene_count$UMICount,decreasing=T),]
use_genes <- get_nonzero_genes(gene_lps) 
gene_lps_bcnorm <- normalize_barcode_sums_to_median(gene_lps[use_genes,])
gene_lps_log <- log_gene_bc_matrix(gene_lps_bcnorm,base=10)
gene_lps_tsne<-data.frame(gene_lps_res$tsne,row.names = "Barcode")
gene_lps_tsne_dist<-as.matrix(dist(gene_lps_tsne[rownames(gene_lps_res)[6:nrow(gene_lps_res)]]))
gene_lps_tsne_dist<-max(gene_lps_tsne_dist)-gene_lps_tsne_dist
lps_hc<-hclust(dist(gene_lps_tsne_dist))
lps_cluster_tsne<-data.frame(clust=cutree(lps_hc,k=9)) #same as Naive clusters
lps_cluster_tsne$clust<-factor(lps_cluster_tsne$clust)
lps_cluster_tsne$clust<-factor(lps_cluster_tsne$clust,levels=c(4,1,2,6,7,9,8,5,3),labels=c("LM","LG1","LG2","LG3","LG4","LG5","LG6","LG7","LE"))

##Differential expression in LPS clusters
lps_gene_count$clust_hclust<-factor(lps_cluster_tsne[lps_gene_count$Barcode,"clust"])
select_lps_cell<-subset(lps_gene_count,clust_hclust!="LM"&clust_hclust!="LE")
select_lps_exp<-gene_exp_lps[,select_lps_cell]
select_lps_exp<-select_lps_exp[apply(select_lps_exp,1,sum)>0,]

select_lps_exp_tpm<-t(t(select_lps_exp)/apply(select_lps_exp,2,sum)*1000000)
select_lps_exp_tpm_log<-log2(select_lps_exp_tpm+1)
select_lps_gene_tpm_sig<-data.frame(p.value=apply(select_lps_exp_tpm_log,1,function(a){kruskal.test(a~select_lps_cell$clust_hclust)$p.value}))
select_lps_gene_tpm_fdr<-apply(select_lps_gene_tpm_sig,2,function(a){p.adjust(a,method = "fdr")})
select_lps_gene_tpm_fdr<-data.frame(select_lps_gene_tpm_fdr)
select_lps_exp_tpm_mean<-t(apply(select_lps_exp_tpm_log,1,function(a){tapply(a,select_lps_cell$clust_hclust,mean)}))
scale_select<-select_lps_exp_tpm_mean[rownames(subset(select_lps_gene_tpm_fdr,p.value<0.0001)),]

scale_select<-t(scale(t(scale_select)))
lps_gene_clust<-data.frame(clust=(cutree(hclust(dist(scale_select)),h=4.8))) #same as Naive clusters
lps_gene_clust$clust<-factor(lps_gene_clust$clust,labels = c("G1 - 4448 genes","G2 - 1654 genes","G3 - 91 genes","G4 - 374 genes"))

select_lps_exp_tpm_pos<-t(apply(select_lps_exp_tpm_log,1,function(b){tapply(b,select_lps_cell$clust_hclust,function(a){length(a[a>0])})}))
select_lps_exp_tpm_pos<-t(t(select_lps_exp_tpm_pos)/c(table(select_lps_cell$clust_hclust)))



##Differential expression between Naive and LPS clusters
gene_nl<-load_cellranger_matrix("./naive_lps")
gene_nl_res<-load_cellranger_analysis_results("./naive_lps")
naive_cluster<-read.delim("naive_cell_clust.tsv",sep=" ")
lps_cluster<-read.delim("lps_cell_clust.tsv",sep=" ")
rownames(lps_cluster)<-sub("-1","-2",rownames(lps_cluster))
cell_clust_sep<-rbind(naive_cluster,lps_cluster)
cell_clust_tsne<-gene_nl_res$tsne 
cell_clust_tsne<-subset(cell_clust_tsne,Barcode %in% rownames(cell_clust_sep))
cell_clust_tsne$Barcode<-as.character(cell_clust_tsne$Barcode)
cell_clust_tsne$clust_raw<-cell_clust_sep[cell_clust_tsne$Barcode,]

select_nl_exp_pos_count<-data.frame(t(apply(select_nl_exp,1,function(b){tapply(b[cell_clust_tsne$Barcode], cell_clust_tsne$clust_raw, function(a){length(a[a>0])})})))
select_nl_exp_neg_count<-data.frame(t(apply(select_nl_exp,1,function(b){tapply(b[cell_clust_tsne$Barcode], cell_clust_tsne$clust_raw, function(a){length(a[a==0])})})))
select_nl_exp_pos_perc<-select_nl_exp_pos_count/(select_nl_exp_pos_count+select_nl_exp_neg_count)

diff_exp_pos_test<-function(clust1,clust2){ 
select_nl_exp_pos_to_test<-cbind(select_nl_exp_pos_count[,c(clust1,clust2)],select_nl_exp_neg_count[,c(clust1,clust2)])
colnames(select_nl_exp_pos_to_test)<-c("pos_c1","pos_c2","neg_c1","neg_c2")
select_nl_exp_pos_to_test<-subset(select_nl_exp_pos_to_test,pos_c1>0&pos_c2>0)
select_nl_exp_pos_diff<-data.frame(t(apply(select_nl_exp_pos_to_test,1,function(a){f_res<-fisher.test(matrix(a,2));return(c(p.value=f_res$p.value,odds.ratio=as.numeric(f_res$estimate)))})))
select_nl_exp_pos_diff$q.value<-p.adjust(select_nl_exp_pos_diff$p.value,method = "fdr")
select_nl_exp_pos_diff$pattern="common"
select_nl_exp_pos_diff[select_nl_exp_pos_diff$q.value<0.05&select_nl_exp_pos_diff$odds.ratio>4,"pattern"]="LPS_down"
select_nl_exp_pos_diff[select_nl_exp_pos_diff$q.value<0.05&select_nl_exp_pos_diff$odds.ratio<0.25,"pattern"]="LPS_up"
select_nl_exp_pos_diff$clust=paste(clust1,clust2,sep="_")
select_nl_exp_pos_diff$gene=rownames(select_nl_exp_pos_diff)
return(select_nl_exp_pos_diff)
}

ng_group<-c("NG1","NG2","NG3","NG4","NG5","NG6","NG7")
lg_group<-c("LG1","LG2","LG3","LG4","LG5","LG6","LG7")

diff_exp_pos_res_matrix<-data.frame()

for (i in 1:7){
  for (j in 1:7){
    diff_exp_pos_res_matrix <- rbind (diff_exp_pos_res_matrix,diff_exp_pos_test(ng_group[i],lg_group[j]))
  }
}

diff_exp_matrix_table<-as.data.frame(table(diff_exp_pos_res_matrix[,c("clust","pattern")]))
diff_exp_matrix_table_prop<- prop.table(diff_exp_matrix_table,1)
common_prop_table<-data.frame(common_prop=diff_exp_matrix_table_prop[,1],clust=rownames(diff_exp_matrix_table_prop))
common_prop_table<-cbind(common_prop_table,matrix(unlist(strsplit(as.character(common_prop_table$clust),split="_")),ncol=2,byrow=T))
colnames(common_prop_table) [3:4] <-c("clust1","clust2")

diff_exp_pos_res_matrix2<-data.frame()
for (i in 1:7){

    diff_exp_pos_res_matrix2 <- rbind (diff_exp_pos_res_matrix2,diff_exp_pos_test(ng_group[i],lg_group[i]))
  
}


##pseudo-time analysis
library(monocle)
fd_nl <- fData(gene_nl)
colnames(fd_nl)[2] <- "gene_short_name"
cds_nl <- newCellDataSet(exprs(gene_nl),
          phenoData = new("AnnotatedDataFrame", data = pData(gene_nl)),
          featureData = new("AnnotatedDataFrame", data = fd_nl),
          lowerDetectionLimit = 0.5,
          expressionFamily = negbinomial.size())
cds_nl <- estimateSizeFactors(cds_nl)
cds_nl <- estimateDispersions(cds_nl)

cds_nl <- detectGenes(cds_nl, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds_nl),num_cells_expressed >= 10))
pData(cds_nl)$Total_mRNAs <- Matrix::colSums(exprs(cds_nl))
cds_nl <- cds_nl[,pData(cds_nl)$Total_mRNAs < 1e6]
pData(cds_nl)$clust_raw<-cell_clust_sep[pData(cds_nl)$barcode,"clust"]
cds_granu <- cds_nl[,!pData(cds_nl)$clust_raw %in% c("NM","LM","NE","LE")]
cds_granu <- cds_granu[,!is.na(pData(cds_granu)$clust_raw)]
diff_test_res_clust <- differentialGeneTest(cds_granu[expressed_genes,],
              fullModelFormulaStr = "~clust_raw")
ordering_genes_clust <- row.names (subset(diff_test_res_clust, qval < 0.01))
cds_granu_clust <- setOrderingFilter(cds_granu, ordering_genes_clust)
cds_granu_clust <- reduceDimension(cds_granu_clust, max_components = 2,method = 'DDRTree')
cds_granu_clust <- orderCells(cds_granu_clust)


###Figures

##Figure1
#Figure1 A
ggplot(naive_UMI_cell,aes(x=CellCum,y=UMI_per_cell,group=1,col=group))+geom_line(size=1.5)+scale_x_log10()+scale_y_log10()+theme_light()+scale_color_discrete("",limits=c("Burn-in","Cells","Backgroud"))
#Figure1 B
ggplot(naive_gene_count)+geom_violin(aes(y=geneCount,x="Genes/Cell"),width=0.5)+geom_violin(aes(y=UMICount,x="UMI count/Cell"),width=0.5)+geom_boxplot(aes(y=geneCount,x="Genes/Cell"),width=0.05)+geom_boxplot(aes(y=UMICount,x="UMI count/Cell"),width=0.05)+scale_y_log10()+theme_light()
#Figure1 C
marker_gene<-c("Itgam","Ly6g","Cd48","Alas2")
visualize_gene_markers(gene_naive_log,marker_gene,gene_naive_res$tsne[c("TSNE.1","TSNE.2")])
#Figure1 D
ggplot(naive_cluster_tsne_type,aes(x=TSNE.1,y=TSNE.2,col=CellType ))+geom_point()

##Figure2
#Figure2 A
ggplot(naive_gene_count,aes(x=TSNE.1,y=TSNE.2,col=as.factor(clust_hclust)))+geom_point()+theme_classic()
#Figure2 B
pheatmap(gene_naive_tsne_dist,show_rownames = F,show_colnames = F,annotation_row = naive_cluster_tsne,annotation_col = naive_cluster_tsne)
#Figure2 D
select_genes<-gene_ids[gene_ids$gene_symbol %in% c("Ly6g","Cxcr4","Cxcr2","Vcam1","Cebpb","Cebpe","Cebpz","Spi1","Myb"),"ENSG_id"]
mt_select_gene_mt<-melt(select_naive_exp_tpm_mean[select_genes,])
colnames(mt_select_gene_mt) <-c("ENSG_id","clust","mean_exp")
mt_select_gene_mt<-merge(mt_select_gene_mt,gene_ids)
ggplot(mt_select_gene_mt,aes(x=clust,y=mean_exp,col=gene_symbol))+geom_line(aes(group=gene_symbol))+facet_grid(gene_symbol~.)+xlim(levels(mt_select_gene_mt$clust)[2:8])

##Figure3
#Figure3 A
random_select_50cell<-rbind(head(subset(naive_cluster_tsne,clust=="NG1"),50),
head(subset(naive_cluster_tsne2,clust=="NG2"),50),
head(subset(naive_cluster_tsne2,clust=="NG3"),50),
head(subset(naive_cluster_tsne2,clust=="NG4"),50),
head(subset(naive_cluster_tsne2,clust=="NG5"),50),
head(subset(naive_cluster_tsne2,clust=="NG6"),50),
head(subset(naive_cluster_tsne2,clust=="NG7"),50)) 
granu_random_select<-select_naive_exp_tpm_log[,rownames(random_select_50cell)]
granu_random_select<-granu_random_select[apply(granu_random_select,1,mad)>0,] 
random_select_50gene<-subset(naive_gene_clust,rownames(naive_gene_clust) %in% rownames(granu_random_select))
random_select_50gene$geneClust<-factor(random_select_50gene$geneClust,labels  = c("G1","G2","G3","G4"))

random_select_50gene2<-rbind(subset(random_select_50gene,geneClust=="G1"),
                    head(subset(random_select_50gene,geneClust=="G2"),29),
                    head(subset(random_select_50gene,geneClust=="G3"),8),
                    head(subset(random_select_50gene,geneClust=="G4"),11))  
to_heat<-t(scale(t(granu_random_select[rownames(random_select_50gene2)[order(random_select_50gene2$geneClust)],])))
to_heat[to_heat< -2] = -2
to_heat[to_heat> 2] = 2
pheatmap(to_heat,cluster_cols = F,cluster_rows = F,show_rownames = F,show_colnames = F,annotation_col = random_select_50cell,annotation_row =random_select_50gene2 ,color = colorpanel(100,low="blue",mid="white",high="red"))
#Figure3 B
mt_scale_select<-melt(scale_select)
colnames(mt_scale_select)<-c("Gene","CellClust","scale_exp")
mt_scale_select$GeneClust<-naive_gene_clust[mt_scale_select$Gene,"clust"]
ggplot(mt_scale_select)+geom_violin(data=mt_scale_select,aes(x=CellClust,y=scale_exp,fill=GeneClust))+facet_wrap(~GeneClust)+theme(legend.position = "none") +scale_fill_discrete("gene Clust",limits=levels(naive_gene_clust$geneClust)[c(3,2,4,1)]) 
 mt_scale_select_mean<-melt(tapply(mt_scale_select$scale_exp,mt_scale_select[,c("CellClust","GeneClust")],mean))
 ggplot()+geom_violin(data=mt_scale_select,aes(x=CellClust,y=scale_exp,fill=GeneClust))+geom_line(data=mt_scale_select_mean,aes(x=CellClust,y=value ,group=GeneClust))+facet_grid(GeneClust~.)+theme(legend.position = "none") +scale_fill_discrete("gene Clust",limits=levels(naive_gene_clust$geneClust)[c(3,2,4,1)])+scale_x_discrete("Cell Clust",labels=c("NG1","NG2","NG3","NG4","NG5","NG6","NG7"))+theme_classic()
#Figure3 C
gene_enrichement_granylocyte<-read.delim("granu_naive_tpm_clust.go_bp.top5.sig.chart",header=F) 
ggplot(gene_enrichement_granylocyte,aes(x=paste(V1,V5),y=-log10(V4),fill=as.factor(V5)))+geom_bar(stat="identity",width = 0.5)+coord_flip()+ylab("-log10(p-value)")+scale_fill_discrete("Gene Cluster")+scale_x_discrete("",limits=rev(as.character(paste(gene_enrichement_granylocyte$V1,gene_enrichement_granylocyte$V5))),labels=rev(as.character(gene_enrichement_granylocyte$V1)))

##Figure4 
#Figure4 A
ggplot(lps_gene_count)+geom_violin(aes(y=geneCount,x="Genes/Cell"),width=0.5)+geom_violin(aes(y=UMICount,x="UMI count/Cell"),width=0.5)+geom_boxplot(aes(y=geneCount,x="Genes/Cell"),width=0.05)+geom_boxplot(aes(y=UMICount,x="UMI count/Cell"),width=0.05)+scale_y_log10()+theme_light()
#Figure4 B
ggplot(naive_gene_count,aes(x=TSNE.1,y=TSNE.2,col=as.factor(clust_hclust)))+geom_point()+theme_classic()
#Figure4 C
select_genes<-gene_ids[gene_ids$gene_symbol %in% c("Ly6g","Cxcr4","Cxcr2","Vcam1","Cebpb","Cebpe","Cebpz","Spi1","Myb"),"ENSG_id"]
mt_select_gene_mt<-melt(select_lps_exp_tpm_mean[select_genes,])
colnames(mt_select_gene_mt) <-c("ENSG_id","clust","mean_exp")
mt_select_gene_mt<-merge(mt_select_gene_mt,gene_ids)
ggplot(mt_select_gene_mt,aes(x=clust,y=mean_exp,col=gene_symbol))+geom_line(aes(group=gene_symbol))+facet_grid(gene_symbol~.)+xlim(levels(mt_select_gene_mt$clust)[2:8])
#Figure4 D
mt_scale_select_lps<-melt(scale_select_lps)
colnames(mt_scale_select_lps)<-c("Gene","CellClust","scale_exp")
mt_scale_select_lps$GeneClust<-naive_gene_clust[mt_scale_select_lps$Gene,"clust"]
mt_scale_select_lps_mean<-melt(tapply(mt_scale_select_lps$scale_exp,mt_scale_select_lps[,c("CellClust","GeneClust")],mean))
ggplot()+geom_violin(data=mt_scale_select_lps,aes(x=CellClust,y=scale_exp,fill=GeneClust))+geom_line(data=mt_scale_select_lps_mean,aes(x=CellClust,y=value ,group=GeneClust))+facet_grid(GeneClust~.)+theme(legend.position = "none") +scale_x_discrete("Cell Clust",labels=c("LG1","LG2","LG3","LG4","LG5","LG6","LG7"))+theme_classic()
#Figure4 E
gene_enrichement_granylocyte_lps<-read.delim("granu_lps_tpm_clust.go_bp.top5.sig.chart",header=F) 
ggplot(gene_enrichement_granylocyte_lps,aes(x=paste(V1,V5),y=-log10(V4),fill=as.factor(V5)))+geom_bar(stat="identity",width = 0.5)+coord_flip()+ylab("-log10(p-value)")+scale_fill_discrete("Gene Cluster")+scale_x_discrete("",limits=rev(as.character(paste(gene_enrichement_granylocyte_lps$V1,gene_enrichement_granylocyte_lps$V5))),labels=rev(as.character(gene_enrichement_granylocyte_lps$V1)))

##Figure5
#Figure5 A
ggplot(cell_clust_tsne,aes(x=TSNE.1,y=TSNE.2,col=clust_raw))+geom_point()
#Figure5 B
pheatmap( acast(common_prop_table,clust1~clust2,value.var = "common_prop"),cluster_rows = F,cluster_cols = F,width = 120,height = 100)
#Figure5 C
plot_cell_trajectory(cds_granu_clust, color_by = "clust_color",cell_size = 0.5)
#Figure5 D
plot_ng_lg_pos<-function(clus1,clus2){
to_plot<-select_nl_exp_pos_perc[,c(clus1,clus2)]
colnames(to_plot)<-c("g1","g2")
to_plot<-subset(to_plot,g1>0&g2>0)
to_plot$color<-subset(diff_exp_pos_res_matrix2,clust==paste(clus1,clus2,sep="_"))[rownames(to_plot),"pattern"]
to_plot$color<-factor(to_plot$color,levels=c("LPS_up","LPS_down",NA))
ggplot(to_plot,aes(x=g1,y=g2,col=color))+geom_point()+geom_point(size=0.5)+geom_abline(slope=1,intercept = 0,col="blue")+xlab(clus1)+ylab(clus2)+theme_classic()+theme(legend.position = "none")
}
plot_ng_lg_pos("NG1","LG1")
plot_ng_lg_pos("NG2","LG2")
plot_ng_lg_pos("NG3","LG3")
plot_ng_lg_pos("NG4","LG4")
plot_ng_lg_pos("NG5","LG5")
plot_ng_lg_pos("NG6","LG6")
plot_ng_lg_pos("NG7","LG7")

##Figure6
#Figure6 B
diff_exp_pos_ng_enrich<-read.delim("diff_exp_pos_ng.clust.go_bp.top5.sig.chart",header=F)
colnames(diff_exp_pos_ng_enrich) <-c("GOTERM","q.value","type")
diff_exp_pos_ng_enrich<-merge(diff_exp_pos_ng_enrich,unique(diff_exp_pos_NG[,c("clust","pattern","type")]))
diff_exp_pos_ng_up<-rownames(diff_exp_pos_ng_enrich_matrix[order(apply(diff_exp_pos_ng_enrich_matrix,1,mean)),])
ggplot(subset(diff_exp_pos_ng_enrich,pattern=="LPS_up"),aes(x=GOTERM,y=-log10(q.value),fill=pattern))+geom_bar(stat="identity",width=0.5)+facet_grid(pattern~clust,scale="free_y",space="free_y")+coord_flip()+xlim(diff_exp_pos_ng_up)
#Figure6 C
diff_exp_ng_up<-diff_exp_ng_acast[apply(diff_exp_ng_acast,1,function(a){length(colnames(diff_exp_ng_acast_c)[a %in% c("LPS_up")])})>0,]
diff_exp_ng_fill<-melt(as.matrix(diff_exp_ng_up))
colnames(diff_exp_ng_fill) <- c("gene","clust","pattern")
diff_exp_ng_fill$pattern <-factor(diff_exp_ng_fill$pattern,levels = c("LPS_up","LPS_down","common","not_express"))
ggplot(diff_exp_ng_fill,aes(x=clust,stratum =pattern,alluvium =gene, fill=pattern, label=pattern))+geom_flow()+geom_stratum()+theme_classic()+xlab("clust")

