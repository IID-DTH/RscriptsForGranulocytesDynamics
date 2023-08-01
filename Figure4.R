

########bulk RNA-seq
bg_all<-ballgown(dataDir = "./ballgown",samplePattern = "mRNA")
sample_bulk<-read.delim("sampleSheet_bulk.txt")
rownames(sample_bulk)<-sample_bulk$sampleID
sample_bulk$sampleID<-as.character(sample_bulk$sampleID)
sample_bulk$fpkm<-paste0("FPKM.",sample_bulk$sampleID)
pData(bg_all)<-sample_bulk[sampleNames(bg_all),]
gexpr_all<-gexpr(bg_all)
gexpr_all<-gexpr_all[,sample_bulk$fpkm]
colnames(gexpr_all)<-sample_bulk$sampleID

gexpr_all<-gexpr_all[apply(gexpr_all,1,max)>1,]
log_gexpr_all<-log(gexpr_all+1)

stat_results = stattest(bg, feature='gene', meas='FPKM', 
  covariate='group')
sample_bulk$Condition=factor(sample_bulk$Condition,levels = c("NAIVE","LPS","CLP","CLP-LPS"))


##Figure 4A
my_statest<-function(bg_subset,log=T,group){
  expr_subset=gexpr(bg_subset)
  if(log==T){
    expr_subset=log(expr_subset+1)
  }
exp_t_test_p<-function(a){
  t.test(a~pData(bg_subset)[,group])$p.value
}

p_value<-apply(expr_subset,1,exp_t_test_p)
exp_t_test_est<-function(a){
  t.test(a~pData(bg_subset)[,group])$estimate
}
est_sub<-apply(expr_subset,1,exp_t_test_est)

stat_res=data.frame(t(rbind(p_value,est_sub)))
stat_res$lgfc=stat_res[,3]-stat_res[,2]
stat_res$fdr=p.adjust(stat_res$p_value,method = "fdr")
return(stat_res)
}

stat_res_lclg1=my_statest(bg_lcl_g1,log = T,group="Condition")
stat_res_lclg1$cell="G1"
stat_res_lclg1$condition="LPS vs CLP-LPS"
stat_res_lclg1$lgfc=-stat_res_lclg1$lgfc
stat_res_lclg1$geneID=rownames(stat_res_lclg1)
stat_res_lclg2=my_statest(bg_lcl_g2,log = T,group="Condition")
stat_res_lclg2$cell="G2"
stat_res_lclg2$condition="LPS vs CLP-LPS"
stat_res_lclg2$lgfc=-stat_res_lclg2$lgfc
stat_res_lclg2$geneID=rownames(stat_res_lclg2)
stat_res_lclg3=my_statest(bg_lcl_g3,log = T,group="Condition")
stat_res_lclg3$cell="G3"
stat_res_lclg3$condition="LPS vs CLP-LPS"
stat_res_lclg3$lgfc=-stat_res_lclg3$lgfc
stat_res_lclg3$geneID=rownames(stat_res_lclg3)
stat_res_lcl=rbind(stat_res_lclg1,stat_res_lclg2,stat_res_lclg3)
stat_res_lcl$geneCond="none"

stat_res_lcl[stat_res_lcl$lgfc> log(2)&stat_res_lcl$p_value<0.1,"geneCond"]="Up"
stat_res_lcl[stat_res_lcl$lgfc< -log(2)&stat_res_lcl$p_value<0.1,"geneCond"]="Down"
ggplot(stat_res_lcl,aes(x=lgfc,y=-log10(p_value),col=geneCond))+geom_point()+facet_grid(condition~cell)+geom_hline(yintercept = -log10(0.05))+geom_vline(xintercept = c(log(2),log(0.5)))+scale_color_manual(values=brewer.pal(9,"Set1")[c(2,9,1)])+theme_light()

##Figure 4B
log_gexpr_lcl=t(scale(t(log_gexpr_lcl)))
log_gexpr_G1_lcl=log_gexpr_lcl[,c(1,2)]
log_gexpr_G1_lcl=log_gexpr_G1_lcl[order(log_gexpr_G1_lcl[,1]),]
 pheatmap(log_gexpr_G1_lcl,cluster_rows =F,cluster_cols  = F,show_rownames = F,color= colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
 
 log_gexpr_G2_lcl=log_gexpr_lcl[,c(3,4)]
log_gexpr_G2_lcl=log_gexpr_G2_lcl[order(log_gexpr_G2_lcl[,1]),]
 pheatmap(log_gexpr_G2_lcl,cluster_rows =F,cluster_cols  = F,show_rownames = F,color= colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
 
 log_gexpr_G3_lcl=log_gexpr_lcl[,c(5,6)]
log_gexpr_G3_lcl=log_gexpr_G3_lcl[order(log_gexpr_G3_lcl[,1]),]
 pheatmap(log_gexpr_G3_lcl,cluster_rows =F,cluster_cols  = F,show_rownames = F,color= colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

##Figure 4C
go_bp1=read.delim("GO_BP1_david_enrich.txt")
p1= ggplot(subset(go_bp1,CellName=="G1"),aes(x=Term,y=-log10(PValue)))+geom_bar(stat="identity")+coord_flip()+theme_classic()+xlim(rev(as.character(subset(go_bp1,CellName=="G1")$Term)))
p2= ggplot(subset(go_bp1,CellName=="G2"),aes(x=Term,y=-log10(PValue)))+geom_bar(stat="identity")+coord_flip()+theme_classic()+xlim(rev(as.character(subset(go_bp1,CellName=="G2")$Term)))
p3= ggplot(subset(go_bp1,CellName=="G3"),aes(x=Term,y=-log10(PValue)))+geom_bar(stat="identity")+coord_flip()+theme_classic()+xlim(rev(as.character(subset(go_bp1,CellName=="G3")$Term)))

##Figure 4D
entrez_G3=data.frame(entrizID=idConverter(ids=subset(stat_res_lcl,cell=="G3"&geneCond!="none")$geneID,srcSpecies="MUSMU",destSpecies="MUSMU",srcIDType="ENSEMBL",destIDType="EG"))
kk_G3=enrichKEGG(gene=entrez_G3$entrizID,organism  = 'mmu',pvalueCutoff  = 0.05)
dotplot(kk_G3)

##Figure 4E-G
log_gexpr_all_scale=t(scale(t(log_gexpr_all)))

list_tlr4=read.delim("tlr4_pathway.txt")
list_glu=read.delim("glycy_pathway.txt")
list_ppp=read.delim("ppp_pathway.txt")
list_fatty=read.delim("fatty_acid_pathway.txt")

list_tlr4=merge(list_tlr4,geneIDmap,by.x="Gene",by.y="geneName")
list_tlr4=list_tlr4[list_tlr4$ID %in% rownames(log_gexpr_all),]

list_glu=merge(list_glu,geneIDmap,by.x="Gene",by.y="geneName")
list_glu=list_glu[list_glu$ID %in% rownames(log_gexpr_all),]

list_ppp=merge(list_ppp,geneIDmap,by.x="Gene",by.y="geneName")
list_ppp=list_ppp[list_ppp$ID %in% rownames(log_gexpr_all),]

list_fatty=merge(list_fatty,geneIDmap,by.x="Gene",by.y="geneName")
list_fatty=list_fatty[list_fatty$ID %in% rownames(log_gexpr_all),]


exp_weighted_score=data.frame(
  tlr4_score=apply(log_gexpr_all[list_tlr4$ID,]*list_tlr4$weight,2,sum),
  glu_score=apply(log_gexpr_all[list_glu$ID,]*list_glu$weight,2,sum),
  ppp_score=apply(log_gexpr_all[list_ppp$ID,]*list_ppp$weight,2,sum),
  fatty_score=apply(log_gexpr_all[list_fatty$ID,]*list_fatty$weight,2,sum))

exp_weighted_score$sampleID=rownames(exp_weighted_score)
exp_weighted_score=merge(exp_weighted_score,sample_bulk)

p1=ggplot(subset(exp_weighted_score,CellName!="M"),aes(x=CellName,fill=Condition,y=tlr4_score))+geom_boxplot()+theme_classic()+ggtitle("TLR4 pathway: 53 genes")
p2=ggplot(subset(exp_weighted_score,CellName!="M"),aes(x=CellName,fill=Condition,y=glu_score))+geom_boxplot()+theme_classic()+ggtitle("Glycolysis and pyruvate oxidatite pathway: 17 genes")
p3=ggplot(subset(exp_weighted_score,CellName!="M"),aes(x=CellName,fill=Condition,y=ppp_score))+geom_boxplot()+theme_classic()+ggtitle("PPP pathway: 7 genes")
p4=ggplot(subset(exp_weighted_score,CellName!="M"),aes(x=CellName,fill=Condition,y=fatty_score))+geom_boxplot()+theme_classic()+ggtitle("Fatty acid biosynthesis pathway: 5 genes")

grid.arrange(p1,p2,p3,p4)



