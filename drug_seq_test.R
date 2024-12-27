
#Drug-Seq  分析
#参考文献：DRUG-seq Provides Unbiased Biological Activity Readouts for Neuroscience Drug Discovery
#数据集：GSE176150，GSM5357044，VH02001612_S9，SRR14730301
#github上有一些脚本： https://github.com/Novartis/DRUG-seq
#重复出一个样本的上游定量分析流程
setwd("D:\\Cell_panel_screening\\DRUG-Seq\\GSM5357044")
#360 wells, 5w+ genes
m=load("GSM5357044_flowcell_4000_UMI_decode_VH02001612.RData")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

#各样本文库大小
colSums(UMI_decode) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=.))+
  geom_histogram()+
  theme_bw(base_size = 20)

meta=fread("cb1c00920_si_001/File_S1_metadata.csv") %>% 
  filter(plate_barcode=="VH02001612")

#3个sectors，8个 dose，所以每个药物共有24例数据
meta %>% count(compound)

library(limma)
library(edgeR)
library(umap)
library(plotly)

tmp=UMI_decode[rowSums(UMI_decode)>0,]
d<-DGEList(counts=tmp)
d <- calcNormFactors(d, method="TMM")
mat_ump_wbatchE<-log2(edgeR::cpm(d,  normalized.lib.sizes=T,log=F)+1)#scale到一个相对可比的范围    

# set.seed(87654)
# pca_umap <- function(counts=NULL,center_pca=T,metric="euclidean", scale_pca=T,n_neighbors=15,n_pca=20){
#   pca<-prcomp(t(counts), scale = scale_pca, center = center_pca)
#   rownames(pca$x)  <-  colnames(counts)
#   umap  <- umap(pca$x[,1:n_pca],method = "naive",metric=metric, n_neighbors=n_neighbors)
#   rownames(umap$layout)  <-  colnames(counts)
#   return(umap)
# }
counts=mat_ump_wbatchE
pca<-prcomp(t(counts), scale = T, center = T)
rownames(pca$x)  <-  colnames(counts)
umap  <- umap(pca$x[,1:20],method = "naive",n_neighbors=15)

#umap结果可视化
colnames(umap$layout)<-c('umap1','umap2')
color=meta %>% select(Sample,plate_barcode,treatment,batch,treatment_dose,compound) %>% as.data.frame()
rownames(color)=meta %>% mutate(ID=paste0(plate_barcode,"_",well_index)) %>% pull(ID)

tmp=cbind(
  umap$layout,
  color[rownames(umap$layout),]
)

p=tmp %>% plot_ly(x = ~umap1, y = ~umap2,
                color = ~compound,type = "scatter",mode = 'markers', 
                size = ~treatment_dose,fill = ~'',
                text = ~compound)
htmlwidgets::saveWidget(as_widget(p), "tmp_500plotly.html")

#tmp_labels=tmp %>% arrange(desc(treatment_dose)) %>% distinct(.,compound,.keep_all = T)
tmp_labels=tmp %>% arrange(desc(umap2)) %>%  distinct(.,Sample,.keep_all = T)
library(RColorBrewer)
compounds=unique(tmp$compound)
compounds_colors=colorRampPalette(brewer.pal(8, "Set2"))(length(compounds))
names(compounds_colors)=compounds

p=tmp %>% 
  ggplot(aes(x=umap1,y=umap2))+
  geom_point(
    aes(colour=compound,size=treatment_dose)
  )+
  scale_colour_manual(
    values = c(compounds_colors)
  )+
  ggrepel::geom_text_repel(
    aes(x=umap1,y=umap2,label=compound),
    max.overlaps = 100,
    data=tmp_labels
  )+
  facet_wrap(~treatment_dose,nrow=2)+
  theme_bw(base_size = 25)

ggsave('umap1.png',p,width = 20,height = 10)



#heatmap非监督聚类，看一看不同药物之间的clusters
summary(
  apply(mat_ump_wbatchE,
        1,
        sd)
)
x1=mat_ump_wbatchE[ apply(mat_ump_wbatchE,
                          1,
                          sd) >0.9,]

library(pheatmap)
#选择样本进行可视化
select_samples=meta %>% filter(treatment_dose>3) %>% mutate(ID=paste0(plate_barcode,"_",well_index)) %>% 
  select(c(ID,compound,treatment_dose))

annotation_col=data.frame(
  Compounds=select_samples$compound,
  Dose=select_samples$treatment_dose
)
rownames(annotation_col)=select_samples$ID

pheatmap(
  x1[,select_samples$ID],
  cluster_cols = T,
  cluster_rows = F,
  annotation_col = annotation_col,
  annotation_colors = list(
    Compounds=compounds_colors
  ),
  show_rownames = F,
  show_colnames = F,
  scale = 'none',
  filename = 'heatmap.png'
)


#差异表达分析
#findmarkers或许可以借用
#或者使用edgeR，14compounds 分别与DMOS进行比较，得到一些significant genes
library(Seurat)
select_samples=meta %>% filter(treatment_dose==10 | compound=='DMSO') %>% mutate(ID=paste0(plate_barcode,"_",well_index)) %>% 
  select(c(ID,compound,treatment_dose))
select_counts=UMI_decode[,select_samples$ID]

obj=CreateSeuratObject(
  counts = select_counts,
  project= "Drug_seq"
)
obj@meta.data=cbind(obj@meta.data,select_samples)

obj$Group=as.factor(obj$compound)
Idents(obj)='compound'
de.res <- list()
for (i in levels(obj$Group)) {
  print(i)
  if(i!="DMSO"){
    DEG <- FindMarkers(obj, 
                       #第一组哪些样本
                       ident.1 = i, 
                       #第2组哪些样本
                       ident.2 = 'DMSO',
                       test.use = "DESeq2",
                       verbose = FALSE)
    DEG$compound <- i 
    DEG$Gene <- rownames(DEG)
    rownames(DEG) <- NULL
    DEG <- DEG[-1,]
    de.res[[i]] <- DEG
  }
}
out <- do.call(rbind, lapply(de.res, data.frame))
#out$Compare <- paste(levels(combined$Group)[1], "_vs_", levels(combined$Group)[2], sep = "")
#out <- out[,c("CellType", "Gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "Compare")]
openxlsx::write.xlsx(out, "compoundsType_degs_deseq2.xlsx", rowNames = FALSE)
#部分结果可视化
heatmap_markers=out %>%  filter(p_val_adj<0.05 & abs(avg_logFC)>2) %>% 
  group_by(compound) %>% 
  slice_max(.,avg_logFC,n=10) %>% 
  pull(Gene)

obj1=obj
DoHeatmap(object = obj1, features = heatmap_markers,slot='counts')


#对不同药物作用后差异基因进行功能分析
out=readxl::read_xlsx("compoundsType_degs_deseq2.xlsx") %>% 
  dplyr::mutate(GeneSymbol=gsub('\\,grch.*','',Gene)) %>% select(-c(Gene))
#对不同compounds差异基因进行功能分析
getEntrezID=function(geneSymbol){
  library(easyConvert)
  tmp=easyConvert(species = "HUMAN",queryList = geneSymbol,queryType = "SYMBOL")
  # tmp=clusterProfiler::bitr(
  #   geneSymbol,
  #   'SYMBOL',
  #   'ENTREZID',
  #   OrgDb = "org.Hs.eg.db"
  # )
  return(as.character(tmp$ENTREZID))
}

compounds=unique(out$compound)
ups=lapply(
  compounds,
  FUN=function(x){
    gs=out %>% filter(compound==x) %>% filter(p_val_adj<0.05 ) %>% pull(GeneSymbol)
    gs_list=getEntrezID(gs)
    return(gs_list)
  }
)
names(ups)=compounds
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
ck_up <- compareCluster(geneCluster = ups, 
                        fun = "enrichGO",
                        OrgDb="org.Hs.eg.db",pvalueCutoff=0.05,qvalueCutoff  = 0.05)
#将EntrezID，转换为Gene sybmol
ck_up_geneSymbol <- setReadable(ck_up, 'org.Hs.eg.db', 'ENTREZID')
dotplot(ck_up,showCategory = 5)

#进行尝试进行GSEA分析吧
library(msigdbr)
library(GSEABase)
library(easyConvert)
#C2, 包括Chemical and genetic perturbations (CGP) and Canonical pathways
b = msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, entrez_gene)
saveRDS(b,'gsea_c2.rds')

#ontology gene set
C5=msigdbr(species = "Homo sapiens", category = "C5") %>%
  dplyr::select(gs_name, entrez_gene)
saveRDS(C5,'gsea_c5.rds')
b=C5

# oncogenic signature gene sets
C6=msigdbr(species = "Homo sapiens", category = "C6") %>%
  dplyr::select(gs_name, entrez_gene)
saveRDS(C6,'gsea_c6.rds')
b=C6

#C4: computational gene sets：cancer cell, cancer gene, cancer modules
C4=msigdbr(species = "Homo sapiens", category = "C4") %>%
  dplyr::select(gs_name, entrez_gene)
saveRDS(C4,'gsea_c4.rds')
b=C4


gsea_results=lapply(
  compounds,
  FUN=function(x){
    print(x)
    gs=out %>% filter(compound==x) %>% filter(p_val_adj<0.05 ) %>% dplyr::select(c('GeneSymbol','avg_logFC'))
    
    if(nrow(gs)>1){
      a1 = gs$avg_logFC
      a2 = easyConvert(
        species = "HUMAN",
        queryList = gs$GeneSymbol,
        queryType = "SYMBOL"
      )
      # a2=clusterProfiler::bitr(
      #   gs$GeneSymbol,
      #   'SYMBOL',
      #   'ENTREZID',
      #   OrgDb = "org.Hs.eg.db"
      # )
      
      names(a1) = a2$ENTREZID
      b1 = a1[order(a1, decreasing = T)]
      
      gsearesults = clusterProfiler::GSEA(b1, TERM2GENE = b,pvalueCutoff=1)
      return(gsearesults)
    }else{
      gsearesults=list()
      return(gsearesults)
    }
  }
)
names(gsea_results)=compounds
#saveRDS(gsea_results,'gsea_results.rds')
saveRDS(gsea_results,'gsea_c6_results.rds')
#GSEA结果可视化
ps_gseas=list()
ps_gseas=lapply(
  compounds,
  FUN = function(x){
    print(x)
    tmp=gsea_results[[x]]
    if(length(tmp)){
      
      if(nrow(tmp@result)){
        p=tmp@result %>% mutate(Description = stringr::str_to_title(gsub("_", " ", Description))) %>%
          mutate(Description = reorder(Description, setSize)) %>% head(10) %>%
          ggplot(aes(
            x = enrichmentScore,
            y = Description,
            size = setSize,
            color = pvalue
          )) +
          geom_point() +
          scale_size(range = c(1, 5), name = "Size") +
          scale_color_gradient(low = "blue", high = "red") +
          xlab("EScore") + ylab("") + labs(title=x)+
          theme(axis.text.y = element_text(size = 10, color = "black"))+
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 25))+
          theme_bw(base_size = 15)
      }else{
        p=ggplot()
      }
      
    }else{
      p=ggplot()
    }
    return(p)
  }
)


ggsave(
  filename = "GSEA_c6_plots.pdf", 
  plot = gridExtra::marrangeGrob(ps_gseas, nrow=1, ncol=2), 
  width = 12, height = 6
)
#不同的GSEA分析结果中，确实能找到一些clues！！！

#对个别GSEA结果进行可视化
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
gsea_results=readRDS('gsea_c6_results.rds')
enrichplot::gseaplot2(gsea_results$`(S)-Crizotinib`, 
                      geneSetID = 1:3, pvalue_table = TRUE)

enrichplot::gseaplot2(gsea_results$AZD8055, 
                      geneSetID = 1:3, pvalue_table = TRUE)
enrichplot::gseaplot2(gsea_results$Cmp_334, 
                      geneSetID = 1:3, pvalue_table = TRUE)


enrichplot::gseaplot2(gsea_results$Brusatol, 
                      geneSetID = 1:3, pvalue_table = TRUE)
enrichplot::gseaplot2(gsea_results$Homoharringtonine, 
                      geneSetID = 1:3, pvalue_table = TRUE)
enrichplot::gseaplot2(gsea_results$Triptolide, 
                      geneSetID = 1:3, pvalue_table = TRUE)


#尝试一下TSNE,UMAP聚类
select_samples=meta %>% mutate(ID=paste0(plate_barcode,"_",well_index)) %>% 
  select(c(ID,compound,treatment_dose))
select_counts=UMI_decode[,select_samples$ID]
obj=CreateSeuratObject(
  counts = select_counts,
  project= "Drug_seq"
)
obj@meta.data=cbind(obj@meta.data,select_samples)

obj1 <- NormalizeData(object = obj)
obj1 <- FindVariableFeatures(object = obj1)
obj1 <- ScaleData(object = obj1)
obj1 <- RunPCA(object = obj1)
obj1 <- FindNeighbors(object = obj1, dims = 1:30)
obj1 <- FindClusters(object = obj1)
obj1 <- RunUMAP(object = obj1, dims = 1:30)
obj1 <- RunTSNE(object = obj1, dims = 1:30)

Idents(obj1)="compound"
DimPlot(object = obj1, reduction = "umap",label = T)
DimPlot(object = obj1, reduction = "tsne",label = T)

# p=DimPlot(obj1, reduction = "umap", split.by="treatment_dose", label = T,label.size = 5,
#           repel = T) +NoLegend()
p1=DimPlot(obj1, reduction = "umap", label = T,label.size = 5,
          repel = T) +NoLegend()
p2=DimPlot(obj1, reduction = "tsne", label = T,label.size = 5,
           repel = T) +NoLegend()
ps=ggpubr::ggarrange(p1,p2,nrow = 1)
ggsave("umap_tsne_seurat.png",ps, width =12, height = 6)


######################################
#比较我们分析结果与文章中分析结果的一致性
#
#####################################
# writeLines('PATH="D:\\r\\rtools35\\Rtools\\bin;${PATH}"', con = ".Renviron")
# Sys.which("make")
# install.packages("Rcpp")
fread("test/UMI_counts_in_Barcodes.csv") %>% 
  ggplot(aes(x=N))+
  geom_histogram()

fread("test/UMI_duplicates.csv") %>% 
  ggplot(aes(x=Rate))+
  geom_histogram()


test_count=UMI_decode %>% as.data.frame()
colnames(test_count) = gsub("VH02001612_","",colnames(UMI_decode))
#rownames(test_count) =gsub("\\,.*","",rownames(UMI_decode))
tmp1=test_count %>% mutate(ID=gsub("\\,.*","",rownames(UMI_decode))) %>% 
  distinct(.,ID,.keep_all = T)

#对360个样本的表达进行相关性分析
tmp_count=fread("test/sample1_count_matrix_gene_level.csv") %>% select(-c(V1)) %>% as.data.frame()
rownames(tmp_count)=fread("test/sample1_count_matrix_gene_level.csv") %>% pull(V1)
tmp_count[is.na(tmp_count)]=0
p1=colSums(tmp_count) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=.))+
  geom_histogram()+
  theme_bw(base_size = 20)+
  xlim(c(100000,1600000))+
  labs(x='',y='N',title = 'WuXi')

p2=colSums(UMI_decode) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=.))+
  geom_histogram()+
  theme_bw(base_size = 20)+
  xlim(c(100000,1600000))+
  labs(x='',y='N',title = 'GSM5357044')

ps=ggarrange(p1,p2,nrow = 1)
ggsave('tmp.png',ps,width = 10,height = 6)

#使用文库矫正之后的表达数据再计算一次,
#结论： CPM的相关性要明显高于基于counts值计算的相关性
tmp=tmp_count[rowSums(tmp_count)>0,]
d<-DGEList(counts=tmp)
d <- calcNormFactors(d, method="TMM")
wuxi_cpm<-log2(edgeR::cpm(d,  normalized.lib.sizes=T,log=F)+1)#scale到一个相对可比的范围  

test_count=UMI_decode %>% as.data.frame()
tmp=test_count[rowSums(test_count)>0,]
d<-DGEList(counts=tmp)
d <- calcNormFactors(d, method="TMM")
GSM5357044_cpm<-log2(edgeR::cpm(d,  normalized.lib.sizes=T,log=F)+1)

colnames(GSM5357044_cpm) = gsub("VH02001612_","",colnames(GSM5357044_cpm))
#rownames(test_count) =gsub("\\,.*","",rownames(UMI_decode))
tmp1=GSM5357044_cpm %>% as.data.frame() %>% mutate(ID=gsub("\\,.*","",rownames(GSM5357044_cpm))) %>% 
  distinct(.,ID,.keep_all = T)

df=wuxi_cpm %>% as.data.frame() %>% mutate(ID=rownames(.)) %>% 
  inner_join(.,
             tmp1,
             by=c("ID"="ID"),
             suffix = c("_WuXi","_GSM5357044")) 
  
barcodes=colnames(wuxi_cpm)

corx=sapply(
  barcodes,
  FUN=function(x){
    tmp=df %>% select(contains(x))
    res=cor.test(tmp[,1],tmp[,2])
    return(round(res$estimate,3))
  }
)

#总体上看来相关性挺高，应该问题不大，后续下载完整的fastqz文件，再进行整体测试
#1. drug_seq_v1.py
#2. 提取基因水平的表达变化
#3. 输出结果并比较
p1=corx%>% as.data.frame() %>% setnames('Cor') %>% 
  ggplot(aes(x=Cor))+
  geom_histogram()+
  theme_bw(base_size = 20)+
  labs(x='PCC',y='N',title = 'WuXi vs GSM5357044')
ggsave('tmp.png',p1,width = 6,height = 6)
corx%>% as.data.frame() %>% setnames('Cor') %>% 
  mutate(well_index=gsub("\\.cor","",rownames(.))) %>% 
  left_join(.,
            meta) %>% 
  ggplot(aes(x=Cor)) +
  geom_histogram()+
  facet_wrap(~treatment_dose)

#分析结果UMAP的重复性
tmp=tmp_count[rowSums(tmp_count)>0,]
d<-DGEList(counts=tmp)
d <- calcNormFactors(d, method="TMM")
mat_ump_wbatchE<-log2(edgeR::cpm(d,  normalized.lib.sizes=T,log=F)+1)#scale到一个相对可比的范围    

counts=mat_ump_wbatchE
pca<-prcomp(t(counts), scale = T, center = T)
rownames(pca$x)  <-  colnames(counts)
umap  <- umap(pca$x[,1:20],method = "naive",n_neighbors=15)

#umap结果可视化
colnames(umap$layout)<-c('umap1','umap2')
color=meta %>% select(Sample,plate_barcode,treatment,batch,treatment_dose,compound) %>% as.data.frame()
rownames(color)=meta %>% pull(well_index)

tmp=cbind(
  umap$layout,
  color[rownames(umap$layout),]
)

p=tmp %>% plot_ly(x = ~umap1, y = ~umap2,
                  color = ~compound,type = "scatter",mode = 'markers', 
                  size = ~treatment_dose,fill = ~'',
                  text = ~compound)
htmlwidgets::saveWidget(as_widget(p), "tmp_500plotly_wuxi.html")

#tmp_labels=tmp %>% arrange(desc(treatment_dose)) %>% distinct(.,compound,.keep_all = T)
tmp_labels=tmp %>% arrange(desc(umap2)) %>%  distinct(.,Sample,.keep_all = T)
library(RColorBrewer)
compounds=unique(tmp$compound)
compounds_colors=colorRampPalette(brewer.pal(8, "Set2"))(length(compounds))
names(compounds_colors)=compounds

p=tmp %>% 
  ggplot(aes(x=umap1,y=umap2))+
  geom_point(
    aes(colour=compound,size=treatment_dose)
  )+
  scale_colour_manual(
    values = c(compounds_colors)
  )+
  ggrepel::geom_text_repel(
    aes(x=umap1,y=umap2,label=compound),
    max.overlaps = 100,
    data=tmp_labels
  )+
  facet_wrap(~treatment_dose,nrow=2)+
  theme_bw(base_size = 25)

ggsave('umap1_wuxi.png',p,width = 20,height = 10)



#对样本进行分组相关性分析，by compound and dose
compounds
treatment_dose=meta %>% filter(compound!="DMSO") %>% count(treatment_dose) %>% pull(treatment_dose)
res=sapply(
  treatment_dose,
  FUN = function(x){
    abc_dose=x
    res=sapply(
      compounds,
      FUN=function(x,dose){
        well_index=meta %>% filter(compound==x) %>% filter(treatment_dose==dose) %>% pull(well_index)
        if(length(well_index)>1){
          tmp=GSM5357044_cpm %>% as.data.frame() %>% dplyr::select(well_index)
          tmp_cor=cor(tmp)
          return(median(tmp_cor[upper.tri(tmp_cor)]))
        }
      },
      dose=abc_dose
    )
    return(res)
  }
)
res_df=res %>% as.data.frame() %>% setnames(as.character(treatment_dose))

#不同化合物，不同浓度下PCC分布
p=res_df %>% mutate(compound=rownames(.)) %>% 
  tidyr::pivot_longer(!compound,names_to = 'dose',values_to = 'Cor') %>%
  filter(compound!="DMSO") %>% as.data.frame() %>% 
  ggplot(aes(y=compound,x=as.numeric(Cor)))+
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(aes(colour=dose))+
  theme_bw(base_size = 15)+
  labs(x='PCC',y='',title='PCC across replicates by compounds')
ggsave('tmp.png',p,width = 6,height = 6)

#不同浓度下PCC分布
p=res_df %>% mutate(compound=rownames(.)) %>% 
  tidyr::pivot_longer(!compound,names_to = 'dose',values_to = 'Cor') %>%
  filter(compound!="DMSO") %>% as.data.frame() %>% 
  mutate(Dose=forcats::fct_reorder(dose,-desc(as.numeric(dose)))) %>% 
  ggplot(aes(y=Dose,x=as.numeric(Cor)))+
  geom_boxplot()+
  theme_bw(base_size = 15)+
  labs(x='PCC',y='Dose',title='PCC across dose')
ggsave('tmp.png',p,width = 4,height = 6)
 