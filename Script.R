############################################################################################################
# Heat-map
############################################################################################################


library(ggplot2)
library(tidyverse)
library(scales)


# importation des 17756 gènes 

tpm= read.csv("genes_filtres.tsv",sep='\t')

colnames(tpm)[1] ="gene_id"


############################################################################################################
# Ordre des gènes du plus exprimé au moins exprimé
############################################################################################################


# Moyenne TPM par gène (colonnes 2 à 21)
moyennes <- rowMeans(tpm[, 2:21])

# Dictionnaire trié directement
d_tri <- as.list(sort(setNames(moyennes, tpm$gene_id)))

# Noms de gènes triés
genes_tri <- names(d_tri)

# Réordonner tpm selon l'ordre des gènes
tpm_tri <- tpm[match(genes_tri, tpm$gene_id), ]



tpm_long=tpm_tri %>% pivot_longer(cols = -gene_id, names_to = "sample", values_to = "TPM")
tpm_long= tpm_long %>% mutate(log2TPM = log2(TPM))

tpm_long <- tpm_long %>% mutate(gene_id = factor(gene_id, levels = genes_tri))

# Définir l'ordre des colonnes : paires O-T côte à côte

col_order = c("GSM5195615","GSM5195625", "GSM5195616","GSM5195626","GSM5195617","GSM5195627","GSM5195618","GSM5195628","GSM5195619","GSM5195629","GSM5195620","GSM5195630","GSM5195621","GSM5195631","GSM5195622","GSM5195632","GSM5195623","GSM5195633","GSM5195624","GSM5195634")


#Calcul du Pearson par paire tumeur-organoide
rownames(tpm_tri) <- NULL
tpm_wide <- tpm_tri %>% column_to_rownames("gene_id")
paires <- data.frame(
  tumeur = c("GSM5195615","GSM5195616","GSM5195617","GSM5195618","GSM5195619","GSM5195620","GSM5195621","GSM5195622","GSM5195623","GSM5195624"),
  organoid = c("GSM5195625","GSM5195626","GSM5195627","GSM5195628","GSM5195629","GSM5195630","GSM5195631","GSM5195632","GSM5195633","GSM5195634")
)

pearson_values <- mapply(function(tum,org) {
  round(cor(log2(tpm_wide[[tum]]),log2(tpm_wide[[org]]),method = "pearson"),2)
},paires$tumeur,paires$organoid)


#Labels avec Pearson affiché sous la tumeur de chaque paire

labels_simples <- c(
  "GSM5195615"= "T1", "GSM5195625" = "O1", #nouveau
  "GSM5195616" = "T2", "GSM5195626" = "O2", #nouveau
  "GSM5195617" = "T3", "GSM5195627" = "O3", #nouveau
  "GSM5195618" = "T4", "GSM5195628" = "O4", #nouveau
  "GSM5195619" = "T5", "GSM5195629" = "O5", #nouveau
  "GSM5195620" = "T6", "GSM5195630" = "O6", #nouveau
  "GSM5195621" = "T7", "GSM5195631" = "O7", #nouveau
  "GSM5195622" = "T8", "GSM5195632" = "O8", #nouveau
  "GSM5195623" = "T9", "GSM5195633" = "O9", #nouveau
  "GSM5195624" = "T10","GSM5195634" = "O10" #nouveau
)


gene_cible <- "1277"
pos_y <- which(genes_tri == gene_cible)    # indiquer le gène le plus différentiellement exprimé entre tumeurs primaires et organoïdes

#Appliquer l'ordre
tpm_long$sample = factor(tpm_long$sample,levels = col_order)

ggplot(tpm_long, aes(x = sample, y = gene_id, fill = log2TPM)) +
  geom_tile() +
  geom_vline(xintercept=seq(2.5, 18.5,by=2), color="white",linewidth=1.5)+
  scale_fill_gradientn(colours=c("darkblue","cadetblue","skyblue","white","yellow","orange","red"), name = expression(log[2]~TPM)) + 
  annotate("text",x=seq(1.5,19.5,by=2),y=-200,label=paste0("R=", pearson_values),size=3,fontface="italic")+
  annotate("text",x=21,y=-200,label="Pearson R",size=3,fontface="italic",hjust=0)+
  annotate("text",x=21,y=-1400,label="n = 17 756 gènes",size=3,hjust=0)+ 
  annotate("text",x=5.5,y=Inf,label="Traitement-naïve",size=4,fontface="bold",vjust=-0.5)+ 
  annotate("text",x=15.5,y=Inf,label="FOLFIRINOX-treated",size=4,fontface="bold",vjust=-0.5)+
  annotate("text",x=21,y=-2800,label="T : Tumeurs primaires\ncorrespondantes",size=3,hjust = 0,color="red")+ #nouveau
  annotate("text",x=21,y=-4200,label="O : Lignées d'organoïdes",size=3,hjust = 0,color="blue")+
  annotate("text", x = 21.2, y = pos_y,label = "COL1A1", hjust = 0, size = 3, color = "black") +
  annotate("segment",x = 21, xend = 20.6,y = pos_y, yend = pos_y,arrow = arrow(length = unit(0.2, "cm"), type = "closed"),color = "black", linewidth = 0.6) +
  coord_cartesian(clip="off")+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, colour = ifelse(grepl("^T", labels_simples),"red","blue")),
    axis.text.y =  element_blank(), axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold",margin=margin(b=40)),
    plot.margin = margin(t=10,r=40,b=60,l=10)
  ) +
  labs(title = "Heatmaps des gènes regroupés hiérarchiquement issues de deux groupes de patients") +
  scale_x_discrete(labels = labels_simples)









############################################################################################################
# DEG
############################################################################################################



# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with DESeq2
library(DESeq2)

# load counts table from GEO
urld<- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"

path_TPM <- paste(urld, "acc=GSE169321", "file=GSE169321_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl_TPM <- as.matrix(data.table::fread(path_TPM, header=T, colClasses="integer"), rownames="GeneID")

path_raw <- paste(urld, "acc=GSE169321", "file=GSE169321_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl_raw <- as.matrix(data.table::fread(path_raw, header=T, colClasses="integer"), rownames="GeneID")

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID


# pre-filter low count genes
# keep genes with at least N counts > 10, where N = size of smallest group
keep <- rowSums( tbl_TPM >=1 ) >= ncol(tbl_TPM)
tbl_raw <- tbl_raw[keep, ]




#Recherche des DEGs entre les organoïdes et les PDAC n'ayant pas eu de traitement
#################################################################################


# sample selection
gsms_NT <- "00000XXXXX11111XXXXX"
sml_NT <- strsplit(gsms_NT, split="")[[1]]

# filter out excluded samples (marked as "X")
sel_NT <- which(sml_NT != "X")
sml_NT <- sml_NT[sel_NT]
tbl_raw_NT <- tbl_raw[ ,sel_NT]

# group membership for samples
gs_NT <- factor(sml_NT)
groups_NT <- make.names(c("PDAC NT","Org NT"))
levels(gs_NT) <- groups_NT
sample_info_NT <- data.frame(Group = gs_NT, row.names = colnames(tbl_raw_NT))
ds_NT <- DESeqDataSetFromMatrix(countData=tbl_raw_NT, colData=sample_info_NT, design= ~Group)
ds_NT <- DESeq(ds_NT, test="Wald", sfType="poscount")

# extract results for top genes table
r_NT <- results(ds_NT, contrast=c("Group", groups_NT[1], groups_NT[2]), alpha=0.1, pAdjustMethod ="fdr")


# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2.5,1), cex.main=1.2)
plot(r_NT$log2FoldChange, -log10(r_NT$padj), main="Volcano plot des gènes exprimés par les tumeurs primaires et par les organoïdes\nn'ayant pas reçu le traitement",
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r_NT, padj<0.1 & abs(log2FoldChange) >= 1),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.1, sep=""), legend=c("down", "up"), pch=20,col=1:2)


df_DEGs_NT <- as.data.frame(r_NT)
df_DEGs_NT <- df_DEGs_NT[!is.na(df_DEGs_NT$padj) & df_DEGs_NT$padj < 0.1 & 
                           abs(df_DEGs_NT$log2FoldChange) >= 1, ]

DEGs_NT=rownames(df_DEGs_NT)



#Recherche des DEGs entre les organoïdes et les PDAC qui ont eu le traitement
#############################################################################

# sample selection
gsms_FOL <- "XXXXX00000XXXXX11111"
sml_FOL <- strsplit(gsms_FOL, split="")[[1]]

# filter out excluded samples (marked as "X")
sel_FOL <- which(sml_FOL != "X")
sml_FOL <- sml_FOL[sel_FOL]
tbl_raw_FOL <- tbl_raw[ ,sel_FOL]

# group membership for samples
gs_FOL <- factor(sml_FOL)
groups_FOL <- make.names(c("PDAC FOL","Org FOL"))
levels(gs_FOL) <- groups_FOL
sample_info_FOL <- data.frame(Group = gs_FOL, row.names = colnames(tbl_raw_FOL))
ds_FOL <- DESeqDataSetFromMatrix(countData=tbl_raw_FOL, colData=sample_info_FOL, design= ~Group)
ds_FOL <- DESeq(ds_FOL, test="Wald", sfType="poscount")

# extract results for top genes table
r_FOL <- results(ds_FOL, contrast=c("Group", groups_FOL[1], groups_FOL[2]), alpha=0.1, pAdjustMethod ="fdr")

# volcano plot
old.pal <- palette(c("#00BFFF", "#FF3030")) # low-hi colors
par(mar=c(4,4,2.5,1), cex.main=1.2)
plot(r_FOL$log2FoldChange, -log10(r_FOL$padj), main="Volcano plot des gènes exprimés par les tumeurs primaires et par les organoïdes\nayant reçu le traitement",
     xlab="log2FC", ylab="-log10(Padj)", pch=20, cex=0.5)
with(subset(r_FOL, padj<0.1 & abs(log2FoldChange) >= 1),
     points(log2FoldChange, -log10(padj), pch=20, col=(sign(log2FoldChange) + 3)/2, cex=1))
legend("bottomleft", title=paste("Padj<", 0.1, sep=""), legend=c("down", "up"), pch=20,col=1:2)


df_DEGs_FOL <- as.data.frame(r_FOL)
df_DEGs_FOL <- df_DEGs_FOL[!is.na(df_DEGs_FOL$padj) & df_DEGs_FOL$padj < 0.1 & 
                             abs(df_DEGs_FOL$log2FoldChange) >= 1, ]
DEGs_FOL=rownames(df_DEGs_FOL)




#################################################################################################################
# Diagrammes de Venn
#################################################################################################################

library(VennDiagram)


#DEGs entre Organoïdes et les PDAC pour les traités et non traités

y=list(A=DEGs_NT,B=DEGs_FOL)
venn.diagram(
  y,
  filename = "Diagramme de Venn des DEGs entre Organoïdes et les PDAC pour les traités et non traités.png",
  imagetype = "png",
  category.names = c("No Therapy", "FOLFIRINOX"),
  fill = c("blue", "red"),
  alpha = 0.5,
  height = 3000,
  width = 3000,
  cex = 2,
  resolution = 300,
  main = "DEGs entre Organoïdes et PDAC\n(traités vs non traités)",
  main.cex = 2,
  main.fontface = "bold",
  margin = 0.2,
  cat.cex = 1.5,
  main.pos = c(0.5, 0.90),
  cat.dist = 0.05
)
