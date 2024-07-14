library(Seurat)
library(ggplot2)
library(patchwork)
library(DESeq2) 
library(viridis)
library(tidyverse)
library(edgeR)
library(limma)



seurat_young <- readRDS("/home/alsultfh/H1_Project/DATA/Dani/seurat_young_updated.rds")

plotter <- FetchData(seurat_young, vars=c('umap_1', 'umap_2', 'CellType'))



cluster_colors <- c(
	'#A4CEE2',
	'#1F77B5',
	'#2B6E38', 
	'#BBD58D', 
	'#842758', 
	'#E1BB6D', 
	'#C3B1CC', 
	'#673F95', 
	'#010000', 
	'#E99C9A', 
	'#E0A630', 
	'#D33C28')

names(cluster_colors) <- levels(plotter$CellType)

pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig1.pdf')
ggplot(plotter, aes(x=umap_1, y=umap_2, color=CellType))+
geom_point(alpha=0.8, size=0.2) + theme_classic() + scale_color_manual(values=cluster_colors) +
labs(x='UMAP 1', y='UMAP 2')+
theme(legend.position='bottom', legend.direction = "horizontal", legend.background=element_blank()) + 
guides(color = guide_legend(nrow=3, byrow=TRUE, override.aes = list(size=4, alpha=0.9)))
dev.off()



pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig2.pdf')
DotPlot(object = seurat_young, 
          features = unique(c('CRHBP', 'HOPX','CD34','PTPRC','SPINK2', 'MPO','CSF3R','CTSG','PRTN3','MPO', 'ELANE','AZU1','CEBPA','CST7','LYZ','CSTA','IL3RA', 'LSP1', 'IRF8','IL7R', 'SATB1','DNTT','VPREB1','EBF1','CD79A','CD79B','TAL1','GATA1','HBD','HBB','CA1','AHSP','KLF1','RUNX1','HDC','MS4A2','MS4A3','TPSAB1')), 
          group.by = 'CellType',
          col.min = -2.5,
         col.max = 2.5,
         dot.min = 0,
         dot.scale = 6) +
viridis::scale_color_viridis(option='viridis') +
labs(y='Cell types', x='Gene markers')+
theme(axis.text.x = element_text(size = 12, angle=90, hjust=1, vjust=0.5), 
      axis.text.y=element_text(size=12))  + 
      coord_flip()
dev.off()


pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig3.pdf', width=10)
VlnPlot(seurat_young, features='RPLP1', cols=cluster_colors, pt.size =0) + NoLegend() + labs(x='Cell types', y='Normalized expression')
dev.off()


data <- as.data.frame(seurat_young@assays$RNA$data)
h1.matrix <- data[c("H1F0", "HIST1H1A", "HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "H1FX"), ]
hk.matrix <- data[c('RPLP1'), , drop=FALSE]

#devide the Hl matrix by each HK matrix 
h1.matrix_RPLP1 <- h1.matrix / t(hk.matrix)[,1]


corrected_obj <- CreateSeuratObject(h1.matrix_RPLP1)
corrected_obj$cellType <- seurat_young$CellType



pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig4.pdf')
Idents(corrected_obj) <- 'cellType'
DotPlot(corrected_obj, 
        features = c("H1F0", "HIST1H1A", "HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "H1FX"), 
        group.by = "cellType") +
  theme(axis.text.x = element_text(size = 12, angle=90, hjust=1, vjust=0.5), 
      axis.text.y=element_text(size=12)) + 
      labs(y='Cell types', x= 'Normalized expression') +
      viridis::scale_color_viridis(option='viridis') +
      coord_flip()

dev.off()



dge <- DGEList(seurat_young@assays$RNA$counts)
dge <- calcNormFactors(dge)
gene_2_plot <- c("H1F0", "HIST1H1A", "HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "H1FX")
filtered_data <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
filtered_data <- filtered_data[gene_2_plot, ]
filtered_data <- t(filtered_data)
annotation <- FetchData(seurat_young, vars=c('CellType'))

filtered_data <- merge(filtered_data, annotation, by=0)
tmp <- reshape2::melt(filtered_data)

pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig5.pdf', width=10)
ggplot(tmp, aes(x=CellType, y=value, fill=CellType, color=CellType))+
geom_boxplot(outlier.size=1, outlier.alpha=0.1) +
theme_classic() + scale_fill_manual(values=cluster_colors) +
scale_color_manual(values=cluster_colors) +
labs(x='Cell types', y = 'Normalized expression') +
theme(legend.position = 'none', 
axis.text.x = element_text(size = 12, angle=90, hjust=1, vjust=0.5)) +
facet_wrap(~variable, nrow=3)
dev.off()





cts <- AggregateExpression(seurat_young, 
                           group.by = c("CellType", "orig.ident"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)


cts <- cts$RNA
cts_norm <- cpm(cts, log=TRUE)

cts_norm <- cts_norm[gene_2_plot, ]

plotter <- reshape2::melt(cts_norm)
plotter$Var2 <- factor(stringr::str_remove(plotter$Var2, "_mo193|_mo194|_mo196$"), 
                levels=c("HSC", "LMPP", "Cycling_LMPP","CLP", "MEP", "GMP1", "GMP2", "PreB", "Monocytes", "Basophils", "DC","Erythroid"))


pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig6.pdf', width=12)
ggplot(plotter, aes(x= Var2, y= value, fill= Var2))+
geom_boxplot(outlier.size=1, outlier.alpha=0.1) +
theme_bw() + scale_fill_manual(values=cluster_colors) +
scale_color_manual(values=cluster_colors) +
labs(x='Cell Types', y = 'LogCPM expression') +
theme(legend.position = 'none', strip.background = element_rect(fill="white"), 
axis.text.x = element_text(size = 12, angle=90, hjust=1, vjust=0.5)) +
facet_wrap(~Var1, nrow=2)
dev.off()

ann <- data.frame(cell_type=stringr::str_remove(colnames(cts_norm),  "_mo193|_mo194|_mo196$"))
rownames(ann) <- colnames(cts_norm)




pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig6_2.pdf', width=12)
col_ann_HM <- list(cell_type=cluster_colors[
    c("HSC", "LMPP", "CLP", "Cycling_LMPP", "MEP", "GMP1", "GMP2", "PreB", "Monocytes", "Basophils", "DC","Erythroid")
    ])
plt_hm <- cts_norm
order_col <- colnames(plt_hm)[unlist(
    sapply(c("HSC", "LMPP", "CLP", "MEP", "GMP1", "GMP2", "PreB", "Monocytes", "Basophils", "DC","Erythroid"), 
    function(x) { grep(x, colnames(plt_hm)) }))]
pheatmap::pheatmap(plt_hm[, order_col], annotation=ann, annotation_colors=col_ann_HM,
cluster_cols=FALSE, show_colnames=FALSE, scale='row',
cluster_rows=FALSE)
dev.off()


gsea_list <- readRDS('/home/alsultfh/H1_Project/DATA/gsea_list.rds')
gsea_list_all <- readRDS('/home/alsultfh/H1_Project/DATA/gsea_list_all.rds')
terms_2_keep <- c('GOBP_RNA_PROCESSING', 'GOBP_REGULATION_OF_PROGRAMMED_CELL_DEATH', 'GOBP_CELLULAR_RESPONSE_TO_STRESS', 'GOBP_CELL_MORPHOGENESIS', 'GOBP_MITOTIC_CELL_CYCLE', 'GOBP_APOPTOTIC_PROCESS', 'GOBP_CELL_CYCLE', 'GOBP_CHROMATIN_REMODELING', 'GOBP_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_GENE_EXPRESSION', 'GOBP_REGULATION_OF_GROWTH', 'GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE', 'GOBP_REGULATION_OF_PROTEIN_MODIFICATION_PROCESS', 'GOBP_MYELOID_CELL_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_DEVELOPMENTAL_GROWTH', 
'GOBP_PROTEIN_DNA_COMPLEX_ORGANIZATION', 'GOBP_PROTEIN_DNA_COMPLEX_ASSEMBLY', 
'GOBP_NUCLEOSOME_ORGANIZATION', 'GOBP_CHROMATIN_REMODELING', 
'GOBP_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_TELOMERE_ORGANIZATION', 'GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN', 
'GOBP_PROTEIN_LOCALIZATION_TO_CHROMOSOME', 'GOBP_PROTEIN_LOCALIZATION_TO_CENP_A_CONTAINING_CHROMATIN', 'GOBP_PROTEIN_LOCALIZATION_TO_CHROMOSOME_CENTROMERIC_REGION', 'GOBP_MYELOID_CELL_DIFFERENTIATION', 'GOBP_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION', 'GOBP_CHROMOSOME_ORGANIZATION', 'GOBP_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_MEGAKARYOCYTE_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION', 'GOBP_REGULATION_OF_CELL_DEVELOPMENT', 
'GOBP_REGULATION_OF_HEMOPOIESIS', 'GOBP_CELL_MORPHOGENESIS', 'GOBP_MESENCHYME_DEVELOPMENT', 'GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS', 'GOBP_HOMEOSTATIC_PROCESS', 
'GOBP_HEMOPOIESIS', 'GOBP_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION', 'GOBP_SKELETAL_SYSTEM_DEVELOPMENT', 'GOBP_PROTEIN_RNA_COMPLEX_ORGANIZATION', 'GOBP_REGULATION_OF_IMMUNE_RESPONSE')

terms_2_keep <- c(
  # Group 1: Cell Cycle and Chromatin/Chromosome Organization
  'GOBP_MITOTIC_CELL_CYCLE', 
  'GOBP_CELL_CYCLE', 
  'GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE', 
  'GOBP_PROTEIN_DNA_COMPLEX_ORGANIZATION', 
  'GOBP_PROTEIN_DNA_COMPLEX_ASSEMBLY', 
  'GOBP_NUCLEOSOME_ORGANIZATION', 
  'GOBP_CHROMATIN_REMODELING', 
  'GOBP_TELOMERE_ORGANIZATION', 
  'GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN', 
  'GOBP_PROTEIN_LOCALIZATION_TO_CHROMOSOME', 
  'GOBP_PROTEIN_LOCALIZATION_TO_CENP_A_CONTAINING_CHROMATIN', 
  'GOBP_PROTEIN_LOCALIZATION_TO_CHROMOSOME_CENTROMERIC_REGION', 
  'GOBP_CHROMOSOME_ORGANIZATION',

  # Group 2: Cell Differentiation, Growth, and Development
  'GOBP_CELL_MORPHOGENESIS', 
  'GOBP_DEVELOPMENTAL_GROWTH', 
  'GOBP_REGULATION_OF_CELL_DIFFERENTIATION', 
  'GOBP_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION', 
  'GOBP_MYELOID_CELL_DIFFERENTIATION', 
  'GOBP_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION', 
  'GOBP_MEGAKARYOCYTE_DIFFERENTIATION', 
  'GOBP_NEGATIVE_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION', 
  'GOBP_REGULATION_OF_CELL_DEVELOPMENT', 
  'GOBP_REGULATION_OF_HEMOPOIESIS', 
  'GOBP_MESENCHYME_DEVELOPMENT', 
  'GOBP_HOMEOSTATIC_PROCESS', 
  'GOBP_HEMOPOIESIS', 
  'GOBP_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION', 
  'GOBP_SKELETAL_SYSTEM_DEVELOPMENT',

  # Group 3: Response to Stress and Apoptosis
  'GOBP_APOPTOTIC_PROCESS', 
  'GOBP_CELLULAR_RESPONSE_TO_STRESS', 
  'GOBP_REGULATION_OF_PROGRAMMED_CELL_DEATH',

  # Group 4: Regulation of Immune Response and Protein/RNA Processes
  'GOBP_NEGATIVE_REGULATION_OF_GENE_EXPRESSION', 
  'GOBP_REGULATION_OF_PROTEIN_MODIFICATION_PROCESS', 
  'GOBP_RNA_PROCESSING', 
  'GOBP_PROTEIN_RNA_COMPLEX_ORGANIZATION', 
  'GOBP_REGULATION_OF_IMMUNE_RESPONSE', 
  'GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS'
)




plotter <- data.frame(pathway=NULL, padj=NULL, NES=NULL, HL=NULL)
for (HL in names(gsea_list)){
    tmp <- gsea_list[[HL]]
    tmp <- tmp[tmp$pathway %in% terms_2_keep,  c('pathway', 'padj', 'NES')]
    if(nrow(tmp)==0){next}
    tmp$HL <- stringr::str_extract(HL, '[\\dX]')
    plotter <- rbind(plotter, tmp)
}

for (HL in names(gsea_list)){
    tmp <- gsea_list_all[[HL]]
    tmp <- tmp[tmp$pathway %in% terms_2_keep,  c('pathway', 'padj', 'NES')]
    if(nrow(tmp)==0){next}
    tmp$HL <- stringr::str_extract(HL, '[\\dX]')
    plotter <- rbind(plotter, tmp)
}

plotter$pathway <- sub('GOBP_', '', plotter$pathway)
plotter$pathway <- gsub('_', ' ', str_to_title(plotter$pathway))
terms_2_keep_pretty <- sub('GOBP_', '', terms_2_keep)
terms_2_keep_pretty <- gsub('_', ' ', str_to_title(terms_2_keep_pretty))

plotter$HL <- factor(ifelse(plotter$HL == '3', 'HIST1H1D', 
                ifelse(plotter$HL == 'X', 'H1FX',
                ifelse(plotter$HL == '0', 'H1F0',
                ifelse(plotter$HL == '1', 'HIST1H1A','ERROR')))), 
                levels=c('H1F0', 'HIST1H1A', 'HIST1H1D', 'H1FX'))


plotter$pathway <- factor(plotter$pathway, levels=c(terms_2_keep_pretty))
pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig7.pdf', width=5, height= 8.50)
ggplot() + 
geom_point(plotter, mapping= aes(y=pathway, x=HL, size =padj, fill=NES), color='black', shape=21) + 
scale_fill_gradient2(low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(3, 0.5)) +
guides(fill=guide_colorbar(title='Normalized\nEnrichment\nScore'), size=guide_legend(title='Adjusted\nP.value')) + 
labs(y ='GO Biological process pathways') +
theme(
	legend.position='right',
	axis.text.y=element_text(size=8, color='black'),
	axis.title.x=element_blank(),
	axis.text.x=element_text(size = 8, angle=90, vjust=0.5, hjust=1, color='black'),
	legend.key.size = unit(0.3, "cm"),
	legend.title=element_text(size=9, face='bold'), 
	legend.text=element_text(size=7))

dev.off()



interesting_contrasts <- c("logFC_MGUS_vs_old", "logFC_MM_vs_old", "logFC_SMM_vs_old")
results <- readRDS('/home/alsultfh/H1_Project/DATA/Fatima/FgseaResults_ALL_H1_Contrasts.rds')

mgus_gsea <- as.data.frame(results[['logFC_MGUS_vs_old']])
mgus_gsea <- mgus_gsea[c('pathway', 'padj', 'NES')]
mgus_gsea$Contrast <- 'MGUS'

mm_gsea <- as.data.frame(results[['logFC_MM_vs_old']])
mm_gsea <- mm_gsea[c('pathway', 'padj', 'NES')]
mm_gsea$Contrast <- 'MM'

smm_gsea <- as.data.frame(results[['logFC_SMM_vs_old']])
smm_gsea <- smm_gsea[c('pathway', 'padj', 'NES')]
smm_gsea$Contrast <- 'SMM'

all <- rbind(mgus_gsea, mm_gsea, smm_gsea)
all$Contrast <- factor(all$Contrast, levels=c('MGUS', 'SMM', 'MM'))
all$pathway <- sub('SKO0', 'SKO_H1F0', all$pathway)
all$pathway <- sub('SKO1', 'SKO_HIST1H1A', all$pathway)
all$pathway <- sub('SKO2', 'SKO_HIST1H1C', all$pathway)
all$pathway <- sub('SKO3', 'SKO_HIST1H1D', all$pathway)
all$pathway <- sub('SKO4', 'SKO_HIST1H1E', all$pathway)
all$pathway <- sub('SKO5', 'SKO_HIST1H1B', all$pathway)
all$pathway <- sub('SKOX', 'SKO_H1FX', all$pathway)
all$pathway <- factor(all$pathway, levels = rev(c("SKO_H1FX_vs_WT" , "SKO_H1FX_vs_WT_NEG", "SKO_H1FX_vs_WT_POS", "SKO_HIST1H1D_vs_WT", "SKO_HIST1H1D_vs_WT_NEG" , "SKO_HIST1H1D_vs_WT_POS", "SKO_HIST1H1E_vs_WT", "SKO_HIST1H1E_vs_WT_NEG", "SKO_HIST1H1E_vs_WT_POS" , "SKO_HIST1H1B_vs_WT_NEG", "SKO_HIST1H1B_vs_WT_POS", "SKO_HIST1H1B_vs_WT", "SKO_HIST1H1C_vs_WT" , "SKO_HIST1H1A_vs_WT" , "SKO_H1F0_vs_WT", "SKO_H1F0_vs_WT_POS")))




pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig8.pdf', width=5)
ggplot() + 
geom_point(all, mapping= aes(y=pathway, x=Contrast, size =padj, fill=NES), color='black', shape=21) + 
scale_fill_gradient2(low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   
scale_size_continuous(range = c(6, 0.5), breaks=c(0.005, 0.05, 0.1)) +
guides(fill=guide_colorbar(title='Normalized\nEnrichment\nScore'), size=guide_legend(title='Adjusted\nP.value')) + 
labs(y ='Gene sets') +
theme(
	legend.position='right',
	axis.text.y=element_text(size=8, color='black'),
	axis.title.x=element_blank(),
	axis.text.x=element_text(size = 8, angle=90, vjust=0.5, hjust=1, color='black'),
	legend.key.size = unit(0.3, "cm"),
	legend.title=element_text(size=9, face='bold'), 
	legend.text=element_text(size=7))

dev.off()



gsea_list_SC_SKO_complete <- readRDS('/home/alsultfh/H1_Project/DATA/gsea_list_SC_SKO_complete.rds')
terms_2_keep <- c("SKO0_vs_WT"   ,  "SKO0_vs_WT_POS" ,"SKO0_vs_WT_NEG", "SKO1_vs_WT"  ,   "SKO1_vs_WT_POS","SKO1_vs_WT_NEG", "SKO2_vs_WT"  ,  "SKO2_vs_WT_POS", "SKO2_vs_WT_NEG", "SKO3_vs_WT"    , "SKO3_vs_WT_POS",
 "SKO3_vs_WT_NEG", "SKO4_vs_WT" ,    "SKO4_vs_WT_POS" ,"SKO4_vs_WT_NEG", "SKO5_vs_WT" ,    "SKO5_vs_WT_POS" ,"SKO5_vs_WT_NEG" ,"SKOX_vs_WT"  ,   "SKOX_vs_WT_POS" ,"SKOX_vs_WT_NEG")


plotter <- data.frame(pathway=NULL, padj=NULL, NES=NULL, CellType=NULL)
for (CellType in names(gsea_list_SC_SKO_complete)) {
    tmp <- gsea_list_SC_SKO_complete[[CellType]]
    tmp <- tmp[tmp$pathway %in% terms_2_keep,  c('pathway', 'padj', 'NES')]
    if(nrow(tmp)==0)
    {next}
    tmp$CellType <- CellType
    plotter <- rbind(plotter, tmp)
}


plotter$pathway <- sub('SKO3', 'SKO_HIST1H1D', plotter$pathway)
plotter$pathway <- sub('SKO4', 'SKO_HIST1H1E', plotter$pathway)
plotter$pathway <- sub('SKO5', 'SKO_HIST1H1B', plotter$pathway)
plotter$pathway <- sub('SKOX', 'SKO_H1FX', plotter$pathway)



plotter$CellType <- factor(plotter$CellType, levels=rev(c(
"HSC_Vs_LMPP",
"HSC_Vs_CyclingLMPP",
"HSC_Vs_GMP1",
"HSC_Vs_GMP2",
"HSC_Vs_MEP",
"LMPP_Vs_CLP",
"LMPP_Vs_PreB",
"CyclingLMPP_Vs_CLP",
"CyclingLMPP_Vs_PreB",
"CLP_Vs_PreB",
"MEP_Vs_GMP1",
"MEP_Vs_GMP2",
"MEP_Vs_Erythroid",
"GMP1_Vs_Monocytes",
"GMP1_Vs_Basophils",
"GMP1_Vs_DC",
"GMP1_Vs_MEP",
"GMP2_Vs_Monocytes",
"GMP2_Vs_Basophils",
"GMP2_Vs_DC",
"GMP2_Vs_MEP",
"Monocytes_Vs_Basophils",
"Monocytes_Vs_DC",
"Basophils_Vs_DC"))
)





pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Fig9_2.pdf', width=5)
ggplot() + 
geom_point(plotter, mapping= aes(x=pathway, y=CellType, size =padj, fill=NES), color='black', shape=21) + 
scale_fill_gradient2(low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(5, 0.1), breaks=c(0.05, 0.1, 0.2, 0.3, 0.5)) +
labs(x='Genesets', y='Contrasts')+
guides(fill=guide_colorbar(title='Normalized\nEnrichment\nScore'), size=guide_legend(title='Adjusted\nP.value')) + 
theme(
	legend.position='right',
	axis.text.y=element_text(size=8, color='black'),
	axis.title.x=element_blank(),
	axis.text.x=element_text(size = 8, angle=90, vjust=1, hjust=1, color='black'),
	legend.key.size = unit(0.3, "cm"),
	legend.title=element_text(size=9, face='bold'), 
	legend.text=element_text(size=7))

dev.off()
