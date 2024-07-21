##GSEA##

##set working dir


## libraries 
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggsignif)

##loeading the Data 
load("/home/alsultfh/H1_Project/DATA/Fatima/voom_to_dea.RData")
## visualiz Data 
head(gene_anno)

DEA <- read.table("/home/alsultfh/H1_Project/DATA/Fatima/20230615_dea_results.txt", header = TRUE, sep = "\t", row.names = 1, comment.char = "")
#pathway <- read.table("~/Desktop/H1/fgseaRes.txt", header = TRUE, sep = "\t", comment.char = "")
pathway <- fgsea::gmtPathways('/Users/Alsultfh/Desktop/H1/Fatima/c5.go.bp.v2023.2.Hs.symbols.gmt')
#Anodata ??
#merge ...

DEA <- merge(DEA, gene_anno[, c('gene_id', 'gene_name')], by.x=0, by.y='gene_id', all.x=TRUE)

      
      
interesting_contrasts <- c("logFC_MGUS_vs_old", "logFC_MM_vs_old", "logFC_SMM_vs_old", 
                           "logFC_MM_vs_SMM", "logFC_SMM_vs_MGUS", "logFC_MM_vs_MGUS")
results <- list()

for (contrast in interesting_contrasts) {
  print(contrast)
  sig_col <- sub('logFC', 'sig', contrast)
  #tmp <- DEA %>% filter(!!as.name(sig_col) != 0) %>%
    #arrange(!!as.name(contrast))
  tmp <- DEA[DEA[[sig_col]] !=0 , c(contrast, 'gene_name')]
  tmp <- tmp[!duplicated(tmp$gene_name),]
  #We need to sort the input of GSEA by logFC
  tmp <- tmp[order(tmp[[contrast]]), ]
  fgsea_input <- setNames(tmp[[contrast]], c(tmp$gene_name))
  fgseaRes <- fgsea::fgsea(pathways = pathway, stats = fgsea_input, minSize = 10, maxSize = 500)
  results[[contrast]] <- fgseaRes
}

##

saveRDS(results, '/Users/Alsultfh/Desktop/H1/Fatimah/FgseaResults_ALL_Contrasts.rds')

for (contrast in interesting_contrasts) {
  pdf(paste0('/Users/Alsultfh/Desktop/H1/Fatimah/TopBottomGSEA_',contrast,'.pdf'))
  print(contrast)
  fgseaRes <- results[[contrast]] 
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  print(plotGseaTable(pathway[topPathways], fgsea_input, fgseaRes, 
                gseaParam=0.5))
  dev.off()
}


###### Expression analysis 

#load("~/Desktop/H1/Fatimah/voom_to_dea.RData")


#metadata <- read.table("/home/alsultfh/metadata.txt", sep="\t", header = TRUE)
#countData <- read.table("/home/alsultfh/norm_data.txt",sep="\t",header=T)
countData <- readRDS("/home/alsultfh/H1_Project/DATA/norm_data.rds")

#gene_anno <- readRDS("/home/alsultfh/H1_Project/DATA/gene_anno.rds")
countData <- norm_data

gene_2_plot <- c("H1-10", "H1-1"   ,   "H1-2"   ,   "H1-4" , 
                 "H1-3"   ,   "H1-5"  ,   "H1-0")

gene_2_plot %in% gene_anno$gene_name

ids_2_plot <- gene_anno[gene_anno$gene_name %in% gene_2_plot, 'gene_id']

plotter <- countData[ids_2_plot, ]
plotter_mlt <- reshape2::melt(plotter)
plotter_mlt<- setNames(plotter_mlt, c('Gene_ID', 'Sample', 'Expression'))
plotter_mlt$type <- metadata[match(plotter_mlt$Sample, metadata$Sample),'group']

plotter_mlt <- merge(plotter_mlt, gene_anno, by.x='Gene_ID', by.y='gene_id', all.x=TRUE)
plotter_mlt$type <- factor(plotter_mlt$type, levels = c("healthy donors","MGUS","SMM","MM"))
library(ggplot2)

##wilcox 
plotter_mlt <- merge(plotter_mlt, DEA[, c('Row.names', 'p.adj_MGUS_vs_old',  
                                            'p.adj_MM_vs_old', 'p.adj_SMM_vs_old')], 
                     by.x='Gene_ID', by.y='Row.names')

pdf('Gene_expression_per_Type.pdf')
ggplot(plotter_mlt, aes(x=type, y= Expression, fill=type)) + 
  geom_boxplot() + theme_classic() + 
  geom_signif(
    comparisons = list(c("healthy donors", "MGUS"),
                       c("healthy donors","SMM"),
                       c("healthy donors","MM")),
    step_increase = 0.2,
    map_signif_level=TRUE,
    show.legend = TRUE) + ylim(-6, 25) + 
  theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1)) +
  facet_wrap(~gene_name)
dev.off()
### Pval from limma 

plot_gene_boxplot <- function(data, gene, pvalMgus, pvalSMM, pvalMM){
  data <- data[data$gene_name == gene, ]
  p <- ggplot(data, 
              aes(x=type, y= Expression, fill=type)) + 
    geom_boxplot() + theme_classic() + 
    geom_signif(
      comparisons = list(c("healthy donors", "MGUS"),
                         c("healthy donors","SMM"),
                         c("healthy donors","MM")),
                    annotations = c(pvalMgus, 
                                    pvalSMM,
                                    pvalMM),
      textsize = 1.75,
      step_increase = 0.2,
      map_signif_level=TRUE,
      show.legend = TRUE) + ylim(-6, max(data$Expression) +6) + ggtitle(gene) + 
    theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1), 
          axis.title.x = element_blank()) 
    return(p)
}


plot_list = list()
for (gene_name in gtools::mixedsort(unique(plotter_mlt$gene_name), decreasing = TRUE)) {
  print(gene_name)
  plot_list[[gene_name]] <- plot_gene_boxplot(plotter_mlt, gene_name, 
                                              format(signif(unique(plotter_mlt[plotter_mlt$type == 'MGUS' & 
                                                                         plotter_mlt$gene_name == gene_name, 
                                                                       'p.adj_MGUS_vs_old']), 3), scientific=T), 
                                              format(signif(unique(plotter_mlt[plotter_mlt$type == 'SMM' & 
                                                                         plotter_mlt$gene_name == gene_name, 
                                                                       'p.adj_SMM_vs_old']), 3), scientific=T), 
                                              format(signif(unique(plotter_mlt[plotter_mlt$type == 'MM' & 
                                                                         plotter_mlt$gene_name == gene_name, 
                                                                       'p.adj_MM_vs_old']), 3), scientific=T))
}


pdf('Gene_expression_per_Type_2.pdf')
cowplot::plot_grid(plotlist = plot_list, nrow=3, ncol=3)
dev.off()



#Step1:install packages 
# BiocManager::install("ComplexHeatmap")

# Load the ComplexHeatmap package
library(ComplexHeatmap)
plotter <- countData[ids_2_plot, ]
data_matrix <- matrix(rnorm(200), nrow = 8)

rownames(data_matrix) <- paste("Sample")
colnames(data_matrix) <- paste("Feature")

# Create the heatmap for Histone expression in multiple myeloma
pdf('Heatmap.pdf')
annotation_hm <- data.frame(Sample=colnames(plotter))
annotation_hm <- merge(annotation_hm, metadata, by='Sample', all.x=TRUE)
column_ha <- HeatmapAnnotation(Type = annotation_hm$group, 
                               col = list(Type = c("MM" = "red", "healthy donors" = "green", "SMM" = "blue", 'MGUS'='pink')))


plotter_hm <- merge(plotter, gene_anno[, c('gene_id', 'gene_name')], by.x=0, by.y='gene_id', all.x=TRUE)
rownames(plotter_hm) <- plotter_hm$gene_name
plotter_hm <- plotter_hm[, !colnames(plotter_hm) %in% c('Row.names', 'gene_name')] 
#tmp <- t(scale(t(plotter)))
Heatmap(as.matrix(plotter_hm),
        name = "My Heatmap",
        row_names_side = "left",
        column_names_side = "top",
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        show_row_names = TRUE,
        show_column_names = FALSE,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        top_annotation = column_ha,
        col = circlize::colorRamp2(c(min(as.matrix(plotter_hm)), 0, max(as.matrix(plotter_hm))), 
                         c("blue", "white", "red"))
)
dev.off()

###


library(readxl)

head(SKO_H1)
DKO_H1 <- read_excel("~/Desktop/H1/Sumaya/DEG_genes_DKO_WT_contrasts.xls")
##To do 
#transfer name_ID to gene_ID 
#7 coustome pathway for each subtype H1-7
#making custome gmt file 
#biological process 
#H1-depletion and so on .. 

#filter by adj.value +logFC




SKO_H1 <- read_excel("~/Desktop/H1/Sumaya/DEG_genes_SKO_WT_contrasts.xls")
pathway <- fgsea::gmtPathways('/Users/Alsultfh/Desktop/H1/Fatimah/c5.go.bp.v2023.2.Hs.symbols.gmt')

#subset_pathway?

##filtering P.adj and logFC
library(limma)

coefficients <- c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT')

filtered_results <- list()

for(coef in coefficients) {
  top_table <- topTable(fit_cont, coef = coef, n= Inf)
  filtered_table <- top_table[top_table$adj.P.Val < 0.05 & abs(top_table$logFC) > 1, ]
  
  filtered_results[[coef]] <- filtered_table
}

??

topTable(fit_cont, coef = 'SKO0_vs_WT', n= Inf)
topTable(fit_cont, coef = 'SKO1_vs_WT', n= Inf)
topTable(fit_cont, coef = 'SKO2_vs_WT', n= Inf)
topTable(fit_cont, coef = 'SKO3_vs_WT', n= Inf)
topTable(fit_cont, coef = 'SKO4_vs_WT', n= Inf)
topTable(fit_cont, coef = 'SKO5_vs_WT', n= Inf)
topTable(fit_cont, coef = 'SKOX_vs_WT', n= Inf)

##perform custome gmt file and for each list below use one of the subtype 
gs.list <- list(
  "SKO0_vs_WT" = ids2indices(geneIds(getGmt("genesets/disorder.gmt")), id = v$genes$Gene_name),
  "SKO1_vs_WT" = ids2indices(geneIds(getGmt("genesets/behavior.gmt")), id=v$genes$Gene_name),
  "SKO2_vs_WT" = ids2indices(geneIds(getGmt("genesets/head.gmt")), id=v$genes$Gene_name),
  "SKO3_vs_WT" = ids2indices(geneIds(getGmt("genesets/neural.gmt")), id=v$genes$Gene_name),
  "SKO4_vs_WT" = ids2indices(geneIds(getGmt("genesets/postsynapse.gmt")), id=v$genes$Gene_name),
  "SKO5_vs_WT" = ids2indices(geneIds(getGmt("genesets/presynapse.gmt")), id=v$genes$Gene_name),
  "SKOX_vs_WT" = ids2indices(geneIds(getGmt("genesets/synapse.gmt")), id=v$genes$Gene_name))




