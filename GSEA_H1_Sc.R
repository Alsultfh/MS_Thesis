library(limma)
library(fgsea)
library(ggplot2)

everything_gene_anno <- load("/home/alsultfh/H1_Project/DATA/Fatima/voom_to_dea.RData")

seurat_young <- readRDS("/home/alsultfh/H1_Project/DATA/Dani/seurat_young_updated.rds")
fit_cont <- readRDS("/home/alsultfh/H1_Project/DATA/fit_cont.rds")
pathway <- fgsea::gmtPathways('/home/alsultfh/H1_Project/c5.go.bp.v2023.2.Hs.symbols.gmt')


coefficients <- c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT')

gsea_list <- list()
for (coef in coefficients){
    tmp <- topTable(fit_cont, coef=coef, n=Inf, adjust='BH')
    tmp <-  tmp[tmp$adj.P.Val < 0.05 & abs(tmp$logFC) > 1, ]
    tmp$gene_id <-  sub('\\.[0-9]*$', '',rownames(tmp))
    tmp <- merge(tmp , gene_anno[, c('gene_id', 'gene_name')],  by='gene_id', all.x=TRUE)
    tmp <- tmp[!is.na(tmp$gene_name),]
    tmp <- tmp[order(tmp$logFC, decreasing=T),]
    fgsea_input <- setNames(tmp$logFC, c(tmp$gene_name))
    fgseaRes <- as.data.frame(fgsea::fgsea(pathways = pathway, 
                           stats = fgsea_input, minSize = 10, maxSize = 500, 
                           nPermSimple = 1000))
    gsea_list[[coef]] <- fgseaRes[fgseaRes$padj < 0.05, ]
}
saveRDS(gsea_list,'/home/alsultfh/H1_Project/DATA/gsea_list.rds')
#writexl::write_xlsx(gsea_list, "/home/alsultfh/H1_Project/DATA/GSEA_results.xlsx")


gsea_list <- list()
for (coef in coefficients){
    tmp <- topTable(fit_cont, coef=coef, n=Inf, adjust='BH')
    tmp <-  tmp[tmp$adj.P.Val < 0.5 & abs(tmp$logFC) > 1, ]
    tmp$gene_id <-  sub('\\.[0-9]*$', '',rownames(tmp))
    tmp <- merge(tmp , gene_anno[, c('gene_id', 'gene_name')],  by='gene_id', all.x=TRUE)
    tmp <- tmp[!is.na(tmp$gene_name),]
    tmp <- tmp[order(tmp$logFC, decreasing=T),]
    fgsea_input <- setNames(tmp$logFC, c(tmp$gene_name))
    fgseaRes <- as.data.frame(fgsea::fgsea(pathways = pathway, 
                           stats = fgsea_input, minSize = 10, maxSize = 500, 
                           nPermSimple = 1000))
    gsea_list[[coef]] <- fgseaRes[fgseaRes$padj < 0.05, ]
}
saveRDS(gsea_list,'/home/alsultfh/H1_Project/DATA/gsea_list_all.rds')

#writexl::write_xlsx(gsea_list, "/home/alsultfh/H1_Project/DATA/GSEA_results_ALL.xlsx")

terms_2_keep <- c('GOBP_RNA_PROCESSING', 'GOBP_REGULATION_OF_PROGRAMMED_CELL_DEATH', 'GOBP_CELLULAR_RESPONSE_TO_STRESS', 'GOBP_CELL_MORPHOGENESIS', 'GOBP_MITOTIC_CELL_CYCLE', 'GOBP_APOPTOTIC_PROCESS', 'GOBP_CELL_CYCLE', 'GOBP_CHROMATIN_REMODELING', 'GOBP_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_GENE_EXPRESSION', 'GOBP_REGULATION_OF_GROWTH', 'GOBP_REGULATION_OF_MITOTIC_CELL_CYCLE', 'GOBP_REGULATION_OF_PROTEIN_MODIFICATION_PROCESS', 'GOBP_MYELOID_CELL_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION', 'GOBP_DEVELOPMENTAL_GROWTH', 'GOBP_PROTEIN_DNA_COMPLEX_ORGANIZATION', 'GOBP_PROTEIN_DNA_COMPLEX_ASSEMBLY', 'GOBP_NUCLEOSOME_ORGANIZATION', 'GOBP_CHROMATIN_REMODELING', 'GOBP_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_TELOMERE_ORGANIZATION', 'GOBP_PROTEIN_LOCALIZATION_TO_CHROMATIN', 'GOBP_PROTEIN_LOCALIZATION_TO_CHROMOSOME', 'GOBP_PROTEIN_LOCALIZATION_TO_CENP_A_CONTAINING_CHROMATIN', 'GOBP_PROTEIN_LOCALIZATION_TO_CHROMOSOME_CENTROMERIC_REGION', 'GOBP_MYELOID_CELL_DIFFERENTIATION', 'GOBP_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION', 'GOBP_CHROMOSOME_ORGANIZATION', 'GOBP_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION', 'GOBP_MEGAKARYOCYTE_DIFFERENTIATION', 'GOBP_NEGATIVE_REGULATION_OF_MYELOID_CELL_DIFFERENTIATION', 'GOBP_REGULATION_OF_CELL_DEVELOPMENT', 'GOBP_REGULATION_OF_HEMOPOIESIS', 'GOBP_CELL_MORPHOGENESIS', 'GOBP_MESENCHYME_DEVELOPMENT', 'GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS', 'GOBP_HOMEOSTATIC_PROCESS', 'GOBP_HEMOPOIESIS', 'GOBP_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION', 'GOBP_SKELETAL_SYSTEM_DEVELOPMENT', 'GOBP_PROTEIN_RNA_COMPLEX_ORGANIZATION', 'GOBP_REGULATION_OF_IMMUNE_RESPONSE')


gsea_list <- readRDS('/home/alsultfh/H1_Project/DATA/gsea_list.rds')
gsea_list_all <- readRDS('/home/alsultfh/H1_Project/DATA/gsea_list_all.rds')

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



pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/GSEA_HL.pdf')
ggplot() + 
geom_point(plotter, mapping= aes(y=pathway, x=HL, size =padj, fill=NES), color='black', shape=21) + 
scale_fill_gradient2(low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(3, 0.5)) +
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
