library(limma)
library(fgsea)
library(readxl)



fit_cont <- readRDS("~/Desktop/H1/Sumaya/fit_cont.rds")
DEA <- read.table("~/Desktop/H1/Alexsandra/20230615_dea_results.txt", header = TRUE, sep = "\t", row.names = 1, comment.char = "")
voom



everything_gene_anno <- load("~/Desktop/H1/Alexsandra/voom_to_dea.RData")

SKO_H1 <- read_excel("~/Desktop/H1/Sumaya/DEG_genes_SKO_WT_contrasts.xls")
pathway <- fgsea::gmtPathways('/Users/Alsultfh/Desktop/H1/Alexsandra/c5.go.bp.v2023.2.Hs.symbols.gmt')

#subset_pathway?

##filtering P.adj and logFC


coefficients <- c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT')

filtered_results <- list()



# system('rm ~/Desktop/H1/Alexsandra/SKO_patwhays.gmt')#to destroy the file when running it so many time 
# for(coef in coefficients) {
#   top_table <- topTable(fit_cont, coef = coef, n= Inf)
#   filtered_table <- top_table[top_table$adj.P.Val < 0.05 & abs(top_table$logFC) > 1, ]
#   #filtered_results[[coef]] <- filtered_table
#   rownames(filtered_table) <- sub('\\.[0-9]*$', '', c(rownames(filtered_table)))
#   filtered_table <- merge(filtered_table , gene_anno[, c('gene_id', 'gene_name')], by.x=0, by.y='gene_id', all.x=TRUE)
#   filtered_table <- filtered_table[!is.na(filtered_table$gene_name),]
#   name_pathway <- coef
#   description <- 'none'
#   genes_list <- paste(filtered_table$gene_name, collapse = '\t')
#   final_gmt <- paste0(name_pathway, '\t', description ,'\t', genes_list)
#   write(final_gmt,file="~/Desktop/H1/Alexsandra/SKO_patwhays.gmt",append=TRUE)
#   final_gmt <- paste0( coef, '_POS', '\t', 
#                        description ,'\t', 
#                        paste(filtered_table[filtered_table$logFC > 0, 'gene_name'], collapse = '\t'))
#   write(final_gmt,file="~/Desktop/H1/Alexsandra/SKO_patwhays.gmt",append=TRUE)
#   final_gmt <- paste0( coef, '_NEG', '\t', 
#                        description ,'\t', 
#                        paste(filtered_table[filtered_table$logFC < 0, 'gene_name'], collapse = '\t'))
#   write(final_gmt,file="~/Desktop/H1/Alexsandra/SKO_patwhays.gmt",append=TRUE)
# }

## GSEA for each comparison (Healthy ,MM, MGUS, SMM) and H1 linker

custom_pathways <- fgsea::gmtPathways("~/Desktop/H1/Alexsandra/SKO_patwhays.gmt")
DEA <- read.table("~/Desktop/H1/Alexsandra/20230615_dea_results.txt", header = TRUE, sep = "\t", row.names = 1, comment.char = "")

DEA <- merge(DEA, gene_anno[, c('gene_id', 'gene_name')], by.x=0, by.y='gene_id', all.x=TRUE)


interesting_contrasts <- c("logFC_MGUS_vs_old", "logFC_MM_vs_old", "logFC_SMM_vs_old")
                          
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
  fgseaRes <- fgsea::fgsea(pathways = custom_pathways, 
                           stats = fgsea_input, minSize = 10, maxSize = 500, 
                           nPermSimple = 10000000)
  results[[contrast]] <- fgseaRes
}
#write.csv(results,'/Users/Alsultfh/Desktop/H1/Alexsandra/FgseaResults_ALL_H1_Contrasts.csv')

#the split for each knockout pathway by logFC

saveRDS(results, '/Users/Alsultfh/Desktop/H1/Alexsandra/FgseaResults_ALL_H1_Contrasts.rds')

All_H1_Contrasts<- readRDS('/Users/Alsultfh/Desktop/H1/Alexsandra/FgseaResults_ALL_H1_Contrasts.rds')

for (contrast in interesting_contrasts) {
  pdf(paste0('/Users/Alsultfh/Desktop/H1/Alexsandra/TopBottomGSEAH1_',contrast,'.pdf'))
  print(contrast)
  fgseaRes <- results[[contrast]] 
  topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[NES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  print(plotGseaTable(custom_pathways[topPathways], fgsea_input, fgseaRes, 
                      gseaParam=0.5))
  dev.off()
}

##################################### STOP STOP STOP STOP HERE ###################################################


##upset plot(MM,MGUS,SMM VS old)+ SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT'
#SKO_H1
#DEA
library(UpSetR)

coefficients <- c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT')
filtered_results <- list() 

for(coef in coefficients) {
  top_table <- topTable(fit_cont, coef = coef, n= Inf)
  filtered_table <- top_table[top_table$adj.P.Val < 0.05 & abs(top_table$logFC) > 1, ]
  rownames(filtered_table) <- sub('\\.[0-9]*$', '', rownames(filtered_table), perl=T)
  filtered_table <- merge(filtered_table, gene_anno[, c('gene_id', 'gene_name')], by.x=0, by.y='gene_id', all.x=TRUE)
  filtered_table <- filtered_table[!is.na(filtered_table$gene_name),]
  filtered_results[[coef]] <- filtered_table 
}



###

# example of list input (list of named vectors)

listInput <- list(SKO0_vs_WT = ("SKO0_vs_WT"), SKO1_vs_WT= ("SKO1_vs_WT"), 
                  SKO2_vs_WT= ("SKO2_vs_WT"), SKO3_vs_WT= ("SKO3_vs_WT"), 
                  SKO4_vs_WT=("SKO4_vs_WT"), SKO5_vs_WT= ("SKO4_vs_WT"), 
                  SKOX_vs_WT= ("SKOX_vs_WT"))

upset(fromList(listInput), order.by = "freq", keep.order = TRUE)


###

# Load the necessary library
library(readr)

# Read the RDS file
data1 <- readRDS("~/Desktop/H1/Sumaya/fit_cont.rds")

# Write the data to a CSV file

write_csv(data, "output_fit_cont.csv")









library(UpSetR)

# Ensure the data is in the correct format: a binary matrix
# The first column of upset_data should not be included in the UpSet plot
binary_matrix <- as.data.frame(upset_data[, -1], stringsAsFactors=F)

# Create the UpSet plot
upset(fromList(binary_matrix), 
      sets = c("SKO0_vs_WT", "SKO1_vs_WT", "SKO2_vs_WT", "SKO3_vs_WT", "SKO4_vs_WT", "SKO5_vs_WT", "SKOX_vs_WT"),
      number.angles = 30, 
      point.size = 3.5, 
      line.size = 2, 
      main.bar.color = "gray60", 
      matrix.color="gray60")




















#trash code
# remove genes without annotation
#ids <- c("ENSG00000259820", "ENSG00000273576")
#filtered_table <- filtered_table[!rownames(filtered_table)%in%ids,]
#filtered_table <- filtered_table[rownames(filtered_table)%in%biotypes$ensembl_gene_id,]
#rownames(filtered_table) <- biotypes[biotypes$ensembl_gene_id%in%rownames(filtered_table),]$hgnc_symbol

##converting gene ID to gene name in a different way 
library(biomaRt)


mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes_ids <- sub('\\.[0-9]*$', '', c(rownames(fit_cont)))

biotypes = getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position',
                                 'strand','gene_biotype','percentage_gene_gc_content'),
                 filters = 'ensembl_gene_id', values=genes_ids, mart=mart) #19508     8


#Let's check duplicates
which(duplicated(biotypes$ensembl_gene_id)) #no duplicates

sum(!(genes_ids %in% biotypes$ensembl_gene_id))
#there are 2 genes from fit_cont without annotation
#Check what type of ensembl ID are...
ids <- c(genes_ids[!genes_ids %in% biotypes$ensembl_gene_id])
# "ENSG00000259820" "ENSG00000273576"


write.table(biotypes,'gene_annotation.txt', sep='\t', dec='.', quote=FALSE, row.names=FALSE, col.names=TRUE)




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

