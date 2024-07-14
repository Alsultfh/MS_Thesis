library(limma)
library(fgsea)
library(readxl)
library(UpSetR)


fit_cont <- readRDS("/home/alsultfh/H1_Project/DATA/fit_cont.rds")
DEA <- read.table("/home/alsultfh/H1_Project/DATA/Fatima/20230615_dea_results.txt", header = TRUE, sep = "\t", row.names = 1, comment.char = "")
everything_gene_anno <- load("/home/alsultfh/H1_Project/DATA/Fatima/voom_to_dea.RData")



listInput <- list(SKO0_vs_WT = ("SKO0_vs_WT"), SKO1_vs_WT= ("SKO1_vs_WT"), 
                  SKO2_vs_WT= ("SKO2_vs_WT"), SKO3_vs_WT= ("SKO3_vs_WT"), 
                  SKO4_vs_WT=("SKO4_vs_WT"), SKO5_vs_WT= ("SKO4_vs_WT"), 
                  SKOX_vs_WT= ("SKOX_vs_WT"))


upset(fromList(listInput), order.by = "freq", keep.order = TRUE)




library(UpSetR)

coefficients <- c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT')
filtered_results <- list() 

for(coef in coefficients) {
  top_table <- topTable(fit_cont, coef = coef, n= Inf)
  filtered_table <- top_table[top_table$adj.P.Val < 0.05 & abs(top_table$logFC) > 1, ]
  rownames(filtered_table) <- sub('\\.[0-9]*$', '', rownames(filtered_table))
  filtered_table <- merge(filtered_table, gene_anno[, c('gene_id', 'gene_name')], by.x=0, by.y='gene_id', all.x=TRUE)
  filtered_table <- filtered_table[!is.na(filtered_table$gene_name),]
  filtered_results[[coef]] <- filtered_table 
  #coefficients <- c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT')
  # upset(coef, sets = c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT'))
  
}




  
  

selected_data <- DEA[, c("sig_MGUS_vs_old", "sig_MM_vs_old", "sig_SMM_vs_old")]
  
 # Creating the UpSet plot
upset(selected_data, sets = c("sig_MGUS_vs_old", "sig_MM_vs_old", "sig_SMM_vs_old"))

#All genes
selected_data_2 <- selected_data
selected_data_2$gene <- rownames(selected_data_2)
upset_list <- list()
upset_list[["sig_MGUS_vs_old"]] <- selected_data_2[selected_data_2$sig_MGUS_vs_old != 0, 'gene']
upset_list[["sig_MM_vs_old"]]   <- selected_data_2[selected_data_2$sig_MM_vs_old != 0, 'gene']
upset_list[["sig_SMM_vs_old"]]  <- selected_data_2[selected_data_2$sig_SMM_vs_old != 0, 'gene']
upset(fromList(upset_list))

# UPDE genes +SKO 
selected_data_2 <- selected_data
selected_data_2$gene <- rownames(selected_data_2)
upset_list <- list()
upset_list[["sig_MGUS_vs_old"]] <- selected_data_2[selected_data_2$sig_MGUS_vs_old ==1, 'gene']
upset_list[["sig_MM_vs_old"]]   <- selected_data_2[selected_data_2$sig_MM_vs_old   ==1, 'gene']
upset_list[["sig_SMM_vs_old"]]  <- selected_data_2[selected_data_2$sig_SMM_vs_old  ==1, 'gene']
#upset(fromList(upset_list))
for (ct in c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT') ){
  tmp <- topTable(fit_cont, coef = ct, number=Inf)
  tmp <- rownames(tmp[abs(tmp$logFC) >1 & tmp$adj.P.Val < 0.001, ])#0.001 Alex
  tmp <- sub('\\.[0-9]*$', '', tmp, perl=T)
  upset_list[[ct]] <- tmp
}
pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Upset_UpDE.pdf')
upset(fromList(upset_list), sets = names(upset_list), order.by = "freq")
dev.off()

# DOWNDE genes +SKO 
selected_data_2 <- selected_data
selected_data_2$gene <- rownames(selected_data_2)
upset_list <- list()
upset_list[["sig_MGUS_vs_old"]] <- selected_data_2[selected_data_2$sig_MGUS_vs_old ==-1, 'gene']
upset_list[["sig_MM_vs_old"]]   <- selected_data_2[selected_data_2$sig_MM_vs_old   ==-1, 'gene']
upset_list[["sig_SMM_vs_old"]]  <- selected_data_2[selected_data_2$sig_SMM_vs_old  ==-1, 'gene']
#upset(fromList(upset_list))
for (ct in c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT') ){
  tmp <- topTable(fit_cont, coef = ct, number=Inf)
  tmp <- rownames(tmp[abs(tmp$logFC) >1 & tmp$adj.P.Val < 0.001, ])#0.001 Alex
  tmp <- sub('\\.[0-9]*$', '', tmp, perl=T)
  upset_list[[ct]] <- tmp
}
pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Upset_DownDE.pdf')
upset(fromList(upset_list), sets = names(upset_list), order.by = "freq")
dev.off()


#All genes + SKO
selected_data_2 <- selected_data
selected_data_2$gene <- rownames(selected_data_2)
upset_list <- list()
upset_list[["sig_MGUS_vs_old"]] <- selected_data_2[selected_data_2$sig_MGUS_vs_old != 0, 'gene']
upset_list[["sig_MM_vs_old"]]   <- selected_data_2[selected_data_2$sig_MM_vs_old != 0, 'gene']
upset_list[["sig_SMM_vs_old"]]  <- selected_data_2[selected_data_2$sig_SMM_vs_old != 0, 'gene']
# upset(fromList(upset_list))

for (ct in c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT') ){
  tmp <- topTable(fit_cont, coef = ct, number=Inf)
  tmp <- rownames(tmp[abs(tmp$logFC) >1 & tmp$adj.P.Val < 0.001, ])#0.001 Alex
  tmp <- sub('\\.[0-9]*$', '', tmp, perl=T)
  upset_list[[ct]] <- tmp
}
pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/Upset_all_0.001.pdf')
upset(fromList(upset_list), sets = names(upset_list), order.by = "freq")
dev.off()
##checking which genes the are involved in a specific intersection in all genes list 
intersected_genes <- intersect(intersect(intersect(intersect(upset_list[['sig_MGUS_vs_old']],
                  upset_list[['sig_MM_vs_old']]),
                  upset_list[['sig_SMM_vs_old']]),
                  upset_list[['SKO3_vs_WT']]),
                  upset_list[['SKOX_vs_WT']])



#selected_data_fit_cont <- fit_cont[, c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT')]
#selected_data_fit_cont <- as.data.frame(fit_cont$contrasts)

# Creating the UpSet plot
#upset(selected_data_fit_cont, sets = c('SKO0_vs_WT', 'SKO1_vs_WT', 'SKO2_vs_WT', 'SKO3_vs_WT', 'SKO4_vs_WT', 'SKO5_vs_WT', 'SKOX_vs_WT'))






