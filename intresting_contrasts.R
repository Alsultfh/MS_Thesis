library(edgeR)
library(limma)
library(Seurat)


seurat_young <- readRDS("/home/alsultfh/H1_Project/DATA/Dani/seurat_young_updated.rds")
cts_norm <- readRDS("/home/alsultfh/H1_Project/DATA/cts_norm.rds")
fit <- readRDS('/home/alsultfh/H1_Project/DATA/Fatima/fit.rds')
custom_pathways <- fgsea::gmtPathways("/home/alsultfh/H1_Project/DATA/Fatima/SKO_patwhays.gmt")
everything_gene_anno <- load("/home/alsultfh/H1_Project/DATA/Fatima/voom_to_dea.RData")


contrast.matrix <- makeContrasts(
  # HSC Comparisons
  HSC_Vs_LMPP = CellTypeHSC - CellTypeLMPP,
  HSC_Vs_CyclingLMPP = CellTypeHSC - CellTypeCycling,
  HSC_Vs_GMP1 = CellTypeHSC - CellTypeGMP1,
  HSC_Vs_GMP2 = CellTypeHSC - CellTypeGMP2,
  HSC_Vs_MEP = CellTypeHSC - CellTypeMEP,
  # LMPP Comparisons
  LMPP_Vs_CyclingLMPP = CellTypeLMPP - CellTypeCycling,
  LMPP_Vs_CLP = CellTypeLMPP - CellTypeCLP,
  LMPP_Vs_PreB = CellTypeLMPP - CellTypePreB,
 
  # Cycling Comparisons
  CyclingLMPP_Vs_CLP = CellTypeCycling - CellTypeCLP,
  CyclingLMPP_Vs_PreB = CellTypeCycling - CellTypePreB,
  
  # CLP Comparisons
  CLP_Vs_PreB = CellTypeCLP - CellTypePreB,

  # MEP Comparisons
  MEP_Vs_GMP1 = CellTypeMEP - CellTypeGMP1,
  MEP_Vs_GMP2 = CellTypeMEP - CellTypeGMP2,
  MEP_Vs_Erythroid = CellTypeMEP - CellTypeErythroid,

  # GMP1 Comparisons
  GMP1_Vs_GMP2 = CellTypeGMP1 - CellTypeGMP2,
  GMP1_Vs_Monocytes = CellTypeGMP1 - CellTypeMonocytes,
  GMP1_Vs_Basophils = CellTypeGMP1 - CellTypeBasophils,
  GMP1_Vs_DC = CellTypeGMP1 - CellTypeDC,
  GMP1_Vs_MEP = CellTypeGMP1 - CellTypeMEP,
  
  # GMP2 Comparisons
  GMP2_Vs_Monocytes = CellTypeGMP2 - CellTypeMonocytes,
  GMP2_Vs_Basophils = CellTypeGMP2 - CellTypeBasophils,
  GMP2_Vs_DC = CellTypeGMP2 - CellTypeDC,
  GMP2_Vs_MEP = CellTypeGMP2 - CellTypeMEP,

  # Monocytes Comparisons
  Monocytes_Vs_Basophils = CellTypeMonocytes - CellTypeBasophils,
  Monocytes_Vs_DC = CellTypeMonocytes - CellTypeDC,
  # Basophils Comparisons
  Basophils_Vs_DC = CellTypeBasophils - CellTypeDC,

  levels = colnames(coef(fit))
  )
 


# step 2
  tmp <- contrasts.fit(fit, contrast.matrix)#linear model lm 
  tmp <- eBayes(tmp)#empirical Bayes smoothing to the standard errors of the estimated coefficients in the linear model fits

saveRDS(tmp, './Results_DE_updated_Lineage.rds')

results_de <- readRDS('./Results_DE_updated_Lineage.rds')

#26 contrasts


interesting_contrasts <- c('HSC_Vs_LMPP','HSC_Vs_CyclingLMPP','HSC_Vs_GMP1','HSC_Vs_GMP2','HSC_Vs_MEP','LMPP_Vs_CyclingLMPP',
              'LMPP_Vs_CLP','LMPP_Vs_PreB','CyclingLMPP_Vs_CLP', 'CyclingLMPP_Vs_PreB', 'CLP_Vs_PreB',
              'MEP_Vs_GMP1', 'MEP_Vs_GMP2','MEP_Vs_Monocytes','MEP_Vs_Basophils','MEP_Vs_DC','MEP_Vs_Erythroid',
              'GMP1_Vs_GMP2','GMP1_Vs_Monocytes', 'GMP1_Vs_Basophils','GMP1_Vs_DC','GMP1_Vs_MEP',
              'GMP2_Vs_Monocytes','GMP2_Vs_Basophils','GMP2_Vs_DC','GMP2_Vs_MEP','Monocytes_Vs_Basophils',
              'Monocytes_Vs_DC','Basophils_Vs_DC')



results <- list()
for (contrast in colnames(results_de$coefficients)){
  print(contrast)
    tmp <- topTable(results_de, coef = contrast, n=Inf, adjust='BH')
#   tmp <- tmp[!duplicated(tmp$gene_name),]
#   tmp <- tmp[order(tmp[[contrast]]), ]
    tmp <- tmp[tmp$adj.P.Val<0.05 & abs(tmp$logFC) > 1, ]
    tmp <- tmp[order(tmp$logFC, decreasing=T),]
    tmp$gene_name <- rownames(tmp)
    fgsea_input <- setNames(tmp$logFC, c(tmp$gene_name))
    fgseaRes <- fgsea::fgsea(pathways = custom_pathways, 
                            stats = fgsea_input, minSize = 10, maxSize = 500, 
                            nPermSimple = 1000)
    results[[contrast]] <- fgseaRes
}

# saveRDS(results,'/home/alsultfh/H1_Project/DATA/gsea_list_SC_SKO_complete.rds' )
results <- readRDS('/home/alsultfh/H1_Project/DATA/gsea_list_SC_SKO_complete.rds' )