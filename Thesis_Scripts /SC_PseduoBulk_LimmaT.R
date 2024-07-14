library(Seurat)
library(dplyr)
library(patchwork)
library(DESeq2) 
library(viridis)
library(tidyverse)
library(edgeR)
library(limma)
library(ggplot2)
library(ComplexHeatmap)

# Load the RDS file
# seurat_young <- readRDS("/home/alsultfh/H1_Project/DATA/Dani/seurat_young.rds")
# saveRDS(seurat_young, "/home/alsultfh/H1_Project/DATA/Dani/seurat_young_updated.rds")
seurat_young <- readRDS("/home/alsultfh/H1_Project/DATA/Dani/seurat_young_updated.rds")
fit_limma_voom <- readRDS("/home/alsultfh/H1_Project/DATA/Fatima/fit.rds")
cts_norm <- readRDS("/home/alsultfh/H1_Project/DATA/cts_norm.rds")


##Psedubulking before DE analysis 


# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
metadata <- seurat_young@meta.data


DefaultAssay(seurat_young)


cts <- AggregateExpression(seurat_young, 
                           group.by = c("CellType", "orig.ident"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)


cts <- cts$RNA


################################################### Limma trend with pseudoBulk ######################
logCPM <- cpm(cts, log=TRUE, prior.count = 3)

lmfit <- lmFit(logCPM)
lmfit <- eBayes(lmfit, trend=TRUE)

d0 <- DGEList(cts)
d0 <- calcNormFactors(d0)

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

snames <- colnames(cts)
CellType <- sub('_[\\w]+', '', snames, perl = TRUE)
SampleID <- sub('[\\w]+_', '', snames, perl = TRUE)

data_table <- data.frame(sample_name = paste(snames), CellType = CellType, SampleID = SampleID)
group <- interaction(CellType, SampleID)

mm <- model.matrix(~0 + CellType)
rownames(mm) <- data_table$sample_name


#Voom with trend 

y <- voom(d, mm, plot = TRUE)


fit_trend <- lmFit(y, mm)
fit_trend <- eBayes(fit_trend, trend=TRUE)

topTable(fit_trend, coef=3)



saveRDS(object = fit, file = "fit_trend.rds")

#####################################################Limma trend WITHOUT pseduobulk #############


# Assuming 'condition' is the correct column representing groups
dge <- DGEList(seurat_young@assays$RNA$counts)
dge <- calcNormFactors(dge)
cell_types <- paste(seurat_young$CellType)
design <- model.matrix(~0+cell_types)
rownames(design) <- colnames(seurat_young)
y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- lmFit(y, design = design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
tt <- topTable(fit, n = Inf, adjust.method = "fdr")


test <- FetchData(seurat_young, c('CellType'))

#replot with y$E





#Contrasts defined by lineages order 

contrast.matrix <- makeContrasts(
  # HSC Comparisons
  HSC_Vs_LMPP = cell_typesHSC - cell_typesLMPP,
  HSC_Vs_CyclingLMPP = cell_typesHSC - cell_typesCycling_LMPP,
  HSC_Vs_GMP1 = cell_typesHSC - cell_typesGMP1,
  HSC_Vs_GMP2 = cell_typesHSC - cell_typesGMP2,
  
  # LMPP Comparisons
  LMPP_Vs_CyclingLMPP = cell_typesLMPP - cell_typesCycling_LMPP,
  LMPP_Vs_CLP = cell_typesLMPP - cell_typesCLP,
  LMPP_Vs_PreB = cell_typesLMPP - cell_typesPreB,
 
  # Cycling Comparisons
  CyclingLMPP_Vs_CLP = cell_typesCycling_LMPP - cell_typesCLP,
  CyclingLMPP_Vs_PreB = cell_typesCycling_LMPP - cell_typesPreB,
  
  # CLP Comparisons
  CLP_Vs_PreB = cell_typesCLP - cell_typesPreB,

  # MEP Comparisons
  MEP_Vs_GMP1 = cell_typesMEP - cell_typesGMP1,
  MEP_Vs_GMP2 = cell_typesMEP - cell_typesGMP2,
  MEP_Vs_Monocytes = cell_typesMEP - cell_typesMonocytes,
  MEP_Vs_Basophils = cell_typesMEP - cell_typesBasophils,
  MEP_Vs_DC = cell_typesMEP - cell_typesDC,
  MEP_Vs_Erythroid = cell_typesMEP - cell_typesErythroid,

  # GMP1 Comparisons
  GMP1_Vs_GMP2 = cell_typesGMP1 - cell_typesGMP2,
  GMP1_Vs_Monocytes = cell_typesGMP1 - cell_typesMonocytes,
  GMP1_Vs_Basophils = cell_typesGMP1 - cell_typesBasophils,
  GMP1_Vs_DC = cell_typesGMP1 - cell_typesDC,
  GMP1_Vs_Erythroid = cell_typesGMP1 - cell_typesErythroid,
  GMP1_Vs_MEP = cell_typesGMP1 - cell_typesMEP,
  
  # GMP2 Comparisons
  GMP2_Vs_Monocytes = cell_typesGMP2 - cell_typesMonocytes,
  GMP2_Vs_Basophils = cell_typesGMP2 - cell_typesBasophils,
  GMP2_Vs_DC = cell_typesGMP2 - cell_typesDC,
  GMP2_Vs_Erythroid = cell_typesGMP2 - cell_typesErythroid,
  GMP2_Vs_MEP = cell_typesGMP2 - cell_typesMEP,

  # Monocytes Comparisons
  Monocytes_Vs_Basophils = cell_typesMonocytes - cell_typesBasophils,
  Monocytes_Vs_DC = cell_typesMonocytes - cell_typesDC,
  Monocytes_Vs_Erythroid = cell_typesMonocytes - cell_typesErythroid,
  # Basophils Comparisons
  Basophils_Vs_DC = cell_typesBasophils - cell_typesDC,
  Basophils_Vs_Erythroid = cell_typesBasophils - cell_typesErythroid,

levels = colnames(coef(fit))
)
 

# step 2
  fit <- lmFit(y, design = design)
  fit <- contrasts.fit(fit, contrast.matrix)#linear model lm 
  tmp <-eBayes(fit)
  
saveRDS(tmp, 'limma_trend_results.rds')
## include all DEG genes

results_list <- list()  
for (contrast in colnames(contrast.matrix)) {
    print(contrast)
    top_table <- topTable(tmp, coef=contrast, sort.by="P", n=Inf)
    filtered_table <- top_table[top_table$adj.P.Val < 0.05 & abs(top_table$logFC) > 0.2 ,]    
    # Store the result in the list
    results_list[[contrast]] <- filtered_table
}



# Define the list of histone linker symbols
histone_linkers_symbol <- c("H1F0", "HIST1H1A", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H1B", "H1FX")

# Initialize an empty list to store results
results_list <- list() 

# Loop through each contrast in your contrast matrix
for (contrast in colnames(contrast.matrix)) {
    print(contrast)
    # Perform DE analysis and obtain the top table of results
    top_table <- topTable(tmp, coef=contrast, sort.by="P", n=Inf)
    
    # Filter results based on adjusted P-value and log fold change
    filtered_table <- top_table[top_table$adj.P.Val < 0.05 & abs(top_table$logFC) > 0.1,]    
    
    # Further filter the table to include only rows corresponding to histone linker genes
    histone_linker_DE <- filtered_table[rownames(filtered_table) %in% histone_linkers_symbol,]
    
    # Store the filtered result in the list
    results_list[[contrast]] = histone_linker_DE 
}

print(results_list)

library(writexl)
write_xlsx(results_list, "results_list_LimmaT.xlsx")

###LIMMA_T BOX plot Thesis 



gene_2_plot <- c("H1F0", "HIST1H1A", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H1B", "H1FX")

# First, check if the genes to plot are in the rownames and filter cts_norm accordingly
filtered_data <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
filtered_data <- filtered_data[gene_2_plot, ]
filtered_data <- t(filtered_data)
annotation <- FetchData(seurat_young, vars=c('CellType'))

filtered_data <- merge(filtered_data, annotation, by=0)

tmp <- reshape2::melt(filtered_data)
pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/TLimma_T.pdf')
ggplot(tmp, aes(x=CellType, y=value, fill=CellType))+
geom_boxplot(outlier.size=1, outlier.alpha=0.1) +
theme_bw() +
facet_wrap(~variable)
dev.off()

 #First, convert to long format and remove suffixes
long_data <- as.data.frame(filtered_data) %>%
  tibble::rownames_to_column('Gene') %>%
  pivot_longer(cols = -Gene, names_to = "CellType", values_to = "Expression") %>%
  mutate(CellType = str_remove(CellType, "_mo193|_mo194|_mo196$"))
# Check the unique CellTypes after modification to ensure the removal worked
print(unique(long_data$CellType))

#  now filter
desired_order <- c("HSC", "LMPP", "Cycling_LMPP","CLP", "PreB", "MEP","Erythroid", "GMP1", "GMP2", "Monocytes", "Basophils", "DC")
long_data_filtered <- long_data %>%
  filter(CellType %in% desired_order)

# Check if any rows are left after filtering
print(nrow(long_data_filtered))



long_data_filtered$CellType <- factor(long_data_filtered$CellType, levels = desired_order)

# Get unique list of genes from the data
genes <- unique(long_data_filtered$Gene)


# Iterate over each gene to create and save a plot
for(gene in genes) {
  # Filter data for the current gene
  gene_data <- filter(long_data_filtered, Gene == gene)
  
  # Create the plot for the current gene
  p <- ggplot(gene_data, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_boxplot() + geom_point() +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), # Vertical x-axis labels for readability
          strip.text.x = element_text(face = "bold")) + # Bold gene names for clarity
    labs(title = paste("Gene Expression of", gene, "by Cell Type"), 
         x = "Cell Type", 
         y = "Limma T w/o Pseudobulk")
  
  # Save the plot with larger dimensions
  # Define the file name based on the gene, adjusting the path as needed
  file_name <- paste0("/home/alsultfh/H1_Project/DATA/Fatima/Plots/Limma_Trend_plot_", gene, ".png")
  ggsave(file_name, plot = p, width = 10, height = 8, dpi = 300)
}











