
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
library(ggsignif)
library(viridis)

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
cts_norm <- cpm(cts, log=TRUE)

############################################################Limma Voom##########################################################

head(cts)
#Create DGEList object


#2#Preprocessing

d0 <- DGEList(cts)

d0

d0 <- calcNormFactors(d0)
d0
dim(d0)

#Filter low-expressed genes
#19812

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left   #17105 #36

#17105 Number of gene left 

snames <- colnames(cts) # Sample names
snames


CellType <- sub('_[\\w]+', '', snames, perl=T)
SampleID <- sub('[\\w]+_', '', snames, perl=T)

data_table <- data.frame(sample_name = paste(snames), CellType=CellType, SampleID=SampleID)

group <- interaction(CellType, SampleID)
group


pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/MDS.pdf')

plotMDS(d, col = as.numeric(colnames(d$counts)))
 dev.off()
 
 #3##Voom transformation and calculation of variance weights
 
# mm <- model.matrix(~0 + data_table$SampleID + data_table$CellType)
mm <- model.matrix(~0 +  CellType)

rownames(mm) <- data_table$sample_name

pdf('/home/alsultfh/H1_Project/DATA/Fatima/Plots/group_mean.pdf')
y <- voom(d, mm, plot = T)
dev.off()


#4. Fitting linear models in limma

fit <- lmFit(y, mm)
head(coef(fit))

#Contrasts defined by lineages order 

contrast.matrix <- makeContrasts(
  # HSC Comparisons
  HSC_Vs_LMPP = CellTypeHSC - CellTypeLMPP,
  HSC_Vs_CyclingLMPP = CellTypeHSC - CellTypeCycling,
  HSC_Vs_GMP1 = CellTypeHSC - CellTypeGMP1,
  HSC_Vs_GMP2 = CellTypeHSC - CellTypeGMP2,
  
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
  MEP_Vs_Monocytes = CellTypeMEP - CellTypeMonocytes,
  MEP_Vs_Basophils = CellTypeMEP - CellTypeBasophils,
  MEP_Vs_DC = CellTypeMEP - CellTypeDC,
  MEP_Vs_Erythroid = CellTypeMEP - CellTypeErythroid,

  # GMP1 Comparisons
  GMP1_Vs_GMP2 = CellTypeGMP1 - CellTypeGMP2,
  GMP1_Vs_Monocytes = CellTypeGMP1 - CellTypeMonocytes,
  GMP1_Vs_Basophils = CellTypeGMP1 - CellTypeBasophils,
  GMP1_Vs_DC = CellTypeGMP1 - CellTypeDC,
  GMP1_Vs_Erythroid = CellTypeGMP1 - CellTypeErythroid,
  GMP1_Vs_MEP = CellTypeGMP1 - CellTypeMEP,
  
  # GMP2 Comparisons
  GMP2_Vs_Monocytes = CellTypeGMP2 - CellTypeMonocytes,
  GMP2_Vs_Basophils = CellTypeGMP2 - CellTypeBasophils,
  GMP2_Vs_DC = CellTypeGMP2 - CellTypeDC,
  GMP2_Vs_Erythroid = CellTypeGMP2 - CellTypeErythroid,
  GMP2_Vs_MEP = CellTypeGMP2 - CellTypeMEP,

  # Monocytes Comparisons
  Monocytes_Vs_Basophils = CellTypeMonocytes - CellTypeBasophils,
  Monocytes_Vs_DC = CellTypeMonocytes - CellTypeDC,
  Monocytes_Vs_Erythroid = CellTypeMonocytes - CellTypeErythroid,
  # Basophils Comparisons
  Basophils_Vs_DC = CellTypeBasophils - CellTypeDC,
  Basophils_Vs_Erythroid = CellTypeBasophils - CellTypeErythroid,

  levels = colnames(mm)
  )
 

# step 2
  tmp <- contrasts.fit(fit, contrast.matrix)#linear model lm 
  tmp <- eBayes(tmp)#empirical Bayes smoothing to the standard errors of the estimated coefficients in the linear model fits

# 13819    32

#step 3



## All contrast all DE 
saveRDS(tmp, './Results_DE_updated_Lineage.rds')
results_de <- readRDS('./Results_DE_updated_Lineage.rds')
## include all DEG genes

results_list <- list()  
for (contrast in colnames(contrast.matrix)) {
    print(contrast)
    top_table <- topTable(tmp, coef=contrast, sort.by="P", adjust='fdr', n=Inf)
    filtered_table <- top_table[top_table$adj.P.Val < 0.05 & abs(top_table$logFC) > 1 ,]    
    # Store the result in the list
    results_list[[contrast]] <- filtered_table
}

# Define the list of histone linker symbols
histone_linkers_symbol <- c("H1F0", "HIST1H1A", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H1B", "H1FX")
histone_linkers_symbol %in% rownames(cts)

# Initialize an empty list to store results
results_list <- list() 

# Loop through each contrast in your contrast matrix
for (contrast in colnames(contrast.matrix)) {
    print(contrast)
    # Perform DE analysis and obtain the top table of results
    top_table <- topTable(tmp, coef=contrast, sort.by="P", adjust='fdr', n=Inf)
    
    # Filter results based on adjusted P-value and log fold change
    filtered_table <- top_table[top_table$adj.P.Val < 0.1 & abs(top_table$logFC) > 0.5,]    
    
    # Further filter the table to include only rows corresponding to histone linker genes
    histone_linker_DE <- filtered_table[rownames(filtered_table) %in% histone_linkers_symbol,]
    
    # Store the filtered result in the list
    results_list[[contrast]] = histone_linker_DE 
}
pdf('/home/alsultfh/H1_Project/DATA/Fatima/Heatmap.pdf')
Idents(seurat_young) <- 'CellType'
DoHeatmap(seurat_young, features=histone_linkers_symbol)
DimPlot(seurat_young, label=TRUE) + NoLegend()
FeaturePlot(seurat_young, features=histone_linkers_symbol, order=TRUE)
dev.off()

pdf('/home/alsultfh/H1_Project/DATA/Fatima/UmspipeFel.pdf')
plotter <- FetchData(seurat_young, vars=c(histone_linkers_symbol,'CellType' ))
plotter <- reshape2::melt(plotter)
cowplot::plot_grid(
cowplot::plot_grid(NULL, 
      DimPlot(seurat_young, label=TRUE, pt.size=3) + NoLegend() + NoAxes(),
      NULL, nrow=1,
      rel_widths =c(0.2,0.5,0.2)
), 
FeaturePlot(seurat_young, features=c('H1F0', 'HIST1H1E'), pt.size=1, order=TRUE) & NoLegend() & NoAxes(), 
ggplot(plotter[plotter$CellType == 'PreB', ], aes(x=variable, y =value, fill =variable )) + 
geom_violin(scale='width') +scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) + theme_classic() + theme(legend.position='none'),
nrow=3, rel_heights= c(0.3,0.2,0.15))
dev.off()


# Initialize an empty data frame to store results
all_results <- data.frame()

for (contrast in names(results_list)){
  tmp <- results_list[[contrast]]
  tmp$gene <- rownames(tmp)
  results <- tmp[tmp$gene %in% histone_linkers_symbol, ]
  
  if(nrow(results) > 0) {
    results$contrast <- contrast  # Add a column for the contrast
    all_results <- rbind(all_results, results)  # Append to the all_results data frame
  }
}

# Write the results to a CSV file
#write.csv(all_results, "/home/alsultfh/H1_Project/DATA/Fatima/contrast_results_updated.csv", row.names = FALSE)

#contrast_results <- read.csv("/home/alsultfh/H1_Project/DATA/Fatima/contrast_results_updated.csv")


#### Boxplot for histone after LV and pdeudoBulk



#This plot generated using log trasform after performing pseudo-bulk workflow##
cts_norm <- cpm(cts, log=TRUE)

# List of genes to plot
gene_2_plot <- c("H1F0", "HIST1H1A", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H1B", "H1FX")

# First, check if the genes to plot are in the rownames and filter cts_norm accordingly
filtered_data <- cts_norm[rownames(cts_norm) %in% gene_2_plot, ]
 
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


desired_order <- c("HSC", "LMPP", "Cycling_LMPP","CLP", "PreB", "MEP","Erythroid", "GMP1", "GMP2", "Monocytes", "Basophils", "DC")

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
         y = "LogCPM")
  
  # Save the plot with larger dimensions
  # Define the file name based on the gene, adjusting the path as needed
  file_name <- paste0("/home/alsultfh/H1_Project/DATA/Fatima/Plots/Limma_VOOM_plot_", gene, ".png")
  ggsave(file_name, plot = p, width = 10, height = 8, dpi = 300)
}




##signifcant annotation 
