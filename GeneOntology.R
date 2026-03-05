# Install / Load Libraries
install.packages("tibble")
install.packages("tibble")
library(tibble)
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(AnnotationDbi)
library(BiocManager)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pointblank)
library(enrichplot)
library(DOSE)
library(ggplot2)

# Set working directory
setwd('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/GO_Results/FC_0')


# Define a function to run the GO analysis and save output
runGO <- function(sigs,condition_name){
  file_name <- affix_date(filename = condition_name)
  
  rownames(sigs) <- sigs[,1]
  enriched <- rownames(sigs[sigs$log2FoldChange > 0,])
  GO_results <- enrichGO(gene = enriched, "org.Hs.eg.db",keyType = "SYMBOL",ont="BP")
  results <- as.data.frame(GO_results)
  write.csv(results, file = file_name)
  
  return (results)
}

# Define a GO DOWN function
run_GO_down <- function(sigs,condition_name){
  file_name <- affix_date(filename = condition_name)
  
  rownames(sigs) <- sigs[,1]
  enriched <- rownames(sigs[sigs$log2FoldChange < 0,])
  GO_results <- enrichGO(gene = enriched, "org.Hs.eg.db",keyType = "SYMBOL",ont="BP")
  results <- as.data.frame(GO_results)
  write.csv(results, file = file_name)
  
  return (results)
}

# Define a plot function
plot_GO <- function(results,category_count){
  x <- plot(dotplot(results, showCategory = category_count))
  x}

# Load all of the data.
trap_C_60_SB = read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/CaMK2a/C_60_SB.csv')
trap_C_60_SH = read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/CaMK2a/C_60_SH.csv')
trap_C_15_SB = read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/CaMK2a/C_15_SB.csv')
trap_C_15_SH = read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/CaMK2a/C_15_SH.csv')

# Run function on all the data
results_a = runGO(trap_C_60_SB,"Gene_Ontology_60_Shock_Box_CaMK2a_FC_0")
results_b = runGO(trap_C_60_SH,"Gene_Ontology_60_Shock_HC_CaMK2a_FC_0")
results_c = runGO(trap_C_15_SB,"Gene_Ontology_15_Shock_Box_CaMK2a_FC_0")
results_d = runGO(trap_C_15_SH,"Gene_Ontology_15_Shock_HC_CaMK2a_FC_0")

# Run Function on all the down Data
results_a_down = run_GO_down(trap_C_60_SB, "DOWN_Gene_Ontology_60_Shock_Box_CaMK2a_FC_0")
results_b_down = run_GO_down(trap_C_60_SH, "DOWN_Gene_Ontology_60_Shock_HC_CaMK2a_FC_0")
results_c_down = run_GO_down(trap_C_15_SB, "DOWN_Gene_Ontology_15_Shock_Box_CaMK2a_FC_0")
results_d_down = run_GO_down(trap_C_15_SH,"DOWN_Gene_Ontology_15_Shock_HC_CaMK2a_FC_0")

# Plot the Results (Some preliminary plotting here, but I have another script for cleaner plotting)
plot_GO(results_a, 20)
plot_GO(results_b, 20)
plot_GO(results_c, 20)
plot_GO(results_d, 20)
