# Emmanuel Makinde
# Install / Import Libraries
install.packages("remotes")
remotes::install_version("grr", version = "0.9.5", repos = "cran.us.r-project.org")
install.packages("pointblank")
BiocManager::install("orthogene")
install.packages("tidyverse")
library(orthogene)
library(tidyverse)
library(pointblank)

# Set wd
setwd('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/CaMK2A_ALL/Ortholog_Conversions')

# Define a function to perform of orthogene conversion and save the output into a csv
ortholog_convert <- function(genes_list, condition_name){
  

  print(paste0("Converting Orthologs for: ", condition_name,"..."))
  
  ortholog_conversion <- convert_orthologs(
    genes_list,
    gene_input = "V1",
    input_species = "mouse",
    output_species = "Homo Sapiens",
    gene_output = "columns",
  )
  file_name <- affix_date(filename = condition_name)
  write.csv(ortholog_conversion, file = file_name)
  "File Saved."
}

# Pull all the genes from the TRAP-seq analysis
TRAP_data <- read.csv("/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/TRAP_Data/Camk2a_ShockvsBox_60min_SIG.csv",header = F)
mouse_gene_names <- TRAP_data["V1"]

TRAP_CaMK2a_Shock_HC_60min <- read.csv('/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/00 TRAP_Data/Camk2a_ShockvsHc_60min.csv',header = F)
TRAP_CaMK2a_Shock_HC_15min <- read.csv('/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/00 TRAP_Data/ShockvsHC_ALL_Camk2a15min.csv',header = F)
TRAP_CaMK2a_Shock_Box_60min <- read.csv('/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/00 TRAP_Data/Camk2a_ShockvsBox_60min.csv',header = F)
TRAP_CaMK2a_Shock_Box_15min <- read.csv('/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/00 TRAP_Data/ShockvsBox_ALL_Camk2a15min.csv',header = F)

# Run the orthogene conversion on each
ortholog_convert(TRAP_CaMK2a_Shock_HC_60min,"60min_CaMK2a_Shock_HC_Orthogene_Conversion_TOTAL.csv")
ortholog_convert(TRAP_CaMK2a_Shock_HC_15min,"15min_CaMK2a_Shock_HC_Orthogene_Conversion_TOTAL.csv")
ortholog_convert(TRAP_CaMK2a_Shock_Box_60min,"60min_CaMK2a_Shock_Box_Orthogene_Conversion_TOTAL.csv")
ortholog_convert(TRAP_CaMK2a_Shock_Box_15min,"15min_CaMK2a_Shock_Box_Orthogene_Conversion_TOTAL.csv")

### Extra ###
# Pull the genes from the CLIP data (hnRNPD Targets) // This is only needed if we are interested in finding the mouse ortholog of human genes
CLIP_data <- read.csv('/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/CLIP_Data/hnRNPD_HITS_PARCLIP_Paper.csv')
human_gene_names <- CLIP_data["GeneName"]


# Generate a report with all the possible orthologs
orth_human <- orthogene::report_orthologs(target_species = "human",
                                          reference_species = "mouse",
                                          method_all_genes = method,
                                          method_convert_orthologs = method)
ortholog_report <- as.data.frame(orth_human)
write.csv(ortholog_report, file = "ortholog_mouse_to_human_report.csv")
