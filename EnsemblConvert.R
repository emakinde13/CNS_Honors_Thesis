library(readr)
library(dplyr)
library(biomaRt)
library(pointblank)


convert_to_ensembl<- function(gene_symbols, output_file_name){
  # Build a clean vector of symbols for querying
  symbols <- gene_symbols %>%
    pull(ortholog_gene) %>%
    as.character() %>%
    trimws() %>%
    na.omit() %>%
    unique()
  
  # Conversion: Gene Symbol to Ensembl ID
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  symbol_to_ensembl <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters    = "hgnc_symbol",
    values     = symbols,
    mart       = ensembl
  )
  # Keep it 1:1 mapping
  symbol_to_ensembl_1to1 <- symbol_to_ensembl %>%
    distinct(hgnc_symbol, .keep_all = TRUE)
  
  # Join back onto your original table (keeps 10811 rows)
  gene_symbols_mapped <- gene_symbols %>%
    left_join(symbol_to_ensembl_1to1, by = c("ortholog_gene" = "hgnc_symbol"))
  
  final_ensembl_map <- as.data.frame(gene_symbols_mapped)
  
  names(final_ensembl_map)[names(final_ensembl_map) == 'ensembl_gene_id'] <- 'Ensembl_ID'
  names(final_ensembl_map)[names(final_ensembl_map) == 'input_gene'] <- 'GeneSymbol'
  
  file_name <- affix_date(filename = output_file_name)
  write.csv(final_ensembl_map, file = file_name)
}
CaMK2a_Shock_HC_60min_Total <- read_csv(
  "/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/hnRNPD_Project/60min_CaMK2a_Shock_HC_Orthogene_Conversion_TOTAL_2025-12-25.csv"
)
CaMK2a_Shock_BC_60min_Total <- read_csv(
  "/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/hnRNPD_Project/60min_CaMK2a_Shock_Box_Orthogene_Conversion_TOTAL_2025-12-25.csv"
)
CaMK2a_Shock_HC_15min_Total <- read_csv(
  "/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/hnRNPD_Project/15min_CaMK2a_Shock_HC_Orthogene_Conversion_TOTAL_2025-12-25.csv"
)
CaMK2a_Shock_BC_15min_Total <- read_csv(
  "/Users/emmanuelmakinde/Documents/Data_Analysis/hnRNPD_Project/hnRNPD_Project/15min_CaMK2a_Shock_Box_Orthogene_Conversion_TOTAL_2025-12-25.csv"
)
convert_to_ensembl(CaMK2a_Shock_HC_60min_Total, "Ensembl_CaMK2a_Shock_HC_60min_Total.csv")
convert_to_ensembl(CaMK2a_Shock_BC_60min_Total, "Ensembl_CaMK2a_Shock_BC_60min_Total.csv")
convert_to_ensembl(CaMK2a_Shock_HC_15min_Total, "Ensembl_CaMK2a_Shock_HC_15min_Total.csv")
convert_to_ensembl(CaMK2a_Shock_BC_15min_Total, "Ensembl_CaMK2a_Shock_BC_15min_Total.csv")
