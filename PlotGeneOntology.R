library(DOSE)
library(ggplot2)
library(stringr)
library(gridExtra)
library(tidyverse)

# Data
data_translation = data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/GO_Results/Overlaps/Translation_GO_Overlap.csv'))
data_synapses = data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/GO_Results/Overlaps/Synapses_GO_Overlap.csv'))
data_validation <- data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/GO_Results/Overlaps/Validation_GO_Overlap.csv'))
data_synapses_FC_0 <- data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/00 Overlap_Data/GO_Results/Overlaps/Synapses_GO_Overlap_FC_0.csv'))

data_validation$Condition <- factor(data_validation$Condition, levels = c("60’ HC", "60' Box", "15’ HC", "15' Box"))
data_validation$Description = str_wrap(data_validation$Description, width = 30)

# Define a function to create a ggplot object
plot_GO_func <- function(data, title="greetings"){
  data$Condition <- factor(data$Condition, levels = c("60' HC", "60’ Box", "15’ HC", "15' Box")) # Orders the x-axis
  data$Description = str_wrap(data$Description, width = 25) # Controls the width of the GO terms
  
  g <- ggplot(data, mapping = aes(x = Condition, y = Description)) + 
    geom_point(mapping = aes(size = Count, color = p.adjust)) +
    theme_bw(base_size = 15) +
    scale_colour_gradient(limits=c(0.00000000, 0.07), low="red") +
    ylab(NULL) + 
    ggtitle(title) +
    theme(axis.title = element_text(size = 12),  # all titles 
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1,
                                     size = 10, color = "black"),
          axis.text.y = element_text(size = 10),
          panel.border = element_rect(color = "black",
                                      size = .5))
  return (g)
  
}

# Plot and assign to object
gA <- plot_GO_func(data_validation, "RNA Processing")
gB <- plot_GO_func(data_subloc, "Subcellular Localization")
gC <- plot_GO_func(data_translation, title="Protein Synthesis")
gD <- plot_GO_func(data_synapses, title="Synaptic Plasticity")
gE <- plot_GO_func(data_synapses_FC_0, "Synaptic Plasticity: |log2FC| > 0")

# Make the sizes the same
gA$widths <- gB$widths
gA$widths <- gC$widths
gA$widths <- gD$widths

gA$heights <- gB$heights
gA$heights <- gC$heights
gA$heights <- gD$heights

# Arrange on one grid
grid.arrange(g1, gB, gC, gD, ncol=2)
