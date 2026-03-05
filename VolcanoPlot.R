library(tidyverse)
library(RColorBrewer) # for a colourful plot
library(ggrepel) 
library(tibble)
library(ggplot2)

# Load the data
data_60_box <- data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/Volcano_Plots/C_60_SB_All.csv'))
data_60_hc <- data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/Volcano_Plots/C_60_SH_All.csv'))
data_15_box <- data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/Volcano_Plots/C_15_SB_All.csv'))
data_15_hc <- data.frame(read.csv('/Users/emmanuelmakinde/Documents/00_Project_hnRNPD/001_Gene_Ontology/Volcano_Plots/C_15_SH_All.csv'))

# Volcano Plotting Function
volcanoPlotter <- function(data, title="Title"){
data$diffexpressed <- 'NA' # Creating fake labels so ggplot isn't mad at us
data$diffexpressed[data$log2FoldChange > 0 & data$padj < 0.01] <- "UP"
data$diffexpressed[data$log2FoldChange < 0 & data$padj < 0.01] <- "DOWN"

data$delabel <- NA
data$delabel[data$diffexpressed != "NO"] <- data$GeneName[data$diffexpressed != "NO"]

thresh = head(arrange(data, pvalue, 12)$pvalue[12])
data$label[data$pvalue <= thresh] <- (data$GeneName[data$pvalue <=thresh])


plot <- ggplot(data=data, mapping = aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal()+
  scale_color_manual(values = c('blue','grey','red'))+
  theme(text=element_text(size=20), 
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(),
        axis.line.x = element_line(color='black'),
        axis.line.y = element_line(color = "black")) +
  geom_text_repel(mapping=aes(label=label), size = 4) +
  labs(title = title)
  
return (plot)
}

vplot_60_box <- volcanoPlotter(data_60_box, title="60' Min Shock/Box")
vplot_60_hc <- volcanoPlotter(data_60_hc, title="60' Min HC/Box")
vplot_15_box <- volcanoPlotter(data_15_box, title="15' Min Shock/Box")
vplot_15_hc <- volcanoPlotter(data_15_hc, title="15' Min HC/Box")

vplot_60_box$widths <- vplot_60_hc$widths
vplot_60_box$widths <- vplot_15_box$widths
vplot_60_box$widths <- vplot_15_hc$widths

vplot_60_box$heights <- vplot_60_hc$heights
vplot_60_box$heights <- vplot_15_box$heights
vplot_60_box$heights <- vplot_15_hc$heights

grid.arrange(vplot_60_box, vplot_60_hc, vplot_15_box, vplot_15_hc, ncol=2)
