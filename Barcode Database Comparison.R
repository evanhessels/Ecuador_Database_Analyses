library(tidyverse)
library(scales)
library(ggrepel)
library(ggplot2)
set.seed(1)
(data <- read.delim("VS test V2.txt") %>% 
    mutate(Mine=as.numeric(gsub(pattern = "%", replacement = "", Mine)),
           NCBI = as.numeric(gsub(pattern = "%", replacement = "", NCBI))) %>%
    pivot_longer(cols = c(Mine, NCBI), names_to = "Dataset",values_to = "ID") %>% 
    mutate(Taxon=factor(Taxon, levels=c("Species", "Genus", "Family", "Order")),
           IDlab = as.character(paste0(ID, "%")),
           Dataset = factor(Dataset, labels = c("Custom database", "GenBank"))))

# Plot in BW --------------------------------------------------------------------

for (i in 1:25) {

p <- ggplot(data = data, aes(x = Taxon, y = ID / 100, 
                        group = Dataset, 
                        linetype = Dataset, 
                        shape = Dataset, 
                        label = IDlab)) +
  geom_point(size = 2, color = "black") +
  geom_line(linewidth = 1, color = "black") +
  geom_text_repel(color = "black", size = 5, box.padding = 1.25, max.overlaps = 2, force = 50, min.segment.length = 1) +
  theme_bw() +
  labs(x = NULL, y = "Percentage of Haplotypes Identified") +
  theme(
    axis.text = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    legend.position = c(0.8, 0.4)  # replaces inside + inside coords
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "twodash")) +
  scale_shape_manual(values = c(16, 17, 15, 3, 4))  # customize as needed

  ggsave(
    filename = paste0("Database_comparison_run_", sprintf("%02d", i), ".pdf"),
    plot = p,
    dpi = 300,
    width = 10,
    height = 7,
    bg = "white"
  )
}
