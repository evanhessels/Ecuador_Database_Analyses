
# First Read in Sequence Data as a DNAbin and setup data---------------------------------------------

library(spider)
library(ape)
library(tidyverse)

Data <- read.dna("Alignment.fasta", format = "fasta")
barcodeDist <- ape::dist.dna(Data, model = "raw", pairwise.deletion = TRUE)

# Clean Labels ------------------------------------------------------------

(barcodeSpp <- dimnames(Data)[[1]] %>% 
   gsub("cf\\.\\_", "", .) %>% 
   gsub("\\_\\(reversed)", "", .) %>% 
   gsub("sp\\.\\_", "sp\\.", .) %>% 
   gsub("\\_[0-9]+$", "", .) %>% 
   gsub("aff\\.\\_", "", .) %>%
   gsub("PCR", "", .) %>%
   gsub("gr\\.\\_", "", .) %>% 
   gsub("\\_\\_", "\\_", .) %>% 
   gsub("Upper\\_Aguarico", "UpperAguarico", .) %>% 
   strsplit(.,  split  = "_") %>% 
   sapply(., \(x) paste(x[1], x[2], sep = "_")))

# Check ThreshOpt to Determine error at 1% as a test --------------------------------

threshOpt(barcodeDist, barcodeSpp, threshold = 0.0001)

# Conduct the analysis over a range of values to determine the optimum threshold --------

threshVal <- seq(0.00,0.015, by = 0.001)
opt <- lapply(threshVal, function(x) threshOpt(barcodeDist, barcodeSpp, thresh = x))
optMat <- do.call(rbind, opt)

# Convert to Percentage ---------------------------------------------------

optMat[1:20, 1:2] <- apply(optMat[1:20, 1:2]*100, 2, function(x) paste0(x, "%"))

# Visualize the results ---------------------------------------------------

graphics::barplot(t(optMat)[4:5,], names.arg=optMat[,1], xlab=("Threshold values"),
                  ylab="Cumulative error", cex.axis=1, cex.names=1.25, cex.lab=1.25)

graphics::legend(x = 3.5, y = 200, legend = c("False positives", "False negatives"),
                 fill = c("grey75", "grey25"), cex= 1)

threshID(barcodeDist, barcodeSpp, threshold = 0.004)


# Convert to ggplot -------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Construct the data frame with threshold values
df <- as.data.frame(optMat)
df$Threshold <- threshVal  # Add the numeric thresholds as a new column

# Create a long-format data frame with cumulative error
df_long <- df %>%
  select(Threshold, FalsePos = 4, FalseNeg = 5) %>%
  pivot_longer(cols = c(FalsePos, FalseNeg),
               names_to = "ErrorType", values_to = "Count") %>%
  mutate(
    ErrorType = recode(ErrorType,
                       "FalsePos" = "False positives",
                       "FalseNeg" = "False negatives")
  )

# Plot #1
p1 <-ggplot(df_long, aes(x = Threshold, y = Count, fill = ErrorType)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0.004, linetype = "dotted", color = "red", linewidth = 0.5) +
  scale_x_continuous(labels = label_percent(accuracy = 0.1)) +
  scale_fill_manual(values = c("grey75", "grey25")) +
  labs(x = "Threshold values", y = "Cumulative error") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.2, 0.70),             # Inside top-left
    legend.justification = c(0, 1),              # Align legend's top-left corner
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA)  # Add semi-transparent background
  )



# Redo Analysis for Second Plot -------------------------------------------

Data2 <- read.dna("Alignment2.fasta", format = "fasta")
barcodeDist2 <- ape::dist.dna(Data2, model = "raw", pairwise.deletion = TRUE)

# Clean Labels ------------------------------------------------------------

(barcodeSpp2 <- dimnames(Data2)[[1]] %>% 
   gsub("cf\\.\\_", "", .) %>% 
   gsub("\\_\\(reversed)", "", .) %>% 
   gsub("sp\\.\\_", "sp\\.", .) %>% 
   gsub("\\_[0-9]+$", "", .) %>% 
   gsub("aff\\.\\_", "", .) %>%
   gsub("PCR", "", .) %>%
   gsub("gr\\.\\_", "", .) %>% 
   gsub("\\_\\_", "\\_", .) %>% 
   gsub("Upper\\_Aguarico", "UpperAguarico", .) %>% 
   strsplit(.,  split  = "_") %>% 
   sapply(., \(x) paste(x[1], x[2], sep = "_")))

# Check ThreshOpt to Determine error at 1% as a test --------------------------------

threshOpt(barcodeDist2, barcodeSpp2, threshold = 0.0001)

# Conduct the analysis over a range of values to determine the optimum threshold --------

threshVal2 <- seq(0.00,0.015, by = 0.001)
opt2 <- lapply(threshVal2, function(x) threshOpt(barcodeDist2, barcodeSpp2, thresh = x))
optMat2 <- do.call(rbind, opt2)

# Convert to Percentage ---------------------------------------------------

optMat2[1:20, 1:2] <- apply(optMat2[1:20, 1:2]*100, 2, function(x) paste0(x, "%"))

# Visualize the results ---------------------------------------------------

graphics::barplot(t(optMat2)[4:5,], names.arg=optMat2[,1], xlab=("Threshold values"),
                  ylab="Cumulative error", cex.axis=1, cex.names=1.25, cex.lab=1.25)

graphics::legend(x = 3.5, y = 200, legend = c("False positives", "False negatives"),
                 fill = c("grey75", "grey25"), cex= 1)

threshID(barcodeDist2, barcodeSpp2, threshold = 0.004)


# Convert to ggplot -------------------------------------------------------

# Construct the data frame with threshold values
df2 <- as.data.frame(optMat2)
df2$Threshold2 <- threshVal2  # Add the numeric thresholds as a new column

# Create a long-format data frame with cumulative error
df_long2 <- df2 %>%
  select(Threshold2, FalsePos = 4, FalseNeg = 5) %>%
  pivot_longer(cols = c(FalsePos, FalseNeg),
               names_to = "ErrorType", values_to = "Count") %>%
  mutate(
    ErrorType = recode(ErrorType,
                       "FalsePos" = "False positives",
                       "FalseNeg" = "False negatives")
  )

# Plot #1
p2 <-ggplot(df_long2, aes(x = Threshold2, y = Count, fill = ErrorType)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0.009, linetype = "dotted", color = "red", linewidth = 0.5) +
  scale_x_continuous(labels = label_percent(accuracy = 0.1)) +
  scale_fill_manual(values = c("grey75", "grey25")) +
  labs(x = "Threshold values", y = "Cumulative error") +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.2, 0.70),             # Inside top-left
    legend.justification = c(0, 1),              # Align legend's top-left corner
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA)  # Add semi-transparent background
  )


# Combine Graphs into one Figure ------------------------------------------

library(patchwork)

FinalGraph <- (p1 + p2) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") &
    theme(legend.position = "right")

ggsave("Full Barcode + Vert General mini.PDF", width = 10, height = 5, plot = last_plot(), dpi = 300 )


