library(ggplot2)
library(patchwork)
library(magick)
library(stringr)
library(grid)

# Set Image folder = project directory ----------------------------------------

img_dir <- "."   

files <- list.files(img_dir, pattern = "jpg|jpeg|png", full.names = TRUE)

# Standardize image size --------------------------------------------------

target_height <- 800

standardize_image <- function(path) {
  img <- image_read(path)
  img <- image_scale(img, paste0("x", target_height))  
  return(img)
}

std_images <- lapply(files, standardize_image)

# Create panels with black background + white labels ----------------------

labels <- tools::file_path_sans_ext(basename(files))
labels <- str_replace_all(labels, "_", " ")

make_panel <- function(img, label) {
  ggplot() +
    annotation_custom(
      rasterGrob(img),
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = Inf
    ) +
    labs(title = label) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      plot.title = element_text(
        hjust = 0.5, size = 14, face = "italic", color = "white"
      ),
      plot.margin = margin(4, 4, 4, 4)
    )
}

panels <- mapply(make_panel, std_images, labels, SIMPLIFY = FALSE)


# Arrange into a grid -----------------------------------------------------

plate <- wrap_plots(panels, ncol = 4) &
  theme(plot.background = element_rect(fill = "black", color = NA))


# Export ------------------------------------------------------------------


ggsave(
  "Ecuador Fish PLate December 9, 25.tiff",
  plate,
  width = 25, height = 15, units = "cm",
  dpi = 600,
  device = "tiff",
  bg = "black"   
)

