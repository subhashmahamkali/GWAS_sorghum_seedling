library(knitr)
library(patchwork)
library(ggplot2)
library(grid)

pca_rmd <- "/Users/subhashmahamkali/Documents/gwas_sap/profiling/1.phenotype_analysis/4.pca_geo.Rmd"
geo_rmd <- "/Users/subhashmahamkali/Documents/gwas_sap/profiling/01_publication/figure_1/generate.Rmd"
combined_out_pdf <- "/Users/subhashmahamkali/Documents/gwas_sap/graphs/01_publication/figure_1/figure1_geo_pca_two_panel_A4_landscape.pdf"

tmp_dir <- tempdir()
pca_r <- file.path(tmp_dir, "build_geo_pca_pca.R")
geo_r <- file.path(tmp_dir, "build_geo_pca_geo.R")

knitr::purl(pca_rmd, output = pca_r, quiet = TRUE)
knitr::purl(geo_rmd, output = geo_r, quiet = TRUE)

orig_ggsave <- ggplot2::ggsave
ggsave <- function(...) invisible(NULL)

pca_env <- new.env(parent = globalenv())
geo_env <- new.env(parent = globalenv())

source(pca_r, local = pca_env)
source(geo_r, local = geo_env)

rm(ggsave)

pca_plot <- pca_env$pca_plot
p <- geo_env$p

# Use one shared palette across both panels and match point styling.
shared_colors <- c(
  "durra" = "#FF7F00",
  "kafir" = "#E69F00",
  "caudatum" = "#984EA3",
  "guinea" = "#56B4E9",
  "bicolor" = "#A6A6A6",
  "milo/durra-bicolor" = "#F781BF",
  "mixed/bicolor" = "#A6A6A6",
  "Wild sorghum" = "#4DAF4A"
)
shared_shapes <- c(
  "durra" = 16,
  "kafir" = 16,
  "caudatum" = 16,
  "guinea" = 16,
  "milo/durra-bicolor" = 16,
  "mixed/bicolor" = 16,
  "Wild sorghum" = 17
)

# Match map point styling to the PCA plot while leaving the origin star intact.
p$layers[[5]]$aes_params$size <- 2
p$layers[[5]]$aes_params$alpha <- 0.6
pca_plot$layers[[1]]$aes_params$size <- 2
pca_plot$layers[[1]]$aes_params$alpha <- 0.6

geo_panel <- p +
  scale_color_manual(
    values = shared_colors,
    breaks = c("durra", "kafir", "caudatum", "guinea", "bicolor"),
    name = "Race"
  ) +
  labs(title = NULL, subtitle = NULL, tag = "A") +
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98),
    legend.position = "none",
    plot.margin = margin(5, 3, 5, 5)
  )

pca_panel <- pca_plot +
  scale_color_manual(values = shared_colors, name = NULL) +
  scale_shape_manual(values = shared_shapes, name = NULL) +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.02, 0.98),
    panel.grid.major = element_line(color = "#D0D0D0", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.35),
    axis.ticks.length = unit(0.14, "cm"),
    aspect.ratio = 0.78,
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.42, "cm"),
    legend.box.margin = margin(10, 0, 0, 0),
    plot.margin = margin(5, 5, 12, 3)
  )

two_panel <- (geo_panel | pca_panel) + plot_layout(widths = c(1.6, 1))

orig_ggsave(
  combined_out_pdf,
  two_panel,
  width = 29.7,
  height = 21,
  units = "cm",
  bg = "white"
)

cat("Saved:", combined_out_pdf, "\n")
