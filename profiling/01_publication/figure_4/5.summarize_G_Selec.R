#!/usr/bin/env Rscript

# ==============================================================================
# Strict UpSet plots for GWAS loci + selection overlaps
# Consistent-size Illustrator-friendly version
#
# This version:
#   1. keeps ALL information (no filtering of intersections)
#   2. keeps panel size CONSISTENT across Architecture / Developmental / Panicle
#   3. restores previous bar border styling
#   4. slightly enlarges overall figure and reduces font/dot sizes for clarity
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

# ==============================================================================
# PATHS
# ==============================================================================

repo_root <- "/mnt/nrdstor/jyanglab/subhash/git/GWAS_sorghum_seedling"

INFILE <- file.path(
  repo_root,
  "largedata",
  "figure_4",
  "panel_A_upset_v5",
  "all_loci_upset_data.tsv"
)

OUTDIR <- file.path(
  repo_root,
  "largedata",
  "figure_4",
  "panel_B_upset_strict_consistent_size"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# FIGURE SIZE
# ==============================================================================

# Keep ALL intersections
MAX_INTERSECTIONS <- Inf

# Consistent size for all individual panels
IND_WIDTH_CM  <- 8.2
IND_HEIGHT_CM <- 6.6

# Composite size for three equally sized panels
COMPOSITE_WIDTH_CM  <- 25.0
COMPOSITE_HEIGHT_CM <- 7.6

# ==============================================================================
# LABELS / TRAIT GROUPS
# ==============================================================================

SEL_BAL <- "Balancing sel."
SEL_POS <- "Positive sel."

canonical_traits <- list(
  Architectural = c(
    "Branch IL",
    "Plant height",
    "Branch number",
    "Tillers"
  ),
  Developmental = c(
    "Days to flower",
    "Extant leaf #",
    "Flag leaf L",
    "Flag leaf W",
    "Leaf angle SD",
    "Median leaf angle",
    "Stem diam. lower",
    "Stem diam. upper",
    "Third leaf L",
    "Third leaf W"
  ),
  Panicle = c(
    "Panicle grain wt",
    "Panicles/plot",
    "Rachis diam. lower",
    "Rachis diam. upper",
    "Rachis length"
  )
)

# Data uses "Architecture"; title can still show "Architectural"
category_data_name <- c(
  Architectural = "Architecture",
  Developmental = "Developmental",
  Panicle = "Panicle"
)

categories_plot <- c("Architectural", "Developmental", "Panicle")
conditions <- c("NR", "HN", "LN")

condition_label <- c(
  NR = "NR",
  HN = "HN",
  LN = "LN"
)

# ==============================================================================
# COLORS
# ==============================================================================

BAR_COLOR     <- "#B7C1C8"
BAR_BORDER    <- "#4D4D4D"

DOT_ACTIVE    <- "#2C2C2C"
DOT_INACTIVE  <- "#D9D9D9"

BAL_COLOR     <- "#E69F00"
POS_COLOR     <- "#CC3311"

SET_BAR_COLOR <- "#6E8898"

# ==============================================================================
# FONT / SIZE SETTINGS
# ==============================================================================

TITLE_SIZE      <- 8.0
AXIS_TITLE_SIZE <- 5.8
AXIS_TEXT_SIZE  <- 4.8
ROW_TEXT_SIZE   <- 5.0
BAR_TEXT_SIZE   <- 2.2
SET_TEXT_SIZE   <- 1.9
DOT_SIZE        <- 1.9
LINE_WIDTH      <- 0.55

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading loci table:\n", INFILE, "\n\n")

all_loci <- fread(INFILE)

cat("Rows loaded:", nrow(all_loci), "\n")
cat("Columns loaded:", ncol(all_loci), "\n\n")

required_cols <- c("locus_id", "chr", "start", "end", "Condition", "category", SEL_BAL, SEL_POS)
missing_required <- setdiff(required_cols, names(all_loci))

if (length(missing_required) > 0) {
  stop(
    "Missing required columns in all_loci_upset_data.tsv: ",
    paste(missing_required, collapse = ", ")
  )
}

all_possible_sets <- unique(c(unlist(canonical_traits), SEL_BAL, SEL_POS))

for (cc in all_possible_sets) {
  if (cc %in% names(all_loci)) {
    if (!is.logical(all_loci[[cc]])) {
      all_loci[, (cc) := as.logical(get(cc))]
    }
    all_loci[is.na(get(cc)), (cc) := FALSE]
  }
}

# ==============================================================================
# FUNCTION: Build exact strict intersection data
# ==============================================================================

make_exact_intersections <- function(loci_dt, catg_plot, cond, traits_in_cat) {
  
  catg_data <- category_data_name[[catg_plot]]
  
  sub <- copy(loci_dt[category == catg_data & Condition == cond])
  
  if (nrow(sub) == 0) return(NULL)
  
  valid_traits <- traits_in_cat[traits_in_cat %in% names(sub)]
  if (length(valid_traits) == 0) return(NULL)
  
  traits_with_hits <- valid_traits[
    sapply(valid_traits, function(x) any(sub[[x]], na.rm = TRUE))
  ]
  
  if (length(traits_with_hits) == 0) return(NULL)
  
  # Selection included in TRUE intersection definition
  upset_sets <- c(traits_with_hits, SEL_BAL, SEL_POS)
  
  sub[, exact_member := apply(.SD, 1, function(x) {
    active_sets <- upset_sets[as.logical(x)]
    paste(active_sets, collapse = " & ")
  }), .SDcols = upset_sets]
  
  sub <- sub[exact_member != ""]
  if (nrow(sub) == 0) return(NULL)
  
  inter <- sub[, .(total_loci = .N), by = exact_member]
  setorder(inter, -total_loci, exact_member)
  
  # Keep ALL intersections
  if (is.finite(MAX_INTERSECTIONS)) {
    inter <- inter[seq_len(min(.N, MAX_INTERSECTIONS))]
    sub <- sub[exact_member %in% as.character(inter$exact_member)]
  }
  
  inter[, has_balancing := grepl(SEL_BAL, exact_member, fixed = TRUE)]
  inter[, has_positive := grepl(SEL_POS, exact_member, fixed = TRUE)]
  
  inter[, selection_class := fifelse(
    has_balancing & has_positive, "Balancing + Positive",
    fifelse(
      has_balancing, "Balancing only",
      fifelse(has_positive, "Positive only", "No selection")
    )
  )]
  
  inter[, n_trait_sets := sapply(as.character(exact_member), function(x) {
    active <- strsplit(x, " & ", fixed = TRUE)[[1]]
    sum(active %in% traits_with_hits)
  })]
  
  inter[, Condition := cond]
  inter[, category := catg_plot]
  
  setcolorder(
    inter,
    c(
      "Condition",
      "category",
      "exact_member",
      "total_loci",
      "n_trait_sets",
      "has_balancing",
      "has_positive",
      "selection_class"
    )
  )
  
  # Dot matrix data
  dot_list <- list()
  
  for (m in inter$exact_member) {
    active_sets <- strsplit(as.character(m), " & ", fixed = TRUE)[[1]]
    
    for (s in upset_sets) {
      dot_list[[length(dot_list) + 1]] <- data.table(
        exact_member = as.character(m),
        set_name = s,
        active = s %in% active_sets,
        row_type = ifelse(s %in% c(SEL_BAL, SEL_POS), "selection", "trait")
      )
    }
  }
  
  dots <- rbindlist(dot_list)
  
  # Set sizes
  set_sizes <- rbindlist(lapply(upset_sets, function(s) {
    data.table(
      set_name = s,
      size = sum(sub[[s]], na.rm = TRUE),
      row_type = ifelse(s %in% c(SEL_BAL, SEL_POS), "selection", "trait")
    )
  }))
  
  row_levels <- rev(upset_sets)
  
  set_sizes[, set_name := factor(set_name, levels = row_levels)]
  dots[, set_name := factor(set_name, levels = row_levels)]
  
  inter[, exact_member := factor(exact_member, levels = as.character(inter$exact_member))]
  dots[, exact_member := factor(exact_member, levels = levels(inter$exact_member))]
  
  # Connection lines
  conn_list <- list()
  
  for (m in levels(dots$exact_member)) {
    active_rows <- dots[exact_member == m & active == TRUE]$set_name
    
    if (length(active_rows) > 1) {
      y_pos <- as.numeric(factor(active_rows, levels = levels(dots$set_name)))
      
      conn_list[[length(conn_list) + 1]] <- data.table(
        exact_member = m,
        ymin = min(y_pos),
        ymax = max(y_pos)
      )
    }
  }
  
  conn <- if (length(conn_list) > 0) rbindlist(conn_list) else NULL
  
  if (!is.null(conn)) {
    conn[, exact_member := factor(exact_member, levels = levels(dots$exact_member))]
  }
  
  list(
    sub = sub,
    inter = inter,
    dots = dots,
    conn = conn,
    set_sizes = set_sizes,
    upset_sets = upset_sets,
    traits_with_hits = traits_with_hits,
    catg_plot = catg_plot,
    cond = cond
  )
}

# ==============================================================================
# FUNCTION: Make one strict UpSet plot
# ==============================================================================

make_strict_upset_plot <- function(loci_dt, catg_plot, cond, traits_in_cat) {
  
  dat <- make_exact_intersections(
    loci_dt = loci_dt,
    catg_plot = catg_plot,
    cond = cond,
    traits_in_cat = traits_in_cat
  )
  
  if (is.null(dat)) return(NULL)
  
  inter <- dat$inter
  dots <- dat$dots
  conn <- dat$conn
  set_sizes <- dat$set_sizes
  
  # Save exact count table
  table_name <- paste0(
    "strict_counts_",
    tolower(catg_plot),
    "_",
    cond,
    ".tsv"
  )
  
  fwrite(inter, file.path(OUTDIR, table_name), sep = "\t")
  
  # ---------------------------------------------------------------------------
  # Top bar plot: exact intersection size
  # Previous border styling retained
  # ---------------------------------------------------------------------------
  
  p_bars <- ggplot(inter, aes(x = exact_member, y = total_loci)) +
    geom_col(
      fill = BAR_COLOR,
      color = BAR_BORDER,
      width = 0.56,
      linewidth = 0.35
    ) +
    geom_text(
      aes(label = total_loci),
      vjust = -0.35,
      size = BAR_TEXT_SIZE,
      fontface = "bold",
      color = "black"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(y = "Exact\nintersection\nsize") +
    theme_classic(base_size = 7) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank(),
      axis.title.y = element_text(face = "bold", size = AXIS_TITLE_SIZE),
      axis.text.y  = element_text(size = AXIS_TEXT_SIZE),
      plot.margin = margin(2, 3, 0, 3)
    )
  
  # ---------------------------------------------------------------------------
  # Dot matrix
  # ---------------------------------------------------------------------------
  
  dots[, dot_class := fifelse(
    active == TRUE & set_name == SEL_BAL, "balancing",
    fifelse(
      active == TRUE & set_name == SEL_POS, "positive",
      fifelse(active == TRUE, "trait", "inactive")
    )
  )]
  
  dot_color_map <- c(
    "trait" = DOT_ACTIVE,
    "balancing" = BAL_COLOR,
    "positive" = POS_COLOR,
    "inactive" = DOT_INACTIVE
  )
  
  p_dots <- ggplot(dots, aes(x = exact_member, y = set_name)) +
    geom_point(
      aes(color = dot_class),
      size = DOT_SIZE,
      show.legend = FALSE
    ) +
    scale_color_manual(values = dot_color_map)
  
  if (!is.null(conn) && nrow(conn) > 0) {
    p_dots <- p_dots +
      geom_segment(
        data = conn,
        aes(
          x = exact_member,
          xend = exact_member,
          y = ymin,
          yend = ymax
        ),
        inherit.aes = FALSE,
        linewidth = LINE_WIDTH,
        color = DOT_ACTIVE
      )
  }
  
  p_dots <- p_dots +
    geom_hline(
      yintercept = 2.5,
      linetype = "dashed",
      color = "#999999",
      linewidth = 0.22
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(
        size = ROW_TEXT_SIZE,
        face = "bold",
        color = sapply(levels(dots$set_name), function(lbl) {
          if (lbl == SEL_BAL) return(BAL_COLOR)
          if (lbl == SEL_POS) return(POS_COLOR)
          return("black")
        })
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#EDEDED", linewidth = 0.22),
      plot.margin = margin(0, 3, 2, 3)
    )
  
  # ---------------------------------------------------------------------------
  # Left set-size bars
  # Previous border styling retained
  # ---------------------------------------------------------------------------
  
  set_sizes[, fill_col := fifelse(
    set_name == SEL_BAL, BAL_COLOR,
    fifelse(set_name == SEL_POS, POS_COLOR, SET_BAR_COLOR)
  )]
  
  p_sets <- ggplot(set_sizes, aes(x = size, y = set_name)) +
    geom_col(
      aes(fill = fill_col),
      width = 0.58,
      color = BAR_BORDER,
      linewidth = 0.35,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    geom_text(
      data = set_sizes[size > 0],
      aes(label = size),
      hjust = 1.12,
      size = SET_TEXT_SIZE,
      fontface = "bold",
      color = "white"
    ) +
    scale_x_reverse(expand = expansion(mult = c(0.22, 0))) +
    labs(x = "Set size", y = NULL) +
    theme_classic(base_size = 7) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.title.x = element_text(face = "bold", size = AXIS_TITLE_SIZE),
      axis.text.x  = element_text(size = AXIS_TEXT_SIZE),
      plot.margin = margin(0, 0, 2, 3)
    )
  
  list(
    p_bars = p_bars,
    p_dots = p_dots,
    p_sets = p_sets,
    title = paste0(catg_plot, " - ", condition_label[[cond]]),
    inter = inter
  )
}

# ==============================================================================
# FUNCTION: Assemble one panel
# ==============================================================================

assemble_upset <- function(parts) {
  
  title_grob <- textGrob(
    parts$title,
    gp = gpar(fontface = "bold", fontsize = TITLE_SIZE)
  )
  
  left_w  <- 0.23
  right_w <- 0.77
  
  # Give slightly more room to bottom matrix area
  top_h <- 0.36
  bot_h <- 0.64
  
  blank_left <- ggplot() + theme_void()
  
  p_bars2 <- parts$p_bars +
    scale_x_discrete(expand = expansion(add = c(0.35, 0.35))) +
    theme(plot.margin = margin(2, 3, 0, 3))
  
  p_dots2 <- parts$p_dots +
    scale_x_discrete(expand = expansion(add = c(0.35, 0.35))) +
    theme(plot.margin = margin(0, 3, 2, 3))
  
  p_sets2 <- parts$p_sets +
    theme(plot.margin = margin(0, 0, 2, 3))
  
  aligned_right <- cowplot::align_plots(
    p_bars2,
    p_dots2,
    align = "v",
    axis = "lr"
  )
  
  p_bars_aligned <- aligned_right[[1]]
  p_dots_aligned <- aligned_right[[2]]
  
  aligned_left <- cowplot::align_plots(
    blank_left,
    p_sets2,
    align = "v",
    axis = "lr"
  )
  
  blank_left_aligned <- aligned_left[[1]]
  p_sets_aligned <- aligned_left[[2]]
  
  top_row <- plot_grid(
    blank_left_aligned,
    p_bars_aligned,
    nrow = 1,
    rel_widths = c(left_w, right_w),
    align = "h",
    axis = "tb"
  )
  
  bottom_row <- plot_grid(
    p_sets_aligned,
    p_dots_aligned,
    nrow = 1,
    rel_widths = c(left_w, right_w),
    align = "h",
    axis = "tb"
  )
  
  combined <- plot_grid(
    top_row,
    bottom_row,
    ncol = 1,
    rel_heights = c(top_h, bot_h),
    align = "v",
    axis = "lr"
  )
  
  plot_grid(
    title_grob,
    combined,
    ncol = 1,
    rel_heights = c(0.07, 0.93)
  )
}

# ==============================================================================
# FUNCTION: Save plot
# ==============================================================================

save_panel <- function(plot_obj, pdf_path, png_path, width_cm, height_cm) {
  
  ggsave(
    filename = pdf_path,
    plot = plot_obj,
    width = width_cm,
    height = height_cm,
    units = "cm",
    device = grDevices::cairo_pdf,
    limitsize = FALSE
  )
  
  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width_cm,
    height = height_cm,
    units = "cm",
    dpi = 600,
    limitsize = FALSE
  )
}

# ==============================================================================
# GENERATE FIGURES
# ==============================================================================

cat("Generating strict UpSet plots...\n\n")

all_count_tables <- list()
table_idx <- 0L

for (cond in conditions) {
  
  panels_this_cond <- list()
  
  for (catg_plot in categories_plot) {
    
    cat("Processing:", cond, catg_plot, "\n")
    
    parts <- make_strict_upset_plot(
      loci_dt = all_loci,
      catg_plot = catg_plot,
      cond = cond,
      traits_in_cat = canonical_traits[[catg_plot]]
    )
    
    if (is.null(parts)) {
      cat("  Skipping:", cond, catg_plot, "\n")
      next
    }
    
    assembled <- assemble_upset(parts)
    panels_this_cond[[length(panels_this_cond) + 1]] <- assembled
    
    fname_base <- paste0(
      "strict_upset_",
      tolower(catg_plot),
      "_",
      cond
    )
    
    save_panel(
      plot_obj = assembled,
      pdf_path = file.path(OUTDIR, paste0(fname_base, ".pdf")),
      png_path = file.path(OUTDIR, paste0(fname_base, ".png")),
      width_cm = IND_WIDTH_CM,
      height_cm = IND_HEIGHT_CM
    )
    
    cat("  Saved individual:", fname_base, "\n")
    
    table_idx <- table_idx + 1L
    all_count_tables[[table_idx]] <- copy(parts$inter)
  }
  
  if (length(panels_this_cond) > 0) {
    
    legend_grob <- ggplot() +
      annotate("point", x = 1.00, y = 1, color = BAL_COLOR, size = 2.2) +
      annotate(
        "text",
        x = 1.12,
        y = 1,
        label = "Balancing selection",
        hjust = 0,
        size = 2.2,
        fontface = "bold"
      ) +
      annotate("point", x = 2.62, y = 1, color = POS_COLOR, size = 2.2) +
      annotate(
        "text",
        x = 2.74,
        y = 1,
        label = "Positive selection",
        hjust = 0,
        size = 2.2,
        fontface = "bold"
      ) +
      annotate("point", x = 4.00, y = 1, color = DOT_ACTIVE, size = 2.2) +
      annotate(
        "text",
        x = 4.12,
        y = 1,
        label = "Trait hit",
        hjust = 0,
        size = 2.2
      ) +
      annotate("point", x = 4.95, y = 1, color = DOT_INACTIVE, size = 2.2) +
      annotate(
        "text",
        x = 5.07,
        y = 1,
        label = "No hit",
        hjust = 0,
        size = 2.2
      ) +
      xlim(0.8, 6.0) +
      ylim(0.8, 1.2) +
      theme_void()
    
    cond_composite <- plot_grid(
      plot_grid(
        plotlist = panels_this_cond,
        ncol = length(panels_this_cond),
        align = "none",
        rel_widths = c(1, 1, 1)   # keep panel size consistent
      ),
      legend_grob,
      ncol = 1,
      rel_heights = c(0.91, 0.09)
    )
    
    cond_label_file <- switch(
      cond,
      "NR" = "nitrogen_response",
      "HN" = "high_nitrogen",
      "LN" = "low_nitrogen"
    )
    
    fname_base <- paste0("strict_upset_composite_", cond_label_file)
    
    save_panel(
      plot_obj = cond_composite,
      pdf_path = file.path(OUTDIR, paste0(fname_base, ".pdf")),
      png_path = file.path(OUTDIR, paste0(fname_base, ".png")),
      width_cm = COMPOSITE_WIDTH_CM,
      height_cm = COMPOSITE_HEIGHT_CM
    )
    
    cat("  Saved composite:", fname_base, "\n\n")
  }
}

# ==============================================================================
# SAVE COMBINED COUNT TABLE
# ==============================================================================

if (length(all_count_tables) > 0) {
  
  all_counts <- rbindlist(all_count_tables, use.names = TRUE, fill = TRUE)
  
  fwrite(
    all_counts,
    file.path(OUTDIR, "strict_upset_all_exact_intersection_counts.tsv"),
    sep = "\t"
  )
  
  cat("Saved combined exact count table:\n")
  cat(file.path(OUTDIR, "strict_upset_all_exact_intersection_counts.tsv"), "\n\n")
}

cat("All outputs saved in:\n", OUTDIR, "\n")
cat("Done.\n")
