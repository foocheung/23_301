library(Seurat)
library(ggplot2)
library(patchwork)

# load object
obj_panGI <- readRDS("object_final_annotated_panGI_withDay.rds")

# sanity check that Day exists
table(obj_panGI$Day)

# 1) Single UMAP, colored by Day
p_day <- DimPlot(
  obj_panGI,
  reduction = "umap",
  group.by  = "Day",
  pt.size   = 0.3,
  raster    = TRUE
) + ggtitle("UMAP colored by Day")

# 2) Faceted UMAP: split by Day (same layout, one panel per Day)
p_split <- DimPlot(
  obj_panGI,
  reduction = "umap",
  group.by  = "Day",
  split.by  = "Day",
  pt.size   = 0.3,
  raster    = TRUE
) + ggtitle("UMAP split by Day")

# show them
p_day
p_split


##########


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load object
obj_panGI <- readRDS("object_final_annotated_panGI_withDay.rds")

# Check the Treatment groups
unique(obj_panGI$Treatment)

# Ensure Treatment is a factor with the right order
obj_panGI$Treatment <- factor(obj_panGI$Treatment, levels = c("Untreated", "Treated", "Spike"))

# --- 1) UMAP for Untreated ---
p_untreated <- DimPlot(
  obj_panGI,
  reduction = "umap",
  cells = WhichCells(obj_panGI, expression = Treatment == "Untreated"),
  pt.size = 1,
  raster = TRUE
) + ggtitle("UMAP: Untreated")

# --- 2) UMAP for Treated ---
p_treated <- DimPlot(
  obj_panGI,
  reduction = "umap",
  cells = WhichCells(obj_panGI, expression = Treatment == "Treated"),
  pt.size = 1,
  raster = TRUE
) + ggtitle("UMAP: Treated")

# --- 3) UMAP for Spike ---
p_spike <- DimPlot(
  obj_panGI,
  reduction = "umap",
  cells = WhichCells(obj_panGI, expression = Treatment == "Spike"),
  pt.size = 1,split.by = "Lane",
  raster = TRUE
) + ggtitle("UMAP: Spike")

# Show all three vertically (or side-by-side)
p_untreated / p_treated / p_spike



###########

# Spike-only
obj_spike <- subset(obj_panGI, subset = Treatment == "Spike")

table(obj_spike$Day)
table(obj_spike$Lane)

p_spike_day_lane <- DimPlot(
  obj_spike,
  reduction = "umap",
  group.by  = "Day",
  split.by  = "Lane",
  pt.size   = 1,       # much bigger
  raster    = FALSE    # force proper ggplot, no fade
) + ggtitle("Spike: UMAP by Day split by Lane")


p_spike_day_lane


##############

library(Seurat)
library(ggplot2)
library(dplyr)

# Spike-only (you already have this)
obj_spike <- subset(obj_panGI, subset = Treatment == "Spike")

# Make a combined Day_Lane factor
obj_spike$Day_Lane <- interaction(obj_spike$Day, obj_spike$Lane, sep = "_", drop = TRUE)

# One UMAP, all Spike cells overlaid, colored by Day_Lane
p_spike_day_lane_overlay <- DimPlot(
  obj_spike,
  reduction = "umap",
  group.by  = "Day_Lane",  # color by combined Day+Lane
  pt.size   = 1.5,
  raster    = FALSE
) + ggtitle("Spike: UMAP colored by Day_Lane")

p_spike_day_lane_overlay


########################

p_spike_day_lane / p_spike_day_lane_overlay

########################


library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Load object
obj_panGI <- readRDS("object_final_annotated_panGI_withDay.rds")

# Ensure Treatment is a factor with the right order
obj_panGI$Treatment <- factor(obj_panGI$Treatment,
                              levels = c("Untreated", "Treated", "Spike"))

## ------------------------------
## 1) UMAPs by Treatment
## ------------------------------

# Untreated
p_untreated <- DimPlot(
  obj_panGI,
  reduction = "umap",
  cells = WhichCells(obj_panGI, expression = Treatment == "Untreated"),
  pt.size = 1,
  raster = TRUE
) + ggtitle("UMAP: Untreated")

# Treated
p_treated <- DimPlot(
  obj_panGI,
  reduction = "umap",
  cells = WhichCells(obj_panGI, expression = Treatment == "Treated"),
  pt.size = 1,
  raster = TRUE
) + ggtitle("UMAP: Treated")

# Spike (simple)
p_spike <- DimPlot(
  obj_panGI,
  reduction = "umap",
  cells = WhichCells(obj_panGI, expression = Treatment == "Spike"),
  pt.size = 1,
  raster = TRUE
) + ggtitle("UMAP: Spike")

## ------------------------------
## 2) Spike-only QC plots
## ------------------------------

# Spike-only object
obj_spike <- subset(obj_panGI, subset = Treatment == "Spike")

# (a) Split by Lane, color = Day
p_spike_day_lane <- DimPlot(
  obj_spike,
  reduction = "umap",
  group.by  = "Day",
  split.by  = "Lane",
  pt.size   = 1,
  raster    = FALSE
) + ggtitle("Spike: UMAP by Day split by Lane")

# (b) Overlay, color = Day_Lane
obj_spike$Day_Lane <- interaction(obj_spike$Day, obj_spike$Lane,
                                  sep = "_", drop = TRUE)

p_spike_day_lane_overlay <- DimPlot(
  obj_spike,
  reduction = "umap",
  group.by  = "Day_Lane",
  pt.size   = 1.5,
  raster    = FALSE
) + ggtitle("Spike: UMAP colored by Day_Lane (overlay)")

## ------------------------------
## 3) Join ALL plots
## ------------------------------

# Top row: Untreated / Treated / Spike
top_row <- p_untreated | p_treated | p_spike

# Bottom row: Spike QC (split vs overlay)
bottom_row <- p_spike_day_lane | p_spike_day_lane_overlay

# Full layout
full_plot <- (top_row) / (bottom_row) +
  plot_annotation(title = "UMAP overview by Treatment and Spike batch QC")

full_plot
