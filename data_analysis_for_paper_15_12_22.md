---
title: "Impact of ROCKi Treatment on Keratinocytes is Reversible: Single-Cell RNAseq Data Analysis"
author: "George T. Hall"
date: "Compiled on 15 December 2022"
output:
    html_document:
        toc: true
        toc_depth: 2
        toc_float:
            collapsed: false
        number_sections: true
        code_folding: none
        keep_md: yes
---

Copyright University College London 2020-2022. This software is licenced under
the terms of the GNU GENERAL PUBLIC LICENSE Version 3. See COPYING.txt for the
licence details.


```r
knitr::opts_chunk$set(results = "hide", message = FALSE, warning = FALSE,
                      fig.align = "center")
# Tidyverse
library(dplyr)
library(tidyr)

# Figure generation
library(ggplot2)
library(ggrepel)
library(grid)
library(png)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggnewscale)

# Analysis tools
library(Seurat)
library(sctransform)
library(harmony)
library(kBET)
library(SCINA)
library(slingshot)
library(SingleCellExperiment)
library(DAseq)

# Set ggplot theme
theme_update(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
theme_update(axis.text = element_text(size = 10, face = "plain"))
theme_update(axis.title = element_text(size = 10, face = "plain"))

axis_text_theme <- theme(axis.title = element_text(face = "bold",
                                                   color = "black"),
                         axis.text = element_text(face = "bold",
                                                  color = "black"))

blue_color <- "#6E86FF"
red_color <- "#FF8585"

save_images <- TRUE
save_tiff <- function(img, img_name, units = "cm", width = 10, height = 10,
                      res = 300) {
    img_path <- paste0("paper_figs/", img_name)
    tiff(img_path, units = units, width = width, height = height, res = res)
    print(img)
    dev.off()
}
```

# Data pre-processing

## Load data

After sequencing, we generated the gene expression matrix using cellranger. We
used the standard cellranger pipeline consisting of the mkfastq, count, and
aggr steps. We analysed the gene expression data in Seurat using the standard
single-cell RNAseq analysis workflow. We only include genes that appear in at
least 3 cells, and only include cells that express at least 200 genes.


```r
set.seed(10403) # specify default R random seed (for reproducibility)

skin_data <- Read10X(data.dir = "3donor_filtered_feature_bc_matrix")
cells <- list()
cells$all <- CreateSeuratObject(counts = skin_data,
                                project = "3DonorEffectOfRocki",
                                min.cells = 3, min.features = 200,
                                names.delim = "-", names.field = 2)
cells$all <- RenameIdents(object = cells$all, "1" = "1Control-6D ",
                          "2" = "1ROCKi-6D ", "3" = "1Control-6D 6D ",
                          "4" = "1ROCKi-6D 6D ", "5" = "2Control-6D ",
                          "6" = "2ROCKi-6D ", "7" = "2Control-6D 6D ",
                          "8" = "2ROCKi-6D 6D ", "9" = "3Control-6D ",
                          "10" = "3ROCKi-6D ", "11" = "3Control-6D 6D ",
                          "12" = "3ROCKi-6D 6D ")
```

## Filter low quality cells

We removed low quality cells and barcodes corresponding to multiple cells or
empty droplets by excluding cells associated with either an anomalous number of
genes, number of molecules, or percentage of mitochondrial DNA. We used
different thresholds of these measures for each sample to account for
inter-sample variation.


```r
# Compute percentage of mitochondrial DNA in each cell
cells$all[["percent.mt"]] <- PercentageFeatureSet(cells$all, pattern = "^MT-")

vlnplot_limits <- function(data, metric, ylab, limits) {
    df <- data.frame(x = Idents(data), y = data[[metric]][, 1])
    plot <- ggplot() + geom_violin(data = df, aes(x = x, y = y, fill = x)) +
            geom_errorbar(data = limits, aes(x = x, ymin = lower, ymax = upper),
                          color = "red") +
            labs(y = ylab) +
            scale_y_continuous(labels = scales::comma) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none", plot.title = element_blank(),
                  axis.title.x = element_blank()) + axis_text_theme

    return(plot)
}

ncount_df <- data.frame(x = levels(Idents(cells$all)),
                        lower = c(8000, 15000, 8000, 8000,
                                  8000, 8000, 5000, 8000,
                                  8000, 8000, 5000, 5000),
                        upper = c(50000, 62500, 50000, 50000,
                                  50000, 50000, 50000, 50000,
                                  50000, 50000, 40000, 40000))

(ncount_vln <- vlnplot_limits(cells$all, "nCount_RNA",
                              "Number of RNA\nmolecules per cell",
                              ncount_df) + ylim(c(0, 100000)))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(ncount_vln, "fig_s4a.tiff")
}

nfeature_df <- data.frame(x = levels(Idents(cells$all)),
                     lower = c(2500, 3000, 2500, 2500,
                               2500, 2500, 2000, 2500,
                               2500, 2500, 1500, 1500),
                     upper = c(7500, 8500, 7500, 7500,
                               7500, 7500, 7500, 7500,
                               7500, 7500, 6000, 6000))

(nfeature_vln <- vlnplot_limits(cells$all, "nFeature_RNA",
                                "Number of genes\ndetected per cell",
                                nfeature_df))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-3-2.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(nfeature_vln, "fig_s4b.tiff")
}

pct_mito_df <- data.frame(x = levels(Idents(cells$all)),
                          lower = c(rep(2, 8), rep(3, 4)),
                          upper = c(rep(10, 8), rep(12, 4)))

(pct_mito_vln <- vlnplot_limits(cells$all, "percent.mt",
                                paste0("Percentage of mitochrondrial\n",
                                       "DNA detected per cell"),
                                pct_mito_df) + ylim(c(0, 20)))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-3-3.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(pct_mito_vln, "fig_s4c.tiff")
}

pass_qc <- c()
for (j in seq_len(length(cells$all$orig.ident))) {
    i <- as.integer(as.character(cells$all$orig.ident[j]))
    pass_qc <- c(pass_qc, all(cells$all$nFeature_RNA[j] >= nfeature_df$lower[i],
                              cells$all$nFeature_RNA[j] <= nfeature_df$upper[i],
                              cells$all$nCount_RNA[j] >= ncount_df$lower[i],
                              cells$all$nCount_RNA[j] <= ncount_df$upper[i],
                              cells$all$percent.mt[j] >= pct_mito_df$lower[i],
                              cells$all$percent.mt[j] <= pct_mito_df$upper[i]))
}
cells$all <- cells$all[, pass_qc]
```

## Filter outliers

We removed non-keratinocytes by excluding cells expressing at least one count
of markers of melanocytes (MLANA, PMEL, MITF), mesenchymal cells (MTRNR2L6,
MTRNR2L10, MTRNR2L7, MTRNR2L1, PRKAR2B, NR2F1), or fibroblasts (ACTG2, DLK1).


```r
cells$all_with_outliers <- cells$all

# Remove melanocytes
cells$all <- subset(cells$all, subset = MLANA > 1 | PMEL > 1 | MITF > 1,
                    invert = TRUE)
# Remove mesenchymal cells
cells$all <- subset(cells$all,
                    subset = MTRNR2L6 > 1 | MTRNR2L10 > 1 | MTRNR2L7 > 1 |
                             MTRNR2L1 > 1 | PRKAR2B > 1 | NR2F1 > 1,
                    invert = TRUE)
# Remove fibroblasts
cells$all <- subset(cells$all, subset = ACTG2 > 1 | DLK1 > 1, invert = TRUE)
```

## Normalise the data

We normalised the gene expressions using sctransform. 


```r
cells$all <- SCTransform(cells$all)
cells$all_outliers_removed <- cells$all
```

## Regress out cell cycle

We then reduced the effect of cell cycle by assigning scores to each cell
representing the probabilities of it being in different stages of the cell
cycle (using Seurat’s CellCycleScoring function) and then regressing out the
difference between these scores using sctransform. The authors of Seurat
recommend this approach when dealing with cells undergoing differentiation
(https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html -
"Alternative Workflow").


```r
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
cells$all <- CellCycleScoring(cells$all, s.features = s_genes,
                              g2m.features = g2m_genes)
cells$all$CC.difference <- cells$all$S.Score - cells$all$G2M.Score
cells$all <- SCTransform(cells$all, vars.to.regress = c("CC.difference"))
cells$all_cell_cycle_regressed <- cells$all
```

## Run PCA and remove donor effect

We carried out dimensionality reduction using principal components analysis
(PCA), using the 3000 most highly variable genes as features. We used the 50
most significant components since there was no clear point at which the
generated components became less significant. We reduced the inter-donor
variation with Harmony.


```r
# Run PCA
cells$all <- RunPCA(cells$all,
                    features = VariableFeatures(cells$all, assay = "SCT"),
                    npcs = 50)
(elbow_plot <- ElbowPlot(cells$all, ndims = 50) + axis_text_theme)
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(elbow_plot, "fig_s5.tiff")
}

# Remove inter-donor differences using Harmony
donor_idents <- sapply(cells$all$orig.ident,
                       function(x) {
                           if (x %in% c("1", "2", "3", "4")) {
                               "1"
                           } else if (x %in% c("5", "6", "7", "8")) {
                               "2"
                           } else {
                               "3"
                           }
                       })
cells$all$donor_idents <- donor_idents

cells$all <- cells$all %>% RunHarmony("donor_idents", assay.use = "SCT")
cells$all_harmonised <- cells$all

sample_class <- sapply(cells$all$orig.ident,
                       function(x) {
                           if (x %in% c("1", "5", "9")) {
                               "Control-6D"
                           } else if (x %in% c("2", "6", "10")) {
                               "ROCKi-6D"
                           } else if (x %in% c("3", "7", "11")) {
                               "Control-12D"
                           } else {
                               "ROCKi-12D"
                           }
                       })
cells$all$sample_class <- sample_class
cells$all$sample_class <- factor(cells$all$sample_class,
                                           levels = c("Control-6D", "ROCKi-6D",
                                                      "Control-12D",
                                                      "ROCKi-12D"))

treatment_class <- sapply(cells$all$sample_class,
                       function(x) {
                           if (x %in% c("Control-6D", "Control-12D")) {
                               "Control"
                           } else {
                               "ROCKi"
                           }
                       })
cells$all$treatment_class <- treatment_class
```

## Run UMAP

We visualised the cells using UMAP.


```r
cells$d6 <- subset(cells$all,
                   subset = sample_class %in% c("Control-6D", "ROCKi-6D"))
cells$d12 <- subset(cells$all,
                    subset = sample_class %in% c("Control-12D", "ROCKi-12D"))

for (s in c("d6", "d12")) {
    cells[[s]] <- RunUMAP(cells[[s]], dims = 1:50, reduction = "harmony")
}

umap_theme <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                    axis.text.y = element_blank(), axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    plot.title = element_blank())
chart_margin_theme <- theme(plot.margin = margin(0, 1.2, 0, 1.2, "cm"))

labels <- c("Control", "ROCKi-treated")
legend_plot <- ggplot(data.frame(class = as.character(seq_len(length(labels))),
                                 x = seq_len(length(labels))),
                      aes(x = x, y = x, col = class)) +
               geom_point() +
               scale_color_manual(labels = labels,
                                  values = c(blue_color, red_color)) +
               theme(legend.title = element_blank(),
                     legend.key = element_rect(fill = NA, color = NA),
                     legend.margin = margin(c(0, 0, 0, 0)),
                     legend.text = element_text(face = "bold"),
                     legend.background = element_rect(fill = "#00000000",
                                                      size = 0.5,
                                                      linetype = "solid")) +
               guides(colour = guide_legend(override.aes = list(size = 2)))
(umap_legend <- as_ggplot(get_legend(legend_plot, "left")))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(umap_legend, "umap_legend.tiff")
}

plot_treatment_class_umap <- function(cells) {
    p <- DimPlot(cells, group.by = "treatment_class",
                 cols = c(blue_color, red_color),
                 pt.size = -0.1) +
         coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

for (i in list(c("d6", "fig_5a.tiff"), c("d12", "fig_5b.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    umap <- plot_treatment_class_umap(cells[[s]])
    print(umap)
    if (save_images) {
        save_tiff(umap, out_name)
    }
}
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-8-2.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-8-3.png" style="display: block; margin: auto;" />

Visualising the 6D and 12D cell populations using UMAP provides an intuitive
indication that the difference between treated and control cells is greater at
6D than it is at 12D. There is substantially less overlap between the treated
and untreated cells at 6D than at 12D. We now measure the amount of overlap
using DAseq.

## Differential abundance testing using DAseq


```r
create_daseq_umap <- function(daseq_out) {
    daseq_umap <- daseq_out$pred.plot + umap_theme +
                  theme(legend.position = "none")

    return(daseq_umap)
}

for (i in list(c("d6", 6, "Control-6D", "ROCKi-6D", "fig_5c.tiff"),
               c("d12", 12, "Control-12D", "ROCKi-12D", "fig_5d.tiff"))) {
    s <- i[1]
    num_days <- i[2]
    l1 <- i[3]
    l2 <- i[4]
    out_name <- i[5]

    res <- getDAcells(Embeddings(object = cells[[s]], reduction = "harmony"),
                      cell.labels = as.character(cells[[s]]$sample_class),
                      labels.1 = c(l1), labels.2 = c(l2),
                      plot.embedding = Embeddings(object = cells[[s]],
                                                  reduction = "umap"))

    num_da_cells <- length(c(res$da.down, res$da.up))
    proportion_da <- num_da_cells / length(colnames(cells[[s]]))
    print(paste0("Proportion of cells in differentially abundant regions at ",
                 num_days, "D: ", proportion_da))

    daseq_umap <- create_daseq_umap(res)
    print(daseq_umap)
    if (save_images) {
        save_tiff(daseq_umap, out_name)
    }
}
```

```
## Calculating DA score vector.
## Running GLM.
## Test on random labels.
## Setting thresholds based on permutation
## [1] "Proportion of cells in differentially abundant regions at 6D: 0.902621278614636"
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

```
## Calculating DA score vector.
## Running GLM.
## Test on random labels.
## Setting thresholds based on permutation
## [1] "Proportion of cells in differentially abundant regions at 12D: 0.371958285052144"
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-9-2.png" style="display: block; margin: auto;" />

```r
(daseq_legend <- as_ggplot(get_legend(res$pred.plot +
                                      umap_theme +
                                      theme(legend.position = "top",
                                            legend.text = element_blank()))))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-9-3.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(daseq_legend, "daseq_legend.tiff")
}
```

```
## quartz_off_screen 
##                 2
```

# Differential Expression Analysis

We carried out differential expression analysis of treated and control cells at
6D and 12D using the standard Wilcoxon rank sum test in Seurat. 


```r
vol_plot <- function(marks, not_sig_color) {
    marks$p_val_adj_neg_log <- (-1) * log10(marks$p_val_adj)
    marks$p_val_adj_neg_log[is.infinite(marks$p_val_adj_neg_log)] <- 300
    sig_marks_pos <- subset(marks, avg_log2FC > 1)
    sig_marks_neg <- subset(marks, avg_log2FC < -1)

    # Assign classes to genes based on how DE they are. Used to set colours
    de_class <- ifelse(marks$avg_log2FC < -1 & marks$p_val_adj < 0.01,
                       "Downreg",
                       ifelse(marks$avg_log2FC > 1 & marks$p_val_adj < 0.01,
                              "Upreg",
                              ifelse(marks$p_val_adj < 0.01,
                                     "Low-LFC", "Not-Sig")))
    de_class <- factor(de_class, levels = c("Not-Sig", "Low-LFC", "Downreg",
                                            "Upreg"))

    x_sublabel_y <- -50 # Position the sub-labels of the x-axis correctly
    plot <- ggplot(marks, aes(x = avg_log2FC, y = p_val_adj_neg_log,
                              label = rownames(marks))) +
            geom_point(aes(color = de_class, fill = de_class), size = 2,
                       shape = 21, alpha = 0.5) +
            scale_color_manual(values = c("gray", not_sig_color, "red",
                                          "white")) +
            scale_fill_manual(values = c("gray", not_sig_color, "white",
                                         "red")) +
            geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
            geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
            xlim(-2.5, 2.5) + ylim(-log10(0.01), 350) +
            xlab(bquote(bold(Log["2"]~"fold change"))) +
            ylab(bquote(bold(-Log["10"] ~ P))) +
            theme(legend.position = "none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            axis_text_theme +
            coord_cartesian(clip = "off")

    return(plot)
}

Idents(cells$all) <- cells$all$sample_class
not_sig_color <- "#4D4D4D"

for (i in list(c("ROCKi-6D", "Control-6D", "fig_5e.tiff"),
               c("ROCKi-12D", "Control-12D", "fig_5f.tiff"))) {
    id1 <- i[1]
    id2 <- i[2]
    out_name <- i[3]

    marks <- FindMarkers(cells$all, ident.1 = id1, ident.2 = id2)
    vol <- vol_plot(marks, not_sig_color)
    print(vol)

    if (save_images) {
        save_tiff(vol, out_name, height = 5, width = 10)
    }
}
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-10-1.png" width="100%" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-10-2.png" width="100%" style="display: block; margin: auto;" />

```r
labels <- c(bquote(bold(Log["2"]~"fold change ≤ -1")),
            bquote(bold("|"*Log["2"]~"fold change| < 1")),
            bquote(bold(Log["2"]~"fold change ≥ 1")))

legend_plot <- ggplot(data.frame(class = as.character(seq_len(length(labels))),
                                 x = seq_len(length(labels))),
                      aes(x = x, y = x, color = class, fill = class)) +
               geom_point(size = 3, shape = 21, alpha = 0.5) +
               scale_fill_manual(values = c("white", not_sig_color, "red"),
                                 labels = labels) +
               scale_color_manual(values = c("red", not_sig_color, "white"),
                                  labels = labels) +
               theme(legend.title = element_blank(),
                     legend.key = element_rect(fill = NA, color = NA),
                     legend.text = element_text(face = "bold"),
                     legend.margin = margin(unit(c(0, 0, 0, 0), "cm"))) +
               guides(colour = guide_legend(override.aes = list(size = 2)))

(volplot_legend <- as_ggplot(get_legend(legend_plot, "right")))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-10-3.png" width="100%" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(volplot_legend, "volplot_legend.tiff")
}
```

# Cell Type Proportions in Treated / Untreated

We used SCINA[^1] to estimate the cell type proportions in the treated/control
cells. We used cell type markers published by Enzo et al._et al_ [^2] to
classify cells as either basal, transient amplifying (TA), or terminal
differentiated (TD), since these are the classes for which Enzo et al. provide
distinct markers.


```r
signatures <- preprocess.signatures("celltype_markers_basal_ta_td.csv")
print(signatures)
```

```
## $Basal
## [1] "AVLN"   "AURKB"  "CCNA2"  "CKAP2L" "FOXM1"  "HMGB2"  "LMNB1" 
## 
## $TA
## [1] "KRT14" "TP63"  "ITGA6" "ITGB1" "BIRC5"
## 
## $TD
## [1] "SERPINB3" "SFN"      "KRT10"    "TGM1"     "IVL"      "SPINK5"
```

```r
for (s in c("all", "d6", "d12")) {
    preds <- SCINA(as.matrix(GetAssayData(cells[[s]])), signatures)$cell_labels
    preds <- factor(preds, levels = c("Basal", "TA", "TD", "unknown"))
    cells[[s]]$cell_type_pred <- preds
}

cells$ctrl6d <- subset(cells$d6, treatment_class == "Control")
cells$rocki6d <- subset(cells$d6, treatment_class == "ROCKi")
cells$ctrl12d <- subset(cells$d12, treatment_class == "Control")
cells$rocki12d <- subset(cells$d12, treatment_class == "ROCKi")

create_cell_type_umap <- function(cells) {
    # Force factors back into correct order as this inexplicably disappears
    # when passed to a function
    cells$cell_type_pred <- factor(cells$cell_type_pred,
                                   levels = c("Basal", "TA", "TD", "unknown"))
    p <- DimPlot(cells, group.by = "cell_type_pred") +
         scale_color_manual(labels = c("Holoclone-forming",
                                       "Mero- or Paraclone-forming",
                                       "Differentiated",
                                       "Unclassified keratinocytes"),
                            values = c("#80C9EA", "#DD6E79", "#43863E",
                                       "#BBBBBB")) +
         coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

for (i in list(c("ctrl6d", "fig_6a.tiff"), c("rocki6d", "fig_6c.tiff"),
               c("ctrl12d", "fig_6b.tiff"), c("rocki12d", "fig_6d.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    cell_type_umap <- create_cell_type_umap(cells[[s]])
    print(cell_type_umap)
    if (save_images) {
        save_tiff(cell_type_umap, out_name)
    }
}
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-2.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-3.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-4.png" style="display: block; margin: auto;" />

```r
cell_type_prop_graph <- function(df, names) {
    p <- ggplot(df, aes(x = Class, fill = CellType)) +
         geom_bar(position = "fill") + ylab("Proportion of cell type") +
         labs(fill = "Cell type") +
         scale_fill_manual(labels = c("Holoclone-forming",
                                      "Mero- or Paraclone-forming",
                                      "Differentiated",
                                      "Unclassified keratinocytes"),
                           values = c("#80C9EA", "#DD6E79", "#43863E",
                                      "#BBBBBB")) +
         scale_x_discrete(labels = names) +
         chart_margin_theme + axis_text_theme +
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               axis.title.x = element_blank())

    return(p)
}

for (i in list(c("d6", 6, "fig_6e.tiff"), c("d12", 12, "fig_6f.tiff"))) {
    s <- i[1]
    num_days <- i[2]
    out_name <- i[3]

    s1 <- paste0("Control-", num_days, "D")
    s2 <- paste0("ROCKi-", num_days, "D")
    cell_type_df <- data.frame(CellType = cells[[s]]$cell_type_pred,
                               Class = cells[[s]]$sample_class)
    cell_type_df$Class <- factor(cell_type_df$Class,
                                 levels = c(s1, s2))
    prop_graph <- cell_type_prop_graph(cell_type_df,
                                       c(paste0(s1, " "), paste0(s2, " ")))
    prop_graph_no_legend <- prop_graph + theme(legend.position = "none")
    print(prop_graph_no_legend)

    if (save_images) {
        save_tiff(prop_graph_no_legend, out_name, width = 6, height = 6)
    }
}
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-5.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-6.png" style="display: block; margin: auto;" />

We can see that at the 6D point there is a large difference in the proportion
of cell types, whereas this difference reduces at 12D.

## Generating confidence intervals

We now compute 95% confidence intervals around the differences in cell type
proportions at both timepoints.


```r
prop_table_cell_type <- function(x) {
    return(prop.table(table(factor(x, levels = c("Basal", "TA", "TD",
                                                 "unknown")))))
}

calc_cell_type_props <- function(cell_types, treatment_classes) {
    # Calculate the proportion of each cell type for both treatment classes.
    df <- data.frame(type = cell_types, treatment = treatment_classes)

    aggred_means <- aggregate(df$type, list(df$treatment),
                              prop_table_cell_type)
    for (i in 1:2) {
        aggred_means[, 2][i, ] <- aggred_means[, 2][i, ] /
                                  sum(aggred_means[, 2][i, ][1:3])
    }

    proportions <- c(aggred_means[1, ][[2]][1], aggred_means[1, ][[2]][2],
                     aggred_means[1, ][[2]][3], aggred_means[2, ][[2]][1],
                     aggred_means[2, ][[2]][2], aggred_means[2, ][[2]][3])

    return(sum(abs(proportions[1:2] - proportions[4:5])))
}

# Set random seed for reproducibility
set.seed(1234)

for (i in list(c("d6", 6), c("d12", 12))) {
    s <- i[1]
    num_days <- i[2]

    obs_total_prop_diff <- calc_cell_type_props(cells[[s]]$cell_type_pred,
                                                cells[[s]]$treatment_class)
    prop_diffs <- c()
    for (r in 1:10000) {
        # Generate resample indices
        c_idx <- sample(which(cells[[s]]$treatment_class == "Control"),
                          length(which(cells[[s]]$treatment_class == "Control")),
                          replace = TRUE)
        r_idx <- sample(which(cells[[s]]$treatment_class == "ROCKi"),
                          length(which(cells[[s]]$treatment_class == "ROCKi")) - 1,
                          replace = TRUE)
        idx <- c(c_idx, r_idx)

        # Find cell type and treatment class labels for this resample
        cell_types <- cells[[s]]$cell_type_pred[idx]
        treat_classes <- cells[[s]]$treatment_class[idx]

        total_prop_diff <- calc_cell_type_props(cell_types, treat_classes)
        prop_diffs <- c(prop_diffs, total_prop_diff)
    }

    alpha_vals <- c(0.975, 0.025)
    residuals <- prop_diffs - obs_total_prop_diff
    print(paste0("95% confidence interval for total absolute difference in",
                 "cell type proportions at ", num_days, "D:"))
    print(as.numeric(obs_total_prop_diff - quantile(residuals, alpha_vals)))
}
```

```
## [1] "95% confidence interval for total absolute difference incell type proportions at 6D:"
## [1] 0.1808670 0.2551785
## [1] "95% confidence interval for total absolute difference incell type proportions at 12D:"
## [1] -0.009402354  0.022321367
```


```r
create_dotplot <- function(cells, signatures) {
    dp <- DotPlot(cells, group.by = "sample_class_space",
                  features = list("Holoclone" = signatures$Basal,
                                  "Mero-/Paraclone" = signatures$TA,
                                  "Differentiated" = signatures$TD)) +
          axis_text_theme +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text = element_text(size = 10),
                axis.title = element_blank(),
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 10, face = "bold"))

    return(dp)
}

cells$all$sample_class_space <- cells$all$sample_class
levels(cells$all$sample_class_space) <- c("Control-6D ", "ROCKi-6D ",
                                          "Control-6D 6D ", "ROCKi-6D 6D ")

cells$holo <- subset(cells$all, subset = cell_type_pred == "Basal")
cells$ta <- subset(cells$all, subset = cell_type_pred == "TA")
cells$td <- subset(cells$all, subset = cell_type_pred == "TD")

for (i in list(c("holo", "fig_6g.tiff"), c("ta", "fig_s2a.tiff"),
               c("td", "fig_s2b.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    dp <- create_dotplot(cells[[s]], signatures)
    print(dp)

    if (save_images) {
        save_tiff(dp, out_name, width = 20, height = 8)
    }
}
```

# Trajectory Inference and Comparison

We used Slingshot[^3] to infer the differentiation trajectory for each sample
class. We used the UMAP reduction for trajectory inference as we found that
using the reduction produced by Harmony resulted in unrealistically similar
pseudotimes among the majority of the cells for each class, whereas using the
UMAP reduction resulted in a more constant increase in pseudotime, which we
considered more realistic. We used three clusters for each sample class (that
can be thought of as corresponding to basal, TA, and TD cell types) as an input
to Slingshot. We used the cluster likely to contain mostly basal cells as the
starting point for each inferred trajectory.


```r
cells$ctrl6d <- subset(cells$d6, treatment_class == "Control")
cells$rocki6d <- subset(cells$d6, treatment_class == "ROCKi")
cells$ctrl12d <- subset(cells$d12, treatment_class == "Control")
cells$rocki12d <- subset(cells$d12, treatment_class == "ROCKi")

for (s in c("ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    cells[[s]] <- RunUMAP(cells[[s]], dims = 1:50, reduction = "harmony")
}

prepare_cells_for_slingshot <- function(cells) {
    cells <- FindNeighbors(cells, reduction = "harmony", dims = 1:50)
    cells <- FindClusters(cells, resolution = 0.1)
    if (all(cells$sample_class == "Control-6D")) {
        cells <- RenameIdents(object = cells,
                              "0" = "TA", "1" = "Basal", "2" = "TD")
    } else {
        cells <- RenameIdents(object = cells,
                              "0" = "TA", "1" = "TD", "2" = "Basal")
    }

    return(cells)
}

for (s in c("d6", "d12", "ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    cells[[s]] <- prepare_cells_for_slingshot(cells[[s]])
}
```

## Computing Trajectories


```r
# Change UMAP values to make figures all similarly oriented
for (s in c("ctrl12d")) {
    flipped <- (-1) * Embeddings(object = cells[[s]], reduction = "umap")[, 1]
    cells[[s]]@reductions$umap@cell.embeddings[, 1] <- flipped
}
for (s in c("ctrl6d", "rocki6d")) {
    flipped <- (-1) * Embeddings(object = cells[[s]], reduction = "umap")[, 2]
    cells[[s]]@reductions$umap@cell.embeddings[, 2] <- flipped
}

colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)

for (s in c("d6", "d12", "ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    cells[[paste0(s, "_sce")]] <- as.SingleCellExperiment(cells[[s]])
}

sshots <- list()

for (s in c("d6", "d12")) {
    sshots[[s]] <- slingshot(cells[[paste0(s, "_sce")]], Idents(cells[[s]]),
                             reducedDim = "HARMONY", start.clus = "Basal")
}

for (s in c("ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    sshots[[s]] <- slingshot(cells[[paste0(s, "_sce")]], Idents(cells[[s]]),
                             reducedDim = "UMAP", start.clus = "Basal")
}

umap_curves <- list()

for (s in c("d6", "d12")) {
    umap_embed <- reducedDims(cells[[paste0(s, "_sce")]])$UMAP
    umap_curves[[s]] <- slingCurves(embedCurves(sshots[[s]], umap_embed))[[1]]
}

for (s in c("ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    umap_curves[[s]] <- slingCurves(sshots[[s]])[[1]]
}

for (num_days in c(6, 12)) {
    s <- paste0("d", num_days)
    ctrl_idxs <- which(cells[[s]]$treatment_class == "Control")
    rocki_idxs <- which(cells[[s]]$treatment_class == "ROCKi")
    sshots[[paste0("ctrl_", num_days, "d_together")]] <- sshots[[s]][, ctrl_idxs]
    sshots[[paste0("rocki_", num_days, "d_together")]] <- sshots[[s]][, rocki_idxs]
}

xmin_sep <- min(reducedDim(cells$ctrl6d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 1],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 1])
xmax_sep <- max(reducedDim(cells$ctrl6d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 1],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 1])
ymin_sep <- min(reducedDim(cells$ctrl6d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 2],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 2])
ymax_sep <- max(reducedDim(cells$ctrl6d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 2],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 2])

create_traj_umap <- function(sshot_to_plot, umap_curve, cell_types) {
    points_df <- data.frame(u1 = reducedDims(sshot_to_plot)$UMAP[, 1],
                            u2 = reducedDims(sshot_to_plot)$UMAP[, 2],
                            cell_types = cell_types)
    points_df$cell_types <- factor(points_df$cell_types,
                                   levels = c("Basal", "TA", "TD", "unknown"))
    curve_df <- data.frame(u1 = umap_curve$s[umap_curve$ord, ][, 1],
                           u2 = umap_curve$s[umap_curve$ord, ][, 2])

    arrowhead_length <- unit(0.2, "cm")
    arrow_width <- 1

    ptimes <- sshot_to_plot$slingPseudotime_1
    basal_ptimes <- ptimes[which(cell_types == "Basal")]
    ta_ptimes <- ptimes[which(cell_types == "TA")]
    td_ptimes <- ptimes[which(cell_types == "TD")]

    percentile <- 0.97
    basal_ptimes_dist <- abs(quantile(basal_ptimes, percentile) - basal_ptimes)
    basal_percentile_ptime_idx <- which.min(basal_ptimes_dist)
    ta_ptimes_dist <- abs(quantile(ta_ptimes, percentile) - ta_ptimes)
    ta_percentile_ptime_idx <- which.min(ta_ptimes_dist)
    td_ptimes_dist <- abs(quantile(td_ptimes, percentile) - td_ptimes)
    td_percentile_ptime_idx <- which.min(td_ptimes_dist)

    basal_percentile_ptime <- basal_ptimes[basal_percentile_ptime_idx]
    ta_percentile_ptime <- ta_ptimes[ta_percentile_ptime_idx]
    td_percentile_ptime <- td_ptimes[td_percentile_ptime_idx]

    end_idx <- round(length(sshot_to_plot$slingPseudotime_1) * percentile)
    basal_pct_overall_idx <- which(sort(ptimes) == basal_percentile_ptime)
    ta_pct_overall_idx <- which(sort(ptimes) == ta_percentile_ptime)

    p <- ggplot() +
         geom_point(data = points_df, aes(x = u1, y = u2, color = cell_types),
                    size = -0.1) +
         scale_color_manual(labels = c("Holoclone-forming",
                           "Mero- or Paraclone-forming",
                           "Differentiated", "Unknown"),
                           values = c("#80C9EA", "#DD6E79",
                                      "#43863E", "#BBBBBB")) +
         geom_path(data = curve_df[seq_len(end_idx), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[1:2, ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "open", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[seq_len(basal_pct_overall_idx), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[seq_len(ta_pct_overall_idx), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         coord_fixed() + umap_theme

     return(p)
}

for (i in list(c("ctrl6d", "fig_7a.tiff"), c("rocki6d", "fig_7b.tiff"),
               c("ctrl12d", "fig_7c.tiff"), c("rocki12d", "fig_7d.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    traj_umap <- create_traj_umap(sshots[[s]], umap_curves[[s]],
                                   Idents(cells[[s]])) +
                  theme(legend.position = "none")
    print(traj_umap)

    if (save_images) {
        save_tiff(traj_umap, out_name, width = 10, height = 8)
    }
}
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-15-2.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-15-3.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-15-4.png" style="display: block; margin: auto;" />

## Comparing trajectories

We compared trajectories by sampling 2500 pseudotimes uniformly at random for
each sample class, sorting them, and then plotting them against one another. We
also computed the “speed” of the differentiation at each point by calculating
the difference between each neighbouring pair of ordered pseudotimes and
plotting the moving average of 151 pseudotimes, in order to reduce noise. We
take cell type of a cell to be the cluster to which it belongs to improve the
appearance of the trajectory plots when colouring by cell type.


```r
for (s in c("d6", "d12")) {
    cells[[s]]$cell_type <- Idents(cells[[s]])
}

# mov_avg function copied from stackoverflow.com/a/4862334/6914552
mov_avg <- function(x, n = 151) {
    stats::filter(x, rep(1 / n, n), sides = 2)
}

create_traj_df <- function(ptimes, name, cell_types) {
    quantile_points <- seq(0.001, 1, 0.001)
    quantiled_data <- as.numeric(quantile(ptimes, quantile_points))
    sorted_cell_types <- cell_types[order(ptimes)]
    num_points <- length(ptimes)
    cell_indices <- round(num_points * quantile_points)
    cell_types_subset <- sorted_cell_types[cell_indices]

    return(data.frame(pct = quantile_points,
                      ptime = quantiled_data,
                      cell_type = factor(cell_types_subset,
                                         levels = c("Basal", "TA", "TD",
                                                    "unknown")),
                      diff = c(NA, mov_avg(diff(quantiled_data))),
                      sample_class = rep(name, length(quantile_points))))
}

traj_plot_margin_theme <- theme(plot.margin = margin(0.1, 1.2, 0.1, 1.2, "cm"))
diff_comparison_plots_theme <- theme(legend.title = element_blank())

point_size <- 1
create_traj_plot <- function(traj_df) {
    shade_alpha <- 0.2

    basal_upper_lim <- quantile(traj_df[traj_df$cell_type == "Basal", ]$ptime,
                                0.75)
    ta_upper_lim <- quantile(traj_df[traj_df$cell_type == "TA", ]$ptime, 0.95)

    basal_col <- "#80C9EA"
    ta_col <- "#DD6E79"
    td_col <- "#43863E"

    traj_plot <- ggplot(traj_df,
                        aes(x = pct, y = ptime, group = sample_class)) +
                 annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
                          ymax = basal_upper_lim, alpha = shade_alpha,
                          fill = basal_col) +
                 annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = basal_upper_lim, ymax = ta_upper_lim,
                          alpha = shade_alpha, fill = ta_col) +
                 annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = ta_upper_lim, ymax = Inf,
                          alpha = shade_alpha, fill = td_col) +
                 geom_point(aes(color = cell_type), size = point_size) +
                 scale_color_manual(values = c(basal_col, ta_col, td_col)) +
                 new_scale_color() +
                 geom_line(aes(linetype = sample_class), size = 0.4) +
                 scale_linetype_manual(values = c("solid", "dashed")) +
                 scale_x_continuous(labels = scales::percent,
                                    limits = c(0, 1)) +
                 ylim(c(0, max(traj_df$ptime))) +
                 xlab("Percentage through differentation process") +
                 ylab("Pseudotime") + traj_plot_margin_theme +
                 diff_comparison_plots_theme + axis_text_theme

    return(traj_plot)
}

create_speed_plot <- function(traj_df) {
    speed_plot <- ggplot(traj_df,
                         aes(x = pct, y = diff, group = sample_class)) +
                  geom_point(aes(color = cell_type), size = point_size) +
                  scale_color_manual(values = c("#80C9EA", "#DD6E79",
                                                "#43863E")) +
                  new_scale_color() +
                  geom_line(aes(linetype = sample_class), size = 0.4) +
                  scale_linetype_manual(values = c("solid", "dashed")) +
                  scale_x_continuous(labels = scales::percent,
                                     limits = c(0, 1)) +
                  xlab("Percentage through differentation process") +
                  ylab("Differentiation speed") + traj_plot_margin_theme +
                  diff_comparison_plots_theme + axis_text_theme

    return(speed_plot)
}

obs_results_ptime <- list()
obs_results_speed <- list()
for (i in list(c(6, "fig_7e.tiff", "fig_s3a.tiff"),
               c(12, "fig_7f.tiff", "fig_s3b.tiff"))) {
    num_days <- i[1]
    traj_graph_out_name <- i[2]
    speed_graph_out_name <- i[3]
    s <- paste0("d", num_days)
    ctrl_idxs <- which(cells[[s]]$treatment_class == "Control")
    rocki_idxs <- which(cells[[s]]$treatment_class == "ROCKi")
    ctrl_cell_types <- cells[[paste0("d", num_days)]][, ctrl_idxs]$cell_type
    rocki_cell_types <- cells[[paste0("d", num_days)]][, rocki_idxs]$cell_type
    cp <- sshots[[paste0("ctrl_", num_days, "d_together")]]$slingPseudotime_1
    rp <- sshots[[paste0("rocki_", num_days, "d_together")]]$slingPseudotime_1

    traj_df_ctrl <- create_traj_df(cp, paste0("Control-", num_days, "D"),
                                   ctrl_cell_types)
    traj_df_rocki <- create_traj_df(rp, paste0("ROCKi-", num_days, "D"),
                                    rocki_cell_types)
    traj_df <- rbind(traj_df_ctrl, traj_df_rocki)

    ptime_abs_diffs <- abs(traj_df_ctrl$ptime - traj_df_rocki$ptime)
    obs_results_ptime[[paste0("d", num_days)]] <- mean(ptime_abs_diffs)

    speed_abs_diffs <- abs(traj_df_ctrl$diff - traj_df_rocki$diff)
    obs_results_speed[[paste0("d", num_days)]] <- mean(speed_abs_diffs,
                                                      na.rm = TRUE)

    traj_plot <- create_traj_plot(traj_df) + theme(legend.position = "none")
    print(traj_plot)
    speed_plot <- create_speed_plot(traj_df) + theme(legend.position = "none")
    print(speed_plot)

    if (save_images) {
        save_tiff(traj_plot, traj_graph_out_name, width = 10, height = 8)
        save_tiff(speed_plot, speed_graph_out_name, width = 14, height = 8)
    }
}
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-16-1.png" width="100%" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-16-2.png" width="100%" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-16-3.png" width="100%" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-16-4.png" width="100%" style="display: block; margin: auto;" />

```r
legend_plot <- ggplot(traj_df, aes(x = pct, y = ptime, group = sample_class)) +
               geom_point(aes(color = cell_type), size = 3) +
               scale_color_manual(values = c("#80C9EA", "#DD6E79", "#43863E"),
                                  labels = c("Holoclone-forming",
                                             "Mero- or Paraclone-forming",
                                             "Differentiated")) +
               theme(legend.title = element_blank(),
                     legend.text = element_text(face = "bold"))
(dot_colour_legend <- as_ggplot(get_legend(legend_plot, "right")))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-16-5.png" width="100%" style="display: block; margin: auto;" />

```r
legend_plot <- ggplot(traj_df, aes(x = pct, y = ptime, group = sample_class)) +
               geom_line(aes(linetype = sample_class), size = 0.4) +
               scale_linetype_manual(values = c("solid", "dashed"),
                                     labels = c("Control", "ROCKi-treated")) +
               theme(legend.title = element_blank(),
                     legend.text = element_text(face = "bold"))
(linetype_legend <- as_ggplot(get_legend(legend_plot, "right")))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-16-6.png" width="100%" style="display: block; margin: auto;" />

```r
(legends <- ggarrange(linetype_legend, dot_colour_legend, ncol = 2))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-16-7.png" width="100%" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(legends, "fig_7_legends.tiff", width = 10, height = 8)
}
```

## Generating confidence intervals


```r
# Generate bootstrapped data as in "experiment_code_trajectories.R" and save
# the outputs as "trajectory_comparison_results.txt"

# bs_out = "bootstrap out"
bs_out <- read.table("trajectory_comparison_results.txt")
bs_out <- bs_out[, c(3, 4, 5, 6)]
colnames(bs_out) <- c("TrajDiff6D", "TrajDiff12D", "SpeedDiff6D",
                      "SpeedDiff12D")

for (num_days in c(6, 12)) {
    traj_diffs <- bs_out[paste0("TrajDiff", num_days, "D")][, 1]
    obs_result_ptime <- obs_results_ptime[[paste0("d", num_days)]]
    diffs_quants_traj <- quantile((traj_diffs - obs_result_ptime), alpha_vals)
    residuals_traj <- as.numeric(obs_result_ptime - diffs_quants_traj)
    print(paste0("95% confidence interval for difference in pseudotimes at ",
                 num_days, "D: ", paste(residuals_traj, collapse = " ")))

    speed_diffs <- bs_out[paste0("SpeedDiff", num_days, "D")][, 1]
    obs_result_speed <- obs_results_speed[[paste0("d", num_days)]]
    diffs_quants_speed <- quantile((speed_diffs - obs_result_speed), alpha_vals)
    residuals_speed <- as.numeric(obs_result_speed - diffs_quants_speed)
    print(paste0("95% confidence interval for difference in speed at ",
                 num_days, "D: ", paste(residuals_speed, collapse = " ")))
}
```

```
## [1] "95% confidence interval for difference in pseudotimes at 6D: 12.4460321784284 15.5195479079284"
## [1] "95% confidence interval for difference in speed at 6D: 0.0366647847933335 0.0519434975433335"
## [1] "95% confidence interval for difference in pseudotimes at 12D: 1.65452072266176 3.97100868716176"
## [1] "95% confidence interval for difference in speed at 12D: 0.0146296286652594 0.0317649334152594"
```

# Supplementary figures

## Expression of VIM in fibroblasts


```r
cells$all_with_outliers <- cells$all_with_outliers %>%
                           SCTransform() %>%
                           CellCycleScoring(s.features = s_genes,
                                            g2m.features = g2m_genes)
cells$all_with_outliers$CC.difference <- cells$all_with_outliers$S.Score -
                                         cells$all_with_outliers$G2M.Score
cells$all_with_outliers <- cells$all_with_outliers %>%
                           SCTransform(vars.to.regress = c("CC.difference")) %>%
                           RunPCA(features = VariableFeatures(cells$all_with_outliers,
                                                              assay = "SCT"),
                                  npcs = 50)
donor_idents <- sapply(cells$all_with_outliers$orig.ident,
                       function(x) {
                           if (x %in% c("1", "2", "3", "4")) {
                               "1"
                           } else if (x %in% c("5", "6", "7", "8")) {
                               "2"
                           } else {
                               "3"
                           }
                       })
cells$all_with_outliers$donor_idents <- donor_idents
cells$all_with_outliers <- cells$all_with_outliers %>%
                           RunHarmony("donor_idents", assay.use = "SCT") %>%
                           RunUMAP(dims = 1:50, reduction = "harmony")
(p <- FeaturePlot(cells$all_with_outliers, "VIM"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(p, "fig_s4a.tiff", width = 10, height = 8)
}
```

## Co-expression of VIM and KRT14


```r
cells$all <- RunUMAP(cells$all, reduction = "harmony", dims = 1:50)
p <- FeaturePlot(cells$all, features = c("KRT14", "VIM"), blend = TRUE,
                 split.by = "sample_class")
blended_fp <- list()
for (i in 1:4) {
    q <- ggarrange(p[[((i - 1) * 4) + 1]] + umap_theme,
                   p[[((i - 1) * 4) + 2]] + umap_theme,
                   p[[((i - 1) * 4) + 3]] + umap_theme,
                   p[[((i - 1) * 4) + 4]] + theme(plot.title = element_blank()),
                   nrow = 1, widths = c(2, 2, 2, 3))
    blended_fp[[i]] <- q
    print(q)
}
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-19-2.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-19-3.png" style="display: block; margin: auto;" /><img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-19-4.png" style="display: block; margin: auto;" />

```r
save_tiff(ggarrange(plotlist = blended_fp, ncol = 1), "fig_s4b.tiff",
          width = 15, height = 15)
```


# References {#endnotes}

[^1]: Zhang Z, Luo D, Zhong X, et al. _SCINA: A Semi-Supervised Subtyping Algorithm of Single Cells and Bulk Samples_. Genes (Basel). 2019;10(7):531.

[^2]: Enzo, E., Secone Seconetti, A., Forcato, M. et al. _Single-keratinocyte transcriptomic analyses identify different clonal types and proliferative potential mediated by FOXM1 in human epidermal stem cells_. Nat Commun 12, 2505 (2021).

[^3]: Street, K., Risso, D., Fletcher, R. B., et al. _Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics_. BMC Genomics 19(477), 2018.

# Appendix

## All code for this report

<!-- Suggested in https://bookdown.org/yihui/rmarkdown-cookbook/code-appendix.html -->

```r
knitr::opts_chunk$set(results = "hide", message = FALSE, warning = FALSE,
                      fig.align = "center")
# Tidyverse
library(dplyr)
library(tidyr)

# Figure generation
library(ggplot2)
library(ggrepel)
library(grid)
library(png)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(ggnewscale)

# Analysis tools
library(Seurat)
library(sctransform)
library(harmony)
library(kBET)
library(SCINA)
library(slingshot)
library(SingleCellExperiment)
library(DAseq)

# Set ggplot theme
theme_update(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
theme_update(axis.text = element_text(size = 10, face = "plain"))
theme_update(axis.title = element_text(size = 10, face = "plain"))

axis_text_theme <- theme(axis.title = element_text(face = "bold",
                                                   color = "black"),
                         axis.text = element_text(face = "bold",
                                                  color = "black"))

blue_color <- "#6E86FF"
red_color <- "#FF8585"

save_images <- TRUE
save_tiff <- function(img, img_name, units = "cm", width = 10, height = 10,
                      res = 300) {
    img_path <- paste0("paper_figs/", img_name)
    tiff(img_path, units = units, width = width, height = height, res = res)
    print(img)
    dev.off()
}
set.seed(10403) # specify default R random seed (for reproducibility)

skin_data <- Read10X(data.dir = "3donor_filtered_feature_bc_matrix")
cells <- list()
cells$all <- CreateSeuratObject(counts = skin_data,
                                project = "3DonorEffectOfRocki",
                                min.cells = 3, min.features = 200,
                                names.delim = "-", names.field = 2)
cells$all <- RenameIdents(object = cells$all, "1" = "1Control-6D ",
                          "2" = "1ROCKi-6D ", "3" = "1Control-6D 6D ",
                          "4" = "1ROCKi-6D 6D ", "5" = "2Control-6D ",
                          "6" = "2ROCKi-6D ", "7" = "2Control-6D 6D ",
                          "8" = "2ROCKi-6D 6D ", "9" = "3Control-6D ",
                          "10" = "3ROCKi-6D ", "11" = "3Control-6D 6D ",
                          "12" = "3ROCKi-6D 6D ")
# Compute percentage of mitochondrial DNA in each cell
cells$all[["percent.mt"]] <- PercentageFeatureSet(cells$all, pattern = "^MT-")

vlnplot_limits <- function(data, metric, ylab, limits) {
    df <- data.frame(x = Idents(data), y = data[[metric]][, 1])
    plot <- ggplot() + geom_violin(data = df, aes(x = x, y = y, fill = x)) +
            geom_errorbar(data = limits, aes(x = x, ymin = lower, ymax = upper),
                          color = "red") +
            labs(y = ylab) +
            scale_y_continuous(labels = scales::comma) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "none", plot.title = element_blank(),
                  axis.title.x = element_blank()) + axis_text_theme

    return(plot)
}

ncount_df <- data.frame(x = levels(Idents(cells$all)),
                        lower = c(8000, 15000, 8000, 8000,
                                  8000, 8000, 5000, 8000,
                                  8000, 8000, 5000, 5000),
                        upper = c(50000, 62500, 50000, 50000,
                                  50000, 50000, 50000, 50000,
                                  50000, 50000, 40000, 40000))

(ncount_vln <- vlnplot_limits(cells$all, "nCount_RNA",
                              "Number of RNA\nmolecules per cell",
                              ncount_df) + ylim(c(0, 100000)))
if (save_images) {
    save_tiff(ncount_vln, "fig_s4a.tiff")
}

nfeature_df <- data.frame(x = levels(Idents(cells$all)),
                     lower = c(2500, 3000, 2500, 2500,
                               2500, 2500, 2000, 2500,
                               2500, 2500, 1500, 1500),
                     upper = c(7500, 8500, 7500, 7500,
                               7500, 7500, 7500, 7500,
                               7500, 7500, 6000, 6000))

(nfeature_vln <- vlnplot_limits(cells$all, "nFeature_RNA",
                                "Number of genes\ndetected per cell",
                                nfeature_df))
if (save_images) {
    save_tiff(nfeature_vln, "fig_s4b.tiff")
}

pct_mito_df <- data.frame(x = levels(Idents(cells$all)),
                          lower = c(rep(2, 8), rep(3, 4)),
                          upper = c(rep(10, 8), rep(12, 4)))

(pct_mito_vln <- vlnplot_limits(cells$all, "percent.mt",
                                paste0("Percentage of mitochrondrial\n",
                                       "DNA detected per cell"),
                                pct_mito_df) + ylim(c(0, 20)))
if (save_images) {
    save_tiff(pct_mito_vln, "fig_s4c.tiff")
}

pass_qc <- c()
for (j in seq_len(length(cells$all$orig.ident))) {
    i <- as.integer(as.character(cells$all$orig.ident[j]))
    pass_qc <- c(pass_qc, all(cells$all$nFeature_RNA[j] >= nfeature_df$lower[i],
                              cells$all$nFeature_RNA[j] <= nfeature_df$upper[i],
                              cells$all$nCount_RNA[j] >= ncount_df$lower[i],
                              cells$all$nCount_RNA[j] <= ncount_df$upper[i],
                              cells$all$percent.mt[j] >= pct_mito_df$lower[i],
                              cells$all$percent.mt[j] <= pct_mito_df$upper[i]))
}
cells$all <- cells$all[, pass_qc]
cells$all_with_outliers <- cells$all

# Remove melanocytes
cells$all <- subset(cells$all, subset = MLANA > 1 | PMEL > 1 | MITF > 1,
                    invert = TRUE)
# Remove mesenchymal cells
cells$all <- subset(cells$all,
                    subset = MTRNR2L6 > 1 | MTRNR2L10 > 1 | MTRNR2L7 > 1 |
                             MTRNR2L1 > 1 | PRKAR2B > 1 | NR2F1 > 1,
                    invert = TRUE)
# Remove fibroblasts
cells$all <- subset(cells$all, subset = ACTG2 > 1 | DLK1 > 1, invert = TRUE)
cells$all <- SCTransform(cells$all)
cells$all_outliers_removed <- cells$all
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
cells$all <- CellCycleScoring(cells$all, s.features = s_genes,
                              g2m.features = g2m_genes)
cells$all$CC.difference <- cells$all$S.Score - cells$all$G2M.Score
cells$all <- SCTransform(cells$all, vars.to.regress = c("CC.difference"))
cells$all_cell_cycle_regressed <- cells$all
# Run PCA
cells$all <- RunPCA(cells$all,
                    features = VariableFeatures(cells$all, assay = "SCT"),
                    npcs = 50)
(elbow_plot <- ElbowPlot(cells$all, ndims = 50) + axis_text_theme)
if (save_images) {
    save_tiff(elbow_plot, "fig_s5.tiff")
}

# Remove inter-donor differences using Harmony
donor_idents <- sapply(cells$all$orig.ident,
                       function(x) {
                           if (x %in% c("1", "2", "3", "4")) {
                               "1"
                           } else if (x %in% c("5", "6", "7", "8")) {
                               "2"
                           } else {
                               "3"
                           }
                       })
cells$all$donor_idents <- donor_idents

cells$all <- cells$all %>% RunHarmony("donor_idents", assay.use = "SCT")
cells$all_harmonised <- cells$all

sample_class <- sapply(cells$all$orig.ident,
                       function(x) {
                           if (x %in% c("1", "5", "9")) {
                               "Control-6D"
                           } else if (x %in% c("2", "6", "10")) {
                               "ROCKi-6D"
                           } else if (x %in% c("3", "7", "11")) {
                               "Control-12D"
                           } else {
                               "ROCKi-12D"
                           }
                       })
cells$all$sample_class <- sample_class
cells$all$sample_class <- factor(cells$all$sample_class,
                                           levels = c("Control-6D", "ROCKi-6D",
                                                      "Control-12D",
                                                      "ROCKi-12D"))

treatment_class <- sapply(cells$all$sample_class,
                       function(x) {
                           if (x %in% c("Control-6D", "Control-12D")) {
                               "Control"
                           } else {
                               "ROCKi"
                           }
                       })
cells$all$treatment_class <- treatment_class
cells$d6 <- subset(cells$all,
                   subset = sample_class %in% c("Control-6D", "ROCKi-6D"))
cells$d12 <- subset(cells$all,
                    subset = sample_class %in% c("Control-12D", "ROCKi-12D"))

for (s in c("d6", "d12")) {
    cells[[s]] <- RunUMAP(cells[[s]], dims = 1:50, reduction = "harmony")
}

umap_theme <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                    axis.text.y = element_blank(), axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    plot.title = element_blank())
chart_margin_theme <- theme(plot.margin = margin(0, 1.2, 0, 1.2, "cm"))

labels <- c("Control", "ROCKi-treated")
legend_plot <- ggplot(data.frame(class = as.character(seq_len(length(labels))),
                                 x = seq_len(length(labels))),
                      aes(x = x, y = x, col = class)) +
               geom_point() +
               scale_color_manual(labels = labels,
                                  values = c(blue_color, red_color)) +
               theme(legend.title = element_blank(),
                     legend.key = element_rect(fill = NA, color = NA),
                     legend.margin = margin(c(0, 0, 0, 0)),
                     legend.text = element_text(face = "bold"),
                     legend.background = element_rect(fill = "#00000000",
                                                      size = 0.5,
                                                      linetype = "solid")) +
               guides(colour = guide_legend(override.aes = list(size = 2)))
(umap_legend <- as_ggplot(get_legend(legend_plot, "left")))
if (save_images) {
    save_tiff(umap_legend, "umap_legend.tiff")
}

plot_treatment_class_umap <- function(cells) {
    p <- DimPlot(cells, group.by = "treatment_class",
                 cols = c(blue_color, red_color),
                 pt.size = -0.1) +
         coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

for (i in list(c("d6", "fig_5a.tiff"), c("d12", "fig_5b.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    umap <- plot_treatment_class_umap(cells[[s]])
    print(umap)
    if (save_images) {
        save_tiff(umap, out_name)
    }
}
create_daseq_umap <- function(daseq_out) {
    daseq_umap <- daseq_out$pred.plot + umap_theme +
                  theme(legend.position = "none")

    return(daseq_umap)
}

for (i in list(c("d6", 6, "Control-6D", "ROCKi-6D", "fig_5c.tiff"),
               c("d12", 12, "Control-12D", "ROCKi-12D", "fig_5d.tiff"))) {
    s <- i[1]
    num_days <- i[2]
    l1 <- i[3]
    l2 <- i[4]
    out_name <- i[5]

    res <- getDAcells(Embeddings(object = cells[[s]], reduction = "harmony"),
                      cell.labels = as.character(cells[[s]]$sample_class),
                      labels.1 = c(l1), labels.2 = c(l2),
                      plot.embedding = Embeddings(object = cells[[s]],
                                                  reduction = "umap"))

    num_da_cells <- length(c(res$da.down, res$da.up))
    proportion_da <- num_da_cells / length(colnames(cells[[s]]))
    print(paste0("Proportion of cells in differentially abundant regions at ",
                 num_days, "D: ", proportion_da))

    daseq_umap <- create_daseq_umap(res)
    print(daseq_umap)
    if (save_images) {
        save_tiff(daseq_umap, out_name)
    }
}

(daseq_legend <- as_ggplot(get_legend(res$pred.plot +
                                      umap_theme +
                                      theme(legend.position = "top",
                                            legend.text = element_blank()))))

if (save_images) {
    save_tiff(daseq_legend, "daseq_legend.tiff")
}
vol_plot <- function(marks, not_sig_color) {
    marks$p_val_adj_neg_log <- (-1) * log10(marks$p_val_adj)
    marks$p_val_adj_neg_log[is.infinite(marks$p_val_adj_neg_log)] <- 300
    sig_marks_pos <- subset(marks, avg_log2FC > 1)
    sig_marks_neg <- subset(marks, avg_log2FC < -1)

    # Assign classes to genes based on how DE they are. Used to set colours
    de_class <- ifelse(marks$avg_log2FC < -1 & marks$p_val_adj < 0.01,
                       "Downreg",
                       ifelse(marks$avg_log2FC > 1 & marks$p_val_adj < 0.01,
                              "Upreg",
                              ifelse(marks$p_val_adj < 0.01,
                                     "Low-LFC", "Not-Sig")))
    de_class <- factor(de_class, levels = c("Not-Sig", "Low-LFC", "Downreg",
                                            "Upreg"))

    x_sublabel_y <- -50 # Position the sub-labels of the x-axis correctly
    plot <- ggplot(marks, aes(x = avg_log2FC, y = p_val_adj_neg_log,
                              label = rownames(marks))) +
            geom_point(aes(color = de_class, fill = de_class), size = 2,
                       shape = 21, alpha = 0.5) +
            scale_color_manual(values = c("gray", not_sig_color, "red",
                                          "white")) +
            scale_fill_manual(values = c("gray", not_sig_color, "white",
                                         "red")) +
            geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
            geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
            xlim(-2.5, 2.5) + ylim(-log10(0.01), 350) +
            xlab(bquote(bold(Log["2"]~"fold change"))) +
            ylab(bquote(bold(-Log["10"] ~ P))) +
            theme(legend.position = "none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            axis_text_theme +
            coord_cartesian(clip = "off")

    return(plot)
}

Idents(cells$all) <- cells$all$sample_class
not_sig_color <- "#4D4D4D"

for (i in list(c("ROCKi-6D", "Control-6D", "fig_5e.tiff"),
               c("ROCKi-12D", "Control-12D", "fig_5f.tiff"))) {
    id1 <- i[1]
    id2 <- i[2]
    out_name <- i[3]

    marks <- FindMarkers(cells$all, ident.1 = id1, ident.2 = id2)
    vol <- vol_plot(marks, not_sig_color)
    print(vol)

    if (save_images) {
        save_tiff(vol, out_name, height = 5, width = 10)
    }
}

labels <- c(bquote(bold(Log["2"]~"fold change ≤ -1")),
            bquote(bold("|"*Log["2"]~"fold change| < 1")),
            bquote(bold(Log["2"]~"fold change ≥ 1")))

legend_plot <- ggplot(data.frame(class = as.character(seq_len(length(labels))),
                                 x = seq_len(length(labels))),
                      aes(x = x, y = x, color = class, fill = class)) +
               geom_point(size = 3, shape = 21, alpha = 0.5) +
               scale_fill_manual(values = c("white", not_sig_color, "red"),
                                 labels = labels) +
               scale_color_manual(values = c("red", not_sig_color, "white"),
                                  labels = labels) +
               theme(legend.title = element_blank(),
                     legend.key = element_rect(fill = NA, color = NA),
                     legend.text = element_text(face = "bold"),
                     legend.margin = margin(unit(c(0, 0, 0, 0), "cm"))) +
               guides(colour = guide_legend(override.aes = list(size = 2)))

(volplot_legend <- as_ggplot(get_legend(legend_plot, "right")))
if (save_images) {
    save_tiff(volplot_legend, "volplot_legend.tiff")
}
signatures <- preprocess.signatures("celltype_markers_basal_ta_td.csv")
print(signatures)

for (s in c("all", "d6", "d12")) {
    preds <- SCINA(as.matrix(GetAssayData(cells[[s]])), signatures)$cell_labels
    preds <- factor(preds, levels = c("Basal", "TA", "TD", "unknown"))
    cells[[s]]$cell_type_pred <- preds
}

cells$ctrl6d <- subset(cells$d6, treatment_class == "Control")
cells$rocki6d <- subset(cells$d6, treatment_class == "ROCKi")
cells$ctrl12d <- subset(cells$d12, treatment_class == "Control")
cells$rocki12d <- subset(cells$d12, treatment_class == "ROCKi")

create_cell_type_umap <- function(cells) {
    # Force factors back into correct order as this inexplicably disappears
    # when passed to a function
    cells$cell_type_pred <- factor(cells$cell_type_pred,
                                   levels = c("Basal", "TA", "TD", "unknown"))
    p <- DimPlot(cells, group.by = "cell_type_pred") +
         scale_color_manual(labels = c("Holoclone-forming",
                                       "Mero- or Paraclone-forming",
                                       "Differentiated",
                                       "Unclassified keratinocytes"),
                            values = c("#80C9EA", "#DD6E79", "#43863E",
                                       "#BBBBBB")) +
         coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

for (i in list(c("ctrl6d", "fig_6a.tiff"), c("rocki6d", "fig_6c.tiff"),
               c("ctrl12d", "fig_6b.tiff"), c("rocki12d", "fig_6d.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    cell_type_umap <- create_cell_type_umap(cells[[s]])
    print(cell_type_umap)
    if (save_images) {
        save_tiff(cell_type_umap, out_name)
    }
}

cell_type_prop_graph <- function(df, names) {
    p <- ggplot(df, aes(x = Class, fill = CellType)) +
         geom_bar(position = "fill") + ylab("Proportion of cell type") +
         labs(fill = "Cell type") +
         scale_fill_manual(labels = c("Holoclone-forming",
                                      "Mero- or Paraclone-forming",
                                      "Differentiated",
                                      "Unclassified keratinocytes"),
                           values = c("#80C9EA", "#DD6E79", "#43863E",
                                      "#BBBBBB")) +
         scale_x_discrete(labels = names) +
         chart_margin_theme + axis_text_theme +
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               axis.title.x = element_blank())

    return(p)
}

for (i in list(c("d6", 6, "fig_6e.tiff"), c("d12", 12, "fig_6f.tiff"))) {
    s <- i[1]
    num_days <- i[2]
    out_name <- i[3]

    s1 <- paste0("Control-", num_days, "D")
    s2 <- paste0("ROCKi-", num_days, "D")
    cell_type_df <- data.frame(CellType = cells[[s]]$cell_type_pred,
                               Class = cells[[s]]$sample_class)
    cell_type_df$Class <- factor(cell_type_df$Class,
                                 levels = c(s1, s2))
    prop_graph <- cell_type_prop_graph(cell_type_df,
                                       c(paste0(s1, " "), paste0(s2, " ")))
    prop_graph_no_legend <- prop_graph + theme(legend.position = "none")
    print(prop_graph_no_legend)

    if (save_images) {
        save_tiff(prop_graph_no_legend, out_name, width = 6, height = 6)
    }
}
prop_table_cell_type <- function(x) {
    return(prop.table(table(factor(x, levels = c("Basal", "TA", "TD",
                                                 "unknown")))))
}

calc_cell_type_props <- function(cell_types, treatment_classes) {
    # Calculate the proportion of each cell type for both treatment classes.
    df <- data.frame(type = cell_types, treatment = treatment_classes)

    aggred_means <- aggregate(df$type, list(df$treatment),
                              prop_table_cell_type)
    for (i in 1:2) {
        aggred_means[, 2][i, ] <- aggred_means[, 2][i, ] /
                                  sum(aggred_means[, 2][i, ][1:3])
    }

    proportions <- c(aggred_means[1, ][[2]][1], aggred_means[1, ][[2]][2],
                     aggred_means[1, ][[2]][3], aggred_means[2, ][[2]][1],
                     aggred_means[2, ][[2]][2], aggred_means[2, ][[2]][3])

    return(sum(abs(proportions[1:2] - proportions[4:5])))
}

# Set random seed for reproducibility
set.seed(1234)

for (i in list(c("d6", 6), c("d12", 12))) {
    s <- i[1]
    num_days <- i[2]

    obs_total_prop_diff <- calc_cell_type_props(cells[[s]]$cell_type_pred,
                                                cells[[s]]$treatment_class)
    prop_diffs <- c()
    for (r in 1:10000) {
        # Generate resample indices
        c_idx <- sample(which(cells[[s]]$treatment_class == "Control"),
                          length(which(cells[[s]]$treatment_class == "Control")),
                          replace = TRUE)
        r_idx <- sample(which(cells[[s]]$treatment_class == "ROCKi"),
                          length(which(cells[[s]]$treatment_class == "ROCKi")) - 1,
                          replace = TRUE)
        idx <- c(c_idx, r_idx)

        # Find cell type and treatment class labels for this resample
        cell_types <- cells[[s]]$cell_type_pred[idx]
        treat_classes <- cells[[s]]$treatment_class[idx]

        total_prop_diff <- calc_cell_type_props(cell_types, treat_classes)
        prop_diffs <- c(prop_diffs, total_prop_diff)
    }

    alpha_vals <- c(0.975, 0.025)
    residuals <- prop_diffs - obs_total_prop_diff
    print(paste0("95% confidence interval for total absolute difference in",
                 "cell type proportions at ", num_days, "D:"))
    print(as.numeric(obs_total_prop_diff - quantile(residuals, alpha_vals)))
}
create_dotplot <- function(cells, signatures) {
    dp <- DotPlot(cells, group.by = "sample_class_space",
                  features = list("Holoclone" = signatures$Basal,
                                  "Mero-/Paraclone" = signatures$TA,
                                  "Differentiated" = signatures$TD)) +
          axis_text_theme +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                axis.text = element_text(size = 10),
                axis.title = element_blank(),
                legend.title = element_text(size = 10, face = "bold"),
                legend.text = element_text(size = 10, face = "bold"))

    return(dp)
}

cells$all$sample_class_space <- cells$all$sample_class
levels(cells$all$sample_class_space) <- c("Control-6D ", "ROCKi-6D ",
                                          "Control-6D 6D ", "ROCKi-6D 6D ")

cells$holo <- subset(cells$all, subset = cell_type_pred == "Basal")
cells$ta <- subset(cells$all, subset = cell_type_pred == "TA")
cells$td <- subset(cells$all, subset = cell_type_pred == "TD")

for (i in list(c("holo", "fig_6g.tiff"), c("ta", "fig_s2a.tiff"),
               c("td", "fig_s2b.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    dp <- create_dotplot(cells[[s]], signatures)
    print(dp)

    if (save_images) {
        save_tiff(dp, out_name, width = 20, height = 8)
    }
}
cells$ctrl6d <- subset(cells$d6, treatment_class == "Control")
cells$rocki6d <- subset(cells$d6, treatment_class == "ROCKi")
cells$ctrl12d <- subset(cells$d12, treatment_class == "Control")
cells$rocki12d <- subset(cells$d12, treatment_class == "ROCKi")

for (s in c("ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    cells[[s]] <- RunUMAP(cells[[s]], dims = 1:50, reduction = "harmony")
}

prepare_cells_for_slingshot <- function(cells) {
    cells <- FindNeighbors(cells, reduction = "harmony", dims = 1:50)
    cells <- FindClusters(cells, resolution = 0.1)
    if (all(cells$sample_class == "Control-6D")) {
        cells <- RenameIdents(object = cells,
                              "0" = "TA", "1" = "Basal", "2" = "TD")
    } else {
        cells <- RenameIdents(object = cells,
                              "0" = "TA", "1" = "TD", "2" = "Basal")
    }

    return(cells)
}

for (s in c("d6", "d12", "ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    cells[[s]] <- prepare_cells_for_slingshot(cells[[s]])
}
# Change UMAP values to make figures all similarly oriented
for (s in c("ctrl12d")) {
    flipped <- (-1) * Embeddings(object = cells[[s]], reduction = "umap")[, 1]
    cells[[s]]@reductions$umap@cell.embeddings[, 1] <- flipped
}
for (s in c("ctrl6d", "rocki6d")) {
    flipped <- (-1) * Embeddings(object = cells[[s]], reduction = "umap")[, 2]
    cells[[s]]@reductions$umap@cell.embeddings[, 2] <- flipped
}

colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)

for (s in c("d6", "d12", "ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    cells[[paste0(s, "_sce")]] <- as.SingleCellExperiment(cells[[s]])
}

sshots <- list()

for (s in c("d6", "d12")) {
    sshots[[s]] <- slingshot(cells[[paste0(s, "_sce")]], Idents(cells[[s]]),
                             reducedDim = "HARMONY", start.clus = "Basal")
}

for (s in c("ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    sshots[[s]] <- slingshot(cells[[paste0(s, "_sce")]], Idents(cells[[s]]),
                             reducedDim = "UMAP", start.clus = "Basal")
}

umap_curves <- list()

for (s in c("d6", "d12")) {
    umap_embed <- reducedDims(cells[[paste0(s, "_sce")]])$UMAP
    umap_curves[[s]] <- slingCurves(embedCurves(sshots[[s]], umap_embed))[[1]]
}

for (s in c("ctrl6d", "rocki6d", "ctrl12d", "rocki12d")) {
    umap_curves[[s]] <- slingCurves(sshots[[s]])[[1]]
}

for (num_days in c(6, 12)) {
    s <- paste0("d", num_days)
    ctrl_idxs <- which(cells[[s]]$treatment_class == "Control")
    rocki_idxs <- which(cells[[s]]$treatment_class == "ROCKi")
    sshots[[paste0("ctrl_", num_days, "d_together")]] <- sshots[[s]][, ctrl_idxs]
    sshots[[paste0("rocki_", num_days, "d_together")]] <- sshots[[s]][, rocki_idxs]
}

xmin_sep <- min(reducedDim(cells$ctrl6d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 1],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 1])
xmax_sep <- max(reducedDim(cells$ctrl6d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 1],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 1],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 1])
ymin_sep <- min(reducedDim(cells$ctrl6d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 2],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 2])
ymax_sep <- max(reducedDim(cells$ctrl6d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki6d_sce, "UMAP")[, 2],
                reducedDim(cells$ctrl12d_sce, "UMAP")[, 2],
                reducedDim(cells$rocki12d_sce, "UMAP")[, 2])

create_traj_umap <- function(sshot_to_plot, umap_curve, cell_types) {
    points_df <- data.frame(u1 = reducedDims(sshot_to_plot)$UMAP[, 1],
                            u2 = reducedDims(sshot_to_plot)$UMAP[, 2],
                            cell_types = cell_types)
    points_df$cell_types <- factor(points_df$cell_types,
                                   levels = c("Basal", "TA", "TD", "unknown"))
    curve_df <- data.frame(u1 = umap_curve$s[umap_curve$ord, ][, 1],
                           u2 = umap_curve$s[umap_curve$ord, ][, 2])

    arrowhead_length <- unit(0.2, "cm")
    arrow_width <- 1

    ptimes <- sshot_to_plot$slingPseudotime_1
    basal_ptimes <- ptimes[which(cell_types == "Basal")]
    ta_ptimes <- ptimes[which(cell_types == "TA")]
    td_ptimes <- ptimes[which(cell_types == "TD")]

    percentile <- 0.97
    basal_ptimes_dist <- abs(quantile(basal_ptimes, percentile) - basal_ptimes)
    basal_percentile_ptime_idx <- which.min(basal_ptimes_dist)
    ta_ptimes_dist <- abs(quantile(ta_ptimes, percentile) - ta_ptimes)
    ta_percentile_ptime_idx <- which.min(ta_ptimes_dist)
    td_ptimes_dist <- abs(quantile(td_ptimes, percentile) - td_ptimes)
    td_percentile_ptime_idx <- which.min(td_ptimes_dist)

    basal_percentile_ptime <- basal_ptimes[basal_percentile_ptime_idx]
    ta_percentile_ptime <- ta_ptimes[ta_percentile_ptime_idx]
    td_percentile_ptime <- td_ptimes[td_percentile_ptime_idx]

    end_idx <- round(length(sshot_to_plot$slingPseudotime_1) * percentile)
    basal_pct_overall_idx <- which(sort(ptimes) == basal_percentile_ptime)
    ta_pct_overall_idx <- which(sort(ptimes) == ta_percentile_ptime)

    p <- ggplot() +
         geom_point(data = points_df, aes(x = u1, y = u2, color = cell_types),
                    size = -0.1) +
         scale_color_manual(labels = c("Holoclone-forming",
                           "Mero- or Paraclone-forming",
                           "Differentiated", "Unknown"),
                           values = c("#80C9EA", "#DD6E79",
                                      "#43863E", "#BBBBBB")) +
         geom_path(data = curve_df[seq_len(end_idx), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[1:2, ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "open", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[seq_len(basal_pct_overall_idx), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[seq_len(ta_pct_overall_idx), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         coord_fixed() + umap_theme

     return(p)
}

for (i in list(c("ctrl6d", "fig_7a.tiff"), c("rocki6d", "fig_7b.tiff"),
               c("ctrl12d", "fig_7c.tiff"), c("rocki12d", "fig_7d.tiff"))) {
    s <- i[1]
    out_name <- i[2]
    traj_umap <- create_traj_umap(sshots[[s]], umap_curves[[s]],
                                   Idents(cells[[s]])) +
                  theme(legend.position = "none")
    print(traj_umap)

    if (save_images) {
        save_tiff(traj_umap, out_name, width = 10, height = 8)
    }
}
for (s in c("d6", "d12")) {
    cells[[s]]$cell_type <- Idents(cells[[s]])
}

# mov_avg function copied from stackoverflow.com/a/4862334/6914552
mov_avg <- function(x, n = 151) {
    stats::filter(x, rep(1 / n, n), sides = 2)
}

create_traj_df <- function(ptimes, name, cell_types) {
    quantile_points <- seq(0.001, 1, 0.001)
    quantiled_data <- as.numeric(quantile(ptimes, quantile_points))
    sorted_cell_types <- cell_types[order(ptimes)]
    num_points <- length(ptimes)
    cell_indices <- round(num_points * quantile_points)
    cell_types_subset <- sorted_cell_types[cell_indices]

    return(data.frame(pct = quantile_points,
                      ptime = quantiled_data,
                      cell_type = factor(cell_types_subset,
                                         levels = c("Basal", "TA", "TD",
                                                    "unknown")),
                      diff = c(NA, mov_avg(diff(quantiled_data))),
                      sample_class = rep(name, length(quantile_points))))
}

traj_plot_margin_theme <- theme(plot.margin = margin(0.1, 1.2, 0.1, 1.2, "cm"))
diff_comparison_plots_theme <- theme(legend.title = element_blank())

point_size <- 1
create_traj_plot <- function(traj_df) {
    shade_alpha <- 0.2

    basal_upper_lim <- quantile(traj_df[traj_df$cell_type == "Basal", ]$ptime,
                                0.75)
    ta_upper_lim <- quantile(traj_df[traj_df$cell_type == "TA", ]$ptime, 0.95)

    basal_col <- "#80C9EA"
    ta_col <- "#DD6E79"
    td_col <- "#43863E"

    traj_plot <- ggplot(traj_df,
                        aes(x = pct, y = ptime, group = sample_class)) +
                 annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
                          ymax = basal_upper_lim, alpha = shade_alpha,
                          fill = basal_col) +
                 annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = basal_upper_lim, ymax = ta_upper_lim,
                          alpha = shade_alpha, fill = ta_col) +
                 annotate("rect", xmin = -Inf, xmax = Inf,
                          ymin = ta_upper_lim, ymax = Inf,
                          alpha = shade_alpha, fill = td_col) +
                 geom_point(aes(color = cell_type), size = point_size) +
                 scale_color_manual(values = c(basal_col, ta_col, td_col)) +
                 new_scale_color() +
                 geom_line(aes(linetype = sample_class), size = 0.4) +
                 scale_linetype_manual(values = c("solid", "dashed")) +
                 scale_x_continuous(labels = scales::percent,
                                    limits = c(0, 1)) +
                 ylim(c(0, max(traj_df$ptime))) +
                 xlab("Percentage through differentation process") +
                 ylab("Pseudotime") + traj_plot_margin_theme +
                 diff_comparison_plots_theme + axis_text_theme

    return(traj_plot)
}

create_speed_plot <- function(traj_df) {
    speed_plot <- ggplot(traj_df,
                         aes(x = pct, y = diff, group = sample_class)) +
                  geom_point(aes(color = cell_type), size = point_size) +
                  scale_color_manual(values = c("#80C9EA", "#DD6E79",
                                                "#43863E")) +
                  new_scale_color() +
                  geom_line(aes(linetype = sample_class), size = 0.4) +
                  scale_linetype_manual(values = c("solid", "dashed")) +
                  scale_x_continuous(labels = scales::percent,
                                     limits = c(0, 1)) +
                  xlab("Percentage through differentation process") +
                  ylab("Differentiation speed") + traj_plot_margin_theme +
                  diff_comparison_plots_theme + axis_text_theme

    return(speed_plot)
}

obs_results_ptime <- list()
obs_results_speed <- list()
for (i in list(c(6, "fig_7e.tiff", "fig_s3a.tiff"),
               c(12, "fig_7f.tiff", "fig_s3b.tiff"))) {
    num_days <- i[1]
    traj_graph_out_name <- i[2]
    speed_graph_out_name <- i[3]
    s <- paste0("d", num_days)
    ctrl_idxs <- which(cells[[s]]$treatment_class == "Control")
    rocki_idxs <- which(cells[[s]]$treatment_class == "ROCKi")
    ctrl_cell_types <- cells[[paste0("d", num_days)]][, ctrl_idxs]$cell_type
    rocki_cell_types <- cells[[paste0("d", num_days)]][, rocki_idxs]$cell_type
    cp <- sshots[[paste0("ctrl_", num_days, "d_together")]]$slingPseudotime_1
    rp <- sshots[[paste0("rocki_", num_days, "d_together")]]$slingPseudotime_1

    traj_df_ctrl <- create_traj_df(cp, paste0("Control-", num_days, "D"),
                                   ctrl_cell_types)
    traj_df_rocki <- create_traj_df(rp, paste0("ROCKi-", num_days, "D"),
                                    rocki_cell_types)
    traj_df <- rbind(traj_df_ctrl, traj_df_rocki)

    ptime_abs_diffs <- abs(traj_df_ctrl$ptime - traj_df_rocki$ptime)
    obs_results_ptime[[paste0("d", num_days)]] <- mean(ptime_abs_diffs)

    speed_abs_diffs <- abs(traj_df_ctrl$diff - traj_df_rocki$diff)
    obs_results_speed[[paste0("d", num_days)]] <- mean(speed_abs_diffs,
                                                      na.rm = TRUE)

    traj_plot <- create_traj_plot(traj_df) + theme(legend.position = "none")
    print(traj_plot)
    speed_plot <- create_speed_plot(traj_df) + theme(legend.position = "none")
    print(speed_plot)

    if (save_images) {
        save_tiff(traj_plot, traj_graph_out_name, width = 10, height = 8)
        save_tiff(speed_plot, speed_graph_out_name, width = 14, height = 8)
    }
}

legend_plot <- ggplot(traj_df, aes(x = pct, y = ptime, group = sample_class)) +
               geom_point(aes(color = cell_type), size = 3) +
               scale_color_manual(values = c("#80C9EA", "#DD6E79", "#43863E"),
                                  labels = c("Holoclone-forming",
                                             "Mero- or Paraclone-forming",
                                             "Differentiated")) +
               theme(legend.title = element_blank(),
                     legend.text = element_text(face = "bold"))
(dot_colour_legend <- as_ggplot(get_legend(legend_plot, "right")))

legend_plot <- ggplot(traj_df, aes(x = pct, y = ptime, group = sample_class)) +
               geom_line(aes(linetype = sample_class), size = 0.4) +
               scale_linetype_manual(values = c("solid", "dashed"),
                                     labels = c("Control", "ROCKi-treated")) +
               theme(legend.title = element_blank(),
                     legend.text = element_text(face = "bold"))
(linetype_legend <- as_ggplot(get_legend(legend_plot, "right")))

(legends <- ggarrange(linetype_legend, dot_colour_legend, ncol = 2))
if (save_images) {
    save_tiff(legends, "fig_7_legends.tiff", width = 10, height = 8)
}
# Generate bootstrapped data as in "experiment_code_trajectories.R" and save
# the outputs as "trajectory_comparison_results.txt"

# bs_out = "bootstrap out"
bs_out <- read.table("trajectory_comparison_results.txt")
bs_out <- bs_out[, c(3, 4, 5, 6)]
colnames(bs_out) <- c("TrajDiff6D", "TrajDiff12D", "SpeedDiff6D",
                      "SpeedDiff12D")

for (num_days in c(6, 12)) {
    traj_diffs <- bs_out[paste0("TrajDiff", num_days, "D")][, 1]
    obs_result_ptime <- obs_results_ptime[[paste0("d", num_days)]]
    diffs_quants_traj <- quantile((traj_diffs - obs_result_ptime), alpha_vals)
    residuals_traj <- as.numeric(obs_result_ptime - diffs_quants_traj)
    print(paste0("95% confidence interval for difference in pseudotimes at ",
                 num_days, "D: ", paste(residuals_traj, collapse = " ")))

    speed_diffs <- bs_out[paste0("SpeedDiff", num_days, "D")][, 1]
    obs_result_speed <- obs_results_speed[[paste0("d", num_days)]]
    diffs_quants_speed <- quantile((speed_diffs - obs_result_speed), alpha_vals)
    residuals_speed <- as.numeric(obs_result_speed - diffs_quants_speed)
    print(paste0("95% confidence interval for difference in speed at ",
                 num_days, "D: ", paste(residuals_speed, collapse = " ")))
}
cells$all_with_outliers <- cells$all_with_outliers %>%
                           SCTransform() %>%
                           CellCycleScoring(s.features = s_genes,
                                            g2m.features = g2m_genes)
cells$all_with_outliers$CC.difference <- cells$all_with_outliers$S.Score -
                                         cells$all_with_outliers$G2M.Score
cells$all_with_outliers <- cells$all_with_outliers %>%
                           SCTransform(vars.to.regress = c("CC.difference")) %>%
                           RunPCA(features = VariableFeatures(cells$all_with_outliers,
                                                              assay = "SCT"),
                                  npcs = 50)
donor_idents <- sapply(cells$all_with_outliers$orig.ident,
                       function(x) {
                           if (x %in% c("1", "2", "3", "4")) {
                               "1"
                           } else if (x %in% c("5", "6", "7", "8")) {
                               "2"
                           } else {
                               "3"
                           }
                       })
cells$all_with_outliers$donor_idents <- donor_idents
cells$all_with_outliers <- cells$all_with_outliers %>%
                           RunHarmony("donor_idents", assay.use = "SCT") %>%
                           RunUMAP(dims = 1:50, reduction = "harmony")
(p <- FeaturePlot(cells$all_with_outliers, "VIM"))
if (save_images) {
    save_tiff(p, "fig_s4a.tiff", width = 10, height = 8)
}
cells$all <- RunUMAP(cells$all, reduction = "harmony", dims = 1:50)
p <- FeaturePlot(cells$all, features = c("KRT14", "VIM"), blend = TRUE,
                 split.by = "sample_class")
blended_fp <- list()
for (i in 1:4) {
    q <- ggarrange(p[[((i - 1) * 4) + 1]] + umap_theme,
                   p[[((i - 1) * 4) + 2]] + umap_theme,
                   p[[((i - 1) * 4) + 3]] + umap_theme,
                   p[[((i - 1) * 4) + 4]] + theme(plot.title = element_blank()),
                   nrow = 1, widths = c(2, 2, 2, 3))
    blended_fp[[i]] <- q
    print(q)
}
save_tiff(ggarrange(plotlist = blended_fp, ncol = 1), "fig_s4b.tiff",
          width = 15, height = 15)
$(document).ready(function() {
  $('.footnotes ol').appendTo('#endnotes');
  $('.footnotes').remove();
});
installed.packages()[names(sessionInfo()$otherPkgs), "Version"]
```

<!-- Javascript to put footnotes in correct location (https://stackoverflow.com/a/56866117) -->
<script type="text/javascript">
$(document).ready(function() {
  $('.footnotes ol').appendTo('#endnotes');
  $('.footnotes').remove();
});
</script>

## Package versions

<!-- https://stackoverflow.com/a/36777618 -->

```r
installed.packages()[names(sessionInfo()$otherPkgs), "Version"]
```

```
##                DAseq SingleCellExperiment SummarizedExperiment 
##              "1.0.0"             "1.12.0"             "1.20.0" 
##              Biobase        GenomicRanges         GenomeInfoDb 
##             "2.50.0"             "1.42.0"             "1.26.2" 
##              IRanges            S4Vectors         BiocGenerics 
##             "2.24.1"             "0.28.1"             "0.36.0" 
##       MatrixGenerics          matrixStats            slingshot 
##              "1.2.1"             "0.61.0"              "1.8.0" 
##            princurve                SCINA               gplots 
##              "2.1.6"              "1.2.0"              "3.1.1" 
##                 MASS                 kBET              harmony 
##           "7.3-53.1"             "0.99.6"              "0.1.0" 
##                 Rcpp          sctransform         SeuratObject 
##              "1.0.7"              "0.3.2"              "4.0.4" 
##               Seurat           ggnewscale         RColorBrewer 
##              "4.0.5"              "0.4.5"              "1.1-2" 
##            patchwork               ggpubr                  png 
##         "1.1.0.9000"              "0.4.0"              "0.1-7" 
##              ggrepel              ggplot2                tidyr 
##              "0.9.1"              "3.3.5"              "1.1.4" 
##                dplyr            rmarkdown 
##              "1.0.7"                "2.9"
```

<!-- Stop large blank section at end of report (https://stackoverflow.com/a/57381047) -->
<div class="tocify-extend-page" data-unique="tocify-extend-page" style="height: 0;"></div>