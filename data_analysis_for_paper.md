---
title: "Impact of ROCKi Treatment on Keratinocytes is Reversible: Single-Cell RNAseq Data Analysis"
author: "George T. Hall"
date: "Compiled on 30 May 2022"
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
skin_data <- CreateSeuratObject(counts = skin_data,
                                project = "3DonorEffectOfRocki",
                                min.cells = 3, min.features = 200,
                                names.delim = "-", names.field = 2)
skin_data <- RenameIdents(object = skin_data, "1" = "1Control-6D ",
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
skin_data[["percent.mt"]] <- PercentageFeatureSet(skin_data, pattern = "^MT-")

VlnPlot_limits <- function(data, metric, ylab, limits) {
    df <- data.frame(x=Idents(data), y = data@meta.data[,metric])
    plot <- ggplot() + geom_violin(data = df, aes(x=x, y=y, fill=x)) +
            geom_errorbar(data = limits, aes(x = x, ymin = lower, ymax = upper),
                          color = "red") +
            labs(y=ylab) +
            scale_y_continuous(labels = scales::comma) +
            theme(axis.text.x = element_text(angle = 45, hjust=1),
                  legend.position = "none", plot.title = element_blank(),
                  axis.title.x = element_blank()) + axis_text_theme

    return(plot)
}

ncount_df <- data.frame(x = levels(Idents(skin_data)),
                        lower = c(8000, 15000, 8000, 8000,
                                  8000, 8000, 5000, 8000,
                                  8000, 8000, 5000, 5000),
                        upper = c(50000, 62500, 50000, 50000,
                                  50000, 50000, 50000, 50000,
                                  50000, 50000, 40000, 40000))

(ncount_vln <- VlnPlot_limits(skin_data, "nCount_RNA",
                               "Number of RNA\nmolecules per cell",
                               ncount_df) + ylim(c(0, 100000)))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(ncount_vln, "fig_s4a.tiff")
}

nfeature_df <- data.frame(x = levels(Idents(skin_data)),
                     lower = c(2500, 3000, 2500, 2500,
                               2500, 2500, 2000, 2500,
                               2500, 2500, 1500, 1500),
                     upper = c(7500, 8500, 7500, 7500,
                               7500, 7500, 7500, 7500,
                               7500, 7500, 6000, 6000))

(nfeature_vln <- VlnPlot_limits(skin_data, "nFeature_RNA",
                               "Number of genes\ndetected per cell",
                               nfeature_df))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-3-2.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(nfeature_vln, "fig_s4b.tiff")
}

pct_mito_df <- data.frame(x = levels(Idents(skin_data)),
                        lower = c(rep(2, 8), rep(3, 4)),
                        upper = c(rep(10, 8), rep(12, 4)))

(pct_mito_vln <- VlnPlot_limits(skin_data, "percent.mt",
                               "Percentage of mitochrondrial\nDNA detected per cell",
                               pct_mito_df) + ylim(c(0, 20)))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-3-3.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(pct_mito_vln, "fig_s4c.tiff")
}

pass_qc <- c()
for (j in 1:length(skin_data@meta.data$orig.ident)) {
    i <- as.integer(as.character(skin_data@meta.data$orig.ident[j]))
    pass_qc <- c(pass_qc, all(skin_data@meta.data$nFeature_RNA[j] >= nfeature_df$lower[i],
                              skin_data@meta.data$nFeature_RNA[j] <= nfeature_df$upper[i],
                              skin_data@meta.data$nCount_RNA[j] >= ncount_df$lower[i],
                              skin_data@meta.data$nCount_RNA[j] <= ncount_df$upper[i],
                              skin_data@meta.data$percent.mt[j] >= pct_mito_df$lower[i],
                              skin_data@meta.data$percent.mt[j] <= pct_mito_df$upper[i]))
}
skin_data <- skin_data[, pass_qc]
```

## Filter outliers

We removed non-keratinocytes by excluding cells expressing at least one count
of markers of melanocytes (MLANA, PMEL, MITF), mesenchymal cells (MTRNR2L6,
MTRNR2L10, MTRNR2L7, MTRNR2L1, PRKAR2B, NR2F1), or fibroblasts (ACTG2, DLK1).


```r
# Remove melanocytes
skin_data <- subset(skin_data, subset = MLANA > 1 | PMEL > 1 | MITF > 1,
                    invert = T)
# Remove mesenchymal cells
skin_data <- subset(skin_data,
                    subset = MTRNR2L6 > 1 | MTRNR2L10 > 1 | MTRNR2L7 > 1 |
                             MTRNR2L1 > 1 | PRKAR2B > 1 | NR2F1 > 1,
                    invert = T)
# Remove fibroblasts
skin_data <- subset(skin_data, subset = ACTG2 > 1 | DLK1 > 1, invert = T)
```

## Normalise the data

We normalised the gene expressions using sctransform. 


```r
skin_data <- SCTransform(skin_data)
skin_data.outliers_removed <- skin_data
```

## Regress out cell cycle

We then reduced the effect of cell cycle by assigning scores to each cell
representing the probabilities of it being in different stages of the cell
cycle (using Seurat???s CellCycleScoring function) and then regressing out the
difference between these scores using sctransform. The authors of Seurat
recommend this approach when dealing with cells undergoing differentiation
(https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html -
"Alternative Workflow").


```r
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
skin_data <- CellCycleScoring(skin_data, s.features = s.genes,
                              g2m.features = g2m.genes)
skin_data@meta.data$CC.difference <- skin_data@meta.data$S.Score - skin_data@meta.data$G2M.Score
skin_data <- SCTransform(skin_data, vars.to.regress = c("CC.difference"))
skin_data.cell_cycle_regressed <- skin_data
```

## Run PCA and remove donor effect

We carried out dimensionality reduction using principal components analysis
(PCA), using the 3000 most highly variable genes as features. We used the 50
most significant components since there was no clear point at which the
generated components became less significant. We reduced the inter-donor
variation with Harmony.


```r
# Run PCA
skin_data <- RunPCA(skin_data, features = skin_data@assays$SCT@var.features,
                    npcs = 50)
(elbow_plot <- ElbowPlot(skin_data, ndims = 50) + axis_text_theme)
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(elbow_plot, "fig_s5.tiff")
}

# Remove inter-donor differences using Harmony
donor_idents <- sapply(skin_data@meta.data$orig.ident,
                       function(x) if (x %in% c("1","2","3","4")) {"1"}
                                   else if (x %in% c("5","6","7","8")) {"2"}
                                   else {"3"})
skin_data@meta.data["donor_idents"] <- donor_idents

skin_data <- skin_data %>% RunHarmony("donor_idents", assay.use = "SCT")
skin_data.harmonised <- skin_data

sample_class <- sapply(skin_data@meta.data$orig.ident,
                       function(x) if (x %in% c("1","5","9")) {"Control-6D"}
                                   else if (x %in% c("2","6","10")) {"ROCKi-6D"}
                                   else if (x %in% c("3","7","11")) {"Control-12D"}
                                   else {"ROCKi-12D"})
skin_data@meta.data$sample_class <- sample_class
skin_data@meta.data$sample_class <- factor(skin_data@meta.data$sample_class,
                                           levels = c("Control-6D", "ROCKi-6D",
                                                      "Control-12D",
                                                      "ROCKi-12D"))

treatment_class <- sapply(skin_data@meta.data$sample_class,
                       function(x) if (x %in% c("Control-6D","Control-12D")) {"Control"}
                                   else {"ROCKi"})
skin_data@meta.data$treatment_class <- treatment_class
```

## Run UMAP

We visualised the cells using UMAP.


```r
d6_cells <- subset(skin_data,
                   subset = sample_class %in% c("Control-6D", "ROCKi-6D"))
d12_cells <- subset(skin_data,
                    subset = sample_class %in% c("Control-12D", "ROCKi-12D"))

d6_cells <- RunUMAP(d6_cells, dims = 1:50, reduction = "harmony")
d12_cells <- RunUMAP(d12_cells, dims = 1:50, reduction = "harmony")

umap_theme <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                    axis.text.y=element_blank(), axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    plot.title = element_blank())
chart_margin_theme <- theme(plot.margin = margin(0, 1.2, 0, 1.2, "cm"))

labels <- c("Control", "ROCKi-treated")
umap_legend_plot <- ggplot(data.frame(class = as.character(1:length(labels)),
                                      x = 1:length(labels)),
                           aes(x = x, y = x, col = class)) +
                    geom_point() +
                    scale_color_manual(labels = labels,
                                       values = c(blue_color, red_color)) +
                    theme(legend.title = element_blank(),
                          legend.key = element_rect(fill = NA, color = NA),
                          legend.margin = margin(c(0,0,0,0)),
                          legend.text = element_text(face = "bold"),
                          legend.background = element_rect(fill="#00000000",
                                                           size=0.5,
                                                           linetype="solid")) +
                    guides(colour = guide_legend(override.aes = list(size=2)))
umap_legend <- as_ggplot(get_legend(umap_legend_plot, "left"))
if(save_images) {
    save_tiff(umap_legend, "umap_legend.tiff")
}

plot_treatment_class_umap <- function(cells) {
    p <- DimPlot(cells, group.by = "treatment_class",
                 cols = c(blue_color, red_color),
                 pt.size = -0.1) +
         coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

(umap_d6 <- plot_treatment_class_umap(d6_cells))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(umap_d6, "fig_5a.tiff")
}

(umap_d12 <- plot_treatment_class_umap(d12_cells))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-8-2.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(umap_d12, "fig_5b.tiff")
}
```

Visualising the 6D and 12D cell populations using UMAP provides an intuitive
indication that the difference between treated and control cells is greater at
6D than it is at 12D. There is substantially less overlap between the treated
and untreated cells at 6D than at 12D. We now measure the amount of overlap
using DAseq.

## Differential abundance testing using DAseq


```r
daseq_out_6d <- getDAcells(d6_cells@reductions$harmony@cell.embeddings,
                           cell.labels = as.character(d6_cells@meta.data$sample_class),
                           labels.1 = c("Control-6D"),
                           labels.2 = c("ROCKi-6D"),
                           plot.embedding = d6_cells@reductions$umap@cell.embeddings)

daseq_out_12d <- getDAcells(d12_cells@reductions$harmony@cell.embeddings,
                            cell.labels = as.character(d12_cells@meta.data$sample_class),
                            labels.1 = c("Control-12D"),
                            labels.2 = c("ROCKi-12D"),
                            plot.embedding = d12_cells@reductions$umap@cell.embeddings)

daseq_legend <- as_ggplot(get_legend(daseq_out_12d$pred.plot +
                                     umap_theme +
                                     theme(legend.position = "top",
                                           legend.text = element_blank())))

if(save_images) {
    save_tiff(daseq_legend, "daseq_legend.tiff")
}

create_daseq_umap <- function(daseq_out) {
    daseq_umap <- daseq_out$pred.plot + umap_theme + theme(legend.position = "none")

    return(daseq_umap)
}

(daseq_umap_6d <- create_daseq_umap(daseq_out_6d))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(daseq_umap_6d, "fig_5c.tiff")
}

(daseq_umap_12d <- create_daseq_umap(daseq_out_12d))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-9-2.png" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(daseq_umap_12d, "fig_5d.tiff")
}

proportion_da_6d <- length(c(daseq_out_6d$da.down, daseq_out_6d$da.up)) / length(colnames(d6_cells))
proportion_da_12d <- length(c(daseq_out_12d$da.down, daseq_out_12d$da.up)) / length(colnames(d12_cells))
```


```r
print(paste("Proportion of cells in differentially abundant regions at 6D:",
            proportion_da_6d))
```

```
## [1] "Proportion of cells in differentially abundant regions at 6D: 0.902621278614636"
```

```r
print(paste("Proportion of cells in differentially abundant regions at 12D:",
            proportion_da_12d))
```

```
## [1] "Proportion of cells in differentially abundant regions at 12D: 0.371958285052144"
```


# Differential Expression Analysis

We carried out differential expression analysis of treated and control cells at
6D and 12D using the standard Wilcoxon rank sum test in Seurat. 


```r
Idents(skin_data) <- skin_data@meta.data$sample_class
marks_6d <- FindMarkers(skin_data, ident.1 = "ROCKi-6D",
                        ident.2 = "Control-6D")
marks_12d <- FindMarkers(skin_data, ident.1 = "ROCKi-12D",
                         ident.2 = "Control-12D")

not_sig_color <- "#4D4D4D"

VolPlot <- function(marks) {
    marks$p_val_adj_neg_log <- (-1) * log10(marks$p_val_adj)
    marks$p_val_adj_neg_log[is.infinite(marks$p_val_adj_neg_log)] <- 300
    sig_marks_pos <- subset(marks, avg_log2FC > 1)
    sig_marks_neg <- subset(marks, avg_log2FC < -1)

    # Assign classes to genes based on how DE they are. Used to set colours
    de_class <- ifelse(marks$avg_log2FC < -1 & marks$p_val_adj < 0.01, 'Downreg',
                  ifelse(marks$avg_log2FC > 1 & marks$p_val_adj < 0.01, 'Upreg',
                         ifelse(marks$p_val_adj < 0.01, "Low-LFC", "Not-Sig")))
    de_class <- factor(de_class, levels = c("Not-Sig", "Low-LFC",
                                            "Downreg", "Upreg"))

    x_sublabel_y <- -50 # Position the sub-labels of the x-axis correctly
    plot <- ggplot(marks, aes(x = avg_log2FC, y = p_val_adj_neg_log,
                              label = rownames(marks))) +
            geom_point(aes(color = de_class, fill = de_class), size = 2,
                       shape = 21, alpha = 0.5) +
            scale_color_manual(values = c("gray", not_sig_color, "red", "white")) +
            scale_fill_manual(values = c("gray", not_sig_color, "white", "red")) +
            geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
            geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
            xlim(-2.5, 2.5) + ylim(-log10(0.01), 350) +
            xlab(bquote(bold(Log["2"]~"fold change"))) +
            ylab(bquote(bold(-Log["10"] ~ P))) +
            theme(legend.position="none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            axis_text_theme +
            coord_cartesian(clip = "off")

    return (plot)
}

(vol_6d <- VolPlot(marks_6d))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-1.png" width="100%" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(vol_6d, "fig_5e.tiff", height = 5, width = 10)
}

labels <- c(bquote(bold(Log["2"]~"fold change ??? -1")),
            bquote(bold("|"*Log["2"]~"fold change| < 1")),
            bquote(bold(Log["2"]~"fold change ??? 1")))

vol_6d_legend_plot <- ggplot(data.frame(class = as.character(1:length(labels)),
                                        x = 1:length(labels)),
                             aes(x = x, y = x, color = class, fill = class)) +
                      geom_point(size = 3, shape = 21, alpha = 0.5) +
                      scale_fill_manual(values = c("white", not_sig_color, "red"),
                                        labels = labels) +
                      scale_color_manual(values = c("red", not_sig_color, "white"),
                                         labels = labels) +
                      theme(legend.title = element_blank(),
                            legend.key = element_rect(fill = NA, color = NA),
                            legend.text = element_text(face = "bold"),
                            legend.margin=margin(unit(c(0, 0, 0, 0), "cm"))) +
                      guides(colour = guide_legend(override.aes = list(size = 2)))

volplot_legend <- as_ggplot(get_legend(vol_6d_legend_plot, "right"))
if(save_images) {
    save_tiff(volplot_legend, "volplot_legend.tiff")
}

(vol_12d <- VolPlot(marks_12d))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-11-2.png" width="100%" style="display: block; margin: auto;" />

```r
if(save_images) {
    save_tiff(vol_12d, "fig_5f.tiff", height = 5, width = 10)
}
```

# Cell Type Proportions in Treated / Untreated

We used SCINA[^1] to estimate the cell type proportions in the treated/control
cells. We used cell type markers published by Enzo et al._et al_ [^2] to
classify cells as either basal, transient amplifying (TA), or terminal
differentiated (TD), since these are the classes for which Enzo et al. provide
distinct markers.


```r
(signatures <- preprocess.signatures("celltype_markers_basal_ta_td.csv"))

skin_data@meta.data$cell_type_pred <- SCINA(as.matrix(GetAssayData(skin_data)),
                                            signatures)$cell_labels
d6_cells@meta.data$cell_type_pred <- SCINA(as.matrix(GetAssayData(d6_cells)),
                                           signatures)$cell_labels
d12_cells@meta.data$cell_type_pred <- SCINA(as.matrix(GetAssayData(d12_cells)),
                                            signatures)$cell_labels

d6_cells@meta.data$cell_type_pred <- factor(d6_cells@meta.data$cell_type_pred,
                                             levels = c("Basal", "TA", "TD", "unknown"))
d12_cells@meta.data$cell_type_pred <- factor(d12_cells@meta.data$cell_type_pred,
                                             levels = c("Basal", "TA", "TD", "unknown"))

create_cell_type_umap <- function(cells) {
    # Force factors back into correct order as this inexplicably disappears
    # when passed to a function
    cells@meta.data$cell_type_pred <- factor(cells@meta.data$cell_type_pred,
                                             levels = c("Basal", "TA", "TD", "unknown"))
    p <- DimPlot(cells, group.by = "cell_type_pred") +
    scale_color_manual(labels=c("Holoclone-forming",
                               "Mero- or Paraclone-forming",
                               "Differentiated", "Unclassified keratinocytes"),
                      values = c("#80C9EA", "#DD6E79", "#43863E", "#BBBBBB")) +
    coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

ctrl6d_cells <- subset(d6_cells, treatment_class == "Control")
rocki6d_cells <- subset(d6_cells, treatment_class == "ROCKi")
ctrl12d_cells <- subset(d12_cells, treatment_class == "Control")
rocki12d_cells <- subset(d12_cells, treatment_class == "ROCKi")

cell_type_umap_ctrl6d <- create_cell_type_umap(ctrl6d_cells)
cell_type_umap_rocki6d <- create_cell_type_umap(rocki6d_cells)
cell_type_umap_ctrl12d <- create_cell_type_umap(ctrl12d_cells)
cell_type_umap_rocki12d <- create_cell_type_umap(rocki12d_cells)

if (save_images) {
    save_tiff(cell_type_umap_ctrl6d, "fig_6a.tiff")
    save_tiff(cell_type_umap_ctrl12d, "fig_6b.tiff")
    save_tiff(cell_type_umap_rocki6d, "fig_6c.tiff")
    save_tiff(cell_type_umap_rocki12d, "fig_6d.tiff")
}

cell_type_df_6d <- data.frame(CellType = d6_cells@meta.data$cell_type_pred,
                           Class = d6_cells@meta.data$sample_class)
cell_type_df_6d$Class <- factor(cell_type_df_6d$Class,
                                levels = c("Control-6D", "ROCKi-6D",
                                           "Control-12D", "ROCKi-12D"))

create_cell_type_proportion_graph <- function(df, names) {
    p <- ggplot(df, aes(x = Class, fill = CellType)) +
    geom_bar(position = "fill") + ylab("Proportion of cell type") +
    labs(fill = "Cell type") +
    scale_fill_manual(labels = c("Holoclone-forming",
                                 "Mero- or Paraclone-forming",
                                 "Differentiated",
                                 "Unclassified keratinocytes"),
                      values = c("#80C9EA", "#DD6E79", "#43863E", "#BBBBBB")) +
    scale_x_discrete(labels = names) +
    chart_margin_theme + axis_text_theme +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          axis.title.x = element_blank())

    return(p)
}

cell_type_proportion_graph_6d <- create_cell_type_proportion_graph(cell_type_df_6d, c("Control-6D ", "ROCKi-6D "))
(cell_type_proportion_graph_6d_no_legend <- cell_type_proportion_graph_6d + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(cell_type_proportion_graph_6d_no_legend, "fig_6e.tiff",
              width = 6, height = 6)
}

cell_type_df_12d <- data.frame(CellType = d12_cells@meta.data$cell_type_pred,
                           Class = d12_cells@meta.data$sample_class)
cell_type_df_12d$Class <- factor(cell_type_df_12d$Class,
                                 levels = c("Control-6D", "ROCKi-6D",
                                            "Control-12D", "ROCKi-12D"))
cell_type_proportion_graph_12d <- create_cell_type_proportion_graph(cell_type_df_12d, c("Control-6D 6D ", "ROCKi-6D 6D "))
(cell_type_proportion_graph_12d_no_legend <- cell_type_proportion_graph_12d + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-12-2.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(cell_type_proportion_graph_12d_no_legend, "fig_6f.tiff",
              width = 6, height = 6)
}
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
## 
## quartz_off_screen 
##                 2 
## quartz_off_screen 
##                 2 
## quartz_off_screen 
##                 2
```

We can see that at the 6D point there is a large difference in the proportion
of cell types, whereas this difference reduces at 12D.

## Generating confidence intervals

We now compute 95% confidence intervals around the differences in cell type
proportions at both timepoints.


```r
prop_table_cell_type <- function(x) {prop.table(table(factor(x, levels =
                                                             c("Basal", "TA",
                                                               "TD",
                                                               "unknown"))))}

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

obs_total_prop_diff_6d <- calc_cell_type_props(d6_cells@meta.data$cell_type_pred,
                                               d6_cells@meta.data$treatment_class)
obs_total_prop_diff_12d <- calc_cell_type_props(d12_cells@meta.data$cell_type_pred,
                                                d12_cells@meta.data$treatment_class)

prop_diffs_6d <- c() # Store results
prop_diffs_12d <- c() # Store results

for (r in 1:10000) {
    # Generate resample indices
    c6d_idx <- sample(which(d6_cells@meta.data$treatment_class == "Control"),
                      length(which(d6_cells@meta.data$treatment_class == "Control")),
                      replace = TRUE)
    r6d_idx <- sample(which(d6_cells@meta.data$treatment_class == "ROCKi"),
                      length(which(d6_cells@meta.data$treatment_class == "ROCKi"))-1,
                      replace = TRUE)
    d6_idx <- c(c6d_idx, r6d_idx)

    c12d_idx <- sample(which(d12_cells@meta.data$treatment_class == "Control"),
                      length(which(d12_cells@meta.data$treatment_class == "Control")),
                      replace = TRUE)
    r12d_idx <- sample(which(d12_cells@meta.data$treatment_class == "ROCKi"),
                      length(which(d12_cells@meta.data$treatment_class == "ROCKi"))-1,
                      replace = TRUE)
    d12_idx <- c(c12d_idx, r12d_idx)

    # Find cell type and treatment class labels for this resample
    d6_cell_type <- d6_cells@meta.data$cell_type_pred[d6_idx]
    d6_treat_class <- d6_cells@meta.data$treatment_class[d6_idx]

    d12_cell_type <- d12_cells@meta.data$cell_type_pred[d12_idx]
    d12_treat_class <- d12_cells@meta.data$treatment_class[d12_idx]

    total_prop_diff <- calc_cell_type_props(d6_cell_type, d6_treat_class)
    prop_diffs_6d <- c(prop_diffs_6d, total_prop_diff)

    total_prop_diff <- calc_cell_type_props(d12_cell_type, d12_treat_class)
    prop_diffs_12d <- c(prop_diffs_12d, total_prop_diff)
}
```


```r
alpha_vals <- c(0.975, 0.025)
residuals_6d <- prop_diffs_6d - obs_total_prop_diff_6d
print("95% confidence interval for total absolute difference in celltype proportions at 6D:")
```

```
## [1] "95% confidence interval for total absolute difference in celltype proportions at 6D:"
```

```r
print(as.numeric(obs_total_prop_diff_6d - quantile(residuals_6d, alpha_vals)))
```

```
## [1] 0.1796758 0.2539271
```

```r
residuals_12d <- prop_diffs_12d - obs_total_prop_diff_12d
print("95% confidence interval for total absolute difference in celltype proportions at 12D:")
```

```
## [1] "95% confidence interval for total absolute difference in celltype proportions at 12D:"
```

```r
print(as.numeric(obs_total_prop_diff_12d - quantile(residuals_12d, alpha_vals)))
```

```
## [1] -0.009077216  0.022292407
```


```r
signatures <- preprocess.signatures("celltype_markers_basal_ta_td.csv")

create_dotplot <- function(cells) {
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

skin_data@meta.data$sample_class_space <- skin_data@meta.data$sample_class
levels(skin_data@meta.data$sample_class_space) <- c("Control-6D ",
                                                    "ROCKi-6D ",
                                                    "Control-6D 6D ",
                                                    "ROCKi-6D 6D ")

holo_cells <- subset(skin_data, subset = cell_type_pred == "Basal")
ta_cells <- subset(skin_data, subset = cell_type_pred == "TA")
td_cells <- subset(skin_data, subset = cell_type_pred == "TD")

(holo_dp <- create_dotplot(holo_cells))
```

```r
(ta_dp <- create_dotplot(ta_cells))
```

```r
(td_dp <- create_dotplot(td_cells))
```

```r
if (save_images) {
    save_tiff(holo_dp, "fig_6g.tiff", width = 20, height = 8)
    save_tiff(ta_dp, "fig_s2a.tiff", width = 20, height = 8)
    save_tiff(td_dp, "fig_s2b.tiff", width = 20, height = 8)
}
```

```
## quartz_off_screen 
##                 2
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
ctrl6d_cells <- subset(d6_cells, treatment_class == "Control")
rocki6d_cells <- subset(d6_cells, treatment_class == "ROCKi")
ctrl12d_cells <- subset(d12_cells, treatment_class == "Control")
rocki12d_cells <- subset(d12_cells, treatment_class == "ROCKi")

ctrl6d_cells <- RunUMAP(ctrl6d_cells, dims = 1:50, reduction = "harmony")
rocki6d_cells <- RunUMAP(rocki6d_cells, dims = 1:50, reduction = "harmony")
ctrl12d_cells <- RunUMAP(ctrl12d_cells, dims = 1:50, reduction = "harmony")
rocki12d_cells <- RunUMAP(rocki12d_cells, dims = 1:50, reduction = "harmony")


prepare_cells_for_slingshot <- function(cells) {
    cells <- FindNeighbors(cells, reduction = "harmony", dims = 1:50)
    cells <- FindClusters(cells, resolution = 0.1)
    cells <- RenameIdents(object = cells, "0" = "TA", "1" = "TD", "2" = "Basal")
    
    return (cells)
}

d6_cells <- prepare_cells_for_slingshot(d6_cells)
d12_cells <- prepare_cells_for_slingshot(d12_cells)
ctrl6d_cells <- prepare_cells_for_slingshot(ctrl6d_cells)
rocki6d_cells <- prepare_cells_for_slingshot(rocki6d_cells)
ctrl12d_cells <- prepare_cells_for_slingshot(ctrl12d_cells)
rocki12d_cells <- prepare_cells_for_slingshot(rocki12d_cells)
```

## Computing Trajectories


```r
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

d6_cells.sce <- as.SingleCellExperiment(d6_cells)
d12_cells.sce <- as.SingleCellExperiment(d12_cells)
ctrl6d_cells.sce <- as.SingleCellExperiment(ctrl6d_cells)
rocki6d_cells.sce <- as.SingleCellExperiment(rocki6d_cells)
ctrl12d_cells.sce <- as.SingleCellExperiment(ctrl12d_cells)
rocki12d_cells.sce <- as.SingleCellExperiment(rocki12d_cells)

# Change UMAP values to make figures all similarly oriented
reducedDim(rocki6d_cells.sce, "UMAP")[,2] <- (-1) * reducedDims(rocki6d_cells.sce)$UMAP[,2]
reducedDim(ctrl12d_cells.sce, "UMAP")[,1] <- (-1) * reducedDims(ctrl12d_cells.sce)$UMAP[,1]

sshot_6d <- slingshot(d6_cells.sce, d6_cells$seurat_clusters,
                      reducedDim = "HARMONY", start.clus = "2")
sshot_12d <- slingshot(d12_cells.sce, d12_cells$seurat_clusters,
                       reducedDim = "HARMONY", start.clus = "2")

sshot_ctrl6d_separate <- slingshot(ctrl6d_cells.sce,
                                   ctrl6d_cells$seurat_clusters,
                                   reducedDim = "UMAP", start.clus = "2")
sshot_rocki6d_separate <- slingshot(rocki6d_cells.sce,
                                    rocki6d_cells$seurat_clusters,
                                    reducedDim = "UMAP", start.clus = "2")
sshot_ctrl12d_separate <- slingshot(ctrl12d_cells.sce,
                                    ctrl12d_cells$seurat_clusters,
                                    reducedDim = "UMAP", start.clus = "2")
sshot_rocki12d_separate <- slingshot(rocki12d_cells.sce,
                                     rocki12d_cells$seurat_clusters,
                                     reducedDim = "UMAP", start.clus = "2")

umap_curve_6d <- slingCurves(embedCurves(sshot_6d,
                                         reducedDims(d6_cells.sce)$UMAP))[[1]]
umap_curve_12d <- slingCurves(embedCurves(sshot_12d,
                                          reducedDims(d12_cells.sce)$UMAP))[[1]]

umap_curve_ctrl6d <- slingCurves(sshot_ctrl6d_separate)[[1]]
umap_curve_rocki6d <- slingCurves(sshot_rocki6d_separate)[[1]]
umap_curve_ctrl12d <- slingCurves(sshot_ctrl12d_separate)[[1]]
umap_curve_rocki12d <- slingCurves(sshot_rocki12d_separate)[[1]]

ctrl6d_idx <- which(d6_cells@meta.data$treatment_class == "Control")
rocki6d_idx <- which(d6_cells@meta.data$treatment_class == "ROCKi")
ctrl12d_idx <- which(d12_cells@meta.data$treatment_class == "Control")
rocki12d_idx <- which(d12_cells@meta.data$treatment_class == "ROCKi")

sshot_ctrl6d_together <- sshot_6d[,ctrl6d_idx]
sshot_rocki6d_together <- sshot_6d[,rocki6d_idx]
sshot_ctrl12d_together <- sshot_12d[,ctrl12d_idx]
sshot_rocki12d_together <- sshot_12d[,rocki12d_idx]

xmin_sep <- min(reducedDim(ctrl6d_cells.sce, "UMAP")[,1],
                reducedDim(rocki6d_cells.sce, "UMAP")[,1],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,1],
                reducedDim(rocki12d_cells.sce, "UMAP")[,1])
xmax_sep <- max(reducedDim(ctrl6d_cells.sce, "UMAP")[,1],
                reducedDim(rocki6d_cells.sce, "UMAP")[,1],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,1],
                reducedDim(rocki12d_cells.sce, "UMAP")[,1])
ymin_sep <- min(reducedDim(ctrl6d_cells.sce, "UMAP")[,2],
                reducedDim(rocki6d_cells.sce, "UMAP")[,2],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,2],
                reducedDim(rocki12d_cells.sce, "UMAP")[,2])
ymax_sep <- max(reducedDim(ctrl6d_cells.sce, "UMAP")[,2],
                reducedDim(rocki6d_cells.sce, "UMAP")[,2],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,2],
                reducedDim(rocki12d_cells.sce, "UMAP")[,2])

create_traj_umap <- function(sshot_to_plot, umap_curve, cell_types) {
    points_df <- data.frame(u1 = reducedDims(sshot_to_plot)$UMAP[,1],
                            u2 = reducedDims(sshot_to_plot)$UMAP[,2],
                            cell_types=cell_types)
    points_df$cell_types <- factor(points_df$cell_types,
                                   levels = c("Basal", "TA", "TD", "unknown"))
    curve_df <- data.frame(u1 = umap_curve$s[umap_curve$ord, ][,1],
                           u2 = umap_curve$s[umap_curve$ord, ][,2])

    arrowhead_length <- unit(0.2, "cm")
    arrow_width <- 1

    ptimes <- sshot_to_plot$slingPseudotime_1
    basal_ptimes <- ptimes[which(cell_types == "Basal")]
    ta_ptimes <- ptimes[which(cell_types == "TA")]
    td_ptimes <- ptimes[which(cell_types == "TD")]

    percentile <- 0.97
    basal_percentile_ptime_idx <- which(abs(quantile(basal_ptimes, percentile) - basal_ptimes) ==
                                        min(abs(quantile(basal_ptimes, percentile) - basal_ptimes)))
    ta_percentile_ptime_idx <- which(abs(quantile(ta_ptimes, percentile) - ta_ptimes) ==
                                     min(abs(quantile(ta_ptimes, percentile) - ta_ptimes)))
    td_percentile_ptime_idx <- which(abs(quantile(td_ptimes, percentile) - td_ptimes) ==
                                     min(abs(quantile(td_ptimes, percentile) - td_ptimes)))

    basal_percentile_ptime <- basal_ptimes[basal_percentile_ptime_idx]
    ta_percentile_ptime <- ta_ptimes[ta_percentile_ptime_idx]
    td_percentile_ptime <- td_ptimes[td_percentile_ptime_idx]

    p <- ggplot() +
         geom_point(data = points_df, aes(x = u1, y = u2, color = cell_types),
                    size = -0.1) +
         scale_color_manual(labels=c("Holoclone-forming",
                           "Mero- or Paraclone-forming",
                           "Differentiated", "Unknown"),
                           values = c("#80C9EA", "#DD6E79",
                                      "#43863E", "#BBBBBB")) +
        geom_path(data = curve_df[1:round(length(sshot_to_plot$slingPseudotime_1) * percentile), ],
                  aes(x = u1, y = u2),
                  arrow = arrow(type = "closed", length = arrowhead_length),
                  size = arrow_width) +
         geom_path(data = curve_df[1:2,],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "open", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[1:which(sort(ptimes) == basal_percentile_ptime), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[1:which(sort(ptimes) == ta_percentile_ptime), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         coord_fixed() + umap_theme

     return(p)
}

ctrl6d_traj_umap <- create_traj_umap(sshot_ctrl6d_separate, umap_curve_ctrl6d,
                                     Idents(ctrl6d_cells))
rocki6d_traj_umap <- create_traj_umap(sshot_rocki6d_separate, umap_curve_rocki6d,
                                      Idents(rocki6d_cells))
ctrl12d_traj_umap <- create_traj_umap(sshot_ctrl12d_separate, umap_curve_ctrl12d,
                                      Idents(ctrl12d_cells))
rocki12d_traj_umap <- create_traj_umap(sshot_rocki12d_separate, umap_curve_rocki12d,
                                       Idents(rocki12d_cells))

(ctrl6d_traj_umap_no_legend <- ctrl6d_traj_umap + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

```r
(rocki6d_traj_umap_no_legend <- rocki6d_traj_umap + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-17-2.png" style="display: block; margin: auto;" />

```r
(ctrl12d_traj_umap_no_legend <- ctrl12d_traj_umap + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-17-3.png" style="display: block; margin: auto;" />

```r
(rocki12d_traj_umap_no_legend <- rocki12d_traj_umap + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-17-4.png" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(ctrl6d_traj_umap_no_legend, "fig_7a.tiff", width = 10, height = 8)
    save_tiff(rocki6d_traj_umap_no_legend, "fig_7b.tiff", width = 10, height = 8)
    save_tiff(ctrl12d_traj_umap_no_legend, "fig_7c.tiff", width = 10, height = 8)
    save_tiff(rocki12d_traj_umap_no_legend, "fig_7d.tiff", width = 10, height = 8)
}
```

## Comparing trajectories

We compared trajectories by sampling 2500 pseudotimes uniformly at random for
each sample class, sorting them, and then plotting them against one another. We
also computed the ???speed??? of the differentiation at each point by calculating
the difference between each neighbouring pair of ordered pseudotimes and
plotting the moving average of 151 pseudotimes, in order to reduce noise. We
take cell type of a cell to be the cluster to which it belongs to improve the
appearance of the trajectory plots when colouring by cell type.


```r
d6_cells@meta.data$cell_type <- Idents(d6_cells)
d12_cells@meta.data$cell_type <- Idents(d12_cells)

# mov_avg function copied from stackoverflow.com/a/4862334/6914552
mov_avg <- function(x, n = 151){stats::filter(x, rep(1 / n, n), sides = 2)}

create_traj_df <- function (ptimes, name, cell_types) {
    quantile_points <- seq(0.001,1,0.001)
    quantiled_data <- as.numeric(quantile(ptimes, quantile_points))
    sorted_cell_types <- cell_types[order(ptimes)]
    num_points <- length(ptimes)
    cell_indices <- round(num_points * quantile_points)
    cell_types_subset <- sorted_cell_types[cell_indices]

    return(data.frame(pct = quantile_points,
                      ptime = quantiled_data,
                      cell_type = factor(cell_types_subset,
                                         levels = c("Basal", "TA", "TD", "unknown")),
                      diff = c(NA, mov_avg(diff(quantiled_data))),
                      sample_class = rep(name, length(quantile_points))))
}

traj_df_ctrl6d <- create_traj_df(sshot_ctrl6d_together$slingPseudotime_1, "Control-6D",
                                 d6_cells[,ctrl6d_idx]@meta.data$cell_type)
traj_df_rocki6d <- create_traj_df(sshot_rocki6d_together$slingPseudotime_1, "ROCKi-6D",
                                  d6_cells[,rocki6d_idx]@meta.data$cell_type)
traj_df_6d <- rbind(traj_df_ctrl6d, traj_df_rocki6d)

traj_df_ctrl12d <- create_traj_df(sshot_ctrl12d_together$slingPseudotime_1,
                                  "Control-12D",
                                  d12_cells[,ctrl12d_idx]@meta.data$cell_type)
traj_df_rocki12d <- create_traj_df(sshot_rocki12d_together$slingPseudotime_1,
                                   "ROCKi-12D",
                                   d12_cells[,rocki12d_idx]@meta.data$cell_type)
traj_df_12d <- rbind(traj_df_ctrl12d, traj_df_rocki12d)

obs_result_ptime_6d <- mean(abs(traj_df_ctrl6d$ptime - traj_df_rocki6d$ptime))
obs_result_ptime_12d <- mean(abs(traj_df_ctrl12d$ptime - traj_df_rocki12d$ptime))

obs_result_speed_6d <- mean(abs(traj_df_ctrl6d$diff - traj_df_rocki6d$diff),
                            na.rm = TRUE)
obs_result_speed_12d <- mean(abs(traj_df_ctrl12d$diff - traj_df_rocki12d$diff),
                             na.rm = TRUE)

traj_plot_margin_theme <- theme(plot.margin = margin(0.1, 1.2, 0.1, 1.2, "cm"))
diff_comparison_plots_theme <- theme(legend.title = element_blank())

point_size <- 1
create_traj_plot <- function(traj_df) {
    shade_alpha <- 0.2

    basal_upper_lim <- quantile(traj_df[traj_df$cell_type == "Basal",]$ptime, 0.75)
    ta_upper_lim <- quantile(traj_df[traj_df$cell_type == "TA",]$ptime, 0.95)

    basal_color <- "#80C9EA"
    ta_color <- "#DD6E79"
    td_color <- "#43863E"

    traj_plot <- ggplot(traj_df, aes(x = pct, y = ptime, group = sample_class)) +
                    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
                             ymax = basal_upper_lim, alpha = shade_alpha,
                             fill = basal_color) +
                    annotate("rect", xmin = -Inf, xmax = Inf,
                             ymin = basal_upper_lim, ymax = ta_upper_lim,
                             alpha = shade_alpha, fill = ta_color) +
                    annotate("rect", xmin = -Inf, xmax = Inf,
                             ymin = ta_upper_lim, ymax = Inf,
                             alpha = shade_alpha, fill = td_color) +
                    geom_point(aes(color = cell_type), size = point_size) +
                    scale_color_manual(values = c(basal_color, ta_color, td_color)) +
                    new_scale_color() +
                    geom_line(aes(linetype = sample_class), size = 0.4) +
                    scale_linetype_manual(values = c("solid", "dashed")) +
                    scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
                    ylim(c(0, max(traj_df$ptime))) +
                    xlab("Percentage through differentation process") +
                    ylab("Pseudotime") + traj_plot_margin_theme +
                    diff_comparison_plots_theme + axis_text_theme

    return(traj_plot)
}

traj_plot_6d <- create_traj_plot(traj_df_6d)
traj_plot_12d <- create_traj_plot(traj_df_12d)

(traj_plot_6d_no_legend <- traj_plot_6d + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-18-1.png" width="100%" style="display: block; margin: auto;" />

```r
(traj_plot_12d_no_legend <- traj_plot_12d + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-18-2.png" width="100%" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(traj_plot_6d_no_legend, "fig_7e.tiff", width = 10, height = 8)
    save_tiff(traj_plot_12d_no_legend, "fig_7f.tiff", width = 10, height = 8)
}


create_speed_plot <- function(traj_df) {
    speed_plot <- ggplot(traj_df, aes(x = pct, y = diff, group = sample_class)) +
                  geom_point(aes(color = cell_type), size = point_size) +
                  scale_color_manual(values = c("#80C9EA", "#DD6E79", "#43863E")) +
                  new_scale_color() +
                  geom_line(aes(linetype = sample_class), size = 0.4) +
                  scale_linetype_manual(values = c("solid", "dashed")) +
                  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
                  xlab("Percentage through differentation process") +
                  ylab("Differentiation speed") + traj_plot_margin_theme +
                  diff_comparison_plots_theme + axis_text_theme

    return(speed_plot)
}


(speed_plot_6d <- create_speed_plot(traj_df_6d) + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-18-3.png" width="100%" style="display: block; margin: auto;" />

```r
(speed_plot_12d <- create_speed_plot(traj_df_12d) + theme(legend.position = "none"))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-18-4.png" width="100%" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(speed_plot_6d, "fig_s3a.tiff", width = 14, height = 8)
    save_tiff(speed_plot_12d, "fig_s3b.tiff", width = 14, height = 8)
}


for_dot_colour_legend <- ggplot(traj_df_6d,
                                aes(x = pct, y = ptime, group = sample_class)) +
                         geom_point(aes(color = cell_type), size = 3) +
                         scale_color_manual(values = c("#80C9EA", "#DD6E79", "#43863E"),
                                            labels = c("Holoclone-forming",
                                                       "Mero- or Paraclone-forming",
                                                       "Differentiated")) +
                         theme(legend.title = element_blank(),
                               legend.text = element_text(face = "bold"))
dot_colour_legend <- as_ggplot(get_legend(for_dot_colour_legend, "right"))

for_linetype_legend <- ggplot(traj_df_6d,
                                 aes(x = pct, y = ptime, group = sample_class)) +
                          geom_line(aes(linetype = sample_class), size = 0.4) +
                          scale_linetype_manual(values = c("solid", "dashed"),
                                                labels = c("Control", "ROCKi-treated")) +
                          theme(legend.title = element_blank(),
                                legend.text = element_text(face = "bold"))
linetype_legend <- as_ggplot(get_legend(for_linetype_legend, "right"))

(legends <- ggarrange(linetype_legend, dot_colour_legend, ncol = 2))
```

<img src="data_analysis_for_paper_files/figure-html/unnamed-chunk-18-5.png" width="100%" style="display: block; margin: auto;" />

```r
if (save_images) {
    save_tiff(legends, "fig_7_legends.tiff", width = 10, height = 8)
}
```

## Generating confidence intervals


```r
# Generate bootstrapped data as in "experiment_code_trajectories.R" and save
# the outputs as "trajectory_comparison_results.txt"

bootstrap_results <- read.table("trajectory_comparison_results.txt")
bootstrap_results <- bootstrap_results[, c(3, 4, 5, 6)]
colnames(bootstrap_results) <- c("TrajDiff6D", "TrajDiff12D", "SpeedDiff6D",
                                 "SpeedDiff12D")

diffs_quants_traj6d <- quantile((bootstrap_results$TrajDiff6D - obs_result_ptime_6d),
                                alpha_vals)
residuals_traj6d <- as.numeric(obs_result_ptime_6d - diffs_quants_traj6d)
print(paste("95% confidence interval for difference in pseudotimes at 6D:",
            paste(residuals_traj6d, collapse = " ")))
```

```
## [1] "95% confidence interval for difference in pseudotimes at 6D: 12.4460321784284 15.5195479079284"
```

```r
diffs_quants_traj12d <- quantile((bootstrap_results$TrajDiff12D - obs_result_ptime_12d),
                                alpha_vals)
residuals_traj12d <- as.numeric(obs_result_ptime_12d - diffs_quants_traj12d)
print(paste("95% confidence interval for difference in pseudotimes at 12D:",
            paste(residuals_traj12d, collapse = " ")))
```

```
## [1] "95% confidence interval for difference in pseudotimes at 12D: 1.65452072266176 3.97100868716176"
```

```r
diffs_quants_speed6d <- quantile((bootstrap_results$SpeedDiff6D - obs_result_speed_6d),
                                  alpha_vals)
residuals_speed6d <- as.numeric(obs_result_speed_6d - diffs_quants_speed6d)
print(paste("95% confidence interval for difference in speed at 6D:",
            paste(residuals_speed6d, collapse = " ")))
```

```
## [1] "95% confidence interval for difference in speed at 6D: 0.0366647847933335 0.0519434975433335"
```

```r
diffs_quants_speed12d <- quantile((bootstrap_results$SpeedDiff12D - obs_result_speed_12d),
                                  alpha_vals)
residuals_speed12d <- as.numeric(obs_result_speed_12d - diffs_quants_speed12d)
print(paste("95% confidence interval for difference in speed at 12D:",
            paste(residuals_speed12d, collapse = " ")))
```

```
## [1] "95% confidence interval for difference in speed at 12D: 0.0146296286652594 0.0317649334152594"
```


# References {#endnotes}

[^1]: Zhang Z, Luo D, Zhong X, et al. _SCINA: A Semi-Supervised Subtyping Algorithm of Single Cells and Bulk Samples_. Genes (Basel). 2019;10(7):531.

[^2]: Enzo, E., Secone Seconetti, A., Forcato, M. et al. _Single-keratinocyte transcriptomic analyses identify different clonal types and proliferative potential mediated by FOXM1 in human epidermal stem cells_. Nat Commun 12, 2505 (2021).

[^3]: Street, K., Risso, D., Fletcher, R. B., et al. _Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics_. BMC Genomics 19(477), 2018.

# Appendix

## All code for this report

<!-- Suggested in https://bookdown.org/yihui/rmarkdown-cookbook/code-appendix.html -->

```r
knitr::opts_chunk$set(results = 'hide', message = FALSE, warning = FALSE,
                      fig.align = 'center')

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
theme_update(plot.title = element_text(face="bold", size=12, hjust=0.5))
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
skin_data <- CreateSeuratObject(counts = skin_data,
                                project = "3DonorEffectOfRocki",
                                min.cells = 3, min.features = 200,
                                names.delim = "-", names.field = 2)
skin_data <- RenameIdents(object = skin_data, "1" = "1Control-6D ",
                          "2" = "1ROCKi-6D ", "3" = "1Control-6D 6D ",
                          "4" = "1ROCKi-6D 6D ", "5" = "2Control-6D ",
                          "6" = "2ROCKi-6D ", "7" = "2Control-6D 6D ",
                          "8" = "2ROCKi-6D 6D ", "9" = "3Control-6D ",
                          "10" = "3ROCKi-6D ", "11" = "3Control-6D 6D ",
                          "12" = "3ROCKi-6D 6D ")
# Compute percentage of mitochondrial DNA in each cell
skin_data[["percent.mt"]] <- PercentageFeatureSet(skin_data, pattern = "^MT-")

VlnPlot_limits <- function(data, metric, ylab, limits) {
    df <- data.frame(x=Idents(data), y = data@meta.data[,metric])
    plot <- ggplot() + geom_violin(data = df, aes(x=x, y=y, fill=x)) +
            geom_errorbar(data = limits, aes(x = x, ymin = lower, ymax = upper),
                          color = "red") +
            labs(y=ylab) +
            scale_y_continuous(labels = scales::comma) +
            theme(axis.text.x = element_text(angle = 45, hjust=1),
                  legend.position = "none", plot.title = element_blank(),
                  axis.title.x = element_blank()) + axis_text_theme

    return(plot)
}

ncount_df <- data.frame(x = levels(Idents(skin_data)),
                        lower = c(8000, 15000, 8000, 8000,
                                  8000, 8000, 5000, 8000,
                                  8000, 8000, 5000, 5000),
                        upper = c(50000, 62500, 50000, 50000,
                                  50000, 50000, 50000, 50000,
                                  50000, 50000, 40000, 40000))

(ncount_vln <- VlnPlot_limits(skin_data, "nCount_RNA",
                               "Number of RNA\nmolecules per cell",
                               ncount_df) + ylim(c(0, 100000)))
if(save_images) {
    save_tiff(ncount_vln, "fig_s4a.tiff")
}

nfeature_df <- data.frame(x = levels(Idents(skin_data)),
                     lower = c(2500, 3000, 2500, 2500,
                               2500, 2500, 2000, 2500,
                               2500, 2500, 1500, 1500),
                     upper = c(7500, 8500, 7500, 7500,
                               7500, 7500, 7500, 7500,
                               7500, 7500, 6000, 6000))

(nfeature_vln <- VlnPlot_limits(skin_data, "nFeature_RNA",
                               "Number of genes\ndetected per cell",
                               nfeature_df))
if(save_images) {
    save_tiff(nfeature_vln, "fig_s4b.tiff")
}

pct_mito_df <- data.frame(x = levels(Idents(skin_data)),
                        lower = c(rep(2, 8), rep(3, 4)),
                        upper = c(rep(10, 8), rep(12, 4)))

(pct_mito_vln <- VlnPlot_limits(skin_data, "percent.mt",
                               "Percentage of mitochrondrial\nDNA detected per cell",
                               pct_mito_df) + ylim(c(0, 20)))
if(save_images) {
    save_tiff(pct_mito_vln, "fig_s4c.tiff")
}

pass_qc <- c()
for (j in 1:length(skin_data@meta.data$orig.ident)) {
    i <- as.integer(as.character(skin_data@meta.data$orig.ident[j]))
    pass_qc <- c(pass_qc, all(skin_data@meta.data$nFeature_RNA[j] >= nfeature_df$lower[i],
                              skin_data@meta.data$nFeature_RNA[j] <= nfeature_df$upper[i],
                              skin_data@meta.data$nCount_RNA[j] >= ncount_df$lower[i],
                              skin_data@meta.data$nCount_RNA[j] <= ncount_df$upper[i],
                              skin_data@meta.data$percent.mt[j] >= pct_mito_df$lower[i],
                              skin_data@meta.data$percent.mt[j] <= pct_mito_df$upper[i]))
}
skin_data <- skin_data[, pass_qc]
# Remove melanocytes
skin_data <- subset(skin_data, subset = MLANA > 1 | PMEL > 1 | MITF > 1,
                    invert = T)
# Remove mesenchymal cells
skin_data <- subset(skin_data,
                    subset = MTRNR2L6 > 1 | MTRNR2L10 > 1 | MTRNR2L7 > 1 |
                             MTRNR2L1 > 1 | PRKAR2B > 1 | NR2F1 > 1,
                    invert = T)
# Remove fibroblasts
skin_data <- subset(skin_data, subset = ACTG2 > 1 | DLK1 > 1, invert = T)
skin_data <- SCTransform(skin_data)
skin_data.outliers_removed <- skin_data
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
skin_data <- CellCycleScoring(skin_data, s.features = s.genes,
                              g2m.features = g2m.genes)
skin_data@meta.data$CC.difference <- skin_data@meta.data$S.Score - skin_data@meta.data$G2M.Score
skin_data <- SCTransform(skin_data, vars.to.regress = c("CC.difference"))
skin_data.cell_cycle_regressed <- skin_data
# Run PCA
skin_data <- RunPCA(skin_data, features = skin_data@assays$SCT@var.features,
                    npcs = 50)
(elbow_plot <- ElbowPlot(skin_data, ndims = 50) + axis_text_theme)
if(save_images) {
    save_tiff(elbow_plot, "fig_s5.tiff")
}

# Remove inter-donor differences using Harmony
donor_idents <- sapply(skin_data@meta.data$orig.ident,
                       function(x) if (x %in% c("1","2","3","4")) {"1"}
                                   else if (x %in% c("5","6","7","8")) {"2"}
                                   else {"3"})
skin_data@meta.data["donor_idents"] <- donor_idents

skin_data <- skin_data %>% RunHarmony("donor_idents", assay.use = "SCT")
skin_data.harmonised <- skin_data

sample_class <- sapply(skin_data@meta.data$orig.ident,
                       function(x) if (x %in% c("1","5","9")) {"Control-6D"}
                                   else if (x %in% c("2","6","10")) {"ROCKi-6D"}
                                   else if (x %in% c("3","7","11")) {"Control-12D"}
                                   else {"ROCKi-12D"})
skin_data@meta.data$sample_class <- sample_class
skin_data@meta.data$sample_class <- factor(skin_data@meta.data$sample_class,
                                           levels = c("Control-6D", "ROCKi-6D",
                                                      "Control-12D",
                                                      "ROCKi-12D"))

treatment_class <- sapply(skin_data@meta.data$sample_class,
                       function(x) if (x %in% c("Control-6D","Control-12D")) {"Control"}
                                   else {"ROCKi"})
skin_data@meta.data$treatment_class <- treatment_class
d6_cells <- subset(skin_data,
                   subset = sample_class %in% c("Control-6D", "ROCKi-6D"))
d12_cells <- subset(skin_data,
                    subset = sample_class %in% c("Control-12D", "ROCKi-12D"))

d6_cells <- RunUMAP(d6_cells, dims = 1:50, reduction = "harmony")
d12_cells <- RunUMAP(d12_cells, dims = 1:50, reduction = "harmony")

umap_theme <- theme(axis.line = element_blank(), axis.text.x = element_blank(),
                    axis.text.y=element_blank(), axis.ticks = element_blank(),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    plot.title = element_blank())
chart_margin_theme <- theme(plot.margin = margin(0, 1.2, 0, 1.2, "cm"))

labels <- c("Control", "ROCKi-treated")
umap_legend_plot <- ggplot(data.frame(class = as.character(1:length(labels)),
                                      x = 1:length(labels)),
                           aes(x = x, y = x, col = class)) +
                    geom_point() +
                    scale_color_manual(labels = labels,
                                       values = c(blue_color, red_color)) +
                    theme(legend.title = element_blank(),
                          legend.key = element_rect(fill = NA, color = NA),
                          legend.margin = margin(c(0,0,0,0)),
                          legend.text = element_text(face = "bold"),
                          legend.background = element_rect(fill="#00000000",
                                                           size=0.5,
                                                           linetype="solid")) +
                    guides(colour = guide_legend(override.aes = list(size=2)))
umap_legend <- as_ggplot(get_legend(umap_legend_plot, "left"))
if(save_images) {
    save_tiff(umap_legend, "umap_legend.tiff")
}

plot_treatment_class_umap <- function(cells) {
    p <- DimPlot(cells, group.by = "treatment_class",
                 cols = c(blue_color, red_color),
                 pt.size = -0.1) +
         coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

(umap_d6 <- plot_treatment_class_umap(d6_cells))
if(save_images) {
    save_tiff(umap_d6, "fig_5a.tiff")
}

(umap_d12 <- plot_treatment_class_umap(d12_cells))
if(save_images) {
    save_tiff(umap_d12, "fig_5b.tiff")
}
daseq_out_6d <- getDAcells(d6_cells@reductions$harmony@cell.embeddings,
                           cell.labels = as.character(d6_cells@meta.data$sample_class),
                           labels.1 = c("Control-6D"),
                           labels.2 = c("ROCKi-6D"),
                           plot.embedding = d6_cells@reductions$umap@cell.embeddings)

daseq_out_12d <- getDAcells(d12_cells@reductions$harmony@cell.embeddings,
                            cell.labels = as.character(d12_cells@meta.data$sample_class),
                            labels.1 = c("Control-12D"),
                            labels.2 = c("ROCKi-12D"),
                            plot.embedding = d12_cells@reductions$umap@cell.embeddings)

daseq_legend <- as_ggplot(get_legend(daseq_out_12d$pred.plot +
                                     umap_theme +
                                     theme(legend.position = "top",
                                           legend.text = element_blank())))

if(save_images) {
    save_tiff(daseq_legend, "daseq_legend.tiff")
}

create_daseq_umap <- function(daseq_out) {
    daseq_umap <- daseq_out$pred.plot + umap_theme + theme(legend.position = "none")

    return(daseq_umap)
}

(daseq_umap_6d <- create_daseq_umap(daseq_out_6d))
if(save_images) {
    save_tiff(daseq_umap_6d, "fig_5c.tiff")
}

(daseq_umap_12d <- create_daseq_umap(daseq_out_12d))
if(save_images) {
    save_tiff(daseq_umap_12d, "fig_5d.tiff")
}

proportion_da_6d <- length(c(daseq_out_6d$da.down, daseq_out_6d$da.up)) / length(colnames(d6_cells))
proportion_da_12d <- length(c(daseq_out_12d$da.down, daseq_out_12d$da.up)) / length(colnames(d12_cells))
print(paste("Proportion of cells in differentially abundant regions at 6D:",
            proportion_da_6d))
print(paste("Proportion of cells in differentially abundant regions at 12D:",
            proportion_da_12d))
Idents(skin_data) <- skin_data@meta.data$sample_class
marks_6d <- FindMarkers(skin_data, ident.1 = "ROCKi-6D",
                        ident.2 = "Control-6D")
marks_12d <- FindMarkers(skin_data, ident.1 = "ROCKi-12D",
                         ident.2 = "Control-12D")

not_sig_color <- "#4D4D4D"

VolPlot <- function(marks) {
    marks$p_val_adj_neg_log <- (-1) * log10(marks$p_val_adj)
    marks$p_val_adj_neg_log[is.infinite(marks$p_val_adj_neg_log)] <- 300
    sig_marks_pos <- subset(marks, avg_log2FC > 1)
    sig_marks_neg <- subset(marks, avg_log2FC < -1)

    # Assign classes to genes based on how DE they are. Used to set colours
    de_class <- ifelse(marks$avg_log2FC < -1 & marks$p_val_adj < 0.01, 'Downreg',
                  ifelse(marks$avg_log2FC > 1 & marks$p_val_adj < 0.01, 'Upreg',
                         ifelse(marks$p_val_adj < 0.01, "Low-LFC", "Not-Sig")))
    de_class <- factor(de_class, levels = c("Not-Sig", "Low-LFC",
                                            "Downreg", "Upreg"))

    x_sublabel_y <- -50 # Position the sub-labels of the x-axis correctly
    plot <- ggplot(marks, aes(x = avg_log2FC, y = p_val_adj_neg_log,
                              label = rownames(marks))) +
            geom_point(aes(color = de_class, fill = de_class), size = 2,
                       shape = 21, alpha = 0.5) +
            scale_color_manual(values = c("gray", not_sig_color, "red", "white")) +
            scale_fill_manual(values = c("gray", not_sig_color, "white", "red")) +
            geom_vline(xintercept = -1, linetype = "dashed", alpha = 0.5) +
            geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
            xlim(-2.5, 2.5) + ylim(-log10(0.01), 350) +
            xlab(bquote(bold(Log["2"]~"fold change"))) +
            ylab(bquote(bold(-Log["10"] ~ P))) +
            theme(legend.position="none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black")) +
            axis_text_theme +
            coord_cartesian(clip = "off")

    return (plot)
}

(vol_6d <- VolPlot(marks_6d))
if(save_images) {
    save_tiff(vol_6d, "fig_5e.tiff", height = 5, width = 10)
}

labels <- c(bquote(bold(Log["2"]~"fold change ??? -1")),
            bquote(bold("|"*Log["2"]~"fold change| < 1")),
            bquote(bold(Log["2"]~"fold change ??? 1")))

vol_6d_legend_plot <- ggplot(data.frame(class = as.character(1:length(labels)),
                                        x = 1:length(labels)),
                             aes(x = x, y = x, color = class, fill = class)) +
                      geom_point(size = 3, shape = 21, alpha = 0.5) +
                      scale_fill_manual(values = c("white", not_sig_color, "red"),
                                        labels = labels) +
                      scale_color_manual(values = c("red", not_sig_color, "white"),
                                         labels = labels) +
                      theme(legend.title = element_blank(),
                            legend.key = element_rect(fill = NA, color = NA),
                            legend.text = element_text(face = "bold"),
                            legend.margin=margin(unit(c(0, 0, 0, 0), "cm"))) +
                      guides(colour = guide_legend(override.aes = list(size = 2)))

volplot_legend <- as_ggplot(get_legend(vol_6d_legend_plot, "right"))
if(save_images) {
    save_tiff(volplot_legend, "volplot_legend.tiff")
}

(vol_12d <- VolPlot(marks_12d))
if(save_images) {
    save_tiff(vol_12d, "fig_5f.tiff", height = 5, width = 10)
}
(signatures <- preprocess.signatures("celltype_markers_basal_ta_td.csv"))

skin_data@meta.data$cell_type_pred <- SCINA(as.matrix(GetAssayData(skin_data)),
                                            signatures)$cell_labels
d6_cells@meta.data$cell_type_pred <- SCINA(as.matrix(GetAssayData(d6_cells)),
                                           signatures)$cell_labels
d12_cells@meta.data$cell_type_pred <- SCINA(as.matrix(GetAssayData(d12_cells)),
                                            signatures)$cell_labels

d6_cells@meta.data$cell_type_pred <- factor(d6_cells@meta.data$cell_type_pred,
                                             levels = c("Basal", "TA", "TD", "unknown"))
d12_cells@meta.data$cell_type_pred <- factor(d12_cells@meta.data$cell_type_pred,
                                             levels = c("Basal", "TA", "TD", "unknown"))

create_cell_type_umap <- function(cells) {
    # Force factors back into correct order as this inexplicably disappears
    # when passed to a function
    cells@meta.data$cell_type_pred <- factor(cells@meta.data$cell_type_pred,
                                             levels = c("Basal", "TA", "TD", "unknown"))
    p <- DimPlot(cells, group.by = "cell_type_pred") +
    scale_color_manual(labels=c("Holoclone-forming",
                               "Mero- or Paraclone-forming",
                               "Differentiated", "Unclassified keratinocytes"),
                      values = c("#80C9EA", "#DD6E79", "#43863E", "#BBBBBB")) +
    coord_fixed() + umap_theme + theme(legend.position = "none")

    return(p)
}

ctrl6d_cells <- subset(d6_cells, treatment_class == "Control")
rocki6d_cells <- subset(d6_cells, treatment_class == "ROCKi")
ctrl12d_cells <- subset(d12_cells, treatment_class == "Control")
rocki12d_cells <- subset(d12_cells, treatment_class == "ROCKi")

cell_type_umap_ctrl6d <- create_cell_type_umap(ctrl6d_cells)
cell_type_umap_rocki6d <- create_cell_type_umap(rocki6d_cells)
cell_type_umap_ctrl12d <- create_cell_type_umap(ctrl12d_cells)
cell_type_umap_rocki12d <- create_cell_type_umap(rocki12d_cells)

if (save_images) {
    save_tiff(cell_type_umap_ctrl6d, "fig_6a.tiff")
    save_tiff(cell_type_umap_ctrl12d, "fig_6b.tiff")
    save_tiff(cell_type_umap_rocki6d, "fig_6c.tiff")
    save_tiff(cell_type_umap_rocki12d, "fig_6d.tiff")
}

cell_type_df_6d <- data.frame(CellType = d6_cells@meta.data$cell_type_pred,
                           Class = d6_cells@meta.data$sample_class)
cell_type_df_6d$Class <- factor(cell_type_df_6d$Class,
                                levels = c("Control-6D", "ROCKi-6D",
                                           "Control-12D", "ROCKi-12D"))

create_cell_type_proportion_graph <- function(df, names) {
    p <- ggplot(df, aes(x = Class, fill = CellType)) +
    geom_bar(position = "fill") + ylab("Proportion of cell type") +
    labs(fill = "Cell type") +
    scale_fill_manual(labels = c("Holoclone-forming",
                                 "Mero- or Paraclone-forming",
                                 "Differentiated",
                                 "Unclassified keratinocytes"),
                      values = c("#80C9EA", "#DD6E79", "#43863E", "#BBBBBB")) +
    scale_x_discrete(labels = names) +
    chart_margin_theme + axis_text_theme +
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1),
          axis.title.x = element_blank())

    return(p)
}

cell_type_proportion_graph_6d <- create_cell_type_proportion_graph(cell_type_df_6d, c("Control-6D ", "ROCKi-6D "))
(cell_type_proportion_graph_6d_no_legend <- cell_type_proportion_graph_6d + theme(legend.position = "none"))

if (save_images) {
    save_tiff(cell_type_proportion_graph_6d_no_legend, "fig_6e.tiff",
              width = 6, height = 6)
}

cell_type_df_12d <- data.frame(CellType = d12_cells@meta.data$cell_type_pred,
                           Class = d12_cells@meta.data$sample_class)
cell_type_df_12d$Class <- factor(cell_type_df_12d$Class,
                                 levels = c("Control-6D", "ROCKi-6D",
                                            "Control-12D", "ROCKi-12D"))
cell_type_proportion_graph_12d <- create_cell_type_proportion_graph(cell_type_df_12d, c("Control-6D 6D ", "ROCKi-6D 6D "))
(cell_type_proportion_graph_12d_no_legend <- cell_type_proportion_graph_12d + theme(legend.position = "none"))
if (save_images) {
    save_tiff(cell_type_proportion_graph_12d_no_legend, "fig_6f.tiff",
              width = 6, height = 6)
}
prop_table_cell_type <- function(x) {prop.table(table(factor(x, levels =
                                                             c("Basal", "TA",
                                                               "TD",
                                                               "unknown"))))}

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

obs_total_prop_diff_6d <- calc_cell_type_props(d6_cells@meta.data$cell_type_pred,
                                               d6_cells@meta.data$treatment_class)
obs_total_prop_diff_12d <- calc_cell_type_props(d12_cells@meta.data$cell_type_pred,
                                                d12_cells@meta.data$treatment_class)

prop_diffs_6d <- c() # Store results
prop_diffs_12d <- c() # Store results

for (r in 1:10000) {
    # Generate resample indices
    c6d_idx <- sample(which(d6_cells@meta.data$treatment_class == "Control"),
                      length(which(d6_cells@meta.data$treatment_class == "Control")),
                      replace = TRUE)
    r6d_idx <- sample(which(d6_cells@meta.data$treatment_class == "ROCKi"),
                      length(which(d6_cells@meta.data$treatment_class == "ROCKi"))-1,
                      replace = TRUE)
    d6_idx <- c(c6d_idx, r6d_idx)

    c12d_idx <- sample(which(d12_cells@meta.data$treatment_class == "Control"),
                      length(which(d12_cells@meta.data$treatment_class == "Control")),
                      replace = TRUE)
    r12d_idx <- sample(which(d12_cells@meta.data$treatment_class == "ROCKi"),
                      length(which(d12_cells@meta.data$treatment_class == "ROCKi"))-1,
                      replace = TRUE)
    d12_idx <- c(c12d_idx, r12d_idx)

    # Find cell type and treatment class labels for this resample
    d6_cell_type <- d6_cells@meta.data$cell_type_pred[d6_idx]
    d6_treat_class <- d6_cells@meta.data$treatment_class[d6_idx]

    d12_cell_type <- d12_cells@meta.data$cell_type_pred[d12_idx]
    d12_treat_class <- d12_cells@meta.data$treatment_class[d12_idx]

    total_prop_diff <- calc_cell_type_props(d6_cell_type, d6_treat_class)
    prop_diffs_6d <- c(prop_diffs_6d, total_prop_diff)

    total_prop_diff <- calc_cell_type_props(d12_cell_type, d12_treat_class)
    prop_diffs_12d <- c(prop_diffs_12d, total_prop_diff)
}
alpha_vals <- c(0.975, 0.025)
residuals_6d <- prop_diffs_6d - obs_total_prop_diff_6d
print("95% confidence interval for total absolute difference in celltype proportions at 6D:")
print(as.numeric(obs_total_prop_diff_6d - quantile(residuals_6d, alpha_vals)))

residuals_12d <- prop_diffs_12d - obs_total_prop_diff_12d
print("95% confidence interval for total absolute difference in celltype proportions at 12D:")
print(as.numeric(obs_total_prop_diff_12d - quantile(residuals_12d, alpha_vals)))
signatures <- preprocess.signatures("celltype_markers_basal_ta_td.csv")

create_dotplot <- function(cells) {
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

skin_data@meta.data$sample_class_space <- skin_data@meta.data$sample_class
levels(skin_data@meta.data$sample_class_space) <- c("Control-6D ",
                                                    "ROCKi-6D ",
                                                    "Control-6D 6D ",
                                                    "ROCKi-6D 6D ")

holo_cells <- subset(skin_data, subset = cell_type_pred == "Basal")
ta_cells <- subset(skin_data, subset = cell_type_pred == "TA")
td_cells <- subset(skin_data, subset = cell_type_pred == "TD")

(holo_dp <- create_dotplot(holo_cells))
(ta_dp <- create_dotplot(ta_cells))
(td_dp <- create_dotplot(td_cells))

if (save_images) {
    save_tiff(holo_dp, "fig_6g.tiff", width = 20, height = 8)
    save_tiff(ta_dp, "fig_s2a.tiff", width = 20, height = 8)
    save_tiff(td_dp, "fig_s2b.tiff", width = 20, height = 8)
}
ctrl6d_cells <- subset(d6_cells, treatment_class == "Control")
rocki6d_cells <- subset(d6_cells, treatment_class == "ROCKi")
ctrl12d_cells <- subset(d12_cells, treatment_class == "Control")
rocki12d_cells <- subset(d12_cells, treatment_class == "ROCKi")

ctrl6d_cells <- RunUMAP(ctrl6d_cells, dims = 1:50, reduction = "harmony")
rocki6d_cells <- RunUMAP(rocki6d_cells, dims = 1:50, reduction = "harmony")
ctrl12d_cells <- RunUMAP(ctrl12d_cells, dims = 1:50, reduction = "harmony")
rocki12d_cells <- RunUMAP(rocki12d_cells, dims = 1:50, reduction = "harmony")


prepare_cells_for_slingshot <- function(cells) {
    cells <- FindNeighbors(cells, reduction = "harmony", dims = 1:50)
    cells <- FindClusters(cells, resolution = 0.1)
    cells <- RenameIdents(object = cells, "0" = "TA", "1" = "TD", "2" = "Basal")
    
    return (cells)
}

d6_cells <- prepare_cells_for_slingshot(d6_cells)
d12_cells <- prepare_cells_for_slingshot(d12_cells)
ctrl6d_cells <- prepare_cells_for_slingshot(ctrl6d_cells)
rocki6d_cells <- prepare_cells_for_slingshot(rocki6d_cells)
ctrl12d_cells <- prepare_cells_for_slingshot(ctrl12d_cells)
rocki12d_cells <- prepare_cells_for_slingshot(rocki12d_cells)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)

d6_cells.sce <- as.SingleCellExperiment(d6_cells)
d12_cells.sce <- as.SingleCellExperiment(d12_cells)
ctrl6d_cells.sce <- as.SingleCellExperiment(ctrl6d_cells)
rocki6d_cells.sce <- as.SingleCellExperiment(rocki6d_cells)
ctrl12d_cells.sce <- as.SingleCellExperiment(ctrl12d_cells)
rocki12d_cells.sce <- as.SingleCellExperiment(rocki12d_cells)

# Change UMAP values to make figures all similarly oriented
reducedDim(rocki6d_cells.sce, "UMAP")[,2] <- (-1) * reducedDims(rocki6d_cells.sce)$UMAP[,2]
reducedDim(ctrl12d_cells.sce, "UMAP")[,1] <- (-1) * reducedDims(ctrl12d_cells.sce)$UMAP[,1]

sshot_6d <- slingshot(d6_cells.sce, d6_cells$seurat_clusters,
                      reducedDim = "HARMONY", start.clus = "2")
sshot_12d <- slingshot(d12_cells.sce, d12_cells$seurat_clusters,
                       reducedDim = "HARMONY", start.clus = "2")

sshot_ctrl6d_separate <- slingshot(ctrl6d_cells.sce,
                                   ctrl6d_cells$seurat_clusters,
                                   reducedDim = "UMAP", start.clus = "2")
sshot_rocki6d_separate <- slingshot(rocki6d_cells.sce,
                                    rocki6d_cells$seurat_clusters,
                                    reducedDim = "UMAP", start.clus = "2")
sshot_ctrl12d_separate <- slingshot(ctrl12d_cells.sce,
                                    ctrl12d_cells$seurat_clusters,
                                    reducedDim = "UMAP", start.clus = "2")
sshot_rocki12d_separate <- slingshot(rocki12d_cells.sce,
                                     rocki12d_cells$seurat_clusters,
                                     reducedDim = "UMAP", start.clus = "2")

umap_curve_6d <- slingCurves(embedCurves(sshot_6d,
                                         reducedDims(d6_cells.sce)$UMAP))[[1]]
umap_curve_12d <- slingCurves(embedCurves(sshot_12d,
                                          reducedDims(d12_cells.sce)$UMAP))[[1]]

umap_curve_ctrl6d <- slingCurves(sshot_ctrl6d_separate)[[1]]
umap_curve_rocki6d <- slingCurves(sshot_rocki6d_separate)[[1]]
umap_curve_ctrl12d <- slingCurves(sshot_ctrl12d_separate)[[1]]
umap_curve_rocki12d <- slingCurves(sshot_rocki12d_separate)[[1]]

ctrl6d_idx <- which(d6_cells@meta.data$treatment_class == "Control")
rocki6d_idx <- which(d6_cells@meta.data$treatment_class == "ROCKi")
ctrl12d_idx <- which(d12_cells@meta.data$treatment_class == "Control")
rocki12d_idx <- which(d12_cells@meta.data$treatment_class == "ROCKi")

sshot_ctrl6d_together <- sshot_6d[,ctrl6d_idx]
sshot_rocki6d_together <- sshot_6d[,rocki6d_idx]
sshot_ctrl12d_together <- sshot_12d[,ctrl12d_idx]
sshot_rocki12d_together <- sshot_12d[,rocki12d_idx]

xmin_sep <- min(reducedDim(ctrl6d_cells.sce, "UMAP")[,1],
                reducedDim(rocki6d_cells.sce, "UMAP")[,1],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,1],
                reducedDim(rocki12d_cells.sce, "UMAP")[,1])
xmax_sep <- max(reducedDim(ctrl6d_cells.sce, "UMAP")[,1],
                reducedDim(rocki6d_cells.sce, "UMAP")[,1],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,1],
                reducedDim(rocki12d_cells.sce, "UMAP")[,1])
ymin_sep <- min(reducedDim(ctrl6d_cells.sce, "UMAP")[,2],
                reducedDim(rocki6d_cells.sce, "UMAP")[,2],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,2],
                reducedDim(rocki12d_cells.sce, "UMAP")[,2])
ymax_sep <- max(reducedDim(ctrl6d_cells.sce, "UMAP")[,2],
                reducedDim(rocki6d_cells.sce, "UMAP")[,2],
                reducedDim(ctrl12d_cells.sce, "UMAP")[,2],
                reducedDim(rocki12d_cells.sce, "UMAP")[,2])

create_traj_umap <- function(sshot_to_plot, umap_curve, cell_types) {
    points_df <- data.frame(u1 = reducedDims(sshot_to_plot)$UMAP[,1],
                            u2 = reducedDims(sshot_to_plot)$UMAP[,2],
                            cell_types=cell_types)
    points_df$cell_types <- factor(points_df$cell_types,
                                   levels = c("Basal", "TA", "TD", "unknown"))
    curve_df <- data.frame(u1 = umap_curve$s[umap_curve$ord, ][,1],
                           u2 = umap_curve$s[umap_curve$ord, ][,2])

    arrowhead_length <- unit(0.2, "cm")
    arrow_width <- 1

    ptimes <- sshot_to_plot$slingPseudotime_1
    basal_ptimes <- ptimes[which(cell_types == "Basal")]
    ta_ptimes <- ptimes[which(cell_types == "TA")]
    td_ptimes <- ptimes[which(cell_types == "TD")]

    percentile <- 0.97
    basal_percentile_ptime_idx <- which(abs(quantile(basal_ptimes, percentile) - basal_ptimes) ==
                                        min(abs(quantile(basal_ptimes, percentile) - basal_ptimes)))
    ta_percentile_ptime_idx <- which(abs(quantile(ta_ptimes, percentile) - ta_ptimes) ==
                                     min(abs(quantile(ta_ptimes, percentile) - ta_ptimes)))
    td_percentile_ptime_idx <- which(abs(quantile(td_ptimes, percentile) - td_ptimes) ==
                                     min(abs(quantile(td_ptimes, percentile) - td_ptimes)))

    basal_percentile_ptime <- basal_ptimes[basal_percentile_ptime_idx]
    ta_percentile_ptime <- ta_ptimes[ta_percentile_ptime_idx]
    td_percentile_ptime <- td_ptimes[td_percentile_ptime_idx]

    p <- ggplot() +
         geom_point(data = points_df, aes(x = u1, y = u2, color = cell_types),
                    size = -0.1) +
         scale_color_manual(labels=c("Holoclone-forming",
                           "Mero- or Paraclone-forming",
                           "Differentiated", "Unknown"),
                           values = c("#80C9EA", "#DD6E79",
                                      "#43863E", "#BBBBBB")) +
        geom_path(data = curve_df[1:round(length(sshot_to_plot$slingPseudotime_1) * percentile), ],
                  aes(x = u1, y = u2),
                  arrow = arrow(type = "closed", length = arrowhead_length),
                  size = arrow_width) +
         geom_path(data = curve_df[1:2,],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "open", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[1:which(sort(ptimes) == basal_percentile_ptime), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         geom_path(data = curve_df[1:which(sort(ptimes) == ta_percentile_ptime), ],
                   aes(x = u1, y = u2),
                   arrow = arrow(type = "closed", length = arrowhead_length),
                   size = arrow_width) +
         coord_fixed() + umap_theme

     return(p)
}

ctrl6d_traj_umap <- create_traj_umap(sshot_ctrl6d_separate, umap_curve_ctrl6d,
                                     Idents(ctrl6d_cells))
rocki6d_traj_umap <- create_traj_umap(sshot_rocki6d_separate, umap_curve_rocki6d,
                                      Idents(rocki6d_cells))
ctrl12d_traj_umap <- create_traj_umap(sshot_ctrl12d_separate, umap_curve_ctrl12d,
                                      Idents(ctrl12d_cells))
rocki12d_traj_umap <- create_traj_umap(sshot_rocki12d_separate, umap_curve_rocki12d,
                                       Idents(rocki12d_cells))

(ctrl6d_traj_umap_no_legend <- ctrl6d_traj_umap + theme(legend.position = "none"))
(rocki6d_traj_umap_no_legend <- rocki6d_traj_umap + theme(legend.position = "none"))
(ctrl12d_traj_umap_no_legend <- ctrl12d_traj_umap + theme(legend.position = "none"))
(rocki12d_traj_umap_no_legend <- rocki12d_traj_umap + theme(legend.position = "none"))

if (save_images) {
    save_tiff(ctrl6d_traj_umap_no_legend, "fig_7a.tiff", width = 10, height = 8)
    save_tiff(rocki6d_traj_umap_no_legend, "fig_7b.tiff", width = 10, height = 8)
    save_tiff(ctrl12d_traj_umap_no_legend, "fig_7c.tiff", width = 10, height = 8)
    save_tiff(rocki12d_traj_umap_no_legend, "fig_7d.tiff", width = 10, height = 8)
}
d6_cells@meta.data$cell_type <- Idents(d6_cells)
d12_cells@meta.data$cell_type <- Idents(d12_cells)

# mov_avg function copied from stackoverflow.com/a/4862334/6914552
mov_avg <- function(x, n = 151){stats::filter(x, rep(1 / n, n), sides = 2)}

create_traj_df <- function (ptimes, name, cell_types) {
    quantile_points <- seq(0.001,1,0.001)
    quantiled_data <- as.numeric(quantile(ptimes, quantile_points))
    sorted_cell_types <- cell_types[order(ptimes)]
    num_points <- length(ptimes)
    cell_indices <- round(num_points * quantile_points)
    cell_types_subset <- sorted_cell_types[cell_indices]

    return(data.frame(pct = quantile_points,
                      ptime = quantiled_data,
                      cell_type = factor(cell_types_subset,
                                         levels = c("Basal", "TA", "TD", "unknown")),
                      diff = c(NA, mov_avg(diff(quantiled_data))),
                      sample_class = rep(name, length(quantile_points))))
}

traj_df_ctrl6d <- create_traj_df(sshot_ctrl6d_together$slingPseudotime_1, "Control-6D",
                                 d6_cells[,ctrl6d_idx]@meta.data$cell_type)
traj_df_rocki6d <- create_traj_df(sshot_rocki6d_together$slingPseudotime_1, "ROCKi-6D",
                                  d6_cells[,rocki6d_idx]@meta.data$cell_type)
traj_df_6d <- rbind(traj_df_ctrl6d, traj_df_rocki6d)

traj_df_ctrl12d <- create_traj_df(sshot_ctrl12d_together$slingPseudotime_1,
                                  "Control-12D",
                                  d12_cells[,ctrl12d_idx]@meta.data$cell_type)
traj_df_rocki12d <- create_traj_df(sshot_rocki12d_together$slingPseudotime_1,
                                   "ROCKi-12D",
                                   d12_cells[,rocki12d_idx]@meta.data$cell_type)
traj_df_12d <- rbind(traj_df_ctrl12d, traj_df_rocki12d)

obs_result_ptime_6d <- mean(abs(traj_df_ctrl6d$ptime - traj_df_rocki6d$ptime))
obs_result_ptime_12d <- mean(abs(traj_df_ctrl12d$ptime - traj_df_rocki12d$ptime))

obs_result_speed_6d <- mean(abs(traj_df_ctrl6d$diff - traj_df_rocki6d$diff),
                            na.rm = TRUE)
obs_result_speed_12d <- mean(abs(traj_df_ctrl12d$diff - traj_df_rocki12d$diff),
                             na.rm = TRUE)

traj_plot_margin_theme <- theme(plot.margin = margin(0.1, 1.2, 0.1, 1.2, "cm"))
diff_comparison_plots_theme <- theme(legend.title = element_blank())

point_size <- 1
create_traj_plot <- function(traj_df) {
    shade_alpha <- 0.2

    basal_upper_lim <- quantile(traj_df[traj_df$cell_type == "Basal",]$ptime, 0.75)
    ta_upper_lim <- quantile(traj_df[traj_df$cell_type == "TA",]$ptime, 0.95)

    basal_color <- "#80C9EA"
    ta_color <- "#DD6E79"
    td_color <- "#43863E"

    traj_plot <- ggplot(traj_df, aes(x = pct, y = ptime, group = sample_class)) +
                    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
                             ymax = basal_upper_lim, alpha = shade_alpha,
                             fill = basal_color) +
                    annotate("rect", xmin = -Inf, xmax = Inf,
                             ymin = basal_upper_lim, ymax = ta_upper_lim,
                             alpha = shade_alpha, fill = ta_color) +
                    annotate("rect", xmin = -Inf, xmax = Inf,
                             ymin = ta_upper_lim, ymax = Inf,
                             alpha = shade_alpha, fill = td_color) +
                    geom_point(aes(color = cell_type), size = point_size) +
                    scale_color_manual(values = c(basal_color, ta_color, td_color)) +
                    new_scale_color() +
                    geom_line(aes(linetype = sample_class), size = 0.4) +
                    scale_linetype_manual(values = c("solid", "dashed")) +
                    scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
                    ylim(c(0, max(traj_df$ptime))) +
                    xlab("Percentage through differentation process") +
                    ylab("Pseudotime") + traj_plot_margin_theme +
                    diff_comparison_plots_theme + axis_text_theme

    return(traj_plot)
}

traj_plot_6d <- create_traj_plot(traj_df_6d)
traj_plot_12d <- create_traj_plot(traj_df_12d)

(traj_plot_6d_no_legend <- traj_plot_6d + theme(legend.position = "none"))
(traj_plot_12d_no_legend <- traj_plot_12d + theme(legend.position = "none"))

if (save_images) {
    save_tiff(traj_plot_6d_no_legend, "fig_7e.tiff", width = 10, height = 8)
    save_tiff(traj_plot_12d_no_legend, "fig_7f.tiff", width = 10, height = 8)
}


create_speed_plot <- function(traj_df) {
    speed_plot <- ggplot(traj_df, aes(x = pct, y = diff, group = sample_class)) +
                  geom_point(aes(color = cell_type), size = point_size) +
                  scale_color_manual(values = c("#80C9EA", "#DD6E79", "#43863E")) +
                  new_scale_color() +
                  geom_line(aes(linetype = sample_class), size = 0.4) +
                  scale_linetype_manual(values = c("solid", "dashed")) +
                  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
                  xlab("Percentage through differentation process") +
                  ylab("Differentiation speed") + traj_plot_margin_theme +
                  diff_comparison_plots_theme + axis_text_theme

    return(speed_plot)
}


(speed_plot_6d <- create_speed_plot(traj_df_6d) + theme(legend.position = "none"))
(speed_plot_12d <- create_speed_plot(traj_df_12d) + theme(legend.position = "none"))

if (save_images) {
    save_tiff(speed_plot_6d, "fig_s3a.tiff", width = 14, height = 8)
    save_tiff(speed_plot_12d, "fig_s3b.tiff", width = 14, height = 8)
}


for_dot_colour_legend <- ggplot(traj_df_6d,
                                aes(x = pct, y = ptime, group = sample_class)) +
                         geom_point(aes(color = cell_type), size = 3) +
                         scale_color_manual(values = c("#80C9EA", "#DD6E79", "#43863E"),
                                            labels = c("Holoclone-forming",
                                                       "Mero- or Paraclone-forming",
                                                       "Differentiated")) +
                         theme(legend.title = element_blank(),
                               legend.text = element_text(face = "bold"))
dot_colour_legend <- as_ggplot(get_legend(for_dot_colour_legend, "right"))

for_linetype_legend <- ggplot(traj_df_6d,
                                 aes(x = pct, y = ptime, group = sample_class)) +
                          geom_line(aes(linetype = sample_class), size = 0.4) +
                          scale_linetype_manual(values = c("solid", "dashed"),
                                                labels = c("Control", "ROCKi-treated")) +
                          theme(legend.title = element_blank(),
                                legend.text = element_text(face = "bold"))
linetype_legend <- as_ggplot(get_legend(for_linetype_legend, "right"))

(legends <- ggarrange(linetype_legend, dot_colour_legend, ncol = 2))
if (save_images) {
    save_tiff(legends, "fig_7_legends.tiff", width = 10, height = 8)
}
# Generate bootstrapped data as in "experiment_code_trajectories.R" and save
# the outputs as "trajectory_comparison_results.txt"

bootstrap_results <- read.table("trajectory_comparison_results.txt")
bootstrap_results <- bootstrap_results[, c(3, 4, 5, 6)]
colnames(bootstrap_results) <- c("TrajDiff6D", "TrajDiff12D", "SpeedDiff6D",
                                 "SpeedDiff12D")

diffs_quants_traj6d <- quantile((bootstrap_results$TrajDiff6D - obs_result_ptime_6d),
                                alpha_vals)
residuals_traj6d <- as.numeric(obs_result_ptime_6d - diffs_quants_traj6d)
print(paste("95% confidence interval for difference in pseudotimes at 6D:",
            paste(residuals_traj6d, collapse = " ")))

diffs_quants_traj12d <- quantile((bootstrap_results$TrajDiff12D - obs_result_ptime_12d),
                                alpha_vals)
residuals_traj12d <- as.numeric(obs_result_ptime_12d - diffs_quants_traj12d)
print(paste("95% confidence interval for difference in pseudotimes at 12D:",
            paste(residuals_traj12d, collapse = " ")))

diffs_quants_speed6d <- quantile((bootstrap_results$SpeedDiff6D - obs_result_speed_6d),
                                  alpha_vals)
residuals_speed6d <- as.numeric(obs_result_speed_6d - diffs_quants_speed6d)
print(paste("95% confidence interval for difference in speed at 6D:",
            paste(residuals_speed6d, collapse = " ")))

diffs_quants_speed12d <- quantile((bootstrap_results$SpeedDiff12D - obs_result_speed_12d),
                                  alpha_vals)
residuals_speed12d <- as.numeric(obs_result_speed_12d - diffs_quants_speed12d)
print(paste("95% confidence interval for difference in speed at 12D:",
            paste(residuals_speed12d, collapse = " ")))
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
