library(foreach)
library(doParallel)

library(Seurat)
library(SingleCellExperiment)
library(slingshot)


quiet_sshot <- function(sce, clusts, reduced_dim, start_clus) {
    suppressMessages(slingshot(sce, clusts, reducedDim = reduced_dim,
                               start.clus = start_clus))
}


mov_avg <- function(x, n = 151) {
    stats::filter(x, rep(1 / n, n), sides = 2)
}


create_traj_df <- function(ptimes, name) {
    quantile_points <- seq(0.001, 1, 0.001)
    quantiled_data <- as.numeric(quantile(ptimes, quantile_points))
    num_points <- length(ptimes)
    cell_indices <- round(num_points * quantile_points)

    return(data.frame(pct = quantile_points,
                      ptime = quantiled_data,
                      diff = c(NA, mov_avg(diff(quantiled_data))),
                      sample_class = rep(name, length(quantile_points))))
}


prepare_cells_for_slingshot <- function(cells) {
    cells <- FindNeighbors(cells, reduction = "harmony", dims = 1:50)
    cells <- FindClusters(cells, resolution = 0.1)
    cells <- RenameIdents(object = cells, "0" = "TA", "1" = "TD", "2" = "Basal")

    return (cells)
}


do_traj_comp <- function(r, idx_list) {
    c6d_idx <- idx_list[[1]]
    r6d_idx <- idx_list[[2]]
    d6_idx <- idx_list[[3]]
    c12d_idx <- idx_list[[4]]
    r12d_idx <- idx_list[[5]]
    d12_idx <- idx_list[[6]]

    d6_harmony <- d6_cells@reductions$harmony@cell.embeddings[d6_idx, ]
    d6_sce <- SingleCellExperiment(t(d6_harmony))
    reducedDim(d6_sce, "HARMONY") <- d6_harmony
    sshot_d6 <- quiet_sshot(d6_sce, d6_cells$seurat_clusters[d6_idx],
                            "HARMONY", "2")

    d12_harmony <- d12_cells@reductions$harmony@cell.embeddings[d12_idx, ]
    d12_sce <- SingleCellExperiment(t(d12_harmony))
    reducedDim(d12_sce, "HARMONY") <- d12_harmony
    sshot_d12 <- quiet_sshot(d12_sce, d12_cells$seurat_clusters[d12_idx],
                                "HARMONY", "2")

    traj_df_ctrl6d <- create_traj_df(sshot_d6$slingPseudotime_1[1:length(c6d_idx)],
                                     "Control-6D")
    traj_df_rocki6d <- create_traj_df(sshot_d6$slingPseudotime_1[(length(c6d_idx) + 1):length(d6_idx)],
                                     "ROCKi-6D")

    traj_df_ctrl12d <- create_traj_df(sshot_d12$slingPseudotime_1[1:length(c12d_idx)],
                                     "Control-12D")
    traj_df_rocki12d <- create_traj_df(sshot_d12$slingPseudotime_1[(length(c12d_idx) + 1):length(d12_idx)],
                                      "ROCKi-12D")

    result_ptime_6d <- mean(abs(traj_df_ctrl6d$ptime - traj_df_rocki6d$ptime))
    result_ptime_12d <- mean(abs(traj_df_ctrl12d$ptime - traj_df_rocki12d$ptime))

    result_speed_6d <- mean(abs(traj_df_ctrl6d$diff - traj_df_rocki6d$diff),
                            na.rm = TRUE)
    result_speed_12d <- mean(abs(traj_df_ctrl12d$diff - traj_df_rocki12d$diff),
                             na.rm = TRUE)

    all_results <- c(result_ptime_6d, result_ptime_12d, result_speed_6d, result_speed_12d)
    print(c(r, all_results))


    return(all_results)
}


skin_data <- readRDS("skin_data_end_pipeline_1458110522.rds")

d6_cells <- subset(skin_data, subset = sample_class %in% c("Control-6D", "ROCKi-6D"))
d12_cells <- subset(skin_data, subset = sample_class %in% c("Control-12D", "ROCKi-12D"))

d6_cells <- prepare_cells_for_slingshot(d6_cells)
d12_cells <- prepare_cells_for_slingshot(d12_cells)

c6d_idx_list <- c()
r6d_idx_list <- c()
d6_idx_list <- c()

c12d_idx_list <- c()
r12d_idx_list <- c()
d12_idx_list <- c()

idx_lists <- c()

num_resamples <- 10000

for (r in 1:num_resamples) {
    c6d_idx <- sample(which(d6_cells@meta.data$treatment_class == "Control"),
                  length(which(d6_cells@meta.data$treatment_class == "Control")),
                  replace = TRUE)
    c6d_idx_list[[r]] <- c6d_idx

    r6d_idx <- sample(head(which(d6_cells@meta.data$treatment_class == "ROCKi"), -1),
                      length(which(d6_cells@meta.data$treatment_class == "ROCKi")) - 1,
                      replace = TRUE)
    r6d_idx_list[[r]] <- r6d_idx

    d6_idx <- c(c6d_idx, r6d_idx)
    d6_idx_list[[r]] <- d6_idx

    c12d_idx <- sample(which(d12_cells@meta.data$treatment_class == "Control"),
                       length(which(d12_cells@meta.data$treatment_class == "Control")),
                       replace = TRUE)
    c12d_idx_list[[r]] <- c12d_idx

    r12d_idx <- sample(head(which(d12_cells@meta.data$treatment_class == "ROCKi"), -1),
                       length(which(d12_cells@meta.data$treatment_class == "ROCKi")) - 1,
                       replace = TRUE)
    r12d_idx_list[[r]] <- r12d_idx

    d12_idx <- c(c12d_idx, r12d_idx)
    d12_idx_list[[r]] <- d12_idx

    idx_lists[[r]] <- list(c6d_idx, r6d_idx, d6_idx, c12d_idx, r12d_idx, d12_idx)
}

num_clusters <- 4
registerDoParallel(num_clusters)
sink("trajectory_comparison_results_parallelised.txt")

foreach (r = 1:num_resamples) %dopar% {
    do_traj_comp(r, idx_lists[[r]])
}

sink()
stopImplicitCluster()

# To generate confidence intervals (based on https://ocw.mit.edu/courses/mathematics/18-05-introduction-to-probability-and-statistics-spring-2014/readings/MIT18_05S14_Reading24.pdf):
# In Python:
# import sys
# for file_name in sys.argv[1:]:
#     with open(file_name) as f:
#         for line in f:
#             splat = line[1:-1].split()
#             rejoin = ",".join(splat[1:])
#             print(rejoin[1:-1])
# and save as "bootstrap_results.csv". Then:
# bootstrap_results <- read.csv("trajectory_comparison_results_all_fixed_5k_cut.csv", header = FALSE,
#                               col.names = c("ResampleIndex", "TrajDiff6D",
#                                             "TrajDiff12D", "SpeedDiff6D",
#                                             "SpeedDiff12D"))
#
# If the observed mean absolute difference of pseudotimes at 6D and 12D are
# result_ptime_6d and result_ptime_12d, resp., and if the mean absolute
# differences of speeds are result_speed_6d and result_speed_12d, resp, then a
# 99% confidence interval for each statistic can be computed as follows:
# as.numeric(result_ptime_6d - quantile((bootstrap_results$TrajDiff6D - result_ptime_6d), c(0.995, 0.005)))
# as.numeric(result_ptime_12d - quantile((bootstrap_results$TrajDiff12D - result_ptime_12d), c(0.995, 0.005)))
# as.numeric(result_speed_6d - quantile((bootstrap_results$SpeedDiff6D - result_speed_6d), c(0.995, 0.005)))
# as.numeric(result_speed_12d - quantile((bootstrap_results$SpeedDiff12D - result_speed_12d), c(0.995, 0.005)))
#
# Note that, based on the code in the Rmarkdown file for the data analysis, the observed results can be calculated as
# result_ptime_6d <- mean(abs(traj_df_ctrl6d$ptime - traj_df_rocki6d$ptime))
#   = 14.30379
# result_ptime_12d <- mean(abs(traj_df_ctrl12d$ptime - traj_df_rocki12d$ptime))
#   = 3.015428
# result_speed_6d <- mean(abs(traj_df_ctrl6d$diff - traj_df_rocki6d$diff), na.rm = TRUE)
#   = 0.0448669
# result_speed_12d <- mean(abs(traj_df_ctrl12d$diff - traj_df_rocki12d$diff), na.rm = TRUE)
#   = 0.02387632
#
# The generated results:
# > as.numeric(result_ptime_6d - quantile((bootstrap_results$TrajDiff6D - result_ptime_6d), c(0.995, 0.005)))
# [1] 11.98288 15.96313
# > as.numeric(result_ptime_12d - quantile((bootstrap_results$TrajDiff12D - result_ptime_12d), c(0.995, 0.005)))
# [1] 1.283433 4.269527
# > as.numeric(result_speed_6d - quantile((bootstrap_results$SpeedDiff6D - result_speed_6d), c(0.995, 0.005)))
# [1] 0.03394174 0.05393807
# > as.numeric(result_speed_12d - quantile((bootstrap_results$SpeedDiff12D - result_speed_12d), c(0.995, 0.005)))
# [1] 0.01244744 0.03357274
