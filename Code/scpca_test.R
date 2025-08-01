# ====== Setup ======
library(here)
library(tidyverse)
library(scPCA)
library(Rtsne)
library(umap)
library(ggpubr)
library(GEOquery)
library(genefilter)
library(cluster)
library(xtable)

set.seed(9765967)

# ====== Load and Prepare Data ======
ges <- getGEO("GSE51808")$GSE51808_series_matrix.txt.gz
ggene_names <- read_tsv(here("C:/download/GPL13158-5065.txt"), skip = 16)
var_filt_ges <- varFilter(ges, var.cutoff = 1 - 500 / nrow(exprs(ges)))
control_label <- which(var_filt_ges$`status:ch1` == "control")
target <- t(exprs(var_filt_ges)[, -control_label])
background <- t(exprs(var_filt_ges)[, control_label])
dengue_class <- var_filt_ges$`status:ch1`[-control_label]

# ====== PCA ======
dengue_pca <- prcomp(target, center = TRUE, scale. = FALSE)
group_mem <- if_else(dengue_class == "convalescent", 1, if_else(dengue_class == "DF", 2, 3))
pca_group_sil <- silhouette(group_mem, dist(dengue_pca$x))
pca_df <- data.frame(PC1 = dengue_pca$x[, 1], PC2 = dengue_pca$x[, 2], class = dengue_class) %>%
  mutate(class = if_else(class == "convalescent",
                         paste0("Conv. (", sprintf("%.3f", round(summary(pca_group_sil)$clus.avg.widths[1], 3)), ")"),
                         if_else(class == "DF",
                                 paste0("DF (", sprintf("%.3f", round(summary(pca_group_sil)$clus.avg.widths[2], 3)), ")"),
                                 paste0("DHF (", sprintf("%.3f", round(summary(pca_group_sil)$clus.avg.widths[3], 3)), ")")
                         )))
pca_p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = class)) +
  geom_point() +
  ggtitle("PCA") +
  scale_colour_viridis_d(name = "Class", begin = 0.1, end = 0.9) +
  theme_minimal()

# ====== cPCA ======
dengue_cpca <- scPCA(target, background, center = TRUE, penalties = 0, n_centers = 3, max_iter = 1000)
cpca_group_sil <- silhouette(group_mem, dist(dengue_cpca$x))
cpca_df <- data.frame(cPC1 = -dengue_cpca$x[, 1], cPC2 = dengue_cpca$x[, 2], class = dengue_class) %>%
  mutate(class = if_else(class == "convalescent",
                         paste0("Conv. (", sprintf("%.3f", round(summary(cpca_group_sil)$clus.avg.widths[1], 3)), ")"),
                         if_else(class == "DF",
                                 paste0("DF (", sprintf("%.3f", round(summary(cpca_group_sil)$clus.avg.widths[2], 3)), ")"),
                                 paste0("DHF (", sprintf("%.3f", round(summary(cpca_group_sil)$clus.avg.widths[3], 3)), ")")
                         )))
cpca_p <- ggplot(cpca_df, aes(x = cPC1, y = cPC2, colour = class)) +
  geom_point() +
  ggtitle("cPCA") +
  scale_colour_viridis_d(name = "Class", begin = 0.1, end = 0.9) +
  theme_minimal()

diff_centers_cpca <- lapply(2:5, function(x) {
  scPCA(target, background, center = TRUE, penalties = 0, n_centers = x, max_iter = 1000)
})

# ====== scPCA ======
load(file = here("analyses/dengue_data/cluster_files/scpca.Rdata"))
scpca_group_sil <- silhouette(group_mem, dist(dengue_scpca$x))
scpca_df <- data.frame(scPC1 = dengue_scpca$x[, 1], scPC2 = dengue_scpca$x[, 2], class = dengue_class) %>%
  mutate(class = if_else(class == "convalescent",
                         paste0("Conv. (", sprintf("%.3f", round(summary(scpca_group_sil)$clus.avg.widths[1], 3)), ")"),
                         if_else(class == "DF",
                                 paste0("DF (", sprintf("%.3f", round(summary(scpca_group_sil)$clus.avg.widths[2], 3)), ")"),
                                 paste0("DHF (", sprintf("%.3f", round(summary(scpca_group_sil)$clus.avg.widths[3], 3)), ")")
                         )))
scpca_p <- ggplot(scpca_df, aes(x = scPC1, y = scPC2, colour = class)) +
  geom_point() +
  ggtitle("scPCA") +
  scale_colour_viridis_d(name = "Class", begin = 0.1, end = 0.9) +
  theme_minimal()

# ====== t-SNE (with PCA) ======
dengue_tsne_pca <- Rtsne(target, perplexity = 8, pca = TRUE, max_iter = 1000, theta = 0)
tsne_pca_group_sil <- silhouette(group_mem, dist(dengue_tsne_pca$Y))
tsne_pca_df <- data.frame(TSNE1 = dengue_tsne_pca$Y[, 1], TSNE2 = dengue_tsne_pca$Y[, 2], class = dengue_class) %>%
  mutate(class = if_else(class == "convalescent",
                         paste0("Conv. (", sprintf("%.3f", round(summary(tsne_pca_group_sil)$clus.avg.widths[1], 3)), ")"),
                         if_else(class == "DF",
                                 paste0("DF (", sprintf("%.3f", round(summary(tsne_pca_group_sil)$clus.avg.widths[2], 3)), ")"),
                                 paste0("DHF (", sprintf("%.3f", round(summary(tsne_pca_group_sil)$clus.avg.widths[3], 3)), ")")
                         )))
tsne_pca_p <- ggplot(tsne_pca_df, aes(x = TSNE1, y = TSNE2, colour = class)) +
  geom_point() +
  ggtitle("t-SNE") +
  scale_colour_viridis_d(name = "Class", begin = 0.1, end = 0.9) +
  theme_minimal()

# ====== t-SNE (no PCA) ======
dengue_tsne <- Rtsne(target, perplexity = 8, pca = FALSE, theta = 0)
tsne_group_sil <- silhouette(group_mem, dist(dengue_tsne$Y))
tsne_df <- data.frame(TSNE1 = dengue_tsne$Y[, 1], TSNE2 = dengue_tsne$Y[, 2], class = dengue_class) %>%
  mutate(class = if_else(class == "convalescent",
                         paste0("Conv. (", sprintf("%.3f", round(summary(tsne_group_sil)$clus.avg.widths[1], 3)), ")"),
                         if_else(class == "DF",
                                 paste0("DF (", sprintf("%.3f", round(summary(tsne_group_sil)$clus.avg.widths[2], 3)), ")"),
                                 paste0("DHF (", sprintf("%.3f", round(summary(tsne_group_sil)$clus.avg.widths[3], 3)), ")")
                         )))
tsne_p <- ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, colour = class)) +
  geom_point() +
  ggtitle("t-SNE") +
  scale_colour_viridis_d(name = "Class", begin = 0.1, end = 0.9) +
  theme_minimal()

# ====== UMAP ======
params <- umap.defaults
params$n_neighbors <- 15
params$min_dist <- 0.2
dengue_umap <- umap(target, config = params)
umap_group_sil <- silhouette(group_mem, dist(dengue_umap$layout))
umap_df <- as.data.frame(dengue_umap$layout) %>%
  mutate(class = dengue_class) %>%
  mutate(class = if_else(class == "convalescent",
                         paste0("Conv. (", sprintf("%.3f", round(summary(umap_group_sil)$clus.avg.widths[1], 3)), ")"),
                         if_else(class == "DF",
                                 paste0("DF (", sprintf("%.3f", round(summary(umap_group_sil)$clus.avg.widths[2], 3)), ")"),
                                 paste0("DHF (", sprintf("%.3f", round(summary(umap_group_sil)$clus.avg.widths[3], 3)), ")")
                         )))
colnames(umap_df) <- c("UMAP1", "UMAP2", "class")
umap_p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = class)) +
  geom_point() +
  ggtitle("UMAP") +
  scale_colour_viridis_d(name = "Class", begin = 0.1, end = 0.9) +
  theme_minimal()

# ====== Save/Plot Results ======
ggsave("manuscript/figures/dengue_pca.png", pca_p, dpi = 300, width = 6, height = 4)
ggsave("manuscript/figures/dengue_cpca.png", cpca_p, dpi = 300, width = 6, height = 4)
ggsave("manuscript/figures/dengue_scpca.png", scpca_p, dpi = 300, width = 6, height = 4)
ggsave("manuscript/figures/dengue_tsne.png", tsne_p, dpi = 300, width = 6, height = 4)
ggsave("manuscript/figures/dengue_umap.png", umap_p, dpi = 300, width = 6, height = 4)
