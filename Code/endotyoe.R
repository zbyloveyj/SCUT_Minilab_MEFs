packages_needed <- c("vegan", "cluster", "factoextra", "ggplot2", "philentropy", "openxlsx", 
                     "topicmodels", "Matrix", "BiotypeR", "reshape2", "pheatmap")
to_install <- packages_needed[!(packages_needed %in% installed.packages()[, "Package"])]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(packages_needed, library, character.only = TRUE))

microbiota <- read.xlsx("C:/rwork/Factor/temp-data/filter-final.xlsx", rowNames = TRUE)
data <- t(microbiota)
data <- data / rowSums(data)

jsd_matrix <- distance(as.matrix(data) + 1e-12, method = "jensen-shannon", test.na = FALSE)
jsd_dist <- as.dist(jsd_matrix)

set.seed(1001)
gap_result_pam <- clusGap(as.matrix(jsd_matrix), FUN = pam, K.max = 10, B = 50)
optimal_k_pam <- which.max(gap_result_pam$Tab[, "gap"])
pam_model <- pam(jsd_dist, k = optimal_k_pam, diss = TRUE)
pam_label <- pam_model$clustering
ordination <- cmdscale(jsd_dist, eig = TRUE, k = 2)
ordination_df <- as.data.frame(ordination$points)
colnames(ordination_df) <- c("Axis1", "Axis2")
ordination_df$cluster <- as.factor(pam_label)
ggplot(ordination_df, aes(x = Axis1, y = Axis2, color = cluster)) + geom_point(size = 3) + theme_minimal()

doc_matrix <- as(as.matrix(data), "CsparseMatrix")
lda_k_range <- 2:8
lda_perp_vec <- c()
for (k_val in lda_k_range) {
  lda_fit <- LDA(doc_matrix, k = k_val, method = "Gibbs", control = list(seed = 999))
  perp_val <- perplexity(lda_fit, doc_matrix)
  lda_perp_vec <- c(lda_perp_vec, perp_val)
}
best_k_lda <- lda_k_range[which.min(lda_perp_vec)]
lda_final <- LDA(doc_matrix, k = best_k_lda, method = "Gibbs", control = list(seed = 999))
lda_topics <- posterior(lda_final)$topics
lda_assign <- apply(lda_topics, 1, which.max)
lda_label <- factor(lda_assign)
top_lda_terms <- terms(lda_final, 15)

res_biotype <- biotyper.data.frame(microbiota, k = NULL)
label_biotype <- res_biotype$cluster

transposed_data <- t(microbiota)
transposed_scaled <- scale(transposed_data)
gap_result_kmeans <- clusGap(transposed_scaled, FUN = kmeans, K.max = 10, B = 50)
best_k_kmeans <- which.max(gap_result_kmeans$Tab[, "gap"])
kmeans_fit <- kmeans(transposed_scaled, centers = best_k_kmeans, nstart = 50)
kmeans_label <- kmeans_fit$cluster

cluster_labels <- data.frame(
  Sample = rownames(data),
  PAM = as.factor(pam_label),
  LDA = as.factor(lda_label),
  Biotype = as.factor(label_biotype),
  KMeans = as.factor(kmeans_label)
)

p1 <- fviz_cluster(list(data = transposed_scaled, cluster = kmeans_label), 
                   ellipse.type = "convex", palette = "jco", ggtheme = theme_minimal())
p2 <- fviz_cluster(list(data = transposed_scaled, cluster = pam_label), 
                   ellipse.type = "convex", palette = "npg", ggtheme = theme_light())
p3 <- fviz_cluster(list(data = transposed_scaled, cluster = lda_assign), 
                   ellipse.type = "t", palette = "aaas", ggtheme = theme_classic())
p4 <- fviz_cluster(list(data = transposed_scaled, cluster = label_biotype), 
                   ellipse.type = "norm", palette = "d3", ggtheme = theme_bw())
pheatmap(log1p(microbiota), cluster_rows = TRUE, cluster_cols = TRUE, 
         annotation_col = data.frame(Cluster = factor(pam_label)))

mat_plot <- cmdscale(dist(scale(data)), k = 2)
df_plot <- data.frame(mat_plot, cluster = factor(pam_label))
ggplot(df_plot, aes(X1, X2, color = cluster)) + geom_point(size = 2.8) + 
  stat_ellipse(type = "t") + scale_color_viridis_d() + theme_minimal()

meta_kmeans <- silhouette(kmeans_label, dist(transposed_scaled))
meta_pam <- silhouette(pam_label, jsd_dist)
meta_lda <- silhouette(lda_assign, dist(transposed_scaled))
meta_bio <- silhouette(label_biotype, dist(transposed_scaled))

ari_mat <- matrix(NA, 4, 4)
methods <- list(PAM = pam_label, LDA = lda_assign, KMeans = kmeans_label, Biotype = label_biotype)
for (i in 1:4) {
  for (j in 1:4) {
    ari_mat[i, j] <- mclust::adjustedRandIndex(methods[[i]], methods[[j]])
  }
}
colnames(ari_mat) <- rownames(ari_mat) <- names(methods)
print(ari_mat)

write.xlsx(cluster_labels, "C:/rwork/Factor/enterotype_labels_all_methods.xlsx", rowNames = FALSE)
write.xlsx(data.frame(ARI = ari_mat), "C:/rwork/Factor/ARI_matrix.xlsx", rowNames = TRUE)
