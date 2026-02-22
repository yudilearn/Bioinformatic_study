library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)

#Analisis GO dan KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#Contoh Menganalisis data dari GEO

gset <- getGEO("GSE10072", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
ex <- exprs(gset)

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm =TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}
group_info <- pData(gset)[["source_name_ch1"]]
groups <- make.names(group_info)
gset$group <- factor(groups)
nama_grup <- levels(gset$group)
print(nama_grup)
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)
grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]
contrast_formula <- paste(grup_kanker, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)
probe_ids <- rownames(topTableResults)

gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#Set warna berdasarkan grup
group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )
#PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)

#UMAP digunakan untuk:
#- Mereduksi ribuan gen menjadi 2 dimensi
#- Melihat pemisahan sampel secara global
#- Alternatif PCA (lebih sensitif ke struktur lokal)

#Transpose matriks ekspresi:
#UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

#Jalankan UMAP
umap_result <- umap(umap_input)

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#vulcano plot
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Kanker Paru")

# 1. Ekstrak daftar gen yang UP-regulated
genes_up <- topTableResults %>%
  filter(logFC > 1 & adj.P.Val < 0.01) %>%
  arrange(desc(logFC)) # Urutkan dari yang kenaikannya paling tinggi

# 2. Ekstrak daftar gen yang DOWN-regulated
genes_down <- topTableResults %>%
  filter(logFC < -1 & adj.P.Val < 0.01) %>%
  arrange(logFC) # Urutkan dari yang penurunannya paling drastis

# 3. Melihat jumlah gen yang didapat
message("Jumlah gen yang Up-regulated: ", nrow(genes_up))
message("Jumlah gen yang Down-regulated: ", nrow(genes_down))

# 4. Menampilkan 10 gen teratas untuk masing-masing kategori
print("Top 10 Gen Up-regulated:")
head(genes_up[, c("SYMBOL", "logFC", "adj.P.Val")], 10)

print("Top 10 Gen Down-regulated:")
head(genes_down[, c("SYMBOL", "logFC", "adj.P.Val")], 10)

# 5. Simpan ke file terpisah jika diperlukan
write.csv(genes_up, "Daftar_Gen_UP.csv", row.names = FALSE)
write.csv(genes_down, "Daftar_Gen_DOWN.csv", row.names = FALSE)

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#additional
# 1. Buat data frame yang berisi info grup (Kanker/Normal)
annotation_col <- data.frame(Group = gset$group)

# 2. Berikan nama baris pada annotation_col agar sama dengan nama kolom di mat_heatmap
# Ini penting agar R tahu sampel mana masuk grup mana
rownames(annotation_col) <- colnames(mat_heatmap)

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi Heatmap
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)


# ==============================================================================
# ANALISIS FUNGSIONAL (UP & DOWN-REGULATED)
# ==============================================================================

# Analisis UP-REGULATED
gene_conv_up <- bitr(genes_up$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ego_up <- enrichGO(gene=gene_conv_up$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", readable=TRUE)
ekegg_up <- enrichKEGG(gene=gene_conv_up$ENTREZID, organism='hsa')

# Analisis DOWN-REGULATED
gene_conv_down <- bitr(genes_down$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
ego_down <- enrichGO(gene=gene_conv_down$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", readable=TRUE)
ekegg_down <- enrichKEGG(gene=gene_conv_down$ENTREZID, organism='hsa')

# ==============================================================================
# 9. VISUALISASI & PENYIMPANAN OTOMATIS (VERSI BERSIH)
# ==============================================================================

# --- PROSES UP-REGULATED ---
if(exists("ego_up") && nrow(as.data.frame(ego_up)) > 0) {
  p_go_up <- barplot(ego_up, showCategory=10) + ggtitle("GO Biological Process: UP")
  print(p_go_up)
  ggsave("Plot_GO_UP.png", plot = p_go_up, width = 8, height = 6, dpi = 300)
}

if(exists("ekegg_up") && nrow(as.data.frame(ekegg_up)) > 0) {
  p_kegg_up <- dotplot(ekegg_up, showCategory=10) + ggtitle("KEGG Pathways: UP")
  print(p_kegg_up)
  ggsave("Plot_KEGG_UP.png", plot = p_kegg_up, width = 8, height = 6, dpi = 300)
}

# --- PROSES DOWN-REGULATED ---
if(exists("ego_down") && nrow(as.data.frame(ego_down)) > 0) {
  p_go_down <- barplot(ego_down, showCategory=10) + ggtitle("GO Biological Process: DOWN")
  print(p_go_down)
  ggsave("Plot_GO_DOWN.png", plot = p_go_down, width = 8, height = 6, dpi = 300)
}

if(exists("ekegg_down") && nrow(as.data.frame(ekegg_down)) > 0) {
  p_kegg_down <- dotplot(ekegg_down, showCategory=10) + ggtitle("KEGG Pathways: DOWN")
  print(p_kegg_down)
  ggsave("Plot_KEGG_DOWN.png", plot = p_kegg_down, width = 8, height = 6, dpi = 300)
}

message("Proses selesai! Cek folder kerja Anda untuk melihat file .png")
  