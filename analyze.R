#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=5L) stop("Expected: manifest results_dir pathways.gmt CASE CONTROL")
manifest <- args[1]; out <- args[2]; gmt <- args[3]; case <- args[4]; control <- args[5]
needed <- c("tximport", "DESeq2", "fgsea", "ggplot2", "ggrepel", "SummarizedExperiment")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse=", "), ". Nothing was installed.")
writeLines(capture.output(sessionInfo()), file.path(out,"logs","R_sessionInfo.txt"))
samples <- read.delim(manifest, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(all(c("Sample","Condition") %in% names(samples)), !anyDuplicated(samples$Sample))
if (case==control || !setequal(unique(samples$Condition),c(case,control))) stop("Need exactly case/control conditions")
if (any(table(samples$Condition)<2)) stop("Need at least two biological replicates per group")
files <- setNames(file.path(out,"04_kallisto",samples$Sample,"abundance.tsv"),samples$Sample)
if (!all(file.exists(files))) stop("Missing kallisto abundance.tsv files")
tx2gene <- read.delim(file.path(out,"00_reference","tx2gene.tsv"), stringsAsFactors=FALSE)
# Strip only terminal numeric transcript versions consistently on both sides.
tx2gene$transcript_id <- sub("\\.[0-9]+$", "", tx2gene$transcript_id)
if (any(vapply(split(tx2gene$gene_id,tx2gene$transcript_id),function(x) length(unique(x))>1L,logical(1)))) stop("Transcript version removal creates conflicting gene mappings")
tx2gene <- tx2gene[!duplicated(tx2gene$transcript_id), ,drop=FALSE]
txi <- tximport::tximport(files, type="kallisto", tx2gene=tx2gene,
                         ignoreTxVersion=TRUE, countsFromAbundance="no", dropInfReps = TRUE)
# Use DESeqDataSetFromTximport below; it handles estimated counts and length offsets.
write.csv(data.frame(gene_id=rownames(txi$counts),txi$counts,check.names=FALSE),
          file.path(out,"04_kallisto","estimated_gene_counts.csv"),row.names=FALSE)
samples$Condition <- factor(samples$Condition, levels=c(control, case))
rownames(samples) <- samples$Sample
dds <- DESeq2::DESeqDataSetFromTximport(txi, colData=samples, design=~Condition)
dds <- dds[rowSums(DESeq2::counts(dds) >= 10) >= 2, ]
if (nrow(dds) == 0L) stop("No genes passed count filter; check strandedness, kallisto logs and transcript mapping")
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast=c("Condition", case, control))
resdf <- data.frame(gene_id=rownames(res), as.data.frame(res), check.names=FALSE)
resdf <- resdf[order(resdf$padj, na.last=TRUE), , drop=FALSE]

write.csv(data.frame(gene_id=rownames(dds), DESeq2::counts(dds, normalized=TRUE), check.names=FALSE),
          file.path(out,"05_deseq2","normalized_counts.csv"), row.names=FALSE)
saveRDS(dds, file.path(out,"05_deseq2","deseq2_object.rds"))
pdf(file.path(out,"05_deseq2","MA_plot.pdf")); DESeq2::plotMA(res, ylim=c(-5,5)); dev.off()
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE)
pca <- prcomp(t(SummarizedExperiment::assay(vsd)))
write.csv(data.frame(sample=samples$Sample, condition=samples$Condition, PC1=pca$x[,1], PC2=pca$x[,2]),
          file.path(out,"05_deseq2","PCA_coordinates.csv"), row.names=FALSE)
pdf(file.path(out,"05_deseq2","PCA.pdf"))
plot(pca$x[,1:2], col=as.integer(samples$Condition), pch=19, xlab="PC1", ylab="PC2")
text(pca$x[,1:2], labels=samples$Sample, pos=3, cex=.7)
legend("topright", legend=levels(samples$Condition), col=seq_along(levels(samples$Condition)), pch=19)
dev.off()

# Use symbols extracted from the same GTF as the transcript-to-gene map.
gene_table <- read.delim(file.path(out, "00_reference", "gene_symbols.tsv"),
                        check.names=FALSE, stringsAsFactors=FALSE)
gene_id <- gene_table$gene_id
gene_symbol <- gene_table$gene_name
gene_map <- setNames(gene_symbol, sub("\\.[0-9]+$", "", gene_id))
gene_map <- gene_map[!duplicated(names(gene_map))]
gene_names <- unname(gene_map[sub("\\.[0-9]+$", "", resdf$gene_id)])
missing_names <- is.na(gene_names) | !nzchar(gene_names)
gene_names[missing_names] <- resdf$gene_id[missing_names]
output_df <- resdf
output_df$gene_id <- gene_names
names(output_df)[names(output_df) == "gene_id"] <- "gene_name"
write.csv(output_df, file.path(out, "05_deseq2", "all_genes.csv"), row.names = FALSE)
write.csv(subset(output_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >=1),
          file.path(out, "05_deseq2", "significant_genes.csv"), row.names = FALSE)
# Sample correlations on variance-stabilized expression.
expr <- SummarizedExperiment::assay(vsd)
if (any(!is.finite(expr))) stop("Nonfinite transformed expression")
pdf(file.path(out,"05_deseq2","sample_correlation.pdf"), width=8,height=8)
heatmap(cor(expr), symm=TRUE, scale="none", margins=c(8,8),
        main="Sample correlation (VST)", col=colorRampPalette(c("white","steelblue"))(100))
dev.off()
# Top 30 significant genes by adjusted p-value, shown as row z-scores.
sel <- head(resdf$gene_id[!is.na(resdf$padj) & resdf$padj<0.05 & abs(resdf$log2FoldChange)>=1],30)
if (length(sel)>=2L) {
  mat <- expr[sel,,drop=FALSE]
  mat <- mat[apply(mat,1,sd)>0,,drop=FALSE]
  if (nrow(mat)>=2L) {
    symbols <- unname(gene_map[sub("\\.[0-9]+$","",rownames(mat))])
    bad <- is.na(symbols) | !nzchar(symbols)
    symbols[bad] <- rownames(mat)[bad]
    rownames(mat) <- make.unique(symbols)
    pdf(file.path(out,"05_deseq2","top_genes_heatmap.pdf"),width=9,height=10)
    heatmap(mat,scale="row",margins=c(8,8),
            col=colorRampPalette(c("#4575B4","white","#D73027"))(101),
            main="Top significant genes: row-scaled VST")
    dev.off()
  }
}
write.csv(data.frame(gene_id=resdf$gene_id,gene_name=gene_names,
                    mapped_symbol=!missing_names),
          file.path(out,"05_deseq2","gene_name_mapping.csv"),row.names=FALSE)
# ---- Volcano plot ----

# Use the results containing gene names.
volcano_df <- output_df

# Keep rows with usable fold changes and adjusted p-values.
volcano_df <- volcano_df[
  is.finite(volcano_df$log2FoldChange) &
    is.finite(volcano_df$padj) &
    volcano_df$padj >= 0,
  ,
  drop = FALSE
]

if (nrow(volcano_df) == 0) {
  stop("No genes have usable values for the volcano plot")
}

# Calculate the y-axis; avoid infinity when padj is zero.
volcano_df$neg_log10_padj <- -log10(
  pmax(volcano_df$padj, .Machine$double.xmin)
)

# Classify genes using your existing thresholds.
volcano_df$direction <- "Below thresholds"

volcano_df$direction[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange >= 1
] <- "Upregulated"

volcano_df$direction[
  volcano_df$padj < 0.05 &
    volcano_df$log2FoldChange <= -1
] <- "Downregulated"

# Find the top 5 significant genes in each direction.
up <- subset(volcano_df, direction == "Upregulated")
down <- subset(volcano_df, direction == "Downregulated")

up <- up[order(-up$log2FoldChange, up$padj), , drop = FALSE]
down <- down[order(down$log2FoldChange, down$padj), , drop = FALSE]

top_up <- head(up, 5)
top_down <- head(down, 5)
label_df <- rbind(top_up, top_down)

# Assign labels only to the selected rows.
volcano_df$label <- ""
selected <- rownames(volcano_df) %in% rownames(label_df)
volcano_df$label[selected] <- volcano_df$gene_name[selected]

# Create the plot.
volcano_plot <- ggplot2::ggplot(
  volcano_df,
  ggplot2::aes(x = log2FoldChange, y = neg_log10_padj)
) +
  ggplot2::geom_point(
    ggplot2::aes(color = direction),
    alpha = 0.65,
    size = 1.6
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Upregulated" = "#D73027",
      "Downregulated" = "#4575B4",
      "Below thresholds" = "grey75"
    )
  ) +
  ggplot2::geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey50"
  ) +
  ggplot2::geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    seed = 42,
    max.overlaps = Inf,
    box.padding = 0.5,
    min.segment.length = 0,
    size = 3.5
  ) +
  ggplot2::labs(
    title = paste(case, "vs", control),
    subtitle = "Top 5 significant genes per direction by log2 fold change",
    x = "Log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")

# Save PNG and PDF.
ggplot2::ggsave(
  file.path(out, "05_deseq2", "volcano_plot.png"),
  plot = volcano_plot,
  width = 10, height = 8, dpi = 300
)

ggplot2::ggsave(
  file.path(out, "05_deseq2", "volcano_plot.pdf"),
  plot = volcano_plot,
  width = 10, height = 8
)

# Save the labeled genes.
write.csv(
  label_df[, c(
    "gene_name", "log2FoldChange", "padj", "direction"
  ), drop = FALSE],
  file.path(out, "05_deseq2", "volcano_top_genes.csv"),
  row.names = FALSE
)
rankdf <- resdf[is.finite(resdf$stat), c("gene_id", "stat")]
rankdf$symbol <- unname(gene_map[sub("\\.[0-9]+$", "", rankdf$gene_id)])
rankdf <- rankdf[!is.na(rankdf$symbol) & nzchar(rankdf$symbol), ]
rankdf <- rankdf[order(abs(rankdf$stat), decreasing=TRUE), ]
rankdf <- rankdf[!duplicated(rankdf$symbol), ]
if (nrow(rankdf) < 15L) stop("Too few gene symbols matched the GTF; check reference compatibility")
ranks <- setNames(rankdf$stat, rankdf$symbol)
ranks <- sort(ranks, decreasing=TRUE)
write.csv(data.frame(symbol=names(ranks), wald_stat=unname(ranks)),
          file.path(out,"06_gsea","ranked_genes.csv"), row.names=FALSE)
pathways <- fgsea::gmtPathways(gmt)
set.seed(42)
enrichment <- fgsea::fgsea(pathways=pathways, stats=ranks, minSize=15, maxSize=500)
enrichment <- as.data.frame(enrichment)
# Save running enrichment curves for the top six pathways by adjusted p-value.
if (nrow(enrichment)>0L) {
  valid <- enrichment[is.finite(enrichment$padj), ,drop=FALSE]
  valid <- valid[order(valid$padj), ,drop=FALSE]
  if (nrow(valid)>0L) {
    pdf(file.path(out,"06_gsea","top_pathway_enrichment.pdf"),width=9,height=6)
    for (pname in head(valid$pathway,6)) {
      print(fgsea::plotEnrichment(pathways[[pname]], ranks) +
              ggplot2::labs(title=pname,subtitle=paste(case,"vs",control)))
    }
    dev.off()
    plotdf <- head(valid,20)
    plotdf$pathway <- factor(plotdf$pathway,levels=rev(plotdf$pathway))
    gsea_plot <- ggplot2::ggplot(plotdf,ggplot2::aes(x=NES,y=pathway,fill=NES>0)) +
      ggplot2::geom_col() + ggplot2::scale_fill_manual(values=c("#4575B4","#D73027"),guide="none") +
      ggplot2::labs(title="Top gene sets by adjusted p-value",x="NES",y=NULL) +
      ggplot2::theme_classic(base_size=10)
    ggplot2::ggsave(file.path(out,"06_gsea","top_pathways.pdf"),gsea_plot,width=11,height=8)
  }
}
if (nrow(enrichment)) {
  enrichment$leadingEdge <- vapply(enrichment$leadingEdge, paste, collapse=";", FUN.VALUE=character(1))
  enrichment <- enrichment[order(enrichment$padj), ]
}
write.csv(enrichment, file.path(out,"06_gsea","enrichment.csv"), row.names=FALSE)
cat("Done:", nrow(dds), "tested genes;", length(ranks), "genes ranked by symbol;", nrow(enrichment), "pathways tested.\n")

if (nrow(enrichment)>0L) {
  write.csv(subset(enrichment,!is.na(padj) & padj<0.05),
            file.path(out,"06_gsea","significant_pathways.csv"),row.names=FALSE)
}
