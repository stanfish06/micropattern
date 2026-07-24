#!/usr/bin/env Rscript

# Prepare provisional, batch-corrected RNA clusters for ArchR pseudobulk
# peak calling. These groups are technical peak-calling strata, not final
# biological annotations.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(harmony)
  library(Matrix)
  library(Seurat)
  library(SingleCellExperiment)
  library(zellkonverter)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "."
root <- normalizePath(root, mustWork = TRUE)

input_h5ad <- file.path(
  root,
  "data/GSE245998/processed/gse245998_rna_qc.h5ad"
)
output_csv <- file.path(
  root,
  "data/GSE245998/processed/gse245998_peak_calling_groups.csv"
)
summary_csv <- file.path(
  root,
  "data/GSE245998/processed/gse245998_peak_calling_group_summary.csv"
)
marker_csv <- file.path(
  root,
  "data/GSE245998/processed/gse245998_peak_calling_marker_means.csv"
)

set.seed(245998)

message("Reading retained RNA counts: ", input_h5ad)
sce <- readH5AD(input_h5ad, use_hdf5 = FALSE)
if (!"X" %in% assayNames(sce)) {
  stop("Expected count assay 'X' in the RNA H5AD.", call. = FALSE)
}

counts <- assay(sce, "X")
if (!inherits(counts, "sparseMatrix")) {
  counts <- as(counts, "dgCMatrix")
}
if (any(counts@x < 0) || any(counts@x != round(counts@x))) {
  stop("RNA assay X is not a non-negative integer count matrix.", call. = FALSE)
}

metadata <- as.data.frame(colData(sce))
rownames(metadata) <- colnames(sce)
if (!all(c("batch", "sample", "genotype") %in% colnames(metadata))) {
  stop("RNA metadata lacks batch, sample, or genotype.", call. = FALSE)
}
metadata$batch <- factor(metadata$batch)

message("Creating Seurat object for ", ncol(counts), " retained nuclei.")
obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  project = "GSE245998"
)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)
obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  verbose = FALSE
)
obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 50,
  verbose = FALSE
)
obj <- RunHarmony(
  obj,
  group.by.vars = "batch",
  reduction.use = "pca",
  dims.use = 1:40,
  verbose = FALSE
)
obj <- FindNeighbors(
  obj,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)
obj <- FindClusters(
  obj,
  resolution = 0.8,
  algorithm = 1,
  random.seed = 245998,
  verbose = FALSE
)

clusters <- as.character(Idents(obj))
cluster_sizes <- sort(table(clusters))
if (any(cluster_sizes < 40)) {
  stop(
    "RNA clustering produced peak-calling groups smaller than 40 cells:\n",
    paste(
      names(cluster_sizes[cluster_sizes < 40]),
      cluster_sizes[cluster_sizes < 40],
      sep = "=",
      collapse = "\n"
    ),
    "\nReduce the clustering resolution before ArchR peak calling.",
    call. = FALSE
  )
}

group_names <- paste0("rna_cluster_", clusters)
groups <- data.frame(
  cell_name = colnames(obj),
  peak_calling_group = group_names,
  stringsAsFactors = FALSE
)
write.csv(groups, output_csv, row.names = FALSE)

group_summary <- aggregate(
  rep(1L, nrow(obj[[]])),
  by = list(
    peak_calling_group = group_names,
    sample = as.character(obj$sample),
    genotype = as.character(obj$genotype),
    batch = as.character(obj$batch)
  ),
  FUN = sum
)
colnames(group_summary)[ncol(group_summary)] <- "n_cells"
group_summary <- group_summary[
  order(group_summary$peak_calling_group, group_summary$sample),
]
write.csv(group_summary, summary_csv, row.names = FALSE)

pgc_markers <- c(
  "NANOS3",
  "TFAP2C",
  "PRDM1",
  "SOX17",
  "KIT",
  "DPPA3",
  "PRDM14",
  "NANOG",
  "POU5F1"
)
markers_present <- intersect(pgc_markers, rownames(obj))
if (length(markers_present) == 0) {
  warning("None of the requested PGC markers were present.")
} else {
  normalized <- LayerData(obj, assay = "RNA", layer = "data")
  marker_matrix <- normalized[markers_present, , drop = FALSE]
  marker_means <- vapply(
    split(seq_len(ncol(obj)), group_names),
    function(indices) Matrix::rowMeans(marker_matrix[, indices, drop = FALSE]),
    numeric(length(markers_present))
  )
  marker_means <- as.data.frame(t(marker_means), check.names = FALSE)
  marker_means$peak_calling_group <- rownames(marker_means)
  marker_means$n_cells <- as.integer(table(group_names)[
    marker_means$peak_calling_group
  ])
  marker_means <- marker_means[
    ,
    c("peak_calling_group", "n_cells", markers_present),
    drop = FALSE
  ]
  write.csv(marker_means, marker_csv, row.names = FALSE)
}

message("Peak-calling groups:")
print(sort(table(groups$peak_calling_group)))
message("Wrote: ", output_csv)
message("Wrote: ", summary_csv)
message("Wrote: ", marker_csv)
