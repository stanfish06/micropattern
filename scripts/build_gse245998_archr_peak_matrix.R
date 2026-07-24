#!/usr/bin/env Rscript

# Build a shared GSE245998 cell-by-consensus-peak matrix with ArchR.
#
# Usage:
#   Rscript scripts/build_gse245998_archr_peak_matrix.R \
#     /home/stan/Git/micropattern \
#     data/GSE245998/processed/gse245998_peak_calling_groups.csv \
#     data/GSE245998/processed/archr_peak_matrix
#
# The group CSV must contain:
#   cell_name,peak_calling_group
# where cell_name uses the same "sample#10x-barcode" convention as the QC
# AnnData, for example eb01#AAACAGCCAATAAGCA-1.
#
# Direct requirements (tested with R 4.5.1):
#   ArchR, Matrix, GenomicRanges, S4Vectors, rhdf5,
#   SummarizedExperiment, SingleCellExperiment, zellkonverter
# External executable:
#   MACS2 (on PATH or supplied through MACS2_PATH)

options(stringsAsFactors = FALSE)

required_packages <- c(
  "ArchR",
  "GenomicRanges",
  "Matrix",
  "rhdf5",
  "S4Vectors",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "zellkonverter"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# ArchR's bundled genome annotation lookup uses data() without an explicit
# package name, so the package must be attached rather than only namespace-loaded.
suppressPackageStartupMessages(library(ArchR))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "."
root <- normalizePath(root, mustWork = TRUE)

resolve_path <- function(path) {
  if (grepl("^/", path)) {
    return(path)
  }
  file.path(root, path)
}

group_csv <- if (length(args) >= 2) {
  resolve_path(args[[2]])
} else {
  file.path(
    root,
    "data/GSE245998/processed/gse245998_peak_calling_groups.csv"
  )
}
output_dir <- if (length(args) >= 3) {
  resolve_path(args[[3]])
} else {
  file.path(root, "data/GSE245998/processed/archr_peak_matrix")
}

threads <- as.integer(Sys.getenv("ARCHR_THREADS", "12"))
force <- tolower(Sys.getenv("ARCHR_FORCE", "false")) %in% c(
  "1", "true", "yes"
)
if (is.na(threads) || threads < 1) {
  stop("ARCHR_THREADS must be a positive integer.", call. = FALSE)
}

samples <- c("eb01", "eb02", "eb03", "eb06", "eb07", "eb08")
fragment_accessions <- c(
  eb01 = "GSM7853056",
  eb02 = "GSM7853057",
  eb03 = "GSM7853058",
  eb06 = "GSM7853059",
  eb07 = "GSM7853060",
  eb08 = "GSM7853061"
)
fragment_files <- file.path(
  root,
  "data/GSE245998/raw",
  paste0(
    unname(fragment_accessions[samples]),
    "_",
    samples,
    "_atac_fragments.tsv.gz"
  )
)
names(fragment_files) <- samples
barcode_files <- file.path(
  root,
  "data/GSE245998/processed/barcodes",
  paste0(samples, "_retained_barcodes.txt")
)
names(barcode_files) <- samples

required_files <- c(fragment_files, paste0(fragment_files, ".tbi"), barcode_files)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}
if (!file.exists(group_csv)) {
  stop(
    "Peak-calling group CSV is missing: ",
    group_csv,
    "\nCreate it from the retained RNA AnnData before running ArchR.",
    call. = FALSE
  )
}

valid_barcodes <- lapply(barcode_files, readLines, warn = FALSE)
valid_barcodes <- lapply(valid_barcodes, unique)
expected_cell_names <- unlist(
  Map(
    function(sample, barcodes) paste0(sample, "#", barcodes),
    names(valid_barcodes),
    valid_barcodes
  ),
  use.names = FALSE
)
if (anyDuplicated(expected_cell_names)) {
  stop("Retained sample-prefixed cell names are not unique.", call. = FALSE)
}

groups <- read.csv(group_csv, check.names = FALSE)
required_group_columns <- c("cell_name", "peak_calling_group")
missing_group_columns <- setdiff(required_group_columns, colnames(groups))
if (length(missing_group_columns) > 0) {
  stop(
    "Group CSV is missing columns: ",
    paste(missing_group_columns, collapse = ", "),
    call. = FALSE
  )
}
groups <- groups[, required_group_columns]
groups$cell_name <- as.character(groups$cell_name)
groups$peak_calling_group <- as.character(groups$peak_calling_group)
if (anyDuplicated(groups$cell_name)) {
  stop("Group CSV contains duplicate cell_name values.", call. = FALSE)
}
missing_groups <- setdiff(expected_cell_names, groups$cell_name)
extra_groups <- setdiff(groups$cell_name, expected_cell_names)
if (length(missing_groups) > 0 || length(extra_groups) > 0) {
  stop(
    "Group CSV cell mismatch. Missing retained cells: ",
    length(missing_groups),
    "; unexpected cells: ",
    length(extra_groups),
    call. = FALSE
  )
}
if (anyNA(groups$peak_calling_group) || any(groups$peak_calling_group == "")) {
  stop("Every retained cell must have a peak_calling_group.", call. = FALSE)
}
group_sizes <- sort(table(groups$peak_calling_group))
if (any(group_sizes < 40)) {
  stop(
    "Every peak-calling group must contain at least 40 cells. Small groups:\n",
    paste(
      names(group_sizes[group_sizes < 40]),
      group_sizes[group_sizes < 40],
      sep = "=",
      collapse = "\n"
    ),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
arrow_dir <- file.path(output_dir, "ArrowFiles")
project_dir <- file.path(output_dir, "ArchRProject")
run_work_dir <- file.path(output_dir, "ArchRRun")
dir.create(arrow_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run_work_dir, recursive = TRUE, showWarnings = FALSE)

initial_working_directory <- getwd()
on.exit(setwd(initial_working_directory), add = TRUE)
setwd(run_work_dir)

message("Project root: ", root)
message("Output directory: ", output_dir)
message("Retained cells: ", length(expected_cell_names))
message("Peak-calling groups:")
print(group_sizes)

ArchR::addArchRGenome("hg38")
ArchR::addArchRThreads(threads = threads)

project_arrow_files <- file.path(
  project_dir,
  "ArrowFiles",
  paste0(samples, ".arrow")
)
coverage_dir <- file.path(
  project_dir,
  "GroupCoverages",
  "peak_calling_group"
)
resume_from_coverages <- all(file.exists(project_arrow_files)) &&
  dir.exists(coverage_dir) &&
  length(list.files(coverage_dir, pattern = "\\.coverage\\.h5$")) > 0

if (resume_from_coverages) {
  message("Resuming from existing ArchR Arrow and group-coverage files.")
  arrow_files <- normalizePath(project_arrow_files, mustWork = TRUE)
} else {
  old_working_directory <- getwd()
  setwd(arrow_dir)
  arrow_files <- ArchR::createArrowFiles(
    inputFiles = fragment_files,
    sampleNames = samples,
    outputNames = samples,
    validBarcodes = valid_barcodes,
    minTSS = 0,
    minFrags = 0,
    maxFrags = 1e9,
    addTileMat = TRUE,
    addGeneScoreMat = FALSE,
    force = force,
    threads = threads
  )
  arrow_files <- normalizePath(arrow_files, mustWork = TRUE)
  setwd(old_working_directory)
}

proj <- ArchR::ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = project_dir,
  copyArrows = !resume_from_coverages
)

archr_cells <- ArchR::getCellNames(proj)
if (!setequal(archr_cells, expected_cell_names)) {
  stop(
    "ArchR retained-cell mismatch. Expected ",
    length(expected_cell_names),
    " cells but ArchR contains ",
    length(archr_cells),
    ".",
    call. = FALSE
  )
}
group_lookup <- setNames(groups$peak_calling_group, groups$cell_name)
proj <- ArchR::addCellColData(
  ArchRProj = proj,
  data = unname(group_lookup[archr_cells]),
  cells = archr_cells,
  name = "peak_calling_group",
  force = TRUE
)

if (resume_from_coverages) {
  message("Registering existing sample-aware group coverages.")
  coverage_files <- sort(list.files(
    coverage_dir,
    pattern = "\\.insertions\\.coverage\\.h5$",
    full.names = TRUE
  ))
  coverage_names <- sub(
    "\\.insertions\\.coverage\\.h5$",
    "",
    basename(coverage_files)
  )
  coverage_groups <- sub("\\._\\..*$", "", coverage_names)
  coverage_replicates <- sub("^.*\\._\\.", "", coverage_names)
  coverage_cells <- lapply(
    coverage_files,
    rhdf5::h5read,
    name = "Coverage/Info/CellNames"
  )
  if (!all(unlist(coverage_cells, use.names = FALSE) %in% archr_cells)) {
    stop(
      "Existing group coverage files contain cells outside the ArchR project.",
      call. = FALSE
    )
  }
  coverage_group_names <- unique(coverage_groups)
  coverage_cell_groups <- lapply(coverage_group_names, function(group_name) {
    indices <- which(coverage_groups == group_name)
    do.call(
      S4Vectors::SimpleList,
      setNames(coverage_cells[indices], coverage_replicates[indices])
    )
  })
  coverage_cell_groups <- do.call(
    S4Vectors::SimpleList,
    setNames(coverage_cell_groups, coverage_group_names)
  )
  coverage_params <- S4Vectors::SimpleList(
    groupBy = "peak_calling_group",
    minCells = 40L,
    maxCells = 500L,
    minReplicates = 2L,
    sampleRatio = 0.8,
    kmerLength = 6L,
    cellGroups = coverage_cell_groups
  )
  coverage_metadata <- S4Vectors::DataFrame(
    Group = coverage_groups,
    Name = coverage_names,
    File = normalizePath(coverage_files, mustWork = TRUE),
    nCells = vapply(coverage_cells, length, integer(1)),
    nInsertions = rep(NA_real_, length(coverage_files))
  )
  if (is.null(proj@projectMetadata$GroupCoverages)) {
    proj@projectMetadata$GroupCoverages <- S4Vectors::SimpleList()
  }
  proj@projectMetadata$GroupCoverages[["peak_calling_group"]] <-
    S4Vectors::SimpleList(
      Params = coverage_params,
      coverageMetadata = coverage_metadata
    )
  message("Registered ", length(coverage_files), " existing coverages.")
} else {
  proj <- ArchR::addGroupCoverages(
    ArchRProj = proj,
    groupBy = "peak_calling_group",
    force = force,
    threads = threads
  )
}
coverage_checkpoint <- file.path(
  output_dir,
  "archr_project_after_coverages.rds"
)
saveRDS(proj, coverage_checkpoint, compress = FALSE)
message("Saved coverage checkpoint: ", coverage_checkpoint)

# The Nix shell also contains Python 3.13 tools, while MACS2 is packaged with
# Python 3.11. Prevent the shell PYTHONPATH from leaking into MACS2.
Sys.unsetenv("PYTHONPATH")
macs2_path <- Sys.getenv("MACS2_PATH", "")
if (!nzchar(macs2_path)) {
  macs2_path <- ArchR::findMacs2()
}
if (!nzchar(macs2_path) || !file.exists(macs2_path)) {
  stop(
    "MACS2 was not found. Put macs2 on PATH or set MACS2_PATH.",
    call. = FALSE
  )
}
message("Using MACS2: ", macs2_path)

proj <- ArchR::addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "peak_calling_group",
  pathToMacs2 = macs2_path,
  force = force,
  threads = threads
)
proj <- ArchR::addPeakMatrix(
  ArchRProj = proj,
  binarize = FALSE,
  # PeakMatrix is a derived object tied to the newly called consensus peak set.
  # Overwrite any empty or stale group left by an interrupted prior attempt.
  force = TRUE,
  threads = threads
)
proj <- ArchR::saveArchRProject(
  ArchRProj = proj,
  outputDirectory = project_dir,
  load = TRUE
)

peak_se <- ArchR::getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  binarize = FALSE,
  threads = threads
)
if (ncol(peak_se) != length(expected_cell_names)) {
  stop(
    "Exported PeakMatrix contains ",
    ncol(peak_se),
    " cells; expected ",
    length(expected_cell_names),
    ".",
    call. = FALSE
  )
}

peak_ranges <- SummarizedExperiment::rowRanges(peak_se)
peak_names <- paste0(
  as.character(GenomicRanges::seqnames(peak_ranges)),
  ":",
  GenomicRanges::start(peak_ranges) - 1L,
  "-",
  GenomicRanges::end(peak_ranges)
)
if (anyDuplicated(peak_names)) {
  stop("Consensus peak names are not unique.", call. = FALSE)
}

assay_name <- SummarizedExperiment::assayNames(peak_se)[[1]]
peak_counts <- SummarizedExperiment::assay(peak_se, assay_name)
peak_sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = peak_counts),
  rowData = S4Vectors::DataFrame(
    peak = peak_names,
    chrom = as.character(GenomicRanges::seqnames(peak_ranges)),
    start = GenomicRanges::start(peak_ranges) - 1L,
    end = GenomicRanges::end(peak_ranges),
    row.names = peak_names
  ),
  colData = S4Vectors::DataFrame(
    cell_name = colnames(peak_se),
    sample = sub("#.*$", "", colnames(peak_se)),
    barcode = sub("^[^#]+#", "", colnames(peak_se)),
    row.names = colnames(peak_se)
  )
)

peak_rds <- file.path(output_dir, "gse245998_archr_peak_matrix.rds")
peak_h5ad <- file.path(output_dir, "gse245998_archr_peak_matrix.h5ad")
peak_bed <- file.path(output_dir, "gse245998_consensus_peaks.bed")
saveRDS(peak_sce, peak_rds, compress = FALSE)
write.table(
  data.frame(
    chrom = as.character(GenomicRanges::seqnames(peak_ranges)),
    start = GenomicRanges::start(peak_ranges) - 1L,
    end = GenomicRanges::end(peak_ranges),
    name = peak_names
  ),
  file = peak_bed,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
zellkonverter::writeH5AD(
  peak_sce,
  file = peak_h5ad,
  X_name = "counts",
  compression = "gzip",
  verbose = TRUE
)

summary <- list(
  n_cells = ncol(peak_sce),
  n_peaks = nrow(peak_sce),
  n_nonzero = Matrix::nnzero(peak_counts),
  peak_matrix_binarized = FALSE,
  coordinate_convention = "0-based half-open in H5AD var and BED",
  group_csv = normalizePath(group_csv),
  macs2_path = macs2_path,
  threads = threads
)
dput(
  summary,
  file = file.path(output_dir, "peak_matrix_summary.R")
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "sessionInfo.txt")
)

message("Peak matrix complete:")
message("  cells: ", summary$n_cells)
message("  peaks: ", summary$n_peaks)
message("  nonzero entries: ", summary$n_nonzero)
message("  H5AD: ", peak_h5ad)
message("  BED: ", peak_bed)
