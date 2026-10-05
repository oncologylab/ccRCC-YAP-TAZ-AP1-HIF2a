args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop("Usage: Rscript normalize_csaw.R samples.tsv blacklist.bed bin_width output.tsv")
if (!requireNamespace("csaw", quietly=TRUE) || !requireNamespace("rtracklayer", quietly=TRUE)) stop("Install Bioconductor csaw and rtracklayer")
anno <- read.delim(args[1], stringsAsFactors=FALSE)
if (!all(c("sample", "bam", "normalization_group") %in% names(anno))) stop("Required columns: sample, bam, normalization_group")
if (!nrow(anno) || anyNA(anno[,c("sample","bam","normalization_group")]) || anyDuplicated(anno$sample)) stop("Invalid sample identifiers/groups")
if (!all(file.exists(anno$bam)) || !file.exists(args[2])) stop("Missing BAM or blacklist")
width <- as.integer(args[3])
if (is.na(width) || width < 1) stop("Invalid bin width")
if (file.exists(args[4])) stop("Refusing to overwrite output")
blacklist <- rtracklayer::import(args[2])
param <- csaw::readParam(minq=20, max.frag=800, discard=blacklist, pe="both")
anno$composition_factor <- NA_real_
for (group in unique(anno$normalization_group)) {
    selected <- which(anno$normalization_group == group)
    if (length(selected) < 2) stop("Each normalization group requires at least two libraries")
    bins <- csaw::windowCounts(anno$bam[selected], bin=TRUE, width=width, param=param)
    anno$composition_factor[selected] <- csaw::normFactors(bins, se.out=FALSE)
}
anno$bin_width <- width
write.table(anno, args[4], sep="\t", quote=FALSE, row.names=FALSE)
