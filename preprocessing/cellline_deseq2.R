args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop("Usage: Rscript cellline_deseq2.R counts.tsv samples.tsv annotations.tsv output_directory")
if (!requireNamespace("DESeq2", quietly=TRUE)) stop("Install Bioconductor DESeq2")
counts <- as.matrix(read.delim(args[1], row.names=1, check.names=FALSE))
samples <- read.delim(args[2], stringsAsFactors=FALSE)
annotations <- read.delim(args[3], stringsAsFactors=FALSE)
if (!all(c("sample", "condition") %in% names(samples))) stop("Required columns: sample, condition")
if (!all(c("ENS_ID", "V2") %in% names(annotations))) stop("Required columns: ENS_ID, V2")
conditions <- c("PAR", "shYT", "shHIF2A", "shJUNs")
if (!setequal(samples$condition, conditions) || any(table(samples$condition) != 3)) stop("Three replicates per condition required")
if (anyDuplicated(samples$sample) || !setequal(samples$sample, colnames(counts))) stop("Sample identifiers differ")
if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) stop("Nonnegative integer counts required")
rownames(counts) <- sub("\\..*$", "", rownames(counts))
if (anyDuplicated(rownames(counts)) || anyDuplicated(annotations$ENS_ID)) stop("Duplicate gene identifiers")
samples <- samples[order(match(samples$condition, conditions)), , drop=FALSE]
counts <- counts[, samples$sample, drop=FALSE]
keep <- rowSums(counts[, samples$condition == "PAR", drop=FALSE]) > 15 |
        rowSums(counts[, samples$condition == "shHIF2A", drop=FALSE]) > 15 |
        rowSums(counts[, samples$condition == "shJUNs", drop=FALSE]) > 15
counts <- counts[keep, , drop=FALSE] + 2
rownames(samples) <- samples$sample
samples$condition <- factor(samples$condition, levels=conditions)
if (file.exists(args[4])) stop("Choose a new output directory")
dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts, colData=samples, design=~condition)
dds <- DESeq2::DESeq(dds, test="Wald", fitType="parametric", betaPrior=FALSE,
                    minReplicatesForReplace=Inf, modelMatrixType="standard", parallel=FALSE)
dir.create(args[4], recursive=TRUE)
normalized <- DESeq2::counts(dds, normalized=TRUE)
write.table(normalized, file.path(args[4], "normalized_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)
for (condition in conditions[-1]) {
    result <- DESeq2::results(dds, contrast=c("condition", condition, "PAR"), lfcThreshold=0,
                             cooksCutoff=FALSE, independentFiltering=FALSE, pAdjustMethod="BH", test="Wald")
    result <- as.data.frame(result)
    result$ENS_ID <- rownames(result)
    result$V2 <- annotations$V2[match(result$ENS_ID, annotations$ENS_ID)]
    name <- if (condition == "shHIF2A") "shHIF2a" else condition
    write.table(result, file.path(args[4], paste0("DESeq2_", name, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
}
saveRDS(dds, file.path(args[4], "DESeq2_fit.rds"))
capture.output(sessionInfo(), file=file.path(args[4], "sessionInfo.txt"))
