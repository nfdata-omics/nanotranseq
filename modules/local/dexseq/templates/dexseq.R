#!/usr/bin/env Rscript

library(DEXSeq)

# Nextflow template variables
design_file <- "${design_file}"
gtf_file <- "${gtf}"

# Read design file
samples <- read.csv(design_file, stringsAsFactors=FALSE, header=TRUE)
colnames(samples) <- c("sample_id", "condition")

# Helper to read and merge count files
files <- list.files(pattern = "transcript_counts.tsv\$")
merged_file_pattern <- "all_samples.*transcript_counts.tsv\$"
merged_files <- list.files(pattern = merged_file_pattern)

countData <- NULL

if (length(merged_files) > 0) {
    message("Detected merged count file: ", merged_files[1])
    d_merged <- read.table(merged_files[1], header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

    # Rename columns to match design
    colnames(d_merged)[1] <- "feature_id"
    colnames(d_merged)[2] <- "gene_id"

    # Ensure all samples are present
    missing_cols <- setdiff(samples\$sample_id, colnames(d_merged))
    if (length(missing_cols) > 0) {
        stop("Missing columns for samples: ", paste(missing_cols, collapse=", "))
    }

    # Format for DEXSeq
    # DEXSeq usually works on exons, but we are using it for transcripts here
    # The countData should be a matrix
    countData <- as.matrix(d_merged[, samples\$sample_id])
    rownames(countData) <- d_merged\$feature_id

    geneIDs <- as.character(d_merged\$gene_id)
    featureIDs <- as.character(d_merged\$feature_id)

} else {
    # Individual files logic (simplified)
    # Assumes files match sample IDs
    message("Reading individual files...")
    # ... implementation similar to DRIMSeq ...
    # For now assuming merged file as per current pipeline state
    stop("Individual file support not fully implemented in DEXSeq template yet. Please provide merged file.")
}

# DEXSeq Analysis
tryCatch({

    dxd <- DEXSeqDataSet(
        countData = countData,
        sampleData = samples,
        design = ~ sample + exon + condition:exon,
        featureID = featureIDs,
        groupID = geneIDs
    )

    # Normalization
    dxd <- estimateSizeFactors(dxd)

    # Dispersion estimation
    dxd <- estimateDispersions(dxd)

    # Testing
    dxd <- testForDEU(dxd)

    # Fold changes
    dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")

    # Results
    dxr <- DEXSeqResults(dxd)

    # Convert to standard format for IsoformSwitchAnalyzeR
    # It expects: isoform_id, pvalue, adj_pvalue (padj)

    res_df <- as.data.frame(dxr)

    # We need to map back to isoform/gene structure
    # DEXSeq results columns: groupID (gene), featureID (isoform)

    # Output gene results (aggregated)
    # DEXSeq gives per-exon (per-isoform) results. Gene-level aggregation is usually q-value of the gene.
    # IsoformSwitchAnalyzeR uses the smallest q-value of any isoform as the gene q-value usually.

    res_isoform <- res_df[, c("groupID", "featureID", "pvalue", "padj")]
    colnames(res_isoform) <- c("gene_id", "feature_id", "pvalue", "adj_pvalue")

    # Create gene level results (min padj per gene)
    res_gene <- aggregate(adj_pvalue ~ gene_id, data=res_isoform, FUN=min)
    colnames(res_gene) <- c("gene_id", "adj_pvalue")

    write.csv(res_gene, "dexseq_gene_results.csv", row.names=FALSE)
    write.csv(res_isoform, "dexseq_transcript_results.csv", row.names=FALSE)
    saveRDS(dxd, "dexseq_object.rds")

    # Plotting
    pdf("dexseq_plots.pdf")
    plotDispEsts(dxd)
    plotMA(dxr, cex=0.8)
    dev.off()

    # Versions
    r_version <- paste(R.version\$major, R.version\$minor, sep=".")
    dexseq_version <- as.character(packageVersion('DEXSeq'))
    versions_text <- paste0(
        '"', "${task.process}", '":\n',
        '    r-base: ', r_version, '\n',
        '    bioconductor-dexseq: ', dexseq_version, '\n'
    )
    writeLines(versions_text, "versions.yml")

}, error = function(e) {
    message("DEXSeq failed: ", e\$message)
    write.csv(data.frame(), "dexseq_gene_results.csv")
    write.csv(data.frame(), "dexseq_transcript_results.csv")

    # Versions
    r_version <- paste(R.version\$major, R.version\$minor, sep=".")
    dexseq_version <- as.character(packageVersion('DEXSeq'))
    versions_text <- paste0(
        '"', "${task.process}", '":\n',
        '    r-base: ', r_version, '\n',
        '    bioconductor-dexseq: ', dexseq_version, '\n'
    )
    writeLines(versions_text, "versions.yml")
})
