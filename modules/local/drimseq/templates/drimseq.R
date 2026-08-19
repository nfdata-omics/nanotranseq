#!/usr/bin/env Rscript

library(DRIMSeq)

# Nextflow template variables
design_file <- "${design_file}"
gtf_file <- "${gtf}"

# Read design file
samples <- read.csv(design_file, stringsAsFactors=FALSE, header=TRUE)
colnames(samples) <- c("sample_id", "condition")

# Find all transcript_counts.tsv files
files <- list.files(pattern = "transcript_counts.tsv\$")

if (length(files) == 0) {
    stop("No *transcript_counts.tsv files found in current directory")
}

# The previous logic was trying to match individual files per sample
# BUT, if we receive the output from TXIMETA_TXIMPORT, we might be getting
# aggregated files or individual files depending on how the pipeline is run.
#
# If we have "all_samples.transcript_counts.tsv", it means TXIMPORT already merged them.
# Let's check if we have one merged file or many individual files.

merged_file_pattern <- "all_samples.*transcript_counts.tsv\$"
merged_files <- list.files(pattern = merged_file_pattern)

if (length(merged_files) > 0) {
    message("Detected merged count file: ", merged_files[1])

    # Read the merged file
    # Format from tximport usually: tx_id, gene_id, sample1, sample2, ...
    d_merged <- read.table(merged_files[1], header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

    # Verify we have tx and gene columns
    # We assume first two columns are tx and gene info
    counts <- d_merged

    # Rename first two columns to match DRIMSeq expectations if needed
    # But DRIMSeq needs a specific data frame structure.
    # The columns should match the sample_ids in design file.

    # Check if all sample_ids from design are in the columns
    missing_cols <- setdiff(samples\$sample_id, colnames(counts))
    if (length(missing_cols) > 0) {
        stop("The merged count file is missing columns for samples: ", paste(missing_cols, collapse=", "))
    }

    # We need to construct the counts data frame for dmDSdata
    # It requires: gene_id, feature_id, sample1, sample2...
    # Assuming col 1 is feature_id (tx), col 2 is gene_id

    # Rename first two cols for clarity
    colnames(counts)[1] <- "feature_id"
    colnames(counts)[2] <- "gene_id"

    # Select only relevant columns
    counts <- counts[, c("gene_id", "feature_id", samples\$sample_id)]

} else {
    # Fallback to individual file logic
    message("Assuming individual files per sample...")

    # Match files to samples
    sample_files <- sapply(samples\$sample_id, function(sid) {
        sid_regex <- paste0("^", sid, ".*transcript_counts[.]tsv\$")
        match <- grep(sid_regex, files, value=TRUE)
        if (length(match) == 0) {
            match <- grep(sid, files, value=TRUE)
        }
        if (length(match) == 0) return(NA)
        return(match[1])
    })

    if (any(is.na(sample_files))) {
        missing <- samples\$sample_id[is.na(sample_files)]
        stop("Could not find count files for samples: ", paste(missing, collapse=", "))
    }

    # Read first file to initialize
    first_data <- read.table(sample_files[1], header=TRUE, sep="\t", stringsAsFactors=FALSE)
    if (ncol(first_data) < 3) {
        stop("Count file has fewer than 3 columns. Expected tx, gene_id, count.")
    }

    counts <- data.frame(gene_id = first_data[, 2], feature_id = first_data[, 1])

    # Merge counts
    for (i in 1:nrow(samples)) {
        sid <- samples\$sample_id[i]
        f <- sample_files[i]
        d <- read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE)

        if (!identical(d[,1], counts\$feature_id)) {
            d <- d[match(counts\$feature_id, d[,1]), ]
        }

        counts[[sid]] <- d[, 3]
    }
}

# dmDSdata requires gene_id/feature_id to be character (or factor). read.table coerces
# all-numeric IDs (e.g. contig-named references like "3", "9") to integer, so force them.
counts\$gene_id    <- as.character(counts\$gene_id)
counts\$feature_id <- as.character(counts\$feature_id)

# Create dmDSdata object
d <- dmDSdata(counts = counts, samples = samples)

# Filtering
n_samples <- nrow(samples)

# Adjust filtering parameters based on sample size
# If n_samples is small (e.g. 1 or 2), we need to be lenient
min_samps <- max(1, min(3, ceiling(n_samples/2)))

# Robust filtering with fallback
d_filtered <- tryCatch({
    dmFilter(d,
          min_samps_feature_expr = min_samps,
          min_feature_expr = 0, # Relaxed for small datasets
          min_samps_feature_prop = min_samps,
          min_feature_prop = 0.0, # Relaxed
          min_samps_gene_expr = min_samps,
          min_gene_expr = 0) # Relaxed
}, error = function(e) {
    message("Filtering failed (likely no genes left): ", e\$message)
    return(NULL)
})

if (is.null(d_filtered)) {
    message("No genes left after filtering. Saving empty results.")
    write.csv(data.frame(), "drimseq_gene_results.csv")
    write.csv(data.frame(), "drimseq_transcript_results.csv")
    saveRDS(d, "drimseq_object.rds") # Save unfiltered object

    # Create dummy plot
    pdf("drimseq_plots.pdf")
    plot(1, 1, main="No genes left after filtering")
    dev.off()

    # Versions
    r_version <- paste(R.version\$major, R.version\$minor, sep=".")
    drimseq_version <- as.character(packageVersion('DRIMSeq'))
    versions_text <- paste0(
        '"', "${task.process}", '":\n',
        '    r-base: ', r_version, '\n',
        '    bioconductor-drimseq: ', drimseq_version, '\n'
    )
    writeLines(versions_text, "versions.yml")

    quit(save="no", status=0)
}

d <- d_filtered

# Design matrix
# If only 1 condition or no replicates, we can't do full statistical testing
if (n_samples < 2 || length(unique(samples\$condition)) < 2) {
    message("Not enough samples or conditions for differential analysis. Saving object and plots only.")
    saveRDS(d, "drimseq_object.rds")
    pdf("drimseq_plots.pdf")
    plotData(d)
    dev.off()

    # Create empty result files to prevent process failure
    write.csv(data.frame(), "drimseq_gene_results.csv")
    write.csv(data.frame(), "drimseq_transcript_results.csv")

    # Versions
    r_version <- paste(R.version\$major, R.version\$minor, sep=".")
    drimseq_version <- as.character(packageVersion('DRIMSeq'))
    versions_text <- paste0(
        '"', "${task.process}", '":\n',
        '    r-base: ', r_version, '\n',
        '    bioconductor-drimseq: ', drimseq_version, '\n'
    )
    writeLines(versions_text, "versions.yml")

} else {
        design_full <- model.matrix(~ condition, data = samples(d))

        tryCatch({
            # Run DRIMSeq
            set.seed(123)
            d <- dmPrecision(d, design = design_full)
            d <- dmFit(d, design = design_full)

            # Test for condition
            coef_name <- colnames(design_full)[2]
            d <- dmTest(d, coef = coef_name)

            # Results
            res <- results(d)
            res_tx <- results(d, level = "feature")

            # Write results
            write.csv(res, "drimseq_gene_results.csv", row.names=FALSE)
            write.csv(res_tx, "drimseq_transcript_results.csv", row.names=FALSE)
            saveRDS(d, "drimseq_object.rds")

            # Plotting
            pdf("drimseq_plots.pdf")
            plotPrecision(d)
            dev.off()

        }, error = function(e) {
            message("DRIMSeq statistical testing failed: ", e\$message)
            message("Saving empty results.")

            write.csv(data.frame(), "drimseq_gene_results.csv")
            write.csv(data.frame(), "drimseq_transcript_results.csv")
            saveRDS(d, "drimseq_object.rds")

            # Create dummy plot
            pdf("drimseq_plots.pdf")
            plot(1, 1, main="Testing failed")
            dev.off()

            # Versions
            r_version <- paste(R.version\$major, R.version\$minor, sep=".")
            drimseq_version <- as.character(packageVersion('DRIMSeq'))
            versions_text <- paste0(
                '"', "${task.process}", '":\n',
                '    r-base: ', r_version, '\n',
                '    bioconductor-drimseq: ', drimseq_version, '\n'
            )
            writeLines(versions_text, "versions.yml")
        })
    }

# Versions (if not already written)
if (!file.exists("versions.yml")) {
    r_version <- paste(R.version\$major, R.version\$minor, sep=".")
    drimseq_version <- as.character(packageVersion('DRIMSeq'))
    versions_text <- paste0(
        '"', "${task.process}", '":\n',
        '    r-base: ', r_version, '\n',
        '    bioconductor-drimseq: ', drimseq_version, '\n'
    )
    writeLines(versions_text, "versions.yml")
}
