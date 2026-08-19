#!/usr/bin/env Rscript --vanilla

# Script for importing and processing transcript-level quantifications.
# Written by Karla Ruiz inspired by Lorena's Pantano code.
# Local copy: patched to support oarfish (*.quant) quantifications for long-read DTU.

# Loading required libraries
library(SummarizedExperiment)
library(tximport)

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(y) scan(text=y, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(y) y[2]), names = lapply(args_vals, function(y) y[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Build a table from a SummarizedExperiment object

build_table <- function(se.obj, slot) {
    data.frame(cbind(rowData(se.obj)[,1:2], assays(se.obj)[[slot]]), check.names = FALSE)
}

#' Write a table to a file from a SummarizedExperiment object with given parameters

write_se_table <- function(params, prefix) {
    file_name <- paste0(prefix, ".", params\$suffix)
    write.table(build_table(params\$obj, params\$slot), file_name,
                sep="\t", quote=FALSE, row.names = FALSE)
}

#' Read Transcript Metadata from a Given Path

read_transcript_info <- function(tinfo_path, tx_col, gene_id_col, gene_name_col){
    info <- file.info(tinfo_path)
    if (info\$size == 0) {
        stop("tx2gene file is empty")
    }

    # Read file with actual header to handle variable column counts correctly
    raw_info <- read.csv(tinfo_path, sep="\t", header = TRUE, check.names = FALSE)

    # Select columns by name, falling back to position if name not found
    transcript_info <- data.frame(
        tx = raw_info[[if (tx_col %in% colnames(raw_info)) tx_col else 1]],
        gene_id = raw_info[[if (gene_id_col %in% colnames(raw_info)) gene_id_col else 2]],
        gene_name = raw_info[[if (gene_name_col %in% colnames(raw_info)) gene_name_col else 3]],
        check.names = FALSE
    )

    extra <- setdiff(rownames(txi[[1]]), as.character(transcript_info[["tx"]]))
    transcript_info <- rbind(transcript_info, data.frame(tx=extra, gene_id=extra, gene_name=extra, check.names = FALSE))
    transcript_info <- transcript_info[match(rownames(txi[[1]]), transcript_info[["tx"]]), ]
    rownames(transcript_info) <- transcript_info[["tx"]]

    list(transcript = transcript_info,
        gene = unique(transcript_info[,2:3]),
        tx2gene = transcript_info[,1:2])
}

#' Create a SummarizedExperiment Object

create_summarized_experiment <- function(counts, abundance, length, col_data, row_data) {
    SummarizedExperiment(assays = list(counts = counts, abundance = abundance, length = length),
        colData = col_data,
        rowData = row_data)
}

################################################
################################################
## Main script starts here                    ##
################################################
################################################

# Set defaults and classes
opt <- list(
    tx_col = "transcript_id",
    gene_id_col = "gene_id",
    gene_name_col = "gene_name"
)

# Apply parameter overrides from ext.args
args_opt <- parse_args('$task.ext.args')
for (ao in names(args_opt)) {
    if (ao %in% names(opt)) {
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Define file pattern per quantification type.
# oarfish writes flat <sample>.quant files; salmon/kallisto live in per-sample dirs.
pattern <- switch('$quant_type',
    kallisto = "abundance.tsv",
    oarfish  = "\\\\.quant\$",
    "quant.sf")
fns <- list.files('quants', pattern = pattern, recursive = T, full.names = T)

if ('$quant_type' == 'oarfish') {
    names <- sub('\\\\.quant\$', '', basename(fns))
} else {
    names <- basename(dirname(fns))
}
names(fns) <- names

dropInfReps <- '$quant_type' == "oarfish"

# Import transcript-level quantifications
txi <- tximport(fns, type = '$quant_type', txOut = TRUE, dropInfReps = dropInfReps)

# Read transcript and sample data
transcript_info <- read_transcript_info('$tx2gene', opt\$tx_col, opt\$gene_id_col, opt\$gene_name_col)

# Make coldata just to appease the summarizedexperiment
coldata <- data.frame(files = fns, names = names, check.names = FALSE)
rownames(coldata) <- coldata[["names"]]

# Create initial SummarizedExperiment object
se <- create_summarized_experiment(txi[["counts"]], txi[["abundance"]], txi[["length"]],
    DataFrame(coldata), transcript_info\$transcript)

# Setting parameters for writing tables
params <- list(
    list(obj = se, slot = "abundance", suffix = "transcript_tpm.tsv"),
    list(obj = se, slot = "counts", suffix = "transcript_counts.tsv"),
    list(obj = se, slot = "length", suffix = "transcript_lengths.tsv")
)

# Process gene-level data if tx2gene mapping is available
if ("tx2gene" %in% names(transcript_info) && !is.null(transcript_info\$tx2gene)) {
    tx2gene <- transcript_info\$tx2gene
    gi <- summarizeToGene(txi, tx2gene = tx2gene)
    gi.ls <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
    gi.s <- summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "scaledTPM")

    gene_info <- transcript_info\$gene[match(rownames(gi[[1]]), transcript_info\$gene[["gene_id"]]),]
    rownames(gene_info) <- NULL
    col_data_frame <- DataFrame(coldata)

    # Create gene-level SummarizedExperiment objects
    gse <- create_summarized_experiment(gi[["counts"]], gi[["abundance"]], gi[["length"]],
        col_data_frame, gene_info)
    gse.ls <- create_summarized_experiment(gi.ls[["counts"]], gi.ls[["abundance"]], gi.ls[["length"]],
        col_data_frame, gene_info)
    gse.s <- create_summarized_experiment(gi.s[["counts"]], gi.s[["abundance"]], gi.s[["length"]],
        col_data_frame, gene_info)

    params <- c(params, list(
        list(obj = gse, slot = "length", suffix = "gene_lengths.tsv"),
        list(obj = gse, slot = "abundance", suffix = "gene_tpm.tsv"),
        list(obj = gse, slot = "counts", suffix = "gene_counts.tsv"),
        list(obj = gse.ls, slot = "counts", suffix = "gene_counts_length_scaled.tsv"),
        list(obj = gse.s, slot = "counts", suffix = "gene_counts_scaled.tsv")
    ))
}

# Writing tables for each set of parameters

prefix <- ''
if ('$task.ext.prefix' != 'null'){
    prefix = '$task.ext.prefix'
} else if ('$meta.id' != 'null'){
    prefix = '$meta.id'
}

done <- lapply(params, write_se_table, prefix)

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste(prefix, "R_sessionInfo.log", sep = '.'))
citation("tximeta")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
tximeta.version <- as.character(packageVersion('tximeta'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-tximeta:', tximeta.version)
    ),
'versions.yml')
