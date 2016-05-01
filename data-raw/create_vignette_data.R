############################################################################
### Script to save simulated data set (100 genes) for regsplice vignette ###
### author: Lukas Weber                                                  ###
############################################################################

# data from Simulation5_Charlotte

# from paper:
# Soneson et al. (2016), "Isoform prefiltering improves performance of count-based
# methods for analysis of differential transcript usage", Genome Biology.
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0862-3

# original data files containing reads (FASTQ and BAM) available on ArrayExpress:
# http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3766/

# count files stored on taupo server:
# /home/Shared/data/alt_splicing_simulations/Simulation5_Charlotte/hsapiens/with_diffexpression/non_null_simulation/2_counts/dexseq_nomerge



#################
### load data ###
#################

DIR_COUNTS <- "../../data/Charlotte_simulation5/counts"
DIR_TRUTH <- "../../data/Charlotte_simulation5/truth"

files_counts <- list.files(DIR_COUNTS, full.names = TRUE)
file_truth <- list.files(DIR_TRUTH, full.names = TRUE)


# load counts

load_counts <- function(f) read.table(f, header = FALSE, sep = "\t", col.names = c("exon", "count"))

counts <- lapply(files_counts, load_counts)  # takes a few seconds
names(counts) <- paste0("sample", 1:6)

# rows with ambiguous, empty, low quality, or not aligned reads begin with underscore
remove_rows <- function(d) d[!grepl("^_.*$", d[, 1]), ]
counts <- lapply(counts, remove_rows)

# check exons are the same for each sample
dim(counts[[1]])
counts_dims <- lapply(counts, function(u) dim(u))
all(sapply(counts_dims, identical, dim(counts[[1]])))

head(counts[[1]]$exon)
exon_labels <- lapply(counts, function(u) u$exon)
all(sapply(exon_labels, identical, exon_labels[[1]]))

# collapse counts for all samples into one data frame
counts_only <- sapply(counts, function(d) d$count)
counts <- data.frame(exon = as.character(counts[[1]]$exon), counts_only, stringsAsFactors = FALSE)

dim(counts)
head(counts)
tail(counts)
str(counts)


# load truth labels

truth <- read.table(file_truth, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
truth <- truth[, 1:2]

dim(truth)
head(truth)
tail(truth)
str(truth)

table(truth$ds_status)



##############################
### select first 100 genes ###
##############################

n_sub <- 100

genes_sub <- truth$gene[1:n_sub]
genes_sub

genes <- sapply(strsplit(counts$exon, ":"), function(g) g[[1]])


# counts

counts_sub <- counts[genes %in% genes_sub, ]

dim(counts_sub)
head(counts_sub)
tail(counts_sub)


# truth labels

truth_sub <- truth[1:n_sub, ]

dim(truth_sub)
head(truth_sub)
tail(truth_sub)

table(truth_sub$ds_status)


# check ratios of DS to non-DS genes are approximately the same

table(truth$ds_status)[["1"]] / sum(table(truth$ds_status))
table(truth_sub$ds_status)[["1"]] / sum(table(truth_sub$ds_status))



#######################
### save data files ###
#######################

# save in inst/extdata/ directory, so can access with system.file() in vignette

SAVE_DIR <- "../inst/extdata"

write.table(counts_sub, file = file.path(SAVE_DIR, "counts.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(truth_sub, file = file.path(SAVE_DIR, "truth.txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)


