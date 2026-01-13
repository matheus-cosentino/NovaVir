# workflow/scripts/plot_rarefaction.R

# Setup Logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(biomformat)
library(stringr)
library(GUniFrac)
library(vegan)

# --- Inputs & Outputs ---
biom_file <- snakemake@input[["biom"]]
pdf_out   <- snakemake@output[["pdf"]]
table_out <- snakemake@output[["table"]]

# --- Read Data ---
message("Reading BIOM file: ", biom_file)
new <- read_biom(biom_file)
df <- as(biom_data(new), "matrix")
df <- as.data.frame(df)

# Adjust OTU ID (Last col -> First col logic)
df$OTU_id <- row.names(df)
# Reorder: Put OTU_id first, then samples
df_max <- df[, c(ncol(df), 1:(ncol(df)-1))]
row.names(df_max) <- NULL

# --- Taxonomy Parsing ---
META <- observation_metadata(new)

# Dynamic Taxonomy Cleaning
if(!is.null(META)) {
    # Find columns that look like taxonomy
    tax_cols <- grep("taxonomy", names(META), value=TRUE, ignore.case=TRUE)
    
    # If found, clean them
    if(length(tax_cols) > 0) {
        clean_tax <- function(x) {
            x <- str_remove(as.character(x), pattern = ".+ ?\\_") # Remove prefix (k__, p__, etc.)
            x[is.na(x)] <- ""
            return(x)
        }
        
        # Apply cleaning
        META[tax_cols] <- lapply(META[tax_cols], clean_tax)
        
        # Paste columns together separated by ;
        if(length(tax_cols) > 1) {
            taxon <- do.call(paste, c(META[tax_cols], sep=";"))
        } else {
            taxon <- META[[tax_cols]]
        }
        df_max$taxonomy <- taxon
    } else {
        df_max$taxonomy <- "Unclassified"
    }
} else {
    df_max$taxonomy <- "Unclassified"
}

# Save Raw Table
write.table(df_max, file = table_out, sep = "\t", quote = F, row.names = TRUE, col.names = NA)

# --- Normalization & Rarefaction Logic ---
# Prepare numeric matrix for analysis (Remove ID and Taxonomy)
numeric_cols <- names(df_max)[!names(df_max) %in% c("OTU_id", "taxonomy")]

#avoid errors within 1 sample runs
otu_table <- df_max[, numeric_cols, drop = FALSE]
rownames(otu_table) <- df_max$OTU_id

# Convert to numeric
otu_table[] <- lapply(otu_table, as.numeric)
rownames(otu_table) <- df_max$OTU_id

# Clean empty rows
otu_table <- otu_table[rowSums(otu_table) > 0, , drop=FALSE]

# --- Plotting ---
message("Generating Rarefaction Curve...")
pdf(file = pdf_out)

# Calculate Curves
# Transpose: vegan rarecurve expects Rows=Samples
rarefactionCurve <- rarecurve(data.frame(t(otu_table)),
                              step = 20,
                              col = "black",
                              lty = "solid",
                              label = FALSE,
                              xlab = "Total Reads",
                              ylab = "Richness (Families/OTUs)",
                              main = "Rarefaction Curves (All Samples)")

# Labeling logic (Top 21 or all)
labelCutoff <- min(21, ncol(otu_table))

# Calculate slopes/angles for labeling
slope <- vector()
SampleID <- vector()

for (i in seq_along(rarefactionCurve)) {
  curve_i <- rarefactionCurve[[i]]
  richness <- 1000
  
  if (length(curve_i) > 6) {
    y_end <- curve_i[length(curve_i)]
    y_prev <- curve_i[length(curve_i) - 5]
    x_end <- attr(curve_i, "Subsample")[length(curve_i)]
    x_prev <- attr(curve_i, "Subsample")[length(curve_i) - 5]
    
    richness <- (y_end - y_prev) / (x_end - x_prev)
  }
  slope <- c(slope, richness)
  SampleID <- c(SampleID, names(otu_table)[i])
}

curvedf <- data.frame(SampleID, slope)
ordered_vector <- order(as.numeric(curvedf$slope), decreasing = TRUE)

# Add red lines for top subsampled
for (i in 1:labelCutoff) {
  idx <- ordered_vector[i]
  N <- attr(rarefactionCurve[[idx]], "Subsample")
  lines(N, rarefactionCurve[[idx]], col="red")
  text(max(N), max(rarefactionCurve[[idx]]), SampleID[idx], cex=0.6, pos=4)
}

dev.off()
message("Done.")