#!/usr/bin/env Rscript

# --- 1. Load Libraries ---
# Ensure all required packages are loaded.
if (!requireNamespace("argparse", quietly = TRUE)) stop("Package 'argparse' is required.")
if (!requireNamespace("readr", quietly = TRUE)) stop("Package 'readr' is required.")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
if (!requireNamespace("ape", quietly = TRUE)) stop("Package 'ape' is required.")
if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

library(argparse)
library(readr)
library(dplyr)
library(ape)
library(stringr)
library(ggplot2)

# --- 2. Define Functions ---
# Copy the PlotPP function from your Rmd for plotting RdRp motifs.
PlotPP <- function(ps){
  # Bind Local Variables
  segstrt <- segend <- segment <- NULL
  
  # plot variables
  pp.col <- c("#440154", "gray60", "#3b528b", "gray60", "#fde725")
  vbump = 0.85
  
  # Query
  ps.query <- data.frame( segment = c("query"),
                          segstrt = 1,
                          segend = ps$qlen)
  
  # Palmprint
  ps.pp <- data.frame(segment = unlist(ps$order_labels),
                         segstrt = unlist(ps$order_starts),
                         segend = unlist(ps$order_ends),
                         color= unlist(ps$order_colors))
  
  ppPlot <- ggplot() +
    geom_segment(data = ps.query, aes(x = segstrt, xend = segend, y = 0, yend = 0), linetype = 1, size = 4, color = "gray") +
    xlab("query position (aa)") +
    geom_segment(data = ps.pp, aes(x = segstrt, xend = segend, y = 0, yend = 0), linetype = 1, size = 10, color = ps.pp$color) +
    geom_text(data = ps.pp, aes(x = ((segstrt + segend)/2) , y = vbump, label = segment, hjust = "center"), color = ps.pp$color) +
    ggtitle(label = paste0(">", ps$query)) +
    ylim(c(-1, 1)) + xlim(c(0, NA)) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(), legend.position="none",
          panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), plot.background=element_rect(fill = "white", colour = "white"))
  
  return(ppPlot)
}

# --- NEW Function for Single Plot Generation ---
GenerateSinglePlot <- function(duskmatter_path, contig_label, sample_name, output_path) {
    cat("INFO: Generating single RdRp motif PNG...\n")
    
    duskmatter_df <- read.table(duskmatter_path, header = TRUE, sep = "\t")
    
    # Filter for the specific contig label
    fev <- duskmatter_df %>% filter(Label == contig_label)
    
    if (nrow(fev) == 0) {
        stop(paste("Contig label", contig_label, "not found in", duskmatter_path))
    }
    
    # Re-implement the data preparation logic from the Rmd
    mA <- 11; mB <- 14; mC <- 7
    pp.col.canonical <- c("#440154", "gray60", "#3b528b", "gray60", "#fde725")
    fev$query <- fev$Label
    fev$qlen <- nchar(fev$aaseq)
    fev$order_labels <- list(NA); fev$order_colors <- list(NA); fev$order_starts <- list(NA); fev$order_ends <- list(NA)

    ps <- fev[1, ] # Process only the first (and only) row
    if (ps$pssm_ABC == "ABC") {
        pA <- ps$posA; pB <- ps$posB; pC <- ps$posC
        labels <- c("A", "v1", "B", "v2", "C"); colors <- pp.col.canonical
        starts <- c(pA, pA + mA + 1, pB, pB + mB + 1, pC); ends <- c(pA + mA, pB, pB + mB, pC, pC + mC - 1)
    } else if (ps$pssm_ABC == "CAB") {
        pC <- ps$posC; pA <- ps$posA; pB <- ps$posB
        labels <- c("C", "v1", "A", "v2", "B"); colors <- c(pp.col.canonical[5], pp.col.canonical[2], pp.col.canonical[1], pp.col.canonical[4], pp.col.canonical[3])
        starts <- c(pC, pC + mC + 1, pA, pA + mA + 1, pB); ends <- c(pC + mC, pA, pA + mA, pB, pB + mB)
    } else { return(NULL) } # Skip if no valid motif

    fev[1, ]$order_labels <- list(labels); fev[1, ]$order_colors <- list(colors); fev[1, ]$order_starts <- list(starts); fev[1, ]$order_ends <- list(ends)
    
    pp.plot <- PlotPP(fev[1, , drop = FALSE])
    ggsave(filename = output_path, plot = pp.plot, width = 8, height = 6, dpi = 300)
    cat(paste("  - Wrote", output_path, "\n"))
}

# --- 3. Parse Command-Line Arguments ---
parser <- ArgumentParser(description = "Generate dynamic FASTA and PNG outputs for DiscoVir report.")
parser$add_argument("--sample_name", required = TRUE, help = "Sample identifier.")
parser$add_argument("--fasta_path", required = TRUE, help = "Path to the main contigs FASTA file.")
parser$add_argument("--diamond_path", required = TRUE, help = "Path to the main diamond annotation file.")
parser$add_argument("--duskmatter_path", required = TRUE, help = "Path to the Duskmatter RdRp annotation file.")
parser$add_argument("--viral_summary_path", required = TRUE, help = "Path to the viral contigs summary TSV file.")
parser$add_argument("--output_dir", required = TRUE, help = "Directory to save output files.")

args <- parser$parse_args()

# This block makes the script runnable from the command line for debugging,
# but the new Snakemake rules will call the functions directly.
if (sys.nframe() == 0) {
    # --- 4. Read Input Files ---
    cat("INFO: Reading input files...\n")
    fasta <- read.FASTA(args$fasta_path)
    names(fasta) <- str_remove_all(names(fasta), pattern = " .+")
    
    # This is the file that tells us which families to process
    viral_summary_metrics <- read.table(args$viral_summary_path, sep = "\t", header = TRUE, check.names = FALSE)
    family_names <- colnames(viral_summary_metrics)[-1] # Get family names from columns, excluding 'Metric'
    
    # We also need the main diamond file to link contigs to families
    diamond_data <- read.table(args$diamond_path, sep = "\t", fill = TRUE, header = TRUE)
    virus_data <- diamond_data %>% filter(Overall == "Viruses")
    
    duskmatter_df <- read.table(args$duskmatter_path, header = TRUE, sep = "\t")
    
    # --- 5. Generate Per-Family FASTA Files ---
    cat("INFO: Generating per-family FASTA files...\n")
    for (family in family_names) {
        # Find which contigs belong to this family
        qseqids_for_family <- virus_data$qseqid[which(virus_data$Family == family)]
        
        if (length(qseqids_for_family) > 0) {
            index <- which(names(fasta) %in% qseqids_for_family)
            fasta_subset <- fasta[index]
            
            # Define output file path and write the FASTA
            output_fasta_path <- file.path(args$output_dir, paste0(args$sample_name, "_", family, ".fasta"))
            write.FASTA(fasta_subset, file = output_fasta_path)
            cat(paste("  - Wrote", output_fasta_path, "\n"))
        }
    }
    
    # --- 6. Generate Per-Contig RdRp Motif PNGs ---
    cat("INFO: Generating RdRp motif PNG files...\n")
    
    # Re-implement the data preparation logic from the Rmd
    mA <- 11; mB <- 14; mC <- 7
    pp.col.canonical <- c("#440154", "gray60", "#3b528b", "gray60", "#fde725")
    fev <- duskmatter_df
    fev$query <- fev$Label
    fev$qlen <- nchar(fev$aaseq)
    fev$order_labels <- list(NA); fev$order_colors <- list(NA); fev$order_starts <- list(NA); fev$order_ends <- list(NA)
    
    for (i in 1:nrow(fev)) {
      ps <- fev[i, ]
      if (ps$pssm_ABC == "ABC") {
        pA <- ps$posA; pB <- ps$posB; pC <- ps$posC
        labels <- c("A", "v1", "B", "v2", "C"); colors <- pp.col.canonical
        starts <- c(pA, pA + mA + 1, pB, pB + mB + 1, pC); ends <- c(pA + mA, pB, pB + mB, pC, pC + mC - 1)
      } else if (ps$pssm_ABC == "CAB") {
        pC <- ps$posC; pA <- ps$posA; pB <- ps$posB
        labels <- c("C", "v1", "A", "v2", "B"); colors <- c(pp.col.canonical[5], pp.col.canonical[2], pp.col.canonical[1], pp.col.canonical[4], pp.col.canonical[3])
        starts <- c(pC, pC + mC + 1, pA, pA + mA + 1, pB); ends <- c(pC + mC, pA, pA + mA, pB, pB + mB)
      } else { next }
      
      fev[i, ]$order_labels <- list(labels); fev[i, ]$order_colors <- list(colors); fev[i, ]$order_starts <- list(starts); fev[i, ]$order_ends <- list(ends)
      
      # Now plot and save
      pp.plot <- PlotPP(fev[i, , drop = FALSE])
      output_png_path <- file.path(args$output_dir, paste0(args$sample_name, "_", fev$query[i], "RdRp_Motifs.png"))
      ggsave(filename = output_png_path, plot = pp.plot, width = 8, height = 6, dpi = 300)
      cat(paste("  - Wrote", output_png_path, "\n"))
    }
    
    cat("SUCCESS: Dynamic file generation complete.\n")
}
