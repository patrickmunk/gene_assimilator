#!/usr/bin/env Rscript

##### GeneAssimiloR #####

# This works even in windows when R and Rscript is in the path
# Rscript GeneAssimilatoR.R --dbdir databases --outputdir panresdb --prefix pan

# Opt-parse arguments from the command line
library("optparse")

option_list = list(
  make_option(c("-d", "--dbdir"), type="character", default=NULL,
              help="Directory with multi-fasta gene database files", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="Directory for GeneAssimilator output", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="Prefix for the new assimilated db. Used in files and genes", metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### MAIN PROGRAM START ###

# Testing. Manual set CLI arguments - Overwrites CLI input
#opt$dbdir = "databases"
#opt$outputdir = "geneAssimilatorOut"
#opt$prefix = "pan"

# Allowed fasta input db file extensions
allowed_exts = c("fa", "fna", "fsa", "fasta")

# Identifies the sequences in a number of fasta input seqs
# Renames the fasta seqs so they conform and save the associations between old and new file to disk
# Outputs rewritten fasta
source("scripts/1_gene_memberships.R")

# Cluster analysis of rewritten files
source("scripts/2_cluster_analysis.R")

# Now run some basic descriptive plots for the results of gene assimilation
source("scripts/3_plotting.R")

print("Finished! Thanks for using GeneAssimilatoR")