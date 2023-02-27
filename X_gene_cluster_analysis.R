# Run uclust and process the results
require("tidyverse")

# System call to uclust
dir.create(file.path(".", "usearch_clustering"), showWarnings = FALSE)

#shell("scripts/run_uclust.sh")
# Make this more elegant. Distribute executable? Cross-OS? Not hardcoded params

# Get unique sequences from the input sequences
# Construct arguments for usearch fastx_uniques
us.arg.in.path = file.path(opt$outputdir, mergedDBDirName, "combined_database.fa")
us.arg.out.path = file.path(opt$outputdir, mergedDBDirName, "combined_database.uniq.fa")
us.arg.threads = 8
#us.combined.args = c("-fastx_uniques", us.arg.in.path, "-fastaout", us.arg.out.path, "-sizeout -relabel Uniq -threads", us.arg.threads)
us.combined.args = c("-fastx_uniques", us.arg.in.path, "-fastaout", us.arg.out.path, "-threads", us.arg.threads) # better names?

#system2("dependencies/usearch", "-fastx_uniques combined_database.fa -fastaout combined_database.uniq.fa -sizeout -relabel Uniq -threads 8")
system2("dependencies/usearch", us.combined.args)

# Save another version of unique genes with user-specified gene abbreviations
userFileOutputName = paste(opt$prefix, ".fa", sep = "")
userFastaPrefix = paste(toupper(opt$prefix), "_", sep = "")
us.arg.out.path = file.path(opt$outputdir, mergedDBDirName, userFileOutputName)
us.combined.args = c("-fastx_uniques", us.arg.in.path, "-fastaout", us.arg.out.path, "-relabel", userFastaPrefix, "-threads", us.arg.threads) # better names?
system2("dependencies/usearch", us.combined.args)

# usearch -fastx_uniques input.fasta -fastaout uniques.fasta -sizeout -relabel Uniq


# Works. test
# TODO: change input sequence so something else
# Change so its not hardcoded too
us.arg.in.path = us.arg.out.path
userFileOutputName = paste(opt$prefix, ".nr90.fa", sep = "")
us.arg.out.path = file.path(opt$outputdir, mergedDBDirName, userFileOutputName)
us.combined.args = c("-cluster_fast", us.arg.in.path, "-id 0.9 -query_cov 0.9 -target_cov 0.9 -centroids", us.arg.out.path, "-uc usearch.uc90.tsv", "-threads", us.arg.threads)

system2("dependencies/usearch", us.combined.args)


uclustData = read_delim(file = "overview/giga_res.uc90.tsv", col_names = F)
colnames(uclustData) = c("linetype", "clustnum", "len", "idpercent", "notsure", 
                         "notsure2", "represent_len", "notsure3", "sequence", 
                         "represent_seq")

uclustData = uclustData %>% 
  mutate(chosenSeq = case_when(represent_seq == "*" ~ sequence,
                                                         represent_seq != "*" ~ represent_seq))


# Summary data per gene cluster
uclustSummary = uclustData %>%
  group_by(clustnum, chosenSeq) %>%
  add_count(clustnum) %>% 
  filter(linetype == "H") %>%
  mutate(idpercent = as.numeric(idpercent)) %>%
  mutate(represent_len = as.numeric(represent_len)) %>%
  group_by(clustnum, chosenSeq, n) %>%
  summarise(maxLen = max(len),
            meanLen = mean(len), 
            minLen = min(len),
            maxID = max(idpercent),
            meanID = mean(idpercent),
            minID = min(idpercent)) %>%
  arrange(-n)

# Export the uclust summary
write_tsv(uclustSummary, "overview/gene_cluster_summary.tsv")

# Test with writing a cluster of genes to a file
uclustData %>% filter(clustnum == 23) %>% pull(sequence)
