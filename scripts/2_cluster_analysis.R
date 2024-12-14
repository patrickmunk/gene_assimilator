# Run uclust and process the results
require("tidyverse")

# System call to uclust
dir.create(file.path(".", "usearch_clustering"), showWarnings = FALSE)
dir.create(file.path(opt$outputdir, "msa"), showWarnings = FALSE)

#shell("scripts/run_uclust.sh")
# Make this more elegant. Distribute executable? Cross-OS? Not hardcoded params

# The usearch to use. Either one auto-installed in user directory or one in user path
#usearch.path = "dependencies/usearch"
usearch.path = "usearch"

# Get unique sequences from the input sequences
# Construct arguments for usearch fastx_uniques
us.arg.in.path = file.path(opt$outputdir, mergedDBDirName, "combined_database.fa")
us.arg.out.path = file.path(opt$outputdir, mergedDBDirName, "combined_database.uniq.fa")
us.arg.out.tbl = file.path(opt$outputdir, mergedDBDirName, "combined_database.uniq.tsv")
us.arg.threads = 8 # Shared for all usearch calls

us.combined.args = c("-fastx_uniques", us.arg.in.path, "-fastaout", us.arg.out.path, "-threads", us.arg.threads, "-tabbedout", us.arg.out.tbl, "-strand both") # better names?

system2(usearch.path, us.combined.args)


# Save another version of unique genes with user-specified gene abbreviations
userFileOutputName = paste(opt$prefix, ".fa", sep = "")
userFastaPrefix = paste(tolower(opt$prefix), "_", sep = "")
us2.arg.out.path = file.path(opt$outputdir, mergedDBDirName, userFileOutputName)
us2.arg.out.tbl = file.path(opt$outputdir, mergedDBDirName, paste(userFastaPrefix, "uniq.tsv", sep = ""))
#us2.combined.args = c("-fastx_uniques", us.arg.in.path, "-fastaout", us2.arg.out.path, "-relabel", userFastaPrefix, "-threads", us.arg.threads) # better names?
us2.combined.args = c("-fastx_uniques", us.arg.in.path, "-fastaout", us2.arg.out.path, "-relabel", userFastaPrefix, "-threads", us.arg.threads, "-tabbedout", us2.arg.out.tbl) # better names?
system2(usearch.path, us2.combined.args)


# This is perhaps not needed now: we can use one of usearch outputs instead
usearchIdenticalClust = read_delim(file = us2.arg.out.tbl, col_names = F) %>%
  select(1,2)
colnames(usearchIdenticalClust) = c("shortname", "userGeneName")

# Now cluster each remaining user-named gene at 90% ID/cov
us3.arg.in.path = us2.arg.out.path
us3.arg.out.path = file.path(opt$outputdir, mergedDBDirName, paste(opt$prefix, ".nr90.fa", sep = ""))
us3.arg.out2.path =  file.path(opt$outputdir, mergedDBDirName, paste(opt$prefix, ".uc90.tsv", sep = ""))
us3.arg.out3.path = file.path(opt$outputdir, "msa", paste(opt$prefix, ".msa_", sep = ""))
us.combined.args = c("-cluster_fast", us3.arg.in.path, "-id 0.9 -query_cov 0.9 -target_cov 0.9 -centroids", us3.arg.out.path, "-uc ", us3.arg.out2.path, "-threads", us.arg.threads, "-msaout", us3.arg.out3.path)
system2(usearch.path, us.combined.args)

# Read in Usearch cluster data table
uclustData = read_delim(file = us3.arg.out2.path, col_names = F)
colnames(uclustData) = c("linetype", "clustnum", "len", "idpercent", "notsure", 
                         "notsure2", "represent_len", "notsure3", "userGeneName", 
                         "represent_seq")

uclustData = uclustData %>% 
  mutate(chosenSeq = case_when(represent_seq == "*" ~ userGeneName,
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
overviewClustFileSummaryPath = file.path(opt$outputdir, overviewDirName, paste(opt$prefix, "cluster_summary.tsv", sep = "_"))
write_tsv(uclustSummary, overviewClustFileSummaryPath)

# Make a master table with gene info
masterTable = usearchIdenticalClust %>%
  unique() %>%
  left_join(GeneOverviewTable, by = "shortname") %>%
  left_join(select(uclustData, "userGeneName", "chosenSeq") %>% unique) %>%
  arrange(userGeneName, chosenSeq, shortname, fa_name) %>%
  relocate(userGeneName, gene_len)

# Export the Master table
overviewMasterTblPath = file.path(opt$outputdir, overviewDirName, paste(opt$prefix, "master_gene_tbl.tsv", sep = "_"))
write_tsv(masterTable, overviewMasterTblPath)
