# Creates some visual overview of input database overview

# Determine which exact sequences are shared between input databases

seqDBoverlap = masterTable %>%
  # Change between userGeneName and chosenSeq - exact vs id clusters 
  select(userGeneName, database) %>% 
  #select(chosenSeq, database) %>% 
  table %>%
  as_tibble()

# Identify the number of times each gene is found
genePrevalence = masterTable %>% 
  select(userGeneName) %>% 
  #select(chosenSeq) %>% 
  table %>% 
  as_tibble() %>%
  arrange(-n)
colnames(genePrevalence) = c("userGeneName", "prevalence")

genePrevalenceOrder = genePrevalence$userGeneName

genePrevalence = genePrevalence %>%
  mutate(userGeneNameFact = factor(userGeneName, levels = genePrevalenceOrder))


#seqDBoverlap %>% group_by(userGeneName) %>% summarise(dbcount = sum(n)) %>% ggplot(aes(userGeneName, dbcount)) + geom_tile()

# Heatmap of gene belonging
seqDBoverlap %>% 
  group_by(userGeneName) %>% 
  left_join(genePrevalence, by = "userGeneName") %>% 
  ggplot(aes(userGeneNameFact, database, fill = as.factor(n))) + 
  geom_tile() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_brewer()
ggsave(file.path(opt$outputdir, overviewDirName, "exact_seq_db_overlap.png"),
       width = 14, height = 6)

# Generalize to a 