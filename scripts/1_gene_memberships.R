# Read and edit fasta file databases
require("tidyverse")
require("seqinr")

# If user has supplied an exclusion list - read that in
if (class(opt$exclude) == "character") {
  exclusion = read_tsv(opt$exclude) %>% 
    pull(1)
}

# Create overview of sequences in each input database
inputMultiFastas = list.files(opt$dbdir)
userFileExts = tools::file_ext(inputMultiFastas)
#inputMultiFastas = inputMultiFastas[grep(fileExtension, inputMultiFastas)]
inputMultiFastas = inputMultiFastas[userFileExts %in% allowed_exts]
inputMultiFastasPath = paste(opt$dbdir, inputMultiFastas, sep ="/")
fileAbbreviations = substr(inputMultiFastas, 1, 3) %>% toupper()

SelectRelevantFiles = function(search_path, good_file_exts) {
  input_files = list.files(search_path)
  good_files = input_files[input_files %>% tools::file_ext() %in% allowed_exts]
  print(paste("Out of", length(input_files), "files,", length(good_files), "are sound like fasta files and are used."))
  return(good_files)
}

# Hardcoded version
# files2include = SelectRelevantFiles("databases", allowed_exts)
files2include = SelectRelevantFiles(opt$dbdir, allowed_exts)


MakeDirs = function(start_path, make_dirs) {
  # Start by creating the output dir requested by user if it does not exist
  print(paste("creating the directory:", start_path))
  ifelse(!dir.exists(file.path(file.path(start_path))), 
         dir.create(file.path(file.path(start_path))), FALSE)
  # Now also make the required sub-directories
  for (i in make_dirs) {
    dir_i_path = file.path(start_path, i)
    print(paste("creating the directory:", i))
    ifelse(!dir.exists(file.path(file.path(dir_i_path))), 
           dir.create(file.path(file.path(dir_i_path))), FALSE)
  }
}

# Make the required output dirs
overviewDirName = "overview"
renamedDirName = "renamed_dbs"
mergedDBDirName = "merged_dbs"
dirs2make = c(overviewDirName, renamedDirName, mergedDBDirName)

MakeDirs(opt$outputdir, dirs2make)

# Renamed fasta files name addition
renamedFileString = "_rn"
renamedFileExt = "fa"

# Information relevant to the files supplied by the user
userGeneDbs = list()
userGeneDbs$inputMultiFastas = files2include
userGeneDbs$userGeneDbExts = tools::file_ext(userGeneDbs$inputMultiFastas)
userGeneDbs$fileAbbreviations = substr(userGeneDbs$inputMultiFastas, 1, 3) %>% toupper()
userGeneDbs$inputMultiFastaPaths = file.path(opt$dbdir, userGeneDbs$inputMultiFastas)
userGeneDbs$inputFastaBasenames = tools::file_path_sans_ext(inputMultiFastas)
userGeneDbs$outputMultiFastaFiles = paste(userGeneDbs$inputFastaBasenames, renamedFileString, ".", "fa", sep = "")
userGeneDbs$outputMultiFastaPaths = file.path(opt$outputdir, renamedDirName, userGeneDbs$outputMultiFastaFiles)

# Function that creates overview table from collection of fasta files
# Version that also writes a new renamed version of each file
MakeOverviewTable = function(input_fasta_files, db_shortnames, output_fasta_files, database_names) {
  geneOverview = data.frame()
  for (i in 1:length(input_fasta_files)) {
    print(paste("Reading", input_fasta_files[i]))
    fasta_i = read.fasta(input_fasta_files[i],
                         whole.header = T,
                         as.string = T,
                         forceDNAtolower = F)
    # Look for genes to exclude if user has supplied list of exclusion
    if (class(opt$exclude) == "character") {
      fasta_i = fasta_i[!names(fasta_i) %in% exclusion]
    }
    # Make dataset
    fasta_i_headers = names(fasta_i)
    database_i_name = database_names[i]
    fasta_i_data = data.frame(shortname = paste(db_shortnames[i], 1:length(fasta_i_headers), sep = "_"),
                              database = database_i_name,
                              fa_header = fasta_i_headers,
                              gene_len = getLength(fasta_i))
    # Add newly added gene data to the overview table
    geneOverview = rbind(geneOverview, fasta_i_data)
    # Write renamed fasta to file
    print(paste("Writing", output_fasta_files[i]))
    write.fasta(fasta_i, fasta_i_data$shortname, 
                file.out = output_fasta_files[i], 
                nbchar = 60, as.string = T)
    
  }
  geneOverview$fa_name = str_split(geneOverview$fa_header, " ", simplify = T)[,1]
  #geneOverview$fa_descript = str_split(geneOverview$fa_header, " ", simplify = T)[,2]
  return(geneOverview)
}

# Make gene overview and export renamed fasta files
GeneOverviewTable = MakeOverviewTable(input_fasta_files = userGeneDbs$inputMultiFastaPaths, 
                  output_fasta_files = userGeneDbs$outputMultiFastaPaths, 
                  db_shortnames = userGeneDbs$fileAbbreviations, 
                  database_names = userGeneDbs$inputFastaBasenames) %>%
  as_tibble()

write_tsv(x = GeneOverviewTable, 
          file = file.path(opt$outputdir, 
                           overviewDirName, 
                           "input_gene_database_overview.tsv"))

# Also make a single combined file with all the renamed references
# In the future it should be combined with the function above so we dont
# have to read in the fastas 2 times

# New renamed files of individual databases
#renamedFilePaths = paste(paste(inputDir, "_ren", sep = ""), "/", 
#      tools::file_path_sans_ext(inputMultiFastas), 
#      "_ren.fa", sep = "")

# Version that writes capital base letters
#MakeCombinedFasta = function(renamed_fasta_dir, out_fasta_file) {
#  input_fasta_paths = list.files(renamed_fasta_dir)
#  renamedFastas = paste(renamed_fasta_dir, input_fasta_paths, sep = "/")
#  outFile <- file(out_fasta_file, "w")
#  for (i in renamedFastas){
#    x = readLines(i)
#    y = toupper(x)
#    writeLines(y, outFile)
#  }
#  close(outFile)
#}

# Version that writes capital base letters - better version
MakeCombinedFasta = function(renamed_fasta_dir, out_file_path) {
  input_fasta_paths = list.files(renamed_fasta_dir)
  renamedFastas = paste(renamed_fasta_dir, input_fasta_paths, sep = "/")
  outFile <- file(out_file_path, "w")
  for (i in renamedFastas){
    x = readLines(i)
    y = toupper(x)
    writeLines(y, outFile)
  }
  close(outFile)
}

MakeCombinedFasta(file.path(opt$outputdir, renamedDirName), file.path(opt$outputdir, mergedDBDirName, "combined_database.fa"))
