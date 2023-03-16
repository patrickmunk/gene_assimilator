# The program needs the uclust algorithm, but for licencing reasons,
# We cannot distribute the software with this program

rPackageList = c("R.utils")

webUsearchUrl = "https://drive5.com/downloads/usearch11.0.667_win32.exe.gz"

# Check if the software is already installed
#usearchExePath = Sys.which("ipconfig") %>% as.character()

# If not, check if the software is already in the dependencies/ directory
#if (nchar(usearchExePath) == 1) {
#  alreadyDownloaded = sum(file.exists("dependencies/usearch", "dependencies/usearch.exe")) > 0
#}

#if (alreadyDownloaded == TRUE) {
#  usearchExePath = "dependencies/usearch"
#} else {
#  download.file(webUsearchUrl, "dependencies/test", mode = "wb")
#  install.packages("R.utils")
#  R.utils::gunzip("dependencies/test", "dependencies/test.exe")
#}

# If not, then download the binary and put in there
#https://drive5.com/downloads/usearch11.0.667_i86linux32.exe
#https://drive5.com/downloads/usearch11.0.667_i86linux32.gz

# Ensure required R packages
InstallRequiredRPackages = function(required_packages) {
  new_packages = required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages)
}

# Ensure that our dependency is installed
EnsureDependency = function(web_url, dependency_name) {
  # First check if its in system path
  dependency_path = Sys.which(dependency_name) %>% as.character()
  if (nchar(dependency_path) > 0) {
    paste("Is already in sys path? Yes here: ", dependency_path)
    return(dependency_path)
  }
  if (nchar(dependency_path) == 0) {
    # Check if potential binary is already in dependency dir
    potential_binary_files = file.path("dependencies", 
                                       paste(dependency_name, c("", ".exe"), sep = ""))
    is_downloaded = sum(file.exists(potential_binary_files)) > 0
    print(paste("Is downloaded?", is_downloaded))
  } 
  if (is_downloaded == F) {
    download.file(webUsearchUrl, "dependencies/temp.gz", mode = "wb")
    # Install R package dependency for gzip
    InstallRequiredRPackages("R.utils")
    out_binary_name = file.path("dependencies", dependency_name)
    # Add exe extension on Windows systems
    if (Sys.info()[1] == "Windows") {
      out_binary_name = paste(out_binary_name, ".exe", sep ="")
    }
    R.utils::gunzip("dependencies/temp.gz", file.path(out_binary_name))
  }
}

EnsureDependency(webUsearchUrl, "usearch")
