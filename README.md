# Gene Assimilator

We are often in a situation where we have several collections of genes that we want to compare and combine in a new gene collection.
This project tries to adress that by collecting data on and clustering gene catalogs and outputting combined catalogs and homology-reduced versions.

In order to use the program, you need to gather a collection of multi-fasta files in a directory.
These will be the ones analyzed, clustered and merged into an output database with annotation files.
The three first letter in each input fasta file will be used as a prefix for some of the output.

## Dependencies
You need to have usearch in your path (https://www.drive5.com/usearch/download.html).
We also need a couple of R packages: optparse, tidyverse, seqinr.
The pipeline was found to both work in R 4.1.1 on Windows 10 and R 4.2.0 on Linux (CentOS).

## How to run the program
Rscript GeneAssimilatoR.R --dbdir <myInputDir> --outputdir <myOutputDir> --prefix <myOutputPrefix>

So something like this should work with the example files included:
Rscript GeneAssimilatoR.R --dbdir example_databases/ --outputdir testOutDir/ --prefix test
