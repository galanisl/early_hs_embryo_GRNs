# Function definitions ----------------------------------------------------

#' Download FASTQ files from ENA
#' 
#' Given a table with sample metadata and a target directory, this function
#' downloads FASTQ files from the FTP addresses in the fastq_ftp column.
#' 
#' @param sdata tibble; A tibble with sample metadata. It must contain a column
#' with the library layout (SINGLE or PAIRED), a column with the FTP address
#' from which the FASTQ files are downloaded and a column of sample IDs.
#' @param target_dir character; The relative path to the directory where the 
#' FASTQ files will be downloaded.
#' 
#' @return The same value as the last system call.
#' 
#' @author Gregorio Alanis-Lobato \email{gregorio.alanis@crick.ac.uk}
#' 
dw_fastq <- function(sdata, target_dir){
  if(sdata$library_layout[1] == "SINGLE"){
    for(i in 1:nrow(sdata)){
      # wget is called with options -nv (no verbose, just print error messages and 
      # basic info), -nc (no clobber, only download new files), -w 5 (wait 5s between
      # downloads), -L and -P (download to target directory)
      system(paste0("wget -L ftp://", sdata$fastq_ftp[i], " -nv -nc -w 5 -P ", target_dir))
    }
  }else{
    for(i in 1:nrow(sdata)){
      fastq <- strsplit(sdata$fastq_ftp[i], ";")[[1]]
      system(paste0("wget -L ftp://", fastq[1], " -nv -nc -w 5 -P ", target_dir))
      system(paste0("wget -L ftp://", fastq[2], " -nv -nc -w 5 -P ", target_dir))
    }
  }
}

#' Merge FASTQ files from the same library
#' 
#' Given a table with sample metadata and a target directory, this function
#' merges FASTQ files from the same library based on the sample_accession column.
#' 
#' @param sdata tibble; A tibble with sample metadata. It must contain a column
#' with the library layout (SINGLE or PAIRED), a column with the FTP address
#' from which the FASTQ files are downloaded and a column of sample IDs.
#' @param fastq_dir character; The relative path to the directory where the 
#' FASTQ files were downloaded.
#' @param target_dir character; The relative path to the directory where the 
#' merged FASTQ files will be created.
#' 
#' @return The same value as the last system call.
#' 
#' @author Gregorio Alanis-Lobato \email{gregorio.alanis@crick.ac.uk}
#' 
merge_fastq <- function(sdata, fastq_dir, target_dir){
  sam_ids <- unique(sdata$sample_accession)
  if(sdata$library_layout[1] == "SINGLE"){
    fastq_files <- paste0(fastq_dir, "/", sapply(strsplit(sdata$fastq_ftp, "/"), 
                          function(x) x[length(x)]))
    for(i in seq_along(sam_ids)){
      fq <- fastq_files[sdata$sample_accession == sam_ids[i]]
      fq <- paste(fq, collapse = " ")
      system(paste0("cat ", fq, " > ", 
                    target_dir, "/", sam_ids[i], ".fastq.gz"))
    }
  }else{
    end1 <- sapply(strsplit(sdata$fastq_ftp, ";"), `[`, 1)
    end1 <- paste0(fastq_dir, "/", sapply(strsplit(end1, "/"), 
                                          function(x) x[length(x)]))
    end2 <- sapply(strsplit(sdata$fastq_ftp, ";"), `[`, 2)
    end2 <- paste0(fastq_dir, "/", sapply(strsplit(end2, "/"), 
                                          function(x) x[length(x)]))
    for(i in seq_along(sam_ids)){
      fq_end1 <- paste(end1[sdata$sample_accession == sam_ids[i]], collapse = " ")
      fq_end2 <- paste(end2[sdata$sample_accession == sam_ids[i]], collapse = " ")
      
      system(paste0("cat ", fq_end1, " > ", 
                    target_dir, "/", sam_ids[i], "_R1.fastq.gz"))
      system(paste0("cat ", fq_end2, " > ", 
                    target_dir, "/", sam_ids[i], "_R2.fastq.gz"))
    }
  }
}


# Argument parsing --------------------------------------------------------

library(optparse)

option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL, 
              help = paste("A CSV file with the sample metadata. The following columns are mandatory:",
                           "\t\t\tlibrary_layout: SINGLE or PAIRED.",
                           "\t\t\tfastq_ftp: FTP address to FASTQ file(s).",
                           "\t\t\tsample_accession: Sample IDs that allow to merge FASTQ files from the same library.",
                           sep = "\n")),
  make_option(c("-d", "--dir"), type = "character", default = "./fastq", 
              help = "The destination directory for the FASTQ files [default = fastq].")
)

opt_parser <-  OptionParser(usage = "Rscript dw_fastq.R [options]", 
                            option_list=option_list)
opt <-  parse_args(opt_parser)

if(is.null(opt$file)){
  print_help(opt_parser)
  stop("The tab-separated file with sample metadata is mandatory.", call. = FALSE)
}

library(readr)
library(dplyr)

sample_data <- read_csv(opt$file) # The table with the sample metadata
target_path <- opt$dir            # The destination directory for the FASTQ files 


# Download the FASTQ files ------------------------------------------------

print("Downloading FASTQ files...")
# Check if dir exists in the current directory. If not, create it. Then 
# download the FASTQ files.
if(!dir.exists(target_path)){
  dir.create(target_path)  
}
dw_fastq(sample_data, target_path)
print("Done downloading FASTQ files...")


# Merge FASTQ files from the same library ---------------------------------
# Check if merging is necessary
gr_sample_data <- sample_data %>% 
  group_by(sample_accession) %>% 
  summarise(runs = n())

if(any(gr_sample_data$runs > 1)){
  print("Merging FASTQ files by library...")
  # Check if dir exists in the current directory. If not, create it. Then 
  # merge FASTQ files from the same library.
  if(!dir.exists(paste0(target_path, "/merged"))){
    dir.create(paste0(target_path, "/merged"))  
  }
  merge_fastq(sample_data, target_path, paste0(target_path, "/merged"))
  print("Done merging FASTQ files...")
}
