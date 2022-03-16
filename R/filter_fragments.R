#' Filter DNA fragments in BAM files by length and position
#'
#' This function takes a BAM file as input and returns a filtered BAM file
#' by size or/and position, or viceversa.
#'
#'
#' @param bam Path to the input file.
#' @param bin_path Path to fastQC executable. Default path tools/samtools/samtools.
#' @param verbose Enables progress messages. Default False.
#' @param min_frag_size Minimum fragment size to keep. Default 1.
#' @param max_frag_size Maximum fragment size to keep. Default 100000000.
#' @param chr Chromosome to keep. Only single chromosomes. Default None.
#' @param start_pos Starting position to search within a chromosome. Default 1
#' @param end_pos Last position to search within a chromosome. Default None
#' @param bed Path to a BED input file with multiple genome positions. Default None
#' @param threads Number of threads. Default 3
#' @export


filter_fragments=function(bin_path="tools/samtools/samtools",bam="",bed="",min_frag_size=35,
  max_frag_size=180, position="",verbose=FALSE,threads=3){
  options(scipen=999)

  sample_name=ULPwgs::get_sample_name(bam)
  chrs <- get_chr_names_in_bam(bin_path = bin_path, bam = bam, verbose = verbose)
  if (bed=="" & position==""){
      if (verbose) {
        print(paste(bin_path,"view -h -b",bam," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
        min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
        ")'", bin_path, "view -h >",paste0(sample_name,"_",min_frag_size,"_",max_frag_size,".bam")))
      }
        system(paste(bin_path,"view -h -b",bam," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
          min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
          ")'",bin_path, "view -h >",paste0(sample_name,"_",min_frag_size,"_",max_frag_size,".bam")))

    ULPwgs::sort_and_index(bin_path,file=paste0(sample_name,"_",min_frag_size,"-",max_frag_size,".bam"),
    threads=threads,verbose=verbose)

  }else if (bed!="" & position==""){

    chr_check <- system(paste(bin_path, " view ", bam, " | head -n 1 | awk -F \"\t\" '{print $3}'"), intern = TRUE)

    ref_data <- read.table(bed, comment.char = "",stringsAsFactors=FALSE,colClasses="character")

    if (!grepl("chr", chr_check)) {
      ref_data[, 1] <- gsub("chr", "", ref_data[, 1])
    }

    dat <- parallel::mclapply(1:nrow(ref_data), FUN = function(x) {
      if (verbose) {
        print(paste(bin_path,"view -h -b",bam, paste0(ref_data[x,1],":",ref_data[x,2],
        "-",ref_data[x,3])," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
        min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
        ")'", bin_path, "view -h >",paste0(ref_data[x,1],":",ref_data[x,2],"-",ref_data[x,3],".bam")))
            }
      tryCatch(
        {
          dat <- system(paste(bin_path,"view -h -b",bam,paste0(ref_data[x,1],":",ref_data[x,2],
          "-",ref_data[x,3])," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
          min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
          ")' ", bin_path, "view -h >",paste0(sample_name,"_",min_frag_size,
          "-",max_frag_size,"_",ref_data[x,1],":",ref_data[x,2],"-",ref_data[x,3],".bam")))
              },
        error = function(e) {
          return(NULL)
        }
      )
    }, mc.cores = threads)

  }else if (position!=""& bed==""){
      if (verbose) {
        print(paste(bin_path,"view -h",bam,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
        min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
        ")' ", bin_path, "view -h >",paste0(sample_name,"_",min_frag_size,"-",max_frag_size,"_",position,".bam")))
      }
    system(paste(bin_path,"view -h",bam,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
    min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
    ")' ", bin_path, "view -h >",paste0(sample_name,"_",min_frag_size,"-",max_frag_size,"_",position,".bam")))



  }else{
    return("ERROR. Multiple input methods selected")
  }


}


#' Merge BAM files in directory
#'
#' This function takes a BAM file and merges it with others found within a
#' directory. This function is still Work In Progress (WIP).
#'
#' @param out_bam Name of merged BAM.
#' @param bam_dir Path to directory with BAM files to merge.
#' @param verbose Enables progress messages. Default False.
#' @export


merge_bam=function(bin_path="tools/samtools/samtools",out_bam="",bam_dir="",verbose=FALSE){
    if(verbose){
      print(paste(bin_path,"merge",paste0(out_bam,".MERGED.bam"), paste0(bam_dir,"/*.bam")))
    }
    system(paste(bin_path,"merge",paste0(out_bam,".MERGED.bam"), paste0(bam_dir,"/*.bam")))
  }
