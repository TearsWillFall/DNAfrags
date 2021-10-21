#' Function to collect chromosome named
#'
#' This function takes a BAM file and collects the chr names from the bam
#' file.
#'
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param bam Path to directory with BAM files to merge.
#' @param verbose Enables progress messages. Default False.
#' @export

get_chr_names_in_bam=function(bin_path="tools/samtools/samtools",bam="",verbose=FALSE){
  if(verbose){
    print(paste0(bin_path," view -H ",bam," | grep @SQ"))
  }
  chr=read.table(system(paste0(bin_path," view -H ",bam," | grep @SQ"),intern=TRUE))
  chr=chr[2]
  chr=unlist(lapply(chr,FUN=function(x){strsplit(x,":")[[1]][2]}))
  return(chr)
}
