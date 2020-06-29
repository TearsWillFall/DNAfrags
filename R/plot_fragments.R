#' Plot DNA fragment length distribution
#'
#' This function takes a BAM file as input and plots the distribution of DNA fragments by their length.
#'
#'
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param verbose Enables progress messages. Default False.
#' @param file Path to BAM file.
#' @export


plot_fragments=function(bin_path="tools/samtools/samtools",file="",verbose=FALSE){
  sample_name=get_sample_name(file)
  if(verbose){
    print(paste(paste0("./",bin_path),"view",file," | awk '{sub(\"^-\", \"\", $9); print $9}' >",paste0(sample_name,"_fragment_length.txt")))
  }
  system(paste(paste0("./",bin_path),"view",file," | awk '{sub(\"^-\", \"\", $9); print $9}' >",paste0(sample_name,"_fragment_length.txt")))
}
