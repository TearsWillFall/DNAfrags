#' Filter DNA fragments in BAM files by length and position
#'
#' This function takes a BAM file as input and returns a filtered BAM file
#' by size or/and position, or viceversa.
#'
#'
#' @param file Path to the input file.
#' @param bin_path Path to fastQC executable. Default path tools/samtools/samtools.
#' @param verbose Enables progress messages. Default False.
#' @param min_frag_size Minimum fragment size to keep. Default 1.
#' @param max_frag_size Maximum fragment size to keep. Default 10000000.
#' @param chr Chromosome to keep. Only single chromosomes. Default None.
#' @param start_pos Starting position to search within a chromosome. Default 1
#' @param end_pos Last position to search within a chromosome. Default None
#' @param bed Path to a BED input file with multiple genome positions. Default None
#' @export


filter_fragments=function(bin_path="tools/samtools/samtools",file="",bed="",min_frag_size=0,max_frag_size=100000000,chr="",start_pos=1,end_pos="",verbose=FALSE){
  options(scipen=999)
  position=""
  if(bed==""){
    if (!chr==""){
        position=paste0("chr",chr)
      if (!start_pos==1){
        position=paste0(position,":",start_pos)
        if (!end_pos==""){
          position=paste0(position,"-",end_pos)
          }

      }
      if (verbose){
        print(paste("./",bin_path,"view -h",file,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=",max_frag_size,")' | \ samtools view -b >",paste0(get_sample_name(file),"_",min_frag_size,"-",max_frag_size,"_",position,".bam")))
      }
      system(paste("./",bin_path,"view -h",file,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=",max_frag_size,")' | \ samtools view -b >",paste0(get_sample_name(file),"_",min_frag_size,"-",max_frag_size,"_",position,".bam")))
    }else{
      if (verbose){
        print(paste("./",bin_path,"view -h",file,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=",max_frag_size,")' | \ samtools view -b >",paste0(get_sample_name(file),"_",min_frag_size,"-",max_frag_size,".bam")))

      }
      system(paste("./",bin_path,"view -h",file,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=",max_frag_size,")' | \ samtools view -b >",paste0(get_sample_name(file),"_",min_frag_size,"-",max_frag_size,".bam")))
    }

  }else{
    if (verbose){
      print(paste("./",bin_path,"view -h -L",bed,file," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=",max_frag_size,")' | \ samtools view -b >",paste0(get_sample_name(file),"_",min_frag_size,"-",max_frag_size,"_multiple-pos.bam")))

    }
    system(paste("./",bin_path,"view -h -L",bed,file," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=",max_frag_size,")' | \ samtools view -b >",paste0(get_sample_name(file),"_",min_frag_size,"-",max_frag_size,"_multiple-pos.bam")))
  }


  }



}
