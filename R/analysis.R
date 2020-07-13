#' Trims reads/fragments in fastqc files
#'
#' This function takes a fastqc sequence and trims the reads returning a fastqc file with the trimmed reads.
#'
#'
#' @param verbose Enables progress messages. Default False.
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param quality Min quality of reads to be trimmed.
#' @param first_base First base to keep. Trims everything before this base.
#' @param last_base Last base to keep. Trims everything after this base.
#' @export

trim_fragments=function(bin_path="~/tools/bamUtil/bam",quality=33,first_base="",last_base="",file_R1="",file_R2=""){
  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(file_R1)

  if (!file_R2==""){
  sample_name=intersect_sample_name(file_path=file_R1,file_path2=file_R2)
  output_dir=paste0(output_dir,sep,sample_name,"_trimmed_",first_base,"-",last_base)
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
    fun=function(x,y){
      if(verbose){
        print(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",x,"-o",paste0(get_sample_name(x),"_",y,".fastq.gz")))
      }
      system(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",x,"-o",paste0(get_sample_name(x),"_",y,".fastq.gz")))
    }
    mapply(x=c(file_R1,file_R2),y=c(1,2),FUN=fun)

  }else{
  output_dir=paste0(output_dir,sep,sample_name,"_trimmed_",first_base,"-",last_base)
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  if(verbose){
    print(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",file_R1,"-o",paste0(get_sample_name(file_R1),".fastq.gz")))
  }
  system(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",file_R1,"-o",paste0(get_sample_name(file_R1),".fastq.gz")))
  }
}
