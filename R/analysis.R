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
#' @param output_dir Path to the output directory.
#' @export



trim_fragments=function(bin_path="tools/fastx_toolkit/bin/fastx_trimmer",quality=33,first_base="",last_base="",file_R1="",file_R2="",output_dir="",verbose=FALSE){
  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=ULPwgs::get_sample_name(file_R1)

  if (!file_R2==""){
  sample_name=ULPwgs::intersect_sample_name(file_path=file_R1,file_path2=file_R2)
  output_dir=paste0(output_dir,sep,sample_name,"_trimmed_",first_base,"-",last_base)
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
    fun=function(x,verbose,bin_path,quality,first_base,last_base,output_dir){
      if(verbose){
        print(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",x,"-o",paste0(output_dir,"/",paste0(ULPwgs::get_sample_name(x),".fastq.gz"))))
      }
      system(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",x,"-o",paste0(output_dir,"/",paste0(ULPwgs::get_sample_name(x),".fastq.gz"))))
    }
  mapply(FUN=fun,x=c(file_R1,file_R2),MoreArgs=list(verbose,bin_path,quality,first_base,last_base,output_dir))

  }else{
  sample_name=unlist(strsplit(sample_name),"_")[1]
  output_dir=paste0(output_dir,sep,,"_trimmed_",first_base,"-",last_base)
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  if(verbose){
    print(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",file_R1,"-o",paste0(output_dir,"/",sample_name,".fastq.gz")))
  }
  system(paste(bin_path,paste0("-Q",quality),"-f",first_base,"-l",last_base,"-z -i",file_R1,"-o",paste0(output_dir,"/",sample_name,".fastq.gz")))
  }
}
