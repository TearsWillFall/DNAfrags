bam="/home/osvaldas/Work/All/SRR11742859_SORTED.RMDUP.SORTED.BAM/SRR11742859.SORTED.RMDUP.SORTED.bam"
get_sample_name=function(file_path=""){
  sample_name=unlist(strsplit(basename(file_path),"\\."))[1]
  return(sample_name)
}

plot_fragments=function(path_bin="tools/samtools/samtools",bam=""){
  sample_name=get_sample_name(bam)
  system(paste(paste0("./",path_bin),"view | awk '{sub(\"^-\", \"\", $9);}'|sort | uniq -c >",paste0(sample_name,"_fragment_length.txt")))
}
