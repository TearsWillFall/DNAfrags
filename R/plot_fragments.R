#' Plot DNA fragment length distribution
#'
#' This function takes a BAM file as input and plots the length distribution of DNA fragments.
#'
#'
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param verbose Enables progress messages. Default False.
#' @param min_frag_length Minimum fragment length to plot. Default 2.
#' @param max_frag_length Max fragment length to plot. When not provided it uses the rounded value of MEDIAN + DEVIATIONS*MEDIAN ABSOLUTE DEVIATION as cut point.
#' @param deviations MEDIAN + DEVIATIONS*MEDIAN ABSOLUTE DEVIATION. Default 10
#' @param width_span Window of width span at which a peak is greater than all other elements around it. Default 3
#' @param min_frgl_maximum Minimum fragment length at which to plot maximums peaks. Default 2
#' @param max_frgl_maximum Maximum fragment length at which to plot maximums peaks. Default 167
#' @param min_maximum_distance Minimum distance between local maximum peaks to plot. Default 10
#' @param max_maximum_distance Maximum distance between local maximum peaks  to plot. Default 12
#' @param file Path to BAM file.
#' @export



plot_fragments=function(bin_path="tools/samtools/samtools",file="",verbose=FALSE,min_frag_length=2,max_frag_length="",deviations=10,width_span=3,min_frgl_maximum=2,max_frgl_maximum=167,min_maximum_distance=10,max_maximum_distance=12){
  options(scipen=999)
  sample_name=get_sample_name(file)
  if(verbose){
    print(paste(paste0("./",bin_path),"view -F 4 -f 2",file," | awk '{sub(\"^-\", \"\", $9); print $9}' >",paste0(sample_name,"_fragment_length.txt")))
  }
  system(paste(paste0("./",bin_path),"view -F 4 -f 2",file," | awk '{sub(\"^-\", \"\", $9); print $9}' >",paste0(sample_name,"_fragment_length.txt")))

  ### TODO add verbose
  pdf(file=paste0(sample_name,".pdf"))
  data=read.table(paste0(sample_name,"_fragment_length.txt"))
  med=median(data$V1)
  mads=mad(data$V1)
  med.mad=med*deviations+mads

  if (max_frag_length==""){
    max_frag_length=round(med.mad)
  }

  sub=data.frame(frags=data[data$V1>=min_frag_length &data$V1<=max_frag_length,])
  cnt=dplyr::count(sub)
  local_maximums=cnt[ggpmisc:::find_peaks(cnt$freq,span=width_span),]
  local_maximums=local_maximums[local_maximums$frags>=min_frgl_maximum & local_maximums$frags<=max_frgl_maximum,]
  maximums_distance=as.data.frame(t(combn(local_maximums$frags,2)))
  maximums_distance=cbind(maximums_distance,dif=maximums_distance$V2-maximums_distance$V1)
  sorted_maximums=maximums_distance[maximums_distance$dif>=min_maximum_distance & maximums_distance$dif<=max_maximum_distance,]
  tmp_maximums=sorted_maximums

  best_solution=c()
  while (!dim(tmp_maximums)[1]==0){
    solution=c()
    solution=append(solution,list(data.frame(pos=row.names(tmp_maximums[1,]),V1=tmp_maximums[1,]$V1,V2=tmp_maximums[1,]$V2)))
    tmp_maximums2=tmp_maximums
    while(!dim(tmp_maximums2)[1]==0){
      if (solution[[length(solution)]]$V2==tmp_maximums2[1,]$V1){
        solution=append(solution,list(data.frame(pos=row.names(tmp_maximums2[1,]),V1=tmp_maximums2[1,]$V1,V2=tmp_maximums2[1,]$V2)))
      }
      tmp_maximums2=tmp_maximums2[-1,]
    }
    if (length(solution)>length(best_solution)){
      best_solution=solution
    }
    tmp_maximums=tmp_maximums[-1,]
  }


  best_solution=Reduce(function(x, y) merge(x, y, all=TRUE), best_solution)
  print(best_solution)
  ggplot2::ggplot(cnt, aes(x =frags,y=freq)) +
  ggplot2::geom_line(size=2) +
  ggplot2::scale_x_continuous(breaks=seq(min_frag_length,max_frag_size,30))+
  ggplot2::geom_vline(data=local_maximums[local_maximums$frags %in% unique(c(best_solution$V1,best_solution$V2,local_maximums[length(local_maximums$frags),])),],aes(xintercept=frags),lty="dashed",size=0.8)+
  ggplot2::geom_vline(xintercept=167,col = "green",lty="dashed",size=1.2)+
  ggplot2::ggtitle("Fragment length distribution") +
  ggplot2::xlab("Fragment length (Pb)") +
  ggplot2::ylab("Counts")+
  ggplot2::theme_classic()
  dev.off()

}
