#' Plot DNA fragment length distribution
#'
#' This function takes a file with the fragment length as input and plots the length distribution of DNA fragments.
#' The input file must not have headers and must contain the length of a single fragment per line.
#'
#'
#' @param verbose Enables progress messages. Default False.
#' @param min_frag_length Minimum fragment length to plot. Default 2.
#' @param max_frag_length Max fragment length to plot. When not provided it uses the rounded value of MEDIAN + DEVIATIONS*MEDIAN ABSOLUTE DEVIATION as cut point.
#' @param deviations MEDIAN + DEVIATIONS*MEDIAN ABSOLUTE DEVIATION. Default 10
#' @param width_span Window of width span at which a peak is greater than all other elements around it. Default 3
#' @param min_frgl_maximum Minimum fragment length at which to plot maximums peaks. Default 2
#' @param min_maximum_distance Minimum distance between local maximum peaks to plot. Default 10
#' @param max_maximum_distance Maximum distance between local maximum peaks  to plot. Default 12
#' @param vline Fragment length where to draw a vertical line. Default none
#' @param file Path to file with fragment length.
#' @export



plot_fragments_length <- function(file = "", verbose = FALSE, min_frag_length = 2, max_frag_length = "", deviations = 10, width_span = 3, min_frgl_maximum = 2, min_maximum_distance = 10, max_maximum_distance = 11, vline = "") {
  options(scipen = 999, warn = -1)


  ### TODO add verbose
  sample_name <- ULPwgs::get_sample_name(gsub("_fragment_length*", "", file))

  data <- read.table(file)


  med <- median(data$V1)
  mads <- mad(data$V1)
  med.mad <- med + mads * deviations

  if (max_frag_length == "") {
    max_frag_length <- round(med.mad)
  }

  data.cnt <- plyr:::count(data)
  sub <- data.frame(frags = data.cnt[data.cnt$V1 >= min_frag_length & data.cnt$V1 <= max_frag_length, ]$V1, freq = data.cnt[data.cnt$V1 >= min_frag_length & data.cnt$V1 <= max_frag_length, ]$freq)
  cnt <- sub
  local_maximums <- cnt[ggpmisc:::find_peaks(cnt$freq, span = width_span), ]


  filtered_local_maximums <- local_maximums[local_maximums$frags >= min_frgl_maximum & local_maximums$frags <= data.cnt[data.cnt$freq == max(data.cnt[!data.cnt$V1 == 0, ]$freq), ]$V1, ]
  maximums_distance <- as.data.frame(t(combn(filtered_local_maximums$frags, 2)))
  maximums_distance <- cbind(maximums_distance, dif = maximums_distance$V2 - maximums_distance$V1)
  sorted_maximums <- maximums_distance[maximums_distance$dif >= min_maximum_distance & maximums_distance$dif <= max_maximum_distance, ]
  tmp_maximums <- sorted_maximums

  best_solution <- c()
  while (!dim(tmp_maximums)[1] == 0) {
    solution <- c()
    solution <- append(solution, list(data.frame(start = tmp_maximums[1, ]$V1, end = tmp_maximums[1, ]$V2)))
    tmp_maximums2 <- tmp_maximums
    while (!dim(tmp_maximums2)[1] == 0) {
      if (solution[[length(solution)]]$end == tmp_maximums2[1, ]$V1) {
        solution <- append(solution, list(data.frame(start = tmp_maximums2[1, ]$V1, end = tmp_maximums2[1, ]$V2)))
      }
      tmp_maximums2 <- tmp_maximums2[-1, ]
    }
    if (length(solution) > length(best_solution)) {
      best_solution <- solution
    }
    tmp_maximums <- tmp_maximums[-1, ]
  }


  best_solution <- Reduce(function(x, y) merge(x, y, all = TRUE), best_solution)
  if (!length(best_solution) == 0) {
    best_solution <- best_solution[order(best_solution$start), ]
    best_solution <- cbind(best_solution, dif = best_solution$end - best_solution$start)
    best_solution <- rbind(best_solution, c(best_solution[length(best_solution$dif), ]$end, filtered_local_maximums[length(filtered_local_maximums$frags), ]$frags, filtered_local_maximums[length(filtered_local_maximums$frags), ]$frags - best_solution[length(best_solution$dif), ]$end))
  } else {
    print(paste("No optimal local maximums found for min_frgl_maximum:", min_frgl_maximum, ",min_maximum_distance:", min_maximum_distance, " and max_maximum_distance:", max_maximum_distance))
  }

  ## Generate logs

  log_file <- paste0(sample_name, "_fragment_length_distribution.txt")
  cat(paste(Sys.time(), "\n\n"), file = log_file, append = FALSE)
  cat(paste("## PARAM \n"), file = log_file, append = TRUE)
  param <- data.frame(file = file, verbose = verbose, vline = vline, min_frag_length = min_frag_length, max_frag_length = max_frag_length, deviations = deviations, width_span = width_span, min_frgl_maximum = min_frgl_maximum, min_maximum_distance = min_maximum_distance, max_maximum_distance = max_maximum_distance)
  write.table(param, file = log_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("\n"), file = log_file, append = TRUE)
  info <- data.frame(Median = med, Median_Absolute_Deviation = mads, Mode = data.cnt[data.cnt$freq == max(data.cnt$freq), ]$V1, Min_Fragment_Size = min(data), Max_Fragment_Size = max(data), Mean = mean(as.numeric(data$V1)), Standard_Deviation = sd(as.numeric(data$V1)))
  cat(paste("## STATS \n"), file = log_file, append = TRUE)
  write.table(info, file = log_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("\n"), file = log_file, append = TRUE)
  cat(paste("## PLOT_DATA \n"), file = log_file, append = TRUE)
  write.table(cnt, file = log_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("\n"), file = log_file, append = TRUE)
  cat(paste("## ALL_MAXIMUMS \n"), file = log_file, append = TRUE)
  write.table(local_maximums, file = log_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(paste("\n"), file = log_file, append = TRUE)
  cat(paste("## BEST_LOCAL_MAXIMUMS \n"), file = log_file, append = TRUE)
  write.table(best_solution, file = log_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

  ## Generate plot

  pdf(file = paste0(sample_name, "_fragment_length_distribution.pdf"))
  p <- ggplot2:::ggplot(cnt, ggplot2::aes(x = frags, y = freq)) +
    ggplot2::geom_line(size = 2) +
    ggplot2::scale_x_continuous(breaks = seq(min_frag_length, max_frag_length, 30)) +
    ggplot2::geom_vline(data = local_maximums[local_maximums$frags %in% best_solution$start, ], ggplot2::aes(xintercept = frags), lty = "dashed", size = 0.8) +
    ggplot2::ggtitle(paste("Fragment length distribution for sample", sample_name)) +
    ggplot2::xlab("Fragment length (Pb)") +
    ggplot2::ylab("Counts") +
    ggplot2::theme_classic()
  p <- p + ggplot2::geom_vline(xintercept = max(best_solution$end), col = "red", lty = "dashed", size = 1.2) + ggplot2::annotate(geom = "text", label = max(best_solution$end), x = max(best_solution$end), y = 0.75 * local_maximums[local_maximums$frags == max(best_solution$end), ]$freq, hjust = -1, col = "red")
  if (!vline == "") {
    p <- p + ggplot2::geom_vline(xintercept = vline, col = "green", lty = "dashed", size = 1.2)
  }
  print(p)
  dev.off()
}


#' Extract fragment length from a BAM file.
#'
#' This function takes a BAM file as input and outputs the fragment length of all reads in a TXT file.
#'
#'
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param verbose Enables progress messages. Default False.
#' @param remove_unmapped Removes unmapped reads and/or mates. Only in BAM files. Default TRUE
#' @param bam Path to BAM file.
#' @param threads Number of threads to use.
#' @param bed Path to BED file
#' @param tmp_dir Path to tmp file dir
#' @export

get_fragments_length <- function(bin_path = "tools/samtools/samtools", bam = "", bed="", remove_unmapped = TRUE,
                                 verbose = FALSE, threads = 1,tmp_dir="/tmp") {
  sample_name <- ULPwgs::get_sample_name(bam)
  flags <- ""
  if (remove_unmapped) {
    flags <- "-F 4 -f 2"
  }

  tmp_dir=paste0(" -T ",tmp_dir)
  chr_check <- system(paste(bin_path, " view ",tmp_dir, bam, " | head -n 1 | awk -F \"\t\" '{print $3}'"), intern = TRUE)


  if (bed != "") {
    ref_data <- read.table(bed, comment.char = "")

    if (!grepl("chr", chr_check)) {
      ref_data[, 1] <- gsub("chr", "", ref_data[, 1])
    }
  }




  chrs <- get_chr_names_in_bam(bin_path = bin_path, bam = bam, verbose = verbose)
  if (bed==""){
    dat <- parallel::mclapply(chrs, FUN = function(x) {
      if (verbose) {
        print(paste(bin_path, " view ", tmp_dir, flags, bam, x, paste0(" | awk '{sub(\"^-\", \"\", $9); print $9}' |sort |uniq -c | sort -k2 -V | awk '{print \"", sample_name, "\"\"\\t\"\"", x, "\"\"\\t\"$2\"\\t\"$1}' ")))
      }
      tryCatch(
        {
          dat <- read.table(text = system(paste(bin_path, " view ",tmp_dir, flags, bam, x, paste0(" | awk '{sub(\"^-\", \"\", $9); print $9}' |sort |uniq -c | sort -k2 -V | awk '{print \"", sample_name, "\"\"\\t\"\"", x, "\"\"\\t\"$2\"\\t\"$1}'")), intern = TRUE))
        },
        error = function(e) {
          return(NULL)
        }
      )
    }, mc.cores = threads)

    dat <- dplyr::bind_rows(dat)
    names(dat) <- c("SAMPLE", "REGION", "SIZE", "COUNT")
    dat_tmp <- dat %>%
      dplyr::group_by(SAMPLE, SIZE) %>%
      dplyr::summarise(COUNT = sum(COUNT))
    dat_tmp$REGION <- "GENOME"
    dat <- dplyr::bind_rows(dat, dat_tmp)
    write.table(file = paste0(sample_name, "_fragment_length.txt"), dat, quote = FALSE, col.names = TRUE, row.names = FALSE)
    return(dat)
}else{
  dat <- parallel::mclapply(1:nrow(ref_data), FUN = function(x) {
    if (verbose) {
      print(paste(bin_path, " view ",tmp_dir, flags, bam, paste0(ref_data[x,1],":",ref_data[x,2],"-",ref_data[x,3]), paste0(" | awk '{sub(\"^-\", \"\", $9); print $9}' |sort |uniq -c | sort -k2 -V | awk '{print \"", sample_name, "\"\"\\t\"\"", x, "\"\"\\t\"$2\"\\t\"$1}' ")))
    }
    tryCatch(
      {
        dat <- read.table(text = system(paste(bin_path, " view ",tmp_dir, flags, bam, x, paste0(" | awk '{sub(\"^-\", \"\", $9); print $9}' |sort |uniq -c | sort -k2 -V | awk '{print \"", sample_name, "\"\"\\t\"\"", x, "\"\"\\t\"$2\"\\t\"$1}'")), intern = TRUE))
      },
      error = function(e) {
        return(NULL)
      }
    )
    if (!is.null(ref_data[x,4])){
      dat$GENE=ref_data[x,4]
    }else{
      dat$GENE=""
    }
  }, mc.cores = threads)

}
dat <- dplyr::bind_rows(dat)
names(dat) <- c("SAMPLE", "REGION", "SIZE", "COUNT","GENE")
write.table(file = paste0(sample_name, "_fragment_length.txt"), dat, quote = FALSE, col.names = TRUE, row.names = FALSE)
return(dat)


}


#' Extract fragment length from a BAM file.
#'
#' This function takes a BAM file and a BED file as an input and outputs the fragment length distributions of all the regions in a TXT file.
#'
#'
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param verbose Enables progress messages. Default False.
#' @param bed Path to BED file.
#' @param bam Path to BAM file.
#' @param max_frag_length Maximum fragment length to keep. Default 1000.
#' @param mapq Minimum MapQ of the reads to keep. Default 10.
#' @param start Downstream distance from (start+end)/2 position in bed file
#' @param end Upstream distance from (start+end)/2 position in bed file
#' @param start_bin Downstream distance from (start+end)/2 position in bed file for bin
#' @param end_bin Upstream distance from (start+end)/2 position in bed file for bin
#' @param threads Number of cores. Default 1
#' @param tmp_dir Temporary directory
#' @param output_dir Directory to output results.
#' @param mode Mode in which to output the fragment length. Default 0. Mode 0: Returns fragment length for R1 and R2 reads. Mode 1: Returns fragment legnth for R1 reads. Mode 2: Returns fragment length for R2 reads.
#' @export



get_fragment_length_bed <- function(bin_path = "tools/samtools/samtools", bam = "", bed = "", max_frag_length = 1000, mapq = 10, threads = 1, output_dir = "", verbose = FALSE, mode = 0, start = NULL, end = NULL, start_bin = NULL, end_bin = NULL, tmp_dir = "/tmp") {
  sample_name <- ULPwgs::get_sample_name(bam)

  awk_file_filter <- system.file("shell", "filter.awk", package = "DNAfrags")
  awk_file_stats <- system.file("shell", "stats.awk", package = "DNAfrags")

  chr_check <- system(paste(bin_path, "view", bam, " | head -n 1 | awk -F \"\t\" '{print $3}'"), intern = TRUE)


  if (bed != "") {
    ref_data <- read.table(bed, comment.char = "")

    if (!grepl("chr", chr_check)) {
      ref_data[, 1] <- gsub("chr", "", ref_data[, 1])
    }
  }
  if (!is.null(start) & !is.null(end) & !is.null(end_bin) & !is.null(start_bin)) {
    ctrl_bin <- ref_data
    ctrl_bin$V4 <- as.integer((as.numeric(ctrl_bin$V2) + as.numeric(ctrl_bin$V3)) / 2)
    ctrl_bin$V2 <- ctrl_bin$V4 - start_bin
    ctrl_bin$V3 <- ctrl_bin$V4 + end_bin
    left_flank <- ctrl_bin
    right_flank <- ctrl_bin
    left_flank$V2 <- left_flank$V4 - start
    left_flank$V3 <- left_flank$V4 - start_bin
    right_flank$V3 <- right_flank$V4 + end
    right_flank$V2 <- right_flank$V4 + end_bin
    ref_data <- rbind(left_flank, ctrl_bin, right_flank)
  }

  data <- data.frame(
    chr = ref_data[, 1], r_start = (as.numeric(ref_data[, 2]) + 1),
    r_end = (as.numeric(ref_data[, 3]) + 1)
  ) %>%
    dplyr::mutate(f_start = ifelse((r_start - max_frag_length) < 1, 1,
      r_start - max_frag_length
    ), f_end = (r_end + max_frag_length), r_id = as.character(ref_data[, 4])) %>%
    dplyr::filter(!grepl("_", chr))
  FUN <- function(x, bin_path, bam, mapq, awk_file_filter, awk_file_stats, max_frag_length, verbose, mode) {
    region_data <- t(x)

    position <- paste0(region_data[1], ":", as.numeric(region_data[4]), "-", as.numeric(region_data[5]))

    if (mode == 0) {
      func <- paste0(
        "{ ", bin_path, " view ", bam, " -T ", tmp_dir, " -f 99 ", position, " | awk -v MIN_MAPQ=", mapq,
        " -v MAX_FRAGMENT_LEN=", max_frag_length, " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]),
        " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", region_data[6], " -f ", awk_file_filter, " ; ", bin_path, " view ", bam, " -T ", tmp_dir, " -f 163 ", position, " | awk -v MIN_MAPQ=", mapq,
        " -v MAX_FRAGMENT_LEN=", max_frag_length, " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]),
        " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", region_data[6], " -f ", awk_file_filter, "; } | sort -k9 -n | awk -v MIN_MAPQ=", mapq,
        " -v MAX_FRAGMENT_LEN=", max_frag_length, " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]),
        " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", region_data[6], " -f ", awk_file_stats
      )
      if (verbose) {
        print(func)
      }
      fragment_data <- read.csv(text = system(func, intern = TRUE), header = FALSE, sep = "\t")
    } else if (mode == 1) {
      func <- paste0(
        "{ ", bin_path, " view ", bam, " -T ", tmp_dir, " -f 99 ", position, " | awk -v MIN_MAPQ=", mapq,
        " -v MAX_FRAGMENT_LEN=", max_frag_length, " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]),
        " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", region_data[6], " -f ", awk_file_filter, "; } | sort -k9 -n | awk -v MIN_MAPQ=", mapq,
        " -v MAX_FRAGMENT_LEN=", max_frag_length, " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]),
        " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", region_data[6], " -f ", awk_file_stats
      )
      if (verbose) {
        print(func)
      }
      fragment_data <- read.csv(text = system(func, intern = TRUE), header = FALSE, sep = "\t")
    } else if (mode == 2) {
      func <- paste0(
        "{ ", bin_path, " view ", bam, " -T ", tmp_dir, " -f 163 ", position, " | awk -v MIN_MAPQ=", mapq,
        " -v MAX_FRAGMENT_LEN=", max_frag_length, " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]),
        " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", region_data[6], " -f ", awk_file_filter, "; } | sort -k9 -n | awk -v MIN_MAPQ=", mapq,
        " -v MAX_FRAGMENT_LEN=", max_frag_length, " -v CHR=", region_data[1], " -v R_START=", as.numeric(region_data[2]),
        " -v R_END=", as.numeric(region_data[3]), " -v R_ID=", region_data[6], " -f ", awk_file_stats
      )
      if (verbose) {
        print(func)
      }
      fragment_data <- read.csv(text = system(func, intern = TRUE), header = FALSE, sep = "\t")
    } else {
      (
        stop(paste("Mode:", mode, "is not a valid mode"))
      )
    }

    names(fragment_data) <- c("Region_ID", "Chr", "Region_Start", "Region_End", "Number_of_Reads", "Frag_len_med", "Frag_len_avg", "Frag_len_sd", "Frag_len_distr", "Motif_dist")
    fragment_data$Chr <- as.character(fragment_data$Chr)
    return(fragment_data)
  }

  tictoc::tic("Analysis time: ")
  cl <- parallel::makeCluster(threads)
  df_list <- pbapply::pbapply(
    X = data, 1, FUN = FUN, bin_path = bin_path, bam = bam, mapq = mapq, awk_file_filter = awk_file_filter,
    max_frag_length = max_frag_length, awk_file_stats = awk_file_stats, verbose = verbose, cl = cl, mode = mode
  )
  on.exit(parallel::stopCluster(cl))

  df <- dplyr::bind_rows(df_list) %>% dplyr::arrange(Chr, Region_Start, Region_End)
  df$bed <- ULPwgs::get_sample_name(bed)
  sep <- "/"
  if (output_dir == "") {
    sep <- ""
  }

  out_file <- paste0(output_dir, sep, sample_name, ".fragment_length_regions.txt")

  write.table(df, quote = FALSE, row.names = FALSE, out_file)
  tictoc::toc()
}
