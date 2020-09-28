#/usr/bin/env awk
{
    fl_count[NR] = $9;
    fl_dist[$9":"] = fl_dist[$9":"]+1;
    mot = substr($10, 1, 4);
    motif_dist[substr($10, 1, 4)":"] = motif_dist[mot":"]+1;
}
END {
    fl_median = 0
    fl_average = 0
    fl_sd = 0
    fl_str_dist=""
    for( fl in fl_dist ) {
        if (fl_str_dist != ""){
        fl_str_dist = fl_str_dist"|"fl fl_dist[fl]
        }
        else{
          fl_str_dist = fl_dist[fl]
        }
    }

    motif_str_dist=""
    for( motif in motif_dist ) {
        if (motif_str_dist != ""){
        motif_str_dist  = motif_str_dist"|"motif motif_dist[motif]
        }
        else{
          motif_str_dist = motif_dist[motif]
        }
    }



    if (NR > 1) {
        if ((NR % 2) == 1) {
            fl_median = fl_count[(NR + 1) / 2];
        } else {
            fl_median = (fl_count[NR / 2] + fl_count[(NR / 2) + 1]) / 2.0;
        }
        fl_sum = 0;
        for( i = 1; i <= length( fl_count ); i++ ) {
            fl_sum += fl_count[i];
        }
        fl_average = fl_sum / NR
        fl_sumsd = 0;
        for( i = 1; i <= length( fl_count ); i++ ) {
            fl_sumsd += (fl_count[i] - fl_average)^2;
        }
        fl_sd = (fl_sumsd /(NR - 1))^(1/2)
    } else {
        if (NR == 1) {
            fl_median = fl_count[1]
            fl_average = fl_count[1]
            fl_sd = 0
        } else {
            fl_median = 0
            fl_average = 0
            fl_sd = 0
        }
    }
    printf( R_ID"\t"CHR"\t"R_START"\t"R_END"\t%d\t%d\t%d\t%d\t%s\n", NR , fl_median, fl_average , fl_sd, fl_str_dist, motif_str_dist);
}
