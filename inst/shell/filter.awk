#/usr/bin/env awk
{
    bPrint = 1
    if ($5 < MIN_MAPQ) bPrint = 0
    if ($9 > MAX_FRAGMENT_LEN) bPrint = 0
    if ($9 <= 0) bPrint = 0
    if (($4+$9) < R_START) bPrint = 0
    if ($4 > R_END) bPrint = 0
    if (bPrint == 1)
        print R_ID"\t"CHR"\t"R_START"\t"R_END"\t"$1"\t"$3"\t"$4"\t"($4+$9)"\t"$9";
}
