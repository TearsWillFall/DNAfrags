#!/bin/bash
#$ -l h_rt=48:00:0
#EXTRA -t 1-200
#EXTRA -pe smp 4


###################################################################################
# Submit Myryad
###################################################################################
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=Germline_PT23
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=Plasma_PT23
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=E23_Plasma_BL_d0_BSwg_R2
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=E36_Plasma_BL_d0_BSwg_R2
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=E123_Plasma_PD_d0_BSwg_R2
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=APC21_Plasma_BL_d0_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=APC21_Plasma_PD_d108_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=E23_Plasma_BL_d0_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=E36_Plasma_BL_d0_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=E123_Plasma_PD_d999_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=H1_Plasma_HV_d0_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=H1_Plasma_HV_d0_BSlp_R2
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=V5322_Plasma_BL_d0_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=V5322_Plasma_PD_d10_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=15-638
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=15-690
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA34_334
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA34_334hc_WGSsf_80_01
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA34_334hcfull
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA43_434
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA43_434hc_WGSsf_80_01
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA43_434hcfull
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA43_437
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA43_437hc_WGSsf_80_01
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA43_437hcfull
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA79_36
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA79_36hc_WGSsf_80_01
#Submit: qsub -N fl -l h_rt=48:00:00 -l mem=4G -l tmpfs=10G -wd /lustre/scratch/scratch/regmpcr/fragment_length/data/LOG /lustre/scratch/scratch/regmpcr/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=CA79_36hcfull


###################################################################################
# Submit Gamble
###################################################################################
#Submit: qsub -N fl -l h_rt=48:00:00 -wd /SAN/colcc/stratosphere/fragment_length/data/LOG /SAN/colcc/stratosphere/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=E36_Plasma_BL_d0_BSlp_R1
#Submit: qsub -N fl -l h_rt=48:00:00 -wd /SAN/colcc/stratosphere/fragment_length/data/LOG /SAN/colcc/stratosphere/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=Germline_PT60
#Submit: qsub -N fl -l h_rt=48:00:00 -wd /SAN/colcc/stratosphere/fragment_length/data/LOG /SAN/colcc/stratosphere/bash-qsub/pipeline_fl.sh -e -r=hg19_W500kbs -s=Germline_PT78


###################################################################################
# Set analysis parameters
###################################################################################
SERVER="GAMBLE"
SECURITY="NA:NA"
LIST_SAMPLES=""
SAMPLE="MISSING"
REGIONS="MISSING"
EXECUTE="NO"
MAX_FRAGMENT_LEN=1000
MIN_MAPQ=10
LIST_STEPS_PRIN="START S1.01 END"
LIST_STEPS_EXEC=""
N_CORES=4
for i in "$@";do
	case $i in
		-e|--exec) EXECUTE="YES";;
		-c|--cluster-name) SERVER="${i#*=}";;
		-u=*|--user=*) SECURITY="${i#*=}";;
		-s=*|--sample=*) SAMPLE="${i#*=}";;
		-r=*|--region=*) REGIONS="${i#*=}";;
		-p=*|--process-steps=*) LIST_STEPS_PRIN="${i#*=}";;
		-l=*|--max-fragment-len=*) MAX_FRAGMENT_LEN="${i#*=}";;
		-q=*|--min-mapq=*) MIN_MAPQ="${i#*=}";;
	*) echo "UNKNOWN PARAM: $i";;
	esac
done
if [[ $EXECUTE =~ "YES" ]];then LIST_STEPS_EXEC=$LIST_STEPS_PRIN; else LIST_STEPS_EXEC="";fi


###################################################################################
# Load requred modules in MYRIAD
###################################################################################
if [[ $SERVER =~ "MYRIAD" ]];then
	source /etc/profile.d/modules.sh
	module unload compilers/intel/2018/update3
	module load gcc-libs/4.9.2
	module load compilers/gnu/4.9.2
	module load openblas/0.2.14/gnu-4.9.2
	module load fftw/3.3.6-pl2/gnu-4.9.2
	module load ghostscript/9.19/gnu-4.9.2
	module load texinfo/6.6/gnu-4.9.2
	module load texlive/2019
	module load gsl/2.4/gnu-4.9.2
	module load hdf/5-1.8.15/gnu-4.9.2
	module load netcdf/4.3.3.1/gnu-4.9.2
	module load jags/4.2.0/gnu.4.9.2-openblas
	module load screen/4.2.1
	module load xorg-utils/X11R7.7
	module load java/1.8.0_92
	module load julia/1.1.0
	module load python/3.7.4
	module load r/3.6.0-openblas/gnu-4.9.2
	module load git/2.19.1
	module load delly/0.7.8-bindist
	module load fastqc/0.11.8
	module load samtools/1.9/gnu-4.9.2
	module load gatk/4.0.8.0
	module load picard-tools/2.18.9
	module load skewer/0.2.2
	module load htslib/1.7
	module load fastqc/0.11.8
	module load bwa/0.7.12/gnu-4.9.2
	module load picard-tools/2.18.9
	module load bamtools/2.4.0/gnu-4.9.2
	module load bedtools/2.25.0
	#module load platypus/3e72641
	module load bcftools/2.1/gnu-4.9.2
	module load gatk/4.0.8.0
	cd /home/regmpcr/Scratch/bash-qsub
	export PYTHONPATH=/home/regmpcr/Scratch/python/3.7:$PYTHONPATH
	export R_HOME=/shared/ucl/apps/R/R-3.6.0-OpenBLAS/bin/R
	export R_LIB_PATHS=/home/regmpcr/Scratch/R/3.6:$R_LIB_PATHS
	export R_LIBS_USER=/home/regmpcr/Scratch/R/3.6:$R_LIBS_USER
	export R_MAX_NUM_DLLS=614
	export PATH=/shared/ucl/apps/fastqc/0.11.8/FastQC/fastqc:$PATH
	PATH_DATA=/lustre/scratch/scratch/regmpcr/fragment_length/data/
	genome=/lustre/scratch/scratch/regmpcr/Genomes/hs37d5.fa
fi


###################################################################################
# Folders vars
###################################################################################
PATH_ALIGN_BAM=$PATH_DATA/BAM
PATH_REPOR_TXT=$PATH_DATA/TXT
PATH_REPOR_BED=$PATH_DATA/BED
AWK_FILE_FILTE=$PATH_DATA/../source/filter.awk
AWK_FILE_STATS=$PATH_DATA/../source/stats.awk
AWK_FILE_DETAI=$PATH_DATA/../source/detail.awk
mkdir -p $PATH_ALIGN_BAM
mkdir -p $PATH_REPOR_TXT
mkdir -p $PATH_REPOR_BED
cd $PATH_DATA
BAM_FILE=$PATH_ALIGN_BAM/$SAMPLE.RMDUP.SORTED.bam
BED_FILE=$PATH_REPOR_BED/$REGIONS.bed
echo "$(date) Start analysis"
echo "$(date) Analysis steps: $LIST_STEPS_PRIN"
echo "$(date) Analysis steps to execute: $LIST_STEPS_EXEC"
echo "$(date) Sample to analyise: $SAMPLE"
echo "$(date) BAM to analyise: $BAM_FILE"
echo "$(date) Regions to scan: $REGIONS"
echo "$(date) BED to use: $BED_FILE"


###################################################################################
# Start the analysis
###################################################################################
if ([[ -f $BAM_FILE ]] && [[ -f $BED_FILE ]]);then
	###################################################################################
	# Set output file
	###################################################################################
	OUT_FILE_D="$PATH_REPOR_TXT"/"$SAMPLE"_"$REGIONS"_DETAIL.txt
	OUT_FILE_S="$PATH_REPOR_TXT"/"$SAMPLE"_"$REGIONS"_stats.txt
	if [[ $LIST_STEPS_EXEC =~ " S1.01 " ]];then echo "$(date) File OUT Stats: $OUT_FILE_S";fi
	if [[ $LIST_STEPS_EXEC =~ " S1.02 " ]];then echo "$(date) File OUT Detail: $OUT_FILE_D";fi
	if [[ $LIST_STEPS_EXEC =~ " S1.01 " ]];then echo -e "REGION_ID\tREGION_CHR\tREGION_START\tREGION_END\tFRAGMENT_COUNT\tFRAGMENT_MEDIAN\tFRAGMENT_AVERAGE\tFRAGMENT_SD\tDIST" > $OUT_FILE_S;fi
	if [[ $LIST_STEPS_EXEC =~ " S1.02 " ]];then echo -e "REGION_ID\tREGION_CHR\tREGION_START\tREGION_END\tREAD_ID\tREAD_CHR\tREAD_START\tREAD_END\tFRAGMENT_LEN" > $OUT_FILE_D;fi

	###################################################################################
	# Check chr naming convention
	###################################################################################
	CHR1st=$(samtools view $BAM_FILE | head -n 1 | awk -F "\t" '{print $3}')

	###################################################################################
	# Loop through all the BED regions
	###################################################################################
	while read r;do
		region=$r
		IFS=$'\t' read -r -a region <<< "$region"
		CHR=${region[0]}
		if [[ $CHR1st =~ "chr" ]];then
			if [[ ! $CHR =~ "chr" ]];then
				CHR="chr$CHR"
			fi
		else
			CHR=${CHR/chr/}
		fi
		R_START=$((${region[1]}+1))
		R_END=$((${region[2]}+1))
		F_START=$(($R_START-$MAX_FRAGMENT_LEN))
		F_END=$(($R_END+$MAX_FRAGMENT_LEN))
		R_ID=${region[3]}
		echo "Processing region: $R_ID"
		if [[ $LIST_STEPS_EXEC =~ " S1.01 " ]];then
			{ samtools view -f 99 $BAM_FILE $CHR:$F_START-$F_END | awk \
			-v MIN_MAPQ=$MIN_MAPQ -v MAX_FRAGMENT_LEN=$MAX_FRAGMENT_LEN \
			-v CHR=$CHR -v R_START=$R_START -v R_END=$R_END -v R_ID=$R_ID \
			-f $AWK_FILE_FILTE ; \
			samtools view -f 163 $BAM_FILE $CHR:$F_START-$F_END | awk \
			-v MIN_MAPQ=$MIN_MAPQ -v MAX_FRAGMENT_LEN=$MAX_FRAGMENT_LEN \
			-v CHR=$CHR -v R_START=$R_START -v R_END=$R_END -v R_ID=$R_ID \
			-f $AWK_FILE_FILTE ; } | sort -k9 -n | awk \
			-v CHR=$CHR -v R_START=$R_START -v R_END=$R_END -v R_ID=$R_ID \
			-f $AWK_FILE_STATS >> $OUT_FILE_S
		fi
		if [[ $LIST_STEPS_EXEC =~ " S1.02 " ]];then
			{ samtools view -f 99 $BAM_FILE $CHR:$F_START-$F_END | awk \
			-v MIN_MAPQ=$MIN_MAPQ -v MAX_FRAGMENT_LEN=$MAX_FRAGMENT_LEN \
			-v CHR=$CHR -v R_START=$R_START -v R_END=$R_END -v R_ID=$R_ID \
			-f $AWK_FILE_FILTE ; \
			samtools view -f 163 $BAM_FILE $CHR:$F_START-$F_END | awk \
			-v MIN_MAPQ=$MIN_MAPQ -v MAX_FRAGMENT_LEN=$MAX_FRAGMENT_LEN \
			-v CHR=$CHR -v R_START=$R_START -v R_END=$R_END -v R_ID=$R_ID \
			-f $AWK_FILE_FILTE ; } | awk -f $AWK_FILE_DETAI >> $OUT_FILE_D
		fi
	done < $BED_FILE
	if [[ $LIST_STEPS_EXEC =~ " S1.02 " ]];then
		tar -czvf $OUT_FILE_D.tar.gz $OUT_FILE_D
	fi
fi
echo "$(date) END analysis"
