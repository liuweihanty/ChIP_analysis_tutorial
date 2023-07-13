#!/bin/bash
#This scrip takes the forward and reverse read fastq files of ChIP Seq/CnR, runs FastQC, bwa alignment, and samtool filtering 

#PBS -N run_CnR
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -o /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/logs/run_CnR_seq.out
#PBS -e /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/logs/run_CnR_seq.err


#specify the input fastq files
echo $fq_F
echo $fq_R
echo $fq_location

p_value=0.1
q_value=0.1
run_macs2=false

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -macs2)
            run_macs2=true
            ;;
        -p)
            p_value="$2"
            shift
            ;;
        -q)
            q_value="$2"
            shift
            ;;
        *)
            echo "Unknown option: $key"
            exit 1
            ;;
    esac
    shift
done


#specify the location of reference genome fast file for bwa mem
GENOME_FASTA=/gpfs/data/mcnerney-lab/mcnerney/reference/hg19/hg19_Ordered.fa

#specify the path to the blasklist file
blacklist=/gpfs/data/mcnerney-lab/liuweihan/chip_seq/Human_CD34_HSC_Cut_Run_Jeff/hg19-blacklist.v2.bed

#grab base name of the fastq files
base=`basename $fq_F .fastq.gz`
echo "Sample name is $base"


#specify the number of cores to use for various downstream analysis
cores=8


#create all the output directories
# The -p option means mkdir will create the whole path if it does not exist and refrain from complaining if it does exist

mkdir -p /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa
mkdir -p /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bigwigs


# set up output filenames and locations

bwa_mem_out=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.sam

samtools_q30_in=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.sam
samtools_q30_out=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.bam

samtools_sort_in=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.bam
samtools_sort_out=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.srt.bam

samtools_dedup_in=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.srt.bam
samtools_dedup_out=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.srt.dedup.bam

bedtools_rm_blacklist_in=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.srt.dedup.bam
bedtools_rm_blacklist_out=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.srt.dedup.blkrm.bam

bamCoverage_in=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bwa/${base}.q30.srt.dedup.blkrm.bam
bamCoverage_out=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/bigwigs/${base}.q30.srt.dedup.blkrm.bw

#run the jobs

#genome alignment
echo "Run bwa mem"

cd $fq_location

module load gcc/6.2.0
module load intel/2017
module load bwa/0.7.17


bwa mem -t $cores \
$GENOME_FASTA \
$fq_F \
$fq_R \
> $bwa_mem_out

#filter and clean bam files
echo "Run samtools filter"

module load gcc/6.2.0
module load samtools/1.10

#q30 filtering
samtools view -bSh -q 30 -o $samtools_q30_out $samtools_q30_in

#sort
samtools sort -o $samtools_sort_out $samtools_sort_in

#dedup
samtools rmdup $samtools_dedup_in $samtools_dedup_out

#index and generate bai file
samtools index $samtools_dedup_out


echo "run bedtools to remove blacklist regions from bam files"
module load gcc/6.2.0
module load bedtools/2.29.0

bedtools intersect -abam $bedtools_rm_blacklist_in -b $blacklist -v > $bedtools_rm_blacklist_out
#remove the intermediate index file
rm *.dedup.bam.bai

echo "index the final bam file"
samtools index $bedtools_rm_blacklist_out

#generate bigwig files
#you need deepTools for this, and deepTools could be called on as long as you loaded the correct python version, so you don't need to load deepTools separately.
echo "generating bigwig files"

module load gcc/6.2.0
module load python/3.8.1

bamCoverage -b $bamCoverage_in -o $bamCoverage_out

#macs2 peak calling
if [[ $run_macs2=true ]]; then

    echo "Running macs2 analysis"

    # Check if -p or -q are specified
    if [[ -z $p_value ]] && [[ -z $q_value ]]; then
        echo "Error: Either -p or -q option must be specified with -macs2."
        exit 1
    fi
 
    mkdir -p /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/macs2
    macs2_outdir=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/output/macs2
    module load gcc/6.2.0
    module load intel/2017
    module load python/2.7.13
    module load macs2/2.1.0

    if [[ -n $p_value ]]; then
        echo "-p $p_value"
        macs2 callpeak -t $bedtools_rm_blacklist_out \
                       -f BAMPE \
                       -g hs \
                       -n $base \
                       -p $p_value \
                       -B \
                       --outdir $macs2_outdir

    elif [[ -n $q_value ]]; then
        echo "-q $q_value"
        macs2 callpeak -t $bedtools_rm_blacklist_out \
                       -f BAMPE \
                       -g hs \
                       -n $base \
                       -p $q_value \
                       -B \
                       --outdir $macs2_outdir
    else
        echo "No parameter (-p or -q) specified."
        exit 1
    fi

else
    echo "Skipping the macs2 peak calling analysis. You can run it separately"
fi


