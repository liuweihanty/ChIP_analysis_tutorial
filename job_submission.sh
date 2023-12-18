
#!/bin/bash
#specify your project directory,change for different analysis
project_dir=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR 

#change directory to where the fastq files are
cd $project_dir/input

#this for loop will take the input fastq files and run the scripts for all of them one pair after another

for i in $(ls *R1*.gz)
do
otherfilename="${i/R1/R2}"
echo $i
echo $otherfilename


#here you can specify whether to run MACS2 peak calling and the p value threshold, these two parameters will be passed along to the run_job.sh file
#whether to run macs2, if you include this flag, the problem will run macs2 peak caller, if not, the program will skip macs2.
run_macs2=true
#p value for macs2. Use p value as the significant thrshold and specify it to be 0.1 if you are running IDR.
p_value=0.1  # Adjust as needed

sbatch --job-name=run_ChIP_wrapper --time=12:00:00 \
           -o $project_dir/logs/run_Chip_seq.out \
           -e $project_dir/logs/run_Chip_seq.err \
           --partition=tier2q \
           --nodes=1 \
           --ntasks-per-node=8 \
           --mem-per-cpu=10000 \
           --wrap="sh $project_dir/scripts/run_job.sh $i $otherfilename $run_macs2 $p_value $project_dir"


      
 done
