# Table of contents <br>
 - [Introduction](#introduction)
 - [Demo data](#demo_data)
 - [Analysis workflow](#analysis_workflow)
 - [Step by step analysis](#Step_by_step_analysis)
 - [Broad peakcalling](#broad_peakcalling)

## Introduction <br>
This tutorial walks step-by-step tutorial of analysis pipeline for ChIP-seq/CUT&RUN. In my experience, I found you can generally use the same analysis workflow for the two types of experiment, but there are studies proposing tailored CUT&RUN analysis tools such as [SEACR](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4), you are welcome to experimening orthogonal approaches and becnchmark their performance. 

## Demo data
I will be using an example data set to illustrate this workflow. This is a CUT&RUN experiment on human CD34+ HSPC, probing for SMARCA4 binding with and without knocking out CUX1 using CRISPR/Cas9. The metadata is whown below:

## Analysis workflow
![GitHub Logo](https://github.com/liuweihanty/ChIP_analysis_tutorial/blob/main/figures/ChIP_CnR_workflow_chart.png)

## Step by step analysis
* ### Run fastqc 
  You can do this in either linux or R. [fastqcr](http://www.sthda.com/english/wiki/fastqcr-an-r-package-facilitating-quality-controls-of-sequencing-data-for-large-numbers-of-samples) package provided an easy implementation in R language. You can run this in your local desktop. The fastqcr code is attached in this folder named fastqc.Rmd. Please see the [document](https://github.com/liuweihanty/ChIP_analysis_tutorial/blob/f4982c5fd9c9e25d493fb50f1813dc429562869b/fastqc.Rmd) for details.

* ### Remove illumina library adaptor.
  * #### You will obtain adaptor sequence from the sequencing facility. If not, from fastqc results, in the last section "adaptor sequence" you will see it. Typical Illumina library sequencing adaptor sequences could be seen [here](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314) <br>
  * #### Use [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim the adaptors. See [here](https://cutadapt.readthedocs.io/en/stable/installation.html) for how to install Cutadapt <br>
  * #### Run the code below, swap the AACCGGTT for your actual adaptor sequence <br>
   ```cutadapt -a AACCGGTT -o output.fastq input.fastq```

* ### (Optional) Run fastqcr again to ensure the adaptors are successfully removed
  
* ### Set up your working directory. (the demo example folder names are written in parenthesis)
  * #### create your project folder. </SMARCA4_ChIP/>
  * #### create four sub-folders underneath your project folder
     * ##### </SMARCA4_ChIP/input> trimmed fastqs
     * ##### </SMARCA4_ChIP/output> the analysis output
         * ###### </SMARCA4_ChIP/output/bwa> alignment output
         * ###### </SMARCA4_ChIP/output/macs2> peak callint output
     * ##### </SMARCA4_ChIP/logs> the error and output records files for debugging
     * ##### </SMARCA4_ChIP/scripts> the analysis scripts
           
* ### Run the job
    #### Now that we have the adaptor trimmed fastqs, it's time to proceeed to next steps. In the flow chart above, we finished steps 1 and 2 so far. Step 3 to 6 will be implemented in an automated workflow, which is organized into two bash scripts: <br>
    * #### job_submission.sh: this script specify and select the two fastq files(forward and reverse reads) for each sample, and send the fastqs to the "run.sh" script below
    * #### run.sh:  this scripts takes in the forward and reverse fastqs for each sample from the job_submission.sh file above, and performs steps 3-6 in the flow chart on each sample, in a paralelled fashion( all samples will be simutaneously analyzed), so no matter how many samples you have, you just need to run once.
    #### Now let's take a look inside of each file and I will explain what each code chunk do:
    #### job_submission.sh
   
   
   
   * ### change directory to the folder that contains your job file and job submission file, type in ``` chmod +x * ```
  
## Broad peak calling 



