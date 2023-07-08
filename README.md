# Table of contents <br>
 - [Introduction](#introduction)
 - [Analysis workflow](#analysis_workflow)
 - [Step by step analysis](#Step_by_step_analysis)
 - [Broad peakcalling](#broad_peakcalling)

## Introduction <br>
This tutorial walks step-by-step tutorial of analysis pipeline for ChIP-seq/CUT&RUN. In my experience, I found you can generally use the same analysis workflow for the two types of experiment, but there are studies proposing tailored CUT&RUN analysis tools such as [SEACR](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4), you are welcome to experimening orthogonal approaches and becnchmark their performance. 

## Analysis workflow
![GitHub Logo](https://github.com/liuweihanty/ChIP_analysis_tutorial/blob/main/figures/ChIP_CnR_workflow.png)

## Step by step analysis
* ### Run fastqc
  You can do this in either linux or R. [fastqcr](http://www.sthda.com/english/wiki/fastqcr-an-r-package-facilitating-quality-controls-of-sequencing-data-for-large-numbers-of-samples) package provided an easy implementation in R language. You can run this in your local desktop. The fastqcr code is attached in this folder named fastqc.Rmd. Please see the [document](https://github.com/liuweihanty/ChIP_analysis_tutorial/blob/f4982c5fd9c9e25d493fb50f1813dc429562869b/fastqc.Rmd) for details.

* ### Remove illumina library adaptor.
  1. #### You will obtain adaptor sequence from the sequencing facility. If not, from fastqc results, in the last section "adaptor sequence" you will see it. Typical Illumina library sequencing adaptor sequences could be seen [here](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314) <br>
  2. #### Use[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim the adaptors. See [here](https://cutadapt.readthedocs.io/en/stable/installation.html) for how to install Cutadapt <br>
  3. #### Run the code below, swap the AACCGGTT for your actual adaptor sequence <br>
   ```cutadapt -a AACCGGTT -o output.fastq input.fastq```

* ### (optional) Run fastqcr again to ensure the adaptors are successfully removed
  
* ### set up your working directory. (the demo example folder names are written in parenthesis)
  1. #### create your project folder. <SMARCA4_ChIP>
  2. #### create four sub-folders underneath your project folder
     * ##### </SMARCA4_ChIP/input> trimmed fastqs
     * ##### </SMARCA4_ChIP/output> the analysis output
         * ###### </SMARCA4_ChIP/output/bwa> alignment output
         * ###### </SMARCA4_ChIP/output/macs2> peak callint output
           
    
* ### Trim fastq files for sequencing adaptors
  * #### If you don't know the sequencing adaptor 
* ### modify the job submission file (wrapper file)
* ### change directory to the folder that contains your job file and job submission file, type in ``` chmod +x * ```
  
## Broad peak calling 



