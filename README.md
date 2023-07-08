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
  You can do this in either linux or R. [fastqcr](http://www.sthda.com/english/wiki/fastqcr-an-r-package-facilitating-quality-controls-of-sequencing-data-for-large-numbers-of-samples) package provided an easy implementation in R language. You can run this in your local desktop. The fastqcr code is attached in this folder. Please see inside for details.
  
* ### Remove illumina library adaptor.
   You will obtain sequencing adaptor sequence from the sequencing facility. If not, you can run the fastqc and in the last section "adaptpor sequence" you will see it. Typical Illumina library sequencing adaptor sequences could be seen [here](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314)
* ### set up your working directory
    
* ### Trim fastq files for sequencing adaptors
  * #### If you don't know the sequencing adaptor 
* ### modify the job submission file (wrapper file)
* ### change directory to the folder that contains your job file and job submission file, type in ``` chmod +x * ```
  
## Broad peak calling 



