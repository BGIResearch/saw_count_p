# Change Log

## 2.3.2(2022-09-13)

* add header for saturation file
* update for checking input/output file format

## 2.3.1(2022-07-28)

* remove reading index of input bam file
* generate csi index of output bam file

## 2.3.0(2022-05-24)

* add parameter --multi_map for correcting multi-mapping reads
* add project checkGTF for checking gtf/gff file format
* stat MID count of exon when dumping gene expression file
* exit if no valid gene found in gtf/gff files

## 2.2.2(2022-04-11)

* add parameter --sn for parsing resolution

## 2.2.1(2022-02-24)

* reduce log info when blocking

## 2.2.0(2022-01-20)

* replace gene expression file with h5 file
* add checking path is exists of output files
* fix int overflow when get system memory info

## 2.1.1(2022-01-11)

* fix bug: pseudo blocking

## 2.1.0(2021-12-24)

* set buffer size when read and write bam files
* split data for umi correction when memory is not enough
* add parameter "-m" for setting memory
* catch IO exception and abort
* set intron threshold as 50%
* update stat of antisense reads
* keep the gene and transcript when they have same name
* parse all feature type include "mRNA" in gff files 

## 2.0.2(2021-12-06)

* fix bug: block before join read thread

## 2.0.1(2021-11-26)

* update doc of workflow, umi correction, annotation
* refactor annotation module, umi correction module
* add unittest

## 2.0.0(2021-11-12)

* Speed up the program by treating single gene as process unit
* Remove scRNA mode
* Replace final output of sequencing saturation with intermediate output

## 1.2.1(2021-10-20)

* Support 32bp umi at most

## 1.2.0(2021-09-03)

* Support new input BAM formats that has "Cx/Cy/UR" tags in records
* Add comment for "XF" tag in output BAM header
* Change tag from "XF:Z" to "XF:i"

## 1.1.0(2021-08-11)

* add uniq reads in saturation file

## 1.0.8(2021-07-21)

* change value of x-axis in sequence saturation from mean reads to total reads

## 1.0.7(2021-07-16)

* trim '"' in gtf files
* fix random problem when correcting umi

## 1.0.6(2021-07-13)

* support empty line and space in gff files 

## 1.0.5(2021-06-28)

* support empty line and space in gtf files

## 1.0.4(2021-06-24)

* change the format of expression file, split barcode into x/y coordinates, and add header line

## 1.0.3(2021-05-31)

* fix integer overflow when stat reads

## 1.0.2(2021-03-19)

* change intenel directory name from fixed to random

## 1.0.1(2021-02-04)

* set transcript_name as transcript_id when parsing gtf file if transcript_name is empty

## 1.0.0(2021-01-06)

* the first release version

