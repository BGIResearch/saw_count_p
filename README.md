# Bam2Gem

## Introduction

Bam2Gem implement these steps in RNA-seq pipeline:

1. mapping quality filter
2. deduplication
3. set annotations
4. stat gene expression data

optional functions:

1. umi correction
2. sequencing saturation

## Compile

### Platform & Environment

* centos-7.0+
* gcc-9.1.0
* cmake-3.17.2

### Prerequisites

| Library           | Version | Description                        | Link                                                                           |
| ----------------- | ------- | ---------------------------------- | ------------------------------------------------------------------------------ |
| htslib            | 1.14.0  | process bam/sam data               | https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2  |
| spdlog            | 1.5.0   | logging module                     | https://github.com/gabime/spdlog/archive/v1.5.0.zip                            |
| CLI11             | 1.9.0   | parse command line parameters      | https://github.com/CLIUtils/CLI11/releases/download/v1.9.0/CLI11.hpp           |
| libdeflate        | 1.5     | accelerate bgzf IO                 | https://github.com/ebiggers/libdeflate/archive/v1.5.zip                        |
| ygg               | master  | interval tree for find overlapping | https://github.com/tinloaf/ygg.git                                             |
| doctest           | 2.3.7   | optional for unittest              | https://github.com/onqtam/doctest/archive/2.3.7.tar.gz                         |
| readerwriterqueue | 1.0.5   | rw-queue                           | https://github.com/tigeroses/readerwriterqueue/archive/refs/tags/v1.0.5.tar.gz |
| libx              | 1.1     | common library                     | https://github.com/tigeroses/libx/archive/refs/tags/v1.1.tar.gz                |
| hdf5-1.12.1       | 1.12.1  | hdf5 library                       | https://github.com/HDFGroup/hdf5                                               |

### Compile

1. enter the code path
2. modify the value of `libPath` in file *script/build.sh*
3. run *script/build.sh*
4. the compile result saved as: *install/bin/Bam2Gem*

### UnitTest

Unit testing is already supported.  
Repeat steps in *Compile* and modify step three: `sh ./script/build.sh ON` , 
then the binary executable file saved as: *install/bin/unittest*

## Run

### Usage

```text
$./install/bin/Bam2Gem -h
Bam2Gem: mapping quality filter, deduplication, set annotation, stat gene expression data.
Usage: ./install/bin/Bam2Gem [OPTIONS]

Options:
  -h,--help                             Print this help message and exit
  -I,-i TEXT REQUIRED                   Input bam filename or file list separated by comma
  -O,-o TEXT REQUIRED                   Output bam filename
  -A,-a TEXT:FILE REQUIRED              Input annotation filename
  -S,-s TEXT REQUIRED                   Output summary filename
  -E,-e TEXT REQUIRED                   Output barcode gene expression filename
  --sn TEXT REQUIRED                    STOmics Chip Serial Number
  -Q,-q INT:POSITIVE                    Set mapping quality threshold, default 10
  -C,-c INT:POSITIVE                    Set cpu cores, default detect
  -M,-m INT:POSITIVE                    Set avaliable memory(GB), default detect
  --save_lq                             Save low quality reads, default false
  --save_dup                            Save duplicate reads, default false
  --umi_on                              Enable umi correction, default disable
  --umi_min_num INT:POSITIVE            Minimum umi number for correction, default 5
  --umi_mismatch INT:POSITIVE           Maximum mismatch for umi correction, default 1
  --umi_len INT:POSITIVE                UMI length, default 10
  --sat_file TEXT                       Output sequencing saturation file, default None
  --multi_map                           Enable multi-mapping reads correction, default disable

Bam2Gem version: 2.3.0
```

Required parameters:

* -i filename. Input bam filename
* -o filename. Output bam filename
* -a filename. Input annotation filename
* -s filename. Output summary filename
* -e filename. Output barcode gene expression filename

Optional parameters:

* -q integer. Set mapping quality threshold, default 10
* --save_lq. Save low quality reads(less than paramter of '-q'), default not save
* --save_dup. Save duplicate reads, default not save.
* --umi_on. Open umi correction, default off
* --umi_min_num integer. Minimum umi number for correction, default 5
* --umi_mismatch integer. Maximum mismatch for umi correction, default 1
* --sat_file filename. Output sequencing saturation file, depend on --umi_on
* ...

### Example

#### Normal Mode

```text
$./install/bin/Bam2Gem \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/summary.txt \
 -e batch25/exp.gef
```

Input parameters: `-i batch25/batch25.bam -a batch25/batch25.gtf`  
Ouput parameters: `-o batch25/exp.bam -s batch25/summary.txt -e batch25/exp.gef`

#### Save Mode

```text
$./install/bin/Bam2Gem \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/summary.txt \
 -e batch25/exp.gef \
 --save_lq \
 --save_dup
```

Save reads with low quality or duplicate, also set bam flags with *BAM_FQCFAIL* or *BAM_FDUP*

#### Use Umi Correction

```text
$./install/bin/Bam2Gem \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/summary.txt \
 -e batch25/exp.gef \
 --umi_on \
 --umi_min_num 5 \
 --umi_mismatch 1
```

Use umi correction for deduplicate

#### Sequencing Saturation

```text
$./install/bin/Bam2Gem \
 -i batch25/batch25.bam \
 -o batch25/exp.bam \
 -a batch25/batch25.gtf \
 -s batch25/summary.txt \
 -e batch25/exp.gef \
 --umi_on \
 --sat_file batch25/exp.saturation.txt
```

### Results

* output bam. The bam file contains gene annotation information
* summary. Text file contains stat metrics of filter, deduplication and annotation,example:

  1. filter and deduplication metrics table

  | TOTAL_READS | PASS_FILTER | UNIQUE_READS | FAIL_FILTER_RATE | DUPLICATION_RATE |
  | ----------- | ----------- | ------------ | ---------------- | ---------------- |
  | 9999404     | 6450176     | 4036230      | 35.4944          | 37.4245          |

  2. annotation metrics table

  | TOTAL_READS | MAP    | EXONIC | INTRONIC | INTERGENIC | TRANSCRIPTOME | ANTISENSE |
  | ----------- | ------ | ------ | -------- | ---------- | ------------- | --------- |
  | 949906      | 949906 | 855200 | 354      | 42033      | 850051        | 9833      |
  | 100.0       | 100.0  | 90.0   | 0.0      | 9.9        | 89.5          | 1.0       |

* sequencing saturation(exists if using parameter '--sat_file'), each line separated by space, no header line

  | coorY  | coorX | geneIndex | MIDIndex | readCount |
  | ------ | ----- | --------- | -------- | --------- |
  | 108565 | 84833 | 10539     | 1018895  | 1         |
  | 114929 | 88696 | 10539     | 929853   | 4         |
  
  coorY,coorX: coordinates of dnb, both together represent barcode; type is uint32
  geneIndex: index of gene; type is uint32
  MIDIndex: MID sequence; type is uint64
  readCount: read number of specific MID; type is uint32

* gene expression matrix. It use hdf5 format for storing gene expression data.

  
### Error Code

| code       | desc                        |
| ---------- | --------------------------- |
| SAW-A20001 | user-input parameters error |
| SAW-A20002 | disk IO error               |
| SAW-A20003 | program error               |
