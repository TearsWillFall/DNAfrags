# DNAfrags
DNAfrags is a tool to filter reads in BAM files by fragment length and genome position in R code language
## 1. Description
DNAfrags simplifies the process of filtering sequence reads in BAM files by fragment length and genome position, by internally making use of Samtools command-line tricks to access this information.

**Tools:**
* [samtools](https://github.com/samtools/samtools/):Tools (written in C using htslib) for manipulating next-generation sequencing data

## 2. System Requirements
Tested on fresh Ubuntu and Arch Linux installations. Should be working on other Linux distros too as long as equivalent packages are provided.

In order to be able to download and compile the source files of all the required tools the following programs are required:
* git
* make
* gcc
* autoconf

These tools can and should be installed using the terminal with the following commands:

* **For Ubuntu:**

  ```
  sudo apt install make gcc autoconf
  ```

* **For Arch Linux:**

  ```
  sudo pacman -S make gcc autoconf
  ```
  
Additional dependencies may be needed to succesfully install `devtools` package in R:
  
* **For Ubuntu:**

  ```
  sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev openssl
  ```
  
* **For Arch Linux:**

  ```
  sudo pacman -S build-essential libcurl-gnutls libxml2-dev openssl
  ```
  
## 3. Installation Instructions
  In order to install `DNAfrags` package we will be using R `devtools`:
  ```
  install.packages("devtools")
  devtools::install_github("TearsWillFall/DNAfrags")
  ```
If `devtools` package installation fails check System Requirements section, as you may be missing a dependency.

Once the `DNAfrags` package is installed we can use the function `install_required_tools()` to download and set up all the tools required for the bioinformatic     process. This will create a directory named `tools` in the current working directory with all the tools. **Note: All functions within this package call scripts from the `tools` directory, therefore if this directory is moved, deleted or the current working directory is changed, these functions will fail.**
  

```
DNAfrags::install_required_tools()
```

Or, alternatively.

```
library("DNAfrags")
install_required_tools() 
```
## 4. Usage

**Examples:**

In this example reads in the BAM file called "myBAM.bam" are only filtered by fragment length. Only fragments between lengths 10 and 150 (inclusive) will be kept  for all chromosomes.
```
DNAfrags::filter_fragments(file="myBAM.bam",min_frag_size=10,max_frag_size=150)
```

In this example reads in the BAM file called "myBAM.bam" are filtered by fragment length and a single genomic regions. Only fragments mapped across chromosome 6 between lengths 10 and 150 (inclusive) will be kept.
```
DNAfrags::filter_fragments(file="myBAM.bam",min_frag_size=10,max_frag_size=150,chr=6)
```

In this example reads in the BAM file called "myBAM.bam" are filtered by fragment length and a single genomic regions. Only fragments mapped after the initial 5566 bases in chromosome 6 between lengths 10 and 150 (inclusive) will be kept.
```
DNAfrags::filter_fragments(file="myBAM.bam",min_frag_size=10,max_frag_size=150,chr=6,start_pos=5566)
```

In this example reads in the BAM file called "myBAM.bam" are filtered by fragment length and a single genomic regions. Only fragments mapped between the bases 5566 and 1000000 (inclusive) in chromosome 6 between lengths 10 and 150 (inclusive) will be kept.
```
DNAfrags::filter_fragments(file="myBAM.bam",min_frag_size=10,max_frag_size=150,chr=6,start_pos=5566,end_pos=1000000)
```
**Note:To deal with multiple genomic regions 2 solutions exist**


* **BED solution [Faster]**

In this example reads in the BAM file called "myBAM.bam" are filtered by fragment length in multiple genomic regions. Fragments between lengths 10 and 150 (inclusive) will be kept for specified regions.
```
DNAfrags::filter_fragments(file="myBAM.bam",min_frag_size=10,max_frag_size=150,bed="myBED.bed")
```

* **Merge BAM solution**

In this example reads in the BAM file called "myBAM.bam" are filtered by fragment length in multiple genomic regions. Fragments between lengths 10 and 150 (inclusive) will be kept for specified regions.
```
DNAfrags::filter_fragments(file="myBAM.bam",min_frag_size=10,max_frag_size=150,chr=6,start_pos=5566)
DNAfrags::filter_fragments(file="myBAM.bam",min_frag_size=10,max_frag_size=150,chr=12,start_pos=123123,end_pos=1232131)
```
Then the BAM files for each genomic region are merged together into a single BAM file.
```
DNAfrags::merge_bam(out_bam="MyBAM",bam_dir=".")
```

