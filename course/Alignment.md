Sequence Alignment Tutorial
========================================================
author: "MRC CSC Bioinformatics Core Team"
date:http://mrccsc.github.io/training.html
width: 1440
height: 1100
autosize: true
css:style.css


MRC CSC Training Resources
========================================================
* [Reproducible R](http://mrccsc.github.io/Reproducible-R/)
* [Intermediate R - Data analysis sand Visualisation](http://bioinformatics-core-shared-training.github.io/r-intermediate/)
* [Statistics in R](http://mrccsc.github.io/StatisticsInR/)
* [Genomic File Formats](http://mrccsc.github.io/genomic_formats/)
* [ChIP-seq (short)](http://mrccsc.github.io/ChIPseq_short/)
* [RNA-seq (short)](http://mrccsc.github.io/RNAseq_short/)

Overview
========================================================

* [Introduction to Seqeuncing Technology](#/Intro)
* [Revisiting File Formats](#/FileFormat)
* [Quality Assessment of Sequencing](#/SeqQC)
* [Introduction to Alignment](#/align)
* [Aligning Sequences](#/align2)
* [Sorting and indexing](#/Sorting)
* [Summary & Post-Alignment QC](#/QC)
* [Reading BAM/SAM](#/readbam)
* [Visualisation](#/vis)


Downloading data
========================================================

https://fileexchange.imperial.ac.uk/pickup.php?claimID=T3eei7jDEorDnrCf&claimPasscode=29FSJPJ5PCYPJB8H&emailAddr=gdharmal%40ic.ac.uk

- Download the zip and unzip under "course" directory (this will create a new directory "Data")



Introduction to Seqeuncing Technology
========================================================
id: Intro

Illumina - Sequencing by synthesis

![Illumina](./images/Illumina.png)

#### Image Source: Illumina Inc

Introduction to Seqeuncing Technology
========================================================

![Illumina](./images/Illumina2.png)

Image Source: Illumina Inc


Introduction to Seqeuncing Technology
========================================================

![Multiplexing](./images/Multiplex_Seqeuncing.jpg)


Revisiting File Formats
========================================================
type: section
id: FileFormat


Revisiting File Formats
========================================================

## Genomic File Format Course
<br><br><br>
## http://mrccsc.github.io/genomic_formats/genomicFileFormats.html


Quality Assessment of Sequencing
========================================================
type: section
id: SeqQC


Quality Assessment of Sequencing
========================================================

<br>
### Tools for Alignment & QC
![Multiplexing](./images/Alignment_flow.jpg)

Quality Assessment of Sequencing
========================================================
- Summary
- Quality score distribtion (5' to 3')
- Adapter contamination
- GC Content
- Distribution of A/T/G/C along the sequence


Quality Assessment of Sequencing
========================================================

<br>
### R Packages: ShortRead, Rqc, SeqTools

<br>
### Quality assessment using ShortRead

<br> Key Functions in ShortRead
- <b>FastqStreamer</b> Iterate through FASTQ files in chunks
- <b>readFastq</b> Read an entire FASTQ file into memory
- <b>alphabetFrequency</b> Nucleotide or quality score use per read
- <b>alphabetByCycle</b> Nucleotide or quality score use by cycle
- <b>report</b> Generate a quality assessment report



```r
library(ShortRead)

fastqfile <- system.file(package="ShortRead", "extdata", "E-MTAB-1147","ERR127302_1_subset.fastq.gz")
fq <- readFastq(fastqfile)
```

Quality Assessment of Sequencing
========================================================

<br>
### Exploring sequences and Quality Scores

```r
fq
```

```
class: ShortReadQ
length: 20000 reads; width: 72 cycles
```

```r
head(sread(fq), 2)
```

```
  A DNAStringSet instance of length 2
    width seq
[1]    72 GTCTGCTGTATCTGTGTCGGCTGTCTCGCGG...CAATGAAGGCCTGGAATGTCACTACCCCCAG
[2]    72 CTAGGGCAATCTTTGCAGCAATGAATGCCAA...GTGGCTTTTGAGGCCAGAGCAGACCTTCGGG
```

```r
head(quality(fq), 2)
```

```
class: FastqQuality
quality:
  A BStringSet instance of length 2
    width seq
[1]    72 HHHHHHHHHHHHHHHHHHHHEBDBB?B:BBG...FEFBDBD@DDECEE3>:?;@@@>?=BAB?##
[2]    72 IIIIHIIIGIIIIIIIHIIIIEGBGHIIIIH...IHIIIHIIIIIGIIIEGIIGBGE@DDGGGIG
```

Quality Assessment of Sequencing
========================================================

<br>
### Quality assessment Report

```r
qaSummary <- qa(fastqfile, type="fastq")
browseURL(report(qaSummary))
```



Quality Assessment of Sequencing
========================================================

<br><br>
### Non-R approach: FASTQC

- Handles FASTQ, BAM and SAM Files
- Provies summary table and graphs
- Exports results in HTML
- GUI & ability to run multiple FASTQ in offline mode

### Webpage: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/


Quality Assessment of Sequencing
========================================================
![FASTQC1](./images/FAQSTQC1.png)


Quality Assessment of Sequencing
========================================================
![FASTQC2](./images/FASTQC2.png)


Adapter & Quality Trimming
========================================================

- Seqtk (https://github.com/lh3/seqtk)
- Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
- Cutadpat (https://github.com/marcelm/cutadapt)
- TrimGalore! - Wrapper for Cutadapt (http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)


Time for Exercises!
========================================================
* [Exercises Part1](./AlignmentExercise1.html)
<br><br>

* [Exercises Part1 Solutions](./AlignmentExercise1_Solutions.html)




Sequence Alignment
========================================================
type: section
id: align


Sequence Alignment
========================================================

![Alignment](./images/alignment1.jpg)

### Challenges
- Millions of short reads (36bp to 100bp)
- Identifying unique positions in the genome
- Repetitive regions in the genome
- Orientation of reads relative to genome is not known
- Sequencing reads have error

Sequence Alignment
========================================================

<br>
### How to
- Align reads efficiently in terms of memory and time?
- Account for mismatches and repeat regions 


<br><br>
### Two approaches
- Spaced seed and Extend
    + MAQ
    + BFAST
    + SHRiMP
- Burrows-Wheeler transform (BWT)
    + BWA
    + Bowtie


Spaced Seed & Extend
========================================================

![seed](./images/spacedseed.jpeg)

- Reference genome cut in to small "seeds"
- Pairs of spaced seeds are stored in index (lookup table)
- Seqeucning reads cut into 4 equal parts (read seeds)
- Seed pairs are used as keys to search the lookup table

#### Image: Nature Biotechnology 2009 27(5): 455-457

Burrows-Wheeler transform (BWT)
========================================================

![bwt](./images/bwt.jpeg)

- BWT helps to index entire human genome in less than 2 gb memory
- Aligner matches end of reads against the index (increasing one base at a time)
- If no perfect alignment found,  back up and try again with substitution
- Breakdown approach: First solving simple sub-problem (aligning one base) and use that to solve slighly harder problem (aligning two base)
- ~30 fold faster than spaced seed methods

#### Image: Nature Biotechnology 2009 27(5): 455-457

Spliced read aligners
========================================================

![alignment2](./images/alignment2.jpg)


- Aligning transcript reads to genome (reads overlapping with splice junction)
- With or without known splice junction databases
- Aligners: Tophat, STAR, MapSplice


Aligning Sequences
========================================================
type: section
id: align2


Aligning sequences with Rbowtie
========================================================

<br>
### Required Input
- Genome fasta
- Sequencing reads (fastq)
- User specified parameters 

<br><br><br>
### Data

- ChIPseq of CTCF, bone marrow macrophage (Mus musculus) from ENCODE (https://www.encodeproject.org/experiments/ENCSR000CFJ/)
- Reads from chr1 used for this exercise

Aligning sequences with Rbowtie
========================================================

<br>
### Building Index for the genome

#### Help for building index: bowtie_build_usage()


```r
library(Rbowtie)
bowtie_build_usage()
 [1] "Usage: bowtie2-build-s [options]* <reference_in> <ebwt_outfile_base>"            
 [2] "    reference_in            comma-separated list of files with ref sequences"    
 [3] "    ebwt_outfile_base       write Ebwt data to files with this dir/basename"     
 [4] "Options:"                                                                        
 [5] "    -f                      reference files are Fasta (default)"                 
 [6] "    -c                      reference sequences given on cmd line (as <seq_in>)" 
 [7] "    -C/--color              build a colorspace index"                            
 [8] "    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting"    
 [9] "    -p/--packed             use packed strings internally; slower, uses less mem"
[10] "    --bmax <int>            max bucket sz for blockwise suffix-array builder"    
[11] "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)"    
[12] "    --dcv <int>             diff-cover period for blockwise (default: 1024)"     
[13] "    --nodc                  disable diff-cover (algorithm becomes quadratic)"    
[14] "    -r/--noref              don't build .3/.4.ebwt (packed reference) portion"   
[15] "    -3/--justref            just build .3/.4.ebwt (packed reference) portion"    
[16] "    -o/--offrate <int>      SA is sampled every 2^offRate BWT chars (default: 5)"
[17] "    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)" 
[18] "    --ntoa                  convert Ns in reference to As"                       
[19] "    --seed <int>            seed for random number generator"                    
[20] "    -q/--quiet              verbose output (for debugging)"                      
[21] "    -h/--help               print detailed description of tool and its options"  
[22] "    --usage                 print this usage message"                            
[23] "    --version               print version information and quit"                  
```

Aligning sequences with Rbowtie
========================================================

<br>
### Building Index


```r
refFasta <- list.files(path="Data",pattern="chr1.fa$", full=T)
tmp <- bowtie_build(references=refFasta, outdir="Data", prefix="mm9index_chr1", force=TRUE)
```

Aligning sequences with Rbowtie
========================================================

<br>
#### Help for aligning sequences: bowtie_usage()


```r
bowtie_usage()
 [1] "Usage: "                                                                          
 [2] "bowtie-build-s [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]"      
 [3] ""                                                                                 
 [4] "  <m1>    Comma-separated list of files containing upstream mates (or the"        
 [5] "          sequences themselves, if -c is set) paired with mates in <m2>"          
 [6] "  <m2>    Comma-separated list of files containing downstream mates (or the"      
 [7] "          sequences themselves if -c is set) paired with mates in <m1>"           
 [8] "  <r>     Comma-separated list of files containing Crossbow-style reads.  Can be" 
 [9] "          a mixture of paired and unpaired.  Specify \"-\" for stdin."            
[10] "  <s>     Comma-separated list of files containing unpaired reads, or the"        
[11] "          sequences themselves, if -c is set.  Specify \"-\" for stdin."          
[12] "  <hit>   File to write hits to (default: stdout)"                                
[13] "Input:"                                                                           
[14] "  -q                 query input files are FASTQ .fq/.fastq (default)"            
[15] "  -f                 query input files are (multi-)FASTA .fa/.mfa"                
[16] "  -r                 query input files are raw one-sequence-per-line"             
[17] "  -c                 query sequences given on cmd line (as <mates>, <singles>)"   
[18] "  -C                 reads and index are in colorspace"                           
[19] "  -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C"  
[20] "  --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively"         
[21] "  -s/--skip <int>    skip the first <int> reads/pairs in the input"               
[22] "  -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)"    
[23] "  -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads"                
[24] "  -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads"               
[25] "  --phred33-quals    input quals are Phred+33 (default)"                          
[26] "  --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)"        
[27] "  --solexa-quals     input quals are from GA Pipeline ver. < 1.3"                 
[28] "  --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3"                
[29] "  --integer-quals    qualities are given as space-separated integers (not ASCII)" 
[30] "Alignment:"                                                                       
[31] "  -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities"  
[32] "    or"                                                                           
[33] "  -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)"          
[34] "  -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)" 
[35] "  -l/--seedlen <int> seed length for -n (default: 28)"                            
[36] "  --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)" 
[37] "  -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)"   
[38] "  -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)" 
[39] "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)"    
[40] "  --nofw/--norc      do not align to forward/reverse-complement reference strand" 
[41] "  --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)"  
[42] "  --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)"   
[43] "  -y/--tryhard       try hard to find valid alignments, at the expense of speed"  
[44] "  --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)" 
[45] "Reporting:"                                                                       
[46] "  -k <int>           report up to <int> good alignments per read (default: 1)"    
[47] "  -a/--all           report all alignments per read (much slower than low -k)"    
[48] "  -m <int>           suppress all alignments if > <int> exist (def: no limit)"    
[49] "  -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best" 
[50] "  --best             hits guaranteed best stratum; ties broken by quality"        
[51] "  --strata           hits in sub-optimal strata aren't reported (requires --best)"
[52] "Output:"                                                                          
[53] "  -t/--time          print wall-clock time taken by search phases"                
[54] "  -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)"   
[55] "  --quiet            print nothing but the alignments"                            
[56] "  --refout           write alignments to files refXXXXX.map, 1 map per reference" 
[57] "  --refidx           refer to ref. seqs by 0-based index rather than name"        
[58] "  --al <fname>       write aligned reads/pairs to file(s) <fname>"                
[59] "  --un <fname>       write unaligned reads/pairs to file(s) <fname>"              
[60] "  --max <fname>      write reads/pairs over -m limit to file(s) <fname>"          
[61] "  --suppress <cols>  suppresses given columns (comma-delim'ed) in default output" 
[62] "  --fullref          write entire ref name (default: only up to 1st space)"       
[63] "Colorspace:"                                                                      
[64] "  --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)"    
[65] "     or"                                                                          
[66] "  --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred" 
[67] "  --col-cseq         print aligned colorspace seqs as colors, not decoded bases"  
[68] "  --col-cqual        print original colorspace quals, not decoded quals"          
[69] "  --col-keepends     keep nucleotides at extreme ends of decoded alignment"       
[70] "SAM:"                                                                             
[71] "  -S/--sam           write hits in SAM format"                                    
[72] "  --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments"  
[73] "  --sam-nohead       supppress header lines (starting with @) for SAM output"     
[74] "  --sam-nosq         supppress @SQ header lines for SAM output"                   
[75] "  --sam-RG <text>    add <text> (usually \"lab=value\") to @RG line of SAM header"
[76] "Performance:"                                                                     
[77] "  -o/--offrate <int> override offrate of index; must be >= index's offrate"       
[78] "  -p/--threads <int> number of alignment threads to launch (default: 1)"          
[79] "  --mm               use memory-mapped I/O for index; many 'bowtie's can share"   
[80] "  --shmem            use shared mem for index; many 'bowtie's can share"          
[81] "Other:"                                                                           
[82] "  --seed <int>       seed for random number generator"                            
[83] "  --verbose          verbose output (for debugging)"                              
[84] "  --version          print version information and quit"                          
[85] "  -h/--help          print this usage message"                                    
```


Aligning sequences with Rbowtie
========================================================


```r
bowtie(sequences="Data/ENCFF001LGM_chr1.fastq",
       index=file.path("Data", "mm9index_chr1"), 
       outfile="Data/CTCF_mm9_MF.sam",
       sam=TRUE, 
       best=TRUE, 
       force=TRUE,
       threads=2,
       type="single", m=1
       )
```


```r
# reads processed: 975396
# reads with at least one reported alignment: 967484 (99.19%)
# reads that failed to align: 7340 (0.75%)
# reads with alignments suppressed due to -m: 572 (0.06%)
# Reported 967484 alignments to 1 output stream(s)
```



Spliced Alignment using SpliceMap
========================================================

- To generate spliced alignments (RNA-Seq)
- SpliceMap uses bowtie
- Accepts paramets as named list
- Limitation: reads should be at least 50bp


### Data
- RNA-Seq of Liver (Mus musculus E12.5) from ENCODE: https://www.encodeproject.org/experiments/ENCSR648YEP/
- Single end
- Reads from chr1 used for this exercise


Spliced Alignment using SpliceMap
========================================================

<br>
### Aligning the sequences using SpliceMap


```r
readsFiles <- "Data/ENCFF905YGD_chr1.fastq"
refDir <- "Data/chr1.fa"
samFiles <- "Data/ENCFF905YGD_chr1.sam"
cfg <- list(genome_dir=refDir,
            reads_list1=readsFiles,
            read_format="FASTQ",
            quality_format="phred-33",
            outfile=samFiles,
            temp_path=tempdir(),
            max_intron=400000,
            max_multi_hit=10,
            seed_mismatch=1,
            read_mismatch=2,
            num_chromosome_together=2,
            bowtie_base_dir="Data/mm9index_chr1",
            num_threads=4,
            try_hard="yes",
            selectSingleHit=TRUE)
res <- SpliceMap(cfg)
```



Sorting & Indexing
========================================================
type: section
id: Sorting


Sorting & Indexing the SAM/BAM
========================================================

<br>
### Sorting SAM/BAM

- Most post-processing methods requires reads sorted by aligned position
- Sorted and indexed large files enable random access of required genomic regions

### Tools
- R: Rsamtools
- Non-R: samtools, Picard


Sorting & Indexing the SAM/BAM
========================================================
<br>
### Convert SAM to BAM

- `asBam` converts SAM to BAM, sort and index

```r
library(Rsamtools)
insam <- "Data/CTCF_mm9_MF.sam"
outbam <- "Data/CTCF_mm9_MF"
asBam(insam, outbam, overwrite=T)
```

- Sorting and Indexing a BAM

```r
sortBam(InputBAM, SortedBAM, byQname=F, maxMemory=512)
indexBam(SortedBAM)
```


Few key functions in Rsamtools
========================================================

- <b>scanBamHeader</b>: Read the SAM header section
- <b>countBam</b>: Count number reads (accepts `ScanBamParam` object)
- <b>mergeBam</b>: Merge BAM files
- <b>asSam</b>: Convert BAM to SAM

<br><br>

































```
Error in value[[3L]](cond) : 
  failed to open BamFile: file(s) do not exist:
  'Data/CTCF_mm9_MF.bam'
```
