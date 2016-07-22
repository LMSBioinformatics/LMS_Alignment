Sequence Alignment Tutorial
========================================================
css: Rpress.css
author: "MRC CSC Bioinformatics Core Team"
date:http://mrccsc.github.io/training.html
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
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
* [Reading BAM](#/readbam)
* [Visualisation](#/SeqQC)



Introduction to Seqeuncing Technology
========================================================
id: Intro

Illumina - Sequencing by synthesis

[Figure]


Introduction to Seqeuncing Technology
========================================================

![Multiplexing](./images/Multiplex_Seqeuncing.jpg)


Revisiting File Formats
========================================================
id: FileFormat

## Genomic File Format Course
<br><br><br>
## http://mrccsc.github.io/genomic_formats/genomicFileFormats.html)


Quality Assessment of Sequencing
========================================================
id: SeqQC

- Summary
- Quality score distribtion (5' to 3')
- Adapter contamination
- GC Content
- Distribution of A/T/G/C along the sequence


Quality Assessment of Sequencing
========================================================

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
id: align


![Alignment](./images/alignment1.jpg)

### Challenges
- Millions of short reads (36bp to 100bp)
- Identifying unique positions in the genome
- Repetitive regions in the genome
- Orientation of reads relative to genome is not known
- Sequencing reads have error

Sequence Alignment
========================================================

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

#{.right}
![seed](./images/spacedseed.jpeg)

- Reference genome cut in to small "seeds"
- Pairs of spaced seeds are stored in index (lookup table)
- Seqeucning reads cut into 4 equal parts (read seeds)
- Seed pairs are used as keys to search the lookup table


Burrows-Wheeler transform (BWT)
========================================================

#{.right}
![bwt](./images/bwt.jpeg)

- BWT helps to index entire human genome in less than 2 gb memory
- Aligner matches end of reads against the index (increasing one base at a time)
- If no perfect alignment found,  back up and try again with substitution
- Breakdown approach: First solving simple sub-problem (aligning one base) and use that to solve slighly harder problem (aligning two base)
- ~30 fold faster than spaced seed methods


Spliced read aligners
========================================================

![alignment2](./images/alignment2.jpg)


- Aligning transcript reads to genome (reads overlapping with exon-exon junction)
- With or without known exon-exon junction databases
- Aligners: Tophat, STAR, MapSplice


Aligning Sequences
========================================================
type: section
id: align2


Aligning sequences with Rbowtie
========================================================

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

### Building Index for the genome

#### Help for building index: bowtie_build_usage()


```r
library(Rbowtie)
bowtie_build_usage()
```

```
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

### Building Index


```r
refFasta <- list.files(path="Data",pattern="*.fa$", full=T)
tmp <- bowtie_build(references=refFasta, outdir="Data", prefix="mm9index_chr1_chr2", force=TRUE)
```

Aligning sequences with Rbowtie
========================================================

#### Help for aligning sequences: bowtie_usage()


```r
bowtie_usage()
```

```
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
       index=file.path("Data", "mm9index_chr1_chr2"), 
       outfile="Data/CTCF_mm9_MF.sam", 
       sam=TRUE, 
       best=TRUE, 
       force=TRUE,
       type="single", m=1
       )
```


```r
# reads processed: 2184953
# reads with at least one reported alignment: 2164338 (99.06%)
# reads that failed to align: 18292 (0.84%)
# reads with alignments suppressed due to -m: 2323 (0.11%)
# Reported 2164338 alignments to 1 output stream(s)
```



Spliced Alignment using SpliceMap
========================================================

- To generate spliced alignments (RNA-Seq)
- SpliceMap uses bowtie
- Accepts paramets as named list
- Limitation: reads should be at least 50bp


### Data
- RNA-Seq of NIH3T3 (Mus musculus) from ENCODE: https://www.encodeproject.org/experiments/ENCSR000CLW/
- Single end
- Reads from chr1 used for this exercise


Spliced Alignment using SpliceMap
========================================================

### Aligning the sequences using SpliceMap


```r
readsFiles <- "Data/ENCFF001QSC_chr1.fastq"
refDir <- "Data/chr1.fa"
indexDir <- file.path(tempdir(), "refsIndex")
samFiles <- "Data/ENCFF001QSC_chr1.sam"
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
            bowtie_base_dir="Data/mm9index_chr1_chr2",
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

<br><br>
### Sorting SAM/BAM

- Most post-processing methods requires reads sorted by aligned position
- Sorted and indexed large files enable random access of required genomic regions

### Tools
- R: Rsamtools
- Non-R: samtools, Picard


Sorting & Indexing the SAM/BAM
========================================================
<br><br>
### Convert SAM to BAM

- `asBam` converts SAM to BAM, sort and index

```r
library(Rsamtools)
insam <- "Data/CTCF_mm9_MF.sam"
outbam <- "Data/CTCF_mm9_MF"
asBam(insam, outbam)
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

```r
headerinfo <- scanBamHeader("Data/CTCF_mm9_MF.bam")

chrsinbam <- lapply(scanBamHeader("Data/CTCF_mm9_MF.bam"),function(h) {names(h$targets)})
chrsinbam <- as.vector(unlist(chrsinbam))
chrsinbam
```

```
[1] "chr1" "chr2"
```


Sorting & Indexing using samtools
========================================================


```r
samtools view -b in.sam > out.bam
samtools sort in.bam outputprefix
samtools index sorted.bam
```


Summary & Post-Alignment QC
========================================================
type: section
id: QC


Summary & Post-Alignment QC
========================================================
<br>
### Alignment Summary
- Alignment rate
- Primary and secondary alignments
- Percentage of Duplicate reads
- rRNA/mtDNA reads
- Paired reads: properly paired, unpaired, chimeric pair


```r
param1 <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE))
countBam("Data/CTCF_mm9_MF.bam",param=param1)
```

```
  space start end width            file records nucleotides
1    NA    NA  NA    NA CTCF_mm9_MF.bam 2164338    77916168
```

<br>
### Using samtools

```r
samtools flagstat in.bam
```


ChIP-Seq QC
========================================================

- Reads summary
- Duplicate percentage
- Estimation of Fragment Length
- Cross-Coverage score at the fragment length
- Cross-coverage score at the fragment length over Cross-coverage at the read length
- Percentage of reads wthin peaks
- Percentage of reads wthin Blacklist regions

## R package: ChIPQC (https://bioconductor.org/packages/release/bioc/html/ChIPQC.html)


RNA-Seq QC
========================================================

- Reads Summary
- Uniformity of 5' to 3' gene coverage bias
- Percentage of reads on exons
- Strand specificity
- Correlation between replicates

### Tools: 
- RSeQC http://rseqc.sourceforge.net/
- Picard (CollectRnaSeqMetrics)
- RNA-SeQC 


Reading BAM/SAM
========================================================
type: section
id: readbam


Reading BAM/SAM
========================================================

### Methods for reading BAM/SAM

- <b>readAligned</b> from ShortRead package 
    – Accept multiple formats – BAM, export
    - Reads all files in a directory
    - Reads base call qualities, chromosome, position, and strand
- <b>scanBam</b> from Rsamtools package
    - scanBam reads BAM files into list structure
    - Options to select what fields and which records to import using <b>ScanBamParam</b>
- <b>readGAlignments</b> and <b>readGAlignmentPairs</b> from GenomicAlignments package
    - can accept <b>ScanBamParam</b> object



Reading BAM/SAM
========================================================

<b>ScanBamParam:</b>


```r
# Constructor
ScanBamParam(flag = scanBamFlag(), what = character(0), which)

# Constructor helpers
scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA, 
    hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
    isFirstMateRead = NA, isSecondMateRead = NA, isNotPrimaryRead = NA,
    isNotPassingQualityControls = NA, isDuplicate = NA)
```

Reading BAM/SAM
========================================================

- GenomicRanges package defines the GAlignments class – a specialised class for storing set of genomic alignments (ex: sequencing data) 
- Only BAM support now – future version may include other formats
- The readGAlignments function takes an additional argument, <b>param</b> allowing the user to customise which genomic regions and which fields to read from BAM
-<b>param</b> can be constructed using </b>ScanBamParam</b> function


```r
library(GenomicAlignments)
bamin <- "Data/CTCF_mm9_MF.bam"
SampleAlign <- readGAlignments(bamin)
SampleAlign
```

```
GAlignments object with 2164338 alignments and 0 metadata columns:
            seqnames strand       cigar    qwidth     start       end
               <Rle>  <Rle> <character> <integer> <integer> <integer>
        [1]     chr1      -         36M        36   3010814   3010849
        [2]     chr1      -         36M        36   3010814   3010849
        [3]     chr1      +         36M        36   3010860   3010895
        [4]     chr1      +         36M        36   3010917   3010952
        [5]     chr1      +         36M        36   3010949   3010984
        ...      ...    ...         ...       ...       ...       ...
  [2164334]     chr2      -         36M        36 182009922 182009957
  [2164335]     chr2      -         36M        36 182009922 182009957
  [2164336]     chr2      -         36M        36 182010636 182010671
  [2164337]     chr2      -         36M        36 182012129 182012164
  [2164338]     chr2      -         36M        36 182012130 182012165
                width     njunc
            <integer> <integer>
        [1]        36         0
        [2]        36         0
        [3]        36         0
        [4]        36         0
        [5]        36         0
        ...       ...       ...
  [2164334]        36         0
  [2164335]        36         0
  [2164336]        36         0
  [2164337]        36         0
  [2164338]        36         0
  -------
  seqinfo: 2 sequences from an unspecified genome
```



Reading Sequence alignments (BAM/SAM)
========================================================

We can also customise which regions to read

```r
region <- GRanges("chr1",IRanges(3000000,5000000))
param1 <- ScanBamParam(what=c("rname", "pos", "cigar","qwidth"),which=region, flag=scanBamFlag(isUnmappedQuery=FALSE, isDuplicate=FALSE, isNotPassingQualityControls=FALSE))
SampleAlign1 <- readGAlignments(bamin,param=param1)
```



