# BIOL4315_Lab2 = Quality Check, Processing & Alignment

# install packages
install.packages("BiocManager")
BiocManager::install("Rqc")
BiocManager::install("quasr")
install.packages("here")

# install.packages("Hmisc")

# Aggregating multiple fastqc reports into a data frame
# install.packages("fastqcr")

# Load the package
library(Rqc)
library(fastqcr) 
library(QuasR)
library(ShortRead)
library(here)

# Get the path to the file
folder <- system.file(package="ShortRead", "extdata/E-MTAB-1147")

# Feeds fastq.qz files in "folder" to quality check function
qcRes <- rqc(path = folder, pattern = ".fastq.gz", openBrowser=FALSE, outdir="outputs")


# Sequencing quality per base/cycle
rqcCycleQualityBoxPlot(qcRes)

# Sequence content per base/cycle
rqcCycleBaseCallsLinePlot(qcRes)

# Read frequency plot
rqcReadFrequencyPlot(qcRes)

# Demo QC directory containing zipped FASTQC reports
qc.dir <- system.file("fastqc_results", package = "fastqcr")
qc <- qc_aggregate(qc.dir)
qc

# Inspecting QC problems

# See which modules failed in the most samples
qc_fails(qc, "module")
# Or, see which samples failed the most
qc_fails(qc, "sample")

# Building multi QC reports
qc_report(qc.dir, result.file = "outputs/multi-qc-report" )

# Building one-sample QC reports (+ interpretation)
qc.file <- system.file("fastqc_results", "S1_fastqc.zip", package = "fastqcr")

# View the report rendered by R functions
qc_report(qc.file, result.file = "outputs/one-sample-report",
          interpret = TRUE)

# QC for long reads
# conda install wget
# wget -nc -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR153/006/ERR1539006/ERR1539006.fastq.gz

# Run sequali on downloaded file
# sequali --outdir outputs/sequali_reports/ data/ERR1539006.fastq.gz

# Filtering and trimming short reads
# Obtain a list of fastq file paths
fastqFiles <- system.file(package="ShortRead",
                          "extdata/E-MTAB-1147",
                          c("ERR127302_1_subset.fastq.gz",
                            "ERR127302_2_subset.fastq.gz")
)

# Defined processed fastq file names in the outputs folder
outfiles <- paste0("outputs/", c("processed_1_", "processed_2_"), basename(fastqFiles))

# Process fastq files
# Remove reads that have more than 1 N, (nbases)
# Trim 3 bases from the end of the reads (truncateendbases)
# Remove ACCCGGGA patern if it occurs at the start (lpattern)
# Remove reads shorter than 40 base-pairs (minlength)
preprocessReads(fastqFiles, outfiles, 
                nBases=1,
                truncateEndBases=3,
                Lpattern="ACCCGGGA",
                minLength=40)

# Obtain a list of fastq file paths
fastqFile <- system.file(package="ShortRead",
                         "extdata/E-MTAB-1147",
                         "ERR127302_1_subset.fastq.gz")
# Read fastq file
fq = readFastq(fastqFile)

# Get quality scores per base as a matrix
qPerBase = as(quality(fq), "matrix")

# Get number of bases per read that have quality score below 20
# We use this
qcount = rowSums( qPerBase <= 20) 

# Number of reads where all phred scores >= 20
fq[which(qcount == 0)]

#We can finally write out the filtered fastq file with the ShortRead::writeFastq() function

#mode = 'a' allows you to rewrite files that are already written
ShortRead::writeFastq( fq[which(qcount == 0)], here::here("outputs/Qfiltered3.fastq"), mode = 'a', compress = FALSE) 
#~/Desktop/ASharmaBIOl4315_R/BIOL4315_Lab2/outputs

# Set up streaming with block size 1000
# Every time we call the yield() function 1000 read portion
# Of the file will be read successively.
f <- FastqStreamer(fastqFile,readerBlockSize=1000) 

# We set up a while loop to call yield() function to
# Go through the file
while(length(fq <- yield(f))) {
  
  # remove reads where all quality scores are < 20 
  # get quality scores per base as a matrix
  qPerBase = as(quality(fq), "matrix")
  
  # get number of bases per read that have Q score < 20
  qcount = rowSums( qPerBase <= 20) 
  
  # write fastq file with mode="a", so every new block
  # is written out to the same file
  writeFastq(fq[which(qcount == 0)], 
             here::here("outputs", paste(basename(fastqFile), "Qfiltered3.fastq", sep="_")), 
             mode="a", compress = FALSE)
}

# Take Illumina reads from Q8 and run rqc() function

# Feeds fastq.qz files to quality check function
qcres <- rqc(path = "data/", pattern = "^ERR11203340_.*\\.fastq\\.gz$", pair = c(1,1), openBrowser=FALSE, outdir="outputs")

# Sequencing quality per base/cycle
rqcCycleQualityBoxPlot(qcres)

# Sequence content per base/cycle
rqcCycleBaseCallsLinePlot(qcres)

# Read frequency plot
rqcReadFrequencyPlot(qcres)

