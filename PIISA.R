cat("\014")
writeLines(c("\nWelcome to PIISA, an interactive Pipeline for Iterative and Interactive Sequence analysis!",
             "Reminder that throughout the program you will be prompted for inputs! Type these inputs from the keyboard and press enter to continue.",
             "At any input stage in the program enter q to quit the progam, (Please note that an error message will be displayed upon quitting).",
             "At prompts suggested values are enclosed in parentheses (), and lists of all options are enclosed in square brackets [].",
             "To select a default value simply hit enter without typing anything."))
start <- readline("Press enter to continue: ")
if(start == "q"){stop()}
while(start != "")
{
  start <- readline("Try again, press enter without typing anything: ")
  if(start == "q"){stop()}
}

writeLines(c(" ","Checking for required packges and attempting to download.",
          "If an error occurs at this stage it may be necessary to manually install the required packages."
          ,"Please refer to FAQ in manual for more information if this occurs."))
#Check for all neccessary packages
#install Bioconductor
if(!requireNamespace('BiocManager',quietly = TRUE)){
  install.packages('BiocManager')
}
library(BiocManager)
#Install dada2
if(!requireNamespace('dada2',quietly = TRUE)){
  BiocManager::install("dada2", version = "3.12")  
}
library(dada2); packageVersion("dada2")
#install DECIPHER
if(!requireNamespace('DECIPHER',quietly = TRUE)){
  BiocManager::install("DECIPHER", quietly=TRUE)  
}
library(DECIPHER); packageVersion("DECIPHER")
#Load ggplot
if(!requireNamespace('ggplot2',quietly = TRUE)){
  install.packages('ggplot2')  
}
library(ggplot2); packageVersion("ggplot2")
#Load phyloseq
if(!requireNamespace('phyloseq',quietly = TRUE)){
  BiocManager::install("phyloseq", quietly=TRUE)  
}
library(phyloseq); packageVersion("phyloseq")

writeLines("Finished loading packages\n",)

#Setting working directory to file location
script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir)

#Name of folder for input files
ifolder = "Input"
#Name of folder for output files
ofolder = paste(getwd(),"Output",sep='/')

#Get naming pattern for forward and reverse reads
writeLines(c("Please enter the pattern for forward and reverse files. For example, forward files:", 
  "'Tree550mcrA_R1.fastq.sanger.gz' and 'Well2mcrA_R1.fastq.sanger.gz'", "Have pattern '_R1.fastq.'"))
Forward <- readline("Please enter the pattern for your forward Files: ")
if(Forward == "q"){stop()}
fnFs <- sort(list.files(path=(paste(getwd(),ifolder, sep="/")),pattern=Forward, full.names = TRUE))
while(identical(fnFs,character(0)))
{
  Forward <- readline("Error opening files, please check that files are in 'Input' folder and spelling is corect and try again: ")
  if(Forward == "q"){stop()}
  fnFs <- sort(list.files(path=(paste(getwd(),ifolder, sep="/")),pattern=Forward, full.names = TRUE))
}
Reverse <- readline("Please enter the pattern for your reverse Files: ")
if(Reverse == "q"){stop()}
fnRs <- sort(list.files(path=(paste(getwd(),ifolder, sep="/")), pattern=Reverse, full.names = TRUE))
while(identical(fnRs,character(0)))
{
  Reverse <- readline("Error opening files, please check that files are in 'Input' folder and spelling is corect and try again: ")
  if(Reverse == "q"){stop()}
  fnRs <- sort(list.files(path=(paste(getwd(),ifolder, sep="/")), pattern=Reverse, full.names = TRUE))
}

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Create output folder if it doesn't exist already
dir.create(file.path(ofolder), showWarnings = FALSE)

writeLines("\nTaking inputs for quality scores plot.")
#Get parameters for plots
#Select plot width
w <- readline("Please enter the desired width for the quality scores plot in inches (8): ")
if(w == "q"){stop()} else if(w == ""){w = 8}
while(as.numeric(w) < 0 || as.numeric(w) > 50)
{
  w <- readline("Invalid entry, please enter a number for plot width in inches (8): ")
  if(w == "q"){stop()} else if(w == ""){w = 8}
}
#Select plot height
h <- readline("Please enter the desired height for the quality scores plot in inches (8): ")
if(h == "q"){stop()}else if(h == ""){h = 8}
while(as.numeric(h) < 0 || as.numeric(h) > 50)
{
  h <- readline("Invalid entry, please enter a number for plot height in inches (8): ")
  if(h == "q"){stop()} else if(h == ""){w = 8}
}
#select font family
f <- readline("Please enter the desired font for quality score plots (Helvetica): ")
if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
while(!(f %in% names(pdfFonts())))
{
  disp <- readline("Error, invalid font. Would you like to see a list of valid fonts? [y/n]: ")
  if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
  while(disp != 'y' && disp != 'n')
  {
    disp <- readline("Error, invalid selection. Enter y to see available fonts of n to enter font: ")
    if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
  }
  f <- readline("Please enter the desired font for quality score plots (Helvetica): ")
  if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
}
#Select paper size
p <- readline("Would you like generated pdf's size to be 8.5\"x11\"? (y) [y/n]: ")
if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "default"} else if(p == 'n'){p = "special"}
while(p != "special" && p != "default")
{
  p <- readline("Unexpected entry, please enter y for letter paper, or n for custom sizing: ")
  if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "default"} else if(p == 'n'){p = "special"}
}

#Put plots of quality scores in pdf
writeLines("\nPlotting quality scores and writing to pdf in Output folder")
pdf(file = file.path(paste(ofolder,"Quality_Scores.pdf",sep='/')) ,width = as.numeric(w),height = as.numeric(h), 
    family = f, paper = p)
#Check quality scores of forward reads
print(plotQualityProfile(fnFs))
#Check quality scores of reverse reads
print(plotQualityProfile(fnRs))
dev.off()
writeLines("Done")

# Place filtered files in filtered/ subdirectory
filtFs <- file.path("./Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
writeLines(" ")

#Enter first loop for trimming, dada analysis, merging, and chimera removal
lcv = 'y'
first = TRUE;
count <- 1
while(lcv != 'n')
{
  ftd.folder <- paste("FilterTrimDada_Run",count,sep = "")
  #Create folder if it doesn't exist already
  dir.create(file.path(paste(ofolder,ftd.folder, sep = "/")), showWarnings = FALSE)
  
  #Set windows = TRUE if you are using a windows machine, else windows=FALSE
  if(first)
  {
    win <- readline("Are you using a windows computer? [y/n]: ")
    if(win == "q"){stop()}
    while(win != 'y' && win != 'n')
    {
      win <- readline("Unexpected selection, please enter y if you are on a windows machine, or n if not: ")
      if(win == "q"){stop()}
    }
    if(win == 'y')
    {
      windows = TRUE
    }
    else
    {
      windows = FALSE;
    }
  }
  
  #Select value for trimming low quality scores based on quality score plots
  writeLines(c("","Truncation values for forward and reverse reads dictate where the reads are trimmed on the right.",
             "These values should be based on quality scores. Ensure forward and reverse reads maintain overlap"))
  trim1 <- readline("Please enter a truncation value for forward reads (240): ")
  if(trim1 == "q"){stop()}
  else if(trim1 == ""){trim1 = 240}
  while(as.numeric(trim1) < 0 || as.numeric(trim1) > 300)
  {
    trim1 <- readline("Invalid entry, please enter a number based on quality scores plot: ")
    if(trim1 == "q"){stop()}
    else if(trim1 == ""){trim1 = 240}
  }
  trim2 <- readline("Please enter a truncation value for reverse reads (240): ")
  if(trim2 == "q"){stop()}
  else if(trim2 == ""){trim2 = 240}
  while(as.numeric(trim2) < 0 || as.numeric(trim2) > 300)
  {
    trim2 <- readline("Invalid entry, please enter a number based on quality scores plot: ")
    if(trim2 == "q"){stop()}
    else if(trim2 == ""){trim2 = 240}
  }
  writeLines(c("","Left trim values for forward and reverse reads dictate where the reads are trimmed.",
               "If your primers have not been removed yet enter trim values equal to primer length."))
  trim3 <- readline("Please enter a left trim value for forward reads (0): ")
  if(trim3 == "q"){stop()}
  else if(trim3 == ""){trim3 = 0}
  while(as.numeric(trim3) < 0 || as.numeric(trim3) > 80)
  {
    trim3 <- readline("Invalid entry, please enter 0 if your primers are removed, or the length of the primer in nucleotides if not: ")
    if(trim3 == "q"){stop()}
    else if(trim3 == ""){trim3 = 0}
  }
  trim4 <- readline("Please enter a left trim value for reverse reads (0): ")
  if(trim4 == "q"){stop()}
  else if(trim4 == ""){trim4 = 0}
  while(as.numeric(trim4) < 0 || as.numeric(trim4) > 80)
  {
    trim4 <- readline("Invalid entry, please enter 0 if your primers are removed, or the length of the primer in nucleotides if not: ")
    if(trim4 == "q"){stop()}
    else if(trim4 == ""){trim4 = 0}
  }
  writeLines(c("","Max expected error values dictate how many error we expect for each read.",
               "MaxEE values can be increased for lower quality reads or decreased for higher quality reads."))
  maxEEF <- readline("Please enter the max expected error value for forward reads (2): ")
  if(maxEEF == "q"){stop()}
  else if(maxEEF == ""){maxEEF = 2}
  while(as.numeric(maxEEF) < 0 || as.numeric(maxEEF) > 30)
  {
    maxEEF <- readline("Invalid entry, please enter a positive integer value: ")
    if(maxEEF == "q"){stop()}
    else if(maxEEF == ""){maxEEF = 2}
  }
  maxEER <- readline("Please enter the max expected error value for reverse reads (2): ")
  if(maxEER == "q"){stop()}
  else if(maxEER == ""){maxEER = 2}
  while(as.numeric(maxEER) < 0 || as.numeric(maxEER) > 30)
  {
    maxEER <- readline("Invalid entry, please enter a positive integer value: ")
    if(maxEER == "q"){stop()}
    else if(maxEER == ""){maxEER = 2}
  }

  #Filter and trim data
  writeLines("\nPerforming filtering and trimming (This may take some time)")
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       maxN=0, maxEE=c(as.numeric(maxEEF),as.numeric(maxEER)), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=!windows, truncLen=c(as.numeric(trim1),as.numeric(trim2)),
                       trimLeft = c(as.numeric(trim3), as.numeric(trim4))) 
  writeLines("Done, filtered files written to 'Filtered' directory.")
  
  #Learn errors for F reads
  writeLines("\nLearning error rates and plotting to Output pdf (This may take some time)")
  set.seed(100)
  writeLines("Forward reads")
  errF <- learnErrors(filtFs, multithread=TRUE, verbose=FALSE)
  #Learn errors for R reads - also takes a while
  writeLines("Reverse reads")
  errR <- learnErrors(filtRs, multithread=TRUE, verbose = FALSE)

  writeLines("\nTaking inputs for error plots.")
  #Get parameters for error plots
  #select plot width
  w <- readline("Please enter the desired width for the error plots in inches (8.5): ")
  if(w == "q"){stop()} else if(w == ""){w = 8.5}
  while(as.numeric(w) < 0 || as.numeric(w) > 50)
  {
    w <- readline("Invalid entry, please enter a number for plot width in inches (8.5): ")
    if(w == "q"){stop()} else if(w == ""){w = 8.5}
  }
  #Select plot height
  h <- readline("Please enter the desired height for the error plots in inches (11): ")
  if(h == "q"){stop()}else if(h == ""){h = 11}
  while(as.numeric(h) < 0 || as.numeric(h) > 50)
  {
    h <- readline("Invalid entry, please enter a number for plot height in inches (11): ")
    if(h == "q"){stop()} else if(h == ""){w = 11}
  }
  #Select font
  f <- readline("Please enter the desired font for error plots (Helvetica): ")
  if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
  while(!(f %in% names(pdfFonts())))
  {
    disp <- readline("Error, invalid font. Would you like to see a list of valid fonts? [y/n]: ")
    if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
    while(disp != 'y' && disp != 'n')
    {
      disp <- readline("Error, invalid selection. Enter y to see available fonts of n to enter font: ")
      if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
    }
    f <- readline("Please enter the desired font for error plots (Helvetica): ")
    if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
  }
  #Select paper size (letter or custom)
  p <- readline("Would you like generated pdf's size to be 8.5\"x11\"? (y) [y/n]: ")
  if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "default"} else if(p == 'n'){p = "special"}
  while(p != "special" && p != "default")
  {
    p <- readline("Unexpected entry, please enter y for letter paper, or n for custom sizing: ")
    if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "default"} else if(p == 'n'){p = "special"}
  }
  
  #Print error plots to pdf
  pdf(file = file.path(paste(paste(ofolder,ftd.folder, sep = "/"), paste(paste("Error_Plots",count, sep="_"),"pdf",sep="."), sep = '/')),
      width = as.numeric(w),height = as.numeric(h), family = f, paper = p)
  print(plotErrors(errF, nominalQ=TRUE))
  print(plotErrors(errR, nominalQ=TRUE))
  dev.off()
  writeLines("Done")
  
  #Run dereplication before running dada
  writeLines("\nDereplicating data")
  derepFs <- derepFastq(filtFs, verbose=FALSE)
  derepRs <- derepFastq(filtRs, verbose=FALSE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  writeLines("Done")
  
  #Get type of data pooling from user
  p <- readline("Would you like to pool data for analysis (n) [y/n/p]: ")
  if(p == "q"){stop()}
  else if(p == ""){p = "n"}
  while(p != 'y' && p != 'n' && p != 'p')
  {
    p <- readline("Unexpected selection, please enter n (non-pooled), p (pooled), or p (pseudo-pooled): ")
    if(p == "q"){stop()}
    else if(p == ""){p = "n"}
  }
  if(p == 'y'){
    pooled <= TRUE
  } else if(p == 'n'){
    pooled <- FALSE
  } else if(p == 'p'){
    pooled <- "pseudo"
  }
  writeLines("Running dada2 algorithm on forward and reverse reads")
  
  #Apply the core sample interference algorithm on forward and reverse reads
  dadaFs <- dada(filtFs, err=errF, pool = pooled, multithread=TRUE)
  writeLines(" ")
  dadaRs <- dada(filtRs, err=errR, pool = pooled, multithread=TRUE)
  writeLines("\nDone")
  
  writeLines("\nForward reads dada summary:")
  print(dadaFs[[1]])
  writeLines("\nReverse reads dada summary:")
  print(dadaRs[[1]])
  writeLines("")

  #Get parameters for merging from user
  overlap <- readline("Enter the minimum overlap for merging forward and reverse reads (12): ")
  if(overlap == "q"){stop()}
  else if(overlap == ""){overlap = 12}
  while(as.numeric(overlap) < 0 || as.numeric(overlap) > 80)
  {
    overlap <- readline("Invalid entry, please enter a positive integer value: ")
    if(overlap == "q"){stop()}
    else if(overlap == ""){overlap = 12}
  }
  mismatch <- readline("Enter the maximum allowed mismatch for merging forward and reverse reads (0): ")
  if(mismatch == "q"){stop()}
  else if(mismatch == ""){mismatch = 0}
  while(as.numeric(mismatch) < 0 || as.numeric(mismatch) > 50)
  {
    mismatch <- readline("Invalid entry, please enter a positive integer value: ")
    if(mismatch == "q"){stop()}
    else if(mismatch == ""){mismatch = 0}
  }
  
  #Merge paired ends
  writeLines("Merging forward and reverse reads")
  mergers<-mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = as.numeric(overlap), maxMismatch = as.numeric(mismatch), 
                      returnRejects = FALSE, propagateCol = character(0), justConcatenate = FALSE, verbose = FALSE)
  writeLines("Done\n")
  
  #Make a sequence table and write to csv
  seqtab <- makeSequenceTable(mergers)
  write.csv(seqtab, file = paste(ofolder,"/",ftd.folder,"/","sequences",".csv",sep=''))
  
  writeLines("Removing chimeras")
  
  #Remove chimeras and write sequence table to csv
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  write.csv(seqtab.nochim, file = paste(paste(ofolder,ftd.folder, sep = "/"),"/sequences_nochim",".csv",sep=""))
  writeLines("Done")
  
  #Check the outcome of removing chimeras
  fraction = sum(seqtab.nochim)/sum(seqtab)
  writeLines(paste("The fraction of sequences remaining after removing chimeras is", fraction, sep=" "))
  

  #As a final check, look at the number of reads that made it through the pipeline at each step
  writeLines("Writing summary of analysis to csv")
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sample.names
  write.csv(track,file=paste(paste(ofolder,ftd.folder, sep = "/"),"/","Summary",".csv",sep=''))
  writeLines("Done\n")
  
  writeLines("Writing summary of input parameters to txt file")
  #Output summary of results
  ofile <- paste(paste(paste(ofolder,ftd.folder, sep = "/"),"Parameters.txt", sep = "/"))
  sink(ofile)
  cat(paste("Forward reads truncation value:", trim1, sep=" "))
  cat("\n")
  cat(paste("Reverse reads truncation value:", trim2, sep=" "))
  cat("\n")
  cat(paste("Forward reads left trim value:", trim3, sep=" "))
  cat("\n")
  cat(paste("Reverse reads left trim value:", trim4, sep=" "))
  cat("\n")
  cat(paste("Forward reads maxEE:", maxEEF, sep=" "))
  cat("\n")
  cat(paste("Reverse reads maxEE:", maxEER, sep=" "))
  cat("\n")
  if(p == 'y'){
    cat("Pooled = TRUE\n")
  } else if(p == 'n'){
    cat("Pooled = FALSE\n")
  } else if(p == 'p'){
    cat("Pooled = PSEUDO\n")
  }
  cat("\nForward reads dada summary:\n")
  print(dadaFs[[1]])
  cat("\nReverse reads dada summary:\n")
  print(dadaRs[[1]])
  cat("\n")
  cat(paste("Minimum overlap value for merging:", overlap,sep=" "))
  cat("\n")
  cat(paste("Maximum mismatch value for merging:", mismatch,sep=" "))
  cat("\n")
  cat(paste("Results of removing chimeras:", sum(seqtab.nochim), "non chimeric seqs/", sum(seqtab), "original seqs =",fraction, sep=" "))
  cat("\n")
  sink()
  writeLines("Done")
  
  first = FALSE
  lcv = readline("Would you like to re-run error analysis and dada algorithm step? [y/n]: ")
  if(lcv == "q"){stop()}
  count = count + 1
  while(lcv != 'y' && lcv != 'n')
  {
    lcv <- readline("Unexpected selection, please enter y to run again, or n to move on: ")
    if(lcv == "q"){stop()}
  }
}
lcv = 'y'
writeLines(" ")
count = 1
while(lcv == 'y')
{
  tax.folder <- paste(ofolder,"/","AssnTax_Run",count,sep = "")
  #Create folder if it doesn't exist already
  dir.create(file.path(tax.folder), showWarnings = FALSE)
  #Set seed so it is more reproducible
  set.seed(100)

  #testing new DB making algorithm
  db <- readline("Please enter the name of the database you are using for comparison: ")
  if(db == "q"){stop()}
  while(!file.exists(db))
  {
    db <- readline("Error, unable to open file, please check spelling and try again: ")
    if(db == "q"){stop()}
  }
  
  #Create name for csv file
  csv <- strsplit(db, ".",fixed=-T)[[1]][1]
  csv <- paste(csv,count, sep="")
  #Prompt for assignTaxonomy parameters
  rc <- readline("Would you like to allow reverse compliment classification? (n) [y/n]: ")
  if(rc == "q"){stop()}
  else if(rc == ""){rc = "n"}
  while(rc != 'y' && rc != 'n')
  {
    p <- readline("Unexpected selection, please enter y for RC classification, otherwise n: ")
    if(rc == "q"){stop()}
    else if(rc == ""){rc = "n"}
  }
  mb <- readline("Please enter the minimum bootstrap value (50): ")
  if(mb == "q"){stop()}
  else if(mb == ""){mb = 50}
  while(as.numeric(mb) < 0 || as.numeric(mb) > 100)
  {
    mb <- readline("Error, invalid entry. Please enter a value from 0-100: ") 
    if(mb == "q"){stop()}
    else if(mb == ""){mb = 50}
  }
  writeLines("\nAssigning taxonomy.")
  if(rc == 'y')
  {
    taxa <- assignTaxonomy(seqtab.nochim, db, tryRC=TRUE, minBoot = as.numeric(mb))
  }
  else
  {
    taxa <- assignTaxonomy(seqtab.nochim, db, minBoot = as.numeric(mb))
  }
  writeLines("Done")
  taxa.print <- taxa # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL

  #writing results of assign taxonomy to csv file
  writeLines("\nWriting output of taxonomy analysis to csv")
  write.csv(taxa.print, file = paste(tax.folder,paste(csv,"csv", sep='.'), sep='/'))
  writeLines("Done")
  count = count + 1
  lcv = readline("Would you like to re-run assign taxonomy? [y/n]: ")
  if(lcv == "q"){stop()}
  while(lcv != 'y' && lcv != 'n')
  {
    lcv <- readline("Unexpected selection, please enter y to run again, or n to move on: ")
    if(lcv == "q"){stop()}
  }
}

theme_set(theme_bw())


#Get parameters for plots
writeLines("\nTaking inputs for abundance plots.")
#select plot width
w <- readline("Please enter the desired width for the abundance plots in inches (11): ")
if(w == "q"){stop()} else if(w == ""){w = 11}
while(as.numeric(w) < 0 || as.numeric(w) > 50)
{
  w <- readline("Invalid entry, please enter a number for plot width in inches (11): ")
  if(w == "q"){stop()} else if(w == ""){w = 11}
}
#Select plot height
h <- readline("Please enter the desired height for the abundance plots in inches (8.5): ")
if(h == "q"){stop()}else if(h == ""){h = 8.5}
while(as.numeric(h) < 0 || as.numeric(h) > 50)
{
  h <- readline("Invalid entry, please enter a number for plot height in inches (8.5): ")
  if(h == "q"){stop()} else if(h == ""){w = 8.5}
}
#Select font
f <- readline("Please enter the desired font for abundance plots (Helvetica): ")
if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
while(!(f %in% names(pdfFonts())))
{
  disp <- readline("Error, invalid font. Would you like to see a list of valid fonts? [y/n]: ")
  if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
  while(disp != 'y' && disp != 'n')
  {
    disp <- readline("Error, invalid selection. Enter y to see available fonts of n to enter font: ")
    if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
  }
  f <- readline("Please enter the desired font for abundance plots (Helvetica): ")
  if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
}
#Select paper size (letter or custom)
p <- readline("Would you like generated pdf's size to be landscape 8.5\"x11\"? (y) [y/n]: ")
if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "a4r"} else if(p == 'n'){p = "special"}
while(p != "special" && p != "a4r")
{
  p <- readline("Unexpected entry, please enter y for landscape letter paper, or n for custom sizing: ")
  if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "a4r"} else if(p == 'n'){p = "special"}
}

#Print abundance plots to pdf
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))
ps.bar <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
pdf(file = file.path(paste(ofolder, "Abundance_Plots.pdf",sep="/")),
    width = as.numeric(w),height = as.numeric(h), family = f, paper = p)
print(plot_bar(ps.bar, fill = "Kingdom") + geom_bar(aes(color = Kingdom, fill = Kingdom), colour='black', stat = "identity", position = "stack") +
        labs(x = "", y = "Relative Abundance\n") + scale_fill_brewer(palette = "Paired") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
print(plot_bar(ps.bar, fill = "Phylum") + geom_bar(aes(color = Phylum, fill = Phylum), colour='black', stat = "identity", position = "stack") +
        labs(x = "", y = "Relative Abundance\n") + scale_fill_brewer(palette = "Paired") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
print(plot_bar(ps.bar, fill = "Class") + geom_bar(aes(color = Class, fill = Class), colour='black', stat = "identity", position = "stack") +
        labs(x = "", y = "Relative Abundance\n") + scale_fill_brewer(palette = "Paired") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
print(plot_bar(ps.bar, fill = "Order") + geom_bar(aes(color = Order, fill = Order), colour='black', stat = "identity", position = "stack") +
        labs(x = "", y = "Relative Abundance\n") + scale_fill_brewer(palette = "Paired") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
print(plot_bar(ps.bar, fill = "Family") + geom_bar(aes(color = Family, fill = Family), colour='black', stat = "identity", position = "stack") +
        labs(x = "", y = "Relative Abundance\n") + scale_fill_brewer(palette = "Paired") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
print(plot_bar(ps.bar, fill = "Genus") + geom_bar(aes(color = Genus, fill = Genus), colour='black', stat = "identity", position = "stack") +
        labs(x = "", y = "Relative Abundance\n") + scale_fill_brewer(palette = "Paired") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
print(plot_bar(ps.bar, fill = "Species") + geom_bar(aes(color = Species, fill = Species), colour='black', stat = "identity", position = "stack") +
        labs(x = "", y = "Relative Abundance\n") + scale_fill_brewer(palette = "Paired") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")))
dev.off()


#Get parameters for plots
writeLines("\nTaking inputs for diversity plots.")
#select plot width
w <- readline("Please enter the desired width for the diversity plots in inches (6): ")
if(w == "q"){stop()} else if(w == ""){w = 6}
while(as.numeric(w) < 0 || as.numeric(w) > 50)
{
  w <- readline("Invalid entry, please enter a number for plot width in inches (6): ")
  if(w == "q"){stop()} else if(w == ""){w = 6}
}
#Select plot height
h <- readline("Please enter the desired height for the diversity plots in inches (6): ")
if(h == "q"){stop()}else if(h == ""){h = 6}
while(as.numeric(h) < 0 || as.numeric(h) > 50)
{
  h <- readline("Invalid entry, please enter a number for plot height in inches (6): ")
  if(h == "q"){stop()} else if(h == ""){w = 6}
}
#Select font
f <- readline("Please enter the desired font for diversity plots (Helvetica): ")
if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
while(!(f %in% names(pdfFonts())))
{
  disp <- readline("Error, invalid font. Would you like to see a list of valid fonts? [y/n]: ")
  if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
  while(disp != 'y' && disp != 'n')
  {
    disp <- readline("Error, invalid selection. Enter y to see available fonts of n to enter font: ")
    if(disp == "q"){stop()} else if(disp == 'y'){print(names(pdfFonts()))}
  }
  f <- readline("Please enter the desired font for diversity plots (Helvetica): ")
  if(f == "q"){stop()} else if(f == ""){f = "Helvetica"}
}
#Select paper size (letter or custom)
p <- readline("Would you like generated pdf's size to be 8.5\"x11\"? (y) [y/n]: ")
if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "default"} else if(p == 'n'){p = "special"}
while(p != "special" && p != "default")
{
  p <- readline("Unexpected entry, please enter y for letter paper, or n for custom sizing: ")
  if(p == "q"){stop()} else if(p == "" || p == 'y'){p = "default"} else if(p == 'n'){p = "special"}
}

#Print diversity plots to pdf
# and now we can call the plot_richness() function on our phyloseq object
pdf(file = file.path(paste(ofolder, "Diversity_Plots.pdf",sep="/")),
    width = as.numeric(w),height = as.numeric(h), family = f, paper = p)
print(plot_richness(ps, measures=c("Simpson", "Shannon"))) #add simpson instead of chao1, observed vs chao1
print(plot_richness(ps, measures=c("Observed", "Chao1")))
dev.off()




