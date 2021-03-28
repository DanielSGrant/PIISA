writeLines(c("Checking for required packges and attempting to download.",
          "If an error occurs at this stage it may be necessary to manually install the required packages.",
           "Installation instructions can be found at: https://benjjneb.github.io/dada2/dada-installation.html"))
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
#Install DESeq2
if(!requireNamespace('DESeq2',quietly = TRUE)){
  BiocManager::install("DESeq2", quietly=TRUE)  
}
library(DESeq2);


writeLines("Finished loading packages\n",)

#Setting working directory to file location
#setwd(getSrcDirectory()[1])
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
Reverse <- readline("Please enter the pattern for your reverse Files: ")

#Call in my forward and reverse reads by common names
fnFs <- sort(list.files(path=(paste(getwd(),ifolder, sep="/")),pattern=Forward, full.names = TRUE))
fnRs <- sort(list.files(path=(paste(getwd(),ifolder, sep="/")), pattern=Reverse, full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Create output folder if it doesn't exist already
dir.create(file.path(ofolder), showWarnings = FALSE)

#Put plots of quality scores in pdf
writeLines("Plotting quality scores and writing to pdf in Output folder")
pdf(file = file.path(paste(ofolder,"Quality_Scores.pdf",sep='/')) ,width = 5,height = 5)
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
  #Set windows = TRUE if you are using a windows machine, else windows=FALSE
  if(first)
  {
    win <- readline("Are you using a windows computer? [y/n]: ")
    while(win != 'y' && win != 'n')
    {
      win <- readline("Unexpected selection, please enter y if you are on a windows machine, or n if not: ")
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
  writeLines(c("","Truncation values for forward and reverse reads dictate where the reads are trimmed. Value permitted are between 0-250.",
             "If you want to trim the last 10 cycles, enter 240. Ensure forward and reverse reads maintain overlap"))
  trim1 <- readline("Please enter a right trim value for forward reads (240): ")
  trim2 <- readline("Please enter a right trim value for reverse reads (240): ")
  trim3 <- readline("Please enter a left trim value for forward reads (0): ")
  trim4 <- readline("Please enter a left trim value for reverse reads (0): ")
  maxEEF <- readline("Please enter the max expected error value for forward reads (2): ")
  maxEER <- readline("Please enter the max expected error value for reverse reads (2): ")

  #Filter and trim data
  writeLines("Performing filtering and trimming")
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       maxN=0, maxEE=c(as.numeric(maxEEF),as.numeric(maxEER)), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=!windows, truncLen=c(as.numeric(trim1),as.numeric(trim2)),
                       trimLeft = c(as.numeric(trim3), as.numeric(trim4))) 
  writeLines("Done, filtered files written to 'Filtered' directory.")
  
  #Learn errors for F reads
  writeLines("Learning error rates and plotting to Output pdf (This may take some time)")
  set.seed(100)
  errF <- learnErrors(filtFs, multithread=TRUE, verbose=FALSE)
  #Learn errors for R reads - also takes a while
  errR <- learnErrors(filtRs, multithread=TRUE, verbose = FALSE)

  #Print error plots to pdf
  pdf(file = file.path(paste(ofolder, paste(paste("Error_Plots",count, sep="_"),"pdf",sep="."), sep = '/')),width = 5,height = 5)
  print(plotErrors(errF, nominalQ=TRUE))
  print(plotErrors(errR, nominalQ=TRUE))
  dev.off()
  writeLines("Done\n")
  
  #Run dereplication before running dada
  derepFs <- derepFastq(filtFs, verbose=FALSE)
  derepRs <- derepFastq(filtRs, verbose=FALSE)
  # Name the derep-class objects by the sample names
  names(derepFs) <- sample.names
  names(derepRs) <- sample.names
  
  #Get type of data pooling from user
  p <- readline("Would you like to pool data for analysis (n) [y/n/p]: ")
  while(p != 'y' && p != 'n' && p != 'p')
  {
    p <- readline("Unexpected selection, please enter n (non-pooled), p (pooled), or p (pseudo-pooled): ")
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
  writeLines("Done\n")

  #Get parameters for merging from user
  overlap <- readline("Enter the minimum overlap for merging forward and reverse reads (12): ")
  mismatch <- readline("Enter the maximum allowed mismatch for merging forward and reverse reads (0): ")
  
  #Merge paired ends
  writeLines("Merging forward and reverse reads")
  mergers<-mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = as.numeric(overlap), maxMismatch = as.numeric(mismatch), returnRejects = FALSE, propagateCol = character(0), justConcatenate = FALSE, verbose = FALSE)
  writeLines("Done\n")
  
  #Make a sequence table and write to csv
  seqtab <- makeSequenceTable(mergers)
  write.csv(seqtab, file = paste(ofolder,"/","sequences_",count,".csv",sep=''))
  
  writeLines("Removing chimeras")
  
  #Remove chimeras adn write sequence table to csv
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
  write.csv(seqtab.nochim, file = paste(ofolder,"/sequences_nochim_",count,".csv",sep=""))
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
  write.csv(track,file=paste(ofolder,"/","Summary_",count,".csv",sep=''))
  writeLines("Done\n")
  
  writeLines("Writing summary of input parameters to txt file")
  #Output summary of results
  ofile <- paste(paste(ofolder, paste("Parameters",count,sep="_"), sep = '/'),"txt",sep=".")
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
  count = count + 1
  while(lcv != 'y' && lcv != 'n')
  {
    lcv <- readline("Unexpected selection, please enter y to run again, or n to move on: ")
  }
}
lcv = 'y'
writeLines("\n")
count = 1
while(lcv == 'y')
{
  #Set seed so it is more reproducible
  set.seed(100)

  #testing new DB making algorithm
  db <- readline("Please enter the name of the database you are using for comparison: ")
  csv <- strsplit(db, ".",fixed=-T)[[1]][1]
  csv <- paste(csv,count, sep = "")
  rc <- readline("Would you like to allow reverse compliment classification? (n) [y/n]: ")
  while(rc != 'y' && rc != 'n')
  {
    p <- readline("Unexpected selection, please enter y for RC classification, otherwise n: ")
  }
  mb <- readline("Please enter the minimum bootstrap value (60): ")
  if(rc == 'y')
  {
    taxa <- assignTaxonomy(seqtab.nochim, db, tryRC=TRUE, minBoot = as.numeric(mb))
  }
  else
  {
    taxa <- assignTaxonomy(seqtab.nochim, db, minBoot = as.numeric(mb))
  }
  #is.na <- taxa[,"Phylum"] %in% "NA"
  #taxa.trim <- taxa[!is.na,]
  taxa.print <- taxa.trim # Removing sequence rownames for display only
  rownames(taxa.print) <- NULL

  #writing results of assign taxonomy to csv file
  writeLines("Writing output of taxonomy analysis to csv")
  write.csv(taxa.print, file = paste(ofolder,paste(csv,"csv", sep='.'), sep='/'))
  writeLines("Done")
  count = count + 1
  lcv = readline("Would you like to re-run assign taxonomy? [y/n]: ")
  while(lcv != 'y' && lcv != 'n')
  {
    lcv <- readline("Unexpected selection, please enter y to run again, or n to move on: ")
  }
}

theme_set(theme_bw())
# count table:
#asv_tab <- t(seqtab.nochim)
#count_tab <- asv_tab
#tax_tab <- taxa

#deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~Location)

#count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
#tax_tab_phy <- tax_table(tax_tab)


#ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa))
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
pdf(file = file.path(paste(ofolder, "Bray_NMDS.pdf",sep="/")),width = 5,height = 5)
print(plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS"))
dev.off()

ps.bar <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
pdf(file = file.path(paste(ofolder, "Abundance_Plots.pdf",sep="/")),width = 5,height = 6, paper = "a4r")
print(plot_bar(ps.bar, fill = "Kingdom"))
print(plot_bar(ps.bar, fill = "Phylum"))
print(plot_bar(ps.bar, fill = "Class"))
print(plot_bar(ps.bar, fill = "Order"))
print(plot_bar(ps.bar, fill = "Family"))
print(plot_bar(ps.bar, fill = "Genus"))
print(plot_bar(ps.bar, fill = "Species"))
dev.off()






# and now we can call the plot_richness() function on our phyloseq object
pdf(file = file.path(paste(ofolder, "Diversity_Plots.pdf",sep="/")),width = 5,height = 5)
print(plot_richness(ps, measures=c("Chao1", "Shannon")))
dev.off()




