# SeizurePaper_1
Scripts used for RNA alignment for seizure study

The following scripts were used to convert Ion Torrent uBam files to fastq files, clean with Trim-Galore, and then align to the Human Reference genome (Grch38) using 2-pass STAR and Bowtie2.  Scripts are designed to be run on ubam files from multiple chips. The final Samtools cript merges the data together.  

Scripts are designed to run in a Project folder, with subsequent folders labeled ANALYSIS, DATA, and Scripts.  DATA contains all of the ubam files. ANALYSIS is a by sample series of folders containing intermediates.  All scripts reside in Scripts.  

The scripts are called from the project folder using a list of sample identifiers (usually unique sample IDS at the begining of the file name)

Scripts are called with...
>cat samples.txt | parallel -j n Scripts/<name_of_Script>
  where n is the number of threads depending on your systems available memory. 
 
