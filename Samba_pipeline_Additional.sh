#!/bin/bash
#SBATCH -p serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=23:55:00
#SBATCH --mem=118G
#SBATCH -o Random_sample_SAMBAR.%J.out
#SBATCH -e Random_sample_SAMBAR.%J.err

module load all
module load gencore
module load gencore_variant_detection
module load gencore_variant_detection/2.0
#module load samtools/intel/1.6

############################################################################### File conversion to prepare the input files ##########################################################

#(1) Convert from vcf to PED/MAP using PLINK:#From Venu: 

/scratch/mv83/Software/plink_beta5.3/plink --vcf L_genotype.vcf --recode --out Xenopus_LGenome --allow-extra-chr --no-pheno --no-sex --biallelic-only strict

#(2) Check your files: Does the number of lines in the PED file correspond with the number of individuals? Does the number of lines in the MAP file correspond with the number of SNPs?

#(3) Convert from PED/MAP to RAW/BIM format using PLINK:

/scratch/mv83/Software/plink_beta5.3/plink --file Xenopus_LGenome --chr-set 95 --allow-extra-chr --make-bed --recode A --out Xenopus_LGenome


############################################################################### SambaR Basic Run ##########################################################################################

# activate the environment on you command line before running this shell script

conda activate /scratch/da2451/miniconda/SambaR

R 

source("/scratch/da2451/Xenopus/Boissinotlab/gvcf/Xenopus_laevis/All-data-pop/gvcfs/Analysis_AllData/SambaR/All_SNPs/SambaR/SAMBAR_v1.07.txt")

getpackages()
 

#setwd("C:/path/to/workdir/") #where you stored your RAW/BIM files

#To check whether your input files are indeed present, type: 
list.files()  #If your input files are present, import your data using the ‘importdata’ function. 

#Specify the order and color of populations  (for example geographical order)

popvector <- c("blue", "green", "purple", "yellow")
colvector <- c("deepskyblue","darkgreen","darkviolet","gold") 
poporder <- c("blue", "green", "purple", "yellow")

importdata(inputprefix="Xenopus_LGenome", colourvector=colvector, pop_order=poporder, sumstatsfile=FALSE, depthfile=FALSE, geofile= "geofile.txt", samplefile= "popfile.txt")

#importdata(inputprefix="XenopusSamba2", colourvector=colvector, pop_order=poporder, sumstatsfile=TRUE, depthfile=TRUE, geofile= "geofile.txt") 

#importdata(inputprefix="NonHybridSamba2", colourvector=colvector, pop_order=poporder, sumstatsfile=FALSE, depthfile=FALSE, geofile= "geofile_nonHybrids.txt", samplefile= "popfile_nonHybrids.txt")


filterdata(indmiss=0.25,snpmiss=0.05,min_mac=2,dohefilter=TRUE,min_spacing=500,nchroms=NULL,TsTvfilter=NULL) 

#filterdata(indmiss=0.5,snpmiss=0.1,min_mac=2,nchroms=NULL,TsTvfilter=NULL) 

#filterdata(indmiss=0.5,snpmiss=0.5,min_mac=2,dohefilter=FALSE,nchroms=NULL,TsTvfilter=NULL) 

# Number of chromosomes (Page -26) of the reference genome assembly (optional). If a value is provided, SambaR will generate additional plots showing the number of SNPs per chromosome, and try to determine the X-chromosome. 

exportsambarfiles()  # To create the input files needed for other Additional analysis.

#STRUCTURE Analysis: K is set based on your population assumption 

findstructure(Kmax=6,add_legend=TRUE,legend_pos="topleft",legend_cex=2) 

createmaps()
#create_sambarmaps(K_max=6,radius_ratio=50)

#1# If the columns longitude and latitude are present in the inds dataframe, it will also produce geographical maps. 
calckinship()

calcdistance()

calcdiversity(nrsites=NULL) 

nrow(snps)
nrow(snps[snps$filter,])

############################################################################### SambaR Selection Run ##########################################################################################


#Selection Analysis: selectionanalyses function runs up to 4 selection scans: Fsthet, GWDS, OutFlank, and PCadapt. - Page 43
EXAMPLES:
#To run selection analyses for the metapopulation: 
selectionanalyses(do_meta=TRUE ,export="pdf",do_fsthet=TRUE,my_correction=NULL,export_data=TRUE) 
#To run selection analyses for all pairwise population comparisons: 
selectionanalyses(do_meta=FALSE,do_pairwise=TRUE,export="pdf",do_fsthet=TRUE,export_data=TRUE) 
#To run selection analyses for a pooled comparison between a predefined population division: 
selectionanalyses(do_meta=FALSE,do_pheno=TRUE,onlypooled=TRUE,export=”pdf”,phenolabels=c(“high”,”low”), do_thin=FALSE) 

## (3) ## 
# Sliding window Tajima’s D - Page 46:This function will export (in pdf format) graphs with sliding window Tajima’s D estimates. It will also produce a dataframe called wintaj. 
wintajd(my_chrom=25,winsize=1000000,winstep=200000)


## In R (5) ## page 48 from SambaR manual 

#Create an input file with geographical coordinates, sample names (called ‘individual’), and environmental variables, and store this txt file (for example under the name “evofile$

runRDA(export="pdf",legendpos="topleft",envfile=”evofile.txt",doall=FALSE)

backupdata("mySNPdata")


############################################################################### SambaR Additional Run ##########################################################################################

## Additional Analysis:

## (1) ##
#Linkage Disequilibrium - Page 39 - From the inputfiles directory (if you executed the exportsambarfiles() function), PED and MAP files ending on ‘filter2.number.map’ and ‘filter2.number.ped’- (Note: for this analysis you want to use ‘filter2’, not ‘filter’. 

plink --noweb --cow --allow-extra-chr --file yourpop.filter2.number --r2 --ld-window-kb 250 --ld-window-r2 0 –out yourpop

#Then in R to plot:
To plot these estimates, execute on the R command line:
LD_plot(export=”pdf”, xrange=c(0,1000000),stepsize=100000,addstripchart=FALSE, plotmeans=TRUE)
LDperchrom()

## (2) ##
# Plot PLINK relatedness (kinship) estimates - Page 41:

#SambaR contains a function to plot sample relatedness measures (pi_hat scores) generated by PLINK. For that purpose, you will find, if you have run the exportsambarfiles() function, within the inputfiles directory files called ‘metapop.filter.number.ped’, ‘metapop.filter.number.map’, ‘metapop.allinds.filter.number.ped’ and ‘metapop.allinds.filter.number.map’. 

#Navigate on the command line of your computer to this directory and execute these plink commands:

plink --file metapop.allinds.filter.number –-allow-extra-chr –genome -out plink.kin.allinds 
plink --file metapop.retainedinds.filter.number –-allow-extra-chr –genome -out plink.kin.retainedinds

#This will output to the SambaR input files directory a file called ‘plink.kin.allinds.genome’ and ‘plink.kin.retainedinds.genome’. In R, execute:

plinkrelatedness(export=”pdf”,per_pop=FALSE)

#Output files - This function will export several plots to the Kinship subdirectory: 
# - Relatedness.plink.between.pdf 
# - Relatedness.plink.within.pdf
# - Relatedness.plink.persample.pdf
# - Relatedness.plink.matrix.pdf

## In R (3) ## Page 46
# Sliding window Tajima’s D - Page 46:This function will export (in pdf format) graphs with sliding window Tajima’s D estimates. It will also produce a dataframe called wintaj. 
wintajd(my_chrom=25,winsize=1000000,winstep=200000)

## (4) ##  Page - 47
# Find and plot genes close to outlier SNPs:
# If you used reference mapping rather than denovo mapping, Sambar will have exported to the selection directory txt.files in BED format. These files contain lists with outlier SNPs, and can be used to find nearby genes using the software BEDTOOLS2 . To do so, you need the have the annotation file (in gff format) of the genome which you used as reference.
# a) On Unix command line execute the following command to convert this annotation file from gff format to BED format:

cut -f1,4,5,9 reindeer_coding_gene_annotations.gff > genes.bed

# b) Sort both BED files:

bedtools2/bin/sortBed -i reindeergenes.bed > genes.sorted.bed
bedtools2/bin/sortBed -i outliers.bed > outliers.sorted.bed

#To find the 10 closest features (including overlapping features), and to report the distance between SNP and feature (-D a flag), execute:
bedtools2/bin/closestBed -a outliers.sorted.bed -b genes.sorted.bed -D a -k 10 > putative.txt

#To subsequently select features which are within 200kb distance from an outlier SNP: 
awk '$9 <= 200000' putative.txt | cut -f1-8 | uniq -f7 | cut -f1,6,7,8 > putative. 200kb.bed 

# It might be that the gene names (4th column of BED file) are uninformative, and that you have to perform a blast with the gene sequence to find the actual gene name. 

# To plot the positions of the outlier SNPs and the adjacent genes, copypaste the BED-file to the selection subdirectory and execute:

multiplotscaffold(my_bed=”putativegenes.within200kb.genenames.bed",background_pop=NULL,doexport=TRUE,x_range=NULL,y_loc=c(0.4,0.65))

# This will export pdf files of all contigs/chromosomes containing outlier SNPs to the selection subdirectory. With the background_pop flag you can determine which population is shown on the background. With the x_range flag you can determine the region of the contig/chromosome you want to display. With the y_loc flag you can edit the locations of the gene names.  


## In R (5) ## page 48

#Create an input file with geographical coordinates, sample names (called ‘individual’), and environmental variables, and store this txt file (for example under the name “evofile.txt”) in the directory which also contains the PED and MAP input files. Next execute:

runRDA(export="pdf",legendpos="topleft",envfile=”evofile.txt",doall=FALSE)

## (6) ## Page 54
# Circosplot with Bayesass migration rates:

# SambaR contains functions to create input files for Bayesass3-SNPs, and to subsequently create a circos plot of the migration rates calculated by the software Bayesass3-SNPs. Execute:  
exportsambarfiles()
#Afterwards, you wil find in the inputfiles directory a file called ‘Bayesassinput.immanc.txt’. 

BA3-SNPS-Ubuntu64 --file Bayesassinput.immanc.txt --loci 15000 -u -t -s 10 -i 1000000 -b 100000
#This will run Bayesass with the defaults settings of 1 million iterations, a burn-in of 100000 generations, a seed of 10, and with delta (A, F and M) values of 0.10000000000000001, and for 15000 loci. Set the number of loci equal to:
nrow(snps[snps$filter,])

## See page 54 for run details

#To visualize migration rates in a circos plot execute the following SambaR function:
plotmigration("bayesassmatrix.txt",export=”pdf”,addlabels=TRUE) 
#The circosplot will be exported to the Divergence subdirectory.

## (7) ##: TreeMix plots - Page 55

#To generate the required input format for the software TreeMix (developed by the Pritchard Lab), run the command:
exporttreemix(snpsfilter=snps$filter,exportname="Treemixinput.snpsfilter.txt")

#In a Unix environment, run the software Treemix using the following commands (assuming the outgroup population is called ‘OutPop’:
gzip treemixinput.snpsfilter.txt
treemix -i Treemixinput.snpsfilter.txt.gz -m 0 -root OutPop -o Treemixout.0
treemix -i Treemixinput.snpsfilter.txt.gz -m 1 -root OutPop -o Treemixout.1
treemix -i Treemixinput.snpsfilter.txt.gz -m 2 -root OutPop -o Treemixout.2

#Afterwards transfer all output files to your SambaR working directory, and run in R the commands:
plottreemix(prefix="Treemixout.0",export="pdf",myxmin=0,plotname="Treemixoutput.0")
plottreemix(prefix="Treemixout.1",export="pdf",myxmin=0,plotname="Treemixoutput.1")
plottreemix(prefix="Treemixout.2",export="pdf",myxmin=0,plotname="Treemixoutput.2")

#Output plots will be exported to the working directory.

