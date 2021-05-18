###################################
# Imputation Pipeline
###################################

## Requisites:

#~~~~~~ Programs ~~~~~~

# * R 4.0.3
# 
# * plink 1.9
# 
# * vcftools-0.1.16
# 
# * python-2.7.18
# 
# * bcftools-1.12
# 
# * perl-5.26.2
# 
# * perl-vcftools-vcf-0.1.16
# 
# * p7zip-15.09

#~~~~~~ Data ~~~~~~

# * Data must be in plink format (bed + bim + fam) in the data folder.

#~~~~~~ Oncoarray Ancestry Informative Marker SNPs ~~~~~~

# * AIMS_SNPs.txt file must be in the data folder.

#~~~~~~ liftover ~~~~~~

# <https://github.com/sritchie73/liftOverPlink/blob/master/liftOverPlink.py>

# * liftover/liftOverPlink.py 

# * liftover/liftOver

# * liftover/hg18ToHg38.over.chain.gz if assemby is hg18


#~~~~~~ bim check ~~~~~~

# <https://www.well.ox.ac.uk/~wrayner/tools/>

# * bimcheck/HRC-1000G-check-bim.pl

# * bimcheck/CreateTOPMed.pl

# * bimcheck/PASS.Variants.TOPMed_freeze5_hg38_chr20_dbSNP.tab

# To obtain PASS.Variants.TOPMed_freeze5_hg38_chr20_dbSNP.tab file:

# 1.- Download it

# system("curl 'https://bravo.sph.umich.edu/freeze5/hg38/download/all' -H 'Accept-Encoding: gzip, deflate, br' -H --compressed > ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz")

# 2.- Once downloaded the VCF can be converted to an HRC formatted reference legend using the code here: CreateTOPMed.zip

# system("perl ./CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz")

# !!! 10 hours of execution!!!

# 3.- FOR THIS EXAMPLE WE EXTRACT CHR20

# system("zcat /mnt/typhon/data/references/SNPs/TOPMed/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz | awk '$1 == 20' | gzip -c > /mnt/hydra/ubs/shared/users/Anna/GENRISK/bimcheck/PASS.Variants.TOPMed_freeze5_hg38_chr20_dbSNP.tab.gz")

#############################################################

# set working directory
setwd("./Practical")

######################################################################################
# QUALITY CONTROL
######################################################################################

# * dataDir: data directory with the genotyped plink files
dataDir<-"./data/"

# * plinkFile: name of the plink files with the genotyped data (plinkFile.bed; plinkFile.bim and plinkFile.fam)
plinkFile<-"ExampleData"

rmarkdown::render("QC.Rmd", params = list(
  dataDir = dataDir,
  plinkFile = plinkFile
))

outputDir<-"./QC/"

# Choose threshold of ancestry PCA
pca <- read.table(paste0(outputDir, "pca.eigenvec"), header = FALSE)
outpop <- pca[pca[,3]<(-0.1),]

# plot
plot(pca[,3], pca[,4], pch = 20, xlab = "PCA 1", ylab = "PCA 2", main = plinkFile)
points(outpop[,3], outpop[,4], col="red", pch = 20)

# save plot
png(paste0(outputDir, "Ancestry.png"), res = 200, 1500, 1500)
plot(pca[,3], pca[,4], pch = 20, xlab = "PCA 1", ylab = "PCA 2", main = plinkFile)
points(outpop[,3], outpop[,4], col="red", pch = 20)
dev.off()

# Mark samples
samples <- read.table(paste0(outputDir, "marked_samples.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(samples) <- samples[,2]
samples$ExcludePCA <- "No"
samples[as.character(outpop[,2]),"ExcludePCA"] <- "Yes"

# Save file
write.table(samples, file = paste0(outputDir, "marked_samples.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n")
cat(paste0("There are ",sum(samples$ExcludePCA=="Yes")," samples outside the ancestry cluster of the data."))
cat("\n")


if(sum(samples$ExcludePCA=="Yes")>0)  print(samples[samples$ExcludePCA=="Yes",c(1,2,11)])


## In this study, the marked samples are:
marked_samples<-samples[apply(samples=="Yes",1,sum)>0,c(2,7:11)]
marked_samples


###############################################################################################
# REMOVE SAMPLES IF NEEDED
###############################################################################################

# read marked samples file

# select samples to filter out
exclude <- samples[samples$ExcludeMissings == "Yes" | samples$ExcludeSex == "Yes" | samples$ExcludeHeteroz == "Yes"  | samples$ExcludeDups != "No",]

# choose wich related samples to mantain
exclude<-exclude[-5,]

# save the first two columns
write.table(exclude[,1:2], file = paste0(outputDir, "samples2exclude.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


# filter out those samples using plink
system(paste0("plink --bfile ", outputDir, plinkFile, "_QC --remove ", outputDir,
              "samples2exclude.txt --make-bed --out ", outputDir, plinkFile, "_QC_sampfilt"))

# remove some files
unlink(paste0(outputDir, plinkFile, "_QC.*"))
unlink(paste0(outputDir, plinkFile, "_QC_sampfilt.hh"))
unlink(paste0(outputDir, plinkFile, "_QC_sampfilt.log"))



###############################################################################################
# FILTER CHR20 FOR THIS EXAMPLE
###############################################################################################

system(paste0("plink --bfile ", outputDir, plinkFile, "_QC_sampfilt --chr 20 --make-bed --out ", outputDir, plinkFile, "_chr20"))
unlink(paste0(outputDir, plinkFile, "_chr20.log"))


################################################################################################
## PREPARE DATA FOR IMPUTATION
################################################################################################

# * dataDir: data directory with the genotyped and quality controled plink files
dataDir<-outputDir

# * plinkFile: name of the plink files with the genotyped and quality controled data (plinkFile.bed; plinkFile.bim and plinkFile.fam)
plinkFile <- paste0(plinkFile,"_chr20")

# * assembly: genome version of your data. Accepted values are hg18 or hg19.
assembly <- "hg18"

rmarkdown::render("Data_preparation.Rmd", params = list(
  dataDir = dataDir,
  plinkFile = plinkFile,
  assembly = assembly
))


################################################################################################
## IMPUTATION
################################################################################################

# a) Upload VCF files to TOPMed Imputation Server
# 
# <https://imputation.biodatacatalyst.nhlbi.nih.gov/>
#   
# * Minimac 4
# * Reference Panel: TOPMed r2
# * Built: GRCh38/hg38
# * rsq Filter: 0.3
# * Phasing: Eagle v2.4 (unphased input)
# * QC Frequency Check: vs.TOPMed Panel
# * Mode: Quality Control and Imputation
# * check AES 256 encryption (info file inside vcf)
# 
# b) Download results from imputation to ImputationResults folder
# 
# c) Check logs and qcreport.html
# 
# d) Unzip data


dir.create("./ImputationResults/",showWarnings=FALSE)
setwd("./ImputationResults/")

password<-"%bY6A5nOTCtopV"

system(paste0("7z x -p",password," chr_20.zip"))



##########################################################################################################################
##########################################################################################################################
# check info file

library(data.table)
infofile<-fread("chr20.info.gz",sep="\t",stringsAsFactors=FALSE,header=TRUE)
dim(infofile)
# 391605      13

head(infofile)

summary(infofile$Rsq)

table(infofile$Genotyped)
# Genotyped   Imputed 
#     33813    357792


######################################################################################################
######################################################################################################
# extract some random SNPs

idx<-sample(1:nrow(infofile), size=20, replace=FALSE)

# save snps
write.table(infofile[idx,1], file = "snps.txt", sep = "\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

# extraction
system("plink --vcf chr20.dose.vcf.gz --extract snps.txt  --keep-allele-order --make-bed --out chr20_filt")


