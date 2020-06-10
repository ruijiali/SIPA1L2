setwd("D:/tan/ad_bcf_plink/skat")
library(SKAT)

#setID
library(stringr)
df <- read.delim("Michi_filter_final_dup.avinput.variant_function", header=FALSE)
newdf <- df[df$V1 == "exonic"|df$V1 == "splicing",]

dim <- rownames(newdf)
bim <- read.delim("Michi_filter_final_dup.bim", header=FALSE)
bim2 <- bim[dim, ]
bim2$V1 <- newdf$V2

write.table(bim2[,c(1,2)],file="Michi_filter_final_dup.SetID",quote=FALSE,row.names=FALSE,col.names=FALSE)

#generate SSD file
File.Bed <- "./Michi_filter_final_dup.bed"
File.Bim <- "./Michi_filter_final_dup.bim"
File.Fam <- "./Michi_filter_final_dup.fam"
File.Cov <- "./Michi_filter_final_dup.cov"
File.SetID <- "./Michi_filter_final_dup.SetID"
File.Info <- "./Michi_filter_final_dup.info"
File.SSD <- "./Michi_filter_final_dup.SSD"

Generate_SSD_SetID(File.Bed,File.Bim, File.Fam, File.SetID,File.SSD,File.Info)

FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=TRUE, cov_header=FALSE)

SSD.INFO <- Open_SSD(File.SSD,File.Info)

# adjustment is needed
obj<-SKAT_Null_Model(Phenotype ~ COV1 + COV2 + COV3, out_type="D", data=FAM,n.Resampling=1000
                     , type.Resampling="bootstrap")

#obj<-SKAT_Null_Model(Phenotype ~ COV1 + COV2 + COV3 + COV4 , out_type="D", data=FAM,n.Resampling=1000
, type.Resampling="bootstrap")

out.skat<-SKATBinary.SSD.All(SSD.INFO, obj, method="SKAT",method.bin="Hybrid", weights.beta=c(1,25))


#Resampling_FWER_1(out.skat$results$P.value, out.skat$P.value.Resampling, FWER=0.05)
Resampling_FWER(out.skat,FWER=0.05)

QQPlot_Adj(out.skat$results$P.value, out.skat$results$MAP,Is.unadjsted=TRUE)


#EffectiveNumberTest

Get_EffectiveNumberTest(out.skat$results$MAP, alpha=0.05, Is.MidP=TRUE)

#update SKAT
install.packages('githubinstall')
library(githubinstall)
install.packages("devtools")
library(devtools)
install_github("leeshawn/SKAT")