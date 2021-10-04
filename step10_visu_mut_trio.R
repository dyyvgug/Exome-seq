#=======================================================================
# 2021-9-10.Author:Dong Yingying.Visually describe mutation.
#=======================================================================
if (!require("maftools"))
  BiocManager::install("maftools")
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ann_array = list.files(getwd())
ann_array
vcf_list=read.table("trio_snp_indel_VQSR.hg38_multianno.txt",sep = "\t", header = T,quote = "")
fre_list=read.table("trio_snp_indel_VQSR.hg38_gnomad_exome_filtered",sep = "\t",header = F,quote = "")
rare_list=read.table("trio_snp_indel_VQSR.hg38_gnomad_exome_dropped",sep = "\t",header = F,quote = "")
fre_list<-fre_list[,-c(6:16)]
rare_list<-rare_list[,-c(1,8:18)]
names(fre_list)=c("Chr","Start","End","Ref","Alt")
names(rare_list)=c("MAF","Chr","Start","End","Ref","Alt")

fre_mer<-semi_join(vcf_list,fre_list,by="Start")
rare_mer<-semi_join(vcf_list,rare_list,by="Start")
write.table(fre_mer,file = "fre_list.txt",quote = F,sep = '\t',row.names = F)
write.table(rare_mer,file = "rare_list.txt",quote = F,sep = '\t',row.names = F)

all_maf_file <- annovarToMaf(annovar ="trio_snp_indel_VQSR.hg38_multianno.txt",
                                refBuild = "hg38",sep = "\t")
write.table(all_maf_file,file="all_list.maf",quote= F,sep="\t",row.names=F)

fre_maf_file <- annovarToMaf(annovar ="fre_list.txt",
                              refBuild = "hg38",sep = "\t")
write.table(fre_maf_file,file="fre_list.maf",quote= F,sep="\t",row.names=F)

rare_maf_file <- annovarToMaf(annovar ="rare_list.txt",
                                refBuild = "hg38",sep = "\t")
write.table(rare_maf_file,file="rare_list.maf",quote= F,sep="\t",row.names=F)

maf_array<-c("all","fre","rare")
for (i in maf_array) {
  var_maf = read.maf(maf =paste0(i,"_list.maf"))
  maf_sum = getSampleSummary(var_maf)
  write.table(maf_sum,file=paste0(i,"_summary.maf"),quote= F,sep="\t")
  getGeneSummary(var_maf)
  
  pdf(paste0(i,"_sum.pdf"),width = 9,height = 5.7)
  plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
  dev.off()
  }



