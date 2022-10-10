#=======================================================================
# 2021-9-10.Modify:2021-10-8.Author:Dong Yingying.Visually describe mutation.
#=======================================================================
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("maftools"))
  BiocManager::install("maftools")

library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ann_array = list.files(getwd())
ann_array
name="trio_only_pass"
#name="trio_only_pass"
type1="_fre" # MAF>0.05
type2="_rare" # MAF<0.01
type3="_lowFRE" # 0.01<MAF<0.05

vcf_list=read.table(paste0("./anno/",name,".hg38_multianno.txt"),sep = "\t", header = T,quote = "")
fre_list=read.table(paste0(name,type1,".hg38_gnomad_exome_dropped"),sep = "\t",header = F,quote = "")
lowFRE_list=read.table(paste0(name,type1,".hg38_gnomad_exome_filtered"),sep = "\t",header = F,quote = "")

rare_list=read.table(paste0(name,type2,".hg38_gnomad_exome_dropped"),sep = "\t",header = F,quote = "")
lowRARE_list=read.table(paste0(name,type2,".hg38_gnomad_exome_filtered"),sep = "\t",header = F,quote = "")

fre_list<-fre_list[,-c(1,8:18)]
lowFRE_list<-lowFRE_list[,-c(6:16)]
rare_list<-rare_list[,-c(1,8:18)]
lowRARE_list<-lowRARE_list[,-c(6:16)]

names(lowFRE_list)=c("Chr","Start","End","Ref","Alt")
names(lowRARE_list)=c("Chr","Start","End","Ref","Alt")
names(fre_list)=c("MAF","Chr","Start","End","Ref","Alt")
names(rare_list)=c("MAF","Chr","Start","End","Ref","Alt")

fre_mer<-semi_join(vcf_list,fre_list,by="Start")
rare_mer<-semi_join(vcf_list,rare_list,by="Start")
write.table(fre_mer,file = paste0(name,type1,"_list.txt"),quote = F,sep = '\t',row.names = F)
write.table(rare_mer,file = paste0(name,type2,"_list.txt"),quote = F,sep = '\t',row.names = F)
low_fre<-semi_join(lowFRE_list,lowRARE_list,by="Start")
low_mer<-semi_join(vcf_list,low_fre,by="Start")
write.table(low_mer,file = paste0(name,type3,"_list.txt"),quote = F,sep = '\t',row.names = F)

all_maf_file <- annovarToMaf(annovar =paste0("./anno/",name,".hg38_multianno.txt"),
                             refBuild = "hg38",sep = "\t")
write.table(all_maf_file,file=paste0(name,"_all_list.maf"),quote= F,sep="\t",row.names=F)

txt_array<-c(type1,type2,type3)
for (i in txt_array) {
  
  maf_file<-annovarToMaf(annovar =paste0(name,i,"_list.txt"),
                         refBuild = "hg38",sep = "\t")
  write.table(maf_file,file=paste0(name,i,"_list.maf"),quote= F,sep="\t",row.names=F)
}

maf_array<-c("_all",type1,type2,type3)
for (j in maf_array) {
  
  var_maf = read.maf(maf =paste0(name,j,"_list.maf"))
  maf_sum = getSampleSummary(var_maf)
  write.table(maf_sum,file=paste0(name,j,"_summary.maf"),quote= F,sep="\t")
  getGeneSummary(var_maf)
  
  pdf(paste0(name,j,"_sum.pdf"),width = 9,height = 5.7)
  plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
  dev.off()
}

library(dplyr)

# ======================================================
#  Filtering VCF based on gene panel list
# =====================================================
name1="trio_snp_indel_VQSR"
name2="trio_only_pass"

gene<-read.table("kidney_gene_panel_only_name.txt",sep = "\t", header = F,quote = "")
#vq<-read.table(paste0("../anno/",name1,".hg38_multianno.txt"),sep = '\t',header = T,quote = "")
pass<-read.table(paste0("./anno/",name2,".hg38_multianno.txt"),sep = "\t", header = T,quote = "")
# keep only exonic, filter intronic,UTR3
pass_ex<-pass[grep("exonic",pass$Func.refGene),]
pass_af<-read.table("trio_only_pass_rare_list.txt",sep = "\t", header = T,quote = "")
pass_af_ex<-pass_af[grep("exonic",pass_af$Func.refGene),]

names(gene)="Gene.refGene"
pass_mer=semi_join(pass,gene,by="Gene.refGene")
pass_ex_mer=semi_join(pass_ex,gene,by="Gene.refGene")
pass_af_mer=semi_join(pass_af_ex,gene,by="Gene.refGene")

pass<-pass[,-c(6:21)]
pass_ex<-pass_ex[,-c(6:21)]
pass_ex_mer<-pass_ex_mer[,-c(6:21)]


#write.table(pass,file = paste0(name1,"_rem.vcf"),quote= F,sep="\t",row.names=F,col.names = F)
write.table(pass_ex_mer,file = "pass_ex_pan_rem.vcf",quote= F,sep="\t",row.names=F,col.names = T)
write.table(vq_ex_mer,file = "vq_ex_pan_rem.vcf",quote= F,sep="\t",row.names=F,col.names = T)

write.table(pass_af_mer,file = "pass_rare_ex_pan_rem.vcf",quote= F,sep="\t",row.names=F,col.names = T)
write.table(vq_af_mer,file = "vq_rare_ex_pan_rem.vcf",quote= F,sep="\t",row.names=F,col.names = T)

# only gene panel,rare, only pass, maf
var.annovar.maf <- annovarToMaf(annovar = "pass_rare_ex_pan_rem.vcf",
                                refBuild = "hg38",
                                sep = "\t")
write.table(var.annovar.maf,file="pass_rare_ex_pan_rem.maf",quote= F,sep="\t",row.names=F)

var_maf = read.maf(maf ="pass_rare_ex_pan_rem.maf")
sam_sum = getSampleSummary(var_maf)
write.table(sam_sum,file="pass_rare_ex_pan_rem.maf",quote= F,sep="\t")
getGeneSummary(var_maf)

pdf("pass_rare_ex_pan_rem.pdf",width = 8,height = 5.7)
plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
dev.off()
