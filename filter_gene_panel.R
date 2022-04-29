#=======================================================================
# 2021-10.Author:Dong Yingying.Filtering VCF based on gene panel list.
#=======================================================================

library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ann_array = list.files(getwd())
ann_array
name1="trio_snp_indel_VQSR"
name2="trio_only_pass"


gene<-read.table("kidney_gene_panel_only_name.txt",sep = "\t", header = F,quote = "")
vq<-read.table(paste0(name1,".hg38_multianno.txt"),sep = '\t',header = T,quote = "")
pass<-read.table(paste0(name2,".hg38_multianno.txt"),sep = "\t", header = T,quote = "")
# keep only exonic, filter intronic,UTR3
vq_ex<-vq[grep("exonic",vq$Func.refGene),]
pass_ex<-pass[grep("exonic",pass$Func.refGene),]

names(gene)="Gene.refGene"
vq_mer=semi_join(vq,gene,by="Gene.refGene")
pass_mer=semi_join(pass,gene,by="Gene.refGene")
vq_ex_mer=semi_join(vq_ex,gene,by="Gene.refGene")
pass_ex_mer=semi_join(pass_ex,gene,by="Gene.refGene")


pass<-pass[,-c(6:21)]
pass_ex<-pass_ex[,-c(6:21)]
pass_ex_mer<-pass_ex_mer[,-c(6:21)]
vq<-vq[,-c(6:21)]
vq_mer<-vq_mer[,-c(6:21)]
vq_ex_mer<-vq_ex_mer[,-c(6:21)]


#write.table(pass,file = paste0(name1,"_rem.vcf"),quote= F,sep="\t",row.names=F,col.names = F)
write.table(pass_ex_mer,file = "pass_ex_pan_rem.vcf",quote= F,sep="\t",row.names=F,col.names = F)
write.table(vq_ex_mer,file = "vq_ex_pan_rem.vcf",quote= F,sep="\t",row.names=F,col.names = F)

df_vq<-data.frame(vq_ex_mer$Start,vq_ex_mer$Ref,vq_ex_mer$Alt,vq_ex_mer$Gene.refGene)
df_pass<-data.frame(pass_ex_mer$Start,pass_ex_mer$Ref,pass_ex_mer$Alt,pass_ex_mer$Gene.refGene)
write.table(df_pass,file = "only_pass_for_Varcopp.txt",quote= F,sep="\t",row.names=F,col.names = T)
write.table(df_vq,file = "vqsr_for_varcopp.txt",quote= F,sep="\t",row.names=F,col.names = T)

var_vq<-read.table("vqsr_aadiff.txt",sep = '\t',header = T)
df<-cbind(vq_ex_mer,var_vq$Flex_diff)
