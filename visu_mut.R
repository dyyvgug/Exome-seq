#=======================================================================
# 2021-9-10.Author:Dong Yingying.Visually describe mutation.
#=======================================================================
if (!require("maftools"))
  BiocManager::install("maftools")

library(maftools)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
ann_array = list.files(getwd())
ann_array

var.annovar.maf <- annovarToMaf(annovar = "merge.txt",
                                refBuild = "hg38",
                                tsbCol = "Tumor_Sample_Barcode",
                                sep = "\t")
write.table(var.annovar.maf,file="var_annovar.maf",quote= F,sep="\t",row.names=F)

var_maf = read.maf(maf ="var_annovar.maf")
sam_sum = getSampleSummary(var_maf)
write.table(sam_sum,file="sample_summary.maf",quote= F,sep="\t")
getGeneSummary(var_maf)

pdf("maf_sum.pdf",width = 7,height = 4)
plotmafSummary(maf = var_maf, rmOutlier = TRUE, addStat = 'median')
dev.off()

oncoplot(maf = var_maf, top = 30, fontSize = 12 ,showTumorSampleBarcodes = F )
laml.titv = titv(maf = var_maf, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
pdf("mut_int.pdf",width = 7,height = 4)
somaticInteractions(maf = var_maf, top = 25, pvalue = c(0.05, 0.1))
dev.off()
