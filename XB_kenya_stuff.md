# Samples

I am trying to nail down what samples we have mtDNA data for, and what samples we used for NGS.  

BenF put the RAdseq data for the laevis genome and Austin genome here:
`/2/scratch/evanslab/2019_RADseq_KenyaXBXL_GhanaEastfamily/plate1/genotyped/`

and these are the file names:
`mpileup_raw_wildBorealis_AustinGenome.vcf.gz` and `mpileup_raw_wildBorealis_allChrs.vcf.gz`

BenF says they should be a basic bam alignment and samtools genotypic scheme, probably throwing out anything with a map quality < 20, but no other filters. 


# Iqtree

On evanslab: `/2/scratch/evanslab/2019_RADseq_KenyaXBXL_GhanaEastfamily/mtDNA/`
```
iqtree -s Kenya_borealis_only16s.nex -m TEST -nt 1 -pre Kenya_borealis_only16s.nex_
```
```
iqtree -s Kenya_borealis_only16s.nex -m TPM2u+I -bb 1000
```
with mom and dad de novo:
```
iqtree -s Kenya_borealis_only16s_withMomDad.nex -m TEST -nt 1 -pre Kenya_borealis_only16s_withMomDad.nex_
```
```
iqtree -s Kenya_borealis_only16s_withMomDad.nex -m HKY+I -bb 1000
```

# PCA
I'm interested in comparing male and female-only PCAs in the SL and nonSL regions of chr8L. I'm working here (and locally for pca):
```
/2/scratch/evanslab/2019_RADseq_KenyaXBXL_GhanaEastfamily/plate1/genotyped
```

First I extracted only male and only female seqs from the SL and nonSL regions of chr8L:
```
bcftools view mpileup_raw_justBorealis_allChrs.vcf.gz --regions chr8L:1-57000000 --samples Fem_Lukhome_BJE4441_Xb.fq.gz,Fem_Lukhome_BJE4444_Xb.fq.gz,Fem_Lukhome_BJE4445_Xb.fq.gz,Fem_Lukhome_BJE4446_Xb.fq.gz,Fem_Kiminini_BJE4429_Xb.fq.gz,Fem_Kiminini_BJE4433_Xb.fq.gz,Fem_Chesuwe_BJE4479_Xb.fq.gz,Fem_Cheplaskei_BJE4459_Xb.fq.gz,Fem_Cheplaskei_BJE4461_Xb.fq.gz,Fem_Cheplaskei_BJE4470_Xb.fq.gz,Fem_Eldoret_BJE4471_Xb.fq.gz,Fem_Eldoret_BJE4472_Xb.fq.gz,Fem_Eldoret_BJE4474_Xb.fq.gz,Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Fem_NMK_BJE4562_Xb.fq.gz,Fem_NMK_BJE4564_Xb.fq.gz,Fem_NMK_BJE4567_Xb.fq.gz,Fem_NMK_BJE4568_Xb.fq.gz,Fem_Wundanyi_BJE4515_Xb.fq.gz,Fem_Wundanyi_BJE4516_Xb.fq.gz,Fem_Wundanyi_BJE4534_Xb.fq.gz,Fem_Wundanyi_BJE4535_Xb.fq.gz,Fem_Wundanyi_BJE4541_Xb.fq.gz -O z -o mpileup_raw_justBorealis_allChrs_XB_females_chr8L_SL_only.vcf.gz


bcftools view mpileup_raw_justBorealis_allChrs.vcf.gz --regions chr8L:1-57000000 --samples Mal_Lukhome_BJE4442_Xb.fq.gz,Mal_Lukhome_BJE4443_Xb.fq.gz,Mal_Lukhome_BJE4447_Xb.fq.gz,Mal_Kiminini_BJE4430_Xb.fq.gz,Mal_Kiminini_BJE4431_Xb.fq.gz,Mal_Kiminini_BJE4432_Xb.fq.gz,Mal_Chesuwe_BJE4477_Xb.fq.gz,Mal_Chesuwe_BJE4478_Xb.fq.gz,Mal_Chesuwe_BJE4480_Xb.fq.gz,Mal_Kisumu_BJE4391_Xb.fq.gz,Mal_Cheplaskei_BJE4460_Xb.fq.gz,Mal_Cheplaskei_BJE4462_Xb.fq.gz,Mal_Cheplaskei_BJE4465_Xb.fq.gz,Mal_Cheplaskei_BJE4469_Xb.fq.gz,Mal_Eldoret_BJE4473_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz,Mal_NMK_BJE4563_Xb.fq.gz,Mal_Wundanyi_BJE4536_Xb.fq.gz,Mal_Wundanyi_BJE4537_Xb.fq.gz,Mal_Wundanyi_BJE4538_Xb.fq.gz,Mal_Wundanyi_BJE4539_Xb.fq.gz,Mal_Wundanyi_BJE4540_Xb.fq.gz -O z -o mpileup_raw_justBorealis_allChrs_XB_males_chr8L_SL_only.vcf.gz



bcftools view mpileup_raw_justBorealis_allChrs.vcf.gz --regions chr8L:57000001-120101378 --samples Fem_Lukhome_BJE4441_Xb.fq.gz,Fem_Lukhome_BJE4444_Xb.fq.gz,Fem_Lukhome_BJE4445_Xb.fq.gz,Fem_Lukhome_BJE4446_Xb.fq.gz,Fem_Kiminini_BJE4429_Xb.fq.gz,Fem_Kiminini_BJE4433_Xb.fq.gz,Fem_Chesuwe_BJE4479_Xb.fq.gz,Fem_Cheplaskei_BJE4459_Xb.fq.gz,Fem_Cheplaskei_BJE4461_Xb.fq.gz,Fem_Cheplaskei_BJE4470_Xb.fq.gz,Fem_Eldoret_BJE4471_Xb.fq.gz,Fem_Eldoret_BJE4472_Xb.fq.gz,Fem_Eldoret_BJE4474_Xb.fq.gz,Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Fem_NMK_BJE4562_Xb.fq.gz,Fem_NMK_BJE4564_Xb.fq.gz,Fem_NMK_BJE4567_Xb.fq.gz,Fem_NMK_BJE4568_Xb.fq.gz,Fem_Wundanyi_BJE4515_Xb.fq.gz,Fem_Wundanyi_BJE4516_Xb.fq.gz,Fem_Wundanyi_BJE4534_Xb.fq.gz,Fem_Wundanyi_BJE4535_Xb.fq.gz,Fem_Wundanyi_BJE4541_Xb.fq.gz -O z -o mpileup_raw_justBorealis_allChrs_XB_females_chr8L_nonSL_only.vcf.gz


bcftools view mpileup_raw_justBorealis_allChrs.vcf.gz --regions chr8L:57000001-120101378 --samples Mal_Lukhome_BJE4442_Xb.fq.gz,Mal_Lukhome_BJE4443_Xb.fq.gz,Mal_Lukhome_BJE4447_Xb.fq.gz,Mal_Kiminini_BJE4430_Xb.fq.gz,Mal_Kiminini_BJE4431_Xb.fq.gz,Mal_Kiminini_BJE4432_Xb.fq.gz,Mal_Chesuwe_BJE4477_Xb.fq.gz,Mal_Chesuwe_BJE4478_Xb.fq.gz,Mal_Chesuwe_BJE4480_Xb.fq.gz,Mal_Kisumu_BJE4391_Xb.fq.gz,Mal_Cheplaskei_BJE4460_Xb.fq.gz,Mal_Cheplaskei_BJE4462_Xb.fq.gz,Mal_Cheplaskei_BJE4465_Xb.fq.gz,Mal_Cheplaskei_BJE4469_Xb.fq.gz,Mal_Eldoret_BJE4473_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz,Mal_NMK_BJE4563_Xb.fq.gz,Mal_Wundanyi_BJE4536_Xb.fq.gz,Mal_Wundanyi_BJE4537_Xb.fq.gz,Mal_Wundanyi_BJE4538_Xb.fq.gz,Mal_Wundanyi_BJE4539_Xb.fq.gz,Mal_Wundanyi_BJE4540_Xb.fq.gz -O z -o mpileup_raw_justBorealis_allChrs_XB_males_chr8L_nonSL_only.vcf.gz
```
Then I made the PCAs using an R script with no LD pruning because I want to capture linked variants on the W:
```R
#https://github.com/zhengxwen/SNPRelate/issues/13
#http://corearray.sourceforge.net/tutorials/SNPRelate/
#https://github.com/zhengxwen/SNPRelate/wiki/Preparing-Data

#GenotypeVCFs_noBSQR_filtered_aDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_xDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_aDNA_only_no_lowcoverage_individuals.vcf

library("devtools")
library(gdsfmt)
library(SNPRelate)
setwd("/Users/Shared/Previously Relocated Items/Security/projects/XB_sex_determining_gene/XB_chr8L_SL_within_sex_comparisons")
vcf.fn <- "mpileup_raw_justBorealis_allChrs_XB_females_chr8L_SL_only.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "test5.gds", method="copy.num.of.ref") #,ignore.chr.prefix = "chr")
genofile = snpgdsOpen("test5.gds", readonly=FALSE)
samp.annot<-data.frame(pop.group = c("west","west","west","west","west","west","west","west","west","west","west","west","west",
                                     "east","east","east","east","east","east","east","east","east","east","east","east","east"))
add.gdsn(genofile, "sample.annot", samp.annot)
snpgdsSummary("test5.gds")


# LD prinnung
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, missing.rate=0.35, verbose = TRUE)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

# or without LD prunning
#pca <- snpgdsPCA(genofile, num.thread=2)
#pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))


tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#text(tab$EV2, tab$EV1,labels=tab$sample.id, cex= 0.4)

library(ggplot2)
#ggplot(...)+...+ theme(axis.text.x = element_text(angle=60, hjust=1))
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)

pdf("PCA_males_chr8L_SLonly.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("west","west","west","west","west","west","west","west","west","west","west","west","west",
                 "Njoro","Njoro","Njoro","Njoro",
                 "Nairobi","Nairobi","Nairobi","Nairobi",
                 "Wundyani","Wundyani","Wundyani","Wundyani","Wundyani")
tab$samp.color <- c("red","red","red","red","red","red","red","red","red","red","red","red","red",
                    "steel blue","steel blue","steel blue","steel blue",
                    "blue","blue","blue","blue",
                    "dark blue","dark blue","dark blue","dark blue","dark blue")
tab$samp.fieldid <- c("BJE4441","BJE4444","BJE4445","BJE4446","BJE4429","BJE4433","BJE4479","BJE4459","BJE4461","BJE4470","BJE4471","BJE4472","BJE4474",
                      "BJE4363","BJE4364","BJE4367","BJE4368",
                      "BJE4562","BJE4564","BJE4567","BJE4568",
                      "BJE4515","BJE4516","BJE4534","BJE4535","BJE4541")

d<- ggplot(data=tab, aes(x=-EV1,y=EV2, label = samp.fieldid, color = samp.color)) +
  # label axis 
  labs(x=expression("-Eigenvector 1"), y=expression("Eigenvector 2")) +
  # legend details
  scale_colour_manual(name="", values = c("red","red","red","red","red","red","red","red","red","red","red","red","red",
                                                 "steel blue","steel blue","steel blue","steel blue",
                                                 "blue","blue","blue","blue",
                                                 "dark blue","dark blue","dark blue","dark blue","dark blue"),
                                      breaks=c("red","red","red","red","red","red","red","red","red","red","red","red","red",
                                               "steel blue","steel blue","steel blue","steel blue",
                                               "blue","blue","blue","blue",
                                               "dark blue","dark blue","dark blue","dark blue","dark blue"),
                                      labels=c("west","west","west","west","west","west","west","west","west","west","west","west","west",
                                               "Njoro","Njoro","Njoro","Njoro",
                                               "Nairobi","Nairobi","Nairobi","Nairobi",
                                               "Wundyani","Wundyani","Wundyani","Wundyani","Wundyani"))+
  # add points and fieldID labels
  geom_text_repel(aes(-EV1,EV2, label=(samp.fieldid)), size=3) + geom_point(size=3) + 
  # change to cleaner theme
  theme_classic(base_size = 16) +
  # make it clean
  theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  # italicize species names
  #theme(legend.text = element_text(face="italic"))+ 
  # make the text bigger
  theme(text = element_text(size=10)) +
  # move the legend
  #theme(legend.position = c(.2, .2)) +
  # add a title
  #ggtitle("PCA of Sulawesi Samples") + 
  # add space around axis
    theme(axis.text.x = element_text(margin=margin(10,10,10,10,"pt")),
        axis.text.y = element_text(margin=margin(10,10,10,10,"pt"))) +
  # remove boxes around legend symbols
  theme(legend.key = element_blank()) #+
#  annotate(geom = "text", x = 0.18, y = 0.25, label = "Nairobi", color = "blue", angle = 50, size=5)+
#  annotate(geom = "text", x = .25, y = -0.2, label = "Wundyani", color = "blue", angle = 0, size=5)+
#  annotate(geom = "text", x = 0.05, y = -0.18, label = "Njoro", color = "blue", angle = 0, size=5)+
#  annotate(geom = "text", x = -0.2, y = 0.25, label = "West Kenya", color = "red", angle = 0, size=5)
d
dev.off()


closefn.gds(genofile)

```
