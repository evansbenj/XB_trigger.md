# Analyzing tad expression data.

I wrote an R script to analyze Xue's expression data from XB tads:
```R
setwd("/Users/Shared/Relocated\ Items/Security/projects/XB_sex_determining_gene/RNAseq/reexpressiondata")
XB_data <- read.table("bor_tad_borGenome_chr8l_combine.csv", 
                      header = TRUE, sep = ",")
head(XB_data)
ls(XB_data)
hist(XB_data$start, breaks = seq(0,135000000,500000))
library(ggplot2)
g <- ggplot(XB_data, aes(x = start, y = (1-st48_FDR))) +
  geom_point() + geom_smooth(method = "loess")
g

XB_scaffold_data <- read.table("bor_tad_borGenome_Scaffolds_combine.csv", 
                      header = TRUE, sep = ",")



str(XB_data)
dim(XB_data)
rowMeans(XB_data[9:20])

# make it print numbers normally
options(scipen=999)

# plot the difference between the F and M mean expression values
g <- ggplot(XB_data, aes(x = start, y = (rowMeans(XB_data[c(9:12,17:24)]) -
                                           rowMeans(XB_data[c(13:16,25:29)])))) +
  geom_point() + geom_smooth(method = "loess")+ theme_classic()
g

g <- ggplot(XB_scaffold_data , aes(x = start, y = (rowMeans(XB_scaffold_data[c(9:12,17:24)]) 
                         - rowMeans(XB_scaffold_data[c(13:16,25:29)])))) +
  geom_point() + geom_smooth(method = "loess")+ theme_classic()
g

# plot the abs of the difference between the F and M mean expression values
g <- ggplot(XB_data, aes(x = start, y = (abs(rowMeans(XB_data[c(9:12,17:24)]) 
                                             - rowMeans(XB_data[c(13:16,25:29)]))))) +
  geom_point() + geom_smooth(method = "loess") + theme_classic()
g

# sort in reverse order according to the abs of the difference in the 
# mean expression values
reordered_XBdata <- XB_data[order(-abs(rowMeans(XB_data[c(9:12,17:24)]) - 
                                        rowMeans(XB_data[c(13:16,25:29)]))),]
head(reordered_XBdata)

# which ones have an abs difference greater than some threshold?

# first plot the distribution of differences
# limit the y axis to a small value so we can focus on the rare ones
hist(abs(rowMeans(XB_data[c(9:12,17:24)]) - 
            rowMeans(XB_data[c(13:16,25:29)])), breaks = seq(0,2500,10), ylim = c(0,20), xlim = c(400,2500))

# a reasonable threshold is around 400
# make a vector out of the absmean
absmean <- abs(rowMeans(XB_data[c(9:12,17:24)]) - 
                 rowMeans(XB_data[c(13:16,25:29)]))
mean_females <- rowMeans(XB_data[c(9:12,17:24)])
mean_males <- rowMeans(XB_data[c(13:16,25:29)])  


absmean_scaffolds <- abs(rowMeans(XB_scaffold_data[c(9:12,17:24)]) - 
                 rowMeans(XB_scaffold_data[c(13:16,25:29)]))
mean_females_scaffolds <- rowMeans(XB_scaffold_data[c(9:12,17:24)])
mean_males_scaffolds <- rowMeans(XB_scaffold_data[c(13:16,25:29)])  


# add this as a column to XB data
XB_data <- cbind(XB_data, absmean, mean_females, mean_males)
XB_scaffold_data <- cbind(XB_scaffold_data, absmean_scaffolds, mean_females_scaffolds, mean_males_scaffolds)

hist(XB_data$absmean, breaks = seq(0,2500,50), ylim = c(0,20), xlim = c(0,2500))
hist(XB_scaffold_data$absmean, breaks = seq(0,2500,50), ylim = c(0,20), xlim = c(0,2500))


# get the values over some threshold
dim(XB_data[XB_data$absmean > 400, ])
dim(XB_scaffold_data[XB_scaffold_data$absmean > 400, ])
# there are only 16 of these genes mapped, and 4 on scaffolds
# how many map to < 60,000,000?
# if we decrease to 300, we have 18 and 5 respectively
dim(XB_data[(XB_data$absmean > 300) & (XB_data$start < 60000000), ])
# only 13 mapped.
dim(XB_data[(XB_scaffold_data$absmean > 300) & (XB_scaffold_data$start < 60000000), ])
# and still 5 on scaffolds

########
########
########
########
########
########
# so use these 18 and 8 genes because they have major expression level differences
# between males and females
XB_data$trans_id[(XB_data$absmean > 300) & (XB_data$start < 60000000)]
#[1] TRINITY_DN2052_c0_g1_i4 TRINITY_DN2052_c0_g1_i5 TRINITY_DN2052_c0_g1_i1 TRINITY_DN2052_c0_g1_i2
#[5] TRINITY_DN2052_c0_g1_i3 TRINITY_DN2060_c0_g1_i3 TRINITY_DN2060_c0_g1_i2 TRINITY_DN7245_c0_g2_i3
#[9] TRINITY_DN2484_c0_g2_i2 TRINITY_DN2734_c1_g1_i1 TRINITY_DN3153_c0_g3_i1 TRINITY_DN3120_c0_g1_i2
#[13] TRINITY_DN3120_c0_g1_i6 TRINITY_DN2850_c0_g2_i1 TRINITY_DN46_c0_g1_i1   TRINITY_DN46_c0_g1_i6  
#[17] TRINITY_DN7190_c0_g1_i2 TRINITY_DN2734_c1_g1_i2 
# only 11 genes really
XB_scaffold_data$trans_id[XB_scaffold_data$absmean > 300]
#[1] TRINITY_DN10401_c0_g1_i7 TRINITY_DN2660_c0_g1_i1  TRINITY_DN35008_c0_g3_i1 TRINITY_DN1860_c0_g1_i6 
#[5] TRINITY_DN16646_c0_g1_i1
# only 5 genes 
########
########
########
########
########
########
XB_data[(XB_data$trans_id == "TRINITY_DN12136_c1_g1_i1") | (XB_data$trans_id == "TRINITY_DN10625_c1_g1_i1"),]

head(XB_data[(XB_data$absmean > 400) & (XB_data$start < 60000000), ])

XB_data$trans_id[(XB_data$absmean > 400) & (XB_data$start < 60000000)]
highabsmeandifference_on_XBchr8L_SL_region <- as.vector(XB_data$trans_id[(XB_data$absmean > 400) & (XB_data$start < 60000000)])

# we can also select genes with zero expression in males
nomaleexression_on_XBchr8L_SL_region <- as.vector((XB_data$trans_id[(XB_data$mean_males == 0) & (XB_data$start < 60000000)]))
# there are 207 of these. 
nomaleexression_on_scaffolds_SL_region <- as.vector((XB_scaffold_data$trans_id[(XB_scaffold_data$mean_males == 0)]))
# there are 39 on scaffolds
# but some of these also have very low expression in females and are possibly errors
# so let's select ones that are also expressed at least an average of 1 read per female
nomaleexression_and_somefemaleexpression_on_XBchr8L_SL_region <- as.vector((XB_data$trans_id
                                                   [(XB_data$mean_males == 0) &
                                                       (XB_data$mean_females >= 60) 
                                                     & (XB_data$start < 60000000)]))

# there are only 10 of these that map to chr8L
# [1] "TRINITY_DN62490_c0_g1_i2" "TRINITY_DN6717_c0_g1_i4"  "TRINITY_DN17402_c0_g1_i2" "TRINITY_DN3120_c0_g1_i2" 
# [5] "TRINITY_DN6425_c0_g1_i4"  "TRINITY_DN24786_c0_g1_i2" "TRINITY_DN2835_c0_g1_i1"  "TRINITY_DN248_c1_g1_i1"  
# [9] "TRINITY_DN248_c1_g1_i3"


str(XB_data)
XB_data[XB_data$trans_id %in% nomaleexression_and_somefemaleexpression_on_XBchr8L_SL_region, ]


nomaleexression_and_somefemaleexpression_on_scaffolds_SL_region <- as.vector(XB_scaffold_data$trans_id
                    [(XB_scaffold_data$mean_males == 0) &  (XB_scaffold_data$mean_females >= 24) ])

# there are none of these that map to scaffolds


# we can also select genes with very low expression in males 
# such as an average of less than 1 reads per individual
almost_nomaleexression_on_XBchr8L_SL_region <- as.vector((XB_data$trans_id[(XB_data$mean_males < 10) &
                                                                            (XB_data$start < 60000000)]))
# but there are many of these 3438
# how many of them are expressed a bit more in females?
almost_nomaleexression_andsomeinfemalees_on_XBchr8L_SL_region <- as.vector((XB_data$trans_id[(XB_data$mean_males < 10) &
      (XB_data$mean_females > 60) & (XB_data$start < 60000000)]))
# only 19 of these
# [1] "TRINITY_DN62490_c0_g1_i2" "TRINITY_DN6717_c0_g1_i4"  "TRINITY_DN2060_c0_g1_i3"  "TRINITY_DN5168_c0_g2_i2" 
# [5] "TRINITY_DN2734_c1_g1_i1"  "TRINITY_DN42088_c0_g1_i1" "TRINITY_DN17402_c0_g1_i2" "TRINITY_DN8455_c0_g1_i2" 
# [9] "TRINITY_DN8455_c0_g1_i5"  "TRINITY_DN3120_c0_g1_i2"  "TRINITY_DN3120_c0_g1_i6"  "TRINITY_DN6425_c0_g1_i4" 
# [13] "TRINITY_DN24786_c0_g1_i2" "TRINITY_DN2835_c0_g1_i1"  "TRINITY_DN2252_c4_g1_i1"  "TRINITY_DN46_c0_g1_i3"   
# [17] "TRINITY_DN204_c0_g1_i1"   "TRINITY_DN248_c1_g1_i1"   "TRINITY_DN248_c1_g1_i3"  
plot(sort(XB_data[XB_data$trans_id %in% almost_nomaleexression_andsomeinfemalees_on_XBchr8L_SL_region, ]$start))

# check if the more relaxed vector includes all of the strict vector
nomaleexression_and_somefemaleexpression_on_XBchr8L_SL_region %in% almost_nomaleexression_andsomeinfemalees_on_XBchr8L_SL_region


########
########
########
########
########
########
# ok looks good, so use these 43 genes:
almost_nomaleexression_andsomeinfemalees_on_XBchr8L_SL_region
#[1] "TRINITY_DN1286_c1_g1_i3"  "TRINITY_DN8282_c1_g1_i3"  "TRINITY_DN62490_c0_g1_i2" "TRINITY_DN62490_c0_g1_i5"
#[5] "TRINITY_DN3790_c0_g2_i5"  "TRINITY_DN6717_c0_g1_i4"  "TRINITY_DN2060_c0_g2_i10" "TRINITY_DN2060_c0_g2_i8" 
#[9] "TRINITY_DN2060_c0_g2_i6"  "TRINITY_DN2060_c0_g1_i3"  "TRINITY_DN7292_c0_g1_i3"  "TRINITY_DN7293_c0_g1_i3" 
#[13] "TRINITY_DN19793_c0_g1_i1" "TRINITY_DN56299_c0_g1_i2" "TRINITY_DN5130_c0_g1_i2"  "TRINITY_DN5168_c0_g2_i2" 
#[17] "TRINITY_DN2734_c1_g1_i1"  "TRINITY_DN42088_c0_g1_i1" "TRINITY_DN10632_c0_g1_i2" "TRINITY_DN15913_c1_g1_i5"
#[21] "TRINITY_DN24807_c0_g1_i1" "TRINITY_DN17402_c0_g1_i2" "TRINITY_DN8455_c0_g1_i2"  "TRINITY_DN8455_c0_g1_i5" 
#[25] "TRINITY_DN12908_c0_g1_i2" "TRINITY_DN5865_c0_g1_i13" "TRINITY_DN3120_c0_g1_i2"  "TRINITY_DN3120_c0_g1_i6" 
#[29] "TRINITY_DN1663_c3_g2_i2"  "TRINITY_DN6425_c0_g1_i4"  "TRINITY_DN6603_c0_g1_i9"  "TRINITY_DN8654_c1_g1_i2" 
#[33] "TRINITY_DN24786_c0_g1_i2" "TRINITY_DN2835_c0_g1_i1"  "TRINITY_DN3463_c1_g2_i10" "TRINITY_DN3417_c5_g1_i3" 
#[37] "TRINITY_DN29114_c0_g1_i3" "TRINITY_DN2252_c4_g1_i1"  "TRINITY_DN46_c0_g1_i2"    "TRINITY_DN46_c0_g1_i3"   
#[41] "TRINITY_DN204_c0_g1_i1"   "TRINITY_DN248_c1_g1_i1"   "TRINITY_DN248_c1_g1_i3"   "TRINITY_DN2674_c0_g1_i14"
# actually ~35
########
########
########
########
########
########

almost_nomaleexression_andsomeinfemalees_on_scaffolds_SL_region <- as.vector((XB_scaffold_data$trans_id[(XB_scaffold_data$mean_males < 10) &
       (XB_scaffold_data$mean_females > 36)]))
# 1 of these
# [1] "TRINITY_DN54315_c0_g1_i3"





# are some individual males more likely than others to have expression
# when the others do not?  We need to count up the number of times each
# male has expression when all others do not.
# A starting point is to subset the data where the average of a logical of whether
# male expression is zero is 1/9 (i.e., one T,and 8 Falses)
# mean(c(T,F,F,F,F,F,F,F,F)) = 0.1111111
male_one_expression_flag <- NULL
# first make a vector with this value

# I am having trouble doing this so I am going to do it 
# a slow way

male_expression_logicals <- cbind((XB_data$male_stage46_XBO12 == 0),(XB_data$male_stage46_XBO16 == 0),
  (XB_data$male_stage46_XBO23 == 0), (XB_data$male_stage46_XBO26 == 0),
  (XB_data$male_stage48_XBO17 == 0), (XB_data$male_stage48_XBO28 == 0),
 (XB_data$male_stage48_XBO31 == 0), (XB_data$male_stage48_XBO32 == 0), 
  (XB_data$male_stage48_XBO36 == 0))


female_expression_logicals <- cbind((XB_data$female_stage46_XBO15 == 0),(XB_data$female_stage46_XBO24 == 0),
                                  (XB_data$female_stage46_XBO27 == 0), (XB_data$female_stage46_XBO8 == 0),
                                  (XB_data$female_stage48_XBO19 == 0), (XB_data$female_stage48_XBO20 == 0),
                                  (XB_data$female_stage48_XBO21 == 0), (XB_data$female_stage48_XBO29 == 0), 
                                  (XB_data$female_stage48_XBO30 == 0), (XB_data$female_stage48_XBO33 == 0),
                                  (XB_data$female_stage48_XBO34 == 0), (XB_data$female_stage48_XBO35 == 0))

head(male_expression_logicals)
male_one_expression_flag <- rowMeans(male_expression_logicals)
XB_data <- cbind(XB_data,male_one_expression_flag)

# the other thing we can now try is to identify female specific genes but allowing for one male
# to have any level of expression. This can be done using the XB_data$male_one_expression_flag
dim(XB_data[(XB_data$male_one_expression_flag > 0)&(XB_data$male_one_expression_flag < 0.1112)&(XB_data$start < 60000000), ])
# there are 340 genes with a only one male expressed

# how many of these also have females with a reasonable expression level?
dim(XB_data[(XB_data$male_one_expression_flag > 0)&(XB_data$male_one_expression_flag < 0.1112)&(XB_data$start < 60000000)& (XB_data$mean_females > 12), ])
# only 55


# lets look at the column sums to see if any one individual is weird
colSums(male_expression_logicals[(XB_data$male_one_expression_flag > 0)&(XB_data$male_one_expression_flag < 0.1112),])
# [1]  99 165 100 125  48  53  42  45  52
# there are more in stage 46 where only 1 male  has expression
# than in stage 48, but the levels are roughly even


head(female_expression_logicals)
female_one_expression_flag <- rowMeans(female_expression_logicals)
XB_data <- cbind(XB_data,female_one_expression_flag)
dim(XB_data[(XB_data$female_one_expression_flag > 0)&(XB_data$female_one_expression_flag < 0.1112), ])
# there are 654 of these genes with a only one male expressed

# lets look at the column sums to see if any one individual is weird
colSums(female_expression_logicals[(XB_data$female_one_expression_flag > 0)&(XB_data$female_one_expression_flag < 0.1112),])
# [1]  62 170  32  59 114  49  29  21  31  33  15  39
# there might be some weirdness with
# female_stage46_XBO24 and female_stage48_XBO19
```

I then used perl one liners to extract the seqs from Xue's transcriptome (on evanslab: /home/evanslab/borealis_tadpole_transcriptome/transcriptome):
```perl
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' females_high_males_almost_zero borealis_tad_goand_transcriptome.fasta > XB_probe_seqs.fa
```
because 3 were missing, I did this too:
```perl
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(TRINITY_DN8455_c0_g1_i2 TRINITY_DN6425_c0_g1_i4 TRINITY_DN46_c0_g1_i3)}print if $c' borealis_tad_goand_transcriptome.fasta  >> XB_probe_seqs.fa
```



