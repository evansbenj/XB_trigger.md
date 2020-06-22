# W or Z

I'd like to try to get genotypes from RNAseq data from the mother and father and see if these can be used, at least for differentially expressed transcripts, to try to infer contributions of the W-and Z- alleles in females (and maybe males too).

Trimmed RNAseq adult liver data are here:
```
/home/evanslab/borealis_adultFamily_RNAseq/trimmed
```

and the liver transcriptome assembly is here:
```
/home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta
```
Trimmed data from tads are here:
```
/home/evanslab/borealis_tadpole_transcriptome/data/trimmed/2_trimmed_data
```
and the transcriptome assembly is here:
```
/home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta
```


map reads (dad) in this directory `/home/evanslab/borealis_adultFamily_RNAseq/trimmed`:
```
bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE3896_dad_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE3896_dad_liver_R2_paired.fastq.gz'> BJE3896_dad_liver.sam

bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE3897_mom_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE3897_mom_liver_R2_paired.fastq.gz' > BJE3897_mom_liver.sam

bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE3929_boy_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE3929_boy_liver_R2_paired.fastq.gz' > BJE3929_boy_liver.sam

bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4017_boy_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4017_boy_liver_R2_paired.fastq.gz' > BJE4017_boy_liver.sam

bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4039_boy_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4039_boy_liver_R2_paired.fastq.gz' > BJE4039_boy_liver.sam

bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4009_girl_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4009_girl_liver_R2_paired.fastq.gz' > BJE4009_girl_liver.sam

bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4072_girl_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4072_girl_liver_R2_paired.fastq.gz' > BJE4072_girl_liver.sam

bwa mem /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4082_girl_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/trimmed/BJE4082_girl_liver_R2_paired.fastq.gz' > BJE4082_girl_liver.sam

```


make bam files
```
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE3896_mom_liver.bam BJE3897_mom_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE3896_dad_liver.bam BJE3896_dad_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4009_girl_liver.bam BJE4009_girl_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE3929_boy_liver.bam BJE3929_boy_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4072_girl_liver.bam BJE4072_girl_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4017_boy_liver.bam BJE4017_boy_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4082_girl_liver.bam BJE4082_girl_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4039_boy_liver.bam BJE4039_boy_liver.sam
```


delete sam files
```
rm -f BJE3896_mom_liver.sam
```

sort bam files
```
samtools sort BJE3896_mom_liver.bam -o BJE3897_mom_liver_sorted.bam
samtools sort BJE4009_girl_liver.bam -o BJE4009_girl_liver_sorted.bam
samtools sort BJE3929_boy_liver.bam -o BJE3929_boy_liver_sorted.bam
samtools sort BJE3896_dad_liver.bam -o BJE3896_dad_liver_sorted.bam
samtools sort BJE4072_girl_liver.bam -o BJE4072_girl_liver_sorted.bam
samtools sort BJE4017_boy_liver.bam -o BJE4017_boy_liver_sorted.bam
samtools sort BJE4082_girl_liver.bam -o BJE4082_girl_liver_sorted.bam
samtools sort BJE4039_boy_liver.bam -o BJE4039_boy_liver_sorted.bam
```
make bai index for bam files
```
samtools index BJE3896_mom_liver_sorted.bam
```
call genotypees do light filter and outout all calls, including homoz ref
```
Use samtools and bcftools to call genotypes and filter

samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE3897_mom_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE3897_mom_liver_sorted.vcf.gz

samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE3896_dad_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE3896_dad_liver_sorted.vcf.gz

samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE3929_boy_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE3929_boy_liver_sorted.vcf.gz

samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE4009_girl_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE4009_girl_liver_sorted.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE4017_boy_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE4017_boy_liver_sorted.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE4039_boy_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE4039_boy_liver_sorted.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE4072_girl_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE4072_girl_liver_sorted.vcf.gz
samtools mpileup -d8000 -ugf /home/evanslab/borealis_adultFamily_RNAseq/borealis_liver_transcriptome/build_transcriptome/borealis_adult_liver_transcriptome_trinityout.fasta -t DP,AD BJE4082_girl_liver_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o BJE4082_girl_liver_sorted.vcf.gz

```
also for tad RNAseq:
```
samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO12_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO12_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO19_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO19_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO24_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO24_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO29_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO29_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO33_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO33_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO8_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO8_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO15_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO15_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO20_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO20_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO26_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO26_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO30_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO30_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO34_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO34_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO16_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO16_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO21_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO21_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO27_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO27_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO31_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO31_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO35_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO35_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO17_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO17_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO23_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO23_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO28_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO28_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO32_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO32_sorted.bam.vcf.gz


samtools mpileup -d8000 -ugf /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -t DP,AD XBO36_sorted.bam | bcftools call --gvcf 5 -V indels --format-fields GQ -m -O z | bcftools filter -e 'QUAL < 20 | DP < 4 | MQ < 20' -O z -o XBO36_sorted.bam.vcf.gz

```


make index for vcf files
```
tabix -p vcf BJE4009_girl_oviduct_liver_sorted.bam.vcf.gz
```
merge vcf files
```
bcftools merge BJE3897_mom_liver_sorted.bam.vcf.gz BJE3896_dad_testis_liver_sorted.bam.vcf.gz BJE4082_girl_oviduct_liver_sorted.bam.vcf.gz BJE4072_girl_oviduct_liver_sorted.bam.vcf.gz BJE4009_girl_oviduct_liver_sorted.bam.vcf.gz BJE4039_boy_testis_liver_sorted.bam.vcf.gz BJE4017_boy_testis_liver_sorted.bam.vcf.gz BJE3929_boy_testis_liver_sorted.bam.vcf.gz -Oz -o Merged.vcf.gz
```
```
bcftools merge BJE3897_mom_liver_sorted.vcf.gz BJE4009_girl_liver_sorted.vcf.gz BJE4072_girl_liver_sorted.vcf.gz BJE4082_girl_liver_sorted.vcf.gz BJE3896_dad_liver_sorted.vcf.gz BJE3929_boy_liver_sorted.vcf.gz BJE4017_boy_liver_sorted.vcf.gz BJE4039_boy_liver_sorted.vcf.gz -Oz -o XB_liver_RNAseq_Merged.vcf.gz
```

extract allele depth information
```
vcftools --gzvcf Merged.vcf.gz --extract-FORMAT-info AD
```
After renaming, the output file is called `Merged.vcf.gz_out.AD.FORMAT` and it has exactly what I need - allele depth for all individuals that have genotypes for each transcript for each SNP.

Then I went back to R to get the transcript IDs of significantly female biased and male biased transcripts:
```R
borTad_laevisGenome_edgeR_tpm_SL_femalebiased <- borTad_laevisGenome_edgeR_tpm_combine_st46 %>% filter((FDR <= 0.05)& (logFC < -2)) %>% 
  filter(chromosome == "chr8L") %>%
  filter(start<=57000000) 
borTad_laevisGenome_edgeR_tpm_SL_femalebiased$trans_id
write.csv(borTad_laevisGenome_edgeR_tpm_SL_femalebiased$trans_id,
          "borTad_laevisGenome_edgeR_tpm_SL_femalebiased$trans_id.csv", row.names = F)
```

```R
# get the transcript_ids of female_biased trnascripts in the SL region
borTad_laevisGenome_edgeR_tpm_SL_malebiased <- borTad_laevisGenome_edgeR_tpm_combine_st46 %>% filter((FDR <= 0.05)& (logFC > -2)) %>% 
  filter(chromosome == "chr8L") %>%
  filter(start<=57000000) 
borTad_laevisGenome_edgeR_tpm_SL_malebiased$trans_id
write.csv(borTad_laevisGenome_edgeR_tpm_SL_malebiased$trans_id,
          "borTad_laevisGenome_edgeR_tpm_SL_malebiased$trans_id.csv", row.names = F)
```
And used this list to grep the allele depth data of the individual transcripts from the files generated above with vcftools

```
grep -f SL_trans_IDs_sig_Sex_biased Merged.vcf.gz_out.AD.FORMAT > SL_allelic_depth_sig_sex_biased.AD.FORMAT
```

# SNPs in RNAseq data
OK using the approach above I mapped RNAseq data to the de novo transcriptome and called homoz and heteroz sites.  For XB I hypothesize that the female biased expression is because the Z chromosome evolved from a degenrate Y chromosome. If this is true, I predict that there should be very low variation on female biased transcripts (or really any transcript) on the SL region compared to the non-SL region.  

So I need to calculate pi for all transcripts and then test whether it is lower for female biased SL than for female biased nonSL. The reason this is expected is because a degenrate Z from a Y chr would mean that SL genes are mostly from the W.

This can be done for liver, and both XB tads.

Tracing SNPs in the liver RNAseq, which includes a mother, father, and offsping, doesn't help much: a SNP in the mother could be on the W- or Z-.  If it is in the mother and on the W and co-expressed with a Z-linked allele, we expect the W allele to be more highly expressed in the mother and daughers. When I looked at allele depth, most maternal variants (relative to the reference) seemed to be not expressed in males, making it difficult to tell whether this was female-biased regulation or degeneration.

```
zcat Merged_st46.vcf.gz | vcf-to-tab > out.tab
```
in this directory: `/home/evanslab/borealis_adultFamily_RNAseq/trimmed`:

females
```
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_liver_RNAseq_Merged.vcf.gz.tab 11111111 3_4_2_4_7_8 XB_liver_RNAseq_Merged_females.poly XB_liver_RNAseq_Merged_females.poly_by_windows
```
males
```
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_liver_RNAseq_Merged.vcf.gz.tab 11111111 3_4_1_3_5_6 XB_liver_RNAseq_Merged_males.poly XB_liver_RNAseq_Merged_males.poly_by_windows
```
both
```
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_liver_RNAseq_Merged.vcf.gz.tab 11111111 3_4_1_2_3_4_5_6_7_8 XB_liver_RNAseq_Merged_males_and_females.poly XB_liver_RNAseq_Merged_males_and_females.poly_by_windows
```

st46
```
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_st46_RNAseq_Merged.vcf.gz.tab 11111111 3_4_2_5_7_8 XB_st46_RNAseq_Merged_females.poly XB_st46_RNAseq_Merged_females.poly_by_windows

./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_st46_RNAseq_Merged.vcf.gz.tab 11111111 3_4_1_3_4_6 XB_st46_RNAseq_Merged_males.poly XB_st46_RNAseq_Merged_males.poly_by_windows

./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_st46_RNAseq_Merged.vcf.gz.tab 11111111 3_4_1_2_3_4_5_6_7_8 XB_st46_RNAseq_Merged_females_male_and_females.poly XB_st46_RNAseq_Merged_male_and_females.poly_by_windows
```

st48
```
./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_st48_RNAseq_Merged.vcf.gz.tab 1111111111111 3_4_2_3_4_6_7_10_11_12 XB_st48_RNAseq_Merged_females.poly XB_st48_RNAseq_Merged_females.poly_by_windows

./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_st48_RNAseq_Merged.vcf.gz.tab 1111111111111 3_4_1_5_8_9_13 XB_st48_RNAseq_Merged_males.poly XB_st48_RNAseq_Merged_males.poly_by_windows

./Boot_from_tab_diverge_poly_2018_allowmissingdata_transcripts.pl XB_st48_RNAseq_Merged.vcf.gz.tab 1111111111111 3_4_1_2_3_4_5_6_7_8_9_10_11_12_13 XB_st48_RNAseq_Merged_females_male_and_females.poly XB_st48_RNAseq_Merged_male_and_females.poly_by_windows
```

# Angsd

I decided to use angsd to look at polymorphism of individual male and female XB from different places.

The bam files mapped to XL_v9.2 are here:
```
/2/scratch/evanslab/2019_RADseq_KenyaXBXL_GhanaEastfamily/plate1/alinged
```
I indexed the bam file using `samtools faidx file.bam` and I indexed the reference genome using `bwa index XL9_2.fa.gz`.
As directed by the manual I added `-fold 1` because we don't have ancestral states.

Here are the commands I will try:
```
/home/evanslab/bin/angsd/angsd -i sort_Fem_Chesuwe_BJE4479_Xb.fq.gz.bam -doSaf 1 -anc XL9_2.fa.gz -fold 1 -r chr8L: -GL 1 -P 24 -out out 
/home/evanslab/bin/angsd/misc/realSFS out.saf.idx -fold 1 -r chr8L: -P 24 > out.sfs
#use -fold 1 in the above command if you dont have ancestral state.
/home/evanslab/bin/angsd/misc/realSFS saf2theta out.saf.idx -fold 1 -r chr8L: -outname out -sfs out.sfs
#Estimate for every Chromosome/scaffold
/home/evanslab/bin/angsd/misc/thetaStat -r chr8L: do_stat out.thetas.idx
#Do a sliding window analysis based on the output from the make_bed command.
/home/evanslab/bin/angsd/misc/thetaStat do_stat out.thetas.idx -win 50000 -step 10000 -r chr8L: -outnames theta.thetasWindow.gz
```

For tads:
```
/home/evanslab/bin/angsd/angsd -i /home/evanslab/borealis_tadpole_transcriptome/data/trimmed/2_trimmed_data/XBO12_sorted.bam -doSaf 1 -anc /home/evanslab/borealis_tadpole_transcriptome/data/transcriptome/borealis_tad_goand_transcriptome.fasta -fold 1 -GL 1 -P 24 -out out 
/home/evanslab/bin/angsd/misc/realSFS out.saf.idx -fold 1 -P 24 > out.sfs
#use -fold 1 in the above command if you dont have ancestral state.
/home/evanslab/bin/angsd/misc/realSFS saf2theta out.saf.idx -fold 1 -outname out -sfs out.sfs
#Estimate for every Chromosome/scaffold
/home/evanslab/bin/angsd/misc/thetaStat do_stat out.thetas.idx
```

# My polymorphism script

in this directory 
```
/2/scratch/evanslab/2019_RADseq_KenyaXBXL_GhanaEastfamily/plate1/genotyped
```

```
Each sex separately
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_1_2_3_4_6_7_8_15_16_24_25_26_27 Xb_fem_west.poly Xb_fem_west.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_44_45_46_47_48_49_50_55_62_63_64_67_73_74_75 Xb_mal_west.poly Xb_mal_west.poly_by_windows 
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_28_29_30_31_39_40_41_42_43 Xb_fem_east.poly Xb_fem_east.poly_by_windows 
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_76_84_85_86_87_88 Xb_mal_east.poly Xb_mal_east.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_32_33_34_35 Xb_fem_Njoro.poly Xb_fem_Njoro.poly_by_windows 
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_77_78_79_80 Xb_mal_Njoro.poly Xb_mal_Njoro.poly_by_windows

both sexes together
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_1_2_3_4_6_7_8_15_16_24_25_26_27_44_45_46_47_48_49_50_55_62_63_64_67_73_74_75 Xb_fem_and_male_west.poly Xb_fem_and_male_west.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_28_29_30_31_39_40_41_42_43_76_84_85_86_87_88 Xb_fem_and_male_east.poly Xb_fem_and_male_east.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_32_33_34_35_77_78_79_80 Xb_fem_and_male_Njoro.poly Xb_fem_and_male_Njoro.poly_by_windows

just one sex 
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_1_2_3_4_6_7_8_15_16_24_25_26_27_28_29_30_31_39_40_41_42_43_32_33_34_35 Xb_fem_all.poly Xb_fem_all.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_44_45_46_47_48_49_50_55_62_63_64_67_73_74_75_76_77_78_79_80_84_85_86_87_88 Xb_mal_all.poly Xb_mal_all.poly_by_windows 

just from eldoret
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_6 Xb_fem_BJE4471.poly Xb_fem_BJE4471.poly_by_windows
./Boot_from_tab_diverge_poly_2018_allowmissingdata.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L.tab 11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111 3_4_55 Xb_fem_BJE4473.poly Xb_fem_BJE4473.poly_by_windows
```
