# W or Z

I'd like to try to get genotypes from RNAseq data from the mother and father and see if these can be used, at least for differentially expressed transcripts, to try to infer contributions of the W-and Z- alleles in females (and maybe males too).

Trimmed RNAseq adult liver data are here:
```
/home/evanslab/borealis_adultFamily_RNAseq/data/trimmed
```

and the transcriptome assembly is here:
```
/home/evanslab/borealis_adultFamily_RNAseq/data/transcriptome/combine_de_count_transcript/borealis_transcriptome_trinityOut.fasta
```

map reads (dad):
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3896_dad_liver_R1_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3896_dad_testis_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3896_dad_liver_R2_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3896_dad_testis_R2_paired.fastq.gz' > BJE3896_dad_testis_liver.sam
```
mom:
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3897_mom_liver_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3897_mom_liver_R2_paired.fastq.gz' > BJE3896_mom_liver.sam
```
son1:
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3929_boy_liver_R1_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3929_boy_testis_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3929_boy_liver_R2_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE3929_boy_testis_R2_paired.fastq.gz' > BJE3929_boy_testis_liver.sam
```
son2:
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4017_boy_liver_R1_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4017_boy_testis_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4017_boy_liver_R2_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4017_boy_testis_R2_paired.fastq.gz' > BJE4017_boy_testis_liver.sam
```
son3
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4039_boy_liver_R1_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4039_boy_testis_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4039_boy_liver_R2_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4039_boy_testis_R2_paired.fastq.gz' > BJE4039_boy_testis_liver.sam
```
daughter1
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4009_girl_liver_R1_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4009_girl_oviduct_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4009_girl_liver_R2_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4009_girl_oviduct_R2_paired.fastq.gz' > BJE4009_girl_oviduct_liver.sam
```
daughter2
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4072_girl_liver_R1_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4072_girl_oviduct_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4072_girl_liver_R2_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4072_girl_oviduct_R2_paired.fastq.gz' > BJE4072_girl_oviduct_liver.sam
```
daughter3
```
bwa mem borealis_transcriptome_trinityOut.fasta '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4082_girl_liver_R1_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4082_girl_oviduct_R1_paired.fastq.gz' '<zcat /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4082_girl_liver_R2_paired.fastq.gz /home/evanslab/borealis_adultFamily_RNAseq/data/trimmed/BJE4082_girl_oviduct_R2_paired.fastq.gz' > BJE4082_girl_oviduct_liver.sam
```


make bam files
```
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE3896_mom_liver.bam BJE3896_mom_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE3896_dad_testis_liver.bam BJE3896_dad_testis_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4009_girl_oviduct_liver.bam BJE4009_girl_oviduct_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE3929_boy_testis_liver.bam BJE3929_boy_testis_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4072_girl_oviduct_liver.bam BJE4072_girl_oviduct_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4017_boy_testis_liver.bam BJE4017_boy_testis_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4082_girl_oviduct_liver.bam BJE4082_girl_oviduct_liver.sam
samtools view -bt borealis_transcriptome_trinityOut.fasta -o BJE4039_boy_testis_liver.bam BJE4039_boy_testis_liver.sam
```


delete sam files
```
rm -f BJE3896_mom_liver.sam
```

sort bam files
```
samtools sort BJE3896_mom_liver.bam -o BJE3896_mom_liver_sorted.bam
samtools sort BJE4009_girl_oviduct_liver.bam -o BJE4009_girl_oviduct_liver_sorted.bam
samtools sort BJE3929_boy_testis_liver.bam -o BJE3929_boy_testis_liver_sorted.bam
samtools sort BJE3896_dad_testis_liver.bam -o BJE3896_dad_testis_liver_sorted.bam
samtools sort BJE4072_girl_oviduct_liver.bam -o BJE4072_girl_oviduct_liver_sorted.bam
samtools sort BJE4017_boy_testis_liver.bam -o BJE4017_boy_testis_liver_sorted.bam
samtools sort BJE4082_girl_oviduct_liver.bam -o BJE4082_girl_oviduct_liver_sorted.bam
samtools sort BJE4039_boy_testis_liver.bam -o BJE4039_boy_testis_liver_sorted.bam
```
