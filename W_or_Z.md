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
