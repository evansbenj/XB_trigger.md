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

