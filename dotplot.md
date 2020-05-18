# Extract chromosome fasta files

In this directory on info: `/4/ben`
```
zcat XENTR_10.0_genome.fasta.gz | perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(Chr8)}print if $c' > XT_v10_chr8.fasta
```

I think they are softmasked, so I hard masked them like this (using infile editing, so I made a copy first):
```
sed -i.bak -e 's/a/N/g' -e 's/c/N/g' -e 's/t/N/g' -e 's/g/N/g' XL_v9_1_chr8L_hardmasked.fasta
```

# Gepard (program to make dotplots)

I'm working in this directory: `/4/ben/dotplot`

```
java -Xmx8192M -cp Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 /4/ben/XL_v9.1_repeat_masked/XL_v9_1_chr8L_repeatmasked.fasta -seq2 /4/ben/XL_v9.1_repeat_masked/XL_v9_1_chr8S_repeatmasked.fasta  -matrix /4/ben/dotplot/gepard/resources/matrices/edna.mat -outfile XLchr8L_8S.png -word 20

java -Xmx8192M -cp Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 /4/ben/XL_v9.1_repeat_masked/XL_v9_1_chr8L_repeatmasked.fasta -seq2 /4/ben/XT_v10_repeat_masked/XT_v10_chr8_repeatmasked.fasta -matrix /4/ben/dotplot/gepard/resources/matrices/edna.mat -outfile XLchr8L_XTv10_8.png

java -Xmx8192M -cp Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 /4/ben/XL_v9.1_repeat_masked/XL_v9_1_chr8L_repeatmasked.fasta -seq2 /4/ben/XT_v9.0_repeat_masked/XT_v9_0_chr8_repeatmask.fasta -matrix /4/ben/dotplot/gepard/resources/matrices/edna.mat -outfile XLchr8L_XTv9_8.png

java -Xmx8192M -cp Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 /4/ben/XL_v9.1_repeat_masked/XL_v9_1_chr8S_repeatmasked.fasta -seq2 /4/ben/XT_v9.0_repeat_masked/XT_v9_0_chr8_repeatmask.fasta -matrix /4/ben/dotplot/gepard/resources/matrices/edna.mat -outfile XLchr8S_XTv9_8_word40.png -word 40

java -Xmx8192M -cp Gepard-1.40.jar org.gepard.client.cmdline.CommandLine -seq1 /4/ben/XL_v9.1_repeat_masked/XL_v9_1_chr8S_repeatmasked.fasta -seq2 /4/ben/XT_v10_repeat_masked/XT_v10_chr8_repeatmasked.fasta -matrix /4/ben/dotplot/gepard/resources/matrices/edna.mat -outfile XLchr8S_XTv9_8_word40.png -word 40
```
