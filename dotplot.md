# Extract chromosome fasta files
```
zcat XENTR_10.0_genome.fasta.gz | perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(Chr8)}print if $c' > XT_v10_chr8.fasta
```
