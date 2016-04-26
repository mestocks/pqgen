# Population and Quantitative Genetic Tools

### pqtheta
```bash
... | pq-theta nsam [fcol]
```

###### dna2codon

Convert nucleotide sequences into codons
```bash
... | dna2codon

Input:
chr    start    end    name    score    +|-    nuc.1    nuc.2 ... nuc.n

Output:
chr    codon.start    codon.end    name    score    +|-    codon.1    codon.2 ... codon.n    
```

Example input:
```bash
chr1      676151  676152  transcript1       0.65   +       A       A
chr1      676152  676153  transcript1       0.65   +       T       T
chr1      676153  676154  transcript1       0.65   +       G       G
chr1      676154  676155  transcript1       0.65   +       T       T
chr1      676155  676156  transcript1       0.65   +       C       C
chr1      676156  676157  transcript1       0.65   +       G       G
chr1      676157  676158  transcript1       0.65   +       A       A
chr1      676158  676159  transcript1       0.65   +       C       C
chr1      676159  676160  transcript1       0.65   +       G       A
```
Example output:
```bash
chr1      676151  676154  transcript1       0.65   +       ATG     ATG
chr1      676154  676157  transcript1       0.65   +       TCG     TCG
chr1      676157  676160  transcript1       0.65   +       ACG     ACA
```
