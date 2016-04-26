# Population and Quantitative Genetic Tools

*Theta-based nucleotide diversity statistics*

```bash
... | pq-theta nsam [fcol]
```

*Annotation based stats*

```bash
... | dna2codon
```

Take input in bed format, in the form:
```bash
chr    start    end    name    score    strand    nucleotide
```

Output is for example:

```bash
chr1    0    3    transcipt1    1e10    +    ATG
chr1    3    6    transcipt1    1e10    +    TCG
chr1    6    9    transcipt1    1e10    +    ACG
```
