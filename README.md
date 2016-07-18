# Population and Quantitative Genetic Tools

## Quick install

**pq-genetics** makes use of the library **librawk**. To install this dependency, go [here](https://github.com/mspopgen/librawk).

Once **librawk** is install, download and unpack the latest version of **pq-genetics** (*X*, *Y* and *Z* refer to version number):
```bash
wget https://github.com/mspopgen/pq-genetics/archive/vX.Y.Z.tar.gz
tar -zxvf vX.Y.Z.tar.gz
```
Then change into the **pq-genetics** directory, compile the source and install it:
```bash
cd pq-genetics-X.Y.Z
make
make install
```
This will install the compiled code into ```~/.local/``` and put configuration files into ```~/.config/pqgen/```. To call the commands from any directory, add ```~/.local/bin``` to your ```PATH``` environmental variable and ```~/.local/lib``` to your ```LD_LIBRARY_PATH```.

## Usage

##### pq-theta
*Calculate site frequency based stats*

```bash
... | pq-theta nsam [OPTIONS]

Input:
chr    start    end    name    nref    nalt

Output:
chr    region.start    region.end    name    nsam    nsites    seg.sites    thetaW    thetaPi    tajD

OPTIONS

  -b <bool>
    takes values 0 or 1, indicating whether theta values should be 
    give per base pair or summed over the entire region. [1]
  
  -f <int>
    column number (1-indexed) of the factor over which 
    the stats should be calculated. The default is to output 
    stats per chromosome, but the fourth name column could 
    be used instead to calculate over some group of features. [1]

```

##### pq-dna2codon
*Convert nucleotide sequences into codons*

```bash
... | pq-dna2codon

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

##### pq-codon2pnds
*Counts synonymous and non-synonymous sites from codons*

```bash
... | pq-codon2pnds [OPTIONS]

OPTIONS

  -a <STR>
     space delimited file where each line gives a "<codon> <amino_acid>"
     key-value pair determining which amino acid a codon translates to. ["share/codon2aa"]
  
  -f <int>
     column number (1-indexed) of the factor over which 
     the stats should be calculated. The default is to output 
     stats per chromosome, but the fourth name column could 
     be used instead to calculate over some group of features. [1]

  -s <STR>
     space delimited file where each line gives a "<codon> <syn_sites>"
     key-value pair determining how many possible synonymous mutations there are
     for each codon. ["share/codon2syn"]

Input:
chr    codon.start    codon.end    name    score    +|-    codon.REF    codon.1 ... codon.n

Output:
chr    codon.start    codon.end    name    score    +|-    ncodons    nvcodons    nsites.syn    nsites.nsyn    fix.syn    fix.nsyn    [poly.syn    poly.nsyn]

Example input:
chr1      676151  676154  transcript1       0.65   +       ATG     ATG
chr1      676154  676157  transcript1       0.65   +       TCG     TCG
chr1      676157  676160  transcript1       0.65   +       ACG     ACA

Example output:
chr1      676151  676154  transcript1       0.65   +       3	3	4.333333	4.666667	1	0
```
