# Population and Quantitative GENetic toolbox

**pqgen** calculates population and quantitative genetic statistics from genomic data. The focus is on speed and ease of use. Statistics are calculated from the numeric genotype format (GT) used in the Variant Call Format (VCF), making it simple to go from variant calling to population genetic statistics in a few commands. The input format is flexible, but by default uses an extended bed format that makes it trivial to integrate bedtools into the analysis pipeline if, for example, you would like to calculate popgen stats for coding vs non-coding regions, or for synonymous vs non-synonymous sites.

## Preparing an input file

**pqgen** by default expects an input file in an extended bed format:
```bash
  	chrom    start    end    name    GT.1    GT.2  ...  GT.n
```
where GT.1 to GT.n are the genotypes of each individual given in GT format. For example:
```bash
       Contig0    <start>    1    <name>    0/0    0/0    0/1    0/0    0/1
       Contig0    <start>    2    <name>    0/0    0/0    0/0    0/0    0/0
       Contig0    <start>    3    <name>    0/1    0/1    0/0    0/0    0/1
```
Note that the start position is not required by **pqgen** to calculate stats but is included to maintain the bed format. This format can be generated from a vcf file using bcftools:
```bash
    bcftools query -f "%CHROM\t<start>\t%POS\t<name>[\t%GT]\n" aln.vcf > aln.gt.bed
```
Here, we are writing the output to a new file aln.gt.bed.

## Calculating statistics

To calculate some basic site frequency based stats such as the number of segregating sites and Tajima's D use the theta command:
```bash
$ pqgen theta aln.gt.bed
Contig0	0	1372378	Contig0	10	1367408	15844	0.004096	0.004072	-0.029052
Contig1	0	865787	Contig1	10	862828	12965	0.005312	0.005468	0.148427
Contig2	0	786857	Contig2	10	783622	15404	0.006949	0.007175	0.163521
```
By default, statistics will be calculated per chromosome (or contig in this case), but this can be changed (see options below).

All commands follow the same basic usage:

```bash
  pqgen [--help] [--version] <command> [OPTIONS] [FILE]
```
where `<command>` can be anyone of the following:
```bash
  het	       Calculates heterozygosity.
  sfs	       Outputs the folded site frequency spectrum.
  theta        Calculates site frequency based statistics.
```
and `OPTIONS` can include any of the following:

```bash
Options:

  -c <int>
     Column number (1-indexed) giving the chromosome/contig/scaffold ID.

  -d <char>
     Delimiter separating each column. Must be enclosed in single quotes, spaces
     must be escaped (e.g. -d '\ '). Tab delimited by default.

  -f <int>
     Column number (1-indexed) of the factor over which the stats should be calculated.
     By default, statistics are output per chromosome.

  -k <str>
     1-indexed column numbers indicating the columns over which the statistic
     should be calculated. Columns can be separated by commas, and ranges specified
     using dashes. For example, -k 2,5-7 would specify columns 2, 5, 6 and 7.

  -p <int>
     Column number (1-indexed) giving the 1-indexed reference position.
```

## Installation

Download and unpack the latest version of **pqgen** (replacing *X*, *Y* and *Z* with the version number):
```bash
wget https://github.com/mestocks/pqgen/releases/download/vX.Y.Z/pqgen-X.Y.Z.tar.gz
tar -zxvf pqgen-X.Y.Z.tar.gz
```
Then change into the **pqgen** directory, compile the source and install it:
```bash
cd pqgen-X.Y.Z
./configure
make
make install
```

The installation path can be changed using './configure --prefix=<path/to/dir>'. Admin permissions may be required for the 'make install' step.


## Bonus scripts

dna2codon - Collapse consecutive dna sequences into codons.

	Input format:
	```bash
	"Contig0    <pos0>    1    <name>    0.65	+	A    G	..."
	```

dna2div - Calculate divergence based statistics. A maximum of two alleles should be given.

	Input format:
	```bash
	"Contig0    <pos0>    1    <name>    A	G"
	```

codon2pnds - Count the number of silent and replacement substitutions and polymorphisms. The
	     first codon is the outgroup/ancestral codon.

	Input format:
	```bash
	"Contig0    <pos0>    1    <name>    ATG    ATG	..."
	```


