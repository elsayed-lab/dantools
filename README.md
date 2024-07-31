# DanTools

## Contact

Daniel Klimes (<daniel.s.klimes@gmail.com>)

## Overview

DanTools is a tool used for comparing genome sequences between
divergent organisms, such as different species. It does this by
creating a pseudogenome, a modification of a reference genome to match
the sequence of another. In doing so, a list of genetic variants
between given species is produced. Given a GFF file, these variants
can be annotated by their positions relative tofeatures. Additionally,
by keeping track of position shifts from indels, the annotations of
the original reference can be extended to the new species'
genome. This allows for the alignment of RNA-Seq data from divergent
species against a set of genomes with unified annotations. Critically,
DanTools does not require either genome to be well assembled, and can
even accept an RNA-Seq .fastq input.

## Methodology

DanTools in its base form accepts two genomes as input: the base and
the source. In its first step, DanTools breaks the source genome into
sequence "fragments", 50-10,000 base pieces of the genome and its
reverse complement. These fragments are then aligned against the base
genome using HISAT2 with lowered alignment stringency. Freebayes is
used to tabulate variants, which are then used to modify the base
genome. This process of alignment, variant collection, and base
modification is repeated until a satisfactory "modified" genome is
created.

After the final iteration, a new VCF file is generated using optimized
Needleman-Wunsch algorithms. This VCF can optionally be labeled by its
relative position within features and flanking regions. In the final
step, the input GFF (if provided) is then shifted to accomodate indel
mutations that were created in the process of modification.

Because of the flexibility of inputs, DanTools can accept RNA-Seq
reads in fastq format and minimally-assembled genomes.

DanTools additionally provides various helper functions used to
label variants relative to features, translate VCF files into
amino acids, predict translocations, and summarize variant
information (see Usage)

## Installation

### Prerequisites

One may wish to keep this in its own tree, perl provides a number of
mechanisms to do so, here is an example of how one might install
dantools in a self-contained tree with the github-downloaded source in
src/

```{bash, eval=FALSE}
mkdir -p dantools/202408
cd dantools/202408
git clone https://github.com/elsayed-lab/dantools.git
mv dantools src
## tell MakeMaker and Module::Build to use this tree
## All prerequisite perl libraries will go to dantools/202408/lib/perl5
export PERL5LIB="$(pwd)/lib/perl5:${PERL5LIB}"
export PERL_MM_OPT="INSTALL_BASE=$(pwd)"
export PERL_MB_OPT="--install_base $(pwd)"
## If you have App::cpanminus installed
cpanm ExtUtils::AutoInstall inc::Module::Install File::ShareDir::Install
## If not, then a similar invocation of cpan will work,
## you might alway want to check out local::lib
cd src
perl Makefile.PL && make && make install
## Depending on your system, you may need to manually install some dependencies
cpanm Moo::Role Parallel::ForkManager
make && make install
```

DanTools relies on a set of software packages not listed in the
makefile to perform some of its tasks:

- [HISAT2](https://github.com/DaehwanKimLab/hisat2)
- [Freebayes](https://github.com/freebayes/freebayes)
- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)
- [emboss](http://emboss.open-bio.org/html/use/ch02s07.html)

## Usage

Comparisons begin with the dantools pseudogenome command, which
generates the pseudogenome and, more importantly, a VCF file that
describes the variants in the source genome relative to the
base. Invocations differ by source data format:

Source is a genome in fasta format:

```{bash, eval=FALSE}
dantools pseudogen -b base.fasta -s source.fasta
```

Source is an already fragmented genome in fasta format:

```{bash, eval=FALSE}
dantools pseudogen -b base.fasta -s source.fasta --fragment no
```

Source is unpaired RNA-Seq reads in fastq format:

```{bash, eval=FALSE}
dantools pseudogen -b base.fasta --reads-u source.fastq
```

Source is paired-end RNA-Seq reads in fastq format:

```{bash, eval=FALSE}
dantools pseudogen -b base.fasta -1 mate1.fastq -2 mate2.fastq
```

Because many indels were applied during the pseudogenome creation
process, the base genome features need to be shifted in their relative
positions:

```{bash, eval=FALSE}
dantools shift -o shifted.gff -v variants.vcf -f base.gff
```

As an estimate for how well fragments/reads aligned to certain
features, the alignment depth can be summarized over features:

```{bash, eval=FALSE}
dantools summarize-depth -f shifted.gff -d depth.tsv --feature gene
```

Now, any RNA-Seq data aligned against the pseudogenome can be counted with
the shifted GFF file. It is recommended features with low fragment
alignment depth be removed as comparisons there were likely less accurate

One can also perform a variety of variant analyses with DanTools. If one
wants to label the variants according to a GFF file:

```{bash, eval=FALSE}
dantools label -v variants.vcf -f base.gff --features five_prime_UTR,CDS,three_prime_UTR
```

These variants can optionally be translated into amino acid changes
and scored by a mutation scoring matrix using the --translate option.
Both nucleotide and translated outputs can then be summarized by feature
using a set of helper functions:

```{bash, eval=FALSE}
dantools summarize-nuc labeled_nucleotides.tsv
dantools summarize-aa labeled_aa.tsv
```

For alignment of dantools fragment/pseudogen produced genome fragments,
translocation events can also be investigated:

```{bash, eval=FALSE}
dantools transloc alignment.sam
```

A list of DanTools functions can be produced with the simple command:

```{bash, eval=FALSE}
dantools
```

## License

[GPL-3.0](https://github.com/elsayed-lab/dantools/blob/master/LICENSE)
    
