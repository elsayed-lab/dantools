# Dantools

## Contact

Daniel Klimes (<daniel.s.klimes@gmail.com>)

## Overview

Dantools is a tool used for comparing genome sequences between
divergent organisms, such as different species. It does this by
creating a pseudogenome, a modification of a reference genome to match
the sequence of another. In doing so, a list of genetic variants
between given species is produced. Given a GFF file, these variants
can be annotated by their positions relative tofeatures. Additionally,
by keeping track of position shifts from indels, the annotations of
the original reference can be extended to the new species'
genome. This allows for the alignment of RNA-Seq data from divergent
species against a set of genomes with unified annotations. Critically,
Dantools does not require either genome to be well assembled, and can
even accept an RNA-Seq .fastq input.

## Methodology

Dantools in its base form accepts two genomes as input: the base and
the source. In its first step, Dantools breaks the source genome into
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

Because of the flexibility of inputs, Dantools can accept RNA-Seq
reads in fastq format and minimally-assembled genomes.

In future iterations of Dantools, functionality will be added for:
-Determining genome rearrangements
-Translating nucleic acid to amino acid changes

## Installation

### Prerequisites

One may wish to keep this in its own tree, perl provides a number of
mechanisms to do so, here is an example of how one might install
dantools in a self-contained tree with the github-downloaded source in
src/

```{bash, eval=FALSE}
mkdir -p dantools/202406
cd dantools/202406
git clone https://github.com/elsayed-lab/dantools.git
mv dantools src
## tell MakeMaker and Module::Build to use this tree
## All prerequisite perl libraries will go to dantools/202406/lib/perl5
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

Dantools relies on a set of software packages not listed in the
makefile to perform some of its tasks:

- [HISAT2](https://github.com/DaehwanKimLab/hisat2)
- [Freebayes](https://github.com/freebayes/freebayes)
- [samtools](https://github.com/samtools/samtools)
- [bcftools](https://github.com/samtools/bcftools)

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
dantools shift -v variants.vcf -f base.gff
```

Now, RNA-Seq data aligned against the pseudogenome can be counted with
the shifted GFF file.

If one wants to label the variants according to a GFF file:

```{bash, eval=FALSE}
dantools label -v variants.vcf -f base.gff --features five_prime_UTR,CDS,three_prime_UTR
```

These variants can optionally be translated into amino acid changes
and scored by a mutation scoring matrix.

Additional functionalities exist and are listed with the simple command:

```{bash, eval=FALSE}
dantools
```

## License

[GPL-3.0](https://github.com/elsayed-lab/dantools/blob/master/LICENSE)
    
