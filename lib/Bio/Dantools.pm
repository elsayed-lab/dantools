package Bio::Dantools;
use strict;
use autodie;
use diagnostics;
#use lib "$FindBin::Bin/../lib"; #Don't need this now that it's all one file
use warnings;

use FileHandle;
use File::Basename;
use File::Copy;
use File::Path qw"make_path rmtree";
use FindBin;
use List::Util qw"max min sum";
use List::MoreUtils qw"uniq";
use Parallel::ForkManager;
use POSIX qw"floor ceil";
use Text::CSV_XS::TSV;
use Data::Dumper;

use Bio::DB::Fasta;
use Bio::Matrix::IO;
use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::CodonTable;

## The Makefile.PL has a version_from stanza saying the version is in this file.
## I do not see it, so decided to make an executive decision:
my $VERSION = '202410';

# Set up my interrupt trap
BEGIN {
    $SIG{'INT'} = sub {
        exit(1);
    };
}

=head1 NAME

    Bio::Dantools - Tools to compare disparate genomes.

=head1 SYNOPSIS

    Bio::Dantools.pm provides functions which accept disparate inputs (high-throughput
    sequencing data, fasta genomes, gff annotations), ensure that they are ready for
    alignment, and passes them to hisat2 and freebayes for an iterative alignment
    process which results in new genome/gff files that rephrases the source genome
    in terms of the query input data.

    Each of the following commands produces a new pseudogenome rephrased by the base genome and
    a vcf file describing the differences between source and base:

    1.  If you have 2 fasta genomes:

    > dantools pseudogen -r reference.fasta -q query.fasta

    2.  If you have a base fastq genome and an already fragmented fastq genome:

    > dantools pseudogen -r reference.fasta -s query.fasta --fragment no

    3.  If you have unpaired RNA/DNA sequencing reads:

    > dantools pseudogen -r reference.fasta --reads-u source.fastq

    (note to self, make sure this will work with gz/bz2/xz compressed input)

    4.  Paired-end sequencing data:

    > dantools pseudogen -r reference.fasta -1 r1.fastq -2 r2.fastq

    (note to self, make sure this can happily handle different sequence orientations
     and pass the args along to hisat)

    Upon completion, a new gff file may be created which includes new positions for the features
    that represent any indels introduced in the pseudogen process:

    > dantools shift -v variants.vcf -f reference.gff

    (note to self, make sure this has sane defaults for the -v and -f args; I think in most
     cases dantools should be able to correctly guess them and/or figure them out from the
     pseudogen results)

    In addition, one may label observed variants via the input gff:

    > dantools label -v variants.vcf -f reference.gff --features five_prime_UTR,CDS

=head1 METHODS

=head2 C<pseudogen>

    This serves as the main pseudogen worker.  It uses the aligner
    and postprocess scripts to iterate the (re)creation of the
    pseudogenome until the variants observed by freebayes fall below the
    threshold defined by 'min_variants'.

=over

=item C<Arguments>

  Note to Daniel: I am totally bs-ing some of these based on my
  likely incorrect assumptions.

  min_variants: Threshold for falling out of the iterator loop --
    I am a little confused why it is used as '$var_count' in one
    place and directly in another; I need to look a little
    carefully there and see what I am missing.
  fai: Input fasta index file.
  base_idx: hisat2 index of the base genome; e.g. what we are moving
    this species/sample toward or what we are changing to reflect
    this species/sample, depending on how you think.
  base: Current working directory of the aligner over each iteration.
  gff: Genome Feature File describing the base genome.
  keepers: Set of file types to keep upon completion.
  output_name: I think basename of the output, I need to reread to be sure.
  threads: passed to hisat2, number of CPUs to dedicate to the aligner and/or freebayes.
  var_fraction: passed to the min-alternate-fraction argument of freebayes.
  outdir: Output directory for the toys.
  scoremin: Passed along to hisat2.
  bin_size: Passed along to Dantools::bin_maker()
  source: Used for input fasta data by hisat2.
  readsu: Used for input unpaired fastq reads by hisat2.
  reads1: Used for paired fastq reads by hisat2, R1.
  reads2: Ibid, R2.
  input_type: single word describing the various inputs one might use:
    fasta, fastq_u, fastq_p.
  fragment: Used to decide whether to invoke Bio::Dantools::fragment()
    on a large set of assembled chromosomes/contigs. (e.g. a
    genome fasta file).
  lengths: String representing the various lengths to fragment an
    input genome into.
  overlap: Define desired input fasta fragment overlap percentage.
  min_length: Passed along to the fragmenter.

=back

=cut

sub pseudogen {
    my %args = @_;
    chdir($args{'outdir'});

    my $it = 0;
    my $itp = -1;

    my $var_count = $args{'min_variants'} + 1;
    my $fai = $args{'fai'};
    my $reference_idx = $args{'reference_idx'};
    my $og_reference_idx = $reference_idx;
    my $reference = $args{'reference'}; #This changes throughout so it needs a variable
    my $og_reference = $reference;
    my $gff = $args{'gff'};
    my $keepers = $args{'keepers'};
    my $output_name = $args{'output_name'};
    my $threads = $args{'threads'};
    my $var_fraction = $args{'var_fraction'};
    my $var_depth = $args{'var_depth'};
    my $variant_caller = $args{'variant_caller'};
    my $continue = $args{'continue'};
    my $outdir = $args{'outdir'};
    my $min_variants = $args{'min_variants'};
    my $nocap = $args{'nocap'};
    my $scoremin = $args{'scoremin'};
    my $bin_size = $args{'bin_size'};
    my $query = $args{'query'};
    my $readsu = $args{'readsu'};
    my $reads1 = $args{'reads1'};
    my $reads2 = $args{'reads2'};
    my $input_type = $args{'input_type'};
    my $no_fragment = $args{'no_fragment'};
    my $lengths = $args{'lengths'};
    my $overlap = $args{'overlap'};
    my $min_length = $args{'min_length'};
    my $input_name;

    #Begin by fragmenting my reads if a fasta is given
    if (("$input_type" eq 'fasta') & (! $no_fragment)) {
        unless ($continue && -e 'fragments.fasta') {
        ## recommendation from atb: have everything return something
            Bio::Dantools::fragment(input => "$query",
                                    output => "fragments.fasta",
                                    lengths => "$lengths",
                                    overlap => "$overlap",
                                    min_length => "$min_length",
                                    logfile => "fragment_log.txt",
                                    threads => "$threads",
                                );
        }
    }

    #Creating my metadata file
    my $meta = FileHandle->new("> metadata.tsv");
    my $aln_scaled = 0;
    if (("$input_type" eq 'fasta')) {
        print $meta "Reference\tQuery\tIteration\tAlignment\tScale_Alignment\tVariants\n";
    }
    else {
        print $meta "Base\tSource\tIteration\tAlignment\tVariants\n";
    }

    ## recommendation from atb: perl is smrtr than bash and you need
    ## not protect yourself in this fashion.
    until ((("$var_count" < "$min_variants") | ( "$aln_scaled" > 100)) & ("$it" != 0)) {
        #First do some things different for iteration 0
        if (! -d "it${it}") {
            mkdir("it${it}");
        }

        if (! $it == 0) {
            unless ($continue && -s "it${it}/${output_name}_it${it}.fasta") {
                my $consensus = Bio::Dantools::consensus(input_fasta => $reference,
                                                         output_fasta => "it${it}/${output_name}_it${it}.fasta",
                                                         input_vcf => "it${itp}/variants.vcf");
            }

            unless ($continue && -s "it${it}/indexes/${output_name}_it${it}.1.ht2") {
                if (! -d "it${it}/indexes") {
                    mkdir("it${it}/indexes");
                }
                my $indexes = qx"hisat2-build -q -p ${threads} it${it}/${output_name}_it${it}.fasta it${it}/indexes/${output_name}_it${it}";
            }

            unless ($continue && -s "it${it}/${output_name}_it${it}.fasta.fai") {
                my $fai = qx"samtools faidx it${it}/${output_name}_it${it}.fasta";
            }
            $reference_idx = "it${it}/indexes/${output_name}_it${it}";
            $reference = "it${it}/$output_name" . "_it${it}.fasta";
            $fai = "$reference" . ".fai";

        } else {
            #Means this is iteration 0
            if ($fai eq '') {
                my $faidx = qx"samtools faidx -o original_reference.fasta.fai ${reference}";
                $fai = 'original_reference.fasta.fai';
            }

            if ($reference_idx eq '') {
                if (! -d "original_indexes") {
                    mkdir("original_indexes");
                }
                my $hisat_idx=qx"hisat2-build -q -f ${reference} original_indexes/original_indexes";
                $reference_idx = 'original_indexes/original_indexes';
            }
        }

        #Now onto alignment, which is the same for every iteration
        if (("$input_type" eq 'fasta') & (! $no_fragment)) {
            unless ($continue && -s "it${it}/MAPPING.stderr") {
                my $align=qx"hisat2 -x ${reference_idx} -f -p ${threads} --no-softclip --no-spliced-alignment --score-min ${scoremin} -k 1 --no-unal --norc -U fragments.fasta -S it${it}/raw.sam 2>it${it}/MAPPING.stderr";
            }
            $input_name = basename("$query", ('.fasta'));
        }  elsif (("$input_type" eq 'fasta') & ($no_fragment)) {
            unless ($continue && -s "it${it}/MAPPING.stderr") {
                my $align=qx"hisat2 -x ${reference_idx} -f -p ${threads} --no-softclip --no-spliced-alignment --score-min ${scoremin} -k 1 --no-unal --norc -U ${query} -S it${it}/raw.sam 2>it${it}/MAPPING.stderr";
            }
            $input_name = basename("$query", ('.fasta.'));
        } elsif (("$input_type" eq 'fastq_u')) {
            unless ($continue && -s "it${it}/MAPPING.stderr") {
                my $align=qx"hisat2 -x ${reference_idx} -q -p ${threads} --no-softclip --pen-canintronlen G,-8,7.5 --pen-noncanintronlen G,-8,7.5 --score-min ${scoremin} -k 1 --no-unal -U $readsu -S it${it}/raw.sam 2>it${it}/MAPPING.stderr";
            }
        $input_name = basename("$readsu", ('.fastq'));
        } elsif (("$input_type" eq 'fastq_p')) {
            unless ($continue && -s "it${it}/MAPPING.stderr") {
                my $align=qx"hisat2 -x ${reference_idx} -q -p ${threads} --no-softclip --pen-canintronlen G,-8,7.5 --pen-noncanintronlen G,-8,7.5 --score-min ${scoremin} -k 1 --no-unal -1 $reads1 -2 $reads2 -S it${it}/raw.sam 2>it${it}/MAPPING.stderr";
            }
            $input_name = basename("$reads1", ('.fastq'));
        }

        my $sort = qx"samtools sort -l 9 -O BAM -@ ${threads} -o it${it}/sorted.bam it${it}/raw.sam 2>&1";
        my $index = qx"samtools index it${it}/sorted.bam";
        if (-e "it${it}/raw.sam") {
            unlink("it${it}/raw.sam");
        }

        #Now to call variants
        unless ($continue && -s "it${it}/variants.bcf") {
            if ($variant_caller eq 'freebayes') {
                if ($threads > 1) {
                    my $regions=qx"fasta_generate_regions.py ${fai} 100000 > it${it}/freebayes_regions.tmp";
                    my $freebayes=qx"freebayes-parallel it${it}/freebayes_regions.tmp ${threads} -f ${reference} --min-alternate-fraction ${var_fraction} --min-alternate-count ${var_depth} it${it}/sorted.bam > it${it}/variants.vcf";
                } else {
                    my $freebayes=qx"freebayes -f ${reference} --min-alternate-fraction ${var_fraction} --min-alternate-count ${var_depth} it${it}/sorted.bam > it${it}/variants.vcf";
                }
            } elsif ($variant_caller eq 'bcftools') {
                my $filter = 'FORMAT/AD[0:1] >= ' . ${var_depth} . ' && (FORMAT/AD[0:1] / (FORMAT/AD[0:0] + FORMAT/AD[0:1])) >= ' . ${var_fraction};

                my $bcftools=qx"bcftools mpileup -a FORMAT/AD,FORMAT/DP -d 250 -f $reference it${it}/sorted.bam 2>/dev/null | bcftools call -O v -m -v --ploidy 2 | bcftools filter -i \"$filter\" > it${it}/variants.vcf "
            }
            my $view = qx"bcftools view -O bcf -o it${it}/variants.bcf it${it}/variants.vcf";
            my $bcf_idx = qx"bcftools index it${it}/variants.bcf";
        }

        #The above should have just produced everything I need to
        #continue through with my pipeline. That includes the
        #it"$it"/variants.vcf file, the MAPPING.stderr file, and the
        #new genome. The next step is to shift my bins accordingly
        if ("$it" == 0) {
            Bio::Dantools::bin_maker(
                input_fasta => "$og_reference",
                output => "it0/bins.tsv",
                bin_size => "$bin_size",
            );
        } else {
            Bio::Dantools::bin_shifter(
                input_vcf => "it$itp/variants.vcf",
                input_bins => "it$itp/bins.tsv",
                output => "it$it/bins.tsv",
            );
        }

        #Now the processing is done, I just need to do some
        #house-keeping and loop-maintaining

        #Counting up metadata
        my $aln_stderr = FileHandle->new("< it$it/MAPPING.stderr");
        my $aln;
      FIND_ALN: while (<$aln_stderr>) {
            if (index($_, "overall") != -1) {
                my @line = split(/ /, $_);
                $aln = $line[0];
                $aln =~ tr /%//d;
            }
        }

        my $vcf = FileHandle->new("< it$it/variants.vcf");
        $var_count = 0;
      COUNT_VARS: while (<$vcf>) {
            next COUNT_VARS if ($_ =~ /^#/);
            $var_count++;
        }
        if (("$input_type" eq 'fasta')) {
            $aln_scaled = $aln * 2;
            print $meta basename($og_reference, ('.fasta')) . "\t" . "$input_name\t" . "$it\t" . "$aln\t" . "$aln_scaled\t" . "$var_count\n";
            if ("$nocap" == 1) {
                $aln_scaled = 0;
            }
        }  else {
            print $meta basename($og_reference, ('.fasta')) . "\t" . "$input_name\t" . "$it\t" . "$aln\t" . "$var_count\n";
            $aln_scaled = 0;
        }

        #Determine what to delete, if anything
        if (("$it" != 0)) {
            eval {
                if ("$keepers" eq "none") {
                    rmtree("it$itp");
                } elsif ("$keepers" eq "genomes") {
                    unlink(("it$itp/variants.bcf.gz",
                            "it$itp/variants.bcf.gz.csi",
                            "it$itp/variants.vcf",
                            "it$itp/sorted.bam",
                            "it$itp/sorted.bam.bai",
                            "it$itp/bins.tsv"));
                } elsif ("$keepers" eq "alignments") {
                    unlink(("it$itp/variants.bcf.gz",
                            "it$itp/variants.bcf.gz.csi",
                            "it$itp/variants.vcf"));
                    #Genome files don't exist in it0
                    if ("$itp" != 0) {
                        unlink(("it$itp/$output_name" . "_it$itp.fasta",
                                "it$itp/$output_name" . "_it$itp.fasta.fai",
                                "it$itp/bins.tsv"));
                        rmtree("it$itp/indexes");
                    }
                } elsif ("$keepers" eq "variants") {
                    unlink(("it$itp/sorted.bam",
                            "it$itp/sorted.bam.bai"));
                    #Genome files don't exist in it0
                    if ("$itp" != 0) {
                        unlink(("it$itp/$output_name" . "_it$itp.fasta",
                                "it$itp/$output_name" . "_it$itp.fasta.fai",
                                "it$itp/bins.tsv"));
                        rmtree("it$itp/indexes");
                    }
                }
            };
            if ($@) {
                print "Caught error: $@";
            }
        } else {
            copy("it0/bins.tsv", "original_bins.tsv");
        }

        $it++;
        $itp++;
    }

    #Ok and now I need to put the important things together in an
    #output directory
    if (! -d 'output') {
        if ( -e 'output') {
            die "output exists as file, can't overwrite";
        }
        else {
            mkdir 'output';
            mkdir 'output/indexes';
        }
    } elsif (! -d 'output/indexes') {
        if (-e 'output/indexes') {
            die "output/indexes exists as a file, can't overwrite";
        }
        else {
            mkdir('output/indexes');
        }
    }
    copy("$reference", "output/$output_name.fasta");
    copy("it$itp/bins.tsv", "output/bins.tsv");
    copy("it$itp/sorted.bam", "output/alignment.bam");

    #Redo the "catching stuff to delete" based on keepers:
    if ("$keepers" eq "none") {
        rmtree("it$itp");
    }
    elsif ("$keepers" eq "genomes") {
        unlink(("it$itp/variants.bcf.gz",
                "it$itp/variants.bcf.gz.csi",
                "it$itp/variants.vcf",
                "it$itp/sorted.bam",
                "it$itp/sorted.bam.bai",
                "it$itp/bins.tsv"));
    }
    elsif ("$keepers" eq "alignments") {
        unlink(("it$itp/variants.bcf.gz",
                "it$itp/variants.bcf.gz.csi",
                "it$itp/variants.vcf"));
        #Genome files don't exist in it0
        if ("$itp" != 0) {
            unlink(("it$itp/$output_name" . "_it$itp.fasta",
                    "it$itp/$output_name" . "_it$itp.fasta.fai",
                    "it$itp/bins.tsv"));
            rmtree("it$itp/indexes");
        }
    }
    elsif ("$keepers" eq "variants") {
        unlink(("it$itp/sorted.bam",
                "it$itp/sorted.bam.bai"));
        #Genome files don't exist in it0
        if ("$itp" != 0) {
            unlink(("it$itp/$output_name" . "_it$itp.fasta",
                    "it$itp/$output_name" . "_it$itp.fasta.fai",
                    "it$itp/bins.tsv"));
            rmtree("it$itp/indexes");
        }
    }

    #The below is some setup for my needler, which I may just in its
    #own subroutine

    #First I need to break up each genome into its contigs
    Bio::Dantools::contig_split(
        fasta => "$og_reference",
        outdir => "ref",
    );
    Bio::Dantools::contig_split(
        fasta => "output/$output_name.fasta",
        outdir => "alt",
    );

    my @contigs;
    my @breaks;
    my @ref_ends;
    my @alt_ends;
    my $prev_contig = '_balrog_'; #Set it to something that won't be matched
    ## (you shall not pass, flame of udun)
    my $pos=0;
    my $fh = FileHandle->new("original_bins.tsv");
    while (<$fh>) {
        chomp;
        my @line = split(/\t/, $_);
        push @contigs, $line[0];
        push @ref_ends, $line[1];
        if ("$line[0]" ne "$prev_contig") {
            push @breaks, $pos;
        }
        $prev_contig = $line[0];
        $pos++;
    }
    close($fh);
    $fh = FileHandle->new("output/bins.tsv");
    while (<$fh>) {
        chomp;
        my @line = split(/\t/, $_);
        push @alt_ends, $line[1];
    }

    my $hisat_idx = qx"hisat2-build -q -p ${threads} output/${output_name}.fasta output/indexes/${output_name}";

    unless ($continue && -s "output/strict_aln.stderr") {
            if (("$input_type" eq 'fasta') & (! $no_fragment)) {
                my $align=qx"hisat2 -x output/indexes/${output_name} -f -p ${threads} --no-softclip --no-spliced-alignment -k 1 --no-unal --norc -U fragments.fasta -S output/strict_aln.sam 2>output/strict_aln.stderr";
            } elsif (("$input_type" eq 'fasta') & ($no_fragment)) {
                my $align=qx"hisat2 -x output/indexes/${output_name} -f -p ${threads} --no-softclip --no-spliced-alignment -k 1 --no-unal --norc -U ${query} -S output/strict_aln.sam 2>output/stric_aln.stderr";
            } elsif (("$input_type" eq 'fastq_u')) {
                my $align=qx"hisat2 -x output/indexes/${output_name} -q -p ${threads} --no-softclip -k 1 --no-unal -U $readsu -S output/strict_aln.sam 2>output/strict_aln.stderr";
            } elsif (("$input_type" eq 'fastq_p')) {
                my $align=qx"hisat2 -x output/indexes/${output_name} -q -p ${threads} --no-softclip -k 1 --no-unal -1 $reads1 -2 $reads2 -S output/strict_aln.sam 2>output/strict_aln.stderr";
            }
        }

    my $sort = qx"samtools sort -l 9 -O BAM -@ ${threads} -o output/strict_aln.bam output/strict_aln.sam";
    my $index = qx"samtools index output/strict_aln.bam";
    unlink("output/strict_aln.sam");
    my $depth = qx"samtools depth -aa -@ ${threads} -o output/depth.tsv output/strict_aln.bam";

    if (! -d "needle_files") {
        mkdir('needle_files');
    }

    my $pm = Parallel::ForkManager->new(${threads});
  ALNS: for my $idx (0..$#ref_ends) {
        my $pid = $pm->start and next ALNS;
        if (grep { $_ eq $idx } @breaks) {
            #This means I am at the "start" of a contig
            my $stretch = qx"stretcher -asequence alt/$contigs[$idx] -bsequence ref/$contigs[$idx] -sbegin1 1 -send1 $alt_ends[$idx] -sbegin2 1 -send2 $ref_ends[$idx] -aformat3 fasta -sid1 alt -sid2 ref -outfile needle_files/${idx} -gapopen 50 -gapextend 15 2>>needle_files/errors.txt";
        } else {
            my $start1 = $alt_ends[$idx - 1] + 1;
            my $start2 = $ref_ends[$idx - 1] + 1;

            my $stretch = qx"stretcher -asequence alt/$contigs[$idx] -bsequence ref/$contigs[$idx] -sbegin1 $start1 -send1 $alt_ends[$idx] -sbegin2 $start2 -send2 $ref_ends[$idx] -aformat3 fasta -sid1 alt -sid2 ref -outfile needle_files/${idx} -gapopen 50 -gapextend 15 2>>needle_files/errors.txt";
        }
        $pm->finish;
    }
    $pm->wait_all_children;

    #Now there should be a bunch of alignment files which I can go through

    my $file_idx = 0;

    my $out = FileHandle->new("> needle_files/combined_bases.tsv");
    while ("$file_idx" < scalar(@contigs)) {
        #Note in my original I use Bio::AlignIO, maybe this won't work
        my $raw = Bio::SeqIO->new(
            -file => "needle_files/$file_idx",
            -format => 'Fasta');
        my $refseq;
        my @refseq;
        my $altseq;
        my @altseq;
        while (my $seq = $raw->next_seq) {
            if ($seq->id eq 'ref') {
                $refseq = $seq->seq;
                @refseq = split //, $refseq;
            }
            elsif ($seq->id eq 'alt') {
                $altseq = $seq->seq;
                @altseq = split //, $altseq;
            }
        }

        my $idx = 0;
        my $lim = scalar @refseq;

        while ($idx < $lim) {
            print $out "$contigs[$file_idx]\t$refseq[$idx]\t$altseq[$idx]\n";
            $idx++;
        }
        $file_idx++;
    }
    close($out);

    # Now I have a list of every base from each needle alignment next
    # to what it was in the alignment. I could and probably should make
    # this more efficient by skipping the printing to a file and just
    # building a hash of arrays.
    Bio::Dantools::vcf_maker(
        input => "needle_files/combined_bases.tsv",
        output => "output/$output_name.vcf",
        tmpdir => "output",
        depths => "output/depth.tsv",
        sample => "${output_name}",
    );
    if ("$keepers" ne "all") {
        rmtree(['ref', 'alt', 'needle_files', 'original_indexes'], 0, 0);
        for my $file (('original_bins.tsv', 'original_reference.fasta.fai', 'fragments.fasta', 'fragment_log.txt')) {
            if (-e "$file") {
                unlink("$file");
            }
        }
    }

    #This used to automatically generate a new GFF and label the vcf,
    #but I don't feel that is necessary for the general pseudogen. If
    #people want to, they can always call dantools label and dantools
    #shift separately
}
;

#The below function takes a VCF and FASTA file and applies variants
#from the former to the latter. Importantly for me, it does not "skip"
#variants in overlap, but sequentially applies them. This is important
#because it keeps the bin-shifting logic correct
sub consensus {
    my %args = @_;
    my $seqio_in = Bio::SeqIO->new(
        -file => $args{'input_fasta'},
        -format => 'Fasta');
    my $seqio_out = Bio::SeqIO->new(
        -format => 'Fasta',
        -file => "> $args{'output_fasta'}");

    #Read in my VCF:
    my $input_vcf = $args{'input_vcf'};
    my $vcf_fh = FileHandle->new($input_vcf);
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1,});
    my %vcf_hash;

    open($vcf_fh, "<:encoding(utf8)", "$input_vcf") or die "Can't open vcf_fh";

    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'metadata');

    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if ($row->{'contig'} =~ /^#/);
        my $alt = (split(/\,/, $row->{'alternate'}))[0];

        my %hash = (
            pos => $row->{'position'},
            ref => $row->{'reference'},
            alt => $alt
        );
        if (! exists($vcf_hash{$row->{'contig'}})) {
            @vcf_hash{$row->{'contig'}} = ();
        }
        push(@{ $vcf_hash{$row->{'contig'}} }, \%hash);
    }

  CONTIGS: while (my $chrom = $seqio_in->next_seq) {
        #Now to add all of the metadata
        my @sub = sort { $a->{'pos'} <=> $b->{'pos'} } @{ $vcf_hash{ $chrom->{'primary_id'} } };
        my $idx = 0;
        my $contig_length = $chrom->{'primary_seq'}->{'length'};

        my $shift = 0; #how much to adjust positions to account for indels
      VARS: while (my $entry = $sub[$idx]) {
            my $rlength = length($entry->{'ref'});
            substr($chrom->{'primary_seq'}->{'seq'}, $entry->{'pos'} - 1 + $shift, $rlength) = $entry->{'alt'};
            $shift += length($entry->{'alt'}) - length($entry->{'ref'});
            $idx++;
        }

        $seqio_out->write_seq($chrom);
    }

}

#The below can break my genome into pieces
sub fragment {
    my %args = @_;
    my $seqio_out;
    if ($args{'output'} ne 'NO_OUTPUT_PROVIDED') {
        $seqio_out = Bio::SeqIO->new(-format => 'fasta', -file => "> $args{'output'}");
    }
    else {
        $seqio_out = Bio::SeqIO->new(-format => 'fasta', -fh => \*STDOUT);
    }

    my $seqio_in = Bio::SeqIO->new(-format => 'fasta', -file => $args{'input'});

    my @sizes = split(/\,/, $args{'lengths'});
    my $log = FileHandle->new("> $args{'logfile'}");
    my $num_written = 0;
    my $pm = Parallel::ForkManager->new($args{'threads'});

  SEQS: while (my $seq = $seqio_in->next_seq) {
        my $pid = $pm->start and next SEQS;
        my $rev = $seq->revcom;
        my $start = 1;
        my $len = $seq->length;
        my $id = $seq->id;

        for my $size (@sizes) {
            my $chunk_start = 1;
            my $chunk = 1;
            my $increment = ceil($size * ((100 - $args{'overlap'}) / 100));

          CHUNKS: while ($chunk_start < $len) {
                my $chunk_end = $chunk_start + $size - 1;
                $chunk_end = $len if ($chunk_end > $len);
                last CHUNKS if ($chunk_start >= $len);

                my $subseq=$seq->subseq($chunk_start, $chunk_end);
                my $revcom_subseq = $rev->subseq($chunk_start, $chunk_end);
                my $chunk_seq = Bio::Seq->new(-seq => $subseq,
                                              -display_id => qq"${id}_fwd_${size}_${chunk_start}_${chunk_end}");
                if ($chunk_seq->length >= $args{'min_length'}) {
                    my $written = $seqio_out->write_seq($chunk_seq);
                }

                my $rev_chunk_seq = Bio::Seq->new(-seq => $revcom_subseq,
                                                  -display_id => qq"${id}_rev_${size}_${chunk_start}_${chunk_end}");
                if ($rev_chunk_seq->length >= $args{'min_length'}) {
                    my $rev_written =$seqio_out->write_seq($rev_chunk_seq);
                }
                $num_written++;
                $chunk_start = ($chunk_start + $increment);
            }                   #End iterating over every chunk
            print $log "Wrote $num_written fragments of size $size for $id. \n";
        }                       #End iterating over every size
        $pm->finish;
    }
    $log->close;
    $pm->wait_all_children;
}

#The below is used in iteration 0 and printed to it0. It is a very
#simple subroutine which will take an input fasta sequence, bin size,
#and a desired output and then defines those bins.
sub bin_maker {
    my %args = @_;
    my $fasta = Bio::SeqIO->new(-file => "$args{'input_fasta'}",
                                -format => 'fasta');
    my $out = FileHandle->new("> $args{'output'}");
    my $seq;
  SEQS: while($seq = $fasta->next_seq()) {
        my $length = $seq->length;
        my $end = $args{'bin_size'};
        my $contig = $seq->primary_id;
      BINS: while(1) {
            if ($end >= $length) {
                print $out "$contig" . "\t" . "$length\n";
                last BINS;
            }
            else {
                print $out "$contig" . "\t" . "$end\n";
            }
            $end = $end + $args{'bin_size'};
        }
    }
}

#This takes any input bins of the format contig \t position and will
#shift the bins' position according to indels in the vcf file. It
#takes an input vcf, an output file, input bins,
sub bin_shifter {
    my %args = @_;
    my $vcf_fh = FileHandle->new("$args{'input_vcf'}");
    my @vcf;
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($vcf_fh, "<:encoding(utf8)", "$args{'input_vcf'}") or die "Can't open vcf_fh";

    my $contig = '';
    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'metadata');

    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if $row->{'contig'} =~ /^#/;
        my $shift = length((split(/\,/, $row->{'alternate'}))[0]) - length($row->{'reference'});

        next if $shift == 0;

        my %hash = (
            seq_id => $row->{'contig'},
            pos => $row->{'position'},
            ref => $row->{'reference'},
            alt => $row->{'alternate'},
            shift => $shift
        );
        push(@vcf, \%hash);
    }
    close($vcf_fh);
    #print Dumper(@vcf);

    #Importantly I don't sort this. That's because if I do, it's going
    #to sort it in a different way from my original. In dantools
    #pseudogen, the vcf and bins should always be sorted the same anyway
    $contig = '';
    my @contigs;
    my $shiftsum;
    my $vcf_index = 0;
    my $shift = 0;
    my $vcf_row;
    my $prev_row;
    my $vcf_contig;
    my $bins_fh = FileHandle->new("< $args{'input_bins'}");
    my $bins_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($bins_fh, "<:encoding(utf8)", "$args{'input_bins'}") or die "Can't open bins_fh";
    $bins_tsv->column_names('contig', 'end');
    my $out = FileHandle->new("> $args{'output'}");
    my @bins_db;
    if (scalar @vcf == 0) {
        copy("$args{'input_bins'}", "$args{'output'}");
    }
    else {
      BINS: while(my $entry = $bins_tsv->getline_hr($bins_fh)) {
            if ($contig ne $entry->{'contig'}) {
                $shiftsum = 0;
                $contig = $entry->{'contig'};
                push @contigs, $contig;
            }
            my $end = $entry->{'end'};

          SHIFT: while ($vcf_row = $vcf["$vcf_index"]) {
                $vcf_contig = $vcf_row->{'seq_id'};
                if ($contig ne $vcf_contig) {
                    if (! grep /$vcf_contig/, @contigs) {
                        #Means this variant is on next contig;
                        last SHIFT;
                    } else {
                        #Variant is still on previous contig
                        $vcf_index++;
                        $shiftsum = 0;
                        next SHIFT;
                    }
                }
                if ($end < $vcf_row->{'pos'}) { last SHIFT };

                $shift = $vcf_row->{'shift'};
                $shiftsum = $shiftsum + $shift;
                $vcf_index++;
            }
            #Some logic to sheck if deletion covers the feature
            $prev_row = $vcf[$vcf_index - 1];
            my $prev_shift = $prev_row->{'shift'};
            if (($prev_row->{'pos'} - $prev_shift > $end) && ($prev_row->{'seq_id'} eq $contig) && ($vcf_index != 0)) {
                print $out "$contig" . "\t" . ($prev_row->{'pos'} + $shiftsum - $prev_shift) . "\n";
            }
            else {
                my $new_out = $end + $shiftsum;
                print $out "$contig" . "\t" . "$new_out\n";
            }
        }
    }
}

#This function splits a given fasta file into individual files named
#by the seq_id of each. This is important for the needler

sub contig_split {
    my %args = @_;
    my $input_fasta = $args{'fasta'};
    my $outdir = $args{'outdir'};

    mkdir $outdir  if (defined($outdir) && ! -d $outdir);
    my $fasta = Bio::SeqIO->new(-file => "$input_fasta", -format => 'Fasta');
    my $file;

  SEQS: while (my $seq = $fasta->next_seq()) {
        my $contig = $seq->primary_id;
        my $name = "$contig";

        $file = "$outdir/$contig";
        my $out = Bio::SeqIO->new(-file => "> $file", -format => 'Fasta');

        $out->write_seq($seq);
    }
}

#The next subroutine seeks to take a combined_bases.tsv file and
#convert it into a vcf file

sub vcf_maker {
    my %args = @_;
    my $input = $args{'input'};
    my $output = $args{'output'};
    my $tmpdir = $args{'tmpdir'};
    my $depthfile = $args{'depths'};
    my $sample = $args{'sample'};

    if (! -d "$tmpdir") {
        mkdir("$tmpdir") if (! -d $tmpdir);
    }
    ;

    my $vcf = FileHandle->new("> $tmpdir/tmp_variants.vcf");
    my $headfile = FileHandle->new("> $output");

    #I need to read in the data from my depths file, which should
    #follow base-by-base the "alt" column of my
    #combined_bases.tsv. I've checked this and this is true
    my @depths;
    my $depth_fh = FileHandle->new("$depthfile");
    my $depth_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($depth_fh, "<:encoding(utf8)", "$depthfile");
    $depth_tsv->column_names('chromosome', 'pos', 'depth');
    while (my $row = $depth_tsv->getline_hr($depth_fh)) {
        push(@depths, $row->{'depth'});
    }
    ;
    #Note that I don't have any real metadata yet, but I may in the future
    my $header = "##fileformat=VCFv4.2\n";
    print $headfile $header;

    my $raw_fh = FileHandle->new("$input");
    my $raw_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($raw_fh, "<:encoding(utf8)", "$input") or die "Can't open input file";
    $raw_tsv->column_names('chromosome', 'ref', 'alt');

    #The things I'll be building
    my $ref = '';               #String to build ref column of vcf
    my $alt = '';               #String to build alt column of vcf
    my $chrom = '_balrog_';
    my $relpos = 1;     #track position on reference chromosome/contig
    #The below depthidx is how I track depth. Note it starts at -1
    #because when I print something I'm always "on the next
    #loop". Yes, I've tested this
    my $depthidx = -1;
    my $string;
    my $pos;
    my $length;
  VARS: while(my $row = $raw_tsv->getline_hr($raw_fh)) {
        #Check if the chromosome has changed
        next if(! defined($row->{'alt'}));
        if ($chrom ne $row->{'chromosome'}) {
            if ($chrom ne '_balrog_') {
                $length = $relpos - 1;
                print $headfile "##contig=<ID=" . "$chrom" . ",length=" . "$length" . ">\n";
            }
            ;
            $relpos = 1;        #Reset relative chromosome position
            $chrom = $row->{'chromosome'};
        }

        my $rowref = $row->{'ref'};
        my $rowalt = $row->{'alt'};
        if ($rowref eq $rowalt) {
            #I need to adjust using this catch-all
            if ((length($ref) != 0) & ("$ref" ne "$alt")) {
                $pos = $relpos - $ref =~ tr/-//c;
                $ref =~ tr/-//d;
                $alt =~ tr/-//d;
                #I need to double check after removing '-'
                if ("$ref" ne "$alt") {
                    print $vcf
                        $chrom, "\t",
                        $pos, "\t",
                        '.', "\t",
                        $ref, "\t",
                        $alt, "\t",
                        '.', "\t",
                        '.', "\t",
                        "DP=" . $depths[$depthidx], "\t",
                        "GT:DP", "\t",
                        "1/1:" . $depths[$depthidx], "\n";
                }
            }
            $ref = "$rowref";
            $alt = "$rowalt";
            $relpos++;
            $depthidx++;
            next VARS;
        }

        if ($rowref eq '-') {
            $ref = "$ref" . "$rowref";
            $alt = "$alt" . "$rowalt";
            $depthidx++;
            next VARS;
        }
        elsif ($rowalt eq '-') {
            $ref = "$ref" . "$rowref";
            $alt = "$alt" . "$rowalt";
            $relpos++;
            next VARS;
        }
        #Make the logic increment depth properly

        #If neither of the above were evaluated, this must be an SNP:
        if ((length($ref) != 0) & ("$ref" ne "$alt")) {
            #I'm sure there is a more efficient way of doing all
            #this. Right now I'm checking that ref ne alt twice. This is
            #because for shifting, ref and alt need to contain - for
            #indel. However, after deleting them the ref and alt may by
            #chance become the same, requiring me to check again
            $pos = $relpos - $ref =~ tr/-//c;
            $ref =~ tr/-//d;
            $alt =~ tr/-//d;
            if ("$ref" ne "$alt") {
                print $vcf
                    $chrom, "\t",
                    $pos, "\t",
                    '.', "\t",
                    $ref, "\t",
                    $alt, "\t",
                    '.', "\t",
                    '.', "\t",
                    "DP=" . $depths[$depthidx], "\t",
                    "GT:DP", "\t",
                    "1/1:" . $depths[$depthidx], "\n";
            }
        }
        $ref = $rowref;
        $alt = $rowalt;
        $relpos++;
        $depthidx++;
        next VARS;
    }
    close($vcf);
    #Print last contig that wasn't caught:
    $length = $relpos - 1;
    print $headfile "##contig=<ID=" . "$chrom" . ",length=" . "$length" . ">\n";
    print $headfile '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">' . "\n";
    print $headfile '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "\n";
    print $headfile '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">', "\n";

    print $headfile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample}\n";

    #Put the things from my vcf file into my final product
    my $tmp_vcf = FileHandle->new("< $tmpdir/tmp_variants.vcf");
    while (<$tmp_vcf>) {
        print $headfile $_;
    }
    ;
    close($headfile);
    unlink("$tmpdir/tmp_variants.vcf");
}

#The below function seeks to take an input sam file and give a sense
#of which regions in the alignment fragments equal which in the
#reference it was aligned to. By its nature, this would only work with
#alignments from dantools fragment-produced fragments, but it should
#be cool (in my opinion).

#First define a baby function that, given a start position and CIGAR
#string, gives the end position of the alignment
sub calc_end_pos {
    my ($start_pos, $cigar) = @_;
    my $end_pos = $start_pos;

    while ($cigar =~ /(\d+)([MIDNSHPX=])/g) {
        my ($length, $op) = ($1, $2);
        if ($op eq 'M' || $op eq 'D' || $op eq 'N') {
            $end_pos += $length;
        }
    }
    return $end_pos;
}

#Another baby function which, given an array of hashes for one contig,
#parses which to keep and how to do so
sub parse_outputs {
    my ($array_ref, $min_depth, $min_qbase, $min_length, $max_diff, $all) = @_;
    my @contig_hashes = @$array_ref;
    my @output;
    @contig_hashes = grep {$_->{'depth'} > $min_depth &&
                               $_->{'qbases'} > $min_qbase &&
                               $_->{'rlength'} > $min_length &&
                               abs(($_->{'rlength'} - $_->{'qlength'}) / $_->{'rlength'}) < $max_diff
                           } @contig_hashes;

    @contig_hashes = sort {
        $a->{'rcontig'} cmp $b->{'rcontig'} ||
        $a->{'rstart'} <=> $b->{'rstart'}
    } @contig_hashes;

    #Now iterate through each and consider special cases
    my $idx = 0;
    my $idx_limit = scalar(@contig_hashes);
    #First I need to flip each feature with inverted directionality
    #(happens) with minus strand alignments


    if ($all) {
        return (@contig_hashes);
    }

  SUBSEQS: while ($idx < $idx_limit) {
        #I need to resort and index every iteration
        $idx_limit = scalar(@contig_hashes);
        @contig_hashes = sort {
            $a->{'rcontig'} cmp $b->{'rcontig'} ||
                $a->{'rstart'} <=> $b->{'rstart'}
            } @contig_hashes;

        my $first = $contig_hashes[$idx];
        my $second = $contig_hashes[$idx + 1];

        if (! defined($second)) {
            push(@output, $first);
            $idx++;
        } elsif ($first->{'rend'} < $second->{'rstart'}) {
            #means they don't overlap
            push(@output, $first);
            $idx++;
        } elsif ($first->{'rend'} < $second->{'rend'}) {
            #Overlapping but not fully, take weighted mean based on depth
            my $overlap = $first->{'rend'} - $second->{'rstart'} + 1;

            my $first_weight = $first->{'depth'} / ($first->{'depth'} + $second->{'depth'});
            my $second_weight = $second->{'depth'} / ($first->{'depth'} + $second->{'depth'});

            my $first_scale = int(0.5 + $first_weight * $overlap);
            my $second_scale = $overlap - $first_scale;

            $first->{'rend'} = $first->{'rend'} - $first_scale;
            $second->{'rstart'} = $second->{'rstart'} + $second_scale;

            $first->{'rlength'} = $first->{'rend'} - $first->{'rstart'} + 1;
            if ($first->{'rlength'} > 0) {
                $first->{'depth'} = $first->{'qbases'} / $first->{'rlength'};
            } else {
                splice(@contig_hashes, $idx, 1);
            }

            $second->{'rlength'} = $second->{'rend'} - $second->{'rstart'} + 1;
            if ($second->{'rlength'} > 0) {
                $second->{'depth'} = $second->{'qbases'} / $second->{'rlength'};
            } else {
                splice(@contig_hashes, $idx, 1);
            }
        } else {
            #In this case, the second feature is nested
            #within the first, so I make a depth-based
            #determination of what to do. First though,
            #I need to extract ALL contained elements
            if ($first->{'depth'} > $second->{'depth'}) {
                #Remove second subseq
                splice(@contig_hashes, $idx + 1, 1);
            } else {
                #Nest second subseq within first, so I need to add a
                #feature. I also need to update the query positions
                #(via estimate).
                my %third = %{$first};

                $first->{'rend'} = $second->{'rstart'} - 1;
                $third{'rstart'} = $second->{'rend'} + 1;

                #Estimate depth from relative length
                $first->{'rlength'} = $first->{'rend'} - $first->{'rstart'} + 1;
                $third{'rlength'} = $third{'rend'} - $third{'rstart'} + 1;
                #If the shifting inverted the feature (negative
                #length)
                my $total_rlength = $first->{'rlength'} + $third{'rlength'};
                if ($first->{'rlength'} < 1) {
                    splice(@contig_hashes, $idx, 1);
                } else {
                    $first->{'depth'} = $first->{'qbases'} / $first->{'rlength'};
                    $first->{'qbases'} = int($first->{'qbases'} * ($first->{'rlength'} / $total_rlength))
                }

                if ($third{'rlength'} > 0) {
                    $third{'depth'} = int($third{'qbases'} * ($third{'rlength'} / $total_rlength));
                    $third{'depth'} = $third{'qbases'} / $third{'rlength'};
                    push(@contig_hashes, \%third);
                }

            }
        }
    }
    return @output;
}

#My final baby function. Takes an input transloc array of hashes and a
#VCF file. Shifts the translocations relative to variants by changing
#REF positions.

sub transloc_shifter {
    my ($array_ref, $vcf_file) = @_;
    my @contig_hashes = @$array_ref;

    #Read in my VCF
    my $vcf_fh = FileHandle->new($vcf_file);
    my @vcf;
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($vcf_fh, "<:encoding(utf8)", $vcf_file) or die "Could not open VCF $vcf_file\n";
    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'info');
    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if ($row->{'contig'} =~ /^#/);
        next if (! defined($row->{'info'}));
        my $shift = length((split(/\,/, $row->{'alternate'}))[0]) - length($row->{'reference'});
        next if $shift == 0;
        my %hash = (
            seq_id => $row->{'contig'},
            pos => $row->{'position'},
            strand => $row->{'strand'},
            ref => $row->{'ref'},
            alt => $row->{'alt'},
            shift => $shift
        );

        push(@vcf, \%hash);
    }
    close($vcf_fh);
    if (scalar(@vcf) == 0) {
        print STDERR "WARNING: No indels in VCF file, not adjusting output";
        return @contig_hashes;
    }
    @vcf = sort {
        $a->{'seq_id'} cmp $b->{'seq_id'} || $a->{'pos'} <=> $b->{'pos'}
    } @vcf;

    #Index my features to maintain sort order
    my $idx = 0;
    for my $ref (@contig_hashes) {
        $ref->{'idx'} = $idx;
        $idx++;
    }

    #Because I can't assume no overlap (may have passed --all), I need
    #to sort starts and ends separately
    my @starts = @contig_hashes;
    my @ends = @starts;
    @starts = sort {
        $a->{'rcontig'} cmp $b->{'rcontig'} || $a->{'rstart'} <=> $b->{'rstart'}
    } @starts;
    @ends = sort {
        $a->{'rcontig'} cmp $b->{'rcontig'} || $a->{'rstart'} <=> $b->{'rstart'}
    } @ends;

    my @shifted_starts;
    my @shifted_ends;

    my @seen_contigs;
    my $entry;
    my $entry_start;
    my $contig;
    my $entry_idx = 0;
    my $vcf_idx = 0;
    my $shiftsum = 0;
    my $shift = 0;
    my $varpos = 0;
    my $var_contig = '';
    my $vcf_row;
  STARTS: while ($entry = $starts[$entry_idx]) {
        $contig = $entry->{'rcontig'};
        $entry_start = $entry->{'rstart'};
        if (! grep /$contig/, @seen_contigs) {
            $shiftsum = 0;
            push(@seen_contigs, $contig);
            $shift = 0;
        }
      START_SHIFT: while ($vcf_row = $vcf[$vcf_idx]) {
            my $tmp_contig = $vcf_row->{'seq_id'};
            if ($contig ne $tmp_contig) {
                if (! grep /$tmp_contig/, @seen_contigs) {
                    last START_SHIFT;
                } else {
                    $vcf_idx++;
                    next START_SHIFT;
                }
            }
            if ($entry_start < $vcf_row->{'pos'} + $shiftsum) {
                #Note how I added $shiftsum. This is to keep track of
                #where I am on ALT, since that is what my entries are
                #relative to (different from GFF shifter)
                last START_SHIFT;
            }
            $shift = $vcf_row->{'shift'};
            $varpos = $vcf_row->{'pos'} + $shiftsum;
            $shiftsum = $shiftsum + $shift;
            $var_contig = $vcf_row->{'seq_id'};
            $vcf_idx++;
        }
        if (($varpos - $shift > $entry_start ) && ($var_contig eq $contig) && ($idx != 0)) {
            $entry->{'rstart'} = $varpos - $shiftsum + $shift;
        } else {
            $entry->{'rstart'} = $entry_start - $shiftsum; #subtract because I'm backshifting
        }
        push(@shifted_starts, $entry);
        $entry_idx++;
    }
    #Repeat for end positions
    $shiftsum = 0;
    $shift = 0;
    $varpos = 0;
    $var_contig = '';
    @seen_contigs = '';
    $entry_idx = 0;
    $vcf_idx = 0;
    my $entry_end;
  ENDS: while ($entry = $ends[$entry_idx]) {
        $contig = $entry->{'rcontig'};
        $entry_end = $entry->{'rend'};
        if (! grep /$contig/, @seen_contigs) {
            $shiftsum = 0;
            push @seen_contigs, $contig;
            $shift = 0;
        }
      END_SHIFT: while ($vcf_row = $vcf[$vcf_idx]) {
            my $tmp_contig = $vcf_row->{'seq_id'};
            if ($contig ne $tmp_contig) {
                if (! grep /$tmp_contig/, @seen_contigs) {
                    #Means this variant is on the next contig
                    last END_SHIFT;
                }
                else {
                    #Means variant is on previous chromosome
                    $vcf_idx++;
                    next END_SHIFT;
                }
            }
            if ($entry_end < $vcf_row->{'pos'} + $shiftsum) {
                last END_SHIFT;
            }
            #If we get here, it means the variant is on the right
            #chromosome and upstream of our feature
            $shift = $vcf_row->{'shift'};
            $varpos = $vcf_row->{'pos'} + $shiftsum; #kept for later logic
            $shiftsum = $shiftsum + $shift;
            $var_contig = $vcf_row->{'seq_id'};
            $vcf_idx++;
        }
        #Now I should have collected the proper shiftsum:
        if (($varpos - $shift > $entry_end) & ($var_contig eq $contig) & ($idx != 0)) {
            $entry->{'rend'} = $varpos - $shiftsum + $shift;
        }
        else {
            $entry->{'rend'} = $entry_end - $shiftsum; #subtract because I'm backshifting
        }
        push(@shifted_ends, $entry);
        $entry_idx++;
    }
    @shifted_starts = sort {
        $a->{'idx'} <=> $b->{'idx'}
    } @shifted_starts;
    @shifted_ends = sort {
        $a->{'idx'} <=> $b->{'idx'}
    } @shifted_ends;

    my @output;
    my $idx_limit = scalar(@shifted_starts);
    $idx = 0;
    while ($idx < $idx_limit) {
        my $tmp = $shifted_starts[$idx];
        $tmp->{'rend'} = $shifted_ends[$idx]->{'rend'};
        push(@output, $tmp);
        $idx++;
    }
    return @output;
}

#My final subroutine for this translocation calculator will seek to
#merge to input hashes:

sub merge_hashes {
     my ($array_ref) = @_;
    my @hashes = @$array_ref;
    my $sum_aln = 0;
     my $sum_qbases = 0;
     my $rstart = $hashes[0]->{'rstart'};
     my $rend = $hashes[0]->{'rend'};
     my $qstart = $hashes[0]->{'qstart'};
     my $qend = $hashes[0]->{'qend'};
     for my $ref (@hashes) {
         $rstart = $ref->{'rstart'} if ($ref->{'rstart'} < $rstart);
         $rend = $ref->{'rend'} if ($ref->{'rend'} > $rend);
         $qstart = $ref->{'qstart'} if ($ref->{'qstart'} < $qstart);
         $qend = $ref->{'qend'} if ($ref->{'qend'} > $qend);
         $sum_aln += $ref->{'num_aln'};
         $sum_qbases += $ref->{'qbases'};
     }

     #Create a new hash to push onto data:
     my %hash = (
         qcontig => $hashes[0]->{'qcontig'},
         qstrand => $hashes[0]->{'qstrand'},
         qbases => $sum_qbases,
         qstart => $qstart,
         qend => $qend,
         rcontig => $hashes[0]->{'rcontig'},
         rstart => $rstart,
         rend => $rend,
         num_aln => $sum_aln,
     );

     return (%hash);
}

#Now I can go onto the main function. My idea now is this: have a hash
#of hashes which I add data to. For each query, I built up sections of
#each contig separated by a given gap thresh. Every time one finishes,
#I push it to an array of hashes. After each ref contig, I go through
#this array and "decide" which groups to keep.

sub transloc {
    my %args = @_;
    my $input_sam = $args{'input_sam'};
    my $output_file = $args{'output'};
    my $ref_gap = $args{'ref_gap'};
    my $query_gap = $args{'query_gap'};
    my $max_diff = $args{'max_diff'} / 100;
    my $min_qbase = $args{'min_qbase'};
    my $min_depth = $args{'min_depth'};
    my $min_length = $args{'min_length'};
    my $input_vcf = $args{'input_vcf'};
    my $all = $args{'all'};

    my $output;
    if ($output_file eq 'NO_OUTPUT_PROVIDED') {
        $output = *STDOUT
    }
    else {
        $output = FileHandle->new("> $output_file");
    }

    #Set up my reading of the input file:
    my $input_fh = FileHandle->new("$input_sam");

    #my %data_hash;
    my @data_hashes;
    my @output_lines;
    my $current_contig = '';

  SAM: while (<$input_fh>) {
        my @row = split(/\t/, $_);
        next SAM if ($row[0] =~ /^@/);
        if ($current_contig ne $row[2]) {
            if ($current_contig ne '') {
                #Need to add some metadata to hashes:
                for my $ref (@data_hashes) {
                    $ref->{'rlength'} = $ref->{'rend'} - $ref->{'rstart'} + 1;
                    $ref->{'qlength'} = $ref->{'qend'} - $ref->{'qstart'} + 1;
                    $ref->{'depth'} = $ref->{'qbases'} / $ref->{'rlength'};
                }

                my @new_outputs = parse_outputs(\@data_hashes, $min_depth, $min_qbase, $min_length, $max_diff, $all);
                push(@output_lines, @new_outputs);
            }
            $current_contig = $row[2];
            @data_hashes = ();
        }
        #Because some contig names have _ in them, I can't just split this easily.
        my @arr = split(/_/, $row[0]);
        my @qarr = (join('_', @arr[0..$#arr-4]), @arr[-4..-1]);
        my $end_pos = calc_end_pos($row[3], $row[5]);

        my @data_idx = grep {
            $data_hashes[$_]->{'qcontig'} eq $qarr[0] &&
                $data_hashes[$_]->{'qstrand'} eq $qarr[1] &&
                $data_hashes[$_]->{'qstart'} - $query_gap < $qarr[4] &&
                $data_hashes[$_]->{'qend'} + $query_gap > $qarr[3] &&
                $data_hashes[$_]->{'rstart'} - $ref_gap < $end_pos &&
                $data_hashes[$_]->{'rend'} + $ref_gap > $row[3]
            } 0..$#data_hashes;

        if (scalar(@data_idx) == 1) {
            my $idx = $data_idx[0];
            #This means I have an element for this aln already
            $data_hashes[$idx]->{'rend'} = $end_pos if ($data_hashes[$idx]->{'rend'} < $end_pos);
            $data_hashes[$idx]->{'num_aln'} += 1;
            $data_hashes[$idx]->{'qbases'} += $qarr[2];

            $data_hashes[$idx]->{'qend'} = $qarr[4] if ($qarr[4] > $data_hashes[$idx]->{'qend'});
            $data_hashes[$idx]->{'qstart'} = $qarr[3] if ($qarr[4] < $data_hashes[$idx]->{'qstart'});
        } elsif (scalar(@data_idx) == 0) {
            #This means I need to initialize this element
            my %hash = (
                qcontig => $qarr[0],
                qstrand => $qarr[1],
                qbases => $qarr[2],
                qstart => $qarr[3],
                qend => $qarr[4],
                rcontig => $row[2],
                rstart => $row[3],
                rend => $end_pos,
                num_aln => 1,
            );
            push(@data_hashes, \%hash);
        } else {
            #This means some of my groupings have grown together, so I
            #need to add them together

            my @match = @data_hashes[@data_idx];
            my %hash = merge_hashes(\@match);

            #Remove the combining elements
            @data_idx = sort { $b <=> $a } @data_idx;
            for my $idx (@data_idx) { splice(@data_hashes, $idx, 1) };

            push(@data_hashes, \%hash);
        }
    } #end SAM
    #Catch the final contig:
    if ($current_contig ne '') {
        #Need to add some metadata to hashes:
        for my $ref (@data_hashes) {
            $ref->{'rlength'} = $ref->{'rend'} - $ref->{'rstart'} + 1;
            $ref->{'qlength'} = $ref->{'qend'} - $ref->{'qstart'} + 1;
            $ref->{'depth'} = $ref->{'qbases'} / $ref->{'rlength'};
        }

        my @new_outputs = parse_outputs(\@data_hashes, $min_depth, $min_qbase, $min_length, $max_diff, $all);
        push(@output_lines, @new_outputs);
    }


    #If a vcf file is passed, I need to shift my translocations
    if ($input_vcf) {
        @output_lines = transloc_shifter(\@output_lines, "$input_vcf");
    }

    #Re-calculate lengths since this will change them. I won't do
    #depth, since that metric is supposed to be representative of
    #quality and doesn't need the change I feel
    for my $ref (@output_lines) {
        $ref->{'rlength'} = $ref->{'rend'} - $ref->{'rstart'} + 1;
        $ref->{'qlength'} = $ref->{'qend'} - $ref->{'qstart'} + 1;
    }

    #Re-filter after parsing (I did this during parsing too, but may
    #have changed)
    if (! $all) {
        my $prev_num = scalar(@output_lines);
        @output_lines = grep {
            $_->{'rlength'} >= $min_length &&
                $_->{'depth'} >= $min_depth &&
                abs(($_->{'rlength'} - $_->{'qlength'}) / $_->{'rlength'}) < $max_diff &&
                $_->{'qbases'} > $min_qbase
        } @output_lines;
    }

    print $output
        "#Q_CONTIG\tQ_START\tQ_END\tQ_LENGTH\tQ_STRAND\t",
        "R_CONTIG\tR_START\tR_END\tR_LENGTH\tR_STRAND\t",
        "DEPTH\tQ_BASES\n";
  PRINTER: for my $line (@output_lines) {
        my $strand = '+';
        $strand = '-' if ($line->{'qstrand'} eq 'rev');

        print $output
            $line->{'qcontig'}, "\t",
            $line->{'qstart'}, "\t",
            $line->{'qend'}, "\t",
            $line->{'qlength'}, "\t",
            $strand, "\t",
            $line->{'rcontig'}, "\t",
            $line->{'rstart'}, "\t",
            $line->{'rend'}, "\t",
            $line->{'rlength'}, "\t",
            '+', "\t",
            $line->{'depth'}, "\t",
            $line->{'qbases'}, "\n";
    }
}

#This seeks to shift feature annotations according to the indels
#within a vcf file. To this end, it is passed an input gff, an output
#file, and an input vcf. Right now it is compatible with gff only, but
#it shouldn't be too difficult to modify it to accept other file
#types.

sub gff_shifter {
    my %args = @_;
    my @starts;
    my $gff_fh = FileHandle->new("$args{'gff'}");
    #Read in my VCF file:
    my $vcf_fh = FileHandle->new("$args{'vcf'}");
    my @vcf;
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($vcf_fh, "<:encoding(utf8)", "$args{'vcf'}") or die "Can't open VCF file";
    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'info');

    #Load in my VCF features and add them to a DB:
    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if ($row->{'contig'} =~ /^#/);
        next if (! defined($row->{'info'}));
        my $shift = length((split(/\,/, $row->{'alternate'}))[0]) - length($row->{'reference'});
        next if $shift == 0;
        my %hash = (
            seq_id => $row->{'contig'},
            pos => $row->{'position'},
            strand => $row->{'strand'},
            ref => $row->{'ref'},
            alt => $row->{'alt'},
            shift => $shift
        );

        push(@vcf, \%hash);
    }
    close($vcf_fh);

    #Establish output and check if VCF is empty:
    my $gff_out;
    if ("$args{'output'}" eq 'NO_OUTPUT_PROVIDED') {
        $gff_out = *STDOUT;
    }
    else {
        $gff_out = FileHandle->new("> $args{'output'}");
    }

    if (scalar @vcf == 0) {
        print STDERR "WARNING: No indels in the VCF file, copying\n";
        $gff_fh = FileHandle->new("$args{'gff'}");
        while (<$gff_fh>) {
            print $gff_out $_;
        }
        ;
        exit(1);
    }
    ;

    @vcf = sort {
        $a->{'seq_id'} cmp $b->{'seq_id'} || $a->{'pos'} <=> $b->{'pos'}
    } @vcf;

    my $gff_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($gff_fh, "<:encoding(utf8)", "$args{'gff'}");
    $gff_tsv->column_names('contig', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes');

    my $feat;
    my $gff_id = 0;
    my @heading;
  READ_GFF: while (my $row = $gff_tsv->getline_hr($gff_fh)) {
        if ($row->{'contig'} =~ /^#/) {
            push @heading, $row->{'contig'};
            next READ_GFF;
        }

        #For shifting logic, I will create two structures and sort them independently
        my %hash = (
            seq_id => $row->{'contig'},
            start => $row->{'start'},
            end => $row->{'end'},
            strand => $row->{'strand'},
            type => $row->{'type'},
            score => $row->{'score'},
            source => $row->{'source'},
            phase => $row->{'phase'},
            attributes => $row->{'attributes'},
            gff_id => $gff_id,
        );
        push(@starts, \%hash);
        $gff_id++;
    }
    close($gff_fh);
    #Sorting each DB
    my @ends = @starts;        #These are identical until I sort them:
    @starts = sort {
        $a->{'seq_id'} cmp $b->{'seq_id'} || $a->{'start'} <=> $b->{'start'}
    } @starts;
    @ends = sort {
        $a->{'seq_id'} cmp $b->{'seq_id'} || $a->{'end'} <=> $b->{'end'}
    } @ends;

    my $contig;
    my $idx = 0;
  HEADER: foreach my $line (@heading) {
        if (index($line, 'sequence-region') == -1) {
            print $gff_out $line, "\n";
            next HEADER;
        }
        else {
            $line =~ tr /#//d;
            my @line = split / /, $line;
            my $length = $line[3];
            $contig = $line[1];
            while (my $row = $vcf["$idx"]) {
                last if ($contig ne $row->{'seq_id'});
                $length = $length + $row->{'shift'};
                $idx++;
            }
            print $gff_out "##sequence-region " . "$contig " . "1 " . "$length\n";
        }
    }

    my $gff_row;
    my $gff_index = 0;
    my @shifted_starts;
    my @shifted_ends;

    my $shiftsum = 0;
    my $shift = 0;
    my $varpos = 0;
    my $var_contig = '';
    my @seen_contigs= ('');
    my $vcf_row;
    my $feat_start;
    $idx = 0;
  STARTS: while ($gff_row = $starts["$gff_index"]) {
        $contig = $gff_row->{'seq_id'};
        $feat_start = $gff_row->{'start'};
        if (! grep /$contig/, @seen_contigs) {
            $shiftsum = 0;
            push @seen_contigs, $contig;
            $shift = 0;
        }
      START_SHIFT: while ($vcf_row = $vcf["$idx"]) {
            my $tmp_contig = $vcf_row->{'seq_id'};
            if ($contig ne $tmp_contig) {
                if (! grep /$tmp_contig/, @seen_contigs) {
                    #Means this variant is on the next contig
                    last START_SHIFT;
                }
                else {
                    #Means variant is on previous chromosome
                    $idx++;
                    next START_SHIFT;
                }
            }
            if ($feat_start < $vcf_row->{'pos'}) {
                last START_SHIFT;
            }
            #If we get here, it means the variant is on the right
            #chromosome and upstream of our feature
            $shift = $vcf_row->{'shift'};
            $varpos = $vcf_row->{'pos'}; #kept for later logic
            $shiftsum = $shiftsum + $shift;
            $var_contig = $vcf_row->{'seq_id'};
            $idx++;
        }
        ;
        #Now I should have collected the proper shiftsum:
        if (($varpos - $shift > $feat_start) & ($var_contig eq $contig) & ($idx != 0)) {
            $gff_row->{'start'} = $varpos + $shiftsum - $shift;
        }
        else {
            $gff_row->{'start'} = $shiftsum + $feat_start;
        }
        push(@shifted_starts, $gff_row);
        $gff_index++;
    }

    #Now I need to repeat this with my end positions, which requires
    #resetting my variables:
    $shiftsum = 0;
    $shift = 0;
    $varpos = 0;
    $var_contig = '';
    @seen_contigs = '';
    $gff_index = 0;
    $idx = 0;
    my $feat_end;
  ENDS: while ($gff_row = $ends["$gff_index"]) {
        $contig = $gff_row->{'seq_id'};
        $feat_end = $gff_row->{'end'};
        if (! grep /$contig/, @seen_contigs) {
            $shiftsum = 0;
            push @seen_contigs, $contig;
            $shift = 0;
        }
      END_SHIFT: while ($vcf_row = $vcf["$idx"]) {
            my $tmp_contig = $vcf_row->{'seq_id'};
            if ($contig ne $tmp_contig) {
                if (! grep /$tmp_contig/, @seen_contigs) {
                    #Means this variant is on the next contig
                    last END_SHIFT;
                }
                else {
                    #Means variant is on previous chromosome
                    $idx++;
                    next END_SHIFT;
                }
            }
            if ($feat_end < $vcf_row->{'pos'}) {
                last END_SHIFT;
            }
            #If we get here, it means the variant is on the right
            #chromosome and upstream of our feature
            $shift = $vcf_row->{'shift'};
            $varpos = $vcf_row->{'pos'}; #kept for later logic
            $shiftsum = $shiftsum + $shift;
            $var_contig = $vcf_row->{'seq_id'};
            $idx++;
        }
        #Now I should have collected the proper shiftsum:
        if (($varpos - $shift > $feat_end) & ($var_contig eq $contig) & ($idx != 0)) {
            $gff_row->{'end'} = $varpos + $shiftsum - $shift;
        }
        else {
            $gff_row->{'end'} = $feat_end + $shiftsum;
        }
        push(@shifted_ends, $gff_row);
        $gff_index++;
    }

    #Now I need to sort each database according to its tag value
    #"gff_id":
    @shifted_starts = sort {
        $a->{'gff_id'} <=> $b->{'gff_id'}
    } @shifted_starts;
    @shifted_ends = sort {
        $a->{'gff_id'} <=> $b->{'gff_id'}
    } @shifted_ends;


    #Now these should each be back to their original sorting and
    #identical in all but positions:
    $idx = 0;
    while ($gff_row = $shifted_starts["$idx"]) {
        print $gff_out $gff_row->{'seq_id'} . "\t" .
            $gff_row->{'source'} . "\t" .
            $gff_row->{'type'} . "\t" .
            $gff_row->{'start'} . "\t" .
            $shifted_ends[$idx]->{'end'} . "\t" .
            $gff_row->{'score'} . "\t" .
            $gff_row->{'strand'} . "\t" .
            $gff_row->{'phase'} . "\t" .
            $gff_row->{'attributes'} . "\n";
        $idx++;
    }
}

#This next subroutine seeks to take a given VCF, Fasta, and GFF file
#and from them assign amino acid changes I still corresponding to
#them. I also hope to assign a functnaiol score to the mutation using
#something like the BLOSUM index

sub label {
    my %args = @_;
    my $input_gff = $args{'gff'};
    my $input_vcf = $args{'vcf'};
    my $input_fasta = $args{'fasta'};
    my @feature_types = split(/\,/, $args{'features'}); #For example exon, gene, etc.
    my $child_name = $args{'child_name'}; #For example ID
    my $parent_name = $args{'parent_name'}; #Should just be "Parent", but I'm making sure
    my $add_flanks = $args{'add_flanks'};
    my @flank_lengths = split(/\,/, $args{'flank_lengths'});
    my $flank_feature = $args{'flank_feature'};
    my $flank_parent = $args{'flank_parent'};
    my $threads = $args{'threads'};
    my $output_nuc = $args{'output_nuc'};
    my $translate = $args{'translate'};
    my $coding_feature = $args{'coding_feature'};
    my $codon_table = $args{'codon_table'};
    my $score_matrix = $args{'score_matrix'};
    my $output_aa = $args{'output_aa'};
    my $tmpdir = $args{'tmpdir'};
    my $all_vars = $args{'all_vars'};

    chdir($args{'outdir'});

    #Modify my features if flanks are added:
    if ($add_flanks) {
        push(@feature_types, ('up_flank', 'down_flank'));
    }

    #Develop method to determine which features to add
    my %feature_types = map {$_ => 1 } @feature_types;

    #Add functionality for GTF file:
    my $att_sep = '=';
    if ($input_gff =~ /.gtf$/) {
        $att_sep = ' ';
    }
    #Read in my GFF file:
    my @gff;
    my $gff_fh = FileHandle->new("$input_gff");
    my $gff_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($gff_fh, "<:encoding(utf8)", "$input_gff");

    $gff_tsv->column_names('seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes');

    my $feat;
    my $protpos = 1;
    my $outnuc = FileHandle->new("> ${output_nuc}");
    print $outnuc "#CHROM\tPOS\tREF\tALT\tPARENT\tCHILD\tTYPE\tSTRAND\tPOS_PARENT\tPOS_CHILD\tPOS_CODING\tDEPTH\n";
    while (my $row = $gff_tsv->getline_hr($gff_fh)) {
        next if($row->{'seq_id'} =~ /^#/);
        my $feat_type = $row->{'type'};
        #next if(! exists($feature_types{$feat_type}));
        if (exists($feature_types{$feat_type})) {
            my @tmp = split(/;/, $row->{'attributes'});
            my %attributes;
            for my $item (@tmp) {
                my ($key, $value) = split("$att_sep", $item);
                $attributes{$key} = $value;
            }

            my $parent = $attributes{$parent_name};
            my $child = $attributes{$child_name};
            my %feat = (seq_id => $row->{'seq_id'},
                        start => $row->{'start'},
                        end => $row->{'end'},
                        strand => $row->{'strand'},
                        type => $feat_type,
                        parent => $parent,
                        child => $child
                    );
            if (%feat) {
                push(@gff, \%feat);
            }

        }

        #Now I need to figure out whether to add some flanks. This is
        #separate because in some instances you may want to add flanks
        #to the edges of the protein_coding_gene, and annotate
        #variants by CDS

        if (("$feat_type" eq "$flank_feature") & ($add_flanks)) {
            my @tmp = split(/;/, $row->{'attributes'});
            my %attributes;
            for my $item (@tmp) {
                my ($key, $value) = split("$att_sep", $item);
                $attributes{$key} = $value;
            }
            ;
            my $parent = $attributes{$flank_parent};
            my $child = $attributes{$child_name};
            if ($row->{'strand'} eq '-') {
                if ($flank_lengths[0] != 0) {
                    my %up = (seq_id => $row->{'seq_id'},
                              start => $row->{'end'} + 1,
                              end => $row->{'end'} + $flank_lengths[0],
                              strand => '-',
                              type => 'up_flank',
                              parent => $parent,
                              child => "up_flank_" . $child
                          );
                    push(@gff, \%up);
                }
                if ($flank_lengths[1] != 0) {
                    #Need to make sure the downstream flank doesn't go below zero
                    my $tmp_start = $row->{'start'} - $flank_lengths[1];
                    if ($tmp_start < 1) {
                        $tmp_start = 1;
                    }

                    my $tmp_end = $row->{'start'} - 1;
                    unless ($tmp_start >  $tmp_end) {
                        my %down = (seq_id => $row->{'seq_id'},
                                    start => $tmp_start,
                                    end => $tmp_end,
                                    strand => '-',
                                    type => 'down_flank',
                                    parent => $parent,
                                    child => "down_flank_" . $child
                                );
                        push(@gff, \%down);
                    }
                }
            }
            else {
                #I think any good gff has only + or - strands, but
                #this is general so it will treat anything without "-"
                #as a plus strand item
                if ($flank_lengths[0] != 0) {
                    #Need to make sure upstream flank doesn't go below zero
                    my $tmp_start = $row->{'start'} - $flank_lengths[0];
                    if ($tmp_start < 1) {
                        $tmp_start = 1;
                    }

                    my $tmp_end = $row->{'start'} - 1;
                    unless ($tmp_start > $tmp_end) {
                        my %up = (seq_id => $row->{'seq_id'},
                                  start => $tmp_start,
                                  end => $tmp_end,
                                  strand => $row->{'strand'},
                                  type => 'up_flank',
                                  parent => $parent,
                                  child => "up_flank_" . $child
                              );
                        push(@gff, \%up);
                    }
                }
                if ($flank_lengths[1] != 0) {
                    my %down = (seq_id => $row->{'seq_id'},
                                start => $row->{'end'} + 1,
                                end => $row->{'end'} + $flank_lengths[1],
                                strand => $row->{'strand'},
                                type => 'down_flank',
                                parent => $parent,
                                child => "down_flank_" . $child
                            );
                    push(@gff, \%down);
                }
            }
        }
    }

    if (scalar(@gff) == 0) {
        die "GFF file (${input_gff}) has no features\n";
    }
    if (! defined($gff[0]->{'parent'})) {
        die "ERROR: Could not find parent ${parent_name} in GFF/GTF attributes\n";
    }
    if (! defined($gff[0]->{'child'})) {
        die "ERROR: Could not find child ${child_name} in GFF/GTF attributes\n";
    }

    #I'm going to break this into two GFF files
    my @plus_features = grep { $_->{'strand'} eq '+' } @gff;
    my @minus_features = grep { $_->{'strand'} eq '-' } @gff;

    #Now that I've loaded my GFF I want to load my VCF.
    my $vcf_fh = FileHandle->new("$input_vcf");
    my @vcf;
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $vcf_tsv->column_names('seq_id', 'pos', 'index', 'ref', 'alt', 'qual', 'filter', 'info');
    open($vcf_fh, "<:encoding(utf8)", "$input_vcf") or die "Can't open vcf_fh";

    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if $row->{'seq_id'} =~ /^#/;
        next if (! defined($row->{'info'}));
        my @tmp = split /;/, $row->{'info'};
        my %info;
        for my $item (@tmp) {
            my ($key, $value) = split('=', $item);
            $info{$key} = $value;
        }

        my %var = (
            seq_id => $row->{'seq_id'},
            pos => $row->{'pos'},
            strand => 1,
            ref => $row->{'ref'},
            alt => $row->{'alt'},
            depth => $info{'DP'}
        );
        push(@vcf, \%var);
    }

    #Now I need to determine all of the unique contigs I have:
    my @unique_contigs;
    for my $feat (@vcf) {
        push(@unique_contigs, $feat->{'seq_id'});
    }
    @unique_contigs = sort { $a cmp $b } @unique_contigs;
    my $prev = 'FOOBAR';
    @unique_contigs = grep($_ ne $prev && ($prev = $_), @unique_contigs);

    my %positions; #use this to determine when I need to reset my variant position
    #Then, I'm going to separate the GFF into a plus and minus strand
    #using grep. I will then loop through these individually, keeping
    #track of when I change "parent" gene, and produce a log of the
    #amino acid changes
    #I'm going to make this run in parallel to make it somewhat acceptable
    my $pm = Parallel::ForkManager->new($threads);
    my $death = 0;              #Determines if I "die" or not
  CONTIGS: for my $contig (@unique_contigs) {
        my $pid = $pm->start and next CONTIGS;
        my @plus_sub_features = grep { $_->{'seq_id'} eq $contig } @plus_features;
        my @variants = grep { $_->{'seq_id'} eq $contig } @vcf;
        if (scalar(@variants) == 0) {
            $pm->finish;
        }
        ;
        #I don't think I need to sort either group earlier, I should
        #just sort them now for simplicity
        @variants = sort {
            $a->{'pos'} <=> $b->{'pos'}
        } @variants;
        @plus_sub_features = sort {
            $a->{'parent'} cmp $b->{'parent'} ||
                $a->{'start'} <=> $b->{'start'}
            } @plus_sub_features;
        #Set up features no minus strand:
        my @minus_sub_features = grep { $_->{'seq_id'} eq $contig } @minus_features;
        @minus_sub_features = sort {
            $b->{'parent'} cmp $a->{'parent'} ||
                $b->{'start'} <=> $a->{'start'}
            } @minus_sub_features;


        #Now I need to iterate over each feature, keeping track of the parent
        my @labels;    #Instead of printing everytime, which can cause
        #overlap in prints, save info to an array of hashes
        #and print at the end
        my @labeled_idx; #keep track of which variants I've printed

        if (scalar(@plus_sub_features) !=  0) {
            my $parent = $plus_sub_features[0]->{'parent'}; #Start this up
            my $parent_pos = 0; #track number of bases covered in parent
            my $coding_pos = 0; #track number of bases covered in coding sequence
            my $feat_end = 0; #make sure it doesn't fail in first iteration
            my $feat_start;
            my $var_idx = 0;
            my $var = undef;
            my $skipper = '_balrog_'; #if I have overlap, this determines which parent to skip
          FEATURES: for my $feat (@plus_sub_features) {
                next FEATURES if ($feat->{'parent'} eq $skipper);
                if ($parent ne $feat->{'parent'}) {
                    $parent_pos = 0;
                    $coding_pos = 0;
                    $positions{$feat_end} = $var_idx;
                    $parent = $feat->{'parent'};
                    #Determine when I need to reset the variant search:
                    my @options = grep { $_ < ($feat->{'start'}) } (keys %positions);
                    if (scalar(@options) > 0) {
                        my $recent = max @options;
                        $var_idx = $positions{$recent};
                    }
                    else {
                        $var_idx = 0;
                    }
                    last FEATURES if (! defined($var));
                }
                elsif ($feat_end > $feat->{'start'}) {
                    #This means we are in the same parent and have overlap, not good
                    $skipper = $feat->{'parent'};
                    print "WARNING: Features overlap within same parent, skipping $feat->{'parent'}\n";
                    next FEATURES;
                }

                $feat_start = $feat->{'start'};
                $feat_end = $feat->{'end'};
              VARIANTS: while ($var = $variants[$var_idx]) {
                    my $var_pos = $var->{'pos'};
                    if ($var_pos < $feat_start ) {
                        $var_idx++;
                        next VARIANTS;

                    }
                    elsif ($var_pos > $feat_end) {
                        $parent_pos = $parent_pos + $feat_end - $feat_start + 1;
                        if ($feat->{'type'} eq "$coding_feature") {
                            $coding_pos = $coding_pos + $feat_end - $feat_start + 1;
                        }
                        next FEATURES;
                    }

                    my $relpos = $var_pos - $feat_start + 1;
                    my $cd;
                    if ($feat->{'type'} eq "$coding_feature") {
                        $cd = $relpos + $coding_pos;
                    }
                    else {
                        $cd = 'NA';
                    }

                    my %hash = (
                        contig => $contig,
                        var_pos => $var_pos,
                        refn => $var->{'ref'},
                        altn => $var->{'alt'},
                        parent => $feat->{'parent'},
                        child => $feat->{'child'},
                        type => $feat->{'type'},
                        strand => '+',
                        parent_pos => ($relpos + $parent_pos),
                        child_pos => ($relpos),
                        coding_pos => $cd,
                        depth => $var->{'depth'}
                    );
                    push(@labels, \%hash);
                    push(@labeled_idx, $var_idx);
                    $var_idx++;
                } #End Variants
            } #End Features
        } #end scalar(@plus_sub_features) check

        if (scalar(@minus_sub_features) != 0) {

            #Now I'm going to repeat this with a slightly modified pipeline
            #for the minus strand variants
            #I need to reload features because it differs by strand,
            #variants don't

            @variants = sort {
                $b->{'pos'} <=> $a->{'pos'}
            } @variants;
            my %positions = ();
            my $parent = $minus_sub_features[0]->{'parent'}; #Start this up
            my $parent_pos = 0; #track number of bases covered in parent
            my $coding_pos = 0;
            my $feat_start = 100000000000; #make sure it doesn't fail on first iteration
            my $feat_end;
            my $var_idx = 0;
            my $var = undef;
            my $skipper = '_balrog_';
          FEATURES: for my $feat (@minus_sub_features) {
                #For whatever reason, the above sorting can create an empty
                #entry in my hash, so I need to add this caveat
                next FEATURES if($feat->{'parent'} eq $skipper);
                if ($parent ne $feat->{'parent'}) {
                    $parent_pos = 0;
                    $coding_pos = 0;
                    $positions{$feat_start} = $var_idx; #feat_start b/c minus strand
                    $parent = $feat->{'parent'};
                    #Determine when I need to reset the variant search:
                    #Unlike with plus strand, I am going in reverse
                    #order. This means I am checking relative to the ends of
                    #features
                    my @options = grep { $_ > ($feat->{'end'}) } (keys %positions);
                    if (scalar(@options) > 0) {
                        my $recent = min @options;
                        $var_idx = $positions{$recent};
                    }
                    else {
                        $var_idx = 0;
                    }
                    last FEATURES if (! defined($var));
                }
                elsif ($feat_start < $feat->{'end'}) {
                    $skipper = $feat->{'parent'};
                    print "WARNING: Features overlap within the same parent, skipping $feat->{'parent'}\n";
                    next FEATURES;
                }
                ;
                $feat_start = $feat->{'start'};
                $feat_end = $feat->{'end'};
              VARIANTS: while ($var = $variants[$var_idx]) {
                    my $var_pos = $var->{'pos'}; #Some extra logic for - strand
                    my $refn = $var->{'ref'};
                    $var_pos = $var_pos + length($refn) - 1; #Change position for deletions
                    if ($var_pos > $feat_end) {
                        $var_idx++;
                        next VARIANTS;
                    }
                    elsif ($var_pos < $feat_start) {
                        $parent_pos = $parent_pos + $feat_end - $feat_start + 1;
                        if ($feat->{'type'} eq "$coding_feature") {
                            $coding_pos = $coding_pos + $feat_end - $feat_start + 1;
                        }
                        next FEATURES;
                    }

                    $refn =~ tr/ACGT/TGCA/;
                    my $altn = $var->{'alt'};
                    $altn =~ tr/ACGT/TGCA/;
                    my $relpos = $feat_end - $var_pos + 1;
                    my $cd;
                    if ($feat->{'type'} eq "$coding_feature") {
                        $cd = $relpos + $coding_pos;
                    }
                    else {
                        $cd = 'NA';
                    }
                    my %hash = (
                        contig => $contig,
                        var_pos => $var_pos,
                        refn => $refn,
                        altn => $altn,
                        parent => $feat->{'parent'},
                        child => $feat->{'child'},
                        type => $feat->{'type'},
                        strand => '-',
                        parent_pos => ($relpos + $parent_pos),
                        child_pos => ($relpos),
                        coding_pos => $cd,
                        depth => $var->{'depth'}
                    );
                    push(@labels, \%hash);
                    push(@labeled_idx, $var_idx);
                    $var_idx++;
                }
            }
        }
        #if the user wants all variants printed, I need to add those
        #to my labels with a bunch of NA:
        if ($all_vars) {
            my %labeled_idx = map { $_ => 1 } @labeled_idx;
            @variants = @variants[grep { ! exists($labeled_idx{$_}) } 0..$#variants];
            for my $var (@variants) {
                my %hash = (
                    contig => $var->{'seq_id'},
                    var_pos => $var->{'pos'},
                    refn => $var->{'ref'},
                    altn => $var->{'alt'},
                    parent => 'NA',
                    child => 'NA',
                    type => 'NA',
                    strand => 'NA',
                    parent_pos => 'NA',
                    child_pos => 'NA',
                    coding_pos => 'NA',
                    depth => 'NA'
                );
                push (@labels, \%hash);
            }
        }



        #Now I've built up my @labels and need to print it
        if (scalar(@labels) != 0) {
            @labels = sort { $a->{'var_pos'} <=> $b->{'var_pos'} } @labels;
            my $tmp_fh = FileHandle->new("> ${tmpdir}/$contig");
            for my $feat (@labels) {
                print $tmp_fh
                    $feat->{'contig'}, "\t",
                    $feat->{'var_pos'}, "\t",
                    $feat->{'refn'}, "\t",
                    $feat->{'altn'}, "\t",
                    $feat->{'parent'}, "\t",
                    $feat->{'child'}, "\t",
                    $feat->{'type'}, "\t",
                    $feat->{'strand'}, "\t",
                    $feat->{'parent_pos'}, "\t",
                    $feat->{'child_pos'}, "\t",
                    $feat->{'coding_pos'}, "\t",
                    $feat->{'depth'}, "\n";
            }
            close($tmp_fh);
        }
        $pm->finish;
    }
    if ("$death" == 1) {
        return;
    }
    $pm->wait_all_children;
    #Combine all of the temporary files:
    for my $contig (@unique_contigs) {
        if (-e "${tmpdir}/${contig}") {
            my $fh = FileHandle->new("< ${tmpdir}/${contig}");
            while (<$fh>) {
                print $outnuc "$_";
            }
            close($fh);
        }
    }
    close($outnuc);
    rmtree("$tmpdir");

    if (! $translate) {
        exit 0;
    }
    ;
    my $fasta = Bio::DB::Fasta->new("$input_fasta");
    my $parser = Bio::Matrix::IO->new(-format => 'scoring', -file => "$FindBin::Bin/../share/matrices/${score_matrix}");
    my $matrix = $parser->next_matrix;
    #Ok now I've generated the labeled variants file, which I've
    #decided is most effective to read from and do later
    #analysis. Yes, I could just save the info from the last half, but
    #that isn't compatible with the parallelization

    #I think the first thing I want to do is load in the sequence of
    #each feature into a hash:
    my %feature_seqs;
    my $sequence = '';
    my $parent;
    @plus_features = sort {
        $a->{'parent'} cmp $b->{'parent'} ||
            $a->{'start'} <=> $b->{'start'}
        } @plus_features;
    $parent = $plus_features[0]->{'parent'};
    for my $feat (@plus_features) {
        next if ($feat->{'type'} ne "$coding_feature");
        if ($feat->{'parent'} ne $parent) {
            $feature_seqs{$parent} = $sequence;
            $parent = $feat->{'parent'};
            $sequence = $fasta->seq($feat->{'seq_id'}, $feat->{'start'}, $feat->{'end'});
        }
        else {
            $sequence = ($sequence) . ($fasta->seq($feat->{'seq_id'}, $feat->{'start'}, $feat->{'end'}));
        }
    }
    #Catch the final feature (the above won't catch it):
    $feature_seqs{$parent} = $sequence;
    #Repeat for minus features
    @minus_features = sort {
        $a->{'parent'} cmp $b->{'parent'} ||
            $a->{'start'} <=> $b->{'start'}
        } @minus_features;
    $parent = $minus_features[0]->{'parent'};
    for my $feat (@minus_features) {
        next if ($feat->{'type'} ne "$coding_feature");
        if ($feat->{'parent'} ne $parent) {
            $sequence =~ tr/ACGT/TGCA/;
            $sequence = reverse($sequence);
            $feature_seqs{$parent} = $sequence;
            $parent = $feat->{'parent'};
            $sequence = $fasta->seq($feat->{'seq_id'}, $feat->{'start'}, $feat->{'end'});
        }
        else {
            $sequence = $sequence . $fasta->seq($feat->{'seq_id'}, $feat->{'start'}, $feat->{'end'});
        }
    }
    #Now I need to catch the output final feature
    $feature_seqs{$parent} = $sequence;

    #Now I only need to read in this file labeled by CDS
    my @labeled_vars;
    my $labeled_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    my $labeled_fh = FileHandle->new("${output_nuc}");
    $labeled_tsv->column_names('seq_id', 'pos', 'ref', 'alt', 'parent', 'child', 'type', 'strand', 'rel_parent', 'rel_child', 'rel_coding', 'depth');
    while (my $row = $labeled_tsv->getline_hr($labeled_fh)) {
        next if($row->{'seq_id'} =~ /^#/);
        next if($row->{'type'} ne "$coding_feature");
        my $codon = int(($row->{'rel_coding'} - 1) / 3) + 1;
        my $codon_pos = (($row->{'rel_coding'} - 1) % 3) + 1;
        my %hash = (seq_id => $row->{'seq_id'},
                    ref => $row->{'ref'},
                    alt => $row->{'alt'},
                    parent => $row->{'parent'},
                    rel_parent => $row->{'rel_parent'},
                    rel_coding => $row->{'rel_coding'},
                    codon => $codon,
                    codon_pos => $codon_pos
                );
        push(@labeled_vars, \%hash);
    }
    @labeled_vars = sort {
        $a->{'parent'} cmp $b->{'parent'} ||
            $a->{'rel_coding'} <=> $b->{'rel_coding'}
        } @labeled_vars;

    #Important (but weird). The way the below code works is by
    #accumulating data and only printing when its onto the next
    #variant. This means the last variant isn't printed. As such, I
    #will add a dumby variable to the end of my @labeled_vars that
    #will not get printed
    my %hash = (seq_id => 'foobar', ref => 'foobar', alt => 'foobar', parent => $labeled_vars[0]->{'parent'}, rel_parent => 0, rel_coding => 0, codon => 1, codon_pos => 1);
    push (@labeled_vars, \%hash);

    my $var_idx = 0;
    my $ref_codon = '';
    my $alt_codon;
    my $codon_num = 'FOOBAR';
    my @nvar;
    $parent = 'FOOBAR';
    my $complex_vars = 0;
    my $contig;

    #Define my codon table, I need to implement an option separate for
    #all the different ID values.
    $codon_table = Bio::Tools::CodonTable->new( -id => "$codon_table");

    my $outaa = FileHandle->new("> ${output_aa}");
    print $outaa "#CHROM\tPARENT\tCODON\tREF\tALT\tVARS\tTYPE\tFRAME\tSCORE\tSCALED_SCORE\n";

  CODONS: while (my $var = $labeled_vars[$var_idx]) {
        #If I'm still on the same codon:
        if (($codon_num eq $var->{'codon'}) & ($parent eq $var->{'parent'})) {
            #Ensure this works with multiallelic sites:
            my $alt = (split(/,/, $var->{'alt'}))[0];
            substr($alt_codon, $var->{'codon_pos'} - 1, 1) = $alt;
            @nvar = (@nvar, $var->{'ref'} . $var->{'rel_coding'} . $alt);
            $var_idx++;
            next CODONS;
        }
        #This means I'm on a new codon and should print info from the last one
        if (($var_idx != 0) & (index(lc($ref_codon), 'n') == -1)) {
            if ((length($ref_codon) % 3 == 0) & (length($alt_codon) == 3)) {
                my $ref_aa = '';
                my $idx = 0;
                my $score = 0;
                my $scaled_score = 0;
                my $alt_aa = $codon_table->translate($alt_codon);
                while ($idx < length($ref_codon)) {
                    my $tmp = $codon_table->translate(substr($ref_codon, $idx, 3));
                    $ref_aa = $ref_aa . $tmp;
                    $score = $score + $matrix->entry($tmp, $alt_aa);
                    $scaled_score = $scaled_score + $matrix->entry($tmp, $alt_aa) - $matrix->entry($tmp, $tmp);
                    $idx = $idx + 3;
                }
                my $type = 'point';
                if (length($ref_aa) != 1) {
                    $type = 'del';
                }
                print $outaa
                    $contig, "\t",
                    $parent, "\t",
                    $codon_num, "\t",
                    $ref_aa, "\t",
                    $alt_aa, "\t",
                    join(',', @nvar), "\t",
                    $type, "\t",
                    "in", "\t",
                    $score, "\t",
                    $scaled_score, "\n";

            }
            elsif ((length($alt_codon) % 3 == 0) & (length($ref_codon) == 3)) {
                my $idx = 0;
                my $alt_aa = '';
                my $score = 0;
                my $scaled_score = 0;
                my $ref_aa = $codon_table->translate($ref_codon);
                while ($idx < length($alt_codon)) {
                    my $tmp = $codon_table->translate(substr($alt_codon, $idx, 3));
                    $alt_aa = $alt_aa . $tmp;
                    $score = $score + $matrix->entry($ref_aa, $tmp);
                    $scaled_score = $scaled_score + $matrix->entry($ref_aa, $tmp) - $matrix->entry($ref_aa, $ref_aa);
                    $idx = $idx + 3;
                }
                my $type = 'point';
                if (length($alt_aa) != 1) {
                    $type = 'ins'
                }
                print $outaa
                    $contig, "\t",
                    $parent, "\t",
                    $codon_num, "\t",
                    $ref_aa, "\t",
                    $alt_aa, "\t",
                    join(',', @nvar), "\t",
                    $type, "\t",
                    "in", "\t",
                    $score, "\t",
                    $scaled_score, "\n";
            }
            else {
                #Implies a frameshift
                my $ref_aa = $codon_table->translate(substr($ref_codon, 0, 3));
                my $alt_aa = $codon_table->translate(substr($alt_codon, 0, 3));
                if (length($ref_codon) > length($alt_codon)) {
                    print $outaa
                        $contig, "\t",
                        $parent, "\t",
                        $codon_num, "\t",
                        $ref_aa . '@', "\t",
                        $alt_aa, "\t",
                        join(',', @nvar), "\t",
                        'del', "\t",
                        'out', "\t",
                        "NA", "\t",
                        "NA", "\n";
                }
                elsif (length($ref_codon) < length($alt_codon)) {
                    print $outaa
                        $contig, "\t",
                        $parent, "\t",
                        $codon_num, "\t",
                        $ref_aa, "\t",
                        $alt_aa . '@', "\t",
                        join(',', @nvar), "\t",
                        'ins', "\t",
                        'out', "\t",
                        "NA", "\t",
                        "NA", "\n";
                }
                else {
                    $complex_vars++;
                }
            }
        }
        $codon_num = $var->{'codon'};
        $contig = $var->{'seq_id'};
        $parent = $var->{'parent'};
        $ref_codon = substr($feature_seqs{$parent}, (($codon_num * 3) - 3), 3);
        $alt_codon = $ref_codon;
        my $alt = (split(/,/, $var->{'alt'}))[0];
        substr($ref_codon, $var->{'codon_pos'} - 1, 1) = $var->{'ref'}; #necessary for complex
        substr($alt_codon, $var->{'codon_pos'} - 1, 1) = $alt;
        @nvar = ($var->{'ref'} . $var->{'rel_coding'} . $alt);

        $var_idx++;
    }
    if ($complex_vars != 0) {
        print "WARNING: ${complex_vars} complex variants encountered, omitted\n";
    }
}

#Ok and now for what might be the final tool of dantools before
#publication, dantools summary. This will take either the nucleotide
#or amino acid input from dantools translate's output and create a
#summary about it. This includes things like total mutations, change
#in length, number of frameshifts (if any), etc.

sub summarize_aa {
    my %args = @_;
    my $input_vars = $args{'input_vars'};
    my $input_vars_fh = FileHandle->new("$input_vars");
    my $output;
    my $outscore = $args{'outscore'};
    if ($args{'output'} ne 'NO_OUTPUT_PROVIDED') {
        $output = FileHandle->new("> $args{'output'}");
    }
    else {
        $output = *STDOUT;
    }

    my $total_score = 0;
    my $total_sscore = 0;
    my $num_vars = 0;
    my $num_indels = 0;
    my $num_outframe = 0;
    my $parent = '';
    my %all_data;
    my @parent_keys;

    print $output "#PARENT\tVARS\tSCORE\tSCALED_SCORE\tSILENT\tMISSENSE\tNONSTOP\tNONSENSE\tINDELS\tINFRAME\tOUTFRAME\n";

    my $vars_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $vars_tsv->column_names('seq_id', 'parent', 'codon', 'ref', 'alt', 'vars', 'type', 'frame', 'score', 'scaled_score');
  VARS: while (my $entry = $vars_tsv->getline_hr($input_vars_fh)) {
          next if $entry->{'seq_id'} =~ /^#/;

          if ("$parent" ne $entry->{'parent'}) {
              $parent = $entry->{'parent'};
              if (! exists($all_data{"$parent"})) {
                  $all_data{"$parent"} = {
                      num_vars => 0,
                      total_score => 0,
                      total_sscore => 0,
                      num_indels => 0,
                      num_outframe => 0,
                      num_silent => 0,
                      num_missense => 0,
                      num_nonstop => 0,
                      num_nonsense => 0
                  };
                  push (@parent_keys, $parent);
              }
        }
          $all_data{"$parent"}{'num_vars'}++;
          if ($entry->{'type'} ne 'point') {
            $all_data{"$parent"}{'num_indels'}++;
            if ($entry->{'frame'} eq 'out') {
                $all_data{"$parent"}{'num_outframe'}++;
                $all_data{"$parent"}{'total_score'} += $outscore;
                $all_data{"$parent"}{'total_sscore'} += $outscore;
                next VARS; #need to move because score of @ is NA
            }
        } else {
            #Get some extra stats about other types
            if ($entry->{'ref'} eq $entry->{'alt'}) {
                $all_data{"$parent"}{'num_silent'} += 1;
            } elsif ($entry->{'alt'} eq '*') {
                $all_data{"$parent"}{'num_nonsense'} += 1;
            } elsif ($entry->{'ref'} eq '*') {
                $all_data{"$parent"}{'num_nonstop'} += 1;
            } else {
                $all_data{"$parent"}{'num_missense'} += 1;
            }
        }
          $all_data{"$parent"}{'total_score'} += $entry->{'score'};
          $all_data{"$parent"}{'total_sscore'} += $entry->{'scaled_score'};
      }

    for my $i (@parent_keys) {
        print $output
            $i, "\t",
            $all_data{"$i"}{'num_vars'}, "\t",
            $all_data{"$i"}{'total_score'}, "\t",
            $all_data{"$i"}{'total_sscore'}, "\t",
            $all_data{"$i"}{'num_silent'}, "\t",
            $all_data{"$i"}{'num_missense'}, "\t",
            $all_data{"$i"}{'num_nonstop'}, "\t",
            $all_data{"$i"}{'num_nonsense'}, "\t",
            $all_data{"$i"}{'num_indels'}, "\t",
            $all_data{"$i"}{'num_indels'} - $all_data{"$i"}{'num_outframe'}, "\t",
            $all_data{"$i"}{'num_outframe'}, "\n";
    }
}

#The amino acid summarize is great, but especially for flank sequences
#I think it would be good to have a summarize-nuc function. This will
#be very similar, almost identical to the summarize_aa function

sub summarize_nuc {
    my %args = @_;
    my $input_vars = $args{'input_vars'};
    my $input_vars_fh = FileHandle->new("$input_vars");
    my $output;
    if ($args{'output'} ne 'NO_OUTPUT_PROVIDED') {
        $output = FileHandle->new("> $args{'output'}");
    } else {
        $output = *STDOUT;
    }

    my $feat_vars = 0;
    my $up_flank_vars = 0;
    my $down_flank_vars = 0;
    my $feat_indels = 0;
    my $up_flank_indels = 0;
    my $down_flank_indels = 0;
    my $parent = ''; #can't be matched
    my @parent_keys; #keep track of the order of my keys

    my %all_data;
    my %single_data = (
        feat_vars => 0,
        up_flank_vars => 0,
        down_flank_vars => 0,
        feat_indels => 0,
        up_flank_indels => 0,
        down_flank_indels =>  0
        );

    print $output "#PARENT\tFEAT_VARS\tUP_FLANK_VARS\tDOWN_FLANK_VARS\tFEAT_INDELS\tUP_FLANK_INDELS\tDOWN_FLANK_INDELS\n";

    my $vars_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $vars_tsv->column_names('seq_id', 'pos', 'ref', 'alt', 'parent', 'child', 'type', 'strand', 'pos_parent', 'pos_child', 'pos_coding', 'depth');

  VARS: while (my $entry = $vars_tsv->getline_hr($input_vars_fh)) {
        next if $entry->{'seq_id'} =~ /^#/;
        #Ensure this works with multiallelic sites:
        my $alt = (split(/,/, $entry->{'alt'}))[0];

        if ("$parent" ne $entry->{'parent'}) {
            $parent = $entry->{'parent'};
            if (! exists($all_data{"$parent"})) {
                $all_data{"$parent"} = {
                    feat_vars => 0,
                    up_flank_vars => 0,
                    down_flank_vars => 0,
                    feat_indels => 0,
                    up_flank_indels => 0,
                    down_flank_indels =>  0
                };
                push(@parent_keys, $parent);
            }

        }

        if ($entry->{'type'} eq 'up_flank') {
            $all_data{"$parent"}{'up_flank_vars'}++;
            $all_data{"$parent"}{'up_flank_indels'}++ if (length($entry->{'ref'}) != length($alt));
        } elsif ($entry->{'type'} eq 'down_flank') {
            $all_data{"$parent"}{'down_flank_vars'}++;
            $all_data{"$parent"}{'down_flank_indels'}++ if (length($entry->{'ref'}) != length($alt));
        } else {
            $all_data{"$parent"}{'feat_vars'}++;
            $all_data{"$parent"}{'feat_indels'}++ if (length($entry->{'ref'}) != length($alt));
        }
    }

    #Now to print out the output I've built up
    for my $i (@parent_keys) {
        print $output
            $i, "\t",
            $all_data{"$i"}{'feat_vars'}, "\t",
            $all_data{"$i"}{'up_flank_vars'}, "\t",
            $all_data{"$i"}{'down_flank_vars'}, "\t",
            $all_data{"$i"}{'feat_indels'}, "\t",
            $all_data{"$i"}{'up_flank_indels'}, "\t",
            $all_data{"$i"}{'down_flank_indels'}, "\n",
    }

}

sub summarize_depth {
    my %args = @_;
    my $input_gff = $args{'input_gff'};
    my $input_depth = $args{'input_depth'};
    my $output;
    my $parent_name = $args{'parent'};
    my $feature_type = $args{'feature'};
    if ($args{'output'} ne 'NO_OUTPUT_PROVIDED') {
        $output = FileHandle->new("> $args{'output'}");
    }
    else {
        $output = *STDOUT;
    }

    #Add functionality for GTF file:
    my $att_sep = '=';
    if ($input_gff =~ /.gtf$/) {
        $att_sep = ' ';
    }

    my @gff_db;
    my $gff_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $gff_tsv->column_names('seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes');
    my $gff_fh = FileHandle->new("$input_gff");
    open ($gff_fh, "<:encoding(utf8)", "$input_gff");
    while (my $row = $gff_tsv->getline_hr($gff_fh)) {
        next if ($row->{'seq_id'} =~ /^#/);
        next if ($row->{'type'} ne "$feature_type");
        my @tmp = split(/;/, $row->{'attributes'});
        my %attributes;
        for my $item (@tmp) {
            my ($key, $value) = split("$att_sep", $item);
            $attributes{$key} = $value;
        }
        my %hash = (
            seq_id => $row->{'seq_id'},
            start => $row->{'start'},
            end => $row->{'end'},
            parent => $attributes{$parent_name}
        );
        push (@gff_db, \%hash);
    }
    if (scalar(@gff_db) == 0) {
        die "Nothing in GFF!\n";
    }
    #Now that the gff is loaded in, I will begin reading in the depth
    #file.
    my $parent = '_balrog_';    #something that can't be matched
    my $contig = '_balrog_';
    my $depth_fh = FileHandle->new("$input_depth");
    open ($depth_fh, "<:encoding(utf8)", "$input_depth");
    my $depth_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $depth_tsv->column_names('seq_id', 'pos', 'depth');
    my @sub_features;
    my @sub_depth;
    #This probably isn't the best way of doing this but oh well. I've
    #realize that every index of the depth array corresponds to the
    #nucleotide position (minus 1), because the depth file samtools
    #depth produces is always sorted. That should make this much
    #faster
    print $output "#PARENT\tMIN_DEPTH\tMAX_DEPTH\tMEAN_DEPTH\n";
  DEPTH: while (my $row = $depth_tsv->getline_hr($depth_fh)) {
        if ($contig ne $row->{'seq_id'}) {
            #This means I've gotten to the end of a "section" on my
            #depth file
            if ("$contig" eq '_balrog_') {
                $contig = $row->{'seq_id'};
                redo DEPTH;
            }
            else {
                if (scalar(@sub_depth) != 0) {
                    @sub_features = grep { $_->{'seq_id'} eq $contig } @gff_db;
                    if (scalar(@sub_features) == 0) {
                        $contig = $row->{'seq_id'};
                        next DEPTH;
                    }
                    ;
                    @sub_features = sort { $a->{'parent'} cmp $b->{'parent'} } @sub_features;
                    my $parent = $sub_features[0]->{'parent'};
                    my @feature_depths;

                    for my $feat (@sub_features) {
                        if ("$parent" ne $feat->{'parent'}) {
                            print $output "$parent", "\t",
                                min(@feature_depths), "\t",
                                max(@feature_depths), "\t",
                                (sum(@feature_depths) / scalar(@feature_depths)), "\n";
                            $parent = $feat->{'parent'};
                            @feature_depths = ();
                        }
                        @feature_depths = (@feature_depths, @sub_depth[$feat->{'start'} - 1 .. $feat->{'end'} - 1]);
                    }
                    #Catch the final feature which won't be caught by
                    #the above print statement
                    print $output "$parent", "\t",
                        min(@feature_depths), "\t",
                        max(@feature_depths), "\t",
                        (sum(@feature_depths) / scalar(@feature_depths)), "\n";
                    $contig = $row->{'seq_id'};
                }
            }
        }
        else {
            push (@sub_depth, $row->{'depth'});
        }
    }

    #Now I need to catch the final contig
    if (scalar(@sub_depth) != 0) {
        @sub_features = grep { $_->{'seq_id'} eq $contig } @gff_db;
        if (scalar(@sub_features) == 0) {
            exit;
        }
        ;
        @sub_features = sort { $a->{'parent'} cmp $b->{'parent'} } @sub_features;
        my $parent = $sub_features[0]->{'parent'};
        my @feature_depths;

        for my $feat (@sub_features) {
            if ("$parent" ne $feat->{'parent'}) {
                print $output "$parent", "\t",
                    min(@feature_depths), "\t",
                    max(@feature_depths), "\t",
                    (sum(@feature_depths) / scalar(@feature_depths)), "\n";
                $parent = $feat->{'parent'};
                @feature_depths = ();
            }
            @feature_depths = (@feature_depths, @sub_depth[$feat->{'start'} - 1 .. $feat->{'end'} - 1]);
        }
        print $output "$parent", "\t",
            min(@feature_depths), "\t",
            max(@feature_depths), "\t",
            (sum(@feature_depths) / scalar(@feature_depths)), "\n";
    }
}

1;
