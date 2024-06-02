
#Until I have added more functions to this package, I'll keep
#everything here which should make the "use ..." setup much more
#simple.

package Bio::Dantools;
use strict;
use autodie;
use warnings;
use FindBin;
use FileHandle;
use File::Basename;
use File::Path qw( rmtree );
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::CodonTable;
use Bio::DB::Fasta;
use Bio::DB::SeqFeature;
use Bio::SeqFeature::Generic;
use Bio::Matrix::IO;
use Text::CSV_XS::TSV;
use List::Util qw"max min sum";
use File::Copy;
use File::Path qw(make_path);
use POSIX qw"floor ceil";
use lib "$FindBin::Bin/../lib"; #Don't need this now that it's all one file
use Parallel::ForkManager;

#Set up my interrupt trap
BEGIN {
    $SIG{'INT'} = sub {
    exit(1);
};
}

sub pseudogen {
    #This will serve as the main full pseudogen worker
    my %args = @_;
    my $aligner = "$FindBin::Bin/../lib/aligner.sh";
    my $postprocessor = "$FindBin::Bin/../lib/postprocess.sh";
    chdir($args{'outdir'});

    my $it = 0;
    my $itp = -1;

    my $var_count = $args{'min_variants'} + 1;
    my $fai = $args{'fai'};
    my $base_idx = $args{'base_idx'};
    my $og_base_idx = $base_idx;
    my $base = $args{'base'}; #This changes throughout so it needs a variable
    my $og_base = $base;
    my $gff = $args{'gff'};
    my $keepers = $args{'keepers'};
    my $output_name = $args{'output_name'};
    my $threads = $args{'threads'};
    my $var_fraction = $args{'var_fraction'};
    my $outdir = $args{'outdir'};
    my $min_variants = $args{'min_variants'};
    my $scoremin = $args{'scoremin'};
    my $bin_size = $args{'bin_size'};
    my $source = $args{'source'};
    my $readsu = $args{'readsu'};
    my $reads1 = $args{'reads1'};
    my $reads2 = $args{'reads2'};
    my $input_type = $args{'input_type'};
    my $fragment = $args{'fragment'};
    my $lengths = $args{'lengths'};
    my $overlap = $args{'overlap'};
    my $min_length = $args{'min_length'};
    my $input_name;

    if ("$fai" eq '') { $fai = 'NA' }
    if ("$base_idx" eq '') { $base_idx = 'NA' }

    #Begin by fragmenting my reads if a fasta is given
    if (("$input_type" eq 'fasta') & ("$fragment" eq 'yes')) {
        Bio::Dantools::fragment(input => "$source",
                                output => "fragments.fasta",
                                lengths => "$lengths",
                                overlap => "$overlap",
                                min_length => "$min_length",
                                logfile => "fragment_log.txt",
                                threads => "$threads"
            );
        }

    #Creating my metadata file
    my $meta = FileHandle->new("> metadata.tsv");
    if (("$input_type" eq 'fasta')) {
        print $meta "Base\tSource\tIteration\tAlignment\tScale_Alignment\tVariants\n";
    } else {
        print $meta "Base\tSource\tIteration\tAlignment\tVariants\n";
    }

    until (("$var_count" < "$min_variants") & ("$it" != 0)) {
        #Determine which aligner scheme to run based on input type
        if (("$input_type" eq 'fasta') & ("$fragment" eq 'yes')) {
            system("bash", "$aligner",
                   '-b', "$base",
                   '--fai', "$fai",
                   '--indexes', "$base_idx",
                   '-i', "$it",
                   '--output_name', "$output_name",
                   '--source', "fragments.fasta",
                   '-t', "$threads",
                   '--var_fraction', "$var_fraction",
                   '--input_type', "$input_type",
                   '--scoremin', "$scoremin"
                );
            $input_name = basename("$source", ('.fasta'));
        } elsif (("$input_type" eq 'fasta') & ("$fragment" eq 'no')) {
            system("bash", "$aligner",
                   '-b', "$base",
                   '--fai', "$fai",
                   '--indexes', "$base_idx",
                   '-i', "$it",
                   '--output_name', "$output_name",
                   '--source', "$source",
                   '-t', "$threads",
                   '--var_fraction', "$var_fraction",
                   '--input_type', "$input_type",
                   '--scoremin', "$scoremin"
                );
            $input_name = basename("$source", ('.fasta.'));
        } elsif (("$input_type" eq 'fastq_u')) {
            system("bash", "$aligner",
                   '-b', "$base",
                   '--fai', "$fai",
                   '--indexes', "$base_idx",
                   '-i', "$it",
                   '--output_name', "$output_name",
                   '--readsu', "$readsu",
                   '-t', "$threads",
                   '--var_fraction', "$var_fraction",
                   '--input_type', "$input_type",
                   '--scoremin', "$scoremin"
                );
            $input_name = basename("$readsu", ('.fastq'));
        } elsif (("$input_type" eq 'fastq_p')) {
            system("bash", "$aligner",
                   '-b', "$base",
                   '--fai', "$fai",
                   '--indexes', "$base_idx",
                   '-i', "$it",
                   '--output_name', "$output_name",
                   '--reads1', "$reads1",
                   '--reads2', "$reads2",
                   '-t', "$threads",
                   '--var_fraction', "$var_fraction",
                   '--input_type', "$input_type",
                   '--scoremin', "$scoremin"
                );
            $input_name = basename("$reads1", ('.fastq'));
        }

        #The above should have just produced everything I need to
        #continue through with my pipeline. That includes the
        #it"$it"/variants.vcf file, the MAPPING.stderr file, and the
        #new genome. The next step is to shift my bins accordingly
        if ("$it" == 0) {
            Bio::Dantools::bin_maker(
                input_fasta => "$og_base",
                output => "it0/bins.tsv",
                bin_size => "$bin_size"
                );
        } else {
            Bio::Dantools::bin_shifter(
                input_vcf => "it$itp/variants.vcf",
                input_bins => "it$itp/bins.tsv",
                output => "it$it/bins.tsv"
                );
        };

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
        my $aln_scaled = $aln * 2;

        my $vcf = FileHandle->new("< it$it/variants.vcf");
        $var_count = 0;
      COUNT_VARS: while (<$vcf>) {
          next COUNT_VARS if ($_ =~ /^#/);
          $var_count++;
      }
        if (("$input_type" eq 'fasta')) {
            print $meta basename($og_base, ('.fasta')) . "\t" . "$input_name\t" . "$it\t" . "$aln\t" . "$aln_scaled\t" . "$var_count\n";
      } else {
          print $meta basename($og_base, ('.fasta')) . "\t" . "$input_name\t" . "$it\t" . "$aln\t" . "$var_count\n";
      }

        $base_idx = "it$it/indexes/$output_name" . "_it$it";
        if ("$it" != 0) {
            $base = "it$it/$output_name" . "_it$it.fasta";
            $fai = "$base" . ".fai";
        }

        #Determine what to delete, if anything
        if (("$it" != 0)) {
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
        } else {
            copy("it0/bins.tsv", "original_bins.tsv");
        };

        $it++;
        $itp++;
    }

    #Ok and now I need to put the important things together in an
    #output directory
    if (! -d 'output') {
        if ( -e 'output') {
            die "output exists as file, can't overwrite";
        } else {
            mkdir 'output';
            mkdir 'output/indexes';
        }
    } elsif (! -d 'output/indexes') {
        if (-e 'output/indexes') {
            die "output/indexes exists as a file, can't overwrite";
        } else {
            mkdir 'output/indexes';
        }
    }
    copy("$base", "output/$output_name.fasta");
    copy("it$itp/bins.tsv", "output/bins.tsv");
    copy("it$itp/sorted.bam", "output/alignment.bam");

    #Redo the "catching stuff to delete" based on keepers:
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

    #The below is some setup for my needler, which I may just in its
    #own subroutine

    #First I need to break up each genome into its contigs
    Bio::Dantools::contig_split(
        fasta => "$og_base",
        outdir => "ref"
        );
    Bio::Dantools::contig_split(
        fasta => "output/$output_name.fasta",
        outdir => "alt"
        );

    my @contigs;
    my @breaks;
    my @ref_ends;
    my @alt_ends;
    my $prev_contig = '_balrog_'; #Set it to something that won't be matched
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
    my $contigs = join(',', @contigs);
    my $ref_ends = join(',', @ref_ends);
    my $alt_ends = join(',', @alt_ends);
    my $breaks = join(',', @breaks);

    if (("$input_type" eq 'fasta') & ("$fragment" eq 'yes')) {
        system("bash", "$postprocessor",
               "--ref", "$og_base",
               "--alt", "output/$output_name.fasta",
               "--outdir", "$outdir",
               "--output_name", "$output_name",
               "--contigs", "$contigs",
               "--ref_ends", "$ref_ends",
               "--alt_ends", "$alt_ends",
               "--breaks", "$breaks",
               "--aln", "it$it/alignment.bam",
               "--threads", "$threads",
               "--source", "fragments.fasta",
               "--input_type", "$input_type"
            );
    } elsif (("$input_type" eq 'fasta') & ("$fragment" eq 'no')) {
        system("bash", "$postprocessor",
               "--ref", "$og_base",
               "--alt", "output/$output_name.fasta",
               "--outdir", "$outdir",
               "--output_name", "$output_name",
               "--contigs", "$contigs",
               "--ref_ends", "$ref_ends",
               "--alt_ends", "$alt_ends",
               "--breaks", "$breaks",
               "--aln", "it$it/alignment.bam",
               "--threads", "$threads",
               "--source", "$source",
               "--input_type", "$input_type"
            );
    } elsif (("$input_type" eq 'fastq_u')) {
        system("bash", "$postprocessor",
               "--ref", "$og_base",
               "--alt", "output/$output_name.fasta",
               "--outdir", "$outdir",
               "--output_name", "$output_name",
               "--contigs", "$contigs",
               "--ref_ends", "$ref_ends",
               "--alt_ends", "$alt_ends",
               "--breaks", "$breaks",
               "--aln", "it$it/alignment.bam",
               "--threads", "$threads",
               "--readsu", "$readsu",
               "--input_type", "$input_type"
            );
    } elsif (("$input_type" eq 'fastq_p')) {
        system("bash", "$postprocessor",
               "--ref", "$og_base",
               "--alt", "output/$output_name.fasta",
               "--outdir", "$outdir",
               "--output_name", "$output_name",
               "--contigs", "$contigs",
               "--ref_ends", "$ref_ends",
               "--alt_ends", "$alt_ends",
               "--breaks", "$breaks",
               "--aln", "it$it/alignment.bam",
               "--threads", "$threads",
               "--reads1", "$reads1",
               "--reads2", "$reads2",
               "--input_type", "$input_type"
            );
    }

    my $file_idx = 0;
    if (-e 'needle_files/combined_bases.tsv') {
        unlink('needle_files/combined_bases.tsv');
    }
    my $out = FileHandle->new(">> needle_files/combined_bases.tsv");
    while ("$file_idx" < scalar @contigs) {
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
            } elsif ($seq->id eq 'alt') {
                $altseq = $seq->seq;
                @altseq = split //, $altseq;
            }
        };

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
        sample => "${output_name}"
        );
    if ("$keepers" ne "all") {
        rmtree(['ref', 'alt', 'needle_files', 'original_indexes'], 0, 0);
        unlink(('original_bins.tsv',
                'original_base.fasta.fai',
                'fragments.fasta',
               'fragment_log.txt'));
    }

    #This used to automatically generate a new GFF and label the vcf,
    #but I don't feel that is necessary for the general pseudogen. If
    #people want to, they can always call dantools label and dantools
    #shift separately
};


#The below can break my genome into pieces
sub fragment {
    my %args = @_;
    my $seqio_out;
    if ($args{'output'} ne 'NO_OUTPUT_PROVIDED') {
        $seqio_out = Bio::SeqIO->new(-format => 'fasta', -file => "> $args{'output'}");
    } else {
        $seqio_out = Bio::SeqIO->new(-format => 'fasta', -fh => \*STDOUT);
    };

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
        } #End iterating over every chunk
          print $log "Wrote $num_written fragments of size $size for $id. \n";
      } #End iterating over every size
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
        if($end >= $length) {
            print $out "$contig" . "\t" . "$length\n";
            last BINS;
        } else {
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

    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'metadata');

    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if $row->{'contig'} =~ /^#/;
        my $shift = length($row->{'alternate'}) - length($row->{'reference'});
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

    #Importantly I don't sort this. That's because if I do, it's going
    #to sort it in a different way from my original. In dantools
    #pseudogen, the vcf and bins should always be sorted the same anyway
    my $contig;
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
    if(scalar @vcf == 0) {
        copy("$args{'input_bins'}", "$args{'output'}");
    } else {
      BINS: while(my $entry = $bins_tsv->getline_hr($bins_fh)) {
          $contig = $entry->{'contig'};
          if (! grep /$contig/, @contigs) {
              $shiftsum = 0;
              push @contigs, $contig;
          }
          my $end = $entry->{'end'};

        SHIFT: while ($vcf_row = $vcf["$vcf_index"]) {
            $vcf_contig = $vcf_row->{'seq_id'};
            last SHIFT if(($end < $vcf_row->{'pos'}) | ($contig ne $vcf_contig));
            $shift = $vcf_row->{'shift'};
            $shiftsum = $shiftsum + $shift;
            $vcf_index++;
        }
          #Some logic to sheck if deletion covers the feature
          $prev_row = $vcf["$vcf_index" - 1];
          $shift = $prev_row->{'shift'};
          if (($prev_row->{'pos'} - $shift > $end) & ($prev_row->{'seq_id'} eq $contig) & ($vcf_index != 0)) {
              print $out "$contig" . "\t" . ($prev_row->{'pos'} + $shiftsum - $shift) . "\n";
          } else {
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

    mkdir $outdir  if (defined($outdir));
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
        mkdir("$tmpdir");
    };

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
    };
    #Note that I don't have any real metadata yet, but I may in the future
    my $header = "##fileformat=VCFv4.2\n";
    print $headfile $header;

    my $raw_fh = FileHandle->new("$input");
    my $raw_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($raw_fh, "<:encoding(utf8)", "$input") or die "Can't open input file";
    $raw_tsv->column_names('chromosome', 'ref', 'alt');

    #The things I'll be building
    my $ref = ''; #String to build ref column of vcf
    my $alt = ''; #String to build alt column of vcf
    my $chrom = '_balrog_';
    my $relpos = 1; #track position on reference chromosome/contig
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
      if($chrom ne $row->{'chromosome'}){
          if ($chrom ne '_balrog_') {
              $length = $relpos - 1;
              print $headfile "##contig=<ID=" . "$chrom" . ",length=" . "$length" . ">\n";
          };
          $relpos = 1; #Reset relative chromosome position
          $chrom = $row->{'chromosome'};
      }

      my $rowref = $row->{'ref'};
      my $rowalt = $row->{'alt'};
      if($rowref eq $rowalt) {
          #I need to adjust using this catch-all
          if((length($ref) != 0) & ("$ref" ne "$alt")) {
              $pos = $relpos - $ref =~ tr/-//c;
              $ref =~ tr/-//d;
              $alt =~ tr/-//d;
              #I need to double check after removing '-'
              if("$ref" ne "$alt") {
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
      } elsif ($rowalt eq '-') {
          $ref = "$ref" . "$rowref";
          $alt = "$alt" . "$rowalt";
          $relpos++;
          next VARS;
      }
      #Make the logic increment depth properly

      #If neither of the above were evaluated, this must be an SNP:
      if((length($ref) != 0) & ("$ref" ne "$alt")) {
          #I'm sure there is a more efficient way of doing all
          #this. Right now I'm checking that ref ne alt twice. This is
          #because for shifting, ref and alt need to contain - for
          #indel. However, after deleting them the ref and alt may by
          #chance become the same, requiring me to check again
          $pos = $relpos - $ref =~ tr/-//c;
          $ref =~ tr/-//d;
          $alt =~ tr/-//d;
          if("$ref" ne "$alt") {
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
    };
    close($headfile);
    unlink("$tmpdir/tmp_variants.vcf");
}


#This seeks to shift feature annotations according to the indels
#within a vcf file. To this end, it is passed an input gff, an output
#file, and an input vcf. Right now it is compatible with gff only, but
#it shouldn't be too difficult to modify it to accept other file
#types.

sub gff_shifter {
    my %args = @_;
    my $starts = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $gff_fh = FileHandle->new("$args{'gff'}");
    my $gff_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($gff_fh, "<:encoding(utf8)", "$args{'gff'}");
    $gff_tsv->column_names('contig', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes');

    my $feat;
    my $gff_id = 0;
    my @heading;
  READ_GFF: while (my $row = $gff_tsv->getline_hr($gff_fh)) {
      if($row->{'contig'} =~ /^#/) {
          push @heading, $row->{'contig'};
          next READ_GFF;
      };
      my $score;
      if($row->{'score'} eq '.') {
          $score = 123456;
      };
      #For shifting logic, I will create two structures and sort them independently
      $feat = Bio::SeqFeature::Generic->new(
          -seq_id => $row->{'contig'},
          -start => $row->{'start'},
          -end => $row->{'end'},
          -strand => $row->{'strand'},
          -primary_tag => $row->{'type'},
          -score => $score,
          -source_tag => $row->{'source'},
          -tag => {phase => $row->{'phase'},
                   attributes => $row->{'attributes'},
                   gff_id => $gff_id
          });
      $starts->store($feat);
      $gff_id++;
  }
    close($gff_fh);
    #Sorting each DB
    my @starts_db = $starts->features;
    my @ends_db = @starts_db; #These are identical until I sort them:
    @starts_db = sort {
        $a->seq_id cmp $b->seq_id || $a->start <=> $b->start
    } @starts_db;
    @ends_db = sort {
        $a->seq_id cmp $b->seq_id || $a->end <=> $b->end
    } @ends_db;

    my $vcf_fh = FileHandle->new("$args{'vcf'}");
    my $vcf = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($vcf_fh, "<:encoding(utf8)", "$args{'vcf'}") or die "Can't open VCF file";
    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'info');

    #Load in my VCF features and add them to a DB:
    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if ($row->{'contig'} =~ /^#/);
        next if (! defined($row->{'alternate'}));
        my $shift = length($row->{'alternate'}) - length($row->{'reference'});
        next if $shift == 0;
        my $feat = Bio::SeqFeature::Generic->new(
            -seq_id => $row->{'contig'},
            -start => $row->{'position'},
            -end => $row->{'position'},
            -strand => 1,
            -tag => { ref => $row->{'reference'},
                      alt => $row->{'alernate'},
                      shift => $shift
            }
            );
        $vcf->store($feat);
    }
    close($vcf_fh);

    my @vcf_db = $vcf->features;
    @vcf_db = sort {
        $a->seq_id cmp $b->seq_id || $a->start <=> $b->start
    } @vcf_db;

    my $gff_out;
    if ("$args{'output'}" eq 'NO_OUTPUT_PROVIDED') {
        $gff_out = *STDOUT;
    } else {
        $gff_out = FileHandle->new("> $args{'output'}");
    }

    if (scalar @vcf_db == 0) {
        print STDERR "No indels in the VCF file, copying\n";
        $gff_fh = FileHandle->new("$args{'gff'}");
        while (<$gff_fh>) { print $gff_out $_ };
        exit(1);
    };

    my $contig;
    my $idx = 0;
  HEADER: foreach my $line (@heading) {
      my @shift;
      if (index($line, 'sequence-region') == -1) {
          print $gff_out $line, "\n";
          next HEADER;
      } else {
          $line =~ tr /#//d;
          my @line = split / /, $line;
          my $length = $line[3];
          $contig = $line[1];
          while (my $row = $vcf_db["$idx"]) {
              last if ($contig ne $row->seq_id);
              @shift = $row->get_tag_values('shift');
              $length = $length + $shift[0];
              $idx++;
          };
          print $gff_out "##sequence-region " . "$contig " . "1 " . "$length\n";
      };
  };

    my $gff_row;
    my $gff_index = 0;
    my $shifted_starts = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $shifted_ends = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');

    my $shiftsum = 0;
    my @shift = 0;
    my $varpos = 0;
    my $var_contig = '';
    my @seen_contigs = '';
    my $vcf_row;
    my $feat_start;
    $idx = 0;
  STARTS: while ($gff_row = $starts_db["$gff_index"]) {
      $contig = $gff_row->seq_id;
      $feat_start = $gff_row->start;
      if (! grep /$contig/, @seen_contigs) {
          $shiftsum = 0;
          push @seen_contigs, $contig;
          @shift = 0;
      };
    START_SHIFT: while ($vcf_row = $vcf_db["$idx"]) {
        my $tmp_contig = $vcf_row->seq_id;
        if ($contig ne $tmp_contig) {
            if (! grep /$tmp_contig/, @seen_contigs) {
                #Means this variant is on the next contig
                last START_SHIFT;
            } else {
                #Means variant is on previous chromosome
                $idx++;
                next START_SHIFT;
            }
        }
        if ($feat_start < $vcf_row->start) {
            last START_SHIFT;
        }
        #If we get here, it means the variant is on the right
        #chromosome and upstream of our feature
        @shift = $vcf_row->get_tag_values('shift');
        $varpos = $vcf_row->start; #kept for later logic
        $shiftsum = $shiftsum + $shift[0];
        $var_contig = $vcf_row->seq_id;
        $idx++;
    };
      #Now I should have collected the proper shiftsum:
      if (($varpos - $shift[0] > $feat_start) & ($var_contig eq $contig) & ($idx != 0)) {
          $gff_row->start($varpos + $shiftsum - $shift[0]);
      } else {
          $gff_row->start($feat_start + $shiftsum);
      };
      $shifted_starts->store($gff_row);
      $gff_index++;
  };

    #Now I need to repeat this with my end positions, which requires
    #resetting my variables:
    $shiftsum = 0;
    @shift = 0;
    $varpos = 0;
    $var_contig = '';
    @seen_contigs = '';
    $gff_index = 0;
    $idx = 0;
    my $feat_end;
  ENDS: while ($gff_row = $ends_db["$gff_index"]) {
      $contig = $gff_row->seq_id;
      $feat_end = $gff_row->end;
      if (! grep /$contig/, @seen_contigs) {
          $shiftsum = 0;
          push @seen_contigs, $contig;
          @shift = 0;
      };
    END_SHIFT: while ($vcf_row = $vcf_db["$idx"]) {
        my $tmp_contig = $vcf_row->seq_id;
        if ($contig ne $tmp_contig) {
            if (! grep /$tmp_contig/, @seen_contigs) {
                #Means this variant is on the next contig
                last END_SHIFT;
            } else {
                #Means variant is on previous chromosome
                $idx++;
                next END_SHIFT;
            }
        }
        if ($feat_end < $vcf_row->start) {
            last END_SHIFT;
        }
        #If we get here, it means the variant is on the right
        #chromosome and upstream of our feature
        @shift = $vcf_row->get_tag_values('shift');
        $varpos = $vcf_row->start; #kept for later logic
        $shiftsum = $shiftsum + $shift[0];
        $var_contig = $vcf_row->seq_id;
        $idx++;
    };
      #Now I should have collected the proper shiftsum:
      if (($varpos - $shift[0] > $feat_end) & ($var_contig eq $contig) & ($idx != 0)) {
          $gff_row->end($varpos + $shiftsum - $shift[0]);
      } else {
          $gff_row->end($feat_end + $shiftsum);
      };
      $shifted_ends->store($gff_row);
      $gff_index++;
  };

    #Now I need to sort each database according to its tag value
    #"gff_id":
    @starts_db = sort {
        my @pos1 = $a->get_tag_values('gff_id');
        my @pos2 = $b->get_tag_values('gff_id');
        if ($pos1[0] < $pos2[0]) { return -1; }
        elsif ($pos1[0] == $pos2[0]) { return 0; }
        else { return 1; };
    } @starts_db;
    @ends_db = sort {
        my @pos1 = $a->get_tag_values('gff_id');
        my @pos2 = $b->get_tag_values('gff_id');
        if ($pos1[0] < $pos2[0]) { return -1; }
        elsif ($pos1[0] == $pos2[0]) { return 0; }
        else { return 1; };
    } @ends_db;

    #Now these should each be back to their original sorting and
    #identical in all but positions:
    $idx = 0;
    while ($gff_row = $starts_db["$idx"]) {
        my @phase = $gff_row->get_tag_values('phase');
        my $phase = $phase[0];
        my @attributes = $gff_row->get_tag_values('attributes');
        my $attributes = $attributes[0];
        my $strand = $gff_row->strand;
        my $score = $gff_row->score;
        my $tmp_row = $ends_db["$idx"];
        if ($strand == 1) {
            $strand = '+';
        } else {
            $strand = '-';
        }
        if ($score == 123456) {
            $score = '.';
        }
        my $start = $gff_row->start;
        my $end = $ends_db["$idx"]->end;
        print $gff_out $gff_row->seq_id . "\t" .
            $gff_row->source_tag . "\t" .
            $gff_row->primary_tag . "\t" .
            $start . "\t" .
            $end . "\t" .
            $score . "\t" .
            $strand . "\t" .
            $phase . "\t" .
            $attributes . "\n";
        $idx++;
    }
};

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

    chdir($args{'outdir'});

    #Modify my features if flanks are added:
    if ("$add_flanks" eq 'yes') {
        push(@feature_types, ('up_flank', 'down_flank'));
    };
    #Develop method to determine which features to add
    my %feature_types = map {$_ => 1 } @feature_types;
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
    while(my $row = $gff_tsv->getline_hr($gff_fh)) {
        next if($row->{'seq_id'} =~ /^#/);
        my $feat_type = $row->{'type'};
        #next if(! exists($feature_types{$feat_type}));
        if (exists($feature_types{$feat_type})) {
            my @tmp = split(/;/, $row->{'attributes'});
            my %attributes;
            for my $item (@tmp) {
                my ($key, $value) = split('=', $item);
                $attributes{$key} = $value;
            };
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
            };
        };
        #Now I need to figure out whether to add some flanks. This is
        #separate because in some instances you may want to add flanks
        #to the edges of the protein_coding_gene, and annotate
        #variants by CDS

        if (("$feat_type" eq "$flank_feature") & ("$add_flanks" eq 'yes')) {
            my @tmp = split(/;/, $row->{'attributes'});
            my %attributes;
            for my $item (@tmp) {
                my ($key, $value) = split('=', $item);
                $attributes{$key} = $value;
            };
            my $parent = $attributes{$flank_parent};
            my $child = $attributes{$child_name};
            if ($row->{'strand'} eq '-') {
                if ($flank_lengths[0] != 0) {
                    my %down = (seq_id => $row->{'seq_id'},
                             start => $row->{'start'} - $flank_lengths[1],
                             end => $row->{'start'} - 1,
                             strand => '-',
                             type => 'down_flank',
                             parent => $parent,
                             child => "down_flank_" . $child
                        );
                    push(@gff, \%down);
                }
                if ($flank_lengths[1] != 0) {
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
            } else {
                #I think any good gff has only + or - strands, but
                #this is general so it will treat anything without "-"
                #as a plus strand item
                if ($flank_lengths[0] != 0) {
                    my %up = (seq_id => $row->{'seq_id'},
                             start => $row->{'start'} - $flank_lengths[0],
                             end => $row->{'start'} - 1,
                             strand => $row->{'strand'},
                             type => 'up_flank',
                             parent => $parent,
                             child => "up_flank_" . $child
                        );
                    push(@gff, \%up);
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
            };
        };
    };

    #I'm going to break this into two GFF files
    my @plus_features = grep { $_->{'strand'} eq '+' } @gff;
    my @minus_features = grep { $_->{'strand'} eq '-' } @gff;

    #Now that I've loaded my GFF I want to load my VCF.
    my $vcf_fh = FileHandle->new("$input_vcf");
    #my $vcf = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my @vcf;
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $vcf_tsv->column_names('seq_id', 'pos', 'index', 'ref', 'alt', 'qual', 'filter', 'info');
    open($vcf_fh, "<:encoding(utf8)", "$input_vcf") or die "Can't open vcf_fh";

    while(my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if $row->{'seq_id'} =~ /^#/;
        next if (! defined($row->{'info'}));
        my @tmp = split /;/, $row->{'info'};
        my %info;
        for my $item (@tmp) {
            my ($key, $value) = split('=', $item);
            $info{$key} = $value;
        };
        my %var = (
            seq_id => $row->{'seq_id'},
            pos => $row->{'pos'},
            strand => 1,
            ref => $row->{'ref'},
            alt => $row->{'alt'},
            depth => $info{'DP'}
            );
        push(@vcf, \%var);
    };

    #Now I need to determine all of the unique contigs I have:
    my @unique_contigs;
    for my $feat (@vcf) {
        push(@unique_contigs, $feat->{'seq_id'});
    };
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
    my $death = 0; #Determines if I "die" or not
  CONTIGS: for my $contig (@unique_contigs) {
        my $pid = $pm->start and next CONTIGS;
      my @plus_sub_features = grep { $_->{'seq_id'} =~ /$contig/i } @plus_features;
      my @variants = grep { $_->{'seq_id'} =~ /$contig/i } @vcf;
      if (scalar(@variants) == 0) {
          $pm->finish;
      };
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
        my @minus_sub_features = grep { $_->{'seq_id'} =~ /$contig/i } @minus_features;
      @minus_sub_features = sort {
          $b->{'parent'} cmp $a->{'parent'} ||
              $b->{'start'} <=> $a->{'start'}
          } @minus_sub_features;


        #Now I need to iterate over each feature, keeping track of the parent
        my @labels; #Instead of printing everytime, which can cause
      #overlap in prints, save info to an array of hashes
                  #and print at the end
        if (scalar(@plus_sub_features) !=  0) {
      my $parent = $plus_sub_features[0]->{'parent'}; #Start this up
      my $parent_pos = 0; #track number of bases covered in parent
      my $coding_pos = 0; #track number of bases covered in coding sequence
      my $var_idx = 0;
      my $feat_end = 0; #make sure it doesn't fail in first iteration
      my $feat_start;
      my $var = undef;
    FEATURES: for my $feat (@plus_sub_features) {
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
            } else {
                $var_idx = 0;
            };
            last FEATURES if (! defined($var));
        } elsif ($feat_end > $feat->{'start'}) {
            #This means we are in the same parent and have overlap, not good
                $death = 1;
                print "ERROR: Features overlap within same parent\n";
                last CONTIGS;
        };
        $feat_start = $feat->{'start'};
        $feat_end = $feat->{'end'};
      VARIANTS: while ($var = $variants[$var_idx]) {
          my $var_pos = $var->{'pos'};
          if ($var_pos < $feat_start ) {
              $var_idx++;
              next VARIANTS;
          } elsif ($var_pos > $feat_end) {
              $parent_pos = $parent_pos + $feat_end - $feat_start + 1;
              if ($feat->{'type'} eq "$coding_feature") {
                  $coding_pos = $coding_pos + $feat_end - $feat_start + 1;
              }
              next FEATURES;
          };
          my $relpos = $var_pos - $feat_start + 1;
          my $cd;
          if ($feat->{'type'} eq "$coding_feature") {
              $cd = $relpos + $coding_pos;
          } else {
              $cd = 'NA';
          };
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
          $var_idx++;
      };
    };
  }
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
      my $var_idx = 0;
      my $var = undef;
        my $feat_start = 100000000; #make sure it doesn't fail on first iteration
        my $feat_end;
    FEATURES: for my $feat (@minus_sub_features) {
        #For whatever reason, the above sorting can create an empty
        #entry in my hash, so I need to add this caveat
        next FEATURES if(! defined($feat->{'parent'}));
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
            } else {
                $var_idx = 0;
            };
            last FEATURES if (! defined($var));
        } elsif ($feat_start < $feat->{'end'}) {
            print "ERROR: Features overlap within same parent\n";
                $death = 1;
                last CONTIGS;
        };
        $feat_start = $feat->{'start'};
        $feat_end = $feat->{'end'};
      VARIANTS: while ($var = $variants[$var_idx]) {
          my $var_pos = $var->{'pos'}; #Some extra logic for - strand
          my $refn = $var->{'ref'};
          $var_pos = $var_pos + length($refn) - 1; #Change position for deletions
          if ($var_pos > $feat_end) {
              $var_idx++;
              next VARIANTS;
          } elsif ($var_pos < $feat_start) {
              $parent_pos = $parent_pos + $feat_end - $feat_start + 1;
              if ($feat->{'type'} eq "$coding_feature") {
                  $coding_pos = $coding_pos + $feat_end - $feat_start + 1;
              };
              next FEATURES;
          };
          $refn =~ tr/ACGT/TGCA/;
          my $altn = $var->{'alt'};
          $altn =~ tr/ACGT/TGCA/;
          my $relpos = $feat_end - $var_pos + 1;
          my $cd;
          if ($feat->{'type'} eq "$coding_feature") {
              $cd = $relpos + $coding_pos;
          } else {
              $cd = 'NA';
          };
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
          $var_idx++;
      };
    };
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
      };
      close($tmp_fh);
  }
      $pm->finish;
  };
    if ("$death" == 1) { return };
    $pm->wait_all_children;
    #Combine all of the temporary files:
    for my $contig (@unique_contigs) {
        if (-e "${tmpdir}/${contig}") {
            my $fh = FileHandle->new("< ${tmpdir}/${contig}");
        while (<$fh>) {
            print $outnuc "$_";
        };
            close($fh);
    };
    };
    close($outnuc);
    rmtree("$tmpdir");

    if ("$translate" ne 'yes') { exit 0 };
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
        } else {
            $sequence = ($sequence) . ($fasta->seq($feat->{'seq_id'}, $feat->{'start'}, $feat->{'end'}));
        };
    };
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
        } else {
            $sequence = $sequence . $fasta->seq($feat->{'seq_id'}, $feat->{'start'}, $feat->{'end'});
        };
    };
    #Now I need to catch the output final feature
    $feature_seqs{$parent} = $sequence;

    #Now I only need to read in this file labeled by CDS
    my @labeled_vars;
    my $labeled_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    my $labeled_fh = FileHandle->new("${output_nuc}");
    $labeled_tsv->column_names('seq_id', 'pos', 'ref', 'alt', 'parent', 'child', 'type', 'strand', 'rel_parent', 'rel_child', 'rel_coding', 'depth');
    while (my $row = $labeled_tsv->getline_hr($labeled_fh)) {
        next if($row->{'seq_id'} =~ /^#/);
        next if($row->{'type'} ne "$coding_feature"); #Catch potentially empty line
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
    };
    @labeled_vars = sort {
        $a->{'parent'} cmp $b->{'parent'} ||
            $a->{'rel_coding'} <=> $b->{'rel_coding'}
    } @labeled_vars;
    my $var_idx = 0;
    my $ref_codon = '';
    my $alt_codon;
    my $codon_num = 'FOOBAR';
    my @nvar;
    $parent = 'FOOBAR';
    my $complex_vars = 0;

    #Define my codon table, I need to implement an option separate for
    #all the different ID values.
    $codon_table = Bio::Tools::CodonTable->new( -id => "$codon_table");

    my $outaa = FileHandle->new("> ${output_aa}");
    print $outaa "#CHROM\tPARENT\tCODON\tREF\tALT\tVARS\tTYPE\tFRAME\tSCORE\tSCALED_SCORE\n";
  CODONS: while (my $var = $labeled_vars[$var_idx]) {
      #If I'm still on the same codon:
        if (($codon_num eq $var->{'codon'}) & ($parent eq $var->{'parent'})) {
            substr($alt_codon, $var->{'codon_pos'} - 1, 1) = $var->{'alt'};
            @nvar = (@nvar, $var->{'ref'} . $var->{'rel_coding'} . $var->{'alt'});
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
                };
                my $type = 'point';
                if (length($ref_aa) != 1) {
                    $type = 'del';
                };
                print $outaa
                    $var->{'seq_id'}, "\t",
                    $var->{'parent'}, "\t",
                    $codon_num, "\t",
                    $ref_aa, "\t",
                    $alt_aa, "\t",
                    join(',', @nvar), "\t",
                    $type, "\t",
                    "in", "\t",
                    $score, "\t",
                    $scaled_score, "\n";

            } elsif ((length($alt_codon) % 3 == 0) & (length($ref_codon) == 3)) {
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
                };
                print $outaa
                    $var->{'seq_id'}, "\t",
                    $var->{'parent'}, "\t",
                    $codon_num, "\t",
                    $ref_aa, "\t",
                    $alt_aa, "\t",
                    join(',', @nvar), "\t",
                    $type, "\t",
                    "in", "\t",
                    $score, "\t",
                    $scaled_score, "\n";
            } else {
                #Implies a frameshift
                my $ref_aa = $codon_table->translate(substr($ref_codon, 0, 3));
                my $alt_aa = $codon_table->translate(substr($alt_codon, 0, 3));
                if (length($ref_codon) > length($alt_codon)) {
                    print $outaa
                        $var->{'seq_id'}, "\t",
                        $var->{'parent'}, "\t",
                        $codon_num, "\t",
                        $ref_aa . '@', "\t",
                        $alt_aa, "\t",
                        join(',', @nvar), "\t",
                        'del', "\t",
                        'out', "\t",
                        "NA", "\t",
                        "NA", "\n";
                } elsif(length($ref_codon) < length($alt_codon)) {
                    print $outaa
                        $var->{'seq_id'}, "\t",
                        $parent, "\t",
                        $codon_num, "\t",
                        $ref_aa, "\t",
                        $alt_aa . '@', "\t",
                        join(',', @nvar), "\t",
                        'ins', "\t",
                        'out', "\t",
                        "NA", "\t",
                        "NA", "\n";
                } else {
                    $complex_vars++;
                };
            };
        };
        $codon_num = $var->{'codon'};
        $parent = $var->{'parent'};
        $ref_codon = substr($feature_seqs{$parent}, (($codon_num * 3) - 3), 3);
        $alt_codon = $ref_codon;
        substr($ref_codon, $var->{'codon_pos'} - 1, 1) = $var->{'ref'}; #necessary for complex
        substr($alt_codon, $var->{'codon_pos'} - 1, 1) = $var->{'alt'};
        @nvar = ($var->{'ref'} . $var->{'rel_coding'} . $var->{'alt'});

        $var_idx++;
  };
    if ($complex_vars != 0) {
        print "WARNING: ${complex_vars} complex variants encountered, omitted\n";
    };
};

#Ok and now for what might be the final tool of dantools before
#publication, dantools summary. This will take either the nucleotide
#or amino acid input from dantools translate's output and create a
#summary about it. This includes things like total mutations, change
#in length, number of frameshifts (if any), etc.

sub summarize_aa {
    my %args = @_;
    my $input_gff = $args{'input_gff'};
    my $input_depth = $args{'input_depth'};
    my $input_vars = $args{'input_vars'};
    my $input_vars_fh = FileHandle->new("$input_vars");
    my $output;
    my $parent_name = $args{'parent'};
    my $feature_type = $args{'feature'};
    my $outscore = $args{'outscore'};
    if ($args{'output'} ne 'NO_OUTPUT_PROVIDED') {
        $output = FileHandle->new("> $args{'output'}");
    } else {
        $output = *STDOUT;
    };

    my @vars_db;
    my $vars_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $vars_tsv->column_names('seq_id', 'parent', 'codon', 'ref', 'alt', 'vars', 'type', 'frame', 'score', 'scaled_score');
    while (my $row = $vars_tsv->getline_hr($input_vars_fh)) {
        next if $row->{'seq_id'} =~ /^#/;
        my %hash = (
            seq_id => $row->{'seq_id'},
            parent => $row->{'parent'},
            ref => $row->{'ref'},
            alt => $row->{'alt'},
            type => $row->{'type'},
            frame => $row->{'frame'},
            score => $row->{'score'},
            scaled_score => $row->{'scaled_score'}
        );
        push (@vars_db, \%hash);
    };
    my $total_score = 0;
    my $total_sscore = 0;
    my $num_vars = 0;
    my $num_indels = 0;
    my $num_outframe = 0;
    my $parent = $vars_db[0]->{'parent'};
    my @feat_vars;

    print $output "#PARENT\tNUM_VARS\tSCORE\tSCALED_SCORE\tINDELS\tNUM_INFRAME\tNUM_OUTFRAME\n";

    VARS: for my $entry (@vars_db) {
        if ("$parent" ne $entry->{'parent'}) {
            print $output
                "$parent", "\t",
                "$num_vars", "\t",
                "$total_score", "\t",
                "$total_sscore", "\t",
                "$num_indels", "\t",
                ($num_indels - $num_outframe), "\t",
                "$num_outframe", "\n";

        $total_score = 0;
        $total_sscore = 0;
        $num_vars = 0;
        $num_indels = 0;
        $num_outframe = 0;
        $parent = $entry->{'parent'};
    };
    $num_vars++;
        if ($entry->{'type'} ne 'point') {
            $num_indels++;
            if ($entry->{'frame'} eq 'out') {
                $num_outframe++;
                $total_score = $total_score + $outscore;
                $total_sscore = $total_sscore + $outscore;
                next VARS;
            };
        }
        $total_score = $total_score + $entry->{'score'};
        $total_sscore = $total_sscore + $entry->{'scaled_score'};
    };
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
    } else {
        $output = *STDOUT;
    };

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
            my ($key, $value) = split('=', $item);
            $attributes{$key} = $value;
        };
        my %hash = (
            seq_id => $row->{'seq_id'},
            start => $row->{'start'},
            end => $row->{'end'},
            parent => $attributes{$parent_name}
            );
        push (@gff_db, \%hash);
    };
    if (scalar(@gff_db) == 0) { die "Nothing in GFF!\n" };
    #Now that the gff is loaded in, I will begin reading in the depth
    #file.
    my $parent = '_balrog_'; #something that can't be matched
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
            } else {
                if (scalar(@sub_depth) != 0) {
                    @sub_features = grep { $_->{'seq_id'} eq $contig } @gff_db;
                    if (scalar(@sub_features) == 0) {
                        $contig = $row->{'seq_id'};
                        next DEPTH;
                    };
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
                        };
                        @feature_depths = (@feature_depths, @sub_depth[$feat->{'start'} - 1 .. $feat->{'end'} - 1]);
                    };
                    #Catch the final feature which won't be caught by
                    #the above print statement
                    print $output "$parent", "\t",
                        min(@feature_depths), "\t",
                        max(@feature_depths), "\t",
                        (sum(@feature_depths) / scalar(@feature_depths)), "\n";
                    $contig = $row->{'seq_id'};
                };
            }
        } else {
            push (@sub_depth, $row->{'depth'});
        };
    };

    #Now I need to catch the final contig
    if (scalar(@sub_depth) != 0) {
        @sub_features = grep { $_->{'seq_id'} eq $contig } @gff_db;
        if (scalar(@sub_features) == 0) {
            exit;
        };
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
            };
            @feature_depths = (@feature_depths, @sub_depth[$feat->{'start'} - 1 .. $feat->{'end'} - 1]);
        };
        print $output "$parent", "\t",
            min(@feature_depths), "\t",
            max(@feature_depths), "\t",
            (sum(@feature_depths) / scalar(@feature_depths)), "\n";
    };
};

1;
