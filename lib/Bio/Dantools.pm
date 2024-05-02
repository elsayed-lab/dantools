
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
use Bio::DB::SeqFeature;
use Bio::SeqFeature::Generic;
use Text::CSV_XS::TSV;
use File::Copy;
use POSIX qw"floor ceil";
use lib "$FindBin::Bin/../lib"; #Don't need this now that it's all one file
use Parallel::ForkManager;
use Data::Dumper;

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
                                output_name => "$output_name",
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
                   '--input_type', "$input_type"
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
                   '--input_type', "$input_type"
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
                   '--input_type', "$input_type"
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
                   '--input_type', "$input_type"
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
                output => "it$it/bins.tsv",
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

    #Now I should have some important items like my genome, indexes,
    #and alignment in the output directory. Now all I need to do is
    #take all the files in needle_files/ and put them together in a
    #comprehensive "combined_bases.tsv" file

    my $file_idx = 0;
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
            } elsif ($seq-> id eq 'alt') {
                $altseq = $seq->seq;
                @altseq = split //, $altseq;
            }
        };

        my $idx = 0;
        my $lim = scalar @refseq;
        my $out = FileHandle->new(">> needle_files/combined_bases.tsv");

        while ($idx < $lim) {
            print $out "$contigs[$file_idx]\t$refseq[$idx]\t$altseq[$idx]\n";
            $idx++;
        }
        $file_idx++;
    }

    #Now I have a list of every base from each needle alignment next
    #to what it was in the alignment. I could and probably should make
    #this more efficient by skipping the printing to a file and just
    #building a hash of arrays.
    Bio::Dantools::vcf_maker(
        input => "needle_files/combined_bases.tsv",
        output => "output/$output_name.vcf",
        tmpdir => "output",
        depths => "output/depth.tsv"
        );
    if ("$keepers" ne "all") {
        rmtree(['ref', 'alt', 'needle_files', 'original_indexes'], 0, 0);
        unlink(('original_bins.tsv',
                'original_base.fasta.fai',
                'fragments.fasta',
               'fragment_log.txt'));
    }

    #Now the second to last step, labeling the VCF file:
    if ("$gff" ne '') {
        Bio::Dantools::vcf_labeler(
            gff => "$args{'gff'}",
            vcf => "output/$output_name.vcf",
            output => "output/labeled_variants.tsv",
            add_flanks => "$args{'add_flanks'}",
            flank_lengths => "$args{'flank_lengths'}",
            feature_type => "$args{'feature_type'}",
            feature_name => "$args{'feature_name'}"
            );

    #Ok so now I'm finally almost done. I just need one last step, and
    #that is to take my original gff and my vcf file and modiy the
    #former according to the latter
    Bio::Dantools::gff_shifter(
        input_gff => "$args{'gff'}",
        input_vcf => "output/$args{'output_name'}.vcf",
        output => "output/$args{'output_name'}.gff"
        );
    };
};


#The below can break my genome into pieces
sub fragment {
    my %args = @_;

    my $seqio_in = Bio::SeqIO->new(-format => 'fasta', -file => $args{'input'});
    my $seqio_out = Bio::SeqIO->new(-format => 'fasta', -file => "> $args{'output'}");
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
    my $vcf = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($vcf_fh, "<:encoding(utf8)", "$args{'input_vcf'}") or die "Can't open vcf_fh";

    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'metadata');

    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if $row->{'contig'} =~ /^#/;
        my $shift = length($row->{'alternate'}) - length($row->{'reference'});
        next if $shift == 0;
        my $feat = Bio::SeqFeature::Generic->new(
            -seq_id => $row->{'contig'},
            -start => $row->{'position'},
            -end => $row->{'position'},
            -tag => { ref => $row->{'reference'},
                      alt => $row->{'alernate'},
                      shift => $shift
            }
            );
        $vcf->store($feat);
    }
    close($vcf_fh);

    my @db = $vcf->features;
    @db = sort {
        $a->seq_id cmp $b->seq_id || $a->start <=> $b->start
    } @db;

    my $contig;
    my @contigs;
    my $shiftsum;
    my $vcf_index = 0;
    my @shift = 0;
    my $vcf_row;
    my $prev_row;
    my $vcf_contig;
    my $bins_fh = FileHandle->new("< $args{'input_bins'}");
    my $bins_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($bins_fh, "<:encoding(utf8)", "$args{'input_bins'}") or die "Can't open bins_fh";
    $bins_tsv->column_names('contig', 'end');
    my $out = FileHandle->new("> $args{'output'}");

    if(scalar @db == 0) {
        copy("$args{'input_bins'}", "$args{'output'}");
    } else {
      BINS: while(my $bin_row = $bins_tsv->getline_hr($bins_fh)) {
          $contig = $bin_row->{'contig'};
          if (! grep /$contig/, @contigs) {
              $shiftsum = 0;
              push @contigs, $contig;
              @shift = 0;
          }
          my $end = $bin_row->{'end'};

        SHIFT: while ($vcf_row = $db["$vcf_index"]) {
            $vcf_contig = $vcf_row->seq_id;
            last SHIFT if(($end < $vcf_row->start) | ($contig ne $vcf_contig));
            @shift = $vcf_row->get_tag_values('shift');
            $shiftsum = $shiftsum + $shift[0];
            $vcf_index++;
        }
          #Some logic to sheck if deletion covers the feature
          $prev_row = $db["$vcf_index" - 1];
          @shift = $prev_row->get_tag_values('shift');
          if (($prev_row->start - $shift[0] > $end) & ($prev_row->seq_id eq $contig) & ($vcf_index != 0)) {
              print $out "$contig" . "\t" . $prev_row->start . "\n";
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
    my $chrom = '';
    my $relpos = 1; #track position on reference chromosome/contig
    #The below depthidx is how I track depth. Note it starts at -1
    #because when I print something I'm always "on the next
    #loop". Yes, I've tested this
    my $depthidx = -1;
    my $string;
    my $pos;
    my $length;
    my $depth;

  VARS: while(my $row = $raw_tsv->getline_hr($raw_fh)) {
      #Check if the chromosome has changed
      if($chrom ne $row->{'chromosome'} | ! defined($chrom)){
          if ($chrom ne '') {
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
                  $string = "$chrom\t" . "$pos\t" . ".\t" . "$ref\t" . "$alt\t" . ".\t" . ".\t" . "DP=$depths[$depthidx]\n";
                  print $vcf $string;
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
              $string = "$chrom\t" . "$pos\t" . ".\t" . "$ref\t" . "$alt\t" . ".\t" . ".\t" . "DP=$depths[$depthidx]\n";
              print $vcf $string;
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

    print $headfile "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    #Put the things from my vcf file into my final product
    $vcf = FileHandle->new("< $tmpdir/tmp_variants.vcf");
    while (<$vcf>) {
        print $headfile $_;
    };
    close($vcf);
    unlink("$tmpdir/tmp_variants.vcf");
}

#Now I want to take that vcf file and label it by whatever features I
#want. I also have the option of adding utrs if people want

sub vcf_labeler {
    my %args = @_;
    my $input_gff = $args{'gff'};
    my $input_vcf = $args{'vcf'};
    my $output = $args{'output'};
    my $add_flanks = $args{'add_flanks'};
    my $flank_lengths = $args{'flank_lengths'};
    my $feature_type = $args{'feature_type'};
    my $feature_name = $args{'feature_name'};

    #Parse the flank lengths properly as "up,down"
    my @flank_lengths = split(/\,/, $flank_lengths);

    my $gff = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $gff_fh = FileHandle->new("$input_gff");
    my $gff_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($gff_fh, "<:encoding(utf8)", "$input_gff");

    $gff_tsv->column_names('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes');

    my $feat;

    while(my $row = $gff_tsv->getline_hr($gff_fh)) {
        next if($row->{'seqid'} =~ /^#/);
        next if($row->{'type'} ne "$feature_type");

        #Trying to extract gene ID based on the feature name I was given
        my @tmp = split /;/, $row->{'attributes'};
        my %attributes;
        for my $item (@tmp) {
            my ($key, $value) = split('=', $item);
            $attributes{$key} = $value;
        }
        my $gene_id = $attributes{$feature_name} or die "Feature $feature_name not found in GFF";

        $feat = Bio::SeqFeature::Generic->new(
            -seq_id => $row->{'seqid'},
            -start => $row->{'start'},
            -end => $row->{'end'},
            -strand => $row->{'strand'},
            -primary_tag => $row->{'type'},
            -display_name => $gene_id);

        $gff->store($feat);

        #Put in the UTRs as well, which are strand dependent
        if($add_flanks eq 'yes') {
            if($row->{'strand'} eq '+') {
                $feat = Bio::SeqFeature::Generic->new(
                    -seq_id => $row->{'seqid'},
                    -start => ($row->{'start'} - $flank_lengths[0]),
                    -end => ($row->{'start'} - 1),
                    -strand => $row->{'strand'},
                    -primary_tag => 'up_flank',
                    -display_name => $gene_id);
                $gff->store($feat);
                $feat = Bio::SeqFeature::Generic->new(
                    -seq_id => $row->{'seqid'},
                    -start => ($row->{'end'} + 1),
                    -end => ($row->{'end'} + $flank_lengths[1]),
                    -strand => $row->{'strand'},
                    -primary_tag => 'down_flank',
                    -display_name => $gene_id);
                $gff->store($feat);
            } else {
                $feat = Bio::SeqFeature::Generic->new(
                    -seq_id => $row->{'seqid'},
                    -start => ($row->{'start'} - $flank_lengths[1]),
                    -end => ($row->{'start'} - 1),
                    -strand => $row->{'strand'},
                    -primary_tag => 'down_flank',
                    -display_name => $gene_id);
                $gff->store($feat);
                $feat = Bio::SeqFeature::Generic->new(
                    -seq_id => $row->{'seqid'},
                    -start => ($row->{'end'} + 1),
                    -end => ($row->{'end'} + $flank_lengths[0]),
                    -strand => $row->{'strand'},
                    -primary_tag => 'up_flank',
                    -display_name => $gene_id);
                $gff->store($feat);
            }
        }
    };

    my @gff_db = $gff->features();
    @gff_db = sort {
        $a->seq_id cmp $b->seq_id || $a->start <=> $b->start
    } @gff_db;

    #I've manually checked the above and it works as intended

    #Now that the gff is loaded, I need to load in my vcf file by pretty
    #much the same method

    my $vcf_fh = FileHandle->new("$input_vcf");
    my $vcf = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    $vcf_tsv->column_names('seq_id', 'pos', 'index', 'ref', 'alt', 'qual', 'filter', 'info');
    open($vcf_fh, "<:encoding(utf8)", "$input_vcf") or die "Can't open vcf_fh";

    while(my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if $row->{'seq_id'} =~ /^#/;
        my @tmp = split /;/, $row->{'info'};
        my %info;
        for my $item (@tmp) {
            my ($key, $value) = split('=', $item);
            $info{$key} = $value;
        };
        my $feat = Bio::SeqFeature::Generic->new(
            -seq_id => $row->{'seq_id'},
            -start => $row->{'pos'},
            -end => $row->{'pos'},
            -strand => 1,
            -primary_tag => 'variant',
            -tag => {ref =>$row->{'ref'},
                     alt => $row->{'alt'},
                     depth => $info{'DP'}
            }
            );
        my $ref = $row->{'ref'};
        my $alt = $row->{'alt'};
        if(length($ref) == 0) { $ref = 'N' };
        if(length($alt) == 0) { $alt = 'N' };

        $feat->add_tag_value('ref', $row->{'ref'});
        $feat->add_tag_value('alt', $row->{'alt'});
        $vcf->store($feat);
    };

    my @vcf_db = $vcf->features();
    @vcf_db = sort {
        $a->seq_id cmp $b->seq_id || $a->start <=> $b->start
    } @vcf_db;


    if(scalar @vcf_db == 0) { die "Nothing in the VCF file" };

    my $gff_index = 0;
    my $vcf_index = 0;
    my $gff_row;
    my @chroms;
    my $max = scalar(@gff_db);

    my $out = FileHandle->new("> $output");
    print $out "#" . "CHROM\t" . "POS\t" . "ID\t" . "REF\t" . "ALT\t" . "TYPE\t" . "FEATURE\t" . "NAME\t" . "STRAND\t" . "REL_START\t" . "REL_END\t" . "FRAG_DEPTH\n";

  VARS: while(my $vcf_row = $vcf_db["$vcf_index"]) {
      my $var_pos = $vcf_row->start;
      my @ref = $vcf_row->get_tag_values('ref');
      my $ref = $ref[0];
      my @alt = $vcf_row->get_tag_values('alt');
      my $alt = $alt[0];
      my @depth = $vcf_row->get_tag_values('depth');
      my $depth = $depth[0];
      my $var_chrom = $vcf_row->seq_id;

      #Determine the type of variant
      my $type;
      if(length($ref) == length($alt)) { $type = 'SNP' }
      elsif(length($ref) > length($alt)) { $type = 'del' }
      else { $type = 'ins' };

      #The below should loop until I have the proper feature
      my $feat_chrom;
      my $featured = 'no'; #keep track of if it's in a feature
    FEATS: while($gff_row = $gff_db["$gff_index"]) {
        $featured = 'yes';
        if($gff_index >= $max) {
            $featured = 'no';
            last FEATS;
        }
        $feat_chrom = $gff_row->seq_id;
        if(! grep /$feat_chrom/, @chroms) {
            push @chroms, $feat_chrom;
        } #keep track of chromosomes my gff has seen

        if($feat_chrom ne $var_chrom) {
            #If my variant is on a previous chromosome
            if(grep /$var_chrom/, @chroms) {
                $featured = 'no';
                last FEATS;
                #If my variant is on a later chromosome:
            } else {
                $gff_index++;
                next FEATS;
            }
        };
        #The below will stop the loop once a feature contains the variant
        if($gff_row->end < $var_pos) {
            $featured = 'no';
            $gff_index++;
            next FEATS;
        } elsif($gff_row->start > $var_pos) {
            $featured = 'no';
            last FEATS;
        } else {
            #This means the current gff_row contains the variant
            last FEATS;
        }
    }; #end FEATS
      #I should now have whichever feature I want stored in gff_row, so
      #now I just need some logic to determine the things I care about.

      #First, print a bunch of NA if it isn't in a feature
      if( "$featured" eq 'no') {
          print $out "$var_chrom\t" . "$var_pos\t" . ".\t" . "$ref\t" . "$alt\t" . "$type\t" . "NA\t" . "NA\t" . "NA\t" . "NA\t" . "NA\t" . "$depth\n";
          $vcf_index++;
          next VARS;
      }
      #Next, check the strandedness to determine metrics.

      my $rel_start;
      my $rel_end;
      my $strand = $gff_row->strand;;
      if ($strand == 1) {
          $strand = '+';
      } else {
          $strand = '-';
      }

      if($strand eq '+') {
          $rel_start = $var_pos - $gff_row->start + 1;
          $rel_end = $var_pos - $gff_row->end - 1;
      } else {
          #The - strand stuff has additional logic: I need to take the
          #complement of both strands, and I need to shift my
          #variant position to the "end" of any deletions
          $var_pos = $var_pos + length($ref) - 1; #Shouldn't change individual points
          $ref =~ tr/ACGT/TGCA/;
          $alt =~ tr/ACGT/TGCA/;
          #Ok now that is all sorted out, time to figure out the relative position
          $rel_start = $gff_row->end - $var_pos + 1;
          $rel_end = $gff_row->start - $var_pos - 1;
      }

      #Now I should have all the data I need to print my output
      print $out "$feat_chrom\t" . "$var_pos\t" . ".\t" . "$ref\t" . "$alt\t" . "$type\t" . $gff_row->primary_tag . "\t" . $gff_row->display_name . "\t" . $strand . "\t" . "$rel_start\t" . "$rel_end\t" . "$depth\n";
      $vcf_index++;
  }
}


#This seeks to shift feature annotations according to the indels
#within a vcf file. To this end, it is passed an input gff, an output
#file, and an input vcf. Right now it is compatible with gff only, but
#it shouldn't be too difficult to modify it to accept other file
#types.

sub gff_shifter {
    my %args = @_;
    my $starts = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $gff_fh = FileHandle->new("$args{'input_gff'}");
    my $gff_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($gff_fh, "<:encoding(utf8)", "$args{'input_gff'}");
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

    my $vcf_fh = FileHandle->new("$args{'input_vcf'}");
    my $vcf = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $vcf_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($vcf_fh, "<:encoding(utf8)", "$args{'input_vcf'}") or die "Can't open VCF file";
    $vcf_tsv->column_names('contig', 'position', 'index', 'reference', 'alternate', 'quality', 'filter', 'metadata');

    #Load in my VCF features and add them to a DB:
    while (my $row = $vcf_tsv->getline_hr($vcf_fh)) {
        next if $row->{'contig'} =~ /^#/;
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

    if (scalar @vcf_db == 0) {die "No indels in the VCF file!" };
    my $gff_out = FileHandle->new("> $args{'output'}");

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

1;
