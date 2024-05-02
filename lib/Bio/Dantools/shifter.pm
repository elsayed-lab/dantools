#!/usr/bin/env perl

package Bio::dantools;
use strict;
use autodie;
use warnings;
use FindBin;
use FileHandle;
use lib "$FindBin::Bin/../lib";
use Bio::DB::SeqFeature;

#This file should include the scripts necessary to do any and all
#shifting I need. That includes shifting the bins and the features. I
#also need a subroutine that creates the bins for iteration 0

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
    my $vcf_fh = FileHandle->new("> $args{'input_vcf'}");
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
        print "Nothing in the VCF file, copying bins\n";
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
          if (($vcf_row->start - $shift[0] > $end) & ($vcf_row->seq_id eq $contig) & ($vcf_index != 0)) {
              print $out "$contig" . "\t" . $prev_row->start . "\n";
          } else {
              my $new_out = $end + $shiftsum;
              print $out "$contig" . "\t" . "$new_out\n";
          }
      }
    }
}


#This is the most complicated of the bunch of functions defined
#here. It seeks to shift feature annotations according to the indels
#within a vcf file. To this end, it is passed an input gff, an output
#file, and an input vcf. Right now it is compatible with gff only, but
#it shouldn't be too difficult to modify it to accept other file
#types.

sub gff_shifter {
    my %args = @_;
    my $gff = Bio::DB::SeqFeature::Store->new(-adaptor => 'memory');
    my $gff_fh = FileHandle->new("$args{'input_gff'}");
    my $gff_tsv = Text::CSV_XS::TSV->new({binary => 1, });
    open($gff_fh, "<:encoding(utf8)", "$args{'input_gff'}");
    $gff_tsv->column_names('seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes');

    my $feat;
    my @heading;
  READ_GFF: while (my $row = $gff_tsv->getline_hr($gff_fh)) {

  }
}

1;
