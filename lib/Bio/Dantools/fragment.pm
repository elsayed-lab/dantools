#!/usr/bin/env perl

package Bio::dantools::fragment;
use warnings;
use strict;
use autodie;
use POSIX qw"floor ceil";
use Bio::SeqIO;
use FileHandle;

sub fragment {

    my %args = @_;

    my $seqio_in = Bio::SeqIO->new(-format => 'fasta', -file => $args{'input'});
    my $seqio_out = Bio::SeqIO->new(-format => 'fasta', -file => "> $args{'output'}");
    my @sizes = split(/\,/, $args{'lengths'});
    my $log = FileHandle->new("> $args{'logfile'}");
    my $num_written = 0;

  SEQS: while (my $seq = $seqio_in->next_seq) {
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
  }
    $log->close;
}

1;
