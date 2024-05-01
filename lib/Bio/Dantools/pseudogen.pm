#!/usr/bin/env perl

package Bio::dantools::pseudogen;
use strict;
use autodie;
use warnings;
use FindBin;
use FileHandle;
use lib "$FindBin::Bin/../lib";
use Bio::dantools::fragment;

#use lib "$FindBin::Bin/../lib";

sub pseudogen {
    #This will serve as the main full pseudogen worker
    my %args = @_;
    my $bashfile = "$FindBin::Bin/../lib/aligner.sh";

    #Begin by fragmenting my reads
    Bio::dantools::fragment(input => "$args{'source'}", output => "$args{'outdir'}/fragments.fasta",
                            lengths => "$args{'lengths'}", overlap => "$args{'overlap'}",
                            min_length => "$args{'min_length'}",
                            logfile => "$args{'outdir'}/fragment_log.txt",
                            output_name => "$args{'output_name'}"
        );

    my $it = 0;
    my $itp = -1;

    my $var_counts = $args{'min_variants'} + 1;
    my $fai;
    my $base_idx;

    if ("$args{'fai'}" eq '') { $fai = 'NA' }
    if ("$args{'base_idx'}" eq '') { $base_idx = 'NA' }

    print "$args{'min_variants'}";

    until (("$var_counts" < "$args{'min_variants'}") & ("$it" != 0)) {
        system("sh", "$bashfile", '--outdir', "$args{'outdir'}",
               '-b', "$args{'base'}",
               '--fai', "$fai",
               '--indexes', "$base_idx",
               '-i', "$it",
               '--output_name', "$args{'output_name'}",
               '--fragments', "$args{'outdir'}/fragments.fasta",
               '-t', "$args{'threads'}"
            );
        #The above should have just produced everything I need to
        #continue through with my pipeline. That includes the
        #it"$it"/variants.vcf file, the MAPPING.stderr file, and the
        #new genome. The next step is to shift my bins accordingly


    }


};

1;
