#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use File::Which;
use lib "$FindBin::Bin/../lib";
use Bio::Dantools;
use autodie;
use Cwd qw"abs_path cwd getcwd";
use Getopt::Long;

#Set up my interrupt trap
BEGIN {
    $SIG{'INT'} = sub {
        exit(1);
    };
}


#Figure out which method was passed
my $method = shift;

if (! defined($method)) {
    my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/dantools.help");
    while (<$helpdoc>) { print $_ };
    exit;
} elsif ("$method" eq 'help') {
    my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/dantools.help");
    while (<$helpdoc>) { print $_ };
    exit;
};


#Now I need to figure out given the method what I need to run

#This is the whole pseudogenome worker
if ($method eq 'pseudogen') {
    my $reference; #The genome I'll be aligning against
    my $reference_fai = ''; #.fai indexes for my reference, not required
    my $reference_idx = ''; #.ht2 reference indexes, not required
    my $query = ''; #The input genome
    my $readsu = ''; #input unpaired reads
    my $reads1 = ''; #mate 1 input reads
    my $reads2 = ''; #mate 2 input reads
    my $input_type;
    my $no_fragment = 0;
    my $lengths = '200,10000';
    my $min_length = 20;
    my $overlap = 75;
    my $min_variants = 100; #The threshold of variant # to repeat an iteration
    my $variant_caller = 'freebayes';
    my $continue = 0;
    my $output_name = '';
    my $keepers = 'all';
    my $nocap = 0;
    my $outdir = '';
    my $threads = 1;
    my $bin_size = 10000;
    my $var_fraction = 0.501;
    my $var_depth = 2;
    my $help = 0;
    my $scoremin = 'L,0,-1.0';
    GetOptions(
        "reference|r=s" => \$reference,
        "fai=s" => \$reference_fai,
        "ref-idx=s" => \$reference_idx,
        "query|q=s" => \$query,
        "reads-u|u=s" => \$readsu,
        "reads-1|1=s" => \$reads1,
        "reads-2|2=s" => \$reads2,
        "lengths|l=s" => \$lengths,
        "min-length=s" => \$min_length,
        "overlap=f" => \$overlap,
        "min-var=i" => \$min_variants,
        "no-cap" => \$nocap,
        "output-name=s" => \$output_name,
        "keep|k=s" => \$keepers,
        "outdir=s" => \$outdir,
        "threads|t=i" => \$threads,
        "bin-size=i" => \$bin_size,
        "var-fraction=f" => \$var_fraction,
        "var-depth=i" => \$var_depth,
        "var-caller=s" => \$variant_caller,
        "continue" => \$continue,
        "score-min=s" => \$scoremin,
        "no-fragment" => \$no_fragment,
        "help|h" => \$help
        ) or die "Error in Parsing Arguments";
    #Checking input arguments and setting them manually for some
    if (! defined($reference)) { $help = 1 };

    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/pseudogen.help");
        while (<$helpdoc>) { print $_ };
        exit;
    };

    #The first thing I need to do is check to make sure all of my
    #required binaries like hisat2 are loaded
    my @needed = ('hisat2', 'hisat2-build', 'samtools', 'bcftools', 'stretcher');
    for my $bin (@needed) {
        die "Need binary ${bin}, not in PATH\n" unless(which("$bin"));
    };

    if ($outdir eq '') {
        $outdir = getcwd;
    } elsif ( ! -d $outdir ) {
        make_path($outdir);
    }

    if (! grep /^$keepers/, ('all', 'variants', 'alignments', 'genomes', 'none')) {
        die "ERROR: The keepers argument \"$keepers\" is not supported, supported arguments:\nall, variants, alignments, genomes, none\n"
    }
    if (! -e $reference) {
        die "ERROR: input base file does not exist\n";
    }
    $reference = File::Spec->rel2abs($reference);

    $lengths =~ tr/ //d;

    #Figure out what "source_type" was brought in, which should make
    #downstream tools not always have to check this:
    if ($query ne '') {
        die "Error reading input query file $query" if (! -e $query);
        $input_type = 'fasta';
        $query = File::Spec->rel2abs($query);
        if ("$output_name" eq '') {
            $output_name = basename($query, ('.fasta', '.fastq')) . "_on_" . basename($reference, ('.fasta', '.fastq'));
        }
    } elsif ($readsu ne '') {
        die "Error reading input reads $readsu" if (! -e $readsu);
        $input_type = 'fastq_u'; #unpaired fastq
        $readsu = File::Spec->rel2abs($readsu);
        if ("$output_name" eq '') {
            $output_name = basename($readsu, ('.fasta', '.fastq')) . "_on_" . basename($reference, ('.fasta', '.fastq'));
        }
    } elsif (($reads1 ne '') & ("$reads2" ne '')) {
        die "Error reading input reads $reads1" if (! -e "$reads1");
        die "Error reading input reads $reads2" if (! -e "$reads2");
        $input_type = 'fastq_p'; #paired fastq
        $reads1 = File::Spec->rel2abs($reads1);
        $reads2 = File::Spec->rel2abs($reads2);
        if ($output_name eq '') {
            $output_name = basename($reads1, ('.fasta', '.fastq')) . "_on_" . basename($reference, ('.fasta', '.fastq'));
        }
    }

    #Check to make sure the variant caller I'm passed is reasonable
    my @possible_callers = ( 'freebayes', 'bcftools' );
    if (! grep { $_ eq $variant_caller } @possible_callers) {
        die "ERROR: $variant_caller not in list of possible callers\n"
    }

    if ($variant_caller eq 'freebayes') {
        die "ERROR: Need binary freebayes with --var-caller freebayes\n" unless(which('freebayes-parallel'));
    }


    if (($nocap == 1) & ("$input_type" ne 'fasta')) {
        print STDERR "WARNING: Argument 'no-cap' only applies to FASTA inputs\n";
    };

    Bio::Dantools::pseudogen(reference => "$reference",
                             fai => "$reference_fai",
                             reference_idx => "$reference_idx",
                             query => "$query",
                             readsu => "$readsu",
                             reads1 => "$reads1",
                             reads2 => "$reads2",
                             input_type => "$input_type",
                             no_fragment => "$no_fragment",
                             lengths => "$lengths",
                             min_length => "$min_length",
                             overlap => "$overlap",
                             min_variants => "$min_variants",
                             variant_caller => "$variant_caller",
                             continue => "$continue",
                             nocap => "$nocap",
                             keepers => "$keepers",
                             output_name => "$output_name",
                             outdir => "$outdir",
                             threads => "$threads",
                             bin_size => "$bin_size",
                             var_fraction => "$var_fraction",
                             var_depth => "$var_depth",
                             scoremin => "$scoremin"
        );
} elsif ("$method" eq 'fragment') {
    my $lengths = '200,1000';
    my $overlap = 75;
    my $min_length = 20;
    my $output = 'NO_OUTPUT_PROVIDED';
    my $log = 'fragment_log.txt';
    my $threads = 1;
    my $help = 0;

    GetOptions(
        "lengths|l=s" => \$lengths,
        "overlap|o=s" => \$overlap,
        "output=s" => \$output,
        "log=s" => \$log,
        "threads|t=i" => \$threads,
        "help|h" => \$help
        );

    my $input = $ARGV[0];
    $help = 1 if (! defined($input));

    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/fragment.help");
        while (<$helpdoc>) { print $_ };
        exit;
    };
    if (! -r "$input") {
        die "ERROR: couldn't read $input";
    }

    Bio::Dantools::fragment(input => "$input",
                            output => "$output",
                            lengths => "$lengths",
                            overlap => "$overlap",
                            min_length => "$min_length",
                            logfile => "$log",
                            threads => "$threads"
        );


} elsif ("$method" eq 'phylo') {
    die "Sorry, this method is still in development\n";

} elsif ("$method" eq 'label') {
    my $vcf;
    my $gff;
    my $fasta;
    my $features = '';
    my $child_name = 'ID';
    my $parent_name = 'Parent';
    my $threads = 1;
    my $output_nuc = 'nucleotide_vars.tsv';
    my $add_flanks = 0;
    my $all_vars = 0;
    my $flank_lengths = '200,200';
    my $flank_feature = '';
    my $flank_parent;
    my $add_info = 0; #whether to add the INFO field to the output
    #If they want it translated, they need the following:
    my $translate = 0;
    my $coding_feature = 'CDS';
    my $output_aa = 'amino_vars.tsv';
    my $codon_table = 1; #determine the codon table to use
    my $score_matrix = 'blosum62';
    my $outdir = '';
    my $tmpdir = '';
    my $help = 0;
    GetOptions(
        "vcf|v=s" => \$vcf,
        "gff|f=s" => \$gff,
        "fasta|b=s" => \$fasta,
        "features=s" => \$features,
        "child=s" => \$child_name,
        "parent=s" => \$parent_name,
        "add-flanks" => \$add_flanks,
        "flank-lengths=s" => \$flank_lengths,
        "flank-feature=s" => \$flank_feature,
        "flank-parent=s" => \$flank_parent,
        "add-info" => \$add_info,
        "threads|t=i" => \$threads,
        "outdir=s" => \$outdir,
        "output-nuc=s" => \$output_nuc,
        "output-aa=s" => \$output_aa,
        "translate" => \$translate,
        "coding-feature=s" => \$coding_feature,
        "codon-table=i" => \$codon_table,
        "score_matrix=s" => \$score_matrix,
        "tmpdir=s" => \$tmpdir,
        "all" => \$all_vars,
        "help|h" => \$help
        );
    if (! defined($vcf)) { print "ERROR: no VCF file provided\n"; $help = 1 };
    if (! defined($gff)) { print "ERROR: no GFF/GTF file provided\n"; $help = 1 };
    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/label.help");
        while (<$helpdoc>) { print $_ };
        exit;
    };
    $vcf = File::Spec->rel2abs($vcf);
    $gff = File::Spec->rel2abs($gff);

     if ($translate) {
        if (! defined($fasta)) {
            die "ERROR: translate specified but no FASTA provided\n"
        } else {
            $fasta = File::Spec->rel2abs($fasta);
        }

        if (! -e $fasta) {
            die "ERROR: FASTA file $fasta does not exist\n"        }
        elsif (! -r "$fasta") {
            die "ERROR: FASTA file $fasta is read protected\n"
        }
    };
    $fasta = File::Spec->rel2abs($fasta);
    $features =~ tr/ //d;
    if (! -r "$vcf") { die "ERROR: couldn't read $vcf" }
    elsif (! -r "$gff") { die "ERROR: couldn't read $gff" }

    if ($features eq '') { die "ERROR: no target features provided\n" }

    if ("$outdir" eq '') {
        $outdir = getcwd;
    } elsif ( ! -d "$outdir" ) {
        make_path("$outdir");
    }
    if ("$tmpdir" eq '') {
        $tmpdir = "${outdir}/tmp_danfiles";
    }
    if (! -d "$tmpdir") {
        make_path("$tmpdir")
    };

    #Try to make sure the coding feature is included in the features:
    if ((index($features, $coding_feature) == -1) & ($translate)) {
        die "ERROR: Coding feature not included in features\n";
    }

    #Parse the flank_lengths stuff
    if ($add_flanks) {
        $flank_lengths =~ tr/ //d;
        my @tmp = split(/\,/, "$flank_lengths");
        if (scalar(@tmp) == 1) {
            $flank_lengths = "$tmp[0],$tmp[0]";
        } elsif (scalar(@tmp) > 2) {
            die "Only two flank lengths are supported, e.g. 200,200";
        }
    };

    #If someone doesn't tell me what to add flanks to, I will try to
    #figure it out if they pass only a single feature to annotate
    if (($add_flanks) & ("$flank_feature" eq '')) {
        if (scalar(split(/\,/, "$features")) == 1) {
            $flank_feature = "$features";
        } else {
            die "ERROR: add-flanks specified but no flank feature provided\n";
        };
    };

    if (! defined($flank_parent)) {
        $flank_parent = $parent_name
    };

    Bio::Dantools::label(vcf => "$vcf",
                         gff => "$gff",
                         fasta => "$fasta",
                         features => "$features",
                         child_name => "$child_name",
                         parent_name => "$parent_name",
                         flank_lengths => "$flank_lengths",
                         add_flanks => "$add_flanks",
                         flank_feature => "$flank_feature",
                         flank_parent => "$flank_parent",
                         add_info => "$add_info",
                         threads => "$threads",
                         output_nuc => "$output_nuc",
                         translate => "$translate",
                         coding_feature => "$coding_feature",
                         codon_table => "$codon_table",
                         score_matrix => "$score_matrix",
                         output_aa => "$output_aa",
                         outdir => "$outdir",
                         all_vars => "$all_vars",
                         tmpdir => "$tmpdir"
        );

} elsif ("$method" eq 'transloc') {
    #This function seeks to create a TSV file indicating which regions
    #on one chromosome correspond to those regions on another
    #species'.

    my $input_sam;
    my $output;
    my $ref_gap = 1000;
    my $query_gap = 1000;
    my $min_qbase = 10000;
    my $max_diff = 50; #% difference in length between ref and query
    my $min_depth = 1;
    my $min_length = 1;
    my $input_vcf = 0;
    my $all = 0; #print all found regions passing thresholds?
    my $help = 0;
    GetOptions(
        "sam|s=s" => \$input_sam,
        "output|o=s" => \$output,
        "reference-gap|r=i" => \$ref_gap,
        "query-gap|q=i" => \$query_gap,
        "min-qbase|b=i" => \$min_qbase,
        "min-depth|d=i" => \$min_depth,
        "min-length|l=i" => \$min_length,
        "max-diff=i" => \$max_diff,
        "vcf|v=s" => \$input_vcf,
        "all|a" => \$all,
        "help|h" => \$help
        );

    if (! defined($input_sam)) {
        $input_sam = $ARGV[0];
    }
    if (! defined($input_sam)) {
        $help = 1;
    }
    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/transloc.help");
        while (<$helpdoc>) { print $_ };
        exit;
    };

    if (! -e $input_sam) {
        die "ERROR: Input SAM $input_sam does not exist\n";
    } elsif (! -r $input_sam) {
        die "ERROR: Input SAM $input_sam is read-protected\n";
    }

    if (! defined($output)) {
        $output = 'NO_OUTPUT_PROVIDED'
    }

    if ($input_vcf && ! -e $input_vcf) {
        die "ERROR: Input VCF $input_vcf does not exist\n";
    }

    Bio::Dantools::transloc(
        input_sam => "$input_sam",
        output => "$output",
        ref_gap => "$ref_gap",
        query_gap => "$query_gap",
        max_diff => "$max_diff",
        min_qbase => "$min_qbase",
        min_depth => "$min_depth",
        min_length => "$min_length",
        input_vcf => "$input_vcf",
        all => "$all"
    )

} elsif ("$method" eq 'make-vcf') {
    #Just a function for me when I mess up. Uses a combined bases file
    #to generate a vcf file
    my $input;
    my $output;
    my $tmpdir = 'tmp_danfiles';
    my $depths;
    GetOptions(
        "input|i=s" => \$input,
        "output|o=s" => \$output,
        "tmpdir=s" => \$tmpdir,
        "depths|d=s" => \$depths
    );

    Bio::Dantools::vcf_maker(
        input => "$input",
        output => "$output",
        tmpdir => "$tmpdir",
        depths => "$depths"
    );


} elsif ("$method" eq 'shift') {
    my $gff;
    my $vcf;
    my $output = 'NO_OUTPUT_PROVIDED';
    my $help = 0;
    GetOptions(
        "vcf|v=s" => \$vcf,
        "gff|f=s" => \$gff,
        "output|o=s" => \$output,
        "help|h" => \$help
        );
    if (! defined($vcf)) { $help = 1 };
    if (! defined($gff)) { $help = 1 };
    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/shift.help");
        while (<$helpdoc>) { print $_ };
        exit;
    };

    $vcf = File::Spec->rel2abs($vcf);
    $gff = File::Spec->rel2abs($gff);

    if (! -e "$vcf") { die "ERROR: couldn't read $vcf" }
    elsif (! -r "$gff") { die "ERROR: couldn't read $gff" }

    Bio::Dantools::gff_shifter(
        gff => "$gff",
        vcf => "$vcf",
        output => "$output"
    );

} elsif ("$method" eq 'summarize-aa') {
    my $input_vars;
    my $output = 'NO_OUTPUT_PROVIDED';
    my $outscore = 0; #score assigned to out of frame mutations
    my $help = 0;
    GetOptions(
        "vars|v=s" => \$input_vars,
        "output|o=s" => \$output,
        "outframe-score=s" => \$outscore,
        "help|h" => \$help
    );
    if (! defined($input_vars)) {
        $input_vars = $ARGV[0];
    };
    if (! defined($input_vars)) {
        $help = 1;
    }
    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/summarize-aa.help");
        while (<$helpdoc>) { print $_ };
        exit;
    };
    if (! -e "$input_vars") {
        die "Input variant file does not exist\n";
    } elsif (! -r "$input_vars") {
        die "Input variant file is read-protected\n";
    };

    Bio::Dantools::summarize_aa(
        input_vars => "$input_vars",
        output => "$output",
        outscore => "$outscore"
        );
} elsif ("$method" eq 'summarize-nuc') {
    my $input_vars;
    my $output = 'NO_OUTPUT_PROVIDED';
    my $help = 0;
    GetOptions(
        "vars|v=s" => \$input_vars,
        "output|o=s" => \$output,
        "help|h" => \$help
    );
    if (! defined($input_vars)) {
        $input_vars = $ARGV[0]
    }
    if (! defined($input_vars)) {
        $help = 1;
    }

    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/summarize-nuc.help");
        while (<$helpdoc>) { print$_ };
        exit;
    }

    if (! -e $input_vars) {
        die "Input variant file does not exit\n";
    }  elsif (! -r $input_vars) {
        die "Input variant file is read-protected\n"
    }

    Bio::Dantools::summarize_nuc(
        input_vars => "$input_vars",
        output => "$output"
        );


} elsif ("$method" eq 'summarize-depth') {
    my $input_depth;
    my $input_gff;
    my $feature = 'CDS';
    my $parent = 'Parent';
    my $output = 'NO_OUTPUT_PROVIDED';
    my $help = 0;
    GetOptions(
        "gff|f=s" => \$input_gff,
        "output|o=s" => \$output,
        "depth|d=s" => \$input_depth,
        "feature=s" => \$feature,
        "parent=s" => \$parent,
        "help" => \$help
        );
    if (! defined($input_depth)) { $help = 1 };
    if (! defined($input_gff)) { $help = 1 };

    if ($help) {
        my $helpdoc = FileHandle->new("< $FindBin::Bin/../helpdocs/summarize-depth.help");
        while (<$helpdoc>) { print $_ };
        exit;
    };


    if (! -e "$input_depth") {
        die "Input depth file does not exit\n";
    } elsif (! -r "$input_depth") {
        die "Input depth file is read protected\n";
    }

    if (! -e "$input_gff") {
        die "Input GFF does not exit\n";
    } elsif (! -r "$input_gff") {
        die "Input GFF is read protected\n";
    }

    Bio::Dantools::summarize_depth(
        input_gff => "$input_gff",
        input_depth => "$input_depth",
        output => "$output",
        parent => "$parent",
        feature => "$feature"
        );
} else {
    die "Method $method not found\n";
}
;
