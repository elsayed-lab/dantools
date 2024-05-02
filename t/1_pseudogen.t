# -*-Perl-*-

use strict;
use Test::More qw"no_plan";
use FindBin;
use Cwd;
use File::ShareDir qw"dist_dir";
use File::Which;
use File::Slurp;
use File::Path qw( rmtree );
use String::Diff;
use lib "$FindBin::Bin/../lib";
use Bio::Dantools;
use Test::File::ShareDir::Dist { 'Bio-Dantools' => 'share/' };
my $start_dir = dist_dir('Bio-Dantools');

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

my @needed = ('hisat2', 'hisat2-build', 'samtools', 'bcftools', 'freebayes-parallel', 'stretcher');
for my $bin (@needed) {
    die "Need binary ${bin}, not in PATH\n" unless(which("$bin"));
};

Bio::Dantools::pseudogen(gff => "$start_dir/base.gff",
                         base => "$start_dir/base.fasta",
                         fai => "$start_dir/base.fasta.fai",
                         base_idx => "$start_dir/indexes/base",
                         source => "$start_dir/source.fasta",
                         lengths => '100,200',
                         min_length => 20,
                         overlap => 50,
                         min_variants => 1000,
                         keepers => 'all',
                         output_name => 'test',
                         threads => 1,
                         bin_size => 10000,
                         add_flanks => 'yes',
                         flank_lengths => '50,50',
                         feature_types => 'protein_coding_gene,mRNA',
                         feature_name => 'ID',
                         var_fraction => 0.501,
                         outdir => "$start/test_output",
                         input_type => "fasta",
                         fragment => 'yes',
                         scoremin => 'L,0,-1.50'
    );

#Now I should have all the files ready and I can begin checking them
#for consistency. I think the most imporant test is using bcftools
#consensus with the vcf I generate and my original base to make sure I
#perfectly recreate the fragmented genome

system("bcftools view -O bcf -o output/test.bcf output/test.vcf");

system("bcftools index output/test.bcf");

my $error = qx"bcftools consensus -f $start_dir/base.fasta -o output/remade.fasta output/test.bcf 2>&1";

if (index($error, "Applied") == -1) {
    print "$error\n";
}

my $out = File::Slurp::read_file("output/test.fasta");
my $remade = File::Slurp::read_file("output/remade.fasta");

unless (ok($out eq $remade, 'The vcf file recreates the output')) {
    my ($old, $new) = diff($out, $remade);
    diag("--Output--\n${old}\n--Actual--\n${new}\n");
};

chdir($start);
rmtree("test_output");
