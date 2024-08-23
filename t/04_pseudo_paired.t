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
my $start_dir = dist_dir('Bio-Dantools') . "/test_files";

my $start = getcwd();
my $new = 'test_output';
mkdir($new);
chdir($new);

$SIG{'INT'} = sub {
    chdir($start);
    #rmtree('test_output', verbose => 1, error => \my $err_list);
    exit;
};

my @needed = ('hisat2', 'hisat2-build', 'samtools', 'bcftools', 'freebayes-parallel', 'stretcher');
for my $bin (@needed) {
    die "Need binary ${bin}, not in PATH\n" unless(which("$bin"));
};

Bio::Dantools::pseudogen(base => "$start_dir/base.fasta",
                         fai => "$start_dir/base.fasta.fai",
                         base_idx => "$start_dir/indexes/base",
                         source => "$start_dir/source.fasta",
                         min_variants => 100,
                         keepers => 'all',
                         output_name => 'test',
                         continue => 0,
                         threads => 1,
                         bin_size => 10000,
                         var_fraction => 0.501,
                         var_depth => 2,
                         variant_caller => 'bcftools',
                         outdir => "$start/test_output",
                         input_type => "fastq_p",
                         no_fragment => 0,
                         scoremin => 'L,0,-1.50',
                         reads1 => "$start_dir/paired_1.fastq",
                         reads2 => "$start_dir/paired_2.fastq"
    );

#Now I should have all the files ready and I can begin checking them
#for consistency. I think the most imporant test is using bcftools
#consensus with the vcf I generate and my original base to make sure I
#perfectly recreate the fragmented genome

system("bcftools view -O bcf -o output/test.bcf output/test.vcf");

system("bcftools index output/test.bcf");

#I can use bcftools consensus here because my variants should not
#overlap under any circumstance. This justa adds another layer of
#checking instead of using my own script
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

my $exp_vcf = File::Slurp::read_file("${start_dir}/exp_paired_vars.vcf");
my $act_vcf = File::Slurp::read_file("output/test.vcf");

unless (ok($exp_vcf eq $act_vcf, 'Expected VCF file produced')) {
    my ($expected, $actual) = diff($exp_vcf, $act_vcf);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}

my $exp_depth = File::Slurp::read_file("${start_dir}/exp_paired_depth.tsv");
my $act_depth = File::Slurp::read_file("output/depth.tsv");

unless (ok($exp_depth eq $act_depth, 'Expected depth file produced')) {
    my( $expected, $actual) = diff($exp_depth, $act_depth);
    diag("--Output\n${expected}\n--Actual--\n${actual}\n");
};

chdir($start);
#rmtree("test_output");
