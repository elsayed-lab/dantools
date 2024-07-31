# -*-Perl-*-

use strict;
use Test::More qw"no_plan";
use FindBin;
use Cwd;
use File::ShareDir qw"dist_dir";
use File::Which;
use File::Slurp;
use File::Path qw( rmtree );
use String::Diff qw(diff);
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
    rmtree('test_output', verbose => 1, error => \my $err_list);
    exit;
};

Bio::Dantools::transloc(
    input_sam => "${start_dir}/exp_aln.sam",
    output => "$start/test_output/test_transloc.tsv",
    ref_gap => 10,
    query_gap => 10,
    max_diff => 50,
    min_qbase => 100,
    min_depth => 1,
    min_length => 1,
    input_vcf => "${start_dir}/exp_vars.vcf",
    all => 0
);

my $exp_tsv = File::Slurp::read_file("${start_dir}/exp_transloc.tsv");
my $act_tsv = File::Slurp::read_file("test_transloc.tsv");

unless (ok($exp_tsv eq $act_tsv, 'Expected TSV file produced')) {
    my ($expected, $actual) = diff($exp_tsv, $act_tsv);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}

chdir($start);
rmtree("test_output");
