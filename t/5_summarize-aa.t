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

Bio::Dantools::summarize_aa(
        input_vars => "${start_dir}/exp_aa.tsv",
        output => "test_summary.tsv",
        outscore => 0
        );

my $exp_aa = File::Slurp::read_file("${start_dir}/exp_summarize_aa.tsv");
my $act_aa = File::Slurp::read_file("test_summary.tsv");

unless (ok($exp_aa eq $act_aa, 'Expected amino acid summary produced')) {
    my ($expected,$actual) = diff($exp_aa, $act_aa);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}


chdir($start);
rmtree("test_output");
