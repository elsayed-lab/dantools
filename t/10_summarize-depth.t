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

Bio::Dantools::summarize_depth(
    input_gff => "${start_dir}/exp_shift.gff",
    input_depth => "${start_dir}/exp_fasta_depth.tsv",
    output => "test_summary.tsv",
    parent => 'Parent',
    feature => 'CDS'
);

my $exp_depth = File::Slurp::read_file("${start_dir}/exp_summarize_depth.tsv");
my $act_depth = File::Slurp::read_file("test_summary.tsv");

unless (ok($exp_depth eq $act_depth, 'Expected depth summary produced')) {
    my ($expected,$actual) = diff($exp_depth, $act_depth);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}


chdir($start);
rmtree("test_output");
