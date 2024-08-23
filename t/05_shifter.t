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
    rmtree('test_output', verbose => 1, error => \my $err_list);
    exit;
};

Bio::Dantools::gff_shifter(
    gff => "${start_dir}/base.gff",
    vcf => "${start_dir}/exp_fasta_vars.vcf",
    output => "test.gff"
    );

my $exp_gff = File::Slurp::read_file("${start_dir}/exp_shift.gff");
my $act_gff = File::Slurp::read_file("test.gff");

unless (ok($exp_gff eq $act_gff, 'Expected GFF file produced')) {
    my ($expected, $actual) = diff($exp_gff, $act_gff);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}

chdir($start);
rmtree("test_output");
