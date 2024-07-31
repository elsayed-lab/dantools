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

my $outdir = getcwd();
mkdir('tmp_danfiles');
#The purpose of this test is to examine the output of my dantools
#label command and make sure it is accurate.

Bio::Dantools::label(
    vcf => "${start_dir}/exp_vars.vcf",
    gff => "${start_dir}/base.gff",
    fasta => "${start_dir}/base.fasta",
    features => 'CDS',
    child_name => 'ID',
    parent_name => 'Parent',
    flank_lengths => '50,100',
    add_flanks => 1,
    flank_feature => 'CDS',
    flank_parent => 'Parent',
    threads => 1,
    output_nuc => 'test-nuc.tsv',
    translate => 1,
    coding_feature => 'CDS',
    codon_table => 1,
    score_matrix => 'blosum62',
    output_aa => 'test-aa.tsv',
    outdir => "$outdir",
    tmpdir => "tmp_danfiles"
    );

my $exp_nuc = File::Slurp::read_file("${start_dir}/exp_nuc.tsv");
my $act_nuc = File::Slurp::read_file("test-nuc.tsv");

unless (ok($exp_nuc eq $act_nuc, 'Expected nucleotide labels produced')) {
    my ($expected,$actual) = diff($exp_nuc, $act_nuc);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}

my $exp_aa = File::Slurp::read_file("${start_dir}/exp_aa.tsv");
my $act_aa = File::Slurp::read_file("test-aa.tsv");

unless (ok($exp_aa eq $act_aa, 'Expected amino acid labels produced')) {
    my ($expected,$actual) = diff($exp_aa, $act_aa);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}

chdir($start);
rmtree("test_output");
