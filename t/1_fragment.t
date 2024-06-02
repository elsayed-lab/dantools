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

my @needed = ('hisat2', 'hisat2-build', 'samtools', 'bcftools', 'freebayes-parallel', 'stretcher');
for my $bin (@needed) {
    die "Need binary ${bin}, not in PATH\n" unless(which("$bin"));
};

Bio::Dantools::fragment(input => "${start_dir}/source.fasta",
                        output => "fragments.fasta",
                        lengths => '100,200',
                        overlap => "50",
                        min_length => "20",
                        logfile => 'fragment_log.txt',
                        threads => 1
                    );


#Now I should have all the files ready and I can begin checking them
#for consistency. I think the most imporant test is using bcftools
#consensus with the vcf I generate and my original base to make sure I
#perfectly recreate the fragmented genome

my $exp_frag = File::Slurp::read_file("${start_dir}/exp_frag.fasta");
my $act_frag = File::Slurp::read_file("fragments.fasta");

unless (ok($exp_frag eq $act_frag, 'Expected fragments produced')) {
    my ($expected, $actual) = diff($exp_frag, $act_frag);
    diag("--Output--\n${expected}\n--Actual--\n${actual}\n");
}

chdir($start);
rmtree("test_output");
