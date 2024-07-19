use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;
use inc::Module::Install;
## Using ExtUtils::AutoInstall does the following: When running perl
## Makefile.PL for the first time in a fresh directory, it will figure
## out which (if any) modules are not installed.  Then, when make is
## invoked, those will be installed; when used in concert with the way
## we have modules set up on the cluster, this means that all
## dependencies will be installed inside the root of the module
## installation for dantools, which I think is the best of the various
## available installation worlds.
use ExtUtils::AutoInstall(
    -core => [
        'autodie' => '',
        'Bio::SeqIO' => '',
        'Bio::Seq' => '',
        'Bio::DB::Fasta' => '',
        'Bio::Tools::CodonTable' => '',
        'Bio::Matrix::IO' => '',
        'Cwd' => '',
        'FileHandle' => '',
        'FindBin' => '',
        'File::Basename' => '',
        'File::Copy' => '',
        'File::Path' => '',
        'File::Spec' => '',
        'File::Which' => '',
        'Getopt::Long' => '',
        'List::Util' => '',
        'List::MoreUtils' => '',
        'Parallel::ForkManager' => '',
        'POSIX' => '',
        'Text::CSV_XS::TSV' => '',
        ],);

install_share 'share';
auto_install(force => 1,);

## I also made a few changes in the body of the WriteMakefile, nothing
## substantive I think; with the caveat that I moved your aligner.sh
## and postprocess.sh to share/ and they should end up in bin/ upon
## invoking make install.
WriteMakefile(
    LICENSE => 'perl',
    NAME => 'Bio-Dantools',
    AUTHOR => 'Daniel Klimes <daniel.s.klimes@gmail.com>',
    VERSION_FROM => 'lib/Bio/Dantools.pm',
    ABSTRACT => 'dantools: Utilities to compare disparate genomes',
    EXE_FILES => ['script/dantools', 'script/Dantools', 'script/Dantools.pl', 'script/aligner.sh', 'script/postprocess.sh'],
    META_MERGE => {
        'meta-spec' => { version => 2 },
        resources => {
            repository => {
                type => 'git',
                url => 'https://github.com/elsayed-lab/dantools.git',
            },
            bugtracker => {
                web => 'https://github.com/elsayed-lab/dantools/issues',
            },
            homepage => 'https://github.com/elsayed-lab/dantools.git',
        },
    },
);
