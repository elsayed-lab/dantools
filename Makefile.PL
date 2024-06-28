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
<<<<<<< HEAD
    EXE_FILES => ['script/dantools', 'script/Dantools', 'script/Dantools.pl', 'script/aligner.sh', 'script/postprocess.sh'],
    META_MERGE => {
        'meta-spec' => { version => 2 },
        resources => {
            repository => {
                type => 'git',
                url => 'https://github.com/dklimes/dantools.git',
            },
            bugtracker => {
                web => 'https://github.com/dklimes/dantools/issues',
            },
            homepage => 'https://github.com/dklimes/dantools.git',
        },
    },
);
=======
    PREREQ_PM => {
        'autodie' => '0',
        'Bio::DB::SeqFeature' => '0',
        'Bio::SeqIO' => '0',
        'Bio::Seq' => '0',
        'Bio::SeqFeature::Generic' => '0',
        'Bio::DB::Fasta' => '0',
        'Bio::DB::SeqFeature' => '0',
        'Bio::Tools::CodonTable' => '0',
        'Bio::Matrix::IO' => '0',
        'Cwd' => '0',
        'FileHandle' => '0',
        'FindBin' => '0',
        'File::Basename' => '0',
        'File::Copy' => '0',
        'File::Path' => '0',
	'File::ShareDir' => '0',
	'File::Slurp' => '0',
        'File::Spec' => '0',
        'File::Which' => '0',
        'Getopt::Long' => '0',
        'List::Util' => '0',
        'Parallel::ForkManager' => '0',
        'POSIX' => '0',
	'String::Diff' => '0',
	'Test::More' => '0',
        'Text::CSV_XS::TSV' => '0'

    }
    )
>>>>>>> upstream/master
