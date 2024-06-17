use strict;
use warnings;
use ExtUtils::MakeMaker;
use File::ShareDir::Install;
use inc::Module::Install;

install_share 'share';
auto_install( force => 1,);

WriteMakefile(
    LICENSE => 'perl',
    NAME => 'Bio-Dantools',
    AUTHOR => 'Daniel Klimes <daniel.s.klimes@gmail.com>',
    VERSION_FROM => 'lib/Bio/Dantools.pm',
    ABSTRACT => 'dantools: Utilities to compare disparate genomes',
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
