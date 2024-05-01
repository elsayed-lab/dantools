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
    AUTHOR => 'Daniel Klimes <bio.dantools@gmail.com>',
    VERSION_FROM => 'lib/Bio/dantools.pm',
    ABSTRACT => 'dantools: Utilities to compare disparate genomes',
    PREREQ_PM => {
        'autodie' => '0',
        'Bio::DB::SeqFeature' => '0',
            'Getopt::Long' => '0',
            'FileHandle' => '0',
            'Text::CSV_XS::TSV' => '0',
            'Bio::SeqIO' => '0',
            'Bio::AlignIO' => '0',
            'Bio::SeqFeature::Generic' => '0',
            'Cwd' => '0',
            'POSIX' => '0',
            'FindBin' => '0',
            'File::Copy' => '0',
            'File::Basename' => '0',
            'Parallel::ForkManager' => '0'
    }
    )
