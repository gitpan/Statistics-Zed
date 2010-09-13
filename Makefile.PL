use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Statistics::Zed',
    VERSION_FROM      => 'lib/Statistics/Zed.pm', # finds $VERSION
    PREREQ_PM         => { 
        'Test::More' => 0,
	'Math::Cephes' => 0.44,
	'Statistics::Distributions' => 1.02,
        'Statistics::Descriptive' => 2.6,
        'String::Util' => '0.11',
	},
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Statistics/Zed.pm', # retrieve abstract from module
       AUTHOR         => 'rgarton <rgarton AT cpan DOT org>') : ()),
);
