use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Statistics::Zed',
    license             => 'perl',
    dist_author         => 'Roderick Garton <rgarton@cpan.org>',
    dist_version_from   => 'lib/Statistics/Zed.pm',
    build_requires => {
        'Test::More' => 0,
	'Math::Cephes' => 0.44,
	'Statistics::Distributions' => 1.02,
        'Statistics::Descriptive' => 2.6,
        'String::Util' => '0.11',
        'Statistics::Lite' => '3.2',
    },
    add_to_cleanup      => [ 'Statistics-Zed-0.05-*' ],
);

$builder->create_build_script();
