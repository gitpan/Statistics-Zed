#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Statistics::Zed' );
}

diag( "Testing Statistics::Zed $Statistics::Zed::VERSION, Perl $], $^X" );
