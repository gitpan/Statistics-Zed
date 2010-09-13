use strict;
use warnings;
use Test::More tests => 11;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::Zed') };

my $zed = Statistics::Zed->new(ccorr => 1, tails => 2,);
isa_ok($zed, 'Statistics::Zed');

my %ref = (
	z_value => 1.625,
	p_value => 0.10416,
	p2z => 1.5000,
);

my %res = ();

eval { ($res{'z_value'}, $res{'p_value'}) = $zed->score(
	observed => 12,
      expected => 5,
      error => 4,
);};
ok(!$@);

foreach (qw/z_value p_value/) {
    ok(defined $res{$_} );
    ok(equal($res{$_}, $ref{$_}), "$_  $res{$_} = $ref{$_}");
}

eval { $res{'p2z'} = $zed->p2z(value => .066807, tails => 1);};
ok(!$@);
ok(equal($res{'p2z'}, $ref{'p2z'}), "p2z  $res{'p2z'} = $ref{'p2z'}");

eval { $res{'p2z'} = $zed->p2z(value => .133610, tails => 2);};
ok(!$@);
ok(equal($res{'p2z'}, $ref{'p2z'}), "p2z  $res{'p2z'} = $ref{'p2z'}");


sub equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}