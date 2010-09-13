use strict;
use warnings;
use Test::More tests => 7;
use constant EPS => 1e-2;

BEGIN { use_ok('Statistics::Zed') };

my $zed = Statistics::Zed->new(ccorr => 1, tails => 2,);
isa_ok($zed, 'Statistics::Zed');

my %refdat = (
	z_value => 1.625,
	p_value => 0.10416
);

my %res = ();

eval { ($res{'z_value'}, $res{'p_value'}) = $zed->score(
	observed => 12,
      expected => 5, # or key stdev or variance
      error => 4,
);};
ok(!$@);

foreach (qw/z_value p_value/) {
    ok(defined $res{$_} );
    ok(equal($res{$_}, $refdat{$_}), "$_  $res{$_} = $refdat{$_}");
}

sub equal {
    return 1 if $_[0] + EPS > $_[1] and $_[0] - EPS < $_[1];
    return 0;
}