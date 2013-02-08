package Statistics::Zed;

use 5.008008;
use strict;
use warnings;
use Carp qw(croak);
use vars qw($VERSION);
$VERSION = 0.06;
use Math::Cephes qw(:dists);
use Statistics::Lite qw(:all);
use Statistics::Descriptive;
use Statistics::Distributions;
use String::Util qw(nocontent); 
use Scalar::Util qw(looks_like_number);

=head1 NAME

Statistics::Zed - Basic deviation ratio: observed less expected (with optional continuity correction) divided by root variance 

=head1 SYNOPSIS

 use Statistics::Zed 0.06;

 $zed = Statistics::Zed->new(
    ccorr    => 1, 
    tails    => 2,
    precision_s => 5,
    precision_p => 5,
 );

 ($z_value, $p_value, $observed_deviation, $standard_deviation) = $zed->score( # or 'zscore'
    observed => $obs,
    expected => $exp,
    variance => $variance, # or 'stdev'
 );

 $p_value = $zed->z2p(value => $z_value, tails => 1|2);
 $z_value = $zed->p2z(value => $p_value, tails => 1|2);

=head1 DESCRIPTION

+ Calculates a standard, run-of-the-mill z-score: the ratio of an observed deviation to a standard deviation.

+ Provides wraps to convert z-value to p-value, and convert p-value to z-value.

+ Organizes accumulated observed, expected and variance values over two or more runs ahead of making the calculation.

+ Purpose is to support tests in  L<Statistics::Sequences|Statistics::Sequences>.

=head1 METHODS

=head2 new

 $zed = Statistics::Zed->new();

Returns a Statistics::Zed object. Accepts setting of any of the OPTIONS.

=cut

sub new {
    my $class = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
	my $self = {};
	bless $self, $class;

    # Set default values:
    $self->{'tails'} = 2;
    $self->{'ccorr'} = 0;
    
    if (scalar keys %{$args}) {
        foreach (keys %{$args}) {
            $self->{$_} = $args->{$_};
        }
    }
	return $self;
}

=head2 zscore

 $zval = $zed->zscore(observed => 'number', expected => 'number', variance => 'number (non-zero)')
 ($zval, $pval, $obs_dev, $stdev) = $zed->zscore(observed => 'number', expected => 'number', variance => 'number (non-zero)')

I<Alias>: score

You supply the C<observed> and C<expected> values of your statistic, and the C<variance> or C<stdev>.

Optionally specify a logical value for L<ccorr|ccorr> for performing the continuity-correction to the observed deviation, and a value of either 1 or 2 to specify the L<tails|tails> for reading off the probability associated with the I<z>-value. When passed to this function, values of L<ccorr|ccorr> and L<tails|tails> become the values used in this and subsequent tests.

When called in array context, returns an array consisting of the z-statistic ("standard score"), its probability, the observed deviation (the difference between the observed and expected values of your statistic), and the standard deviation (the square-root of the variance supplied).

If you only want the I<z>-value, then expect a string:

 $z_value = $zed->score(...)

The basic formula is the basic:

=for html <p>&nbsp;&nbsp;&nbsp;<i>Z</i> = ( <i>&times;</i> &ndash; <i><o>X</o></i> ) / SD</p>

where I<X> is the expected value (mean, etc.).

=cut

sub zscore {
    my $self = shift;
	my $args = ref $_[0] ? $_[0] : {@_};
    $self->{$_} = $args->{$_} foreach keys %{$args};
    croak "Need to define observed ($self->{'observed'}) and expected ($self->{'expected'}) values for zscore" if nocontent($self->{'observed'}) || nocontent($self->{'expected'});

    my ($z, $p, $obs_dev, $exp_dev) = ();
    $obs_dev = $self->{'observed'} - $self->{'expected'};
    $obs_dev = _ccorr($obs_dev) if $self->{'ccorr'};
    $exp_dev = _exp_dev($self);
    return 0 if !$exp_dev;

    $z = $obs_dev / $exp_dev;
    $p = $self->z2p(value => $z, tails => $self->{'tails'});
    $z = sprintf('%.' . $self->{'precision_s'} . 'f', $z) if $self->{'precision_s'};
    return wantarray ? ($z, $p, $obs_dev, $exp_dev) : $z;
}
*score = \&zscore;

=head2 z2p

 $p = $zed->z2p($z); # assumes 2-tailed
 $p = $zed->z2p(value => $z); # assumes 2-tailed
 $p = $zed->z2p(value => $z, tails => 1);

I<Alias>: C<p_value>

Send a I<z>-value, get its associated I<p>-value, 2-tailed by default, or depending on what the value of $zed->{'tails'} is, or what is sent as the second argument, if anything. If you send just one value (unkeyed), it is taken as the value.

=cut

sub z2p {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : scalar(@_) > 1 ? {@_} : {value => shift};
    $args->{'tails'} ||= 2;
    return '' if nocontent($args->{'value'}) or !looks_like_number($args->{'value'});
    return 1 if $args->{'value'} == 0;
	my $p = Statistics::Distributions::uprob(abs($args->{'value'}));
    $p *= 2 if $args->{'tails'} == 2;
	$p = Statistics::Distributions::precision_string($p);
    $p = 1 if $p > 1;
    $p = 0 if $p < 0;
    $p = sprintf('%.' . $self->{'precision_p'} . 'f', $p) if $self->{'precision_p'};
    return $p;
}
*p_value = \&z2p;

=head2 p2z

 $z_value = $zed->p2z($p) # the p-value is assumed to be 2-tailed
 $z_value = $zed->p2z(value => $p) # the p-value is assumed to be 2-tailed
 $z_value = $zed->p2z(value => $p, tails => 1) # specify 1-tailed probability

Returns the I<z>-value associated with a I<p>-value using L<Math::Cephes|Math::Cephes> C<ndtri> ("phi"). B<The p-value is assumed to be two-tailed>, and so is firstly (before conversion) divided by 2, e.g., .05 becomes .025 so you get I<z> = 1.96.  As a one-tailed probability, it is then assumed to be a probability of being I<greater> than a certain amount, i.e., of getting a I<z>-value I<greater> than that observed. So the phi function is actually given (1 - I<p>-value) to work on - so you get back the I<z>-value up to that point, before the I<p>-value has to fit in. So .055 comes back as 1.598 (speaking of the top-end of the distribution), and .991 comes back as -2.349 (now going from right to left across the distribution). These assumptions are not the same as found in inversion methods in common spreadsheet packages but seem expected by human users like me of I<z>-values and their I<p>-values. (complain if useful)

=cut

sub p2z {
    my $self = shift;
    my ($p_val, $z_val, $tails) = (); 
    if (scalar(@_) > 1) {
        my $args = ref $_[0] ? $_[0] : {@_};
        $p_val = $args->{'value'};
        $tails = $args->{'tails'} || 2;
    }
    else {
        $p_val = shift;
        $tails = 2;
    }
    #return undef if !_valid_p($p_val);
    #croak __PACKAGE__, "::p2z The value sent as a p-value (" . ($p_val || '') . ") does not appear to be valid" if !_valid_p($p_val);

    # Avoid ndtri errors by first accounting for 0 and 1 ...
    if (!$p_val or $p_val < 0) {
        $z_val = '';
    }
    elsif ($p_val >= 1)  {
        $z_val = 0;
    }
    else {
        $p_val /= 2 if $tails == 2; # p-value has been given as two-tailed - we only use one side
        $z_val = ndtri(1 - $p_val); # 
    }
    return $z_val;
}


# Apply the continuity correction to the deviation, e.g. of (x - MCE) in numerator of z-score, if variance is calculated binomially (as hit/miss, not continuous):
sub _ccorr {
   my $dev = shift;
   if ($dev) {
       my $c = abs($dev) - .5;
       $c *= -1 if $dev < 0;
       return $c;
   }
   else {
       return $dev; # thanks for nothing
   }
}

=head2 dump

Prints to STDOUT a line giving the z_value and p_value.

=cut

sub dump {
    my $self = shift;
    print "Z = $self->{'z_value'}, $self->{'tails'}p = $self->{'p_value'}\n";
}

# Calc expected deviation, ensuring valid denominator:
sub _exp_dev { # ugly elongated thing :
    my $self = shift;
    return ($self->{'stdev'} and $self->{'stdev'} > 0) ? 
      $self->{'stdev'} : 
        ($self->{'variance'} and $self->{'variance'} > 0) ?
            sqrt($self->{'variance'}) :
               0;#croak "Need a standard deviation or variance value for z-testing";##return undef;##
}

=head1 Series testing

Aggregate results from multiple tests using a method to initialise the series, another to add data to the series, and a final one to sum the data, compute an aggregate Zscore from them, with all the values of interest ready to dump. L<Statistics::Descriptive|Statistics::Descriptive> objects are kept of each series' observed, expected and variance values, as sent to L<series_update>, and a separate count of the number of series is kept - $zed->{'series'}->{'count'}.

=head2 series_init

Clears any already accumulated data from previous tests.

=head2 series_update

 $zed->series_update() # use any cached values of "observed", "expected" and "variance"
 $zed->series_update(variance => 'number', expected => 'number', observed => 'number') # supply own in a hash

Called once you have performed a test on a sample. It caches the observed, matchcount_expected and variance values from the test.

=head2 series_test

Sums the observed, matchcount_expected and variance values from all the tests updated to the series since calling L<series_init|series_init>, and produces a I<z>-value from these sums. Returns the same values as C<zscore>: wanting an array, then the I<z>-value, its probability, the observed deviation and the standard deviation; otherwise, the I<z>-value alone. Additionally, the C<$zed->{'series'}> object hash is lumped with the usual values, so you can access them like so:

 print "Series summed runs: 
    expected = ", $zed->{'series'}->{'expected'}, " 
    observed = ", $zed->{'series'}->{'observed'}," 
    z = $zed->{'series'}->{'z_value'}, $zed->{'tails'}p = $zed->{'series'}->{'p_value'}\n";

=cut

sub series_init {
    my ($self, %args) = @_;
    $self->{'series_stat'}->{$_} = Statistics::Descriptive::Sparse->new() foreach qw/observed expected variance samplings z_value/;
    $self->{'series'}->{$_} = 0 foreach qw/observed expected variance samplings z_value p_value count/;
   # $self->{'series'}->{'p_value'} = 1;
  # $self->{'series'}->{'expected'} = 1;
}

sub series_update {
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    if (! defined $args->{'variance'}) {
        if (defined $args->{'stdev'}) {
            $args->{'variance'} = $args->{'stdev'}**2;
        }
        else {
            croak "No variance in series_update for calculating deviation";
        }
    }
	foreach (qw/observed expected variance samplings z_value/) {
    	$self->{'series_stat'}->{$_}->add_data($args->{$_}) if defined $args->{$_};
	}
    $self->{'series'}->{'count'}++;
}

sub series_test {
    my ($self, @args) = (shift, @_);
    return  if ! $self->{'series_stat'}->{'variance'}->count();
   # croak 'No data for series-testing appear to have been loaded; maybe you need to call series_update with some data' if ! $self->{'series_stat'}->{'variance'}->count();
    my $o = $self->{'series_stat'}->{'observed'}->sum();
    my $e = $self->{'series_stat'}->{'expected'}->sum();
    my $v = $self->{'series_stat'}->{'variance'}->sum();
    my ($z, $pz, $obs_dev, $exp_dev) = $self->zscore(
        observed => $o,
        expected => $e,
        variance => $v,
        @args,
    );
    $self->{'series'}->{'observed'} = $o;
    $self->{'series'}->{'expected'} = $e;
    $self->{'series'}->{'z_value'} = $z;# || 0;
    $self->{'series'}->{'p_value'} = $pz;# || 1;
    $self->{'series'}->{'obs_dev'} = $obs_dev;
    $self->{'series'}->{'stdev'} = $exp_dev;
    $self->{'series'}->{'variance'} = $v;
    $self->{'series'}->{'samplings'} = $self->{'series_stat'}->{'samplings'}->sum();
	if ($self->{'series_stat'}->{'z_value'}->count()){
    	$self->{'series'}->{'stouffer_z'} = $self->{'series_stat'}->{'z_value'}->sum() / sqrt($self->{'series_stat'}->{'z_value'}->count());
    	$self->{'series'}->{'stouffer_p'} = $self->p_value($self->{'series'}->{'stouffer_z'});
	}
    return wantarray ? ($z, $pz, $obs_dev, $exp_dev) : $z;
}

=head2 series_str

Returns string: line giving the z_value and p_value for the series.

=cut

sub series_str {
    my ($self) = @_;
    return "Z (N = $self->{'series'}->{'samplings'}) = $self->{'series'}->{'z_value'}, $self->{'tails'}p = $self->{'series'}->{'p_value'}";
}

=head2 series_dump

Prints to STDOUT a line giving the z_value and p_value for the series.

=cut

sub series_dump {
    print series_str(@_), "\n";
}

sub _valid_p {
    my $p = shift;
    return undef if nocontent($p) or !looks_like_number($p);
    return ( ($p !~ /^0?\.\d+$/) && ($p !~ /^\d+?\.[\de-]+/ )) || ($p < 0 || $p > 1) ? $p == 1 ? 1 : 0 : 1;
}

1;
__END__

=head1 OPTIONS

The following can be set in the call to C<new> or C<test> or C<score>.

=head2 ccorr

Apply the continuity correction, Stetigkeitskorrektur. Default = 0.

=head2 tails

Tails from which to assess the association I<p>-value (1 or 2). Default = 2.

=head2 precision_s

Precision of the I<z>-value (the statistic). Default = 2.

=head2 precision_p

Precision of the associated I<p>-value. Default = 0 - you get all decimal values available.

=head1 SEE ALSO

L<Statistics::Distributions|Statistics::Distributions> : C<uprob> function, and C<precision_string> method is used here for reporting probability.

L<Statistics::Sequences|Statistics::Sequences> : for application of this module.

=head1 TO DO/BUGS

Other distributions.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2010 R Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut
