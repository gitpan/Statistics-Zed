package Statistics::Zed;

use 5.008008;
use strict;
use warnings;
use Carp qw(croak);
use vars qw($VERSION);
$VERSION = 0.03;
use Statistics::Lite qw(:all);
use Statistics::Descriptive;
use Statistics::Distributions;
use Math::Cephes qw(:dists);
use String::Util qw(hascontent); 

=head1 NAME

Statistics::Zed - Basic ztest/zscore, with optional continuity correction

=head1 SYNOPSIS

 use Statistics::Zed 0.03;

 my $zed = Statistics::Zed->new(
    ccorr    => 1,
    tails    => 2,
    precision_s => 5,
    precision_p => 5,
 );

 my ($z_value, $p_value, $observed_deviation, $standard_deviation) = 
    $zed->score( # or score
        observed => $obs,
        expected => $exp, # or key stdev or variance
        error => $variance,
 );

 my $deviate = $zed->test( # or ztest
        observed => $obs,
        expected => $exp,
        error => $variance, # or key stdev or variance
        samplings => $samplings,
 );

 my $p_value = $zed->p_value($deviate); # or z2p

=head1 DESCRIPTION

Calculates a standard, run-of-the-mill z-score: the ratio of an observed deviation to a standard deviation; and performs a z-test or returns the z-score. Purpose is simply to support L<Statistics::Sequences|Statistics::Sequences>.

=head1 METHODS

=head2 Interface

=head3 new

 $zed = Statistics::Zed->new();

Returns a Statistics::Zed object. Accepts setting of any of the L<OPTIONS|OPTIONS>.

=cut

#-----------------------------------------------
sub new {
#-----------------------------------------------
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

=head2 Stats

=head3 ztest

 $zed->ztest(observed => 'number', expected => 'number', variance => 'number (non-zero)', samplings => 'positive integer')

I<Alias>: test

You supply the C<observed> and c<expected> values of your statistic, and the C<variance> or C<stdev> or just C<error> of any kind; ... and the number of C<samplings> ("sample-size", I<N>, "trials", etc.).

Optionally specify a logical value for L<ccorr|ccorr> for performing the continuity-correction to the observed deviation, and a value of either 1 or 2 to specify the L<tails|tails> for reading off the probability associated with the I<z>-value. When passed to this function, values of L<ccorr|ccorr> and L<tails|tails> become the values used in this and subsequent tests.

When called in array context, returns an array consisting of the z-statistic ("standard normal deviate"), its probability, the observed deviation (the difference between the observed and expected values of your statistic), and the standard deviation (the square-root of the variance supplied).

If you only want the I<z>-value, then expect a string:

 $z_value = $zed->test(...)

The basic formula is the basic:

=for html <p>&nbsp;&nbsp;&nbsp;<i>Z<sub><o>&times;</o></sub></i> = ( <i><o>&times;</o></i> &ndash; <i>&micro;</i> ) / SD / &not;/<i>n</i></p>

=cut

#-----------------------------------------------
sub ztest {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    $self->{$_} = $args->{$_} foreach keys %{$args};
    my ($z, $p, $obs_dev, $exp_dev) = ();
    croak "Need to define observed ($self->{'observed'}) and expected ($self->{'expected'}) values for ztest" if !hascontent($self->{'observed'}) || !hascontent($self->{'expected'});

    $obs_dev = $self->{'observed'} - $self->{'expected'};
    $obs_dev = ccorr($obs_dev) if $self->{'ccorr'};
    $exp_dev = _exp_dev($self, 1);
    return undef if !$exp_dev;
    $z = $obs_dev / $exp_dev;
 
    $p = $self->z_2_p($z, $self->{'tails'});
    $z = sprintf('%.' . $self->{'precision_s'} . 'f', $z) if $self->{'precision_s'};

    return wantarray ? ($z, $p, $obs_dev, $exp_dev) : $z;
}
*test = \&ztest;

=head3 zscore

 $zed->zscore(observed => 'number', expected => 'number', variance => 'number (non-zero)')

I<Alias>: score

You supply the C<observed> and C<expected> values of your statistic, and the C<variance> or C<stdev> or just C<error> of any kind.

Optionally specify a logical value for L<ccorr|ccorr> for performing the continuity-correction to the observed deviation, and a value of either 1 or 2 to specify the L<tails|tails> for reading off the probability associated with the I<z>-value. When passed to this function, values of L<ccorr|ccorr> and L<tails|tails> become the values used in this and subsequent tests.

When called in array context, returns an array consisting of the z-statistic ("standard score"), its probability, the observed deviation (the difference between the observed and expected values of your statistic), and the standard deviation (the square-root of the variance supplied).

If you only want the I<z>-value, then expect a string:

 $z_value = $zed->score(...)

The basic formula is the basic:

=for html <p>&nbsp;&nbsp;&nbsp;<i>Z</i> = ( <i>&times;</i> &ndash; <i><o>X</o></i> ) / SD</p>

where I<X> is the expected value (mean, etc.).

=cut

#-----------------------------------------------
sub zscore {
#-----------------------------------------------
    my $self = shift;
	my $args = ref $_[0] ? $_[0] : {@_};
    $self->{$_} = $args->{$_} foreach keys %{$args};
    my ($z, $p, $obs_dev, $exp_dev) = ();
    croak "Need to define observed ($self->{'observed'}) and expected ($self->{'expected'}) values for zscore" if !hascontent($self->{'observed'}) || !hascontent($self->{'expected'});

    $obs_dev = $self->{'observed'} - $self->{'expected'};
    $obs_dev = ccorr($obs_dev) if $self->{'ccorr'};
    $exp_dev = _exp_dev($self, 0);
    return undef if !$exp_dev;
    $z = $obs_dev / $exp_dev;
 
    $p = $self->z_2_p($z, $self->{'tails'});
    $z = sprintf('%.' . $self->{'precision_s'} . 'f', $z) if $self->{'precision_s'};
    
    return wantarray ? ($z, $p, $obs_dev, $exp_dev) : $z;
}
*score = \&zscore;

=head3 z2p

 $p = $zed->z2p($z)
 $p = $zed->z2p($z, 1)

I<Alias>: C<z_2_p>, C<p_value>

Send a I<z>-value, get its associated I<p>-value, 2-tailed by default, or depending on what the value of $zed->{'tails'} is, or what is sent as the second argument, if anything. Uses L<Math::Cephes|Math::Cephes> C<ndtr>.

=cut

#-----------------------------------------------
sub z2p {
#-----------------------------------------------
    my ($self, $z, $tails) = @_;
    #my $p = (1 - ndtr(abs($z))); # Math::Cephes function
	my $p = Statistics::Distributions::uprob (abs($z));
	#print "$z\t$p\n"; 
	
    $p *= 2 if $tails and $tails == 2;
	$p = Statistics::Distributions::precision_string($p);
    $p = 1 if $p > 1;
    $p = 0 if $p < 0;
    $p = sprintf('%.' . $self->{'precision_p'} . 'f', $p) if $self->{'precision_p'};
    return $p;
}
*z_2_p = \&z2p;
*p_value = \&z2p;

=head3 p2z

 $z = $zed->p2z($p)
 $z = $zed->p2z($p, 2) # 2-tailed probability

I<Alias>: C<p_2_z>

Returns the I<z>-value associated with a I<p>-value, using L<Math::Cephes|Math::Cephes> C<ndtri> function ("phi"). If the class attribute C<tails> has been set to 2, or if the second argument equals 2 (indicating a 2-tailed probability value), the I<p>-value is firstly divided by 2. The absolute value of C<z> is always returned. However, if C<p> equals zero or is undefined, returns undef.

=cut

#-----------------------------------------------
sub p2z {
#-----------------------------------------------
    my ($self, $p, $tails) = @_;

    my $z;
    
    if (!$p) {
        $z = undef;
    }
    elsif ($p == 1)  {
        $z = 0;
    }
    else {
        $p = $p/2 if ($tails and $tails == 2) || $self->{'tails'} == 2;
        $z = abs(ndtri($p)); # Math::Cephes function
    }

    return $z;
}
*p_2_z = \&p2z;


=head3 z2chi

 $chi = $zed->z2chi($z) 

I<Alias>: C<z_2_chi>

Send a I<z>-value, get back a I<chi>-value (the square of the thing ...).

=cut

#-----------------------------------------------
sub z2chi {
#-----------------------------------------------
    my ($self, $z) = @_;
    return $z**2;
}
*z_2_chi = \&z2chi;

=head3 chi2z

 $z = $zed->chi2z($chi)
 $z = $zed->chi2z($chi, $df) 

I<Alias>: C<chi_2_z>

Send a I<chi>-value, get back a I<z>-value. If no degrees-of-freedom after the chi-value (as above), you get back the square-root of the thing. Otherwise:

=for html <p>&nbsp;&nbsp;&nbsp;<i>Z</i> ~= &not;/(2CHI) &ndash;  &not;/(2DF &ndash; 1)</p>

Note that a negative I<Z>-value returned by this method indicates that the chi-square is I<smaller> than expected for so many dfs; the I<p>-value for this I<Z>-value indicates how much smaller.

=cut

#-----------------------------------------------
sub chi2z {
#-----------------------------------------------
    my ($self, $chi, $df) = @_;
    if ($df) {
        return sqrt(2 * $chi) - sqrt( 2 * $df - 1);
    }
    else {
        return sqrt($chi);
    }
}
*chi_2_z = \&chi2z;

# Apply the continuity correction to the deviation, e.g. of (x - MCE) in numerator of z-score, if variance is calculated binomially (as hit/miss, not continuous):
sub ccorr {
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

=head3 vals2z

 $z_aref = $zed->vals2z($chi)
 $z_aref = $zed->vals2z($chi, $df) 

I<Alias>: C<vals_2_z>

Given an array (reference) of numberical values, returns them as an array of standard scores (mean = 0, SD = 1).

=cut

sub vals2z {
   my $self = shift;
   my $vals = ref($_[0]) ? shift : [$_[0]];
   
   if (ref $vals) {
     my @scores = ();
     my $mean = mean(@{$vals});
     my $var = variance(@{$vals});
     foreach (@{$vals}) {
        push @scores, $var ? ($_ - $mean) / sqrt($var) : 0;
     }
     return\@scores;
   }
}
*vals_2_z = \&vals2z;

=head2 dump

Prints to STDOUT a line giving the z_value and p_value.

=cut

sub dump {
    my ($self, $series) = @_;
    print "Z = $self->{'z_value'}, $self->{'tails'}p = $self->{'p_value'}\n";
}

# Calc expected deviation, ensuring valid denominator:
sub _exp_dev {
    my ($self) = @_;
    my $exp_dev;
    my $stdev = $self->{'stdev'} ? 
      $self->{'stdev'} : 
        $self->{'variance'} ?
            sqrt($self->{'variance'}) :
            $self->{'error'} ?
                $self->{'error'} :
                    croak "Need a stdev or variance value for z-testing";##return undef;##
    if ($self->{'test'}) {
        croak "Need non-zero and positive number of samplings/units" if ! $self->{'samplings'} or $self->{'samplings'} < 0;
        $exp_dev = $stdev / sqrt($self->{'samplings'});
    }
    else {
        $exp_dev = $stdev;
    }
    return $exp_dev;
}

=head2 Series testing

A means to aggregate results from multiple tests is supported. Three methods are presently used to effect this.

=head3 series_init

Clears any already accumulated data from previous tests.

=head3 series_update

 $zed->series_update() # use any cached values of "observed", "expected" and "variance"
 $zed->series_update(variance => 'number', expected => 'number', observed => 'number') # supply own in a hash

Called once you have performed a test on a sample. It caches the observed, expectation and variance values from the test.

=head3 series_test

Sums the observed, expectation and variance values from all the tests updated to the series since calling L<series_init|series_init>, and produces a I<z>-value from these sums. It returns nothing in particular, but the following statement shows how the series values can be accessed.

 print "Series summed runs: 
    expected = ", $zed->{'series'}->{'expected'}, " 
    observed = ", $zed->{'series'}->{'observed'}," 
    z = $zed->{'series'}->{'z_value'}, $zed->{'tails'}p = $zed->{'series'}->{'p_value'}\n";

=cut

#-----------------------------------------------
sub series_init {
#-----------------------------------------------
    my ($self, %args) = @_;
    $self->{'series_stat'}->{$_} = Statistics::Descriptive::Sparse->new() foreach qw/observed expected variance z_value/;
    $self->{'series'}->{$_} = 0 foreach qw/observed expected variance z_value p_value /;
   # $self->{'series'}->{'p_value'} = 1;
  # $self->{'series'}->{'expected'} = 1;
}

#-----------------------------------------------
sub series_update {
#-----------------------------------------------
    my $self = shift;
    my $args = ref $_[0] ? $_[0] : {@_};
    if (! defined $args->{'variance'}) {
        if (defined $args->{'stdev'}) {
            $args->{'variance'} = $args->{'stdev'}**2;
        }
        else {
            croak "No variance offered in series_update for calculating deviation";
        }
    }
 	#print("$_\t$args->{$_}\n") foreach qw/observed expected variance z_value/;
	foreach (qw/observed expected variance z_value/) {
    	$self->{'series_stat'}->{$_}->add_data($args->{$_}) if defined $args->{$_};
	}
}

#-----------------------------------------------
sub series_test {
#-----------------------------------------------
    my ($self) = (shift);
    croak 'No data for series-testing appear to have been loaded; maybe you need to call series_update with some data' if ! $self->{'series_stat'}->{'variance'}->count();
    my $o = $self->{'series_stat'}->{'observed'}->sum();
    my $e = $self->{'series_stat'}->{'expected'}->sum();
    my $v = $self->{'series_stat'}->{'variance'}->sum();
    my ($z, $pz, $obs_dev, $stdev) = $self->zscore(
        observed => $o,
        expected => $e,
        variance => $v);
    $self->{'series'}->{'observed'} = $o;
    $self->{'series'}->{'expected'} = $e;
    $self->{'series'}->{'z_value'} = $z;# || 0;
    $self->{'series'}->{'p_value'} = $pz;# || 1;
    #$self->{'series'}->{'r_value'} = $self->z_2_r($z);
    $self->{'series'}->{'obs_dev'} = $obs_dev;
    $self->{'series'}->{'stdev'} = $stdev;
    $self->{'series'}->{'variance'} = $v;
    $self->{'series'}->{'samplings'} = $self->{'series_stat'}->{'variance'}->count();
	if ($self->{'series_stat'}->{'z_value'}->count()){
    	$self->{'series'}->{'stouffer_z'} = $self->{'series_stat'}->{'z_value'}->sum() / sqrt($self->{'series_stat'}->{'z_value'}->count());
    	$self->{'series'}->{'stouffer_p'} = $self->p_value($self->{'series'}->{'stouffer_z'});
	}
    return $self;
}

=head2 series_dump

Prints to STDOUT a line giving the z_value and p_value for the series.

=cut

sub series_dump {
    my ($self) = @_;
    print "Z (N = $self->{'series'}->{'samplings'}) = $self->{'series'}->{'z_value'}, $self->{'tails'}p = $self->{'series'}->{'p_value'}\n";
    
}

1;
__END__

=head1 OPTIONS

The following can be set in the call to L<new|new> or L<test|test> or L<score|score>.

=head2 ccorr

Apply the continuity correction. Default = 0.

=head2 tails

Tails from which to assess the association I<p>-value (1 or 2). Default = 2.

=head2 precision_s

Precision of the I<z>-value (the statistic). Default = 2.

=head2 precision_p

Precision of the associated I<p>-value. Default = 0 - you get all decimal values available.

=head1 SEE ALSO

L<Math::Cephes|Math::Cephes> : the C<ndtr> and C<ndtri> functions are used here for determining probability or I<z>-values, respectively. (This is only in preference to using Statistics::Distributions so that whether we want a I<p>- or I<z>value, the same class of methods is used.)

L<Statistics::Distributions|Statistics::Distributions> : the C<precision_string> method is used here for reporting probability.

L<Statistics::Sequences|Statistics::Sequences> : for application of this module.

=head1 TO DO/BUGS

Other distributions.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2009 R Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut
