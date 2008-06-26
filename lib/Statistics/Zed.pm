package Statistics::Zed;

use 5.008008;
use strict;
use warnings;
use Carp qw(croak);
use Class::OOorNO qw(coerce_array);
use vars qw($VERSION);
$VERSION = 0.01;

use Statistics::Distributions;

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sub new {
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    my $class = shift;
    my $args = Class::OOorNO::coerce_array(@_);
	my $self = {};
	bless $self, $class;

    # Set default values:
    $self->{'tails'} = 2;
    $self->{'s_precision'} = 2;
    $self->{'p_precision'} = 0;
    $self->{'ccorr'} = 0;
    ##$self->{'distribution'} = 'norm';
    
    if (scalar keys %{$args}) {
        foreach (keys %{$args}) {
            $self->{$_} = $args->{$_};
        }
    }

	return $self;
}

#-----------------------------------------------
sub zscore {
#-----------------------------------------------
    my $self = shift;
    my $args = Class::OOorNO::coerce_array(@_);
    foreach (keys %{$args}) {
        $self->{$_} = $args->{$_};
    }
    my ($z, $p, $obs_dev, $exp_dev) = ();

    $obs_dev = $self->{'observed'} - $self->{'expected'};
    $obs_dev = _continuity_correct($obs_dev) if $self->{'ccorr'} && $obs_dev;
    $exp_dev = _exp_dev($self, 0);
    $z = $obs_dev / $exp_dev;
 
    $p = $self->p_value($z);
    $z = sprintf('%.' . $self->{'s_precision'} . 'f', $z) if $self->{'s_precision'};
    
    return wantarray ? ($z, $p, $obs_dev, $exp_dev) : $z;
}

#-----------------------------------------------
sub ztest {
#-----------------------------------------------
    my $self = shift;
    my $args = Class::OOorNO::coerce_array(@_);
    foreach (keys %{$args}) {
        $self->{$_} = $args->{$_};
    }
    my ($z, $p, $obs_dev, $exp_dev) = ();

    $obs_dev = $self->{'observed'} - $self->{'expected'};
    $obs_dev = _continuity_correct($obs_dev) if $self->{'ccorr'} && $obs_dev;
    $exp_dev = _exp_dev($self, 1);
    $z = $obs_dev / $exp_dev;
 
    $p = $self->p_value($z);
    $z = sprintf('%.' . $self->{'s_precision'} . 'f', $z) if $self->{'s_precision'};
    
    return wantarray ? ($z, $p, $obs_dev, $exp_dev) : $z;
}

#-----------------------------------------------
sub p_value {
#-----------------------------------------------
    my ($self, $z) = @_;
    my $p = Statistics::Distributions::uprob(abs($z));
    $p *= 2 if ! $self->{'tails'} or $self->{'tails'} == 2 || ($self->{'tails'} != 1 && $self->{'tails'} != 2);
	$p = Statistics::Distributions::precision_string($p);
    $p = 1 if $p > 1;
    $p = 0 if $p < 0;
    $p = sprintf('%.' . $self->{'p_precision'} . 'f', $p) if $self->{'p_precision'};
    return $p;
}

#-----------------------------------------------
sub r_2_z {
#-----------------------------------------------
    my ($self, $r) = @_;
    croak 'Null or invalid r_value' if !defined $r or $r > 1 or $r < -1;
    my $z = .5 * ( log( (1 + $r) / (1 - $r) ) );
    return $z;
}

#-----------------------------------------------
sub z_2_r {
#-----------------------------------------------
   my ($self, $z) = @_;
   return defined $z ? (exp(2 * $z) - 1) / (exp(2 * $z) + 1) : undef;
}

#-----------------------------------------------
sub chi_2_z {
#-----------------------------------------------
    my ($self, $chi) = @_;
    return sqrt($chi);
}

# Apply the continuity correction to the deviation, e.g. of (x - MCE) in numerator of z-score, if variance is calculated binomially (as hit/miss, not continuous):
sub _continuity_correct {
   my $zed = shift;
   my $c = abs($zed) - .5;
   $c *= -1 if $zed < 0;
   return $c;
}

sub _exp_dev {
    my ($self) = @_;
    my $exp_dev;
    my $std_dev = $self->{'std_dev'} ? 
      $self->{'std_dev'} : 
        $self->{'variance'} ?
            sqrt($self->{'variance'}) :
            croak "Need a std_dev or variance value";
    if ($self->{'test'}) {
        croak "Need non-zero and positive number of samples/units" if ! $self->{'samples'} or $self->{'samples'} < 0;
        $exp_dev = $std_dev / sqrt($self->{'samples'});
    }
    else {
        $exp_dev = $std_dev;
    }
    return $exp_dev;
}

# --------------------
# Aliases
# --------------------
*test = \&ztest;
*score = \&zscore;

1;
__END__

=head1 NAME

Statistics::Zed - Basic ztest/zscore, Fisher's r-to-z, z-to-r, et al.

=head1 VERSION

This is documentation for Version 0.01 of Statistics::Zed (2008.06.22).

=head1 SYNOPSIS

 use Statistics::Zed;

 my $zed = Statistics::Zed->new(
    ccorr    => 1,
    tails    => 2,
    s_precision => 5,
    p_precision => 5,
 );

 my ($z_value, $p_value, $observed_deviation, $standard_deviation) = 
    $zed->score(
        observed => $obs,
        expected => $exp,
        variance => $variance,
 );

 my $deviate = $zed->test(
        observed => $obs,
        expected => $exp,
        variance => $variance,
        samples => $samples,
 );

 my $p_value = $zed->p_value($deviate);

=head1 DESCRIPTION

Calculates a z-statistic: the ratio of an observed deviation to a standard deviation. Purpose is simply to support L<Statistics::Sequences|Statistics::Sequences>, but with some standalone utilities.

=head1 METHODS

=head2 new

 $zed = Statistics::Zed->new();

Returns a Statistics::Zed object. Accepts setting of any of the L<OPTIONS|OPTIONS>.

=head2 ztest

 $zed->ztest(observed => 'number', expected => 'number', variance => 'number (non-zero)', samples => 'number (non-zero)')

I<Alias>: test

You supply the I<observed> and I<expected> values of your statistic, and the I<variance> (or I<std_dev>), and the number of I<samples> ("sample-size", I<N>, etc.).

Optionally specify a logical value for L<ccorr|ccorr> for performing the continuity-correction to the observed deviation, and a value of either 1 or 2 to specify the L<tails|tails> for reading off the probability associated with the z_value.

When called in array context, returns an array consisting of the z-statistic ("standard normal deviate"), its probability, the observed deviation (the difference between the observed and expected values of your statistic), and the standard deviation (the square-root of the variance supplied).

If you only want the z-value, then expect a string:

 $z_value = $zed->test(...)

The basic formula is the basic:

=for html <p>&nbsp;&nbsp;&nbsp;<i>Z<sub><o>&times;</o></sub></i> = ( <i><o>&times;</o></i> &ndash; <i>&micro;</i> ) / SD / &not;/<i>n</i></p>

=head2 zscore

 $zed->zscore(observed => 'number', expected => 'number', variance => 'number (non-zero)')

I<Alias>: score

You supply the I<observed> and I<expected> values of your statistic, and the I<variance> (or I<std_dev>).

Optionally specify a logical value for L<ccorr|ccorr> for performing the continuity-correction to the observed deviation, and a value of either 1 or 2 to specify the L<tails|tails> for reading off the probability associated with the z_value.

When called in array context, returns an array consisting of the z-statistic ("standard score"), its probability, the observed deviation (the difference between the observed and expected values of your statistic), and the standard deviation (the square-root of the variance supplied).

If you only want the z-value, then expect a string:

 $z_value = $zed->score(...)

The basic formula is the basic:

=for html <p>&nbsp;&nbsp;&nbsp;<i>Z</i> = ( <i>&times;</i> &ndash; <i><o>X</o></i> ) / SD</p>

where I<X> is the expected value (mean, etc.).

=head2 p_value

 $p = $zed->p_value($z)

Send a z-value, get its associated p-value, 2-tailed by default.

=head2 r_2_z

 $z = $zed->r_2_z($r)

Performs the Fisher r-to-z transformation. Send a correlation coefficient - get back a z-value.

=head2 z_2_r

 $r = $zed->z_2_r($z)

Send a z-value - get back a correlation coefficient.

=head2 chi_2_z

 $z = $zed->chi_2_z($chi) 

Send a chi-value, get back a z-value (the square-root of the thing ...).

=head1 OPTIONS

The following can be set in the call to L<new|new> or L<test|test> or L<score|score>.

=head2 ccorr

Apply the continuity correction. Default = 0.

=head2 tails

Tails from which to assess the association p_value (1 or 2). Default = 2.

=head2 s_precision

Precision of the z_value. Default = 2.

=head2 p_precision

Precision of the associated p_value. Default = 0.

=head1 SEE ALSO

L<Statistics::Distributions|Statistics::Distributions> : the C<uprob> and C<precision_string> methods are here used for calculating and reporting probability.

=head1 TO DO/BUGS

Other distributions.

=head1 AUTHOR

Roderick Garton, E<lt>rgarton@utas_DOT_edu_DOT_auE<gt>

=head1 COPYRIGHT/LICENSE/DISCLAIMER

Copyright (C) 2008 Roderick Garton 

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself, either Perl version 5.8.8 or, at your option, any later version of Perl 5 you may have available. 

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=cut

1;
