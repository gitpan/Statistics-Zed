#!perl -T

use Test::More;
eval "use Test::Pod::Coverage 1.04";
plan skip_all => "Test::Pod::Coverage 1.04 required for testing POD coverage" if $@;
all_pod_coverage_ok({trustme => ['z_2_p', 'p_value', 'p_2_z', 'r_2_z', 'z_2_r', 'z_2_chi', 'chi_2_z', 'vals_2_z', 'score', 'test']});

