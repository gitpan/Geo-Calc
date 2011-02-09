#!/usr/bin/perl -T
#
# Test geo calc
#

use strict;
use warnings;

use lib qw(. lib tests);

use Test::More tests => 10;

use_ok 'Geo::Calc';

my $gc = Geo::Calc->new( lat => 40.417875, lon => -3.710205 );

is_deeply( $gc->boundry_box( 1, 1, -6 ), { 'lat_max' => '40.413378', 'lon_max' => '-3.704299', 'lon_min' => '-3.716112', 'lat_min' => '40.422371' }, 'boundry box' );

is( $gc->distance_to( { lat => 40.422371, lon => -3.704298 } ), 0.707106, 'distance' );

is( $gc->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 }, -6 ), 0.707095, 'rhumb distance' );

is( $gc->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 ), 314.004851, 'initial bearing' );

is( $gc->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 } ), 45.006766, 'rhumb bearing' );

is( $gc->final_bearing_to( { lat => 40.422371, lon => -3.704298 } ), 134.004851, 'final bearing' );

is_deeply( $gc->rhumb_destination_point( 30, 1, -6 ), { 'lat' => '40.425663', 'lon' => '-3.704298' }, 'rhumb destination' );

is_deeply( $gc->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 ) , { 'lat' => '40.420123', 'lon' => '-3.707252' }, 'midpoint' );

is_deeply( $gc->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180, -6 ), { 'lat' => '40.417875', 'lon' => '-3.704298' }, 'intersection' );
