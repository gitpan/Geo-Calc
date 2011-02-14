#!/usr/bin/perl -T
#
# Test geo calc
#

use strict;
use warnings;

use lib qw(. lib tests);

use Test::More tests => 12;

use_ok 'Geo::Calc';

my $gc = Geo::Calc->new( lat => 40.417875, lon => -3.710205 );

is_deeply( $gc->boundry_box( 1000, 1000, -6 ), { 'lat_max' => '40.422378', 'lon_max' => '-3.704314', 'lon_min' => '-3.716096', 'lat_min' => '40.413372' }, 'boundry box 1/1 km' );

is_deeply( $gc->boundry_box( 6000, 8000 ), { 'lat_max' => '40.453897', 'lon_max' => '-3.674857', 'lon_min' => '-3.745553', 'lat_min' => '40.381853' }, 'boundry box 6/8 km' );

is( $gc->distance_to( { lat => 40.422371, lon => -3.704298 } ), 0.707106, 'distance' );

is( $gc->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 }, -6 ), 0.707095, 'rhumb distance' );

is( $gc->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 ), 314.004851, 'initial bearing' );

is( $gc->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 } ), 45.006766, 'rhumb bearing' );

is( $gc->final_bearing_to( { lat => 40.422371, lon => -3.704298 } ), 134.004851, 'final bearing' );

is_deeply( $gc->rhumb_destination_point( 30, 1, -6 ), { 'lat' => '40.425663', 'lon' => '-3.704298' }, 'rhumb destination' );

is_deeply( $gc->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 ) , { 'lat' => '40.420123', 'lon' => '-3.707252' }, 'midpoint' );

is_deeply( $gc->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180, -6 ), { 'lat' => '40.417875', 'lon' => '-3.704298' }, 'intersection' );

is_deeply( $gc->distance_at(), { m_lon => '84871.014948', m_lat => '111042.645811' }, 'distance at latitude' );
