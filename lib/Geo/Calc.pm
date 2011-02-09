# Copyrights 2011 by Sorin Pop.
# For other contributors see ChangeLog.
# See the manual pages for details on the licensing terms.

package Geo::Calc;

use vars '$VERSION';
$VERSION = '0.03';

use Moose;
use MooseX::FollowPBP;
use MooseX::Method::Signatures;

use Math::Trig qw(:pi asin acos tan deg2rad rad2deg);
use Math::BigFloat;

=head1 NAME

Geo::Calc - simple geo calculator for points and distances

=head1 SYNOPSIS

 use Geo::Calc;

 my $gc            = Geo::Calc->new( lat => 40.417875, lon => -3.710205 );
 my $distance      = $gc->distance_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $brng          = $gc->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $f_brng        = $gc->final_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $midpoint      = $gc->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $destination   = $gc->destination_point( 90, 1, -6 ); # distance in km
 my $bbox          = $gc->boundry_box( 3, 4, -6 ); # in km
 my $r_distance    = $gc->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $r_brng        = $gc->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
 my $r_destination = $gc->rhumb_destination_point( 30, 1, -6 );
 my $point         = $gc->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180, -6 );

=head1 DESCRIPTION

 C<Geo::Calc> implements a variety of calculations for latitude/longitude points

 All these formulare are for calculations on the basis of a spherical earth
(ignoring ellipsoidal effects) which is accurate enough* for most purposes.

 [ In fact, the earth is very slightly ellipsoidal; using a spherical model
gives errors typically up to 0.3% ].

=head1 Geo::Calc->new()

 $gc = Geo::Calc->new( lat => 40.417875, lon => -3.710205 ); # Somewhere in Madrid
 $gc = Geo::Calc->new( lat => 51.503269, lon => 0 ); # The O2 Arena in London

Creates a new Geo::Calc object from a latitude and longitude. The default
precision is -6 for all functions meaning that 6 deciamls

Returns ref to a Geo::Calc object.

=head2 Parameters

=over 4

=item lat

C<>=> latitude of the point ( required )

=item lon

C<>=> longitude of the point ( required )

=item radius

C<>=> earth radius in km ( defaults to 6371 )

=back

=cut

has 'lat' => (
    is       => 'ro',
    isa      => 'Num',
    required => 1,
);

has 'lon' => (
    is       => 'ro',
    isa      => 'Num',
    required => 1,
);

has 'radius' => (
    is       => 'ro',
    isa      => 'Num',
    default  => '6371',
);

=head1 METHODS

=head2 distance_to

 $gc->distance_to( $point, $precision )
 $gc->distance_to( { lat => 40.422371, lon => -3.704298 }, -6 )

This uses the "haversine" formula to calculate great-circle distances between
the two points - that is, the shortest distance over the earth's surface - 
giving an `as-the-crow-flies` distance between the points (ignoring any hills!)

The haversine formula `remains particularly well-conditioned for numerical
computation even at small distances` - unlike calculations based on the spherical
law of cosines. It was published by R W Sinnott in Sky and Telescope, 1984,
though known about for much longer by navigators. (For the curious, c is the
angular distance in radians, and a is the square of half the chord length between
the points).

Returns with the distance using the precision defined or -6
( -6 = 6 decimals after the dot ( eg 4.000001 ) )

=cut

method distance_to( HashRef[Num] $point!, Int $precision? = -6 ) returns (Num) {
    my $lat1 = deg2rad( $self->get_lat() );
    my $lon1 = deg2rad( $self->get_lon() );
    my $lat2 = deg2rad( $point->{lat} );
    my $lon2 = deg2rad( $point->{lon} );

    my $t = sin( ($lat2 - $lat1)/2 ) ** 2 + ( cos( $lat1 ) ** 2 ) * ( sin( ( $lon2 - $lon1 )/2 ) ** 2 );
    my $d = $self->get_radius * ( 2 * atan2( sqrt($t), sqrt(1-$t) ) );

    return $self->_precision( $d, $precision );
}

=head2 bearing_to

 $gc->bearing_to( $point, $precision );
 $gc->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );

In general, your current heading will vary as you follow a great circle path
(orthodrome); the final heading will differ from the initial heading by varying
degrees according to distance and latitude (if you were to go from say 35N,45E
(Baghdad) to 35N,135E (Osaka), you would start on a heading of 60 and end up on
a heading of 120!).

This formula is for the initial bearing (sometimes referred to as forward
azimuth) which if followed in a straight line along a great-circle arc will take
you from the start point to the end point

Returns the (initial) bearing from this point to the supplied point, in degrees
with the specified pricision

see http://williams.best.vwh.net/avform.htm#Crs

=cut

method bearing_to( HashRef[Num] $point!, Int $precision? = -6 ) returns (Num) {
    my $lat1 = deg2rad( $self->get_lat() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlon = deg2rad( $self->get_lon() - $point->{lon} );

    my $brng = atan2( sin( $dlon ) * cos( $lat2 ), ( cos( $lat1 ) * sin( $lat2 ) ) - ( sin( $lat1 ) * cos( $lat2 ) * cos( $dlon ) ) );

    return $self->_ib_precision( rad2deg( $brng ), $precision );
}

=head2 final_bearing_to

 my $f_brng = $gc->final_bearing_to( $point, $precision );
 my $f_brng = $gc->final_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );

Returns final bearing arriving at supplied destination point from this point;
the final bearing will differ from the initial bearing by varying degrees
according to distance and latitude

=cut

method final_bearing_to( HashRef[Num] $point!, Int $precision? = -6 ) returns (Num) {
    my $lat1 = deg2rad( $self->get_lat() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlon = deg2rad( $self->get_lon() - $point->{lon} );

    my $brng = atan2( sin( $dlon ) * cos( $lat2 ), ( cos( $lat1 ) * sin( $lat2 ) ) - ( sin( $lat1 ) * cos( $lat2 ) * cos( $dlon ) ) );

    return $self->_fb_precision( rad2deg( $brng ), $precision );
}

=head2 midpoint_to

 $gc->midpoint_to( $point, $precision );
 $gc->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 );

Returns the midpoint along a great circle path between the initial point and
the supplied point.

see http://mathforum.org/library/drmath/view/51822.html for derivation

=cut

method midpoint_to( HashRef[Num] $point!, Int $precision? = -6 ) returns (HashRef[Num]) {
    my $lat1 = deg2rad( $self->get_lat() );
    my $lon1 = deg2rad( $self->get_lon() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlon = deg2rad( $point->{lon} - $self->get_lon() );

    my $bx = cos( $lat2 ) * cos( $dlon );
    my $by = cos( $lat2 ) * sin( $dlon );

    my $lat3 = atan2( sin( $lat1 ) + sin ( $lat2 ), sqrt( ( ( cos( $lat1 ) + $bx ) ** 2 ) + ( $by ** 2 ) ) );
    my $lon3 = $lon1 + atan2( $by, cos( $lat1 ) + $bx );
    $lon3 -= pi2 while( $lon3 > pi );
    $lon3 += pi2 while( $lon3 <= -(pi) );

    return { lat => $self->_precision( rad2deg($lat3), $precision ), lon => $self->_precision( rad2deg($lon3), $precision ) };
}

=head2 destination_point

 $gc->destination_point( $bearing, $distance, $precision );
 $gc->destination_point( 90, 1, -6 ); # distance in km

Returns the destination point from this point having travelled the given
distance (in km) on the given initial bearing (bearing may vary before
destination is reached)

see http://williams.best.vwh.net/avform.htm#LL

=cut

method destination_point( Num $brng!, Num $dist!, Int $precision? = -6 ) returns (HashRef[Num]) {
    $dist = $dist / $self->get_radius();
    $brng = deg2rad( $brng );
    my $lat1 = deg2rad( $self->get_lat() );
    my $lon1 = deg2rad( $self->get_lon() );

    my $lat2 = asin( sin( $lat1 ) * cos( $dist ) + cos( $lat1 ) * sin( $dist ) * cos( $brng ) );
    my $lon2 = $lon1 + atan2( sin( $brng ) * sin( $dist ) * cos( $lat1 ), cos( $dist ) - sin( $lat1 ) * sin ( $lat2 ) );

    # Normalize longitude so that its in range -PI to +PI
    $lon2 -= pi2 while( $lon2 > pi );
    $lon2 += pi2 while( $lon2 <= -(pi) );

    return { lat => $self->_precision( rad2deg($lat2), $precision ), lon => $self->_precision( rad2deg($lon2), $precision ) };
}

=head2 boundry_box

 $gc->boundry_box( $width, $height, $precision ); # in km
 $gc->boundry_box( 3, 4, -6 ); # in km

Returns the boundry box min/max having the initial point defined as the center
of the boundry box, given the widht and height

=cut

method boundry_box( Num $width!, Num $height!, Int $precision? = -6 ) returns (HashRef[Num]) {
    $height = $width if( !defined( $height ) );
    my @points = ();
    push @points, $self->destination_point( 315 , sqrt( ( ( $height/2 ) ** 2 ) + ( ( $width/2 ) ** 2 ) ), $precision );
    push @points, $self->destination_point( 135 , sqrt( ( ( $height/2 ) ** 2 ) + ( ( $width/2 ) ** 2 ) ), $precision );

    return {
        lat_min => $points[0]->{lat},
        lon_min => $points[0]->{lon},
        lat_max => $points[1]->{lat},
        lon_max => $points[1]->{lon},
    };
}

=head2 rhumb_distance_to

 $gc->rhumb_distance_to( $point, $precision );
 $gc->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 }, -6 );

Returns the distance from this point to the supplied point, in km, travelling
along a rhumb line.

A 'rhumb line' (or loxodrome) is a path of constant bearing, which crosses all
meridians at the same angle.

Sailors used to (and sometimes still) navigate along rhumb lines since it is
easier to follow a constant compass bearing than to be continually adjusting
the bearing, as is needed to follow a great circle. Rhumb lines are straight
lines on a Mercator Projection map (also helpful for navigation).

Rhumb lines are generally longer than great-circle (orthodrome) routes. For
instance, London to New York is 4% longer along a rhumb line than along a
great circle . important for aviation fuel, but not particularly to sailing
vessels. New York to Beijing . close to the most extreme example possible
(though not sailable!) . is 30% longer along a rhumb line.

see http://williams.best.vwh.net/avform.htm#Rhumb

=cut

method rhumb_distance_to( HashRef[Num] $point!, Int $precision? = -6 ) returns (Num) {
    my $lat1 = deg2rad( $self->get_lat() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlat = deg2rad( $point->{lat} - $self->get_lat() );
    my $dlon = abs( deg2rad( $point->{lon} - $self->get_lon() ) );

    my $dphi = log( tan( $lat2/2 + pip4 ) / tan( $lat1/2 + pip4 ) );
    my $q = ( $dphi != 0 ) ? $dlat/$dphi : cos($lat1);# E-W line gives dPhi=0
    $dlon = pi2 - $dlon if ( $dlon > pi );

    my $dist = sqrt( ( $dlat ** 2 ) + ( $q ** 2 ) * ( $dlon ** 2 ) ) * $self->get_radius();

    return $self->_precision( $dist, $precision );
}

=head2 rhumb_bearing_to

 $gc->rhumb_bearing_to( $point, $precision );
 $gc->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );

Returns the bearing from this point to the supplied point along a rhumb line,
in degrees

=cut

method rhumb_bearing_to( HashRef[Num] $point!, Int $precision? = -6 ) returns (Num) {
    my $lat1 = deg2rad( $self->get_lat() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlon = deg2rad( $point->{lon} - $self->get_lon() );


    my $dphi = log( tan( $lat2/2 + pip4 ) / tan( $lat1/2 + pip4 ) );
    if( abs( $dlon ) > pi ) {
        $dlon = ( $dlon > 0 ) ? -(pi2-$dlon) : (pi2+$dlon);
    }

    return $self->_ib_precision( rad2deg( atan2( $dlon, $dphi ) ), $precision );
}

=head2 rhumb_destination_point

 $gc->rhumb_destination_point( $brng, $distance, $precision );
 $gc->rhumb_destination_point( 30, 1, -6 );

Returns the destination point from this point having travelled the given distance
(in km) on the given bearing along a rhumb line.

=cut

method rhumb_destination_point( Num $brng!, Num $dist!, Int $precision? = -6 ) returns (HashRef[Num]) {
    my $d = $dist / $self->get_radius();
    my $lat1 = deg2rad( $self->get_lat() );
    my $lon1 = deg2rad( $self->get_lon() );
    $brng = deg2rad( $brng );

    my $lat2 = $lat1 + $d * cos( $brng );
    my $dlat = $lat2 - $lat1;
    my $dphi = log( tan( $lat2/2 + pip4 ) / tan( $lat1/2 + pip4 ) );
    my $q = ( $dphi != 0 ) ? $dlat/$dphi : cos($lat1);# E-W line gives dPhi=0
    my $dlon = $d * sin( $brng ) / $q;

    if ( abs( $lat2 ) > pip2 ) {
        $lat2 = ( $lat2 > 0 ) ? pi-$lat2 : -(pi-$lat2);
    }
    my $lon2 = $lon1 + $dlon;
    $lon2 -= pi2 while( $lon2 > pi );
    $lon2 += pi2 while( $lon2 <= -(pi) );

    return { lat => $self->_precision( rad2deg($lat2), $precision ), lon => $self->_precision( rad2deg($lon2), $precision ) };
}


=head2 intersection

 $gc->intersection( $brng1, $point, $brng2, $precision );
 $gc->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180, -6 );

Returns the point of intersection of two paths defined by point and bearing

see http://williams.best.vwh.net/avform.htm#Intersection

=cut

method intersection( Num $brng1!, HashRef[Num] $point!, Num $brng2!, Int $precision? = -6 ) returns (HashRef[Num]) {
    my $lat1 = deg2rad( $self->get_lat() );
    my $lon1 = deg2rad( $self->get_lon() );
    my $lat2 = deg2rad( $point->{lat} );
    my $lon2 = deg2rad( $point->{lon} );
    my $brng13 = deg2rad( $brng1 );
    my $brng23 = deg2rad( $brng2 );
    my $dlat = $lat2 - $lat1;
    my $dlon = $lon2 - $lon1;

    my $dist12 = 2 * asin( sqrt( ( sin( $dlat/2 ) ** 2 ) + cos( $lat1 ) * cos( $lat2 ) * ( sin( $dlon/2 ) ** 2 ) ) );
    return undef if( $dist12 == 0 );

    #initial/final bearings between points
    my $brnga = acos( ( sin( $lat2 ) - sin( $lat1 ) * cos( $dist12 ) ) / ( sin( $dist12 ) * cos( $lat1 ) ) ) || 0;
    my $brngb = acos( ( sin( $lat1 ) - sin( $lat2 ) * cos( $dist12 ) ) / ( sin( $dist12 ) * cos( $lat2 ) ) ) || 0;

    my ( $brng12, $brng21 );
    if( sin( $dlon ) > 0 ) {
        $brng12 = $brnga;
        $brng21 = pi2 - $brngb;
    } else {
        $brng12 = pi2 - $brnga;
        $brng21 = $brngb;
    }

    my $alpha1 = $brng13 - $brng12;
    my $alpha2 = $brng21 - $brng23;
    $alpha1 -= pi2 while( $alpha1 > pi );
    $alpha1 += pi2 while( $alpha1 <= -(pi) );
    $alpha2 -= pi2 while( $alpha2 > pi );
    $alpha2 += pi2 while( $alpha2 <= -(pi) );

    return undef if( ( sin( $alpha1 ) == 0 ) and ( sin( $alpha2 ) == 0 ) ); #infinite intersections
    return undef if( sin( $alpha1 ) * sin( $alpha2 ) < 0 ); #ambiguous intersection

    my $alpha3 = acos( -cos( $alpha1 ) * cos( $alpha2 ) + sin( $alpha1 ) * sin( $alpha2 ) * cos( $dist12 ) );
    my $dist13 = atan2( sin( $dist12 ) * sin( $alpha1 ) * sin( $alpha2 ), cos( $alpha2 ) + cos( $alpha1 ) * cos( $alpha3 ) );
    my $lat3 = asin( sin( $lat1 ) * cos( $dist13 ) + cos( $lat1 ) * sin( $dist13 ) * cos( $brng13 ) );
    my $dlon13 = atan2( sin( $brng13 ) * sin( $dist13 ) * cos( $lat1 ), cos( $dist13 ) - sin( $lat1 ) * sin( $lat3 ) );
    my $lon3 = $lon1 + $dlon13;
    $lon3 -= pi2 while( $lon3 > pi );
    $lon3 += pi2 while( $lon3 <= -(pi) );

    return { lat => $self->_precision( rad2deg($lat3), $precision ), lon => $self->_precision( rad2deg($lon3), $precision ) };
}

sub _precision {
    my ( $self, $number, $precision ) = @_;

    die "Error: Private method called" unless (caller)[0]->isa( ref($self) );

    my $mbf = Math::BigFloat->new( $number );
    $mbf->precision( $precision );

    return $mbf->bstr();
}

sub _ib_precision {
    my ( $self, $brng, $precision ) = @_;

    die "Error: Private method called" unless (caller)[0]->isa( ref($self) );

    my $mbf;
    if( $brng =~ m/\./ ) {
        $mbf = Math::BigFloat->new( ( $brng + 360 ) % 360 .'.'. (split('\.', $brng ) )[1] );
    } else {
        $mbf = Math::BigFloat->new( ( $brng + 360 ) % 360 );
    }
    $mbf->precision( $precision );

    return $mbf->bstr();
}

sub _fb_precision {
    my ( $self, $brng, $precision ) = @_;

    die "Error: Private method called" unless (caller)[0]->isa( ref($self) );

    my $mbf;
    if( $brng =~ m/\./ ) {
        $mbf = Math::BigFloat->new( ( $brng + 180 ) % 360 .'.'. (split('\.', $brng ) )[1] );
    } else {
        $mbf = Math::BigFloat->new( ( $brng + 180 ) % 360 );
    }
    $mbf->precision( $precision );

    return $mbf->bstr();
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;
