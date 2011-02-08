# Copyrights 2011 by Sorin Pop.
# For other contributors see ChangeLog.
# See the manual pages for details on the licensing terms.

package Geo::Calc;

use vars '$VERSION';
$VERSION = '0.01';

use Moose;
use MooseX::FollowPBP;
use MooseX::Method::Signatures;

use Math::Trig qw(:pi asin acos tan deg2rad rad2deg);
use Math::BigFloat;

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

=head2 Geo::Calc - simple calculations for geo coordinates
    my $gc = Geo::Calc->new( lat => $latitude, lon => $longitude );
    my $gc = Geo::Calc->new( lat => $latitude, lon => $longitude, radius => 6371 );
=cut

=head2 distance_to - Returns the distance from this point to the supplied point,
in km (using Haversine formula)

from: Haversine formula - R. W. Sinnott, "Virtues of the Haversine",
        Sky and Telescope, vol 68, no 2, 1984
    my $distance = $gc->distance_to( $point, $precision );
    my $distance = $gc->distance_to( { lat => 40.422371, lon => -3.704298 }, -6 );
=cut
sub distance_to {
    my ( $self, $point, $precision ) = @_;

    $precision ||= -6;
    my $lat1 = deg2rad( $self->get_lat() );
    my $lon1 = deg2rad( $self->get_lon() );
    my $lat2 = deg2rad( $point->{lat} );
    my $lon2 = deg2rad( $point->{lon} );

    my $t = sin( ($lat2 - $lat1)/2 ) ** 2 + ( cos( $lat1 ) ** 2 ) * ( sin( ( $lon2 - $lon1 )/2 ) ** 2 );
    my $d = $self->get_radius * ( 2 * atan2( sqrt($t), sqrt(1-$t) ) );

    return $self->_precision( $d, $precision );
}

=head2 bearing_to - Returns the (initial) bearing from this point to the supplied
point, in degrees
    see http://williams.best.vwh.net/avform.htm#Crs
    my $brng = $gc->bearing_to( $point, $precision );
    my $brng = $gc->bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
=cut
sub bearing_to {
    my ( $self, $point, $precision ) = @_;

    $precision ||= -6;
    my $lat1 = deg2rad( $self->get_lat() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlon = deg2rad( $self->get_lon() - $point->{lon} );

    my $brng = atan2( sin( $dlon ) * cos( $lat2 ), ( cos( $lat1 ) * sin( $lat2 ) ) - ( sin( $lat1 ) * cos( $lat2 ) * cos( $dlon ) ) );

    return $self->_ib_precision( rad2deg( $brng ), $precision );
}

=head2 final_bearing_to - Returns final bearing arriving at supplied destination
point from this point; the final bearing will differ from the initial bearing
by varying degrees according to distance and latitude
    my $f_brng = $gc->final_bearing_to( $point, $precision );
    my $f_brng = $gc->final_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
=cut
sub final_bearing_to {
    my ( $self, $point, $precision ) = @_;

    $precision ||= -6;
    my $lat1 = deg2rad( $self->get_lat() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlon = deg2rad( $self->get_lon() - $point->{lon} );

    my $brng = atan2( sin( $dlon ) * cos( $lat2 ), ( cos( $lat1 ) * sin( $lat2 ) ) - ( sin( $lat1 ) * cos( $lat2 ) * cos( $dlon ) ) );

    return $self->_fb_precision( rad2deg( $brng ), $precision );
}

=head2 midpoint_to - Returns the midpoint between this point and the supplied
point.
    see http://mathforum.org/library/drmath/view/51822.html for derivation
    my $midpoint = $gc->midpoint_to( $point, $precision );
    my $midpoint = $gc->midpoint_to( { lat => 40.422371, lon => -3.704298 }, -6 );
=cut
sub midpoint_to {
    my ( $self, $point, $precision ) = @_;

    $precision ||= -6;
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

=head2 destination_point - Returns the destination point from this point having
travelled the given distance (in km) on the given initial bearing (bearing may
vary before destination is reached)
    see http://williams.best.vwh.net/avform.htm#LL
    my $destination = $gc->destination_point( $bearing, $distance, $precision );
    my $destination = $gc->destination_point( 90, 1, -6 ); # distance in km
=cut
sub destination_point {
    my ( $self, $brng, $dist, $precision ) = @_;

    $precision ||= -6;
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

=head2 boundry_box - Returns the boundry box min/max having the initial point
defined as the center of the boundry box, given the widht and height
    my $bbox = $gc->boundry_box( $width, $height, $precision ); # in km
    my $bbox = $gc->boundry_box( 3, 4, -6 ); # in km
=cut
sub boundry_box {
    my ( $self, $width, $height, $precision ) = @_;

    $precision ||= -6;
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

=head2 rhumb_distance_to - Returns the distance from this point to the supplied
point, in km, travelling along a rhumb line
    see http://williams.best.vwh.net/avform.htm#Rhumb
    my $r_distance = $gc->rhumb_distance_to( $point, $precision );
    my $r_distance = $gc->rhumb_distance_to( { lat => 40.422371, lon => -3.704298 }, -6 );
=cut
sub rhumb_distance_to {
    my ( $self, $point, $precision ) = @_;

    $precision ||= -6;
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

=head2 rhumb_bearing_to - Returns the bearing from this point to the supplied
point along a rhumb line, in degrees
    my $r_brng = $gc->rhumb_bearing_to( $point, $precision );
    my $r_brng = $gc->rhumb_bearing_to( { lat => 40.422371, lon => -3.704298 }, -6 );
=cut
sub rhumb_bearing_to {
    my ( $self, $point, $precision ) = @_;

    $precision ||= -6;
    my $lat1 = deg2rad( $self->get_lat() );
    my $lat2 = deg2rad( $point->{lat} );
    my $dlon = deg2rad( $point->{lon} - $self->get_lon() );


    my $dphi = log( tan( $lat2/2 + pip4 ) / tan( $lat1/2 + pip4 ) );
    if( abs( $dlon ) > pi ) {
        $dlon = ( $dlon > 0 ) ? -(pi2-$dlon) : (pi2+$dlon);
    }

    return $self->_ib_precision( rad2deg( atan2( $dlon, $dphi ) ), $precision );
}

=head2 rhumb_destination_point - Returns the destination point from this point
having travelled the given distance (in km) on the given bearing along a rhumb
line
    my $r_destination = $gc->rhumb_destination_point( $brng, $distance, $precision );
    my $r_destination = $gc->rhumb_destination_point( 30, 1, -6 );
=cut
sub rhumb_destination_point {
    my ( $self, $brng, $dist, $precision ) = @_;

    $precision ||= -6;
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


=head2 intersection - Returns the point of intersection of two paths defined
by point and bearing
    see http://williams.best.vwh.net/avform.htm#Intersection
    my $point = $gc->intersection( $brng1, $point, $brng2, $precision );
    my $point = $gc->intersection( 90, { lat => 40.422371, lon => -3.704298 }, 180, -6 );
=cut
sub intersection {
    my ( $self, $brng1, $point, $brng2, $precision ) = @_;

    $precision ||= -6;
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

1;
