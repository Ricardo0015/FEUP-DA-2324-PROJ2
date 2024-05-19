#include "points.h"

points::points(double latitude, double longitude) : lat(latitude), longit(longitude) {}

double points::haversine(points point) const {
    double lat1 = this->lat, lat2 = point.lat;
    double longit1 = this->longit, longit2 = point.longit;

    if ((lat1 == 0 && longit == 0) || (lat2 == 0 && longit == 0)) return -1;

    double dLat = (lat2 - lat1) *
                  M_PI / 180.0;
    double dLon = (longit2 - longit1) *
                  M_PI / 180.0;

    // convert to radians
    lat1 = (lat1) * M_PI / 180.0;
    lat2 = (lat2) * M_PI / 180.0;
    double formula = pow(sin(dLat / 2), 2) +
               pow(sin(dLon / 2), 2) *
               cos(lat1) * cos(lat2);
    double radious = 6371;
    double c = 2 * asin(sqrt(formula));
    return radious * c;
}

double points::getLatitude() const {
    return lat;
}

void points::setLatitude(double latitude) {
    points::lat = latitude;
}

double points::getLongitude() const {
    return longit;
}

void points::setLongitude(double longitude) {
    points::longit = longitude;
}

bool points::operator==(const points &rhs) const {
    return lat == rhs.lat &&
           longit == rhs.longit;
}

bool points::operator!=(const points &rhs) const {
    return !(rhs == *this);
}