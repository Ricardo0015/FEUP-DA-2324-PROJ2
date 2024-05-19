#include "points.h"

/**
 * @brief Constructor for the points class.
 *
 * @param latitude The latitude of the point.
 * @param longitude The longitude of the point.
 *
 * @complexity O(1)
 */
points::points(double latitude, double longitude) : lat(latitude), longit(longitude) {}

/**
 * @brief Calculates the Haversine distance between this point and another point.
 *
 * @param point The other point to calculate the distance to.
 * @return The Haversine distance in kilometers.
 *
 * @complexity O(1)
 */
double points::haversine(points point) const {
    double lat1 = this->lat, lat2 = point.lat;
    double longit1 = this->longit, longit2 = point.longit;

    if ((lat1 == 0 && longit1 == 0) || (lat2 == 0 && longit2 == 0)) return -1;

    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (longit2 - longit1) * M_PI / 180.0;

    // convert to radians
    lat1 = lat1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;
    double formula = pow(sin(dLat / 2), 2) +
                     pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
    double radius = 6371;
    double c = 2 * asin(sqrt(formula));
    return radius * c;
}

/**
 * @brief Gets the latitude of the point.
 *
 * @return The latitude of the point.
 *
 * @complexity O(1)
 */
double points::getLatitude() const {
    return lat;
}

/**
 * @brief Sets the latitude of the point.
 *
 * @param latitude The new latitude of the point.
 *
 * @complexity O(1)
 */
void points::setLatitude(double latitude) {
    points::lat = latitude;
}

/**
 * @brief Gets the longitude of the point.
 *
 * @return The longitude of the point.
 *
 * @complexity O(1)
 */
double points::getLongitude() const {
    return longit;
}

/**
 * @brief Sets the longitude of the point.
 *
 * @param longitude The new longitude of the point.
 *
 * @complexity O(1)
 */
void points::setLongitude(double longitude) {
    points::longit = longitude;
}

/**
 * @brief Equality operator to compare two points.
 *
 * @param rhs The right-hand side point to compare.
 * @return true if the points are equal, false otherwise.
 *
 * @complexity O(1)
 */
bool points::operator==(const points &rhs) const {
    return lat == rhs.lat && longit == rhs.longit;
}

/**
 * @brief Inequality operator to compare two points.
 *
 * @param rhs The right-hand side point to compare.
 * @return true if the points are not equal, false otherwise.
 *
 * @complexity O(1)
 */
bool points::operator!=(const points &rhs) const {
    return !(rhs == *this);
}