//
// Created by ricardo on 25-04-2024.
//

#ifndef FEUP_DA_2324_PROJ2_POINTS_H
#define FEUP_DA_2324_PROJ2_POINTS_H
#include <cmath>

class points {
private:
    double lat;
    double longit;
public:
    points(double latitude, double longitude);

    double getLatitude() const;

    void setLatitude(double latitude);

    double getLongitude() const;

    void setLongitude(double longitude);

    double haversine(points c) const;

    bool operator==(const points &rhs) const;

    bool operator!=(const points &rhs) const;
};


#endif //FEUP_DA_2324_PROJ2_POINTS_H
