#include <cstddef>
#include <cmath>
#include "../headers/edge.h"
#include "../headers/my_utils.h"

#define _USE_MATH_DEFINES


namespace MyNamespace {
    bool edge::operator==(const edge &other) const {
        return (u == other.u && v == other.v);
    }

    bool edge::intersectsStrictly(const edge &other) const {
        if (u.idx == other.u.idx || v.idx == other.v.idx || u.idx == other.v.idx || v.idx == other.u.idx) {
                return false;
        }


        int o1 = orientation(u.point, v.point, other.u.point);
        int o2 = orientation(u.point, v.point, other.v.point);
        int o3 = orientation(other.u.point, other.v.point, u.point);
        int o4 = orientation(other.u.point, other.v.point, v.point);

        if (o1 != o2 && o3 != o4) {
            return true;
        }

        return false;
    }

    point edge::getIntersectionPoint(const edge &other) const {
        point a = u.point;
        point b = v.point;
        point c = other.u.point;
        point d = other.v.point;

        double A1 = b.y - a.y;
        double B1 = a.x - b.x;
        double C1 = A1 * a.x + B1 * a.y;

        double A2 = d.y - c.y;
        double B2 = c.x - d.x;
        double C2 = A2 * c.x + B2 * c.y;

        double det = A1 * B2 - A2 * B1;

        return point{(B2 * C1 - B1 * C2) / det, (A1 * C2 - A2 * C1) / det};
    }

    point edge::getMidPoint() const {
        point midpoint;
        midpoint.x = (u.point.x + v.point.x) / 2.0;
        midpoint.y = (u.point.y + v.point.y) / 2.0;
        return midpoint;
    }

    double edge::find_angle(const edge &other) const {
        double Ax = v.point.x - u.point.x;
        double Ay = v.point.y - u.point.y;

        double Bx = other.v.point.x - other.u.point.x;
        double By = other.v.point.y - other.u.point.y;

        double dot_poduct = Ax * Bx + Ay * By;

        double magnitude_A = std::sqrt(Ax * Ax + Ay * Ay);
        double magnitude_B = std::sqrt(Bx * Bx + By * By);

        double cos = dot_poduct / (magnitude_A * magnitude_B);

        cos = std::max(-1.0, std::min(1.0, cos));

        double angle_radius = std::acos(cos);

        double angle_degrees = angle_radius * 180.0 / M_PI;

        return angle_degrees;
    }

    point edge::getPointAtDistance(const double &d) const {
        double L = std::sqrt((v.point.x - u.point.x) * (v.point.x - u.point.x) + (v.point.y - u.point.y) * (v.point.y - u.point.y));

        double t = d / L;

        point P;
        P.x = u.point.x + t * (v.point.x - u.point.x);
        P.y = u.point.y + t * (v.point.y - u.point.y);

        return P;
    }
}