#include <cstddef>
#include <cmath>
#include "../headers/my_utils.h"

namespace MyNamespace {

    bool point::operator==(const point &other) const {
        return (std::fabs(x - other.x) < EPSILON && std::fabs(y - other.y) < EPSILON);
    }

    bool point::operator!=(const point &other) const {
        return !((std::fabs(x - other.x) < EPSILON && std::fabs(y - other.y) < EPSILON));
    }

    double point::distanceTo(const point &other) const {
        return std::sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
    }

}