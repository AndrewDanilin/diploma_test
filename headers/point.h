#pragma once

#include "point.h"

namespace MyNamespace {

    struct point {
        double x;
        double y;

        bool operator==(const point &other) const;
        bool operator!=(const point &other) const;
        double distanceTo(const point &other) const;

    };

}