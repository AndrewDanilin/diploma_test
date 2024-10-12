#pragma once

#include <vector>
#include "my_utils.h"
#include "point.h"
#include "edge.h"


namespace MyNamespace {
    struct polygon {
        std::vector<point> points;

        bool containtPointGeometryStrictly(const point &p, const double epsilon = EPSILON) const;
        bool intersectsSegment(const point &p1, const point &p2, const edge &e) const;
        
    };

}
