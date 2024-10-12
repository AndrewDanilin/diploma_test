#pragma once

#include "point.h"
#include "vertex.h"
#include <cstddef>

namespace MyNamespace {

    enum EdgeType {
        E_TRENCH,
        E_DB,
        E_TRENCH_TO_DB,
    };

    struct edge {
        size_t idx;
        
        vertex u;
        vertex v;

        EdgeType edgeType;

        double len;
        double cost;

        bool cross;
        bool inner;

        bool operator==(const edge &other) const;
        bool intersectsStrictly(const edge &other) const;
        point getIntersectionPoint(const edge &other) const;
        point getMidPoint() const;
        double find_angle(const edge &other) const;
        point getPointAtDistance(const double &d) const;
    };

}