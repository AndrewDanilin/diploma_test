#pragma once

#include <cstddef>
#include "point.h"

namespace MyNamespace {

    enum VertexType {
        V_TRENCH,
        V_DB
    };

    struct vertex {
        point point;
        size_t idx;
        VertexType vertexType;

        bool operator==(const vertex &other) const;
    };
    
}