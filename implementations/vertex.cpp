#include "../headers/vertex.h"

namespace MyNamespace {
    bool vertex::operator==(const vertex &other) const {
        return point == other.point && idx == other.idx && vertexType == other.vertexType;
    }
}