#include "../headers/polygon.h"
#include "../headers/edge.h"
#include "../headers/point.h"

namespace MyNamespace {
    bool polygon::containtPointGeometryStrictly(const point &p, const double epsilon) const {
        int n = points.size();
        if (n < 3) return false; 

        point extreme = {1e9, p.y};
        int count = 0;

        for (int i = 0; i < n; i++) {
            point u = points[i];
            point v = points[(i + 1) % n]; 

            if (onSegment(p, u, v, epsilon)) {
                return false; 
            }

            if (isGreater(u.y , v.y, epsilon)) std::swap(u, v); 

            if (isGreater(p.y, u.y, epsilon) && isLessOrEqual(p.y, v.y, epsilon)) {
                double x_intersection = u.x + (p.y - u.y) * (v.x - u.x) / (v.y - u.y);
                
                if (isGreater(x_intersection, p.x, epsilon)) {
                    count++;
                }
            }
        }

        return count % 2 == 1;
    }

    bool polygon::intersectsSegment(const point &p1, const point &p2, const edge &e) const {
        for (size_t i = 0; i < points.size() - 1; ++i) {
            if (edge{0, vertex{points[i]}, vertex{points[i + 1]}}.intersectsStrictly(e)) {
                point intersection = edge{0, vertex{points[i]}, vertex{points[i + 1]}}.getIntersectionPoint(e);
                if (intersection == p1 || intersection == p2) {
                    continue;
                }
                return true;
            }
        }
        return false;
    }
}