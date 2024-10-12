#pragma once

#include "point.h"
#include "edge.h"
#include <string>
#include <cstddef>
#include <vector>

#define _USE_MATH_DEFINES

#define EPSILON 1e-5

namespace MyNamespace {
    bool isGreater(double a, double b, const double epsilon = EPSILON);
    bool isLess(double a, double b, const double epsilon = EPSILON);
    int orientation(point p, point q, point r);
    bool isGreaterOrEqual(double a, double b, const double epsilon = EPSILON);
    bool isLessOrEqual(double a, double b, const double epsilon = EPSILON);
    bool intersects_1(double a, double b, double c, double d);
    double area(const point &a, const point &b, const point &c);
    bool onSegment(const point& p, const point& u, const point& v, const double epsilon = EPSILON);
    std::vector<point> generateArc(const edge &e, double h, int num_segments);
    double lenToCost(const double &len, const double &cost_per_meter);
    bool connectedEdges(const edge &e1, const edge &e2);
    void classifyPolygons(std::vector<std::vector<edge>> &components);
    void logMessage(std::string message);
}
