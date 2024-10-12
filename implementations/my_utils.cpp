#include <cmath>
#include <vector>
#include <ctime>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <string>
#include "../headers/my_utils.h"
#include "../headers/polygon.h"
#include "../headers/point.h"
#include "../headers/edge.h"
#include "../headers/vertex.h"


namespace MyNamespace {

    bool isGreater(double a, double b, const double epsilon) {
        return (a > b + epsilon);
    }


    bool isLess(double a, double b, const double epsilon) {
        return (a > b + epsilon);
    }

    int orientation(point p, point q, point r) {
        int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
        if (val < EPSILON) return 0;
        return (val > EPSILON) ? 1 : 2;
    }

    bool isGreaterOrEqual(double a, double b, const double epsilon) {
        return (a > b + epsilon) || (std::abs(a - b) < epsilon);
    }

    bool isLessOrEqual(double a, double b, const double epsilon) {
        return (a < b - epsilon) || (std::abs(a - b) < epsilon);
    }

    bool intersects_1(double a, double b, double c, double d) {
        if (a > b)  std::swap (a, b);
        if (c > d)  std::swap (c, d);
        return std::max(a,c) <= std::min(b,d);
    }

    
    double area(const point &a, const point &b, const point &c) {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    bool onSegment(const point& p, const point& u, const point& v, const double epsilon) {
        if (p.x <= std::max(u.x, v.x) + epsilon && p.x >= std::min(u.x, v.x) - epsilon &&
            p.y <= std::max(u.y, v.y) + epsilon && p.y >= std::min(u.y, v.y) - epsilon) {
            return std::fabs((v.y - u.y) * (p.x - u.x) - (p.y - u.y) * (v.x - u.x)) < epsilon;
        }
        return false;
    }

    std::vector<point> generateArc(const edge &e, double h, int num_segments) {
        point A = e.u.point;
        point B = e.v.point;

        double mid_x = (A.x + B.x) / 2;
        double mid_y = (A.y + B.y) / 2;

        double L = sqrt((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y));

        double perp_x = -(B.y - A.y);
        double perp_y = (B.x - A.x);

        double perp_length = sqrt(perp_x * perp_x + perp_y * perp_y);
        perp_x /= perp_length;
        perp_y /= perp_length;

        double C_x = mid_x + h * perp_x;
        double C_y = mid_y + h * perp_y;

        double R = sqrt((C_x - A.x) * (C_x - A.x) + (C_y - A.y) * (C_y - A.y));

        double start_angle = atan2(A.y - C_y, A.x - C_x);
        double end_angle = atan2(B.y - C_y, B.x - C_x);

        if (end_angle < start_angle) {
            end_angle += 2 * M_PI;
        }

        std::vector<point> arc_points;
        for (int i = 0; i <= num_segments; ++i) {
            double theta = start_angle + (end_angle - start_angle) * (static_cast<double>(i) / num_segments);
            point P;
            P.x = C_x + R * cos(theta);
            P.y = C_y + R * sin(theta);
            arc_points.push_back(P);
        }

        return arc_points;
    }

    double lenToCost(const double &len, const double &cost_per_meter) {
        return len * cost_per_meter;
    }

    bool connectedEdges(const edge &e1, const edge &e2) {
        return e1.u.idx == e2.u.idx || e1.u.idx == e2.v.idx || e1.v.idx == e2.u.idx || e1.v.idx == e2.v.idx;
    }

    // Определяем какие компоненты связности внешние, какие внутренние
    // Нужно, чтобы правильно отрисовать дуги рёбер
    void classifyPolygons(std::vector<std::vector<edge>> &components) {
        std::vector<polygon> polygons;

        for (const auto &comp : components) {
            polygons.push_back(polygon{});
            for (const auto &e : comp) {
                polygons[polygons.size() - 1].points.push_back(e.u.point);
            }
            polygons[polygons.size() - 1].points.push_back(polygons[polygons.size() - 1].points[0]);
        }

        int n = polygons.size();
        std::vector<int> classification(n, 0); 
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j && polygons[j].containtPointGeometryStrictly(polygons[i].points[0])) {
                    classification[i] = -classification[j]; 
                }
            }
            if (classification[i] == 0) {
                classification[i] = 1; 
            }
        }
        
        for (size_t i = 0; i < n; ++i) {
            bool flag;
            if (classification[i] == 1) {
                flag = false;
            } else {
                flag = true;
            }
            for (auto &e: components[i]) {
                    e.inner = flag;
            }
        }
    }

    void logMessage(std::string message) {
        auto now = std::chrono::system_clock::now();
        std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
        std::tm* localTime = std::localtime(&currentTime);

        std::cout << "[" << std::put_time(localTime, "%H:%M:%S") << "] " << message << std::endl; 
    }
}