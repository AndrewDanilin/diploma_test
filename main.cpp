#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <optional>
#include <queue>
#include <utility>
#include <cmath>
#include "third-party/nlohmann/json.hpp"
#include "headers/graph.h"
#include "headers/my_utils.h"
#include "headers/polygon.h"
#include "headers/point.h"

using namespace MyNamespace;

#define BETA 10.0

struct parameters {
    int cost_trench_per_meter;
    int cost_db_per_meter;
    int cost_trench_to_db; 
    int min_length_db;
    int max_length_db;
    int alpha;
};

int readConfig(const std::string path, parameters &params) {
    logMessage("parsing config");

    std::ifstream in(path);
    nlohmann::json data = nlohmann::json::parse(in);

    if (data["cost_trench_per_meter"] == NULL) {
        std::cout << "Cost trench per meter missed in config file" << std::endl;
        return 1;
    }
    params.cost_trench_per_meter = data["cost_trench_per_meter"];

    if (data["cost_db_per_meter"] == NULL) {
        std::cout << "Cost db per meter missed in config file" << std::endl;
        return 1;
    }
    params.cost_db_per_meter = data["cost_db_per_meter"];

    if (data["cost_trench_to_db"] == NULL) {
        std::cout << "Cost trench to db missed in config file" << std::endl;
        return 1;
    } 
    params.cost_trench_to_db = data["cost_trench_to_db"];

    if (data["min_length_db"] == NULL) {
        std::cout << "Min length db missed in config file" << std::endl;
        return 1;
    }
    params.min_length_db = data["min_length_db"];

    if (data["max_length_db"] == NULL) {
        std::cout << "Max length db missed in config file" << std::endl;
        return 1;
    }
    params.max_length_db = data["max_length_db"];

    if (data["alpha"] == NULL) {
        std::cout << "Alpha missed in config file" << std::endl;
        return 1;
    }
    params.alpha = data["alpha"];

    return 0;
}

std::vector<point> parseFeature(const nlohmann::json &feature) {
    if (feature["geometry"]["type"] != "Polygon") {
        std::cout << "Only polygon implemented" << std::endl;
    }

    nlohmann::json geom = feature["geometry"];
    std::vector<point> points;

    for (const auto &coordinate : geom["coordinates"]) {
        for (const auto &coord : coordinate) {
            points.push_back(point{coord[0], coord[1]});
        }
    }

    return points;
}

std::vector<polygon> readGeoJSON(const std::string path) {
    logMessage("parsing input .geojson file");

    std::ifstream in(path);
    nlohmann::json data = nlohmann::json::parse(in);
    std::vector<polygon> polygons;

    for (const auto &feature : data["features"]) {
        std::vector<point> polygon_coords = parseFeature(feature);
        polygons.push_back(polygon{polygon_coords});
    }

    return polygons;
}

nlohmann::json getJsonFromGraph(const graph &g, const bool useArc = false) {
    nlohmann::json resultJson = nlohmann::json::object();
    nlohmann::json trenchFeatureCollection = nlohmann::json::object();
    trenchFeatureCollection["type"] = "FeatureCollection";
    trenchFeatureCollection["name"] = "trench_graph";  
    trenchFeatureCollection["crs"] = nlohmann::json::object();
    trenchFeatureCollection["crs"]["type"] = "name";
    trenchFeatureCollection["crs"]["properties"] = nlohmann::json::object();
    trenchFeatureCollection["crs"]["properties"]["name"] = "urn:ogc:def:crs:EPSG::3857";
    trenchFeatureCollection["features"] = nlohmann::json::array();

    nlohmann::json dbFeatureCollection = nlohmann::json::object();
    dbFeatureCollection["type"] = "FeatureCollection";
    dbFeatureCollection["name"] = "db_graph";  
    dbFeatureCollection["crs"] = nlohmann::json::object();
    dbFeatureCollection["crs"]["type"] = "name";
    dbFeatureCollection["crs"]["properties"] = nlohmann::json::object();
    dbFeatureCollection["crs"]["properties"]["name"] = "urn:ogc:def:crs:EPSG::3857";
    dbFeatureCollection["features"] = nlohmann::json::array();

    nlohmann::json trenchToDbFeatureCollection = nlohmann::json::object();
    trenchToDbFeatureCollection["type"] = "FeatureCollection";
    trenchToDbFeatureCollection["name"] = "trench_to_db_graph";  
    trenchToDbFeatureCollection["crs"] = nlohmann::json::object();
    trenchToDbFeatureCollection["crs"]["type"] = "name";
    trenchToDbFeatureCollection["crs"]["properties"] = nlohmann::json::object();
    trenchToDbFeatureCollection["crs"]["properties"]["name"] = "urn:ogc:def:crs:EPSG::3857";
    trenchToDbFeatureCollection["features"] = nlohmann::json::array();

    for (const auto &e : g.edges) {
        nlohmann::json new_feat_obj = nlohmann::json::object();
        nlohmann::json new_geom_obj = nlohmann::json::object();

        if (useArc && e.edgeType == E_DB && !e.cross) {
            edge cur_edge = e;
            if (!e.inner) {
                cur_edge = edge{0, e.v, e.u};
            } 
            std::vector<point> arc = generateArc(cur_edge, 5, 20);

            new_geom_obj["coordinates"] = nlohmann::json::array();
            for (const auto &p: arc) {
                new_geom_obj["coordinates"].push_back({p.x, p.y});
            }
        } else {
            new_geom_obj["coordinates"] = nlohmann::json::array({{e.u.point.x, e.u.point.y}, {e.v.point.x, e.v.point.y}});
        }

        new_geom_obj["type"] = "LineString";

        new_feat_obj["geometry"] = new_geom_obj;
        new_feat_obj["properties"] = nlohmann::json::object();
        new_feat_obj["type"] = "Feature";

        nlohmann::json *resJson;
        if (e.edgeType == E_DB) {
            resJson = &dbFeatureCollection;
        } else if (e.edgeType == E_TRENCH) {
            resJson = &trenchFeatureCollection;
        } else {
            resJson = &trenchToDbFeatureCollection;
        }

        (*resJson)["features"].push_back(new_feat_obj); 
    }

    for (const auto &v : g.vs) {
        nlohmann::json new_feat_obj2 = nlohmann::json::object();
        nlohmann::json new_geom_obj2 = nlohmann::json::object();
        new_geom_obj2["coordinates"] = {v.point.x, v.point.y};
        new_geom_obj2["type"] = "Point";

        new_feat_obj2["geometry"] = new_geom_obj2;
        new_feat_obj2["properties"] = nlohmann::json::object();
        new_feat_obj2["properties"]["number"] = v.idx;
        new_feat_obj2["type"] = "Feature";

        nlohmann::json *resJson;
        if (v.vertexType == V_DB) {
            resJson = &dbFeatureCollection;
        } else {
            resJson = &trenchFeatureCollection;
        }

        (*resJson)["features"].push_back(new_feat_obj2);
    }

    resultJson["trench_graph"] = trenchFeatureCollection;
    resultJson["db_graph"] = dbFeatureCollection;
    resultJson["trench_to_db_graph"] = trenchToDbFeatureCollection;
    return resultJson;
}

void writeToJSON(const nlohmann::json &json, const std::string &path) {
    logMessage("writing out .geojson file");
    std::ofstream out(path);
    out << std::setw(4) << json << std::endl;
}


// Строит бызовый граф, в котором вершины - это вершины пилогонов из .geojson
void buildBaseTrenchGraph(graph &g, const polygon &poly, size_t &next_idx, const parameters &params) {

    for (size_t i = 0; i < poly.points.size() - 1; ++i) {
        const point p1 = poly.points[i];
        const point p2 = poly.points[i + 1];
        const double dist = p1.distanceTo(p2);
        const double cost = lenToCost(dist, params.cost_trench_per_meter);

        size_t new_idx = (i == poly.points.size() - 2) ? next_idx - poly.points.size() + 2 : next_idx + 1;

        vertex u = vertex{p1, next_idx, VertexType(V_TRENCH)};
        vertex v = vertex{p2, new_idx, VertexType(V_TRENCH)};
        edge edge_uv = edge{new_idx, u, v, EdgeType(E_TRENCH), dist, cost};
        edge edge_vu = edge{next_idx, v, u, EdgeType(E_TRENCH), dist, cost};

        g.vs.push_back(u);
        g.edges.push_back(edge_uv);

        
        if (g.map.find(next_idx) == g.map.end()) {
            g.map[next_idx] = {};
        }

        g.map[next_idx].push_back(edge_uv);

        if (g.map.find(new_idx) == g.map.end()) {
            g.map[new_idx] = {};
        }

        g.map[new_idx].push_back(edge_vu);

        next_idx++;
    }
}

// Удаляет из графа рёбра, которые лежат на пересечении полигонов
void deleteInvalidEdges(graph &g, std::vector<polygon> polygons) {
    std::vector<vertex> to_remove;

    for (auto &poly: polygons) {
        for (const auto &v: g.vs) {
            if (poly.containtPointGeometryStrictly(v.point)) {
                to_remove.push_back(v);
            }
        }
    }
    g.eraseVertices(to_remove);

    std::vector<edge> to_remove_e;
    for (const auto &e: g.edges) {
        for (const auto &v: g.vs) {
            if (v.idx == e.u.idx || v.idx == e.v.idx) {
                continue;
            }

            if (onSegment(v.point, e.u.point, e.v.point)) {
                to_remove_e.push_back(e);
            }
        }
    }
    g.eraseEdges(to_remove_e);

}

// Обрабатывает пересечения, добавляя новые рёбра и вершины, и удаляя лишние
void handleIntersectingPolygons(graph &g, const std::vector<polygon> polygons, size_t &next_idx, const parameters &params) {
    logMessage("processing intersecting polygons");

    std::vector<edge> to_remove;
    std::vector<edge> to_add;

    for (size_t i = 0; i < g.edges.size() - 1; ++i) {
        for (size_t j = i + 1; j < g.edges.size(); ++j) {
            edge e1 = g.edges[i];
            edge e2 = g.edges[j];
            if (e1.intersectsStrictly(e2)) {
                point intersection_p = e1.getIntersectionPoint(e2);
                to_remove.push_back(e1);
                to_remove.push_back(e2);

                vertex new_v = vertex{intersection_p, next_idx, VertexType(V_TRENCH)};
                g.addVertex(new_v);
                next_idx++;

                double dist1 = e1.u.point.distanceTo(intersection_p);
                double cost1 = lenToCost(dist1, params.cost_trench_per_meter);
                to_add.push_back(edge{new_v.idx, e1.u, new_v, EdgeType(E_TRENCH), dist1, cost1});
                
                double dist2 = e1.v.point.distanceTo(intersection_p);
                double cost2 = lenToCost(dist2, params.cost_trench_per_meter);
                to_add.push_back(edge{e1.v.idx, new_v, e1.v, EdgeType(E_TRENCH), dist2, cost2});

                double dist3 = e2.u.point.distanceTo(intersection_p);
                double cost3 = lenToCost(dist3, params.cost_trench_per_meter);
                to_add.push_back(edge{new_v.idx, e2.u, new_v, EdgeType(E_TRENCH), dist3, cost3});

                double dist4 = e2.v.point.distanceTo(intersection_p);
                double cost4 = lenToCost(dist4, params.cost_trench_per_meter);
                to_add.push_back(edge{e2.v.idx, new_v, e2.v, EdgeType(E_TRENCH), dist4, cost4});
            }
        }
    }

    g.eraseEdges(to_remove);
    g.addEdges(to_add);
    deleteInvalidEdges(g, polygons);
}

// Проверяет угол между рёбрами на условие 90 +- альфа
// isCrossing - булевый флаг, которые говорит, что рёбра находятся друг напротив друга
bool isValidCrossing(const edge &e1, const vertex &v1, const vertex &v2, const parameters &params, const bool &isCrossing) {
    double angle = e1.find_angle(edge{0, v1, v2});

    if (isCrossing && angle >= 90.0 - params.alpha && angle <= 90.0 + params.alpha) {
        return true;
    } else if (!isCrossing && angle < BETA) {
        return true;
    }

    return false;
}


// Ищем полигон, который содержит точку
bool findPolygonContainsPoint(const std::vector<polygon> &polygons, const point &p) {
    for (const auto &poly: polygons) {
        if (poly.containtPointGeometryStrictly(p)) {
            return true;
        }
    }
    return false;
}

// Проверяем, что ребро либо строго лежит вне полигонов, либо пересекает какой-то из полигонов
bool isValidEdge(const std::vector<polygon> &polygons, const edge &e) {

    if (!findPolygonContainsPoint(polygons, e.getMidPoint())) {
        return false;
    }

    for (const auto &poly: polygons) {
        if (poly.intersectsSegment(e.u.point, e.v.point, e)) {
            return false;
        }
    }

    return true;
}

void checkAndAdd(
    const edge &first,
    const vertex &v1,
    const vertex &v2,
    graph &g,
    std::vector<edge> &e_to_add,
    std::vector<vertex> &v_to_add,
    const parameters &params,
    const bool &isCrossing,
    const std::vector<polygon> polygons,
    size_t &next_idx) {

    if (isValidCrossing(first, v1, v2, params, isCrossing)) {
        double dist = v1.point.distanceTo(v2.point);
        if (isGreaterOrEqual(dist, params.min_length_db) && isLessOrEqual(dist, params.max_length_db) &&
            isValidEdge(polygons, edge{0, v1, v2})) {
                double cost = lenToCost(dist, params.cost_db_per_meter);
                if (g.mapToDb.find(v1.idx) == g.mapToDb.end()) {
                    vertex new_db_v1 = vertex{v1.point, next_idx, VertexType(V_DB)};
                    edge new_trench_to_db_e = edge{
                        v1.idx,
                        new_db_v1,
                        v1,
                        EdgeType(E_TRENCH_TO_DB),
                        0, 
                        (double) params.cost_trench_to_db
                    };
                    e_to_add.push_back(new_trench_to_db_e);
                    v_to_add.push_back(new_db_v1);
                    next_idx++;
                    g.mapToDb[v1.idx] = new_db_v1;
                }

                if (g.mapToDb.find(v2.idx) == g.mapToDb.end()) {
                    vertex new_db_v2 = vertex{v2.point, next_idx, VertexType(V_DB)};
                    edge new_trench_to_db_e = edge{
                        v2.idx,
                        new_db_v2,
                        v2,
                        EdgeType(E_TRENCH_TO_DB),
                        0, 
                        (double) params.cost_trench_to_db
                    };
                    e_to_add.push_back(new_trench_to_db_e);
                    v_to_add.push_back(new_db_v2);
                    next_idx++;
                    g.mapToDb[v2.idx] = new_db_v2;
                }

                vertex u1 = g.mapToDb[v1.idx];
                vertex u2 = g.mapToDb[v2.idx];

                edge new_e = edge{u2.idx, u1, u2, EdgeType(E_DB), dist, cost, true};
                e_to_add.push_back(new_e);
        }
    }
}

// Добавляем в граф ГНБ рёбра, которые либо лежат внутри полигона, соединяя близкие вершины, либо которые лежат поперёк
// Поперёк это рёбра, угол между которыеми 180 +- альфа
void addDifficultDBEdges(graph &g, size_t &next_idx, const parameters &params, const std::vector<polygon> polygons) {
    logMessage("start to add difficult db edges");

    std::vector<edge> e_to_add;
    std::vector<vertex> v_to_add;

    for (size_t i = 0; i < g.edges.size(); ++i) {
        edge first = g.edges[i];
        if (first.edgeType != E_TRENCH) {
            continue;
        }

        for (size_t j = i + 1; j < g.edges.size(); ++j) {
            edge second = g.edges[j];
            if (second.edgeType != E_TRENCH) {
                continue;
            }

            double cur_angle = first.find_angle(second);
            bool isCrossing = (cur_angle >= 180.0 - BETA || cur_angle <= BETA); 

            checkAndAdd(first, first.u, second.u, g, e_to_add, v_to_add, params, isCrossing, polygons, next_idx);
            checkAndAdd(first, first.u, second.v, g, e_to_add, v_to_add, params, isCrossing, polygons, next_idx);
            checkAndAdd(second, second.u, first.v, g, e_to_add, v_to_add, params, isCrossing, polygons, next_idx);
            checkAndAdd(second, second.v, first.v, g, e_to_add, v_to_add, params, isCrossing, polygons, next_idx);
        }
    }

    g.addEdges(e_to_add);
    g.addVertices(v_to_add);
}

// Добавляем рёбра в граф
// Сначала сортируем старый граф, так, чтобы рёбра шли по порядку
// components - вектор компонент связности графа
void addDBToGraph(graph &g, const std::vector<polygon> &polygons, size_t &next_idx, const parameters &params) {
    logMessage("start building db graph");

    std::pair<graph, std::vector<std::vector<edge>>> return_get_sorted = g.getSorted();
    graph new_graph = return_get_sorted.first;
    std::vector<std::vector<edge>> components = return_get_sorted.second;

    classifyPolygons(components);

    std::vector<edge> e_to_add;
    std::vector<vertex> v_to_add;

    for (const auto &component : components) {
        for (size_t i = 0; i < component.size(); ++i) {
            edge e1 = component[i];
            for (size_t j = (i % component.size()) + 1; j < component.size(); ++j) {
                edge e2 = component[j];
                double angle = e1.find_angle(e2);
                double dist = e1.u.point.distanceTo(e2.u.point);

                if (isGreaterOrEqual(dist, params.min_length_db) && isLessOrEqual(dist, params.max_length_db)) {
                    if (new_graph.mapToDb.find(e1.u.idx) == new_graph.mapToDb.end()) {
                        vertex new_db_v = vertex{e1.u.point, next_idx, VertexType(V_DB)};
                        edge new_trench_to_db_e = edge{
                            e1.idx,
                            new_db_v,
                            e1.u,
                            EdgeType(E_TRENCH_TO_DB),
                            0,
                            (double) params.cost_trench_to_db
                        };
                        e_to_add.push_back(new_trench_to_db_e);
                        new_graph.mapToDb[e1.u.idx] = new_db_v;
                        v_to_add.push_back(new_db_v);
                        next_idx++;
                    } 

                    if (new_graph.mapToDb.find(e2.u.idx) == new_graph.mapToDb.end()) {
                        vertex new_db_v = vertex{e2.u.point, next_idx, VertexType(V_DB)};
                        edge new_trench_to_db_e = edge{
                            e2.u.idx,
                            new_db_v,
                            e2.u,
                            EdgeType(E_TRENCH_TO_DB),
                            0, 
                            (double) params.cost_trench_to_db
                        };
                        e_to_add.push_back(new_trench_to_db_e);
                        new_graph.mapToDb[e2.u.idx] = new_db_v;
                        v_to_add.push_back(new_db_v);
                        next_idx++;
                    }

                    double cost = lenToCost(dist, params.cost_db_per_meter);
                    vertex v1 = new_graph.mapToDb[e1.u.idx];
                    vertex v2 = new_graph.mapToDb[e2.u.idx];

                    bool flag = (components.size() == 1) ? true : component[i].inner;
                    edge new_db_e = edge{v2.idx, v1, v2, EdgeType(E_DB), dist, cost, false, flag};
                    e_to_add.push_back(new_db_e);
                }

                if (angle > EPSILON) {
                    break;
                }      
            }
        }
    }

    

    new_graph.addVertices(v_to_add);
    new_graph.addEdges(e_to_add);

    logMessage("adding difficult db edges");
    addDifficultDBEdges(new_graph, next_idx, params, polygons);

    g = new_graph;
}


// Добавляем доп вершины в базовый граф
void addNewTrenchVerticesUsingLen(graph &g, size_t &next_idx, const parameters &params) {
    logMessage("adding new trench vertices");

    std::vector<edge> e_to_remove;
    std::vector<edge> e_to_add;
    std::vector<vertex> v_to_add;

    for (const auto &e: g.edges) {
        vertex v1 = e.u;
        vertex v2 = e.v;

        double dist = v1.point.distanceTo(v2.point);
        if (dist > params.min_length_db) {

            e_to_remove.push_back(e);

            while (v1.point.distanceTo(v2.point) > params.min_length_db) {
                point new_p = edge{0, v1, v2}.getPointAtDistance(params.min_length_db);
                vertex new_v = vertex{new_p, next_idx, VertexType(V_TRENCH)};
                double dist = v1.point.distanceTo(new_p);
                double cost = lenToCost(dist, params.cost_trench_per_meter);
                edge new_e = edge{next_idx, v1, new_v, EdgeType(E_TRENCH), dist, cost};
                e_to_add.push_back(new_e);
                v_to_add.push_back(new_v);
                v1 = new_v;
                next_idx++;
            }

            double dist = v1.point.distanceTo(v2.point);
            double cost = lenToCost(dist, params.cost_trench_per_meter);
            edge new_e = edge{v2.idx, v1, v2, EdgeType(E_TRENCH), dist, cost};
            e_to_add.push_back(new_e);
        }
    }

    g.eraseEdges(e_to_remove);
    g.addVertices(v_to_add);
    g.addEdges(e_to_add);
}


void buildGraph(const std::vector<polygon> &polygons, graph &g, const parameters &params) {
    logMessage("start building graph");
    
    size_t next_idx = 0;
    for (const auto &poly : polygons) {
        buildBaseTrenchGraph(g, poly, next_idx, params);
    }

    logMessage("base trench graph builded");

    handleIntersectingPolygons(g, polygons, next_idx, params);

    addNewTrenchVerticesUsingLen(g, next_idx, params);

    addDBToGraph(g, polygons, next_idx, params);

    logMessage("graph builded");
}

int main(int argc, char *argv[]) {
    std::string inputFile;
    std::string outputFile;
    std::string configFile;
    bool useArc = false;

    if (argc < 4) {
        std::cout << "Usage: " << argv[0] << " <input file> <output file> <config file> [--usearc]" << std::endl;
        std::cout << "--usearc generate arcs instead of db straight edges";
        return 1;
    }

    inputFile = argv[1];
    outputFile = argv[2];
    configFile = argv[3];

    if (argc == 5 && std::strcmp(argv[4], "--usearc") == 0) {
        useArc = true;
    }

    parameters params;
    if (readConfig(configFile, params) != 0) {
        std::cerr << "Invalid config file" << std::endl;
        return 1;
    };

    std::vector<polygon> polygons = readGeoJSON(inputFile);

    graph g;
    buildGraph(polygons, g, params);

    nlohmann::json output_json = getJsonFromGraph(g, useArc);
    writeToJSON(output_json, outputFile);

    return 0;
}
