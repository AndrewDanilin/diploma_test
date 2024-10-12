#include <utility>
#include <vector>
#include <algorithm>
#include <map>
#include "../headers/graph.h"


namespace MyNamespace {
    void graph::addEdge(const vertex &v1, const vertex &v2, EdgeType edgeType, const double &dist, const double &cost, const bool cross, const bool inner) {
        edge edge_to_add_first = edge{v2.idx, v1, v2, edgeType, dist, cost, cross, inner};
        edge edge_to_add_second = edge{v1.idx, v2, v1, edgeType, dist, cost, cross, inner};

        if (map.find(v1.idx) == map.end()) {
            map[v1.idx] = {};
        }
        map[v1.idx].push_back(edge_to_add_first);

        if (map.find(v2.idx) == map.end()) {
            map[v2.idx] = {};
        }
        map[v2.idx].push_back(edge_to_add_second);

        edges.push_back(edge_to_add_first);
    }

    void graph::addVertex(const vertex &u) {
        vs.push_back(u);
    }    

    void graph::deleteEdgeFromMap(const edge &e) {
        size_t idx1 = e.u.idx;
        size_t idx2 = e.v.idx;

        for (size_t i = 0; i < map[idx1].size(); ++i) {
            if (map[idx1][i].idx == idx2) {
                map[idx1].erase(std::next(map[idx1].begin(), i));
                if (map[idx1].empty()) {
                    map.erase(idx1);
                }
                break;
            }
        }

        for (size_t i = 0; i < map[idx2].size(); ++i) {
            if (map[idx2][i].idx == idx1) {
                map[idx2].erase(std::next(map[idx2].begin(), i));
                if (map[idx2].empty()) {
                    map.erase(idx2);
                }
                break;
            }
        }
    }

    void graph::deleteFromMapByIndex(const size_t &idx) {
        map.erase(idx);
        std::vector<size_t> to_delete;

        for (auto &m : map) {
            for (size_t i = 0; i < m.second.size(); ++i) {
                if (m.second[i].idx == idx) {
                    m.second.erase(m.second.begin() + i);
                    break;
                }
            }
            if (m.second.empty()) {
                to_delete.push_back(m.first);
            }
        }

        for (const auto &i: to_delete) {
            map.erase(i);
        }

    }

    void graph::eraseEdges(const std::vector<edge> es) {
        edges.erase(std::remove_if(edges.begin(), edges.end(),
        [&es](const edge& e) {
            return std::find(es.begin(), es.end(), e) != es.end();
        }),
        edges.end());

        for (const auto &e: es) {
            deleteEdgeFromMap(e);
        }
    }

    void graph::addEdges(const std::vector<edge> es) {
        for (const auto &e: es) {
            addEdge(e.u, e.v, e.edgeType, e.len, e.cost, e.cross, e.inner);
        }
    }

    void graph::addVertices(const std::vector<vertex> vss) {
        for (const auto &v: vss) {
            addVertex(v);
        }
    }

    void graph::eraseVertices(const std::vector<vertex> vertices) {
        vs.erase(std::remove_if(vs.begin(), vs.end(),
        [&vertices](const vertex& v) {
            return std::find(vertices.begin(), vertices.end(), v) != vertices.end();
        }),
        vs.end());

        for (const auto &v: vertices) {
            deleteFromMapByIndex(v.idx);
        }

        edges.erase(std::remove_if(edges.begin(), edges.end(),
        [&vertices](const edge& e) {
            return std::find(vertices.begin(), vertices.end(), e.u) != vertices.end() ||
                std::find(vertices.begin(), vertices.end(), e.v) != vertices.end();
        }),
        edges.end());
    }

    void graph::dfs(const vertex &v,
        std::vector<bool> &visited,
        graph &new_g,
        size_t next_idx, 
        vertex &first_v, 
        const vertex &f_v, 
        const size_t cnt
    ) {
        visited[v.idx] = true;
        vertex new_v = vertex{v.point, next_idx, VertexType{V_TRENCH}};

        if (cnt == 0) {
            first_v = new_v;
        }

        new_g.addVertex(new_v);

        for (const auto &e: map[v.idx]) {
            if (!visited[e.idx]) {
                vertex next_v = vertex{e.v.point, next_idx + 1, VertexType(V_TRENCH)};
                new_g.addEdge(new_v, next_v, EdgeType(E_TRENCH), e.len, e.cost);
                dfs(e.v, visited, new_g, next_idx + 1, first_v, f_v, cnt + 1);
            } else if (cnt > 2 && e.idx == f_v.idx) {
                new_g.addEdge(new_v, first_v, EdgeType(E_TRENCH), e.len, e.cost);
                return;
            }
        }
    }

    std::pair<graph, std::vector<std::vector<edge>>> graph::getSorted() {
        graph new_g;
        std::vector<std::vector<edge>> components;
        size_t max_idx = 0;
        for (const auto &v: vs) {
            if (v.idx > max_idx) {
                max_idx = v.idx;
            }
        }
        std::vector<bool> visited = std::vector(max_idx + 1, false);
        size_t size_edges = new_g.edges.size();
        size_t next_idx = 0;
        for (const auto &v: vs) {
            if (!visited[v.idx]) {
                vertex c_v = vertex{};
                dfs(v, visited, new_g, size_edges + next_idx, c_v, v, 0);
                std::vector<edge> cur_component = {new_g.edges.begin() + size_edges, new_g.edges.end()};
                components.push_back(cur_component);
                size_edges = new_g.edges.size();
            }
        }

        return {new_g, components};
    }
}