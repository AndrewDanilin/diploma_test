#pragma once

#include <cstddef>
#include <vector>
#include "vertex.h"
#include "edge.h"

namespace MyNamespace {

    struct graph {
        std::vector<vertex> vs;
        std::vector<edge> edges;
        std::map<size_t, std::vector<edge>> map;
        std::map<size_t, vertex> mapToDb;

        void addEdge(const vertex &v1,
            const vertex &v2,
            EdgeType edgeType,
            const double &dist,
            const double &cost,
            const bool cross = false,
            const bool inner = false);
          
        void addVertex(const vertex &u);
        void deleteEdgeFromMap(const edge &e);
        void deleteFromMapByIndex(const size_t &idx);
        void eraseEdges(const std::vector<edge> es);
        void addEdges(const std::vector<edge> es);
        void addVertices(const std::vector<vertex> vss);
        void eraseVertices(const std::vector<vertex> vertices);
        void dfs(const vertex &v, 
            std::vector<bool> &visited, 
            graph &new_g, 
            size_t next_idx, 
            vertex &first_v, 
            const vertex &f_v, 
            const size_t cnt
        ); 
        std::pair<graph, std::vector<std::vector<edge>>> getSorted();

    };

}