/*
 * distance table with lazy evaluation, using BFS
 */
#pragma once

#include "graph.hpp"
#include "instance.hpp"
#include "utils.hpp"

struct DistTable {
  const uint V_size;  // number of vertices
  std::vector<std::vector<uint> >
      table;          // distance table, index: agent-id & vertex-id

  std::vector<std::vector<double> >
      double_table;  

  std::vector<std::queue<Vertex*> > OPEN;  // search queue

  inline uint get(uint i, uint v_id);      // agent, vertex-id
  uint get(uint i, Vertex* v);             // agent, vertex
  inline double get_heuristic(uint i, uint v_id);
  double get_heuristic(uint i, Vertex* v);   
  DistTable(const Instance& ins);
  DistTable(const Instance* ins);

  void setup(const Instance* ins);  // initialization
  void reset(const Instance* ins); 
  double distance(uint a, uint b) const;   // NEW: access edge weights
  void set_edge_weight(uint a, uint b, double w); // NEW: modify edge weights
  void increase_edge_weight(uint a, uint b, double w); // NEW: increase edge weights
  double get_edge_weight(uint a, uint b);
  void print_edge_map();
  void clear_traffic() {
    for(int i = 0; i < edge_map->data.size(); i++){
      edge_map->data[i] = 1.0;
    }
    for(int i = 0; i < directed_edge_map->data.size(); i++){
      directed_edge_map->data[i] = 1.0;
    }
  }

  double directed_distance(uint a, uint b) const; 
  void set_directed_edge_weight(uint a, uint b, double w); // NEW: modify edge weights
  void increase_directed_edge_weight(uint a, uint b, double w); // NEW: increase edge weights
  double get_directed_edge_weight(uint a, uint b);

  private:

  // 👇 Private nested struct to manage edge weights
  struct EdgeMap {
    int width, height;
    std::vector<double> data;

    EdgeMap(int w, int h) : width(w), height(h) {
      int total_edges = (w - 1) * h + w * (h - 1);
      data.resize(total_edges, 1.0);
      // std::cout<< total_edges <<std::endl;
    }

    int edge_index(uint a, uint b) const {
      if (a > b) std::swap(a, b);
      int ax = a % width, ay = a / width;
      int bx = b % width, by = b / width;

      if (bx == ax + 1 && by == ay) {
        return ay * (width - 1) + ax; // horizontal edge
      } else if (by == ay + 1 && bx == ax) {
        int base = (width - 1) * height;
        return base + ay * width + ax; // vertical edge
      } else {
        return -1; // not adjacent
      }
    }

    double get(uint a, uint b) const {
      int idx = edge_index(a, b);
      if (idx == -1) return std::numeric_limits<double>::infinity();
      return data[idx];
    }

    void set(uint a, uint b, double w) {
      int idx = edge_index(a, b);
      if (idx != -1) data[idx] = w;
    }


    void increase(uint a, uint b, double w) {
      int idx = edge_index(a, b);
      if (idx != -1){
        data[idx] += w;
      }
    }

  };



  struct DirectedEdgeMap {
    int width, height;
    std::vector<double> data;
  
    DirectedEdgeMap(int w, int h) : width(w), height(h) {
      // Each undirected edge becomes two directed edges
      int undirected_edges = (w - 1) * h + w * (h - 1);
      data.resize(undirected_edges * 2, 1.0); // *2 for both directions
    }
  
    // Returns the base index of the undirected edge, or -1 if not adjacent
    int base_edge_index(uint a, uint b) const {
      int ax = a % width, ay = a / width;
      int bx = b % width, by = b / width;
  
      if (bx == ax + 1 && by == ay) {
        return ay * (width - 1) + ax; // horizontal edge
      } else if (bx == ax - 1 && by == ay) {
        return ay * (width - 1) + bx; // reverse horizontal
      } else if (by == ay + 1 && bx == ax) {
        int base = (width - 1) * height;
        return base + ay * width + ax; // vertical edge
      } else if (by == ay - 1 && bx == ax) {
        int base = (width - 1) * height;
        return base + by * width + ax; // reverse vertical
      } else {
        return -1; // not adjacent
      }
    }
  
    // Return the index in data[] for a directed edge a -> b
    int directed_edge_index(uint a, uint b) const {
      int base = base_edge_index(a, b);
      if (base == -1) return -1;
  
      return 2 * base + ((a < b) ? 0 : 1); // direction matters
    }
  
    double get(uint a, uint b) const {
      int idx = directed_edge_index(a, b);
      if (idx == -1){
        std::cout<< "Invalid edge index: " << a << " -> " << b << std::endl;
        return std::numeric_limits<double>::infinity();
      }
      return data[idx];
    }
  
    void set(uint a, uint b, double w) {
      int idx = directed_edge_index(a, b);
      if (idx != -1) data[idx] = w;
    }
  
    void increase(uint a, uint b, double w) {
      int idx = directed_edge_index(a, b);
      if (idx != -1) data[idx] += w;
    }
  };

  DirectedEdgeMap* directed_edge_map = nullptr;
  EdgeMap* edge_map = nullptr;
};
