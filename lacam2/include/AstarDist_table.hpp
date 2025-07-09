#pragma once

#include "graph.hpp"
#include "instance.hpp"
#include "utils.hpp"
#include "heap.hpp"

struct AstarDistTable {
  const uint V_size;  // number of vertices
  const uint W;
  const std::vector<Vertex*> starts;
  const Graph* G;
  std::vector<std::vector<V_Node> >
      V_Node_table;;  

  std::vector<pqueue<V_Node,cmp_less_f,min_q>> OPEN; 

  int manhattan_dist(const Vertex* v, const Vertex* n);

  inline double get_individual_heuristic(uint i, uint v_id);
  double get_individual_heuristic(uint i, Vertex* v);  
  inline double get_heuristic(uint i, uint v_id);
  double get_heuristic(uint i, Vertex* v);   
  double update_heuristic_table(const Instance* ins, uint i, uint v_index, double g_value);

  AstarDistTable(const Instance& ins);
  AstarDistTable(const Instance* ins);

  void setup(const Instance* ins);  // initialization
  
  void reset(const Instance* ins); 

  void print_edge_map(){
    edge_map->print_edge_map();
  }

  void clear_traffic() {
    for(int i = 0; i < edge_map->data.size(); i++){
      edge_map->data[i] = 1.0;
    }
  }

  void clear_all_traffic() {
    for(int i = 0; i < edge_map_list.size(); i++){
      for(int j = 0; j < edge_map_list[i]->data.size(); j++){
        edge_map_list[i]->data[j] = 1.0;
      }
    }
  }

  void compare_edge_map(){
    for(int i = 0; i < edge_map->data.size(); i++){
      if(edge_map->data[i] != edge_map_list[0]->data[i]){
        std::cout<< "edge map not equal!" <<std::endl;
        std::cout<< edge_map->data[i] << "," << edge_map_list[0]->data[i] <<std::endl;
      }
    }
  }
  void test_dijkstra(uint i, const Instance* ins) ;
  
  void set_edge_weight(uint a, uint b, double w); 

  void increase_edge_weight(uint a, uint b, double w); 

  double get_edge_weight(uint a, uint b);
  
  
  void set_edge_weight(uint agent_id, uint a, uint b, double w);
  
  void increase_edge_weight(uint agent_id, uint a, uint b, double w); 

  double get_edge_weight(uint agent_id, uint a, uint b);

  double decay_local_entries(double decay_factor){
    edge_map->decay_local_entries(decay_factor);
  }

  void copy_global_data_and_clean_local(){
    edge_map->copy_global_data_and_clean_local();
  }

  void sum_global_entries(){
    edge_map->sum_global_entries();
  }

  void increase_local_edge_weight( uint a, uint b, double w) {
    edge_map->increase_local_edge_weight(a, b, w);
  }

  void apply_gaussian_filter(double sigma, const Graph* G) {
    edge_map->apply_graph_2hop_filter(1.0, G);
  }


  private:

  struct EdgeMap {
    int width, height;
    std::vector<double> data;
    std::vector<double> global_data;
    std::vector<double> local_data;


    EdgeMap(int w, int h) : width(w), height(h) {
      int total_edges = (w - 1) * h + w * (h - 1);
      data.resize(total_edges, 1.0);
      global_data.resize(total_edges, 0.0);
      local_data.resize(total_edges, 0.0);
      // std::cout<< total_edges <<std::endl;
    }

    void sum_global_entries(){
      for (size_t i = 0; i < data.size(); ++i) {
        data[i] = global_data[i] + local_data[i];
      }
    }

    void copy_global_data_and_clean_local(){
      global_data = data;
      for (size_t i = 0; i < local_data.size(); ++i) {
        local_data[i] = 0.0;
      }
    }

    void increase_local_edge_weight(uint a, uint b, double w) {
      int idx = edge_index(a, b);
      if (idx != -1){
        local_data[idx] += w;
      }
    }

    // Add to EdgeMap struct:
void apply_graph_2hop_filter(double sigma, const Graph* G) {
    if (sigma <= 0) return;
    
    std::vector<double> filtered_data = data;
    
    // For each vertex in the graph
    for (size_t v = 0; v < G->U.size(); v++) {
        // Get all 2-hop neighbor edges
        std::unordered_map<int, int> neighbor_edges;  // edge_idx -> hop_distance
        
        if(G->U[v] == nullptr) continue;  // skip if vertex is null
        // First hop - direct neighbors
        for (auto&& n1 : G->U[v]->neighbor) {
            int edge_idx = edge_index(v, n1->index);
            if (edge_idx >= 0) {
                neighbor_edges[edge_idx] = 1;
            }
            
            // Second hop - neighbors of neighbors
            for (auto&& n2 : n1->neighbor) {
                int edge_idx2 = edge_index(n1->index, n2->index);
                if (edge_idx2 >= 0 && neighbor_edges.count(edge_idx2) == 0) {
                    neighbor_edges[edge_idx2] = 2;
                }
            }
        }
        
        // Apply Gaussian weights to each edge connected to vertex v
        for (auto&& n : G->U[v]->neighbor) {
            int current_edge = edge_index(v, n->index);
            if (current_edge < 0) continue;
            
            double sum = data[current_edge];  // include self
            double weight_sum = 1.0;          // weight for self
            
            // Apply weights based on hop distance
            for (const auto& [n_edge, hops] : neighbor_edges) {
                if (n_edge != current_edge) {  // skip self
                    double dist = hops;  // use hop count as distance
                    double weight = exp(-(dist * dist) / (2 * sigma * sigma));
                    sum += data[n_edge] * weight;
                    weight_sum += weight;
                }
            }
            filtered_data[current_edge] = sum / weight_sum;
        }
    }
    data = filtered_data;
  }


    // int edge_index(uint a, uint b) const {
    //   // input U_index, width * y + x;
    //   if (a > b) std::swap(a, b);
    //   int ax = a % width, ay = a / width;
    //   int bx = b % width, by = b / width;

    //   if (bx == ax + 1 && by == ay) {
    //     return ay * (width - 1) + ax; // horizontal edge
    //   } else if (by == ay + 1 && bx == ax) {
    //     int base = (width - 1) * height;
    //     return base + ay * width + ax; // vertical edge
    //   } else {
    //     std::cout<< "This should never happened" <<std::endl;
    //     return -1; // not adjacent
    //   }
    // }


    int edge_index(uint a, uint b) const {
      if (a > b) std::swap(a, b);
      int ax = a % width, ay = a / width;
      int bx = b % width, by = b / width;

      if (bx == ax + 1 && by == ay) {  // horizontal edge
          // First all horizontal edges in this row
          return ay * (2*width-1) + ax;
      } else if (by == ay + 1 && bx == ax) {  // vertical edge
          // Then vertical edges in this row
          return ay * (2*width-1) + (width-1) + ax;
      }else {
          std::cout<< "This should never happened" <<std::endl;
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

    void decay_local_entries(double decay_factor){
      for (size_t i = 0; i < data.size(); ++i) {
        local_data[i] *= decay_factor;
      }
    }

    void print_edge_map() const {
      for (size_t i = 0; i < data.size(); ++i) {
        if(data[i] > 1.0){
          std::cout<< "edge: " << i << " weight: " << data[i] << std::endl;
        } 
      }
      std::cout << std::endl;
    }
  };

  EdgeMap* edge_map = nullptr;
  std::vector<EdgeMap*> edge_map_list;
};