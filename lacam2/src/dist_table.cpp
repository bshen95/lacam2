#include "../include/dist_table.hpp"

DistTable::DistTable(const Instance& ins)
    : V_size(ins.G.V.size()), 
    table(ins.N, std::vector<uint>(V_size, V_size*V_size)), 
    double_table(ins.N, std::vector<double>(V_size, V_size*V_size))
{
  setup(&ins);
}

DistTable::DistTable(const Instance* ins)
    : V_size(ins->G.V.size()), 
    table(ins->N, std::vector<uint>(V_size, V_size*V_size)),
    double_table(ins->N, std::vector<double>(V_size, V_size*V_size))
{
  setup(ins);
}

void DistTable::setup(const Instance* ins)
{
  for (size_t i = 0; i < ins->N; ++i) {
    OPEN.push_back(std::queue<Vertex*>());
    auto n = ins->goals[i];
    OPEN[i].push(n);
    table[i][n->id] = 0;
    double_table[i][n->id] = 0;
  }

  if (edge_map) delete edge_map;
  edge_map = new EdgeMap(ins->G.width, ins->G.height);

  if (directed_edge_map) delete directed_edge_map;
  directed_edge_map = new DirectedEdgeMap(ins->G.width, ins->G.height);
}

void DistTable::reset(const Instance* ins)
{
  OPEN.clear();
  for (size_t i = 0; i < ins->N; ++i) {
    std::fill(table[i].begin(), table[i].end(), V_size*V_size);
    std::fill(double_table[i].begin(), double_table[i].end(), V_size*V_size);
    OPEN.push_back(std::queue<Vertex*>());
    auto n = ins->goals[i];
    OPEN[i].push(n);
    table[i][n->id] = 0;
    double_table[i][n->id] = 0;
  }
}

uint DistTable::get(uint i, uint v_id)
{
  if (table[i][v_id] < V_size*V_size) return table[i][v_id];

  /*
   * BFS with lazy evaluation
   * c.f., Reverse Resumable A*
   * https://www.aaai.org/Papers/AIIDE/2005/AIIDE05-020.pdf
   *
   * sidenote:
   * tested RRA* but lazy BFS was much better in performance
   */

  while (!OPEN[i].empty()) {
    auto&& n = OPEN[i].front();
    OPEN[i].pop();
    const int d_n = table[i][n->id];
    for (auto&& m : n->neighbor) {
      const int d_m = table[i][m->id];
      if (d_n + 1 >= d_m) continue;
      table[i][m->id] = d_n + 1;
      OPEN[i].push(m);
    }
    if (n->id == v_id) return d_n;
  }
  // std::cout<<"I am here"<<std::endl;
  return V_size;
}

double DistTable::get_heuristic(uint i, uint v_id)
{
  if (double_table[i][v_id] < V_size*V_size) return double_table[i][v_id];

  /*
   * BFS with lazy evaluation
   * c.f., Reverse Resumable A*
   * https://www.aaai.org/Papers/AIIDE/2005/AIIDE05-020.pdf
   *
   * sidenote:
   * tested RRA* but lazy BFS was much better in performance
   */

  while (!OPEN[i].empty()) {
    auto&& n = OPEN[i].front();
    OPEN[i].pop();
    const double d_n = double_table[i][n->id];
    for (auto&& m : n->neighbor) {
      const double d_m = double_table[i][m->id];
      if (d_n + directed_distance(n->index, m->index) >= d_m) continue;
      double_table[i][m->id] = d_n + directed_distance(n->index, m->index);
      OPEN[i].push(m);
    }
    if (n->id == v_id) return d_n;
  }
  return V_size;
}

double DistTable::distance(uint a, uint b) const {
  return edge_map->get(a, b);
}

double DistTable::get_edge_weight(uint a, uint b){
  return edge_map->get(a, b);
}

void DistTable::set_edge_weight(uint a, uint b, double w) {
  edge_map->set(a, b, w);
}

void DistTable::increase_edge_weight(uint a, uint b, double w) {
  edge_map->increase(a, b, w);
}


double DistTable::directed_distance(uint a, uint b) const {
  // std::cout<<"getting distance:" <<directed_edge_map->get(a, b) <<std::endl;
  return directed_edge_map->get(a, b);
}

void DistTable::print_edge_map() {
  for(auto e : directed_edge_map->data){
    std::cout<< e <<std::endl;
  }
}



double DistTable::get_directed_edge_weight(uint a, uint b){
  return directed_edge_map->get(a, b);
}

void DistTable::set_directed_edge_weight(uint a, uint b, double w) {
  directed_edge_map->set(a, b, w);
}

void DistTable::increase_directed_edge_weight(uint a, uint b, double w) {
  directed_edge_map->increase(a, b, w);
}


double DistTable::get_heuristic(uint i, Vertex* v) { 
  // std::cout<< "getting heuristic:" <<get_heuristic(i, v->id) <<std::endl;
  return get_heuristic(i, v->id);
}
uint DistTable::get(uint i, Vertex* v) { return get(i, v->id); }
