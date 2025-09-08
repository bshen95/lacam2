#include "../include/AstarDist_table.hpp"

AstarDistTable::AstarDistTable(const Instance& ins)
    : V_size(ins.G.V.size()), 
    W(ins.G.width),
    starts(ins.starts),
    V_Node_table(ins.N, std::vector<V_Node>(V_size)),
    G(&ins.G)
{
  setup(&ins);
}

AstarDistTable::AstarDistTable(const Instance* ins)
    : V_size(ins->G.V.size()), 
    W(ins->G.width),
    starts(ins->starts),
    V_Node_table(ins->N, std::vector<V_Node>(V_size)),
    G(&ins->G)
{
  setup(ins);
}

void AstarDistTable::setup(const Instance* ins)
{
  OPEN.resize(ins->N);
  for (size_t i = 0; i < ins->N; ++i) {
    // avoid push back copy;
    // OPEN[i] = pqueue<V_Node, cmp_less_f,min_q>(V_size);
    for (size_t j = 0; j < V_size; ++j) {
      V_Node_table[i][j].expanded = false;
      V_Node_table[i][j].in_queue = false;
      V_Node_table[i][j].generated = false;
      V_Node_table[i][j].v= ins->G.V[j];
    }
    auto n = ins->goals[i];
    OPEN[i].push(&V_Node_table[i][n->id]);
    V_Node_table[i][n->id].g = 0;
    V_Node_table[i][n->id].h = manhattan_dist(n,starts[i]);
    V_Node_table[i][n->id].update_f();
    V_Node_table[i][n->id].in_queue = true;
    V_Node_table[i][n->id].generated = true;
  }

  if (edge_map) delete edge_map;
  edge_map = new EdgeMap(ins->G.width, ins->G.height);
  if(edge_map_list.size() != 0){
    for (int i = 0; i < edge_map_list.size(); i++){
      delete edge_map_list[i];
    }
    edge_map_list.clear();
  }
  edge_map_list.resize(ins->N);
  for (int i = 0; i < ins->N; i++){
    edge_map_list[i] = new EdgeMap(ins->G.width, ins->G.height);
  } 
}

void AstarDistTable::reset(const Instance* ins)
{
  for (size_t i = 0; i < ins->N; ++i) {
    OPEN[i].clear();
    for (size_t j = 0; j < V_size; ++j) {
      V_Node_table[i][j].expanded = false;
      V_Node_table[i][j].generated = false;
    }
    auto n = ins->goals[i];
    OPEN[i].push(&V_Node_table[i][n->id]);
    V_Node_table[i][n->id].g = 0;
    V_Node_table[i][n->id].h = manhattan_dist(n,starts[i]);
    V_Node_table[i][n->id].update_f();
    V_Node_table[i][n->id].generated = true;
  }
}



void AstarDistTable::test_dijkstra(uint i, const Instance* ins) {
  // Priority queue for Dijkstra's algorithm
  using QueueElement = std::pair<double, uint>; // {distance, vertex_id}
  std::priority_queue<QueueElement, std::vector<QueueElement>, std::greater<>> pq;

  // Distance vector initialized to infinity
  std::vector<double> distances(V_size, std::numeric_limits<double>::infinity());
  distances[ins->goals[i]->id] = 0;

  // Push the start node into the queue
  pq.push({0, ins->goals[i]->id});

  while (!pq.empty()) {
      auto [current_dist, current_id] = pq.top();
      pq.pop();

      // Skip if the distance is outdated
      if (current_dist > distances[current_id]) {
          continue;
      }

      // Expand neighbors
      for (auto&& neighbor : G->V[current_id]->neighbor) {
          uint neighbor_id = neighbor->id;
          double edge_weight = get_edge_weight(ins->G.V[current_id]->index, neighbor->index);
          double new_dist = current_dist + edge_weight;

          if (new_dist < distances[neighbor_id]) {
              distances[neighbor_id] = new_dist;
              pq.push({new_dist, neighbor_id});
          }
      }
  }
  // Compare Dijkstra's results with get_individual_heuristic
  for (uint v_id = 0; v_id < V_size; ++v_id) {
      double heuristic = get_heuristic(i, v_id);
      if (distances[v_id] != heuristic) {
          std::cout << "Mismatch at vertex " << v_id
                    << ": Dijkstra = " << distances[v_id]
                    << ", Heuristic = " << heuristic << std::endl;
      }
  }
}

double AstarDistTable::get_individual_heuristic(uint i, uint v_id)
{
  if (V_Node_table[i][v_id].expanded){
    return V_Node_table[i][v_id].g;
  } 
  /*
   * c.f., Reverse Resumable A*
   * https://www.aaai.org/Papers/AIIDE/2005/AIIDE05-020.pdf
   *
   * sidenote:
   * tested RRA* but lazy BFS was much better in performance
   */

  while (!OPEN[i].empty()) {
    auto & current_queue = OPEN[i];
    V_Node* n = OPEN[i].pop();
    const double d_n = n->g;
    V_Node_table[i][n->v->id].expanded = true;
    V_Node_table[i][n->v->id].in_queue = false;
    for (auto &&m : n->v->neighbor) {
      double g = d_n + get_edge_weight(i,n->v->index, m->index);
      double h = manhattan_dist(m, starts[i]);
      if (!V_Node_table[i][m->id].in_queue) {
        V_Node_table[i][m->id].g = g ;
        V_Node_table[i][m->id].h = h;
        V_Node_table[i][m->id].f = g + h;
        V_Node_table[i][m->id].predecessor = n->v->id;
        V_Node_table[i][m->id].in_queue = true;
        OPEN[i].push(&V_Node_table[i][m->id]);
      }else{
        if (g + h < V_Node_table[i][m->id].f){
          V_Node_table[i][m->id].g = g ;
          V_Node_table[i][m->id].h = h;
          V_Node_table[i][m->id].f = g + h;
          V_Node_table[i][m->id].predecessor = n->v->id;
          OPEN[i].decrease_key(&V_Node_table[i][m->id]);
        }
      }
    }
    if (n->v->id == v_id){
      return d_n;
    }
  }
  return V_size;
}


double AstarDistTable::update_heuristic_table(const Instance* ins, uint i, uint v_index, double g_value)
{
  uint v_id = G->U[v_index]->id;
  uint goal_id = ins->goals[i]->id;
  std::vector<uint> path;
  // get shortest path ;
  if (V_Node_table[i][v_id].expanded){
    // already in the table; 
    while(v_id != goal_id ){
      path.push_back(v_id);
      v_id = V_Node_table[i][v_id].predecessor;
    }
    path.push_back(goal_id);
  } else{
    get_individual_heuristic(i, v_id);
    while(v_id != goal_id ){
      path.push_back(v_id);
      v_id = V_Node_table[i][v_id].predecessor;
    }
    path.push_back(goal_id);
  }

  double avg_g_value = g_value / (path.size()- 1);
  for(int j = 0; j < path.size()-1; j++){
    uint v = path[j];
    uint next_v = path[j+1];
    uint edge_index = G->V[v]->index;
    uint next_edge_index = G->V[next_v]->index;
    auto current_weight = get_edge_weight(i, G->V[v]->index, G->V[next_v]->index);
    increase_edge_weight(i, G->V[v]->index, G->V[next_v]->index, 
      avg_g_value);
  }
}

double AstarDistTable::get_heuristic(uint i, uint v_id)
{

  if (V_Node_table[i][v_id].expanded){
    return V_Node_table[i][v_id].g;
  } 
  /*
   * c.f., Reverse Resumable A*
   * https://www.aaai.org/Papers/AIIDE/2005/AIIDE05-020.pdf
   *
   * sidenote:
   * tested RRA* but lazy BFS was much better in performance
   */

  while (!OPEN[i].empty()) {
    auto & current_queue = OPEN[i];
    V_Node* n = OPEN[i].pop();
    const double d_n = n->g;
    V_Node_table[i][n->v->id].expanded = true;
    for (auto &&m : n->v->neighbor) {
      double g = d_n + get_edge_weight(n->v->index, m->index);
      double h = manhattan_dist(m, starts[i]);
      if (!V_Node_table[i][m->id].generated) {
        V_Node_table[i][m->id].g = g ;
        V_Node_table[i][m->id].h = h;
        V_Node_table[i][m->id].f = g + h;
        V_Node_table[i][m->id].predecessor = n->v->id;
        V_Node_table[i][m->id].generated = true;
        OPEN[i].push(&V_Node_table[i][m->id]);
      }else{
        if (g + h < V_Node_table[i][m->id].f){
          V_Node_table[i][m->id].g = g ;
          V_Node_table[i][m->id].h = h;
          V_Node_table[i][m->id].f = g + h;
          V_Node_table[i][m->id].predecessor = n->v->id;
          OPEN[i].decrease_key(&V_Node_table[i][m->id]);
        }
      }
    }
    if (n->v->id == v_id){
      return d_n;
    }
  }
  return V_size;
}

int AstarDistTable::manhattan_dist(const Vertex* v, const Vertex* n){
  return abs(static_cast<int>(v->index / W - n->index / W)) + 
         abs(static_cast<int>(v->index % W - n->index % W));
}

double AstarDistTable::get_edge_weight(uint a, uint b){
  return edge_map->get(a, b);
}

void AstarDistTable::set_edge_weight(uint a, uint b, double w) {
  edge_map->set(a, b, w);
}

void AstarDistTable::increase_edge_weight(uint a, uint b, double w) {
  edge_map->increase(a, b, w);
}

double AstarDistTable::get_edge_weight(uint agent_id, uint a, uint b){
  return edge_map_list[agent_id]->get(a, b);
}

void AstarDistTable::set_edge_weight(uint agent_id, uint a, uint b, double w) {
  edge_map_list[agent_id]->set(a, b, w);
}

void AstarDistTable::increase_edge_weight(uint agent_id, uint a, uint b, double w) {
  edge_map_list[agent_id]->increase(a, b, w);
}


double AstarDistTable::get_heuristic(uint i, Vertex* v) { 
  // std::cout<< "getting heuristic:" <<get_heuristic(i, v->id) <<std::endl;
  return get_heuristic(i, v->id);
}


double AstarDistTable::get_individual_heuristic(uint i, Vertex* v) { 
  // std::cout<< "getting heuristic:" <<get_heuristic(i, v->id) <<std::endl;
  return get_individual_heuristic(i, v->id);
}