#include "../include/planner.hpp"
#include <unordered_set>
#include <map>
LNode::LNode(LNode* parent, uint i, Vertex* v)
    : who(), where(), depth(parent == nullptr ? 0 : parent->depth + 1)
{
  if (parent != nullptr) {
    who = parent->who;
    who.push_back(i);
    where = parent->where;
    where.push_back(v);
  }
}

uint HNode::HNODE_CNT = 0;

// for high-level
HNode::HNode(const Config& _C, DistTable& D, HNode* _parent, const uint _g,
             const uint _h)
    : C(_C),
      parent(_parent),
      neighbor(),
      g(_g),
      h(_h),
      f(g + h),
      priorities(C.size()),
      order(C.size(), 0),
      search_tree(std::queue<LNode*>())
{
  ++HNODE_CNT;

  search_tree.push(new LNode());
  const auto N = C.size();

  // update neighbor
  if (parent != nullptr) parent->neighbor.insert(this);

  // set priorities
  if (parent == nullptr) {
    // initialize
    for (uint i = 0; i < N; ++i) priorities[i] = (float)D.get(i, C[i]) / N;
    // std::random_device rd; // Seed
    // std::mt19937 gen(rd()); // Mersenne Twister generator
    // std::uniform_real_distribution<> dis(0.0, 1.0);
    // for (uint i = 0; i < N; ++i) priorities[i] = dis(gen);
  } else {
    // dynamic priorities, akin to PIBT

    for (size_t i = 0; i < N; ++i) {
      if (D.get(i, C[i]) != 0) {
        priorities[i] = parent->priorities[i] + 1;
      } else {
        priorities[i] = parent->priorities[i] - (int)parent->priorities[i];
      }
    }

    
  }

  //set order
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](uint i, uint j) { return priorities[i] > priorities[j]; });
}

HNode::~HNode()
{
  while (!search_tree.empty()) {
    delete search_tree.front();
    search_tree.pop();
  }
}

MCTNode::MCTNode(MCTNode* _parent, HNode* _HNode, uint _node_id):
      parent( _parent),
      hNode(_HNode),
      visits(0),
      reward(0), 
      node_id(_node_id)
{};


MCTNode::~MCTNode()
{
  delete hNode;
}

Planner::Planner(const Instance* _ins, const Deadline* _deadline,
                 std::mt19937* _MT, const int _verbose,
                 const Objective _objective, const float _restart_rate, const std::string save_tree_file)
    : ins(_ins),
      deadline(_deadline),
      MT(_MT),
      verbose(_verbose),
      objective(_objective),
      RESTART_RATE(_restart_rate),
      N(ins->N),
      V_size(ins->G.size()),
      D(DistTable(ins)),
      PIBT_D(AstarDistTable(ins)),
      traffic_map(&ins->G),
      time_period_traffic_map(time_bucket_size, TrafficMap(&ins->G) ), // Initialize time-period based
      astar_search(&traffic_map, &ins->G),
      guidance_heuristic( &ins->G, &traffic_map, N),
      Q_tables(N, QTable(ins->G.U.size(), 5, -1)),
      loop_cnt(0),
      C_next(N),
      tie_breakers(V_size, 0),
      A(N, nullptr),
      occupied_now(V_size, nullptr),
      occupied_next(V_size, nullptr),
      tree_file(save_tree_file),
      individual_transition_frequency(N)
{

}

Planner::~Planner() {}


void Planner::propagate_order_to_neighbors(HNode* current_node) {
  // Iterate through all neighbors of the current node
  for (auto neighbor : current_node->neighbor) {
    // Update the neighbor's order based on the current node's order
    if(neighbor->order_updated != order_updated_times){
      neighbor->set_priority_and_order(current_node->order);
      neighbor->reordering(N, D);
      while (!neighbor->search_tree.empty()) {
        delete neighbor->search_tree.front();
        neighbor->search_tree.pop();
      }
      neighbor->search_tree.push(new LNode());

      neighbor->order_updated = order_updated_times;
      propagate_order_to_neighbors(neighbor);

    }
  }
}

void Planner::build_dependence_graph(HNode* input_H_goal) {
  if(input_H_goal->h != 0){
    return;
  }
  std::vector<double> agent_cost(N, 0);
  std::vector<double> agent_ratio(N, 0);
  HNode* current = input_H_goal;
  std::vector<std::vector<uint>> solution_nodes = std::vector<std::vector<uint>>(N);
  for (uint i = 0; i < N; ++i) {
    solution_nodes[i].push_back(input_H_goal->C[i]->index);
  }

  while (current->parent != nullptr) {
    get_edge_cost_per_agent(agent_cost,current->C,current->parent->C);
    for(uint i = 0; i < N; ++i) {
      if(solution_nodes[i].size() == 1 
      && current->C[i]->index == input_H_goal->C[i]->index){
        continue; // skip if already at goal
      }
      solution_nodes[i].push_back(current->C[i]->index);
    }
    current = current->parent;
  }

  for (uint i = 0; i < N; ++i) {
    solution_nodes[i].push_back(ins->starts[i]->index);
    std::reverse(solution_nodes[i].begin(), solution_nodes[i].end());
    agent_ratio[i] = agent_cost[i] / (D.get(i, ins->starts[i]) - D.get(i, input_H_goal->C[i]));
  }

  std::vector<uint> ranking(N);
  std::iota(ranking.begin(), ranking.end(), 0);
  std::sort(ranking.begin(), ranking.end(), [&](uint i, uint j) {
      return agent_ratio[i] > agent_ratio[j];
  });

  std::map<std::tuple<int, int>, int> occupancy_map;
  // Suppose you have: std::vector<std::vector<int>> solution_nodes; // [agent][time] = vertex_index
  for (int agent_id = 0; agent_id < N; ++agent_id) {
      for (int t = 0; t < solution_nodes[agent_id].size(); ++t) {
          int vertex_index = solution_nodes[agent_id][t];
          occupancy_map[{vertex_index, t}] = agent_id;
      }
  }

  std::vector<std::set<int>> interacted_agents(N);
  for (auto agent_id : ranking) {
    // process agent_id in order of decreasing agent_ratio
    for (size_t j = 0; j < solution_nodes[agent_id].size() - 1; ++j) {
      Vertex* from_v = ins->G.U[solution_nodes[agent_id][j]];
      Vertex* to_v = ins->G.U[solution_nodes[agent_id][j + 1]];
      int current_time_step = j; // Assuming j starts from 0, so +1 for time step

      // TODO::Due to swap option, sometimes there could be no cached operation. 
      // TODO:: Skip this case for now. 
      if(action_history[agent_id][from_v->index][0] == nullptr) continue;
      for(auto vertex : action_history[agent_id][from_v->index] ){
        if(vertex->index == to_v->index) break; // skip if already at the next vertex
        auto edge = std::make_pair(vertex->index, current_time_step + 1);
        if(occupancy_map.find(edge) != occupancy_map.end()){
          if(occupancy_map[edge] != agent_id){
            interacted_agents[agent_id].insert(occupancy_map[edge]);
          }
        }
      }
    }
  }
  // std::vector<bool> visited(N, false);
  // std::vector<std::vector<int>> depth_clustered_agents;
  // for(auto rank : ranking){
  //   if(visited[rank]) continue; // skip if already visited
  //   depth_clustered_agents.push_back(depth_cluster(interacted_agents, rank, 1, visited));
  // }
  // export_interacted_agents_graph(interacted_agents, "interacted_agents.csv");
  // export_cluster_solution(solution_nodes, depth_clustered_agents, "clustered_solution.csv");
  std::vector<bool> visited(N, false);
  std::vector<int> bfs_order = {};
  for(auto rank : ranking){
    if(visited[rank]) continue;
    bfs_ordering(rank, interacted_agents, visited, bfs_order);
  }
  optimize_traffic_based_on_order(bfs_order, solution_nodes);
  // optimize_traffic_based_on_order(ranking, solution_nodes, interacted_agents);
}


void Planner::bfs_ordering(int node, const std::vector<std::set<int>>& graph, std::vector<bool>& visited, std::vector<int>& result) {
  std::queue<int> q;
  q.push(node);
  while (!q.empty()) {
      int next_node = q.front(); q.pop();
      if(!visited[next_node] ){
        visited[next_node] = true;
        result.push_back(next_node);
        for (int neighbor : graph[next_node]) {
          q.push(neighbor);
        }
      }
  }
}

void Planner::optimize_traffic_based_on_order(std::vector<uint> ranking, 
std::vector<std::vector<uint>>& solution, std::vector<std::set<int>>& interacted_agents){
  std::cout<< "Optimizing traffic based on BFS order..." << std::endl;
  traffic_map.reset();
  traffic_map.initialize_traffic_map(solution);
  // traffic_map.export_traffic_map_csv("old_traffic_map.csv");
  std::vector<std::vector<uint>> revised_path = solution;
  int times = 1; 
  while (times > 0 ){
    std::vector<bool> visited(N, false);
    std::vector<int> bfs_order = {};
    for(auto rank : ranking){
      if(visited[rank]) continue;
      bfs_ordering(rank, interacted_agents, visited, bfs_order);
    }
    for (size_t i = 0; i < bfs_order.size(); ++i) {
      int agent_id = bfs_order[i];
      traffic_map.remove_path(revised_path[agent_id]);
      auto path = astar_search.compute_traffic_path_index(ins->starts[agent_id]->index, ins->goals[agent_id]->index);
      if (path.empty()) {
        std::cout << "No path found for agent " << agent_id << std::endl;
        continue; // No path found, skip this agent
      }
      traffic_map.add_path(path);
      revised_path[agent_id] = path;
      // update the edge weights in the traffic map
    } 
    std::shuffle(ranking.begin(), ranking.end(), *MT);
    times --;
  } 
  guidance_heuristic.set_gudiance_path(revised_path);
}

void Planner::optimize_traffic_based_on_order(std::vector<int>& bfs_order, 
  std::vector<std::vector<uint>>& solution ) {
  // This function optimizes the traffic map based on the current order of agents.
  // It updates the edge weights in the traffic map according to the current order.
  // std::cout<< "Optimizing traffic based on BFS order..." << std::endl;
  traffic_map.reset();
  traffic_map.initialize_traffic_map(solution);
  // traffic_map.export_traffic_map_csv("old_traffic_map.csv");
  std::vector<std::vector<uint>> revised_path = solution;
  std::shuffle(bfs_order.begin(), bfs_order.end(), *MT);
  int times = 10; 
  while (times > 0 ){
    for (size_t i = 0; i < bfs_order.size(); ++i) {
      int agent_id = bfs_order[i];
      traffic_map.remove_path(revised_path[agent_id]);
      auto path = astar_search.compute_traffic_path_index(ins->starts[agent_id]->index, ins->goals[agent_id]->index);
      if (path.empty()) {
        std::cout << "No path found for agent " << agent_id << std::endl;
        continue; // No path found, skip this agent
      }
      traffic_map.add_path(path);
      revised_path[agent_id] = path;
      // update the edge weights in the traffic map
    } 
    // std::shuffle(bfs_order.begin(), bfs_order.end(), *MT);
    times --;
  } 
  // std::cout<< "Finished traffic based on BFS order..." << std::endl;
  // traffic_map.reset();
  // traffic_map.initialize_traffic_map(revised_path);
    // guidance_heuristic.initialized = false;
  guidance_heuristic.set_gudiance_path(revised_path);
  // export_revised_path(revised_path, "revised_path_3.csv");
  // traffic_map.export_traffic_map_csv("revised_traffic_map4.csv");
}

void Planner::export_revised_path(const std::vector<std::vector<uint>>& revised_path, const std::string& filename) {
    std::ofstream fout(filename);
    fout << "agent_id,time_step,vertex_index\n";
    for (size_t agent = 0; agent < revised_path.size(); ++agent) {
        for (size_t t = 0; t < revised_path[agent].size(); ++t) {
            fout << agent << "," << t << "," << revised_path[agent][t] << "\n";
        }
    }
    fout.close();
}


std::vector<int> Planner::depth_cluster(const std::vector<std::set<int>>& graph, int start, int input_depth, std::vector<bool>& visited) {
    std::vector<int> cluster;
    std::queue<std::pair<int, int>> q; // {node, depth}
    q.emplace(start, 0);
    visited[start] = true;

    int visited_count = 0;
    while (!q.empty()) {
        auto [node, depth] = q.front(); q.pop();
        cluster.push_back(node);
        if (depth == input_depth) continue;
        if(visited_count == 50) continue; // limit to 50 nodes
        for (int neighbor : graph[node]) {
            if (!visited[neighbor]) {
              visited[neighbor] = true;
              visited_count++;
              if(visited_count == 50){
                break;
              }
              q.emplace(neighbor, depth + 1);
            }
        }
    }
    return cluster;
}



// Computes the full transitive closure of a directed graph.
// Returns a vector where closure[i] is the set of nodes reachable from i.
std::vector<std::set<int>> Planner::transitiveClosureAll(const std::vector<std::set<int>>& graph) {
    int n = graph.size();
    std::vector<std::set<int>> closure(n);
    std::vector<bool> visited(n);

    // DFS visits all reachable from 'start', marking visited and populating closure[start].
    std::function<void(int,int)> dfs = [&](int start, int u) {
        for (int v : graph[u]) {
            if (!visited[v]) {
                visited[v] = true;
                closure[start].insert(v);
                dfs(start, v);
            }
        }
    };

    // Compute closure for each node
    for (int i = 0; i < n; ++i) {
        std::fill(visited.begin(), visited.end(), false);
        visited[i] = true;
        dfs(i, i);
    }
    return closure;
}





void Planner::export_cluster_solution(
    const std::vector<std::vector<int>>& solution_nodes,
    const std::vector<std::vector<int>>& depth_clustered_agents,
    const std::string& filename)
{
    if (depth_clustered_agents.empty()) return;
    const auto& cluster = depth_clustered_agents[0];

    std::ofstream fout(filename);
    fout << "agent_id,time_step,vertex_index\n";
    for (int agent_id : cluster) {
        for (size_t t = 0; t < solution_nodes[agent_id].size(); ++t) {
            fout << agent_id << "," << t << "," << solution_nodes[agent_id][t] << "\n";
        }
    }
    fout.close();
}

void Planner::export_interacted_agents_graph(const std::vector<std::set<int>>& interacted_agents, const std::string& filename) {
    std::ofstream fout(filename);
    fout << "src,dst\n";
    for (size_t agent = 0; agent < interacted_agents.size(); ++agent) {
        for (int other : interacted_agents[agent]) {
            fout << agent << "," << other << "\n";
        }
    }
    fout.close();
}

void Planner::backpropagate_order(HNode* input_H_goal) {
  HNode* current  = input_H_goal;
  SOLUTION_NODES = std::vector<HNode*>();
  while (current->parent != nullptr) {
    current = current->parent;
  }
  SOLUTION_NODES.push_back(current);
  return;
  // return;
  // order_updated_times ++;
  // // return;
  // std::vector<double> agent_cost(N, 0);
  // std::vector<double> agent_ratio(N, 0);
  // HNode* current = input_H_goal;
  // current->order_updated = order_updated_times;

  // // while (current->parent != nullptr) {
  // //   current->order_updated = order_updated_times;
  // //   current = current->parent;
  // // }

  // current = input_H_goal;
  // while (current->parent != nullptr) {
  //   get_edge_cost_per_agent(agent_cost,current->C,current->parent->C);
  //   SOLUTION_NODES = std::vector<HNode*>();
  //   SOLUTION_NODES.push_back(current);

  //   for (uint i = 0; i < N; ++i) {
  //     if( D.get(i,current->parent->C[i]) == 0){
  //       // a large number 
  //       agent_ratio[i] = 10000;
  //     }else{
  //       agent_ratio[i] = agent_cost[i] / (D.get(i, current->parent->C[i]) -
  //                                         D.get(i, input_H_goal->C[i]));
  //     }
  //   }

  //   // std::uniform_real_distribution<double> noise_dist(-0.1, 0.1); // Adjust range as needed
  //   // for (size_t i = 0; i < agent_ratio.size(); ++i) {
  //   //     agent_ratio[i] = agent_ratio[i] * (1 + noise_dist(*MT));
  //   // }
  //   std::sort(current->parent->order.begin(),current->parent->order.end(), [&](uint i, uint j) {
  //       return (agent_ratio[i] ) < (agent_ratio[j]);
  //   });

  //   // for(auto o :current->parent->order){
  //   //   std::cout<< o << " ";
  //   // }
  //   // std::cout<< std::endl;

  //   // fix order for partent: 
  //   current->parent->set_priority_and_order(current->parent->order);
  //   current->parent->reordering(N,D);
  //   while (!current->parent->search_tree.empty()) {
  //     delete current->parent->search_tree.front();
  //     current->parent->search_tree.pop();
  //   }
  //   current->parent->search_tree.push(new LNode());
  //   current->parent->order_updated = order_updated_times;
  //   // TODO: propagate the order to neighbourhood. 
  //   propagate_order_to_neighbors(current->parent);

  //   // if(makespan > start && makespan < end){
  //     // update the edge weights
  //     // for(int i = 0; i < N; ++i){
  //     //   if(current->parent->C[i]->index != current->C[i]->index){
  //     //     PIBT_D.increase_edge_weight(current->parent->C[i]->index, current->C[i]->index, 1);
  //     //   }
  //     // }
  //   // }

  //   current = current->parent;
  //   // makespan --;
  // }
  // SOLUTION_NODES.push_back(current);
  // PIBT_D.reset(ins);
}


void Planner::increase_visited_node(std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED) {
  // Example usage of PairHash
  std::unordered_map<std::pair<int, int>, int, PairHash> transition_frequency;

  for (auto n : EXPLORED) {
    auto current = n.second;
    for (int i = 0; i < N; ++i) {
        if(current->parent == nullptr || current->C[i] == nullptr){
          std::cout<<" OMMMMMM "<<std::endl;
        }
        if (current->parent != nullptr) {
            if (current->C[i]->index != current->parent->C[i]->index) {
                // Track the transition frequency
                std::pair<int, int> transition = {current->parent->C[i]->index, current->C[i]->index};
                transition_frequency[transition]++;
            }
        }
    }
  }

  std::vector<std::pair<std::pair<int, int>, int>> frequency_vector(transition_frequency.begin(), transition_frequency.end());

    // Sort the vector by frequency in descending order
  std::sort(frequency_vector.begin(), frequency_vector.end(),
            [](const std::pair<std::pair<int, int>, int>& a, const std::pair<std::pair<int, int>, int>& b) {
                return a.second > b.second; // Sort by frequency (value) in descending order
            });

  // Calculate the number of top elements to select (20% of the total size)
  size_t top_20_percent_count = static_cast<size_t>(frequency_vector.size() * 0.1);

  // Print the top 20% of pairs
  
  for (size_t i = 0; i < top_20_percent_count; ++i) {
    std::cout<<"Adding traffic to edge: "<< frequency_vector[i].first.first << " " << 
    frequency_vector[i].first.second << std::endl;
    const auto& pair = frequency_vector[i].first;
    PIBT_D.set_edge_weight(pair.first, pair.second, 1);
  }
  PIBT_D.reset(ins);
}

// void Planner::decay_punishment(){
//   // decay_factor = decay_factor * 0.8;
//   // std::cout<<"decay_factor: " << decay_factor << std::endl;
//   // return; 
//   PIBT_D.clear_traffic();
//   int max_size =  punishment_vector.size();
//   int decay_layer = get_random_int(MT, 0 , max_size - 1);
//   for(auto& pair: punishment_vector[decay_layer]){
//     pair.second = decay_factor * pair.second;
//   }

//   for(int i = 0; i < punishment_vector.size(); i++){
//     for(auto& pair: punishment_vector[i]){
//       PIBT_D.increase_directed_edge_weight(pair.first.first, pair.first.second,
//         pair.second);
//     }
//   }
//   PIBT_D.reset(ins);
// }


// void Planner::decay_punishment(){
//   // decay_factor = decay_factor * 0.8;
//   // std::cout<<"decay_factor: " << decay_factor << std::endl;
//   // return; 
//   PIBT_D.clear_traffic();
//   for(int i = 0; i< decay_counter; i++){
//     for(auto& pair: punishment_vector[i]){
//        PIBT_D.increase_directed_edge_weight(pair.first.first, pair.first.second,
//         pair.second);
//     }
//   }
//   for(int i = decay_counter; i < punishment_vector.size(); i++){
//     for(auto& pair: punishment_vector[i]){
//        PIBT_D.increase_directed_edge_weight(pair.first.first, pair.first.second,
//         decay_factor * pair.second);
//     }
//   }
//   decay_counter -- ;
//   if(decay_counter < 0){
//     decay_counter = punishment_vector.size() - 1;
//     decay_factor = decay_factor * 0.5;
//   }
//   PIBT_D.reset(ins);
// }
 

void Planner::decay_punishment(){
  // decay_factor = decay_factor * 0.8;
  // std::cout<<"decay_factor: " << decay_factor << std::endl;
  // return; 
  PIBT_D.clear_traffic();
  for(auto& pair: punishment_vector[punishment_vector.size()-1]){
    pair.second = 0.5 * pair.second;
  }
  for(int i = 0; i < punishment_vector.size(); i++){
    for(auto& pair: punishment_vector[i]){
       PIBT_D.increase_edge_weight(pair.first.first, pair.first.second, 
      pair.second);
    }
  }
  PIBT_D.reset(ins);
}


// void Planner::set_individual_congestion_map(HNode* H_init){
//   PIBT_D.clear_all_traffic();
//   for (int i = 0; i < N; i++) {
//     for(auto & p : individual_transition_frequency[i]){
//       if (p.first.first == p.first.second) continue;
//       double result = p.second;
//       PIBT_D.increase_edge_weight(i, p.first.first, p.first.second,result);
//     }
//   }
//   PIBT_D.reset(ins);
// }

void Planner::set_individual_congestion_map(HNode* H_init){
  PIBT_D.clear_all_traffic();
  for (int i = 0; i < N; i++) {
    std::unordered_map<std::pair<int, int>, double, PairHash> transition_frequency;
    for(int j  = 0; j < i; j++){
      for(auto & frequency :  individual_transition_frequency[H_init->order[j]]){
        transition_frequency[frequency.first] += frequency.second;
      }
    }
    for(auto& p : transition_frequency){
      if (p.first.first == p.first.second) continue;
      auto increased_weight =  p.second;
      PIBT_D.increase_edge_weight(H_init->order[i], p.first.first, p.first.second, increased_weight);
    }
  }
  PIBT_D.reset(ins);
  // // std::cout<< "start computing map" << std::endl;
  // for (int i = 0; i < N; i++) {
  //   std::unordered_map<std::pair<int, int>, double, PairHash> transition_frequency;
  //   for(int j  = 0; j < N; j++){
  //     auto agent = H_init->order[j];
  //     for(auto & frequency :  individual_transition_frequency[agent]){
  //       transition_frequency[frequency.first] += frequency.second;
  //     }
  //   }
  //   std::vector<double> vertex_congestion(ins->G.U.size(), 0);
  //   std::unordered_map<std::pair<int, int>, double, PairHash> contra_flow;
  //   for (auto p : transition_frequency) {
  //     if (p.first.first == p.first.second) {
  //       vertex_congestion[p.first.second] += p.second;
  //       continue;
  //     }
  //     if( transition_frequency.find({p.first.second, p.first.first}) != transition_frequency.end()){
  //       contra_flow[{std::min(p.first.first, p.first.second), std::max(p.first.first, p.first.second)}] 
  //       = p.second * transition_frequency[{p.first.second,p.first.first}];
  //     }
  //     // Update vertex congestion
  //     vertex_congestion[p.first.second] += p.second;
  //   } 

  //   for(unsigned int e = 0; e  < vertex_congestion.size(); e++){
  //     if(vertex_congestion[e] > 0){
  //       for(auto f : ins->G.U[e]->neighbor){
  //         contra_flow[{std::min(f->index,e), std::max(f->index,e)}] += vertex_congestion[e];  
  //       }
  //     }
  //   }

  //   for(auto& p : contra_flow){
  //     auto increased_weight =  p.second;
  //     PIBT_D.increase_edge_weight(H_init->order[i],p.first.first, p.first.second, increased_weight);
  //   }
  // }
  // std::cout<< "finish computing map" << std::endl;
  // PIBT_D.reset(ins);
}


void Planner::record_each_agent_frequency(HNode* input_H_goal){
  HNode* current = input_H_goal;
  while (current->parent != nullptr) {
    for (int i = 0; i < N; ++i) {
      std::pair<int, int> transition = {current->parent->C[i]->index, current->C[i]->index};
      individual_transition_frequency[i][transition]++;
    }
    current = current->parent;
  }
  reset_congestion_map = true;
}

void Planner::learning_Q_value(HNode* input_H_goal, double is_goal){
  struct PairHash {
    std::size_t operator()(const std::pair<Vertex*, Vertex*>& p) const {
        return std::hash<int>()(p.first->id) ^ (std::hash<int>()(p.second->id) << 1);
    }
  };
  HNode* current = input_H_goal;
  // For each agent, track the cost-to-go (steps to goal)
  std::vector<double> cost_to_go(N, 0.0);
  std::vector<std::unordered_map<std::pair<Vertex*, Vertex*>, int, PairHash>> edge_mapper(N);
  if(!is_goal){
    // use heuristic to fill the value if input is not goal.
    for (int i = 0; i < N; ++i) {
      Vertex* v_curr = current->parent->C[i];
      Vertex* v_next = current->C[i];
      cost_to_go[i] = get_q_value(i, v_curr, v_next);
    }
  }

  std::vector<HNode*> solution_nodes;
  std::vector<std::vector<Vertex*>> each_agent_solution_nodes(N);
  // s -> a -> b -> t,  (s -> a, 2)
  // First, compute the path length for each agent, we'll backtrack and increment cost_to_go for each step
  while (current->parent != nullptr) {
    for (int i = 0; i < N; ++i) {
        if(current->parent->C[i]->index == ins->goals[i]->index && 
           cost_to_go[i] == 0){
          // skip if stay at goal; 
          continue;
        }
        Vertex* v_curr = current->parent->C[i];
        Vertex* v_next = current->C[i];
        cost_to_go[i] += 1.0;
        each_agent_solution_nodes[i].push_back(v_curr);
        // Store the edge and its cost-to-go
        if(edge_mapper[i].find({v_curr, v_next}) == edge_mapper[i].end()) {
          edge_mapper[i][{v_curr, v_next}] = cost_to_go[i];
        }
    }
    solution_nodes.push_back(current);
    current = current->parent;
  }

  for (int i = 0; i < N; ++i) {
    for (auto& p : edge_mapper[i]) {
      // For each edge, update the Q-table.
      if(p.first.first->index == p.first.second->index){
        if(p.first.first->index == ins->goals[i]->index){
          // skip if stay at goal; 
          continue;
        }
      }
      // try propogate Q value k steps
      propogate_q_value_k_steps(i, p.first.first, p.first.second, p.second, 5);
    }
  }
  // for (int i = 0; i < N; ++i) {
  //   // Q_tables[0].print_Q_value_at_state(ins->starts[i]->id);
  //   // std::vector< std::vector<std::tuple<State, int, double>> > q_records(V_size*2);
  //   // std::vector< std::vector<std::tuple<State, int, double>> > updated_q_records(V_size*2);
  //   // for(auto& a : each_agent_solution_nodes[i]){
  //   //   // std::cout<<"Adding solution node for agent: " << i << " at vertex: " << a << std::endl;
  //   //   // Q_tables[i].add_solution_node(a);
  //   //   q_records[a->index] = Q_tables[i].get_q_value_records(a->index);
  //   // }
  //   // Q_tables[i].dump_q_table_csv("q_table_agent_" + std::to_string(i) + ".csv");
  //   for (auto& p : edge_mapper[i]) {
  //     // For each edge, update the Q-table.
  //     if(p.first.first->index == p.first.second->index){
  //       if(p.first.first->index == ins->goals[i]->index){
  //         // skip if stay at goal; 
  //         continue;
  //       }
  //     }
  //     // std::cout<<"Updating Q value for agent: " << i << " from "<< p.first.first->index 
  //             //  << " to " << p.first.second->index << " with cost: " << p.second << std::endl;
  //     // try propogate Q value k steps
  //     propogate_q_value_k_steps(i, p.first.first, p.first.second, p.second, 10);
  //   }
  //   // for(auto& a : each_agent_solution_nodes[i]){
  //   //   // std::cout<<"Adding solution node for agent: " << i << " at vertex: " << a << std::endl;
  //   //   // Q_tables[i].add_solution_node(a);
  //   //   updated_q_records[a->index] = Q_tables[i].get_q_value_records(a->index);
  //   // }
  //   // Q_tables[i].dump_q_table_csv("updated_q_table_agent_" + std::to_string(i) + ".csv");
  //   // for(int j = 0; j < V_size*2; j++){
  //   //   if(q_records[j].size() > 0 && updated_q_records[j].size() > 0){
  //   //     for(auto tuple : q_records[j]){
  //   //       std::cout<<"Q value for agent: " << 0 << " at vertex: " << j << " is: " 
  //   //             <<"action: " << std::get<1>(tuple) 
  //   //             << " and value: " << std::get<2>(tuple)
  //   //       << std::endl;
  //   //     }
  //   //     for(auto tuple : updated_q_records[j]){
  //   //       std::cout<<"Q value for agent: " << 0 << " at vertex: " << j << " is: " 
  //   //             <<"action: " << std::get<1>(tuple) 
  //   //             << " and value: " << std::get<2>(tuple)
  //   //       << std::endl;
  //   //     };
  //   //     Q_tables[i].print_Q_value_at_state(j);
  //   //     bool a = 0;
  //   //   }
  //   // bool b = 0 ;
  //   // }
  // }
  
}
void Planner::increase_traffic_based_on_solution(HNode* input_H_goal) {
    // std::cout<< "Increasing traffic map based on solution" << std::endl;
    guidance_heuristic.initialized = true;
    traffic_map.reset();
    HNode* current = input_H_goal;
    std::vector<std::vector<uint>> solution_nodes = std::vector<std::vector<uint>>(N);
    for (uint i = 0; i < N; ++i) {
      solution_nodes[i].push_back(input_H_goal->C[i]->index);
    }
    while (current->parent != nullptr) {
      for(uint i = 0; i < N; ++i) {
        if(solution_nodes[i].size() == 1 
        && current->C[i]->index == input_H_goal->C[i]->index){
          continue; // skip if already at goal
        }
        solution_nodes[i].push_back(current->C[i]->index);
      }
      current = current->parent;
    }
    for (uint i = 0; i < N; ++i) {
      solution_nodes[i].push_back(ins->starts[i]->index);
      std::reverse(solution_nodes[i].begin(), solution_nodes[i].end());
    }
    for(auto solution : solution_nodes){
      traffic_map.add_incremental_flow_path(solution); 
    }
    traffic_map.record_incremental_flow();
    

    // traffic_map.print_incremental_flow("incremental_flow.csv");
    guidance_heuristic.reset(ins);
    // std::cout<<"finishing increasing traffic map" << std::endl;
}


void Planner::increase_time_dependent_traffic_based_on_solution(HNode* input_H_goal) {
    guidance_heuristic.initialized = true;
    for(uint t = 0; t < time_bucket_size; t++){
      time_period_traffic_map[t].reset();
    }
    int solution_makespan =  input_H_goal->current_make_span +1 ;
    HNode* current = input_H_goal;
    std::vector<std::vector<uint>> solution_nodes = std::vector<std::vector<uint>>(N);
    for (uint i = 0; i < N; ++i) {
      solution_nodes[i].push_back(input_H_goal->C[i]->index);
    }
    while (current->parent != nullptr) {
      for(uint i = 0; i < N; ++i) {
        if(solution_nodes[i].size() == 1 
        && current->C[i]->index == input_H_goal->C[i]->index){
          continue; // skip if already at goal
        }
        solution_nodes[i].push_back(current->C[i]->index);
      }
      current = current->parent;
    }
    for (uint i = 0; i < N; ++i) {
      solution_nodes[i].push_back(ins->starts[i]->index);
      std::reverse(solution_nodes[i].begin(), solution_nodes[i].end());
    }
    if(best_makespan == 0){
      best_makespan = solution_makespan * 2;
      // max_time_period = best_makespan/time_bucket_size;
      // max_time_period = static_cast<int>(std::ceil(static_cast<double>(best_makespan) / time_bucket_size));
      max_time_period = static_cast<int>(std::ceil(static_cast<double>(best_makespan) / time_bucket_size)); 
      if(max_time_period < 1){
        max_time_period = 1;
      }
      // std::cout<<"Initial max time period: " << best_makespan << std::endl;
    }

    for(auto solution : solution_nodes){
      for(uint t = 0; t < time_bucket_size; t++){
        uint start_index = t * max_time_period;
        time_period_traffic_map[t].add_incremental_flow_path_from_time_index(solution,  start_index);
      }
    }
    for(uint t = 0; t < time_bucket_size; t++){
      time_period_traffic_map[t].record_incremental_flow();
    }
    
    current_time_bucket = 0;
    // traffic_map.print_incremental_flow("incremental_flow.csv");
    guidance_heuristic.reset(ins);
    guidance_heuristic.set_traffic_map(&time_period_traffic_map[current_time_bucket]);
    // std::cout<<"finishing increasing traffic map" << std::endl;
}

void Planner::increase_weight_map(HNode* input_H_goal, bool is_goal){
  
  // build_dependence_graph(input_H_goal);


  
  if(is_goal){
    // std::cout<< input_H_goal->current_make_span<< std::endl; 
    // std::cout<<"Increasing weight map based on solution" << std::endl;
    // traffic_map.print_incremental_flow("flow" + std::to_string(solution_count) + ".csv");
    export_solution_from_HNode(input_H_goal, "solution_time_dependent" + std::to_string(solution_count) + ".csv");
    // std::cout<<"Exporting solution to solution" + std::to_string(solution_count) + ".csv" << std::endl;
    solution_count++;
  }else{
    // std::cout<<"Increasing weight map based on partial solution" << std::endl;
  }
  increase_time_dependent_traffic_based_on_solution(input_H_goal);
  // increase_traffic_based_on_solution(input_H_goal);
  
  // learning_Q_value(input_H_goal,is_goal);

  // record_each_agent_frequency(input_H_goal);
  // // // return; 
  // if(is_goal){
  // //   increase_each_agent_cost(input_H_goal);
  // //   // increase_solution_weight(input_H_goal);
    // increase_solution_congestion_cost(input_H_goal);
  // }
    // std::cout<<"Finished increasing weight map" << std::endl;
  // std::cout<<"Function called" << std::endl;
  // if(is_goal){
  //   PIBT_D.copy_global_data_and_clean_local();
  // }
  // increase_solution_congestion_cost(input_H_goal);
  // increase_solution_weight(input_H_goal);
}



void Planner::increase_each_agent_cost(HNode* input_H_goal) {
  std::vector<std::vector<int>> config_c(N); 
  std::vector<std::vector<std::pair<int, double>>> learned_travel_cost(N); 
  // Example usage of PairHash
  HNode* current = input_H_goal;
  for( int i = 0; i < N; i++){
    learned_travel_cost[i].push_back({ins->goals[i]->index, 0});
  } 

  // find cost to go; 
  while (current->parent != nullptr) {
    for (int i = 0; i < N; ++i) {
      config_c[i].push_back(current->C[i]->index);
      if(current->parent->C[i]->index == ins->goals[i]->index && 
        learned_travel_cost[i].back().second == 0){
        // skip if stay at goal; 
        continue;
      }
      if( learned_travel_cost[i].back().first == current->parent->C[i]->index){
        // skip if stay at goal; 
        learned_travel_cost[i].back().second += 1;  
      }else{
        learned_travel_cost[i].push_back({current->parent->C[i]->index, learned_travel_cost[i].back().second + 1});
      }
    }
    current = current->parent;
  }

  for( int i = 0; i < N; i++){
    for(int j = 0; j < learned_travel_cost[i].size(); j++){
      if(learned_travel_cost[i][j].first == ins->goals[i]->index){
        continue;
      }
      auto increased_weight =  learned_travel_cost[i][j].second;
      PIBT_D.update_heuristic_table(ins, i, learned_travel_cost[i][j].first, increased_weight);
    }
  } 
  PIBT_D.reset(ins);
}



// void Planner::increase_each_agent_cost(HNode* input_H_goal) {
//   std::vector<std::pair<int, double>> learned_travel_cost(N); 
//   // Example usage of PairHash
//   HNode* current = input_H_goal;
//   for( int i = 0; i < N; i++){
//     learned_travel_cost[i].first = current->C[i]->index;
//     learned_travel_cost[i].second = 0;
//   } 
//   while (current->parent != nullptr) {
//     for (int i = 0; i < N; ++i) {
//       if(current->parent->C[i]->index == ins->goals[i]->index){
//         // skip if stay at goal; 
//         continue;
//       }
//       std::pair<int, int> transition = {current->parent->C[i]->index, current->C[i]->index};
//       individual_transition_frequency[i][transition]++;
//     }
//     current = current->parent;
//   }
// }

// void Planner::increase_each_agent_cost(HNode* input_H_goal) {
//   // Example usage of PairHash
//   HNode* current = input_H_goal;
//   while (current->parent != nullptr) {
//     for (int i = 0; i < N; ++i) {
//       std::pair<int, int> transition = {current->parent->C[i]->index, current->C[i]->index};
//       individual_transition_frequency[i][transition]++;
//     }
//     current = current->parent;
//   }
// }

void Planner::increase_solution_congestion_cost(HNode* input_H_goal) {
  // Example usage of PairHash
  // double decay = 0.8; 
  // PIBT_D.decay_local_entries( decay);
  std::unordered_map<std::pair<int, int>, double, PairHash> transition_frequency;
  int max_frequent = 0 ;
  HNode* current = input_H_goal;
  while (current->parent != nullptr) {
    for (int i = 0; i < N; ++i) {
      std::pair<int, int> transition = {current->parent->C[i]->index, current->C[i]->index};
      transition_frequency[transition]++;
    }
    current = current->parent;
  }
  std::vector<double> vertex_congestion(ins->G.U.size(), 0);
  std::unordered_map<std::pair<int, int>, double, PairHash> contra_flow;
  for (auto p : transition_frequency) {
    if (p.first.first == p.first.second) {
      vertex_congestion[p.first.second] += p.second;
      continue;
    }
    if( transition_frequency.find({p.first.second, p.first.first}) != transition_frequency.end()){
      contra_flow[{std::min(p.first.first, p.first.second), std::max(p.first.first, p.first.second)}] 
      = p.second * transition_frequency[{p.first.second,p.first.first}];
    }
    // Update vertex congestion
    vertex_congestion[p.first.second] += p.second;
  } 

  for(unsigned int i = 0; i < vertex_congestion.size(); i++){
    if(vertex_congestion[i] > 0){
      for(auto j : ins->G.U[i]->neighbor){
        contra_flow[{std::min(j->index,i), std::max(j->index,i)}] += vertex_congestion[i];  
      }
    }
  }
  for(auto& p : contra_flow){
    auto increased_weight =  p.second;
    // PIBT_D.increase_edge_weight(p.first.first, p.first.second, increased_weight);
    PIBT_D.increase_edge_weight(p.first.first, p.first.second, 
      increased_weight);
  }

  for(auto& p : contra_flow){
    auto traffic = traffic_map.get_incremental_traffic_cost(p.first.first, p.first.second);
    auto pibt_traffic = PIBT_D.get_edge_weight(p.first.first, p.first.second);
    if(traffic != pibt_traffic){
      std::cout<<"Edge: " << p.first.first << " " << p.first.second << std::endl;
      std::cout<<"Mismatch traffic: " << traffic << " vs " << pibt_traffic << std::endl;
    }
    // PIBT_D.increase_edge_weight(p.first.first, p.first.second, increased_weight);
    // PIBT_D.increase_edge_weight(p.first.first, p.first.second, 
      // increased_weight);
  }

  // PIBT_D.reset(ins);
  // for( int i = 0; i < N; i++){
  //   PIBT_D.test_dijkstra(i,ins);
  // }
  // PIBT_D.sum_global_entries();
  // std::cout<<"Starting printing edge map" << std::endl;
  // PIBT_D.print_edge_map();
  // PIBT_D.apply_gaussian_filter();
  // PIBT_D.apply_gaussian_filter(1, &ins->G);
  PIBT_D.reset(ins);
}



void Planner::increase_solution_weight(HNode* input_H_goal) {
  // // return;
  // Example usage of PairHash
  decay_factor = 0.8;
  std::unordered_map<std::pair<int, int>, double, PairHash> transition_frequency;
  int max_frequent = 0 ;
  HNode* current = input_H_goal;
  while (current->parent != nullptr) {
    for (int i = 0; i < N; ++i) {
      std::pair<int, int> transition = {std::min(current->parent->C[i]->index, current->C[i]->index),
                                        std::max(current->parent->C[i]->index, current->C[i]->index)};
      transition_frequency[transition]++;
    }
    current = current->parent;
  }
  std::vector<std::pair<std::pair<int, int>, double>> traffic_vector; 
  for(auto & p : transition_frequency){
    if (p.first.first == p.first.second) {
      for(auto j : ins->G.U[p.first.first]->neighbor){
        PIBT_D.increase_edge_weight(p.first.first, j->index, p.second);
      }
    }else{
      PIBT_D.increase_edge_weight(p.first.first, p.first.second,p.second);  
    }
  }
  // punishment_vector.push_back(traffic_vector);
  // decay_counter = punishment_vector.size() - 1;
  PIBT_D.reset(ins);
}


void Planner::add_punishment(HNode* input_H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED) {
  // Example usage of PairHash
  std::unordered_map<std::pair<int, int>, int, PairHash> transition_frequency;

  for (auto n : EXPLORED) {
    auto current = n.second;
    for (int i = 0; i < N; ++i) {
        if (current->parent != nullptr) {
            if (current->C[i]->index != current->parent->C[i]->index) {
                // Track the transition frequency
                std::pair<int, int> transition = {current->parent->C[i]->index, current->C[i]->index};
                transition_frequency[transition]++;
            }
        }
    }
  }

  std::vector<std::pair<std::pair<int, int>, int>> frequency_vector(transition_frequency.begin(), transition_frequency.end());

    // Sort the vector by frequency in descending order
  std::sort(frequency_vector.begin(), frequency_vector.end(),
            [](const std::pair<std::pair<int, int>, int>& a, const std::pair<std::pair<int, int>, int>& b) {
                return a.second > b.second; // Sort by frequency (value) in descending order
            });

  // Calculate the number of top elements to select (20% of the total size)
  size_t top_20_percent_count = static_cast<size_t>(frequency_vector.size() * 0.3);

  // Print the top 20% of pairs

  for (size_t i = 0; i < top_20_percent_count; ++i) {
    const auto& pair = frequency_vector[i].first;
    PIBT_D.set_edge_weight(pair.first, pair.second, 1);
  }
  PIBT_D.reset(ins);
}





void Planner::pick_restart_nodes(std::stack<HNode*>& OPEN){
  std::uniform_int_distribution<size_t> dist(0, SOLUTION_NODES.size() - 1);
  size_t random_index = dist(*MT);
  OPEN = std::stack<HNode*>();
  // OPEN.push(SOLUTION_NODES[random_index]);
  for(auto n : SOLUTION_NODES){
    if( n->parent == nullptr){
      OPEN.push(n);
    }
  }
  // std::cout<< "Restarting search with " << OPEN.size() << " nodes" << std::endl;
  // bool a = 0 ;
}



void Planner::export_solution_from_HNode(HNode* goal, const std::string& filename) {
    // Backtrack from goal to root, collecting configurations
    std::vector<Config> solution;
    HNode* current = goal;
    while (current != nullptr) {
        solution.push_back(current->C);
        current = current->parent;
    }
    std::reverse(solution.begin(), solution.end());

    // Export as CSV: agent_id, time_step, vertex_index
    std::ofstream fout(filename);
    fout << "agent_id,time_step,vertex_index\n";
    for (size_t t = 0; t < solution.size(); ++t) {
        for (size_t agent = 0; agent < solution[t].size(); ++agent) {
            fout << agent << "," << t << "," << solution[t][agent]->index << "\n";
        }
    }
    fout.close();
}

Solution Planner::backpropagate_solve(std::string& additional_info)
{
  solution_count = 0;
  order_updated_times = 0 ;
  best_makespan = 0;
  node_visit_times = 1;
  uint node_id = 1;
  // setup agents
  for (auto i = 0; i < N; ++i) A[i] = new Agent(i);
  // use action history to record the actions taken by each agent
  action_history = std::vector<std::vector<std::array<Vertex*, 5>>>(
      N, std::vector<std::array<Vertex*, 5>>(
          ins->G.U.size(), std::array<Vertex*, 5>{}
      )
  );

  guidance_heuristic.setup(ins);

  // setup search
  auto OPEN = std::stack<HNode*>();
  auto EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();
  SOLUTION_NODES = std::vector<HNode*>();
  // insert initial node, 'H': high-level node
  auto H_init = new HNode(ins->starts, D, nullptr, 0, get_h_value(ins->starts));
  H_init->setMakeSpan(0);
  punishment_vector.clear();
  // highlevel_node_id ++; 
  // std::cout<< *H_init << std::endl;

  OPEN.push(H_init);
  EXPLORED[H_init->C] = H_init;
  H_init->setNodeID(node_id);
  node_id++;
  std::vector<Config> solution;
  auto C_new = Config(N, nullptr);  // for new configuration
  HNode* H_goal = nullptr;          // to store goal node

  if(verbose == 3){
    std::cout<< "version: 1.4.0\n";
    std::cout<< "events:\n";
  }

  uint restart_cnt = 0;
  while (!OPEN.empty() && !is_expired(deadline)) {
    loop_cnt += 1;
    // if(OPEN.size() == 1){
    //   // set_individual_congestion_map(H_init);
    //   if(reset_congestion_map){
    //     set_individual_congestion_map(H_init);
    //     reset_congestion_map = false;
    //   }
    //   node_visit_times += 1; 
    // } 
    // do not pop here!
    auto H = OPEN.top();  // high-level node

    // low-level search end
    if (H->search_tree.empty()) {
      OPEN.pop();
      std::cout<< "Popping the nodes" << std::endl;
      continue;
    }

    if (OPEN.size() == 1) {
      // set_individual_congestion_map(H_init);
      node_visit_times += 1; 
    }

    // check lower bounds
    // if (H_goal != nullptr && H->f >= H_goal->f) {
    //   OPEN.pop();
    //   continue;
    // }

    if(verbose == 3){
      std::cout<< " - type: expanding" <<std::endl;
      std::cout<< "   id: "<< H->node_id <<std::endl;
      std::cout<< "   pId: " 
            << (H->parent == nullptr ? "0" : std::to_string(H->parent->node_id)) <<std::endl;
      std::cout<< "   f_value: "<< H->f<<std::endl;
    }

    // if (H_goal != nullptr && is_same_config(H->C, ins->goals)) {
    //   // if( H->f >= H_goal->f){
    //   //   increase_weight_map(H,false);
    //   //   pick_restart_nodes(OPEN);
    //   //   // std::cout<< "restarting search" << std::endl;
    //   //   continue;
    //   // }
    //   backpropagate_order(H);
    //   increase_weight_map(H,true);
    //   pick_restart_nodes(OPEN);
    //   // solver_info(1, "main search found solution, cost: ", H->g);
    //   // export_solution_from_HNode(H, "solution_0.csv");
    //   continue;
    // }

    if(H_goal != nullptr  && H->f >= H_goal->f){
      // export_solution_from_HNode(H, "solution_" + std::to_string(restart_cnt) + ".csv");
      restart_cnt++;
      backpropagate_order(H);
      increase_weight_map(H,false);
      pick_restart_nodes(OPEN);
      continue;
    }

    // std::cout<< "   asda   " << std::endl;

    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      H_goal = H;
      solver_info(1, "main search found solution, cost: ", H->g);
      if (objective == OBJ_NONE) break;
      backpropagate_order(H_goal);
      increase_weight_map(H_goal,true);
      pick_restart_nodes(OPEN);
      // export_solution_from_HNode(H, "current_solution.csv");
      continue;
    }
    if(guidance_heuristic.initialized != false){
      if(time_bucket_size != 1){
        uint traffic_index = H->current_make_span  / max_time_period  ;
        if( traffic_index != current_time_bucket){
          // set traffic based on time period;
          current_time_bucket = traffic_index;
          if(current_time_bucket > time_bucket_size -1){
            std::cout<< current_time_bucket << " " << time_bucket_size << std::endl;
            std::cout<< "Warning: current time bucket exceeds limit, setting to max" << std::endl;
          }
          guidance_heuristic.set_traffic_map(&time_period_traffic_map[current_time_bucket]);
          guidance_heuristic.reset(ins);
          // std::cout<< "Updating traffic map to time bucket: " << current_time_bucket << std::endl;
        }
      }
    }
    // create successors at the low-level search
    auto L = H->search_tree.front();
    H->search_tree.pop();
    expand_lowlevel_tree(H, L);


    // create successors at the high-level search
    const auto res = get_new_config(H, L);
    delete L;  // free
    if (!res) continue;

    // create new configuration
    for (auto a : A) C_new[a->id] = a->v_next;

    // check explored list
    const auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      // case found
      // rewrite(H, iter->second, H_goal, OPEN);
      
      rewrite_backpropagate(H, iter->second, H_goal, OPEN, H_init,EXPLORED);
      if(OPEN.size() == 1 ){
        continue;
      }
      // re-insert or random-restart
      auto H_insert = (MT != nullptr && get_random_float(MT) >= RESTART_RATE)
                          ? iter->second
                          : H_init;

      
      // if(H_insert == H_init){
      //   backpropagate_order(H_insert);
      // }
      
      H_insert->set_visit_times(node_visit_times);
      OPEN.push(H_insert);
      // if (H_goal == nullptr || H_insert->f < H_goal->f){
      //   // backpropagate_order(H_insert);
      //   OPEN.push(H_insert);
      //   // record_highlevel_node(H_insert,false,false);
      // }
    } else {
      // insert new search node
      const auto H_new = new HNode(
          C_new, D, H, H->g + get_edge_cost(H->C, C_new), get_h_value(C_new));
      H_new->set_visit_times(node_visit_times);
      EXPLORED[H_new->C] = H_new;
      H_new->setNodeID(node_id);
      node_id++;
      H_new->set_visit_times(node_visit_times);
      H_new->setMakeSpan(H->current_make_span + 1);
      // if (H_goal == nullptr || H_new->f < H_goal->f) 
      OPEN.push(H_new);
      if(verbose == 3){
        std::cout<< " - type: generating" <<std::endl;
        std::cout<< "   id: "<< H_new->node_id <<std::endl;
        std::cout<< "   pId: " 
              << (H_new->parent == nullptr ? "0" : std::to_string(H_new->parent->node_id)) <<std::endl;
        std::cout<< "   f_value: "<< H_new->f<<std::endl;
      }
    }
  }
  
  // backtrack
  if (H_goal != nullptr) {
    auto H = H_goal;
    while (H != nullptr) {
      solution.push_back(H->C);
      H = H->parent;
    }
    std::reverse(solution.begin(), solution.end());
  }


  // print result
  if (H_goal != nullptr && OPEN.empty()) {
    solver_info(1, "solved optimally, objective: ", objective);
  } else if (H_goal != nullptr) {
    solver_info(1, "solved sub-optimally, objective: ", objective);
  } else if (OPEN.empty()) {
    solver_info(1, "no solution");
  } else {
    solver_info(1, "timeout");
  }

  // logging
  additional_info +=
      "optimal=" + std::to_string(H_goal != nullptr && OPEN.empty()) + "\n";
  additional_info += "objective=" + std::to_string(objective) + "\n";
  additional_info += "loop_cnt=" + std::to_string(loop_cnt) + "\n";
  additional_info += "num_node_gen=" + std::to_string(EXPLORED.size()) + "\n";


  // save to tree file 
  if(tree_file != "none"){
    saveTree(tree_file);
  }

  // std::cout<< "version: 1.4.0\n";
  // std::cout<< "events:\n";
  // for(auto e : EXPLORED){
  //   std::cout<< " - type: expanding" <<std::endl;
  //   std::cout<< "   id: "<< e.second->node_id <<std::endl;
  //   std::cout<< "   pId: " 
  //         << (e.second->parent == nullptr ? "null" : std::to_string(e.second->parent->node_id)) <<std::endl;
  //   std::cout<< "   f_value: "<< e.second->f<<std::endl;
  // }

  // memory management
  for (auto a : A) delete a;
  for (auto itr : EXPLORED) delete itr.second;

  return solution;
}





Solution Planner::solve(std::string& additional_info)
{
  uint highlevel_node_id = 0;
  solver_info(1, "start search");

  // setup agents
  for (auto i = 0; i < N; ++i) A[i] = new Agent(i);

  // setup search
  auto OPEN = std::stack<HNode*>();
  auto EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();
  // insert initial node, 'H': high-level node
  auto H_init = new HNode(ins->starts, D, nullptr, 0, get_h_value(ins->starts));
  H_init->setNodeID(highlevel_node_id);
  highlevel_node_id++;
  record_highlevel_node(H_init,false,false);
  // highlevel_node_id ++; 
  // std::cout<< *H_init << std::endl;

  OPEN.push(H_init);
  EXPLORED[H_init->C] = H_init;

  std::vector<Config> solution;
  auto C_new = Config(N, nullptr);  // for new configuration
  HNode* H_goal = nullptr;          // to store goal node


  while (!OPEN.empty() && !is_expired(deadline)) {
    loop_cnt += 1;

    // do not pop here!
    auto H = OPEN.top();  // high-level node
    // std::cout<< *H << std::endl;
    // High_node_info.push_back(Search_Node_Info(*H,true,false));
    record_highlevel_node(H,true,false);
    // low-level search end
    if (H->search_tree.empty()) {
      OPEN.pop();
      continue;
    }

    // check lower bounds
    if (H_goal != nullptr && H->f >= H_goal->f) {
      OPEN.pop();
      continue;
    }

    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      H_goal = H;
      record_highlevel_node(H, true, true);
      solver_info(1, "main search found solution, cost: ", H->g);
      if (objective == OBJ_NONE) break;
      continue;
    }

    // create successors at the low-level search
    auto L = H->search_tree.front();
    H->search_tree.pop();
    expand_lowlevel_tree(H, L);

    // create successors at the high-level search
    const auto res = get_new_config(H, L);
    delete L;  // free
    if (!res) continue;

    // create new configuration
    for (auto a : A) C_new[a->id] = a->v_next;

    // check explored list
    const auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      // case found
      rewrite(H, iter->second, H_goal, OPEN);
      // re-insert or random-restart
      auto H_insert = (MT != nullptr && get_random_float(MT) >= RESTART_RATE)
                          ? iter->second
                          : H_init;
      if (H_goal == nullptr || H_insert->f < H_goal->f){
        OPEN.push(H_insert);
        record_highlevel_node(H_insert,false,false);
      } 
    } else {
      // insert new search node
      const auto H_new = new HNode(
          C_new, D, H, H->g + get_edge_cost(H->C, C_new), get_h_value(C_new));
      H_new ->setNodeID(highlevel_node_id);
      highlevel_node_id ++; 
      record_highlevel_node(H_new,false,false);
      EXPLORED[H_new->C] = H_new;
      if (H_goal == nullptr || H_new->f < H_goal->f) OPEN.push(H_new);
    }
  }
  
  // backtrack
  if (H_goal != nullptr) {
    auto H = H_goal;
    while (H != nullptr) {
      solution.push_back(H->C);
      H = H->parent;
    }
    std::reverse(solution.begin(), solution.end());
  }


  // print result
  if (H_goal != nullptr && OPEN.empty()) {
    solver_info(1, "solved optimally, objective: ", objective);
  } else if (H_goal != nullptr) {
    solver_info(1, "solved sub-optimally, objective: ", objective);
  } else if (OPEN.empty()) {
    solver_info(1, "no solution");
  } else {
    solver_info(1, "timeout");
  }

  // logging
  additional_info +=
      "optimal=" + std::to_string(H_goal != nullptr && OPEN.empty()) + "\n";
  additional_info += "objective=" + std::to_string(objective) + "\n";
  additional_info += "loop_cnt=" + std::to_string(loop_cnt) + "\n";
  additional_info += "num_node_gen=" + std::to_string(EXPLORED.size()) + "\n";


  // save to tree file 
  if(tree_file != "none"){
    saveTree(tree_file);
  }

  // memory management
  for (auto a : A) delete a;
  for (auto itr : EXPLORED) delete itr.second;

  return solution;
}

void Planner::update_ordering(HNode* h_from, HNode* h_to, size_t N, DistTable& D) {
  // dynamic priorities, akin to PIBT
  for (size_t i = 0; i < N; ++i) {
    if (D.get(i, h_to->C[i]) != 0) {
      h_to->priorities[i] = h_from->priorities[i] + 1;
    } else {
      h_to->priorities[i] = h_from->priorities[i] - (int)h_from->priorities[i];
    }
  }
  // set order
  std::iota(h_to->order.begin(), h_to->order.end(), 0);
  std::sort(h_to->order.begin(), h_to->order.end(),
        [&](uint i, uint j) { return h_to->priorities[i] > h_to->priorities[j]; });
}

void Planner::MCT_backpropagate(MCTNode* mct_node, int sample_times,  double reward, int failuare_times){
  mct_node->visits += sample_times;
  mct_node->failuare_times += failuare_times;
  mct_node->reward += reward;
  if(mct_node->parent != nullptr )
  {
    MCT_backpropagate(mct_node->parent, sample_times, reward, failuare_times);
  }
}




MCTNode* Planner::MCT_selection(std::vector<MCTNode*>& node_pool, HNode* goal_node){
  MCTNode* selected_node = nullptr; 
  double max_utc_value = -std::numeric_limits<double>::infinity();
  for(auto node : node_pool){
    if(!node->completed_node){
      // std::cout<<"Node ID: "<< node->node_id <<" ";
      double UCT_value = compute_uct_value(node);
      if(max_utc_value <= UCT_value){
        selected_node = node; 
        max_utc_value = UCT_value;
      }
    }
  }
  // std::cout<<" "<< std::endl; 
  return selected_node;
}



MCTNode* Planner::MCT_selection(std::vector<MCTNode*>& node_pool){
  MCTNode* selected_node = nullptr; 
  // std::cout<<"Pool size: "<< node_pool.size() <<std::endl;
  double max_utc_value = -std::numeric_limits<double>::infinity();
  for(auto node : node_pool){
    if(!node->completed_node){
      // std::cout<<"Node ID: "<< node->node_id <<" ";
      double UCT_value = compute_uct_value(node);
      if(max_utc_value <= UCT_value){
        selected_node = node; 
        max_utc_value = UCT_value;
      }
    }
  }
  // std::cout<<" "<< std::endl; 
  return selected_node;
}

double Planner::compute_uct_value(MCTNode* mcts_node){
  double uct_value = 0;
  if(mcts_node->visits == 0){
    uct_value = std::numeric_limits<double>::infinity();
  } else {
    double parent_visit  =   mcts_node->parent == nullptr ? mcts_node->visits : mcts_node->parent->visits;
    //we use f-value lower the better.
    // set failuare cases with maximal_f_value;
    double exploitation = 1 - minMaxStats.normalize( (mcts_node->reward + mcts_node->failuare_times * minMaxStats.get_maximum() ) / mcts_node->visits); 
    double exploration =  2 * C  * std::sqrt(std::log(parent_visit) / mcts_node->visits);
    uct_value = exploitation + exploration;
  }
  return uct_value;
}



void Planner::print_utc_value(std::vector<MCTNode*>& node_pool){
  std::cout<< "   node_selection:" <<std::endl;
  for(auto mcts_node : node_pool){
    double uct_value = 0;
    double exploitation = 0;
    double exploration = 0;
    if(mcts_node->visits == 0){
      uct_value = std::numeric_limits<double>::infinity();
      exploitation = std::numeric_limits<double>::infinity();
      exploration = std::numeric_limits<double>::infinity();
    } else {
      double parent_visit  =   mcts_node->parent == nullptr ? mcts_node->visits : mcts_node->parent->visits;
      //we use f-value lower the better.
      // set failuare cases with maximal_f_value;
      exploitation = 1 - minMaxStats.normalize( (mcts_node->reward + mcts_node->failuare_times * minMaxStats.get_maximum() ) / mcts_node->visits); 
      exploration = C  * std::sqrt(std::log(parent_visit) / mcts_node->visits);
      uct_value = exploitation + exploration;
    }
    std::cout<< "    - "<< "node_id: "<< mcts_node->node_id<<std::endl;
      // std::cout<< "Exploitation Score: " << exploitation << "; Exploration Score: "<< exploration << "; UTC value: "<<exploitation + exploration <<std::endl;
    std::cout<< "      "<< "reward value: "<< mcts_node->reward <<std::endl;
    std::cout<< "      "<< "exploitation: "<< exploitation <<std::endl;
    std::cout<< "      "<< "exploration: "<< exploration <<std::endl;
    std::cout<< "      "<< "utc value: "<< exploitation + exploration <<std::endl;
  }
}



MCTNode* Planner::MCT_random_successor_generator(Config& C_new, MCTNode* mcts_node, 
  std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, HNode* H_goal){
  MCTNode* MCT_new = nullptr;
  for (int k = 0; k < 10; k++){
    // try 10 times to generate a new successor;
    auto order  = mcts_node->hNode->order;
    if(mcts_node->first_branch){
        order = global_order;
    }else{
      std::shuffle(order.begin(), order.end(), *MT);
    }
    const auto res = get_next_configuration_rollout(mcts_node->hNode->C,  order, nullptr);
    if (!res) {
      continue;
    } 
    for (auto a : A) C_new[a->id] = a->v_next;
    // check explored list 
    auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      rewrite_no_push(mcts_node->hNode, iter->second, H_goal);
      MCT_new = new MCTNode(mcts_node,iter->second, global_MCT_node_id);
      iter->second->set_priority_and_order(order);
      iter->second->reordering(N,D);
      global_MCT_node_id ++;
      MCT_new->setCurrentGValue(mcts_node->curr_g_value + get_edge_cost(mcts_node->hNode->C, C_new));
      mcts_node->first_branch = false;
      return MCT_new;
    }else{
      const auto H_new = new HNode(
        C_new, D, mcts_node->hNode, 
        mcts_node->hNode->g + get_edge_cost(mcts_node->hNode->C, C_new), get_h_value(C_new));
      H_new->setNodeID(global_node_id);
      global_node_id ++;
      H_new->setMakeSpan(mcts_node->hNode->current_make_span + 1);
      H_new->set_priority_and_order(order);
      H_new->reordering(N,D);
      EXPLORED[H_new->C] = H_new;
      MCT_new = new MCTNode(mcts_node,H_new, global_MCT_node_id);
      global_MCT_node_id ++;
      // set MCT node g value;
      MCT_new->setCurrentGValue(mcts_node->curr_g_value + get_edge_cost(mcts_node->hNode->C, C_new));
      mcts_node->first_branch = false;
      return MCT_new;
    }
  }
  MCT_new->completed_node = true;
  return MCT_new;
}

   

MCTNode* Planner::MCT_Learn_order_to_branch(Config& C_new, MCTNode* mcts_node, 
  std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, HNode* H_goal){
  MCTNode* MCT_new = nullptr;
  for (int k = 0; k < 10; k++){
    // try 10 times to generate a new successor;
    auto order  = mcts_node->hNode->order;
    if(mcts_node->first_branch){
        order = global_order;
    }else{
      std::uniform_real_distribution<double> dist(0.0, 0.2); 
      double percentage_preserve = dist(*MT);
      auto preserve_start = static_cast<size_t>(order.size() * (1.0 - percentage_preserve));
      std::sort(order.begin(), order.begin() + preserve_start , [&](uint i, uint j) {
          return mcts_node->agent_ratio[i] < mcts_node->agent_ratio[j]; // Sort in descending order of agent_ratio
      });
    }
    const auto res = get_next_configuration_rollout(mcts_node->hNode->C,  order, nullptr);
    if (!res) {
      continue;
    } 
    for (auto a : A) C_new[a->id] = a->v_next;

    if (mcts_node->generated_configs.find(C_new) != mcts_node->generated_configs.end()) {
      continue; // Skip if the configuration is already generated
    }

   

    // check explored list 
    auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      rewrite_no_push(mcts_node->hNode, iter->second, H_goal);
      MCT_new = new MCTNode(mcts_node,iter->second, global_MCT_node_id);
      iter->second->set_priority_and_order(order);
      iter->second->reordering(N,D);
      global_MCT_node_id ++;
      MCT_new->setCurrentGValue(mcts_node->curr_g_value + get_edge_cost(mcts_node->hNode->C, C_new));
      mcts_node->first_branch = false;
       // Add the new configuration to the set
      mcts_node->generated_configs[iter->second->C] = 0;
      return MCT_new;
    }else{
      const auto H_new = new HNode(
        C_new, D, mcts_node->hNode, 
        mcts_node->hNode->g + get_edge_cost(mcts_node->hNode->C, C_new), get_h_value(C_new));
      H_new->setNodeID(global_node_id);
      global_node_id ++;
      H_new->setMakeSpan(mcts_node->hNode->current_make_span + 1);
      H_new->set_priority_and_order(order);
      H_new->reordering(N,D);
      EXPLORED[H_new->C] = H_new;
      MCT_new = new MCTNode(mcts_node,H_new, global_MCT_node_id);
      global_MCT_node_id ++;
      // set MCT node g value;
      MCT_new->setCurrentGValue(mcts_node->curr_g_value + get_edge_cost(mcts_node->hNode->C, C_new));
      mcts_node->first_branch = false;
      mcts_node->generated_configs[iter->second->C] = 0;
      return MCT_new;
    }
  }
  mcts_node->completed_node = true;
  return MCT_new;
}



MCTNode* Planner::MCT_Learn_order_to_branch(Config& C_new, MCTNode* mcts_node){
  MCTNode* MCT_new = nullptr;
  for (int k = 0; k < 10; k++){
    // try 10 times to generate a new successor;
    auto order  = mcts_node->hNode->order;
    if(mcts_node->first_branch){
        order = global_order;
    }else{
      // std::uniform_real_distribution<double> dist(0.0, 0.2); 
      // double percentage_preserve = dist(*MT);
      // auto preserve_start = static_cast<size_t>(order.size() * (1.0 - percentage_preserve));
      // std::sort(order.begin(), order.begin() + preserve_start , [&](uint i, uint j) {
      //     return mcts_node->agent_ratio[i] < mcts_node->agent_ratio[j]; // Sort in descending order of agent_ratio
      // });

              std::uniform_real_distribution<double> noise_dist(-0.1, 0.1); // Adjust range as needed
        for (size_t i = 0; i < mcts_node->agent_ratio.size(); ++i) {
          mcts_node->agent_ratio[i] = mcts_node->agent_ratio[i] + noise_dist(*MT);
        }
        std::sort(order.begin(), order.end(), [&](uint i, uint j) {
            return (mcts_node->agent_ratio[i] ) < (mcts_node->agent_ratio[j] );
        });
    }
    const auto res = get_next_configuration_rollout(mcts_node->hNode->C,  order, nullptr);
    if (!res) {
      continue;
    } 
    for (auto a : A) C_new[a->id] = a->v_next;

    if (mcts_node->generated_configs.find(C_new) != mcts_node->generated_configs.end()) {
      continue; // Skip if the configuration is already generated
    }
    const auto H_new = new HNode(
      C_new, D, mcts_node->hNode, 
      0, get_h_value(C_new));
    H_new->setNodeID(global_node_id);
    global_node_id ++;
    H_new->setMakeSpan(mcts_node->hNode->current_make_span + 1);
    H_new->set_priority_and_order(order);
    H_new->reordering(N,D);
    MCT_new = new MCTNode(mcts_node,H_new, global_MCT_node_id);
    global_MCT_node_id ++;
    // set MCT node g value;
    MCT_new->setCurrentGValue(mcts_node->curr_g_value + get_edge_cost(mcts_node->hNode->C, C_new));
    mcts_node->first_branch = false;
    mcts_node->generated_configs[H_new->C] = 0;
    return MCT_new;
  }
  mcts_node->completed_node = true;
  return MCT_new;
}



uint Planner::run_completed_lacam(HNode* H_init_node,std::vector<double>& agent_ratio)
{
  // for (auto i = 0; i < N; ++i) A[i] = new Agent(i);
  for (auto i = 0; i < N; ++i){
    A[i]->v_next = nullptr;
    A[i]->v_now = nullptr;
  }
  for (size_t i = 0; i < occupied_next.size(); ++i) {
    occupied_next[i] = nullptr;
    occupied_now[i] = nullptr;
  }

  // setup search
  auto OPEN = std::stack<HNode*>();
  auto EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();
  HNode* H_init = new HNode(
    H_init_node->C,          // Copy the configuration
    D,          // Copy the distance table (if applicable)
    H_init_node->parent,     // Copy the parent pointer
    H_init_node->g,          // Copy the g-value
    H_init_node->h           // Copy the h-value
  );
  H_init->set_priority_and_order(H_init_node->order); // Copy the order
  H_init->reordering(H_init_node->priorities.size(), D); // Reorder priorities


  H_init->parent = nullptr;
  OPEN.push(H_init);
  EXPLORED[H_init->C] = H_init;
  
  auto C_new = Config(N, nullptr);  // for new configuration
  HNode* H_goal = nullptr;          // to store goal node


  while (!OPEN.empty()) {

    // do not pop here!
    auto H = OPEN.top();  // high-level node

    if (H->search_tree.empty()) {
      OPEN.pop();
      continue;
    }

    // check lower bounds
    if (H_goal != nullptr && H->f >= H_goal->f) {
      OPEN.pop();
      continue;
    }

    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      H_goal = H;
      // solver_info(1, "found solution, cost: ", H->g);
      std::fill(agent_ratio.begin(), agent_ratio.end(), 0);
      HNode* current = H_goal;
      while (current != H_init) {
        get_edge_cost_per_agent(agent_ratio,current->C,current->parent->C);
        current = current->parent;
      }
      compute_agent_increase_ratio(agent_ratio, H_init->C, H_goal->C);
    
      uint solution = H_goal->g;
      // solver_info(1, "found solution, cost: ", solution);
      for (auto itr : EXPLORED) delete itr.second;
        // for (auto a : A) delete a;
      return solution;
    }

    // create successors at the low-level search
    auto L = H->search_tree.front();
    H->search_tree.pop();
    expand_lowlevel_tree(H, L);

    // create successors at the high-level search
    const auto res = get_new_config(H, L);
    delete L;  // free
    if (!res) continue;

    // create new configuration
    for (auto a : A) C_new[a->id] = a->v_next;

    // check explored list
    const auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      // case found
      rewrite(H, iter->second, H_goal, OPEN);
      OPEN.push(iter->second);
    } else {
      // insert new search node
      const auto H_new = new HNode(
          C_new, D, H, H->g + get_edge_cost(H->C, C_new), get_h_value(C_new));
      EXPLORED[H_new->C] = H_new;
      OPEN.push(H_new);
    }
  }
}

HNode* Planner::copy_node(HNode* original_node) {
  // Create a new HNode with the same configuration and attributes as the original node
  HNode* copy = new HNode(
      original_node->C,          // Copy the configuration
      D,                         // Copy the distance table
      original_node->parent,     // Copy the parent pointer
      original_node->g,          // Copy the g-value
      original_node->h           // Copy the h-value
  );

  // Copy additional attributes
  copy->set_priority_and_order(original_node->order); // Copy the order
  copy->reordering(original_node->priorities.size(), D); // Reorder priorities

  return copy;
}

uint Planner::run_constrainted_completed_lacam(HNode* H_init_node,std::vector<double>& agent_ratio,
  std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, 
  HNode*& H_goal, std::vector<Config>& solution, std::unordered_set<int>& fixed_agents){
  for (auto i = 0; i < N; ++i){
    A[i]->v_next = nullptr;
    A[i]->v_now = nullptr;
  }
  for (size_t i = 0; i < occupied_next.size(); ++i) {
    occupied_next[i] = nullptr;
    occupied_now[i] = nullptr;
  }

  // setup search
  auto OPEN = std::stack<HNode*>();
  HNode* H_init = new HNode(
    H_init_node->C,          // Copy the configuration
    D,          // Copy the distance table (if applicable)
    H_init_node->parent,     // Copy the parent pointer
    0,    // Copy the g-value
    H_init_node->h           // Copy the h-value
  );
  H_init->set_priority_and_order(H_init_node->order); // Copy the order
  H_init->reordering(H_init_node->priorities.size(), D); // Reorder priorities
  for (auto itr : EXPLORED) delete itr.second;
  EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();

  H_init->parent = nullptr;
  OPEN.push(H_init);
  EXPLORED[H_init->C] = H_init;
  H_init->setMakeSpan(0);
  
  auto C_new = Config(N, nullptr);  // for new configuration
  H_goal = nullptr;          // to store goal node


  while (!OPEN.empty()) {

    // do not pop here!
    auto H = OPEN.top();  // high-level node

    if (H->search_tree.empty()) {
      OPEN.pop();
      continue;
    }

    // check lower bounds
    if (H_goal != nullptr && H->f >= H_goal->f) {
      OPEN.pop();
      continue;
    }

    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      H_goal = H;
      // solver_info(1, "found solution, cost: ", H->g);
      std::fill(agent_ratio.begin(), agent_ratio.end(), 0);
      HNode* current = H_goal;
      while (current != H_init) {
        get_edge_cost_per_agent(agent_ratio,current->C,current->parent->C);
        current = current->parent;
      }
      compute_agent_increase_ratio(agent_ratio, H_init->C, H_goal->C);
    
      uint solution = H_goal->g;
      return solution;
    }

    // create successors at the low-level search
    auto L = H->search_tree.front();
    H->search_tree.pop();
    // expand_lowlevel_tree(H, L);

    std::cout << "           "<< std::endl;
    std::cout << "           "<< std::endl;
    std::cout << "Configuration: [";
    for (size_t i = 0; i < H->C.size(); ++i) {
        if (H->C[i] != nullptr) {
            std::cout << "{ID: " << i << ", Index: " << H->C[i]->index << "}";
        } else {
            std::cout << "Null";
        }
        if (i < H->C.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    std::cout << "           "<< std::endl;
    std::cout << "           "<< std::endl;



    std::cout << "constrainted agents:"<< std::endl;
    std::cout << "           "<< std::endl;
    std::cout << "Configuration: [";
    for (auto i : fixed_agents) {
        if (H->C[i] != nullptr) {
            std::cout << "{ID: " << i << ", Index: " << H->C[i]->index << "}";
        } else {
            std::cout << "Null";
        }
        if (i < H->C.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
    std::cout << "           "<< std::endl;
    std::cout << "           "<< std::endl;

    if( solution.size() >  0 && H->current_make_span < solution.size() -1 ){
      expand_lowlevel_tree_avoid_agent(H, L, fixed_agents);
      auto& next_config = solution[H->current_make_span + 1];
      for (auto a : fixed_agents) {
        if(next_config[a] != ins->goals[a]){
          L->who.push_back(a);
          L->where.push_back(next_config[a]);
        }
      }
      // expand_lowlevel_tree(H, L);
    }else{
      expand_lowlevel_tree(H, L);
    }
    
    // create successors at the high-level search
    const auto res = get_new_config(H, L);
    if (!res) {
      bool a =0;
    }

    delete L;  // free
    if (!res) {
      bool a =0;
      continue;
    }

    // create new configuration
    for (auto a : A) C_new[a->id] = a->v_next;

    if( solution.size() >  0 && H->current_make_span < solution.size() -1 ){
      auto& next_config = solution[H->current_make_span + 1];
      for (auto a : fixed_agents) {
        auto b =  C_new[a]->id;
        auto c = next_config[a]->id;
        if(C_new[a]->id != next_config[a]->id){
          std::cout<<" werid" <<std::endl; 
        }
      }
    }
    // check explored list
    const auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      // case found
      rewrite(H, iter->second, H_goal, OPEN);
      H->current_make_span + 1 < iter->second->current_make_span ? 
      iter->second->setMakeSpan(H->current_make_span + 1 ) : 
        iter->second->setMakeSpan(iter->second->current_make_span);
        
        auto H_insert = (MT != nullptr && get_random_float(MT) >= RESTART_RATE)
        ? iter->second
        : H_init;
        if (H_goal == nullptr || H_insert->f < H_goal->f){
        OPEN.push(H_insert);
        // record_highlevel_node(H_insert,false,false);
        } 
    } else {
      // insert new search node
      const auto H_new = new HNode(
          C_new, D, H, H->g + get_edge_cost(H->C, C_new), get_h_value(C_new));
      H_new->setMakeSpan(H->current_make_span + 1);
      EXPLORED[H_new->C] = H_new;
      OPEN.push(H_new);
    }
  } 


}




uint Planner::run_completed_lacam(HNode* H_init_node,std::vector<double>& agent_ratio,
  std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, HNode*& H_goal
)
{
  // for (auto i = 0; i < N; ++i) A[i] = new Agent(i);
  for (auto i = 0; i < N; ++i){
    A[i]->v_next = nullptr;
    A[i]->v_now = nullptr;
  }
  for (size_t i = 0; i < occupied_next.size(); ++i) {
    occupied_next[i] = nullptr;
    occupied_now[i] = nullptr;
  }

  // setup search
  auto OPEN = std::stack<HNode*>();
  HNode* H_init = new HNode(
    H_init_node->C,          // Copy the configuration
    D,          // Copy the distance table (if applicable)
    H_init_node->parent,     // Copy the parent pointer
    0,    // Copy the g-value
    H_init_node->h           // Copy the h-value
  );
  H_init->set_priority_and_order(H_init_node->order); // Copy the order
  H_init->reordering(H_init_node->priorities.size(), D); // Reorder priorities
  for (auto itr : EXPLORED) delete itr.second;
  EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();

  H_init->parent = nullptr;
  OPEN.push(H_init);
  EXPLORED[H_init->C] = H_init;
  
  auto C_new = Config(N, nullptr);  // for new configuration
  H_goal = nullptr;          // to store goal node


  while (!OPEN.empty()) {

    // do not pop here!
    auto H = OPEN.top();  // high-level node

    if (H->search_tree.empty()) {
      OPEN.pop();
      continue;
    }

    // check lower bounds
    if (H_goal != nullptr && H->f >= H_goal->f) {
      OPEN.pop();
      continue;
    }

    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      H_goal = H;
      // solver_info(1, "found solution, cost: ", H->g);
      std::fill(agent_ratio.begin(), agent_ratio.end(), 0);
      HNode* current = H_goal;
      while (current != H_init) {
        get_edge_cost_per_agent(agent_ratio,current->C,current->parent->C);
        current = current->parent;
      }
      compute_agent_increase_ratio(agent_ratio, H_init->C, H_goal->C);
    
      uint solution = H_goal->g;
      return solution;
    }

    // create successors at the low-level search
    auto L = H->search_tree.front();
    H->search_tree.pop();
    expand_lowlevel_tree(H, L);

    // create successors at the high-level search
    const auto res = get_new_config(H, L);
    delete L;  // free
    if (!res) continue;

    // create new configuration
    for (auto a : A) C_new[a->id] = a->v_next;

    // check explored list
    const auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      // case found
      rewrite(H, iter->second, H_goal, OPEN);
      auto H_insert = (MT != nullptr && get_random_float(MT) >= RESTART_RATE)
                          ? iter->second
                          : H_init;
      if (H_goal == nullptr || H_insert->f < H_goal->f){
        OPEN.push(H_insert);
        // record_highlevel_node(H_insert,false,false);
      } 
      // OPEN.push(iter->second);
    } else {
      // insert new search node
      const auto H_new = new HNode(
          C_new, D, H, H->g + get_edge_cost(H->C, C_new), get_h_value(C_new));
      EXPLORED[H_new->C] = H_new;
      OPEN.push(H_new);
    }
  }
}

MCTNode* Planner::MCT_select_successors_from_pool(MCTNode* mcts_node,  
  std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, 
  std::vector<Config>& solution, HNode*& H_goal)
{

  while(  mcts_node->generated_configs.size() == 0){
    auto results = multiple_simulations(10, H_goal, mcts_node, EXPLORED, solution);
    if(results == std::numeric_limits<uint>::max()){
      mcts_node->completed_node = true;
      return nullptr;
    }
  }
  
  HNode* h_next = nullptr;
  uint best_simulate_results = std::numeric_limits<uint>::max();
  auto best_iter = mcts_node->generated_configs.end();
  
  for (auto iter = mcts_node->generated_configs.begin(); iter != mcts_node->generated_configs.end(); ++iter) {
      if (best_simulate_results > iter->second->simulation_cost) {
          best_simulate_results = iter->second->simulation_cost;
          h_next = iter->second;
          best_iter = iter; // Keep track of the best iterator
      }
  }
  //   std::cout<<" "<<std::endl;
  // std::cout<<" "<<std::endl;
  // std::cout<<" select cost :"<< best_simulate_results<<std::endl;
  // std::cout<<" "<<std::endl;

  // auto agent_ratio = std::vector<double>(N, 0);

  h_next->setNodeID(global_node_id);
  global_node_id ++;
  h_next->setMakeSpan(mcts_node->hNode->current_make_span + 1);
  MCTNode* MCT_new = new MCTNode(mcts_node,h_next, global_MCT_node_id);
  global_MCT_node_id ++;
  MCT_new->setCurrentGValue(mcts_node->curr_g_value + get_edge_cost(mcts_node->hNode->C, h_next->C));
  // Remove the selected iter from generated_configs
  if (best_iter != mcts_node->generated_configs.end()) {
      mcts_node->generated_configs.erase(best_iter);
  }

  return MCT_new;
}


uint Planner::multiple_simulations(int num_simulations, HNode*& H_goal,
   MCTNode* mcts_node, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, std::vector<Config>& solution){
    uint best_simulate_results = std::numeric_limits<uint>::max();

    // num_simulations = 5;
    for (int k = 0; k < num_simulations; k++){
      // try 10 times to generate a new successor;
      auto order = mcts_node->hNode->order;
      if(k > 0){
        std::uniform_real_distribution<double> noise_dist(-0.1, 0.1); // Adjust range as needed
        for (size_t i = 0; i < mcts_node->agent_ratio.size(); ++i) {
          mcts_node->agent_ratio[i] = mcts_node->agent_ratio[i] * (1 + noise_dist(*MT));
        }
        std::sort(order.begin(), order.end(), [&](uint i, uint j) {
            return (mcts_node->agent_ratio[i] ) < (mcts_node->agent_ratio[j] );
        });
      }
      mcts_node->hNode->order = order;
      auto simulate_result = run_completed_lacam(mcts_node->hNode,mcts_node->agent_ratio,EXPLORED,H_goal);

      if(simulate_result == -1){
        continue;
      }
      simulate_result += mcts_node->curr_g_value;
      if(min_f_cost > simulate_result ){
        min_f_cost = simulate_result ;
        best_mct_node = mcts_node;
        solution = std::vector<Config>();
        auto H = H_goal;
        while(H != nullptr ){
          solution.push_back(H->C);
          H = H->parent;
        }
        solver_info(1, "found solution, cost: ", min_f_cost);
      }

      if(best_simulate_results > simulate_result){
        best_simulate_results = simulate_result;
      }
      // if(mcts_node->best_simulation_cost > simulate_result){
      //   mcts_node->best_simulation_cost = simulate_result;
        HNode* Successors = nullptr; 
        auto H = H_goal;
        while(H != nullptr ){
          if(H->parent->parent == nullptr){
            if (mcts_node->generated_configs.find(H->C) != mcts_node->generated_configs.end()) {
              if(mcts_node->generated_configs[H->C]->simulation_cost > simulate_result){
                mcts_node->generated_configs[H->C]->simulation_cost = simulate_result;
                mcts_node->generated_configs[H->C]->order = H->order ; // Copy the order
              }
              std::cout<< "found existing node"<< std::endl;
            }else{
              Successors = new HNode(
                H->C, D, nullptr, 0, H->h);
              Successors->order = H->order;
              Successors->simulation_cost = simulate_result;
              mcts_node->generated_configs[Successors->C] = Successors;
            }
            break; 
          }
          H = H->parent;
        }
      // }
      // if(mcts_node->generated_configs.size() > 200){
      //   // std::cout<< "number of successors:" << mcts_node->generated_configs.size() <<std::endl;
      // }
      // std::cout<< "number of successors:" << mcts_node->generated_configs.size() <<std::endl;
    }
  // std::cout<< "Returned simulation cost "<< best_simulate_results <<std::endl;
    return best_simulate_results;
}


Solution Planner::MCT_multiple_lacam(std::string& additional_info){
  // setup agents
  for (auto i = 0; i < N; ++i) A[i] = new Agent(i);


  auto OPEN = std::vector<MCTNode*>();
  auto H_init = new HNode(ins->starts, D, nullptr, 0, get_h_value(ins->starts));
  H_init->setNodeID(global_node_id);
  global_node_id ++;
  global_order = H_init->order;
  global_MCT_node_id = 1;

  
  auto MCT_init = new MCTNode(nullptr,H_init,global_MCT_node_id);
  global_MCT_node_id ++; 
  H_init->setMakeSpan(0);
  OPEN.push_back(MCT_init);
  MCT_init->initialize_agent_ratio(N);

  minMaxStats = MinMaxStats();
  min_f_cost = std::numeric_limits<uint>::max();
  best_mct_node = nullptr;
  std::vector<Config> solution;
  auto C_new = Config(N, nullptr);  // for new configuration
  auto EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();
  HNode* H_goal = nullptr;


  // for( int j = 0; j < 100000; j ++){
  //   if( j != 0){
  //     std::uniform_real_distribution<double> noise_dist(-0.1, 0.1); // Adjust range as needed
  //     for (size_t i = 0; i < MCT_init->agent_ratio.size(); ++i) {
  //       MCT_init->agent_ratio[i] = MCT_init->agent_ratio[i] * (1 + noise_dist(*MT));
  //     }
  //     std::sort(MCT_init->hNode->order.begin(),MCT_init->hNode->order.end(), [&](uint i, uint j) {
  //         return (MCT_init->agent_ratio[i] ) < (MCT_init->agent_ratio[j] );
  //     });
  //   }
  //   // std::shuffle(MCT_init->hNode->order.begin(), MCT_init->hNode->order.end(), *MT);
  //   auto simulate_result = run_completed_lacam(MCT_init->hNode,MCT_init->agent_ratio
  //     , EXPLORED, H_goal);
  //   if(min_f_cost > simulate_result ){
  //     min_f_cost = simulate_result ;
  //     solver_info(1, "found solution, cost: ", min_f_cost);
  //   }
  //   // solver_info(1, "found solution, cost: ", min_f_cost);
  // }


  for( int j = 0; j < 100000; j ++){
    std::unordered_set<int> fixed_agents = std::unordered_set<int>();
    if( j != 0){
      std::uniform_real_distribution<double> noise_dist(-0.1, 0.1); // Adjust range as needed
      for (size_t i = 0; i < MCT_init->agent_ratio.size(); ++i) {
        MCT_init->agent_ratio[i] = MCT_init->agent_ratio[i] * (1 + noise_dist(*MT));
      }
      std::sort(MCT_init->hNode->order.begin(),MCT_init->hNode->order.end(), [&](uint i, uint j) {
          return (MCT_init->agent_ratio[i] ) < (MCT_init->agent_ratio[j] );
      });

      for( int k = 0; k < N * 0.1; k ++){
        fixed_agents.insert(MCT_init->hNode->order[k]);
      }
    }
    // std::shuffle(MCT_init->hNode->order.begin(), MCT_init->hNode->order.end(), *MT);
    auto simulate_result = run_constrainted_completed_lacam(
    MCT_init->hNode,MCT_init->agent_ratio, EXPLORED, H_goal, solution, fixed_agents);
    // auto simulate_result = run_completed_lacam(MCT_init->hNode,MCT_init->agent_ratio
          // , EXPLORED, H_goal);

    solution = std::vector<Config>();
    auto H = H_goal;
    while(H != nullptr ){
      solution.push_back(H->C);
      H = H->parent;
    }
    std::reverse(solution.begin(), solution.end());

    if(min_f_cost > simulate_result ){
      min_f_cost = simulate_result ;
      solver_info(1, "found solution, cost: ", min_f_cost);
    }
    // solver_info(1, "found solution, cost: ", min_f_cost);
  }
  auto reward = multiple_simulations(10, H_goal, MCT_init, EXPLORED, solution);
  if(reward  == -1){
    MCT_backpropagate(MCT_init,1,0,1);
  }else{
    MCT_backpropagate(MCT_init,1,reward,0);
    minMaxStats.update(reward);
  }

  if(verbose == -1){
    std::cout<< "version: 1.4.0\n";
    std::cout<< "events:\n";
    std::cout<< " - type: expanding" <<std::endl;
    std::cout<< "   id: "<< MCT_init->node_id <<std::endl;
    std::cout<< "   pId: " 
          << (MCT_init->parent == nullptr ? "null" : std::to_string(MCT_init->parent->node_id)) 
          << std::endl;
    std::cout<< "   visits: "<< MCT_init->visits<<std::endl;
    std::cout<< "   reward: "<< MCT_init->reward <<std::endl;
    print_utc_value(OPEN);
  }
  for(int i = 0; i < 10; i++){
    // const auto MCT_new  = MCT_random_successor_generator(C_new, MCT_init, EXPLORED,H_goal);
    const auto MCT_new  = MCT_select_successors_from_pool(MCT_init, EXPLORED, solution,H_goal);
    if(MCT_new != nullptr){
      MCT_new ->initialize_agent_ratio(N);
      auto simulation_cost = multiple_simulations(10, H_goal, MCT_new, EXPLORED, solution);
      // std::cout<<"Simulation cost: "<< simulation_cost <<std::endl;
      if(simulation_cost  == -1){
        MCT_backpropagate(MCT_new,1,0,1);
      }else{
        // simulation_cost += MCT_new->curr_g_value;
        // std::cout<< "simulation cost: "<< simulation_cost <<std::endl;
        MCT_backpropagate(MCT_new,1,simulation_cost,0);
        minMaxStats.update(simulation_cost);
      }
      OPEN.push_back(MCT_new);
      if(verbose == -1){
        std::cout<< " - type: generating" <<std::endl;
        std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
        std::cout<< "   pId: " 
          << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
          << std::endl;
        std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
        std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
      }
    } 
  }
  while (!is_expired(deadline)) {
    loop_cnt += 1;

    // do not pop here!
    auto MCT_NODE =  MCT_selection(OPEN);// high-level node
    if(MCT_NODE == nullptr){
      break;
    }
    // std::cout<<"Selected !" <<std::endl;
    // std::cout << "Selecting node with ID: " << MCT_NODE->node_id << " ;" <<std::endl;
    auto H  = MCT_NODE->hNode; 
    if (H->search_tree.empty()) {
      MCT_NODE -> completed_node = true;
      continue;
    }
    // check lower bounds
    if (MCT_NODE->curr_g_value >= min_f_cost) {
      MCT_NODE -> completed_node = true;
      continue;
    }

    if(verbose == -1){
      std::cout<< " - type: expanding" <<std::endl;
      std::cout<< "   id: "<< MCT_NODE->node_id <<std::endl;
      std::cout<< "   pId: " 
            << (MCT_NODE->parent == nullptr ? "null" : std::to_string(MCT_NODE->parent->node_id)) 
            << std::endl;
      std::cout<< "   visits: "<< MCT_NODE->visits<<std::endl;
      std::cout<< "   reward: "<< MCT_NODE->reward <<std::endl;
      print_utc_value(OPEN);
    }

    // check goal condition
    if ( is_same_config(H->C, ins->goals)) {
      min_f_cost = MCT_NODE->curr_g_value;
      best_mct_node = MCT_NODE;
      MCT_NODE -> completed_node = true;
      solver_info(1, "found solution, cost: ", min_f_cost);
      continue;
    }

    // std::cout<<"I am here !" <<std::endl;
    // create successors at the low-level search
    MCTNode* MCT_new  = MCT_select_successors_from_pool(MCT_NODE, EXPLORED, solution,H_goal);
    if(MCT_new == nullptr){
      MCT_backpropagate(MCT_NODE,1,0,1);
      continue;
    }else{
      MCT_new ->initialize_agent_ratio(N);
      // auto simulation_cost = Lacam_simulator(MCT_new->hNode, H_goal, EXPLORED, 
        // false, node_budget, MCT_NODE->agent_ratio);
      auto simulation_cost = multiple_simulations(10, H_goal, MCT_new, EXPLORED, solution);
      // std::cout<<"Simulation cost: "<< simulation_cost <<std::endl;
      if(simulation_cost  == -1){
        MCT_backpropagate(MCT_new,1,0,1);
      }else{
        // simulation_cost += MCT_new->curr_g_value;
        MCT_backpropagate(MCT_new,1,simulation_cost,0);
        minMaxStats.update(simulation_cost);
      }
      OPEN.push_back(MCT_new); 
      if(verbose == -1){
        std::cout<< " - type: generating" <<std::endl;
        std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
        std::cout<< "   pId: " 
          << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
          << std::endl;
        std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
        std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
      }
    } 
    // }
  }

  if(best_mct_node->parent != nullptr){
    auto H = best_mct_node->parent;
    while(H != nullptr ){
      solution.push_back(H->hNode->C);
      H = H->parent;
    }
  }
  if(solution.size() != 0){
    std::reverse(solution.begin(), solution.end());
  }
  // print result
  if (best_mct_node  != nullptr && OPEN.empty()) {
    solver_info(1, "solved optimally, objective: ", objective);
  } else if (best_mct_node  != nullptr) {
    solver_info(1, "solved sub-optimally, objective: ", objective);
  } else if (OPEN.empty()) {
    solver_info(1, "no solution");
  } else {
    solver_info(1, "timeout");
  }

  // logging
  additional_info +=
      "optimal=" + std::to_string(best_mct_node != nullptr && OPEN.empty()) + "\n";
  additional_info += "objective=" + std::to_string(objective) + "\n";
  additional_info += "loop_cnt=" + std::to_string(loop_cnt) + "\n";
  additional_info += "num_node_gen=" + std::to_string(OPEN.size()) + "\n";


  // save to tree file 
  if(tree_file != "none"){
    saveTree(tree_file);
  }

  // memory management
  for (auto a : A) delete a;
  for (auto itr : OPEN){
    delete itr;
  } 

  return solution;






}


Solution Planner::MCT_vanilla_lacam(std::string& additional_info)
{
  // setup agents
  for (auto i = 0; i < N; ++i) A[i] = new Agent(i);

  // setup search
  auto OPEN = std::vector<MCTNode*>();
  auto H_init = new HNode(ins->starts, D, nullptr, 0, get_h_value(ins->starts));
  H_init->setNodeID(global_node_id);
  global_node_id ++;
  global_order = H_init->order;
  global_MCT_node_id = 1;

  
  auto MCT_init = new MCTNode(nullptr,H_init,global_MCT_node_id);
  global_MCT_node_id ++; 
  H_init->setMakeSpan(0);
  OPEN.push_back(MCT_init);
  MCT_init->initialize_agent_ratio(N);

  minMaxStats = MinMaxStats();
  min_f_cost = std::numeric_limits<uint>::max();
  best_mct_node = nullptr;
  std::vector<Config> solution;
  auto C_new = Config(N, nullptr);  // for new configuration
  auto EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();
  HNode* H_goal = nullptr;
  // const auto H_copy = new HNode(
  //   H_init->C, D, nullptr, H_init->g, H_init->h);
  // H_copy->set_priority_and_order(H_init->order);// Create a non-const copy
  // H_copy->reordering(N, D);
  
  auto r = run_completed_lacam(H_init,MCT_init->agent_ratio,EXPLORED,H_goal);
  // Lacam_simulator(MCT_init->hNode, H_goal, EXPLORED, true, 0,MCT_init->agent_ratio);
  if(r  == -1){
    MCT_backpropagate(MCT_init,1,0,1);
  }else{
    MCT_backpropagate(MCT_init,1,r,0);
    minMaxStats.update(r);
    if(min_f_cost > r){
      min_f_cost = r;
      best_mct_node = MCT_init;
      solution = std::vector<Config>();
      auto H = H_goal;
      while(H != nullptr ){
        solution.push_back(H->C);
        H = H->parent;
      }
      solver_info(1, "found solution, cost: ", min_f_cost);
    }
  }

  // if(verbose == -1){
  //   std::cout<< "version: 1.4.0\n";
  //   std::cout<< "events:\n";
  //   std::cout<< " - type: expanding" <<std::endl;
  //   std::cout<< "   id: "<< MCT_init->node_id <<std::endl;
  //   std::cout<< "   pId: " 
  //         << (MCT_init->parent == nullptr ? "null" : std::to_string(MCT_init->parent->node_id)) 
  //         << std::endl;
  //   std::cout<< "   visits: "<< MCT_init->visits<<std::endl;
  //   std::cout<< "   reward: "<< MCT_init->reward <<std::endl;
  //   print_utc_value(OPEN);
  // }
  
  // uint node_budget = H_goal->current_make_span;
  // for the sake of MCTS, let's force the root not have k branch factors.
  for(int i = 0; i < 10; i++){
    // const auto MCT_new  = MCT_random_successor_generator(C_new, MCT_init, EXPLORED,H_goal);
    const auto MCT_new  = MCT_Learn_order_to_branch(C_new, MCT_init);
    if(MCT_new != nullptr){
      // auto simulation_cost = Lacam_simulator(MCT_new->hNode, H_goal, EXPLORED, 
      //   false, node_budget,MCT_init->agent_ratio);
      auto simulation_cost = run_completed_lacam(MCT_new->hNode,MCT_init->agent_ratio,EXPLORED,H_goal);
      // std::cout<< "Results: "<< results <<std::endl;
      if(simulation_cost  == -1){
        MCT_backpropagate(MCT_new,1,0,1);
      }else{
        simulation_cost += MCT_new->curr_g_value;
        MCT_backpropagate(MCT_new,1,simulation_cost,0);
        minMaxStats.update(simulation_cost);
        if(min_f_cost > simulation_cost){
          min_f_cost = simulation_cost;
          best_mct_node = MCT_new;

          solution = std::vector<Config>();
          auto H = H_goal;
          while(H != nullptr ){
            solution.push_back(H->C);
            H = H->parent;
          }
          solver_info(1, "found solution, cost: ", min_f_cost);
        }
        // solver_info(1, "found solution, cost: ", simulation_cost);
      }
      MCT_new ->initialize_agent_ratio(N);
      OPEN.push_back(MCT_new);

      // if(verbose == -1){
      //   std::cout<< " - type: generating" <<std::endl;
      //   std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
      //   std::cout<< "   pId: " 
      //     << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
      //     << std::endl;
      //   std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
      //   std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
      // }
    } 
  }

  while (!is_expired(deadline)) {
    // if(H_goal != nullptr){
    //   roll_out_cnt = H_goal->current_make_span;
    // }
    // node_budget = H_goal->current_make_span;
    loop_cnt += 1;

    // do not pop here!
    auto MCT_NODE =  MCT_selection(OPEN);// high-level node
    if(MCT_NODE == nullptr){
      break;
    }
    // std::cout << "Selecting node with ID: " << MCT_NODE->node_id << " ;" <<std::endl;
    auto H  = MCT_NODE->hNode; 
    // std::cout << "#Nodes in Queue: " << OPEN.size() <<" #Node ID: "<< MCT_NODE->node_id << " UTC Value: " << MCT_NODE->compute_uct_value() << std::endl ;

    if (H->search_tree.empty()) {
      MCT_NODE -> completed_node = true;
      continue;
    }
    // check lower bounds
    if (MCT_NODE->curr_g_value >= min_f_cost) {
      MCT_NODE -> completed_node = true;
      continue;
    }

    // if(verbose == -1){
    //   std::cout<< " - type: expanding" <<std::endl;
    //   std::cout<< "   id: "<< MCT_NODE->node_id <<std::endl;
    //   std::cout<< "   pId: " 
    //         << (MCT_NODE->parent == nullptr ? "null" : std::to_string(MCT_NODE->parent->node_id)) 
    //         << std::endl;
    //   std::cout<< "   visits: "<< MCT_NODE->visits<<std::endl;
    //   std::cout<< "   reward: "<< MCT_NODE->reward <<std::endl;
    //   print_utc_value(OPEN);
    // }

    // check goal condition
    if ( is_same_config(H->C, ins->goals)) {
      min_f_cost = MCT_NODE->curr_g_value;
      best_mct_node = MCT_NODE;
      MCT_NODE -> completed_node = true;
      solver_info(1, "found solution, cost: ", min_f_cost);
      // std::cout<<" aaaaa"<<std::endl;
      continue;
    }


    // create successors at the low-level search
    MCTNode* MCT_new  = MCT_Learn_order_to_branch(C_new, MCT_NODE);
    if(MCT_new == nullptr){
      MCT_backpropagate(MCT_NODE,1,0,1);
      continue;
    }else{
      // auto simulation_cost = Lacam_simulator(MCT_new->hNode, H_goal, EXPLORED, 
        // false, node_budget, MCT_NODE->agent_ratio);
      auto simulation_cost = run_completed_lacam(MCT_new->hNode,MCT_NODE->agent_ratio,EXPLORED,H_goal);
      if(simulation_cost  == -1){
        MCT_backpropagate(MCT_new,1,0,1);
      }else{
        simulation_cost += MCT_new->curr_g_value;
        MCT_backpropagate(MCT_new,1,simulation_cost,0);
        minMaxStats.update(simulation_cost);
        if(min_f_cost > simulation_cost){
          min_f_cost = simulation_cost;
          best_mct_node = MCT_new;
          solver_info(1, "found solution, cost: ", min_f_cost);
          solution = std::vector<Config>();
          auto H = H_goal;
          while(H != nullptr ){
            solution.push_back(H->C);
            H = H->parent;
          }
        }
      }

      MCT_new ->initialize_agent_ratio(N);
      OPEN.push_back(MCT_new); 
      // if(verbose == -1){
      //   std::cout<< " - type: generating" <<std::endl;
      //   std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
      //   std::cout<< "   pId: " 
      //     << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
      //     << std::endl;
      //   std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
      //   std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
      // }
    } 
    // }
  }
  // for( int i = 1 ; i < 11; i ++){
  //   std::cout<< OPEN[i]->visits << std::endl; 
  // }

  // std::cout<< "SEARCH ENDED" <<std::endl;
  // std::cout<< minMaxStats.minimum<<std::endl;
  // std::cout<< minMaxStats.maximum <<std::endl;
  // std::cout<< "Simulation times: "<< simulation_times <<std::endl;
  // std::cout<< "Walk back times: "<< walk_back_times <<std::endl;
  // backtrack
  if(best_mct_node->parent != nullptr){
    auto H = best_mct_node->parent;
    while(H != nullptr ){
      solution.push_back(H->hNode->C);
      H = H->parent;
    }
  }
  if(solution.size() != 0){
    std::reverse(solution.begin(), solution.end());
  }
  // std::cout<< "rewrite_calls: "<< rewrite_calls <<std::endl; 
  // std::cout<< "sample_times"<< sampleing_times <<std::endl; 
  // print result
  if (best_mct_node  != nullptr && OPEN.empty()) {
    solver_info(1, "solved optimally, objective: ", objective);
  } else if (best_mct_node  != nullptr) {
    solver_info(1, "solved sub-optimally, objective: ", objective);
  } else if (OPEN.empty()) {
    solver_info(1, "no solution");
  } else {
    solver_info(1, "timeout");
  }

  // logging
  additional_info +=
      "optimal=" + std::to_string(best_mct_node != nullptr && OPEN.empty()) + "\n";
  additional_info += "objective=" + std::to_string(objective) + "\n";
  additional_info += "loop_cnt=" + std::to_string(loop_cnt) + "\n";
  additional_info += "num_node_gen=" + std::to_string(OPEN.size()) + "\n";


  // save to tree file 
  if(tree_file != "none"){
    saveTree(tree_file);
  }

  // memory management
  for (auto a : A) delete a;
  for (auto itr : OPEN){
    delete itr;
  } 

  return solution;
}





Solution Planner::MCT_lacam(std::string& additional_info)
{
  // setup agents
  for (auto i = 0; i < N; ++i) A[i] = new Agent(i);

  // setup search
  auto OPEN = std::vector<MCTNode*>();
  auto EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();
  auto H_init = new HNode(ins->starts, D, nullptr, 0, get_h_value(ins->starts));
  H_init->setNodeID(global_node_id);
  global_node_id ++;
  global_order = H_init->order;
  global_MCT_node_id = 1;

  
  auto MCT_init = new MCTNode(nullptr,H_init,global_MCT_node_id);
  global_MCT_node_id ++; 
  H_init->setMakeSpan(0);
  OPEN.push_back(MCT_init);
  MCT_init->initialize_agent_ratio(N);
  EXPLORED[H_init->C] = H_init;
  uint roll_out_cnt =  2 * get_makespan_lower_bound(ins->starts);
  uint makespan = get_makespan_lower_bound(ins->starts);
  minMaxStats = MinMaxStats();
  uint min_f_cost = std::numeric_limits<uint>::max();


  std::vector<Config> solution;
  auto C_new = Config(N, nullptr);  // for new configuration
  HNode* H_goal = nullptr;          // to store goal node
  
  
  // auto r = Lacam_simulator(MCT_init->hNode, H_goal, EXPLORED, true, 0,MCT_init->agent_ratio);
  
  const auto H_copy = new HNode(
    H_init->C, D, nullptr, H_init->g, H_init->h);
  H_copy->set_priority_and_order(H_init->order);// Create a non-const copy
  H_copy->reordering(N, D);
  
  auto r = run_completed_lacam(H_copy,MCT_init->agent_ratio);
  // Lacam_simulator(MCT_init->hNode, H_goal, EXPLORED, true, 0,MCT_init->agent_ratio);
  if(r  == -1){
    MCT_backpropagate(MCT_init,1,0,1);
  }else{
    MCT_backpropagate(MCT_init,1,r,0);
    minMaxStats.update(r);
    if(min_f_cost > r){
      min_f_cost = r;
      // std::cout<< r << std::endl;
      solver_info(1, "found solution, cost: ", min_f_cost);
    }
  }

  if(verbose == -1){
    std::cout<< "version: 1.4.0\n";
    std::cout<< "events:\n";
    std::cout<< " - type: expanding" <<std::endl;
    std::cout<< "   id: "<< MCT_init->node_id <<std::endl;
    std::cout<< "   pId: " 
          << (MCT_init->parent == nullptr ? "null" : std::to_string(MCT_init->parent->node_id)) 
          << std::endl;
    std::cout<< "   visits: "<< MCT_init->visits<<std::endl;
    std::cout<< "   reward: "<< MCT_init->reward <<std::endl;
    print_utc_value(OPEN);
  }
  
  // uint node_budget = H_goal->current_make_span;
  // for the sake of MCTS, let's force the root not have k branch factors.
  for(int i = 0; i < 10; i++){
    // const auto MCT_new  = MCT_random_successor_generator(C_new, MCT_init, EXPLORED,H_goal);
    const auto MCT_new  = MCT_Learn_order_to_branch(C_new, MCT_init, EXPLORED,H_goal);
    if(MCT_new != nullptr){
      // auto simulation_cost = Lacam_simulator(MCT_new->hNode, H_goal, EXPLORED, 
      //   false, node_budget,MCT_init->agent_ratio);
      const auto H_copy = new HNode(
        MCT_new->hNode->C, D, nullptr, 0, MCT_new->hNode->h);
      H_copy->set_priority_and_order(MCT_new->hNode->order);// Create a non-const copy
      H_copy->reordering(N, D);
      auto simulation_cost = run_completed_lacam(H_copy,MCT_init->agent_ratio);
      // std::cout<< "Results: "<< results <<std::endl;
      if(simulation_cost  == -1){
        MCT_backpropagate(MCT_new,1,0,1);
      }else{
        simulation_cost += MCT_new->curr_g_value;
        MCT_backpropagate(MCT_new,1,simulation_cost,0);
        minMaxStats.update(simulation_cost);
        if(min_f_cost > simulation_cost){
          min_f_cost = simulation_cost;
          solver_info(1, "found solution, cost: ", min_f_cost);
        }
        // solver_info(1, "found solution, cost: ", simulation_cost);
      }
      MCT_new ->initialize_agent_ratio(N);
      OPEN.push_back(MCT_new);

      if(verbose == -1){
        std::cout<< " - type: generating" <<std::endl;
        std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
        std::cout<< "   pId: " 
          << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
          << std::endl;
        std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
        std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
      }
    } 
  }



  while (!is_expired(deadline)) {
    // if(H_goal != nullptr){
    //   roll_out_cnt = H_goal->current_make_span;
    // }
    // node_budget = H_goal->current_make_span;
    loop_cnt += 1;

    // do not pop here!
    auto MCT_NODE =  MCT_selection(OPEN,H_goal);// high-level node
    if(MCT_NODE == nullptr){
      break;
    }
    // std::cout << "Selecting node with ID: " << MCT_NODE->node_id << " ;" <<std::endl;
    auto H  = MCT_NODE->hNode; 
    // std::cout << "#Nodes in Queue: " << OPEN.size() <<" #Node ID: "<< MCT_NODE->node_id << " UTC Value: " << MCT_NODE->compute_uct_value() << std::endl ;

    if (H->search_tree.empty()) {
      MCT_NODE -> completed_node = true;
      continue;
    }
    // check lower bounds
    if (H_goal != nullptr && H->f >= H_goal->f) {
      MCT_NODE -> completed_node = true;
      continue;
    }

    if(verbose == -1){
      std::cout<< " - type: expanding" <<std::endl;
      std::cout<< "   id: "<< MCT_NODE->node_id <<std::endl;
      std::cout<< "   pId: " 
            << (MCT_NODE->parent == nullptr ? "null" : std::to_string(MCT_NODE->parent->node_id)) 
            << std::endl;
      std::cout<< "   visits: "<< MCT_NODE->visits<<std::endl;
      std::cout<< "   reward: "<< MCT_NODE->reward <<std::endl;
      print_utc_value(OPEN);
    }
    // std::cout<< "level: "<<H->f <<std::endl;
    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      MCT_NODE -> completed_node = true;
      continue;
    }


    // create successors at the low-level search
    // MCTNode* MCT_new = MCT_random_successor_generator(C_new, MCT_NODE, EXPLORED,H_goal);
    // for(int i = 0; i < 5; i++){
      MCTNode* MCT_new  = MCT_Learn_order_to_branch(C_new, MCT_NODE, EXPLORED,H_goal);
      if(MCT_new == nullptr){
        MCT_backpropagate(MCT_NODE,1,0,1);
        continue;
      }else{
        // auto simulation_cost = Lacam_simulator(MCT_new->hNode, H_goal, EXPLORED, 
          // false, node_budget, MCT_NODE->agent_ratio);
        
        const auto H_copy = new HNode(
          MCT_new->hNode->C, D, nullptr, 0, MCT_new->hNode->h);
        H_copy->set_priority_and_order(MCT_new->hNode->order);// Create a non-const copy
        H_copy->reordering(N, D);
        auto simulation_cost = run_completed_lacam(H_copy,MCT_NODE->agent_ratio);
      
          
        if(simulation_cost  == -1){
          MCT_backpropagate(MCT_new,1,0,1);
        }else{
          simulation_cost += MCT_new->curr_g_value;
          MCT_backpropagate(MCT_new,1,simulation_cost,0);
          minMaxStats.update(simulation_cost);
          if(min_f_cost > simulation_cost){
            min_f_cost = simulation_cost;
            solver_info(1, "found solution, cost: ", min_f_cost);
          }
          // solver_info(1, "found solution, cost: ", simulation_cost);
        }
        MCT_new ->initialize_agent_ratio(N);
        OPEN.push_back(MCT_new); 
        if(verbose == -1){
          std::cout<< " - type: generating" <<std::endl;
          std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
          std::cout<< "   pId: " 
            << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
            << std::endl;
          std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
          std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
        }
      } 
    // }
  }
  // for( int i = 1 ; i < 11; i ++){
  //   std::cout<< OPEN[i]->visits << std::endl; 
  // }

  // std::cout<< "SEARCH ENDED" <<std::endl;
  // std::cout<< minMaxStats.minimum<<std::endl;
  // std::cout<< minMaxStats.maximum <<std::endl;
  // std::cout<< "Simulation times: "<< simulation_times <<std::endl;
  // std::cout<< "Walk back times: "<< walk_back_times <<std::endl;
  // backtrack
  if (H_goal != nullptr) {
    auto H = H_goal;
    while (H != nullptr) {
      solution.push_back(H->C);
      H = H->parent;
    }
    std::reverse(solution.begin(), solution.end());
  }

  // std::cout<< "rewrite_calls: "<< rewrite_calls <<std::endl; 
  // std::cout<< "sample_times"<< sampleing_times <<std::endl; 
  // print result
  if (H_goal != nullptr && OPEN.empty()) {
    solver_info(1, "solved optimally, objective: ", objective);
  } else if (H_goal != nullptr) {
    solver_info(1, "solved sub-optimally, objective: ", objective);
  } else if (OPEN.empty()) {
    solver_info(1, "no solution");
  } else {
    solver_info(1, "timeout");
  }

  // logging
  additional_info +=
      "optimal=" + std::to_string(H_goal != nullptr && OPEN.empty()) + "\n";
  additional_info += "objective=" + std::to_string(objective) + "\n";
  additional_info += "loop_cnt=" + std::to_string(loop_cnt) + "\n";
  additional_info += "num_node_gen=" + std::to_string(EXPLORED.size()) + "\n";


  // save to tree file 
  if(tree_file != "none"){
    saveTree(tree_file);
  }

  // memory management
  for (auto a : A) delete a;
  for (auto itr : EXPLORED) delete itr.second;

  return solution;
}



int Planner::Lacam_simulator(HNode* H_init, 
  HNode*& H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, bool first_run,
   int nodes_budget,std::vector<double>& agent_ratio){
  // This function is used to simulate the search process for nodes_budget number of nodes.
  // we intend to find the path cost from H_init to H_goal. 
  HNode* min_h_node = nullptr;
  // std::cout<< "Simulation started" <<std::endl;
  // std::cout<< nodes_budget <<std::endl;
  uint path_cost = std::numeric_limits<uint>::max();
  // setup search
  auto OPEN = std::stack<HNode*>();
  auto Connection_config = std::unordered_map<Config, HNode*, ConfigHasher>();
  Connection_config[H_init->C] = H_init;
  H_init->cached_g = H_init->g; 

  OPEN.push(H_init);
  auto C_new = Config(N, nullptr);
  uint initial_g_value = H_init->g;

  int number_of_nodes_expanded = 0;
  simulation_times ++;
  while ( first_run || nodes_budget > 0) {
    nodes_budget --;
    if(OPEN.empty()){
      // std::cout<<"OPEN is empty"<<std::endl;
      break;
    }
    // do not pop here!
    auto H = OPEN.top();  // high-level node
    if(min_h_node == nullptr){
        min_h_node = H;
    }else{
      if(min_h_node->h > H->h){
        min_h_node = H ;
      }else if(min_h_node->h == H->h){
        min_h_node->f < H->f ? min_h_node = min_h_node : min_h_node = H;
      }
    }

    // low-level search end
    if (H->search_tree.empty()) {
      OPEN.pop();
      continue;
    }

    // check lower bounds
    if (H_goal != nullptr && H->f >= H_goal->f) {
      OPEN.pop();
      continue;
    }

    // if(H->node_width > 20 && !first_run){
    //   // std::cout<< "Node width is too large" <<std::endl;
    //   OPEN.pop();
    //   continue;
    // }
    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      // first time to run the simulation, we find the goal node.
      H_goal = H;
      solver_info(1, "solver find goal: cost -> ", H->g);

      std::fill(agent_ratio.begin(), agent_ratio.end(), 0);
      HNode* current = H_goal;
      while (current != H_init) {
        get_edge_cost_per_agent(agent_ratio,current->C,current->simulation_parent->C);
        current = current->simulation_parent;
      }
      compute_agent_increase_ratio(agent_ratio, H_init->C, H_goal->C);

      return H->g - initial_g_value;
    }

    // create successors at the low-level search
    auto L = H->search_tree.front();
    H->search_tree.pop();
    expand_lowlevel_tree(H, L);
    H->node_width += 1;
    // create successors at the high-level search
    const auto res = get_new_config(H, L);
    delete L;  // free
    if (!res) {
      // std::cout<< "Simulation Stucked" <<std::endl;
      continue;
    }
    // create new configuration
    for (auto a : A) C_new[a->id] = a->v_next;

    // check explored list
    const auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      // std::cout<< "Node found in the tree" <<std::endl;
      // case found
      int find_goal = rewrite_find_goal(H, iter->second, H_goal, OPEN);
      // iter->second->set_simulation_parent( H );
      // iter->second->set_simulation_parent( H );
      if(find_goal != -1){
        // find better path, return now. 

        std::fill(agent_ratio.begin(), agent_ratio.end(), 0);
        HNode* current = H_goal;
        while (current != H_init) {
          get_edge_cost_per_agent(agent_ratio,current->C,current->simulation_parent->C);
          current = current->simulation_parent;
        }
        compute_agent_increase_ratio(agent_ratio, H_init->C, H_goal->C);

        return find_goal - initial_g_value;
      }else{
        //set cached g value.
        if(Connection_config.find(C_new) != Connection_config.end()){
          if(Connection_config[C_new]->cached_g > H->g + get_edge_cost(H->C, C_new)){
            Connection_config[C_new]->cached_g = H->g + get_edge_cost(H->C, C_new);
          } 
        }else{
          Connection_config[C_new] = iter->second;
          iter->second->cached_g = H->g + get_edge_cost(H->C, C_new);
          // clean the constrain when first time see.
          iter->second->search_tree = std::queue<LNode*>();
          iter->second->search_tree.push(new LNode());
        }
      }
      if (H_goal == nullptr || iter->second->f < H_goal->f){
        update_ordering(H, iter->second, N, D);
        OPEN.push(iter->second);
      } 
    } else {
      // insert new search node
      const auto H_new = new HNode(
          C_new, D, H, H->g + get_edge_cost(H->C, C_new), get_h_value(C_new));
      H_new->setNodeID(global_node_id);
      global_node_id ++;
      EXPLORED[H_new->C] = H_new;
      H_new->set_simulation_parent( H );
      H_new->setMakeSpan(H->current_make_span + 1);
      if (H_goal == nullptr || H_new->f < H_goal->f){
        OPEN.push(H_new);
        number_of_nodes_expanded ++;
      }
    }
  }

  // get_h_value_between_config(agent_h_cost,min_h_node->C,H_init->C);
  if(min_h_node != H_init){
    // std::cout<< " Start back tracking" << std::endl;
    std::fill(agent_ratio.begin(), agent_ratio.end(), 0);
    HNode* current = min_h_node;
    while (current != H_init) {
      get_edge_cost_per_agent(agent_ratio,current->C,current->simulation_parent->C);
      current = current->simulation_parent;
    }
    compute_agent_increase_ratio(agent_ratio, H_init->C, min_h_node->C);
    // std::cout<< " End back tracking" << std::endl;  
  }
  return min_h_node->f - initial_g_value;
  // if(Connection_config.size() == 1){
  //   uint connected_path_cost = Config_A_Star_Search(Connection_config);
  //   return connected_path_cost - initial_g_value;
  // }else if(Connection_config.size() > 1){
  //   uint connected_path_cost = Config_A_Star_Search(Connection_config);
  //   return connected_path_cost - initial_g_value;
  // }
  // if(Connection_config.size() >= 1){
  //   uint connected_path_cost = Config_A_Star_Search(Connection_config);
  //   if (connected_path_cost == H_goal->g + 2*initial_g_value){
  //     walk_back_times ++;
  //   }
  //   return connected_path_cost - initial_g_value;
  // }
  // if(Connection_config.size() == 1){
  //   // only does backward connection.
  //   // best path are from h_init to h_s, and h_s to h_goal.
  //   walk_back_times++;
  //   // return H_goal->g + initial_g_value;
  //   return -1;
  // }else if(Connection_config.size() > 1){
  //   // std::cout<<"Multiple connections"<<std::endl;
  //   uint connected_path_cost = Config_A_Star_Search(Connection_config);
  //   return connected_path_cost - initial_g_value;
  // }
  // simulation failed to find the goal node.
  return -1 ;
}


uint Planner::Config_A_Star_Search(std::unordered_map<Config, HNode*, ConfigHasher>& Connection_config) {
  struct CompareHNode {
    bool operator()(HNode* const& a, HNode* const& b) const {
        return a->cached_g > b->cached_g; // Min-heap based on 'f' value
    }
  };

  // std::cout<< "A* search started" <<std::endl;
  uint min_path_cost = std::numeric_limits<uint>::max();

  // Clean tentative_g_cost
  std::fill(tentative_g_cost.begin(), tentative_g_cost.end(), std::numeric_limits<uint>::max());

  // Resize tentative_g_cost if necessary
  if (global_node_id >= tentative_g_cost.size()) {
      tentative_g_cost.resize(global_node_id + 1, std::numeric_limits<uint>::max());
  }

  // Priority queue for open set
  std::priority_queue<HNode*, std::vector<HNode*>, CompareHNode> open_set;

  // Insert all nodes in Connection_config into the priority queue
  for (auto& entry : Connection_config) {
      HNode* start_node = entry.second;
      tentative_g_cost[start_node->node_id] = start_node->cached_g;
      start_node->cached_g += start_node->h;
      open_set.push(start_node);
  }

  // Perform A* search
  while (!open_set.empty()) {
      HNode* current = open_set.top();
      open_set.pop();

      // Check if we have reached the goal
      if (is_same_config(current->C, ins->goals)) {
          return tentative_g_cost[current->node_id]; // Return the cost to reach the goal
      }

      std::vector<HNode*> successors(current->neighbor.begin(), current->neighbor.end());
      if (current->parent != nullptr) {
          successors.push_back(current->parent); // Include the parent node
      }

      for (HNode* neighbor : successors)  {
          uint tentative_g_cost_value = tentative_g_cost[current->node_id] + get_edge_cost(current, neighbor);
          // If this path to the neighbor is better, update it
          if (tentative_g_cost_value < tentative_g_cost[neighbor->node_id]) {
              tentative_g_cost[neighbor->node_id] = tentative_g_cost_value;
              neighbor->cached_g = tentative_g_cost_value + neighbor->h;
              open_set.push(neighbor);
          }
      }
  }
  // if(min_path_cost == std::numeric_limits<uint>::max()){
  //   std::cout<< "A* search failed" <<std::endl;
  // }
  // std::cout<< "A* search ended" <<std::endl;
  // If the goal is not reachable, return a large value
  return min_path_cost;
}






Solution Planner::MCT_solve(std::string& additional_info)
{

  // setup agents
  for (auto i = 0; i < N; ++i) A[i] = new Agent(i);

  rewrite_calls = 0; 
  sampleing_times = 0;
  // setup search
  auto OPEN = std::vector<MCTNode*>();
  auto EXPLORED = std::unordered_map<Config, HNode*, ConfigHasher>();
  auto H_init = new HNode(ins->starts, D, nullptr, 0, get_h_value(ins->starts));
  uint MCT_node_id = 1 ;
  auto MCT_init = new MCTNode(nullptr,H_init,MCT_node_id);
  uint roll_out_cnt =  2 * get_makespan_lower_bound(ins->starts);

  minMaxStats = MinMaxStats();
  MCT_node_id ++; 
  H_init->setMakeSpan(0);
  OPEN.push_back(MCT_init);
  EXPLORED[H_init->C] = H_init;
  std::vector<Config> solution;
  auto C_new = Config(N, nullptr);  // for new configuration
  HNode* H_goal = nullptr;          // to store goal node
  
  // std::cout<< verbose <<std::endl;
  if(verbose == -1){
    std::cout<< "version: 1.4.0\n";
    std::cout<< "events:\n";
  }
  

  if(verbose == -1){
    std::cout<< " - type: expanding" <<std::endl;
    std::cout<< "   id: "<< MCT_init->node_id <<std::endl;
    std::cout<< "   pId: " 
          << (MCT_init->parent == nullptr ? "null" : std::to_string(MCT_init->parent->node_id)) 
          << std::endl;
    std::cout<< "   visits: "<< MCT_init->visits<<std::endl;
    std::cout<< "   reward: "<< MCT_init->reward <<std::endl;
    print_utc_value(OPEN);
  }
  // for the sake of MCTS, let's force the root not have k branch factors.
  int root_branch = 10; 
  while(root_branch > 0){
    //force generate k branch factors, hopefully it will not repeate.
    auto order  = H_init->order;
    std::shuffle(order.begin(), order.end(), *MT);
    const auto res = get_next_configuration_rollout(H_init->C,  order, nullptr);
    if (!res) {
      continue;
    } 
    for (auto a : A) C_new[a->id] = a->v_next;

    // check explored list
    const auto iter = EXPLORED.find(C_new);
    if (iter != EXPLORED.end()) {
      continue;
    }else{
      const auto H_new = new HNode(
        C_new, D, H_init, H_init->g + get_edge_cost(H_init->C, C_new), get_h_value(C_new));
  
      const auto MCT_new = new MCTNode(MCT_init,H_new,MCT_node_id);
      auto results = downward_rollout_policy(H_new, H_goal, EXPLORED, 1, 1.2*(roll_out_cnt - H_new->current_make_span));
      MCT_backpropagate(MCT_new,rollout_times,results.first,results.second);
      MCT_node_id ++;
      OPEN.push_back(MCT_new);

      if(verbose == -1){
        std::cout<< " - type: generating" <<std::endl;
        std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
        std::cout<< "   pId: " 
          << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
          << std::endl;
        std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
        std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
      }
    }
    root_branch --; 
  }




  while (!is_expired(deadline)) {
    if(H_goal != nullptr){
      roll_out_cnt = H_goal->current_make_span;
    }
    loop_cnt += 1;

    // do not pop here!
    auto MCT_NODE =  MCT_selection(OPEN,H_goal);// high-level node
    if(MCT_NODE == nullptr){
      break;
    }
    if(verbose == -1){
      std::cout<< " - type: expanding" <<std::endl;
      std::cout<< "   id: "<< MCT_NODE->node_id <<std::endl;
      std::cout<< "   pId: " 
            << (MCT_NODE->parent == nullptr ? "null" : std::to_string(MCT_NODE->parent->node_id)) 
            << std::endl;
      std::cout<< "   visits: "<< MCT_NODE->visits<<std::endl;
      std::cout<< "   reward: "<< MCT_NODE->reward <<std::endl;
      print_utc_value(OPEN);
    }
    // std::cout << "Selecting node with ID: " << MCT_NODE->node_id << " ;" <<std::endl;
    auto H  = MCT_NODE->hNode; 
    // std::cout << "#Nodes in Queue: " << OPEN.size() <<" #Node ID: "<< MCT_NODE->node_id << " UTC Value: " << MCT_NODE->compute_uct_value() << std::endl ;

    if (H->search_tree.empty()) {
      MCT_NODE -> completed_node = true;
      continue;
    }
    // check lower bounds
    if (H_goal != nullptr && H->f >= H_goal->f) {
      MCT_NODE -> completed_node = true;
      continue;
    }

    // std::cout<< "level: "<<H->f <<std::endl;
    // check goal condition
    if (H_goal == nullptr && is_same_config(H->C, ins->goals)) {
      MCT_NODE -> completed_node = true;
      H_goal = H;
      solver_info(1, "found solution, cost: ", H->g);
      std::cout << "here" << std::endl; 
      if (objective == OBJ_NONE) break;
      continue;
    }


    int k = 1; 

    while (k > 0 ){
      // if(MCT_NODE->node_id == 64){
      //   std::cout<< "here" <<std::endl;
      // }
      // create successors at the low-level search
      auto L = H->search_tree.front();
      H->search_tree.pop();
      expand_lowlevel_tree(H, L);
      // create successors at the high-level search
      const auto res = get_new_config(H, L);
      delete L;  // free
      if (!res) {
        // set sample failed.
        MCT_backpropagate(MCT_NODE,1,0,1);
        k --; 
        continue;
      } 

      // create new configuration
      for (auto a : A) C_new[a->id] = a->v_next;

      // check explored list
      const auto iter = EXPLORED.find(C_new);
      if (iter != EXPLORED.end()) {
        // case found
        rewrite_rollout(H, iter->second, H_goal);
        auto H_new = (iter->second);
        if (H_goal == nullptr || H_new->f < H_goal->f){
          if(H->current_make_span + 1 < iter->second->current_make_span){
            iter->second ->setMakeSpan(H->current_make_span + 1);
          }
          const auto MCT_new = new MCTNode(MCT_NODE,H_new,MCT_node_id);
          auto results = downward_rollout_policy(H_new, H_goal, EXPLORED, 1, 1.2*(roll_out_cnt - H_new->current_make_span));
          sampleing_times += roll_out_cnt - H_new->current_make_span;
          MCT_backpropagate(MCT_new,rollout_times,results.first,results.second);
          MCT_node_id ++;
          OPEN.push_back(MCT_new);
          if(verbose == -1){
            std::cout<< " - type: generating" <<std::endl;
            std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
            std::cout<< "   pId: " 
              << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
              << std::endl;
            std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
            std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
          }
        } else{
          MCT_backpropagate(MCT_NODE,1,0,1);
        }
        // do a no insert roll out; 
      } else {
        // insert new search node
        const auto H_new = new HNode(
            C_new, D, H, H->g + get_edge_cost(H->C, C_new), get_h_value(C_new));
        H_new ->setMakeSpan(H->current_make_span + 1);
        EXPLORED[H_new->C] = H_new;
        if (H_goal == nullptr || H_new->f < H_goal->f){
          const auto MCT_new = new MCTNode(MCT_NODE,H_new,MCT_node_id);
          auto results = downward_rollout_policy(H_new, H_goal, EXPLORED, 1, 1.2*(roll_out_cnt - H_new->current_make_span));
          sampleing_times += roll_out_cnt - H_new->current_make_span;
          MCT_backpropagate(MCT_new,rollout_times,results.first,results.second);
          MCT_node_id ++;
          OPEN.push_back(MCT_new);
          if(verbose == -1 ){
            std::cout<< " - type: generating" <<std::endl;
            std::cout<< "   id: "<< MCT_new->node_id <<std::endl;
            std::cout<< "   pId: " 
              << (MCT_new->parent == nullptr ? "null" : std::to_string(MCT_new->parent->node_id)) 
              << std::endl;
            std::cout<< "   visits: "<< MCT_new->visits<<std::endl;
            std::cout<< "   reward: "<< MCT_new->reward <<std::endl;
          }
        }else{
          // record sample failed.
          MCT_backpropagate(MCT_NODE,1,0,1);
        }
      }
      k --;
    }
  }
  
  // backtrack
  if (H_goal != nullptr) {
    auto H = H_goal;
    while (H != nullptr) {
      solution.push_back(H->C);
      H = H->parent;
    }
    std::reverse(solution.begin(), solution.end());
  }

  // std::cout<< "rewrite_calls: "<< rewrite_calls <<std::endl; 
  // std::cout<< "sample_times"<< sampleing_times <<std::endl; 
  // print result
  if (H_goal != nullptr && OPEN.empty()) {
    solver_info(1, "solved optimally, objective: ", objective);
  } else if (H_goal != nullptr) {
    solver_info(1, "solved sub-optimally, objective: ", objective);
  } else if (OPEN.empty()) {
    solver_info(1, "no solution");
  } else {
    solver_info(1, "timeout");
  }

  // logging
  additional_info +=
      "optimal=" + std::to_string(H_goal != nullptr && OPEN.empty()) + "\n";
  additional_info += "objective=" + std::to_string(objective) + "\n";
  additional_info += "loop_cnt=" + std::to_string(loop_cnt) + "\n";
  additional_info += "num_node_gen=" + std::to_string(EXPLORED.size()) + "\n";


  // save to tree file 
  if(tree_file != "none"){
    saveTree(tree_file);
  }

  // memory management
  for (auto a : A) delete a;
  for (auto itr : EXPLORED) delete itr.second;

  return solution;
}


bool Planner::get_next_configuration_rollout(const Config& current_config, const std::vector<uint>& order, LNode* L){
  // setup cache
  for (auto a : A) {
    // clear previous cache
    if (a->v_now != nullptr && occupied_now[a->v_now->id] == a) {
      occupied_now[a->v_now->id] = nullptr;
    }
    if (a->v_next != nullptr) {
      occupied_next[a->v_next->id] = nullptr;
      a->v_next = nullptr;
    }

    // set occupied now
    a->v_now = current_config[a->id];
    occupied_now[a->v_now->id] = a;
  }
  // std::cout<< "Start adding constraints: " << std::endl; 
  // add constraints
  if (L != nullptr) {
    for (uint k = 0; k < L->depth; ++k) {
      const auto i = L->who[k];        // agent
      const auto l = L->where[k]->id;  // loc
      // check vertex collision
      if (occupied_next[l] != nullptr) return false;
      // check swap collision
      auto l_pre = current_config[i]->id;
      if (occupied_next[l_pre] != nullptr && occupied_now[l] != nullptr &&
          occupied_next[l_pre]->id == occupied_now[l]->id)
        return false;

      // set occupied_next
      A[i]->v_next = L->where[k];
      occupied_next[l] = A[i];
    }
  } 

  // perform PIBT
  for (auto k : order) {
    auto a = A[k];
    if (a->v_next == nullptr && !funcPIBT(a)) return false;  // planning failure
  }
  return true;
}



void Planner::rewrite_rollout(HNode* H_from, HNode* H_to, HNode* H_goal)
{
    // update neighbors
  H_from->neighbor.insert(H_to);

  // Dijkstra update
  std::queue<HNode*> Q({H_from});  // queue is sufficient
  while (!Q.empty()) {
    auto n_from = Q.front();
    Q.pop();
    for (auto n_to : n_from->neighbor) {
      auto g_val = n_from->g + get_edge_cost(n_from->C, n_to->C);
      if (g_val < n_to->g) {
        if (n_to == H_goal)
        {
          solver_info(1, "rollout cost update: ", n_to->g, " -> ", g_val);
        }
        // for (size_t i = 0; i < N; ++i) {
        //   auto v_i_from = n_from->C[i];
        //   auto v_i_to = n_to->C[i];
        //   // check connectivity
        //   if (v_i_from != v_i_to &&
        //       std::find(v_i_to->neighbor.begin(), v_i_to->neighbor.end(),
        //                 v_i_from) == v_i_to->neighbor.end()) {
            // std::cout<< "WTF !!!here !!!" <<std::endl;
        //   }
        // }
        n_to->g = g_val;
        n_to->f = n_to->g + n_to->h;
        n_to->parent = n_from;
        Q.push(n_to);
        rewrite_calls++;
      }
    }
  }
}

std::pair<double,int> Planner::downward_rollout_policy(HNode* H, HNode*& H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED,  int num_sample, int max_depth){
  // return reward and failure times ;
  // let's return the f-value as reward ! and normolize outside. 
  int original_sample = num_sample;
  int original_max_depth = max_depth;
  double total_sum_of_f_value = 0;
  int success_times = 0;
  // std::vector<std::vector<std::unordered_map<uint, uint>>> sample_array; 
  while (num_sample > 0){
    HNode* current_sample_node = H;
    max_depth = original_max_depth;
    int num_of_failure = 0;
    int num_of_depth = 0; 
    bool sample_success = true;
    // this varaible is used to record the g value of the node that is simulated;
    // in case find existing node;
    int simulated_g_value = H->g;
    // set same order here; 
    auto order  = H->order;
    // std::vector<std::unordered_map<uint, uint>> agent_to_vertex;
    // std::cout<< "Start sampling: "<< std::endl;
    while (max_depth > 0){
      // agent_to_vertex.push_back(computeCellFrequency(current_sample_node->C));
      // check goal condition
      if (is_same_config(current_sample_node->C, ins->goals)) {
        if( H_goal == nullptr || current_sample_node->f < H_goal->f){
          solver_info(1, "rollout-found solution, cost: ", current_sample_node->g);
          H_goal = current_sample_node;
        }
        break;
      }
      // create successors at the high-level search
      LNode* L = nullptr;
      if(current_sample_node == H){
        // only pop dont expand; 
        // add constraint if it is first node;
        L = H->search_tree.front();
      }
      //copy order for now: 
      // shuffle order to sampling; 
      std::shuffle(order.begin(), order.end(), *MT);
      const auto res =  get_next_configuration_rollout(current_sample_node->C,order,L);
      if (!res) {
        // config generated failed. Assume waste one depth.
        // std::cout<< "Sample Failed" << std::endl;
        sample_success = false;
        break;
      }
      // std::cout<< *current_sample_node << std::endl;
    
      // create new configuration
      auto C_new = Config(N, nullptr);  
      for (auto a : A) C_new[a->id] = a->v_next;

      const auto iter = EXPLORED.find(C_new);
      if (iter != EXPLORED.end()) {
        // std::cout<< "found existing node" << std::endl;
        if (iter->second->g  > current_sample_node->g + get_edge_cost(current_sample_node->C, C_new) ) {
          rewrite_rollout(current_sample_node, iter->second, H_goal);
          if(current_sample_node->current_make_span + 1 < iter->second->current_make_span){
            iter->second ->setMakeSpan(current_sample_node->current_make_span + 1);
          }
        }
        simulated_g_value = simulated_g_value + get_edge_cost(current_sample_node->C, C_new);
        current_sample_node = iter->second;
      } else {
        // insert new search node
        const auto H_new = new HNode(
            C_new, D, current_sample_node, current_sample_node->g + get_edge_cost(current_sample_node->C, C_new), get_h_value(C_new));
        H_new ->setMakeSpan(current_sample_node->current_make_span + 1);
        EXPLORED[H_new->C] = H_new;
        simulated_g_value = simulated_g_value + get_edge_cost(current_sample_node->C, C_new);
        current_sample_node = H_new;
      }
      max_depth --; 
      num_of_depth ++; 
    }
    if(sample_success){
      // update the reward !!!
      minMaxStats.update(simulated_g_value + current_sample_node->h);
      total_sum_of_f_value += simulated_g_value + current_sample_node->h;
      success_times++;
    }
    num_sample --;
  }
  return std::make_pair(total_sum_of_f_value/original_sample, original_sample - success_times);
}

// double Planner::downward_rollout_policy(HNode* H, HNode*& H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED,  int num_sample, int max_depth){
//   // std::cout<< "Miximual depth: "<< max_depth << std::endl; 
//   // max_depth =150;
//   int original_sample = num_sample;
//   int original_max_depth = max_depth;
//   double total_sum_of_f_value = 0;
//   uint num_of_goal_hit = 0;
//   // std::cout<<"Start sampling: "<< std::endl; 
//   double success_times = 0;
//   double sampled_f_value = 0;
//   double sum_of_f_ratio = 0;
//   double current_best_f_value = H_goal == nullptr ? N * get_makespan_lower_bound(ins->starts) : H_goal->f;
//   // std::vector<std::vector<std::unordered_map<uint, uint>>> sample_array; 
//   while (num_sample > 0){
//     HNode* current_sample_node = H;
//     max_depth = original_max_depth;
//     int num_of_failure = 0;
//     int num_of_depth = 0; 
//     bool sample_success = true;
//     // this varaible is used to record the g value of the node that is simulated;
//     // in case find existing node;
//     int simulated_g_value = H->g;
//     // set same order here; 
//     auto order  = H->order;
//     // std::vector<std::unordered_map<uint, uint>> agent_to_vertex;
//     // std::cout<< "Start sampling: "<< std::endl;
//     while (max_depth > 0){
//       // agent_to_vertex.push_back(computeCellFrequency(current_sample_node->C));
//       // check goal condition
//       if (is_same_config(current_sample_node->C, ins->goals)) {
//         if( H_goal == nullptr || current_sample_node->f < H_goal->f){
//           solver_info(1, "rollout-found solution, cost: ", current_sample_node->g);
//           H_goal = current_sample_node;
//         }
//         num_of_goal_hit++;
//         break;
//       }
//       // create successors at the high-level search
//       LNode* L = nullptr;
//       if(current_sample_node == H){
//         // only pop dont expand; 
//         // add constraint if it is first node;
//         L = H->search_tree.front();
//       }
//       //copy order for now: 
//       // shuffle order to sampling; 
//       std::shuffle(order.begin(), order.end(), *MT);
//       const auto res =  get_next_configuration_rollout(current_sample_node->C,order,L);
//       if (!res) {
//         // config generated failed. Assume waste one depth.
//         // std::cout<< "Sample Failed" << std::endl;
//         sample_success = false;
//         break;
//       }
//       // std::cout<< *current_sample_node << std::endl;
    
//       // create new configuration
//       auto C_new = Config(N, nullptr);  
//       for (auto a : A) C_new[a->id] = a->v_next;

//       // if(is_same_config(current_sample_node->C, C_new)){
//       //   std::shuffle(order.begin(), order.end(), *MT);
//       //   num_of_depth ++; 
//       //   max_depth --;
//       //   continue;
//       // }
//       // check explored list
//       const auto iter = EXPLORED.find(C_new);
//       if (iter != EXPLORED.end()) {
//         // std::cout<< "found existing node" << std::endl;
//         if (iter->second->g  > current_sample_node->g + get_edge_cost(current_sample_node->C, C_new) ) {
//           rewrite_rollout(current_sample_node, iter->second, H_goal);
//           if(current_sample_node->current_make_span + 1 < iter->second->current_make_span){
//             iter->second ->setMakeSpan(current_sample_node->current_make_span + 1);
//           }
//         }
//         simulated_g_value = simulated_g_value + get_edge_cost(current_sample_node->C, C_new);
//         current_sample_node = iter->second;
//       } else {
//         // insert new search node
//         const auto H_new = new HNode(
//             C_new, D, current_sample_node, current_sample_node->g + get_edge_cost(current_sample_node->C, C_new), get_h_value(C_new));
//         H_new ->setMakeSpan(current_sample_node->current_make_span + 1);
//         EXPLORED[H_new->C] = H_new;
//         current_sample_node = H_new;
//         simulated_g_value = simulated_g_value + get_edge_cost(current_sample_node->C, C_new);
//       }
//       max_depth --; 
//       num_of_depth ++; 
//     }
//     // sample_array.push_back(agent_to_vertex);
//     // std::cout<< "Finish sampling: " << num_of_depth<< std::endl; 
//     // std::cout<< " "<< std::endl; 
//     // std::cout<< " "<< std::endl; 
//     // std::cout<< " "<< std::endl; 
//     // std::cout<< " "<< std::endl; 

//     if(sample_success){
//       sum_of_f_ratio += 1 - (simulated_g_value + current_sample_node->h)/ current_best_f_value;
//       // std::cout<< "f-value: "<< simulated_g_value + current_sample_node->h << std::endl;
//       // std::cout<< "max-f-value: "<< current_best_f_value << std::endl;
//       // std::cout<< "f-value raito:" << std::max(1 - (simulated_g_value + current_sample_node->h)/ current_best_f_value, 0.0) << std::endl;
//     }
//     num_sample --;
//   }

//   // for ( int i = 0 ; i < 9 ; i ++){
//   //   for ( int j = i+1 ; j < 9 ; j ++){
//   //     for( int t = 0; t < 60 ; t ++){
//   //         const auto& map1 = sample_array[i][t];
//   //         const auto& map2 = sample_array[j][t];
//   //         // std::cout<<"  - type: event"<<std::endl;
//   //         // std::cout<<"    agents:"<<std::endl;
//   //         // for(auto entry: map1){
//   //         //   uint cell_index = entry.first;
//   //         //   uint x = cell_index % 32;
//   //         //   uint y = cell_index / 32;
//   //         //   std::cout<<"      - x: "<< x<<std::endl;
//   //         //   std::cout<<"        y: "<< y<<std::endl;
//   //         // }
//   //         // bool a = 0;

//   //                   std::cout<<"  - type: event"<<std::endl;
//   //         std::cout<<"    agents:"<<std::endl;
//   //         for(auto entry: map2){
//   //           uint cell_index = entry.first;
//   //           uint x = cell_index % 32;
//   //           uint y = cell_index / 32;
//   //           std::cout<<"      - x: "<< x<<std::endl;
//   //           std::cout<<"        y: "<< y<<std::endl;
//   //         }
//   //         bool b = 0;
//   //       }
//   //       return 0;
//   //     }


//   // }

//   // for ( int i = 0 ; i < 9 ; i ++){
//   //   for ( int j = i+1 ; j < 9 ; j ++){
//   //   for( int t = 1; t < 2 ; t ++){
//   //           const auto& map1 = sample_array[i][t];
//   //           const auto& map2 = sample_array[j][t];
//   //           int num_of_diff = 0;
//   //           for (const auto& entry : map1) {
//   //               uint cell_index = entry.first;
//   //               uint freq1 = entry.second;
//   //               uint freq2 = map2.count(cell_index) ? map2.at(cell_index) : 0;
//   //               if (freq1 != freq2) {
//   //                   num_of_diff ++;
//   //               } 
//   //               // if (freq1 != freq2) {
//   //               //     std::cout << "Cell Index: " << cell_index << ", Frequency in map " << i << ": " << freq1 << ", Frequency in map " << i+1 << ": " << freq2 << "\n";
//   //               // }
//   //           }
//   //           std::cout << "Number of different cells between map " << i << " and map " << j << ": " << num_of_diff << "\n";
//   //           std::cout << "time step"<< t  << "\n";
//   //         }
//   //   }
//   // }

//   // std::cout<< "Finish all sampling........... " <<std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // std::cout<< " "<< std::endl; 
//   // double f_value_ratio = 1 - total_sum_of_f_value / original_sample / max_f_value ;
//   // double f_value_ratio = 1 - total_sum_of_f_value / original_sample /  max_improvement ;
//   // std::cout<<(double)(min_f_value) <<std::endl;
//   // double f_value_ratio = 1 - total_sum_of_f_value / success_times / (N * get_makespan_lower_bound(ins->starts));
//   // if(H_goal != nullptr){
//   //   if(success_times == 0){
//   //     f_value_ratio = 0 ;
//   //   }else{
//   //     f_value_ratio = 1 - total_sum_of_f_value / success_times / H_goal->f;
//   //   }
//   // }

//   // double ratio = (H_goal->f - H->g) / (total_sum_of_f_value / original_sample - H->g);
//   // std::cout<< ratio << std::endl;
//   // std::cout<< (H_goal->f - H->g) << std::endl;
//   // std::cout<< total_sum_of_f_value / original_sample << std::endl;

//   // std::cout <<" success_ratio: " << success_ratio << " f_value_ratio: "<< f_value_ratio << std::endl;
//   // double reward = std::max((closeness_ratio + congestion_ratio + f_value_ratio)/3, goal_hit_ratio);
//   // double reward = ( (closeness_ratio + congestion_ratio + f_value_ratio)/3 + goal_hit_ratio ) / 2 ;
//   double reward =  sum_of_f_ratio;
//   // double reward =   sum_of_f_ratio / original_sample;
//   // std::cout<<"reward value: " << 2 * success_ratio + 2 * goal_ratio + g_ratio <<"  "<< std::endl;
//   //  + g_increase_ratio;
//   // std::cout<<"reward value: " << 2 * success_ratio + 2 * goal_ratio + g_ratio <<"  "<< std::endl;

//   // std::cout<<"Finishing sampling.       "<< std::endl; 
//   // std::cout<<"Reward: "<< reward<< std::endl; 
//   // std::cout<<"                          "<< std::endl; 
//   // std::cout<<"                          "<< std::endl; 

//   return reward;
// }


// int Planner::get_rollout_increase(){
//   switch (objective)
//   {
//   case OBJ_MAKESPAN:
//     /* code */
//     return 1; 
//     break;
  // case OBJ_SUM_OF_LOSS:

  //   // if (solution.empty()) return 0;
  //   // int c = 0;

  //   // for (size_t i = 0; i < N; ++i) {
  //   //   auto g = A[i];
  //   //   for (size_t t = 1; t < T; ++t) {
  //   //     if (solution[t - 1][i] != g || solution[t][i] != g) ++c;

  //   //   if (A[i]->v_now != g || A[i]->v_next != g) ++c;
  //   // }

  // return c;

  //   break; 
  
  // default:
  //   break;
  // }

// OBJ_NONE, OBJ_MAKESPAN, OBJ_SUM_OF_LOSS

// }


void Planner::rewrite_backpropagate(HNode* H_from, HNode* H_to, HNode* H_goal,
                      std::stack<HNode*>& OPEN, HNode* H_init, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED )
{
  // update neighbors
  H_from->neighbor.insert(H_to);
  // Dijkstra update
  std::queue<HNode*> Q({H_from});  // queue is sufficient
  while (!Q.empty()) {
    auto n_from = Q.front();
    Q.pop();
    for (auto n_to : n_from->neighbor) {
      auto g_val = n_from->g + get_edge_cost(n_from->C, n_to->C);
      if (g_val < n_to->g) {
        if (n_to == H_goal){
          solver_info(1, "rewrite cost update: ", n_to->g, " -> ", g_val);
          n_to->g = g_val;
          n_to->f = n_to->g + n_to->h;
          n_to->parent = n_from;
          backpropagate_order(n_to);
          increase_weight_map(n_to,true);
          pick_restart_nodes(OPEN);
          return;
        }
        n_to->g = g_val;
        n_to->f = n_to->g + n_to->h;
        n_to->parent = n_from;
        n_to->set_visit_times(node_visit_times);
        Q.push(n_to);
        if (H_goal != nullptr && n_to->f < H_goal->f) {
          OPEN.push(n_to);
        }
      }
    }
  }
}

void Planner::rewrite(HNode* H_from, HNode* H_to, HNode* H_goal,
                      std::stack<HNode*>& OPEN)
{
  // update neighbors
  H_from->neighbor.insert(H_to);

  // Dijkstra update
  std::queue<HNode*> Q({H_from});  // queue is sufficient
  while (!Q.empty()) {
    auto n_from = Q.front();
    Q.pop();
    for (auto n_to : n_from->neighbor) {
      auto g_val = n_from->g + get_edge_cost(n_from->C, n_to->C);
      if (g_val < n_to->g) {
        if (n_to == H_goal){
          solver_info(1, "cost update: ", n_to->g, " -> ", g_val);
        }
        n_to->g = g_val;
        n_to->f = n_to->g + n_to->h;
        n_to->parent = n_from;
        Q.push(n_to);
        if (H_goal != nullptr && n_to->f < H_goal->f) {
          OPEN.push(n_to);
        }
      }
    }
  }
}
void Planner::rewrite_no_push(HNode* H_from, HNode* H_to, HNode* H_goal)
{
  // update neighbors
  H_from->neighbor.insert(H_to);
  // Dijkstra update
  std::queue<HNode*> Q({H_from});  // queue is sufficient
  while (!Q.empty()) {
    auto n_from = Q.front();
    Q.pop();
    for (auto n_to : n_from->neighbor) {
      auto g_val = n_from->g + get_edge_cost(n_from->C, n_to->C);
      if (g_val < n_to->g) {
        if (n_to == H_goal){
          solver_info(1, "Simiulation find goal: ", n_to->g, " -> ", g_val);
        }
        n_to->g = g_val;
        n_to->f = n_to->g + n_to->h;
        n_to->parent = n_from;
        n_to->setMakeSpan(n_from->current_make_span + 1);
        Q.push(n_to);
      }
    }
  }
}

int Planner::rewrite_find_goal(HNode* H_from, HNode* H_to, HNode* H_goal,
                      std::stack<HNode*>& OPEN)
{
  // update neighbors
  H_from->neighbor.insert(H_to);
  
  int find_goal = -1;
  // Dijkstra update
  std::queue<HNode*> Q({H_from});  // queue is sufficient
  while (!Q.empty()) {
    auto n_from = Q.front();
    Q.pop();
    for (auto n_to : n_from->neighbor) {
      auto g_val = n_from->g + get_edge_cost(n_from->C, n_to->C);
      if (g_val < n_to->g) {
        if (n_to == H_goal){
          solver_info(1, "Simiulation find goal: ", n_to->g, " -> ", g_val);
          find_goal = H_goal->f;
        }
        n_to->g = g_val;
        n_to->f = n_to->g + n_to->h;
        n_to->parent = n_from;
        n_to->set_simulation_parent(n_from);
        n_to->setMakeSpan(n_from->current_make_span + 1);
        Q.push(n_to);
        if (H_goal != nullptr && n_to->f < H_goal->f) {
          OPEN.push(n_to);
        }
      }
    }
  }
  return find_goal;
}



void Planner::record_highlevel_node(HNode* H_node, bool expanded, bool contain_target){
  if(High_node_info.find(H_node->node_id) != High_node_info.end()) {
    bool H_expanded = High_node_info[H_node->node_id].node_expanded ? true : expanded;
    bool H_contain_target = High_node_info[H_node->node_id].contain_goal_nodes ? true : contain_target;
    High_node_info[H_node->node_id] = Search_Node_Info(*H_node,H_expanded,H_contain_target);
  } else {
    High_node_info[H_node->node_id] = Search_Node_Info(*H_node,expanded,contain_target);
  }
}


void Planner::compute_agent_increase_ratio(std::vector<double>& agent_ratio, 
    const Config& C1, const Config& C2){
    for (uint i = 0; i < N; ++i) {
      agent_ratio[i] = agent_ratio[i] / (D.get(i, C1[i]) - D.get(i, C2[i]));
      int a = D.get(i, C2[i]);
      int b = D.get(i, C1[i]);
      bool c = 0;
    }
}

void Planner::get_edge_cost_per_agent(std::vector<double>& agent_cost, 
  const Config& C1, const Config& C2)
{
  if (objective == OBJ_SUM_OF_LOSS) {
    for (uint i = 0; i < N; ++i) {
      if (C1[i] != ins->goals[i] || C2[i] != ins->goals[i]) {
        agent_cost[i] += 1;
      }
    }
  }else{
      // default: makespan
    for (uint i = 0; i < N; ++i) {
      agent_cost[i] += 1;
    }
  }
}


// void Planner::get_edge_cost_per_agent(std::vector<double>& agent_cost, const Config& C1, const Config& C2)
// {
//   if (objective == OBJ_SUM_OF_LOSS) {
//     for (uint i = 0; i < N; ++i) {
//       if (C1[i] != ins->goals[i] || C2[i] != ins->goals[i]) {
//         agent_cost[i] += 1;
//       }
//     }
//   }else{
//       // default: makespan
//     for (uint i = 0; i < N; ++i) {
//       agent_cost[i] += 1;
//     }
//   }
// }

uint Planner::get_edge_cost(const Config& C1, const Config& C2)
{
  if (objective == OBJ_SUM_OF_LOSS) {
    uint cost = 0;
    for (uint i = 0; i < N; ++i) {
      if (C1[i] != ins->goals[i] || C2[i] != ins->goals[i]) {
        cost += 1;
      }
    }
    return cost;
  }

  // default: makespan
  return 1;
}

uint Planner::get_edge_cost(HNode* H_from, HNode* H_to)
{
  return get_edge_cost(H_from->C, H_to->C);
}

uint Planner::get_h_value(const Config& C)
{
  uint cost = 0;
  if (objective == OBJ_MAKESPAN) {
    for (auto i = 0; i < N; ++i) cost = std::max(cost, D.get(i, C[i]));
  } else if (objective == OBJ_SUM_OF_LOSS) {
    for (auto i = 0; i < N; ++i) cost += D.get(i, C[i]);
  }
  return cost;
}

uint Planner::get_sum_of_distance(const Config& C)
{
  uint cost = 0;
  for (auto i = 0; i < N; ++i) cost += D.get(i, C[i]);
  return cost;
}

int Planner::get_makespan_lower_bound(const Config& C)
{
  uint c = 0;
  for (size_t i = 0; i < N; ++i) {
    c = std::max(c, D.get(i, C[i]));
  }
  return c;
}
void Planner::expand_lowlevel_tree_avoid_agent(HNode* H, LNode* L,std::unordered_set<int>& fixed_agents)
{
  // std::cout<<"Expanding node:"<<std::endl; 
  if (L->depth >= N) return;
  // std::cout<<L->depth<<std::endl;
  while(fixed_agents.find( H->order[L->depth]) != fixed_agents.end()){
    //
    L->depth++;
    if(L->depth >= N) return;
  }
  
  const auto i = H->order[L->depth];
  auto C = H->C[i]->neighbor;
  C.push_back(H->C[i]);
  // randomize
  if (MT != nullptr) std::shuffle(C.begin(), C.end(), *MT);
  // insert
  // Sometime may generate 5 actions include waiting. 
  for (auto v : C) {
    auto n  = new LNode(L, i, v);
    // std::cout<< *n << std::endl;
    H->search_tree.push(n);
  }

}
void Planner::expand_lowlevel_tree(HNode* H, LNode* L)
{
  // std::cout<<"Expanding node:"<<std::endl; 
  if (L->depth >= N) return;
  // std::cout<<L->depth<<std::endl; 
  const auto i = H->order[L->depth];
  auto C = H->C[i]->neighbor;
  C.push_back(H->C[i]);
  // randomize
  if (MT != nullptr) std::shuffle(C.begin(), C.end(), *MT);
  // insert
  // Sometime may generate 5 actions include waiting. 
  for (auto v : C) {
    auto n  = new LNode(L, i, v);
    // std::cout<< *n << std::endl;
    H->search_tree.push(n);
  }

}




// bool Planner::get_new_config_with_path_constraint(HNode* H, LNode* L)
// {
//   // setup cache
//   for (auto a : A) {
//     // clear previous cache
//     if (a->v_now != nullptr && occupied_now[a->v_now->id] == a) {
//       occupied_now[a->v_now->id] = nullptr;
//     }
//     if (a->v_next != nullptr) {
//       occupied_next[a->v_next->id] = nullptr;
//       a->v_next = nullptr;
//     }

//     // set occupied now
//     a->v_now = H->C[a->id];
//     occupied_now[a->v_now->id] = a;
//   }
//   // std::cout<< "Start adding constraints: " << std::endl; 
//   // add constraints
//   auto SIZE  = L->who.size();
//   for (uint k = 0; k < SIZE ; ++k) {
//     const auto i = L->who[k];        // agent
//     const auto l = L->where[k]->id;  // loc
//     // check vertex collision
//     if (occupied_next[l] != nullptr) return false;
//     // check swap collision
//     auto l_pre = H->C[i]->id;
//     if (occupied_next[l_pre] != nullptr && occupied_now[l] != nullptr &&
//         occupied_next[l_pre]->id == occupied_now[l]->id)
//       return false;

//     // set occupied_next
//     A[i]->v_next = L->where[k];
//     occupied_next[l] = A[i];
//   }


//   // perform PIBT
//   for (auto k : H->order) {
//     auto a = A[k];
//     if (a->v_next == nullptr && !funcPIBT(a)) {
//       return false;  // planning failure
//     }
//   }

//   return true;
// }


bool Planner::get_new_config(HNode* H, LNode* L)
{
  // setup cache
  for (auto a : A) {
    // clear previous cache
    if (a->v_now != nullptr && occupied_now[a->v_now->id] == a) {
      occupied_now[a->v_now->id] = nullptr;
    }
    if (a->v_next != nullptr) {
      occupied_next[a->v_next->id] = nullptr;
      a->v_next = nullptr;
    }

    // set occupied now
    a->v_now = H->C[a->id];
    occupied_now[a->v_now->id] = a;
  }
  // std::cout<< "Start adding constraints: " << std::endl; 
  // add constraints
  // if(L->depth != L->who.size()){
  //   std::cout<< "Error: L->depth != L->who.size() "<< std::endl; 
  // }
  for (uint k = 0; k < L->who.size(); ++k) {
    const auto i = L->who[k];        // agent
    const auto l = L->where[k]->id;  // loc
    // check vertex collision
    if (occupied_next[l] != nullptr) {
      // std::cout<<" fall into vertex collision: "<< std::endl;
      return false;
    }
    // check swap collision
    auto l_pre = H->C[i]->id;
    if (occupied_next[l_pre] != nullptr && occupied_now[l] != nullptr &&
        occupied_next[l_pre]->id == occupied_now[l]->id){
      // std::cout<<" fall into swap collision: "<< std::endl;
      return false;
        }

    // set occupied_next
    A[i]->v_next = L->where[k];
    occupied_next[l] = A[i];
  }


  // perform PIBT
  for (auto k : H->order) {
    auto a = A[k];
    if (a->v_next == nullptr && !funcPIBT(a)) {
      return false;  // planning failure
    }
  }

  return true;
}

int Planner::get_action(Vertex* const& v_curr, Vertex* const& v_next, int width) {
    int x1 = v_curr->index % width;
    int y1 = v_curr->index / width;
    int x2 = v_next->index % width;
    int y2 = v_next->index / width;

    if (x1 == x2 && y1 == y2) return 0;         // WAIT
    if (x2 == x1 + 1 && y2 == y1) return 1;     // RIGHT
    if (x2 == x1 - 1 && y2 == y1) return 2;     // LEFT
    if (x2 == x1 && y2 == y1 + 1) return 3;     // DOWN
    if (x2 == x1 && y2 == y1 - 1) return 4;     // UP
}

void  Planner::propogate_q_value_k_steps(int agent_id, Vertex* const& v_curr, Vertex* const& v_next, double cost_to_go, int k_steps)
{ 
  double alpha = 0.2;
  if(k_steps <= 0){
    // std::cout<<"Error: k_steps should be greater than 0"<<std::endl;
    return;
  }
  // cost_to_go -> from v_curr to v_next;
  // for (auto &&m : v_curr->neighbor) {
  //   if( m->id == v_next->id) continue; // skip the next vertex
  //   double cost =  cost_to_go + 2 ;
  //   int action = get_action(v_curr, m, ins->G.width);
  //   double q_table = Q_tables[agent_id](v_curr->id, action);
  //   double delta = cost - Q_tables[agent_id](v_curr->id, action);
  //   Q_tables[agent_id].update(v_curr->id, action, alpha, delta);
  // }
  int action = get_action(v_curr, v_next, ins->G.width);
  // update the Q value for the action to v_next
  double delta = cost_to_go - get_q_value(agent_id, v_curr, v_next);
  Q_tables[agent_id].update(v_curr->index, action, alpha, delta);    

  // update the Q value for the action WAIT
  delta = cost_to_go + 1 - get_q_value(agent_id, v_curr, v_curr);
  Q_tables[agent_id].update(v_curr->index, 0, alpha, delta);    

  for (auto &&m : v_curr->neighbor) {
    if( m->id == v_next->id) continue; // skip the next vertex
    propogate_q_value_k_steps(agent_id, m, v_curr, cost_to_go + 1, k_steps - 1);
  }
}

// void Planner::update_q_value(int agent_id, Vertex* const& v_curr, Vertex* const& v_next, double cost_to_go){
//   int action = get_action(v_curr, v_next, ins->G.width);
//   if(Q_tables[agent_id](v_curr->id,action) == -1){
//     // not found, this cannot happened actually. 
//     // std::cout<<"Error: Q_table not found for agent: "<< agent_id << " at vertex: "<< v_curr->index << " action: "<< action << std::endl;
//     Q_tables[agent_id].set(v_curr->id, action, cost_to_go);
//   }else{
//     // found
//     double q_table = Q_tables[agent_id](v_curr->id, action);
//     double delta = cost_to_go - Q_tables[agent_id](v_curr->id, action);
//     // if(delta < 0){
//     //   std::cout<<"Error: Q_table update with negative delta: "<< delta << " for agent: "<< agent_id << " at vertex: "<< v_curr->index << " action: "<< action << std::endl;
//     //   // delta = 0; // ignore negative delta
//     // }
//     double alpha = 0.4; // learning rate
//     Q_tables[agent_id].update(v_curr->id, action, alpha, delta);
//   }
// }

double Planner::get_q_value(int agent_id, Vertex* const& v_curr, Vertex* const& v_next){
  int action = get_action( v_curr, v_next, ins->G.width);
  if(Q_tables[agent_id](v_curr->index,action) == -1){
    // not found
    double q_value = PIBT_D.get_heuristic(agent_id, v_next);
    Q_tables[agent_id].set(v_curr->index, action, q_value);
    return q_value; 
  }else{ 
    return Q_tables[agent_id](v_curr->index,action);
  }
}


bool Planner::funcPIBT(Agent* ai)
{
  const auto i = ai->id;
  const auto K = ai->v_now->neighbor.size();

  // get candidates for next locations
  for (auto k = 0; k < K; ++k) {
    auto u = ai->v_now->neighbor[k];
    C_next[i][k] = u;
    if (MT != nullptr)
      tie_breakers[u->id] = get_random_float(MT);  // set tie-breaker
  }
  C_next[i][K] = ai->v_now;

  // sort
  // std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
  //           [&](Vertex* const v, Vertex* const u) {
  //             return D.get(i, v) + tie_breakers[v->id] <
  //                    D.get(i, u) + tie_breakers[u->id];
  //           });
  // PIBT_D.compare_edge_map();
  
  // std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
  //           [&](Vertex* const v, Vertex* const u) {
  //             // double b = D.get(i, u);
  //             // double e = D.get(i, v); 
  //             // double f = D.get(i, u);
  //             // return (double)D.get(i, v) + tie_breakers[v->id] <
  //             // (double) D.get(i, u) + tie_breakers[u->id];
  //             return PIBT_D.get_heuristic(i, v) + tie_breakers[v->id] <
  //             PIBT_D.get_heuristic(i, u) + tie_breakers[u->id];
  //             // return PIBT_D.get_individual_heuristic(i, v) + tie_breakers[v->id] <
  //             // PIBT_D.get_individual_heuristic(i, u) + tie_breakers[u->id];
  //           });
  // double selected_a = 10000; 
  // double selected_b = 10000;
  // Vertex* v_selected_a;
  // Vertex* v_selected_b;
  // for (auto k = 0; k < K + 1; ++k) {
  //   auto v = C_next[i][k];
  //   double a = get_q_value(i, ai->v_now,v);
  //   double b = PIBT_D.get_heuristic(i, v); 
  //   std::cout<< a << " " << b << std::endl;
  //   if(a < selected_a){
  //     selected_a = a;
  //     v_selected_a = v;
  //   }
  //   if(b < selected_b){
  //     selected_b = b;
  //     v_selected_b = v;
  //   }
  // }
  // if(v_selected_a != v_selected_b){
  //   std::cout<<"Error: Q value not equal to PIBT heuristic: "<< selected_a << " != " << selected_b << std::endl;
  //   std::cout<<"Agent: "<< i << " Vertex: "<< ai->v_now->id << " Next Vertex A: "<< v_selected_a->id << " Next Vertex B: "<< v_selected_b->id << std::endl;
  // }


  // std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
  //         [&](Vertex* const v, Vertex* const u) {
  //           return get_q_value(i, ai->v_now,v) + tie_breakers[v->id] <
  //           get_q_value(i, ai->v_now,u) + tie_breakers[u->id];
  //         });

  // double selected_a = 10000; 
  // double selected_b = 10000;
  // Vertex* v_selected_a;
  // Vertex* v_selected_b;
  // std::cout<< "Calculating heuritisic:              "<< std::endl;
  // for (auto k = 0; k < K + 1; ++k) {
  //   auto v = C_next[i][k];
  //   double a = get_q_value(i, ai->v_now,v);
  //   double b = PIBT_D.get_heuristic(i, v); 
  //   std::cout<< "Q value:" << a << " PIBT heuristic: " << b << std::endl;
  //   if(a < selected_a){
  //     selected_a = a;
  //     v_selected_a = v;
  //   }
  //   if(b < selected_b){
  //     selected_b = b;
  //     v_selected_b = v;
  //   }
  // }
  // if(v_selected_a != v_selected_b){
  //   std::cout<<"Error: Q value not equal to PIBT heuristic: "<< selected_a << " != " << selected_b << std::endl;
  //   std::cout<<"Agent: "<< i << " Vertex: "<< ai->v_now->id << " Next Vertex A: "<< v_selected_a->id << " Next Vertex B: "<< v_selected_b->id << std::endl;
  // }

  // std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
  //         [&](Vertex* const v, Vertex* const u) {
  //           return get_q_value(i, ai->v_now,v) + tie_breakers[v->id] <
  //           get_q_value(i, ai->v_now,u) + tie_breakers[u->id];
  //         });

  // std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
  //         [&](Vertex* const v, Vertex* const u) {
  //           return PIBT_D.get_heuristic(i, v) + tie_breakers[v->id] <
  //           PIBT_D.get_heuristic(i, u) + tie_breakers[u->id];
  //         });


  if(!guidance_heuristic.initialized){
    std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
              [&](Vertex* const v, Vertex* const u) {
                return D.get(i, v) + tie_breakers[v->id] <
                D.get(i, u) + tie_breakers[u->id];
              });
  }else{
    // std::cout<< "Using guidance heuristic to sort actions for agent: "<< i << std::endl;
      std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
            [&](Vertex* const v, Vertex* const u) {
              return guidance_heuristic.get_Astar_heuristic(i, v->id) + tie_breakers[v->id] <
              guidance_heuristic.get_Astar_heuristic(i, u->id) + tie_breakers[u->id];
            });

      // std::cout<< "hererererer"<<std::endl;
      // std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
      // [&](Vertex* const v, Vertex* const u) {
      //   return guidance_heuristic.get_heuristic(i, v->index) + tie_breakers[v->id] <
      //   guidance_heuristic.get_heuristic(i, u->index) + tie_breakers[u->id];
      // });
  }

  // if(!guidance_heuristic.initialized){
  //   // update the guidance heuristic

  // }
  // std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
  //         [&](Vertex* const v, Vertex* const u) {
  //           return D.get(i, v) + tie_breakers[v->id] <
  //           D.get(i, u) + tie_breakers[u->id];
  //         });
  action_history[ai->id][ai->v_now->index] = C_next[i];  // record the action order;
  

  Agent* swap_agent = swap_possible_and_required(ai);
  if (swap_agent != nullptr)
    std::reverse(C_next[i].begin(), C_next[i].begin() + K + 1);

  // main operation
  for (auto k = 0; k < K + 1; ++k) {
    auto u = C_next[i][k];

    // avoid vertex conflicts
    if (occupied_next[u->id] != nullptr) continue;

    auto& ak = occupied_now[u->id];

    // avoid swap conflicts
    if (ak != nullptr && ak->v_next == ai->v_now) continue;

    // reserve next location
    occupied_next[u->id] = ai;
    ai->v_next = u;

    // priority inheritance
    if (ak != nullptr && ak != ai && ak->v_next == nullptr && !funcPIBT(ak))
      continue;

    // success to plan next one step
    // pull swap_agent when applicable
    if (k == 0 && swap_agent != nullptr && swap_agent->v_next == nullptr &&
        occupied_next[ai->v_now->id] == nullptr) {
      swap_agent->v_next = ai->v_now;
      occupied_next[swap_agent->v_next->id] = swap_agent;
    }
    return true;
  }

  // failed to secure node
  occupied_next[ai->v_now->id] = ai;
  ai->v_next = ai->v_now;
  return false;
}

Agent* Planner::swap_possible_and_required(Agent* ai)
{
  const auto i = ai->id;
  // ai wanna stay at v_now -> no need to swap
  if (C_next[i][0] == ai->v_now) return nullptr;

  // usual swap situation, c.f., case-a, b
  auto aj = occupied_now[C_next[i][0]->id];
  if (aj != nullptr && aj->v_next == nullptr &&
      is_swap_required(ai->id, aj->id, ai->v_now, aj->v_now) &&
      is_swap_possible(aj->v_now, ai->v_now)) {
    return aj;
  }

  // for clear operation, c.f., case-c
  for (auto u : ai->v_now->neighbor) {
    auto ak = occupied_now[u->id];
    if (ak == nullptr || C_next[i][0] == ak->v_now) continue;
    if (is_swap_required(ak->id, ai->id, ai->v_now, C_next[i][0]) &&
        is_swap_possible(C_next[i][0], ai->v_now)) {
      return ak;
    }
  }

  return nullptr;
}

// simulate whether the swap is required
bool Planner::is_swap_required(const uint pusher, const uint puller,
                               Vertex* v_pusher_origin, Vertex* v_puller_origin)
{
  auto v_pusher = v_pusher_origin;
  auto v_puller = v_puller_origin;
  Vertex* tmp = nullptr;
  while (D.get(pusher, v_puller) < D.get(pusher, v_pusher)) {
    auto n = v_puller->neighbor.size();
    // remove agents who need not to move
    for (auto u : v_puller->neighbor) {
      auto a = occupied_now[u->id];
      if (u == v_pusher ||
          (u->neighbor.size() == 1 && a != nullptr && ins->goals[a->id] == u)) {
        --n;
      } else {
        tmp = u;
      }
    }
    if (n >= 2) return false;  // able to swap
    if (n <= 0) break;
    v_pusher = v_puller;
    v_puller = tmp;
  }

  // judge based on distance
  return (D.get(puller, v_pusher) < D.get(puller, v_puller)) &&
         (D.get(pusher, v_pusher) == 0 ||
          D.get(pusher, v_puller) < D.get(pusher, v_pusher));
}

// simulate whether the swap is possible
bool Planner::is_swap_possible(Vertex* v_pusher_origin, Vertex* v_puller_origin)
{
  auto v_pusher = v_pusher_origin;
  auto v_puller = v_puller_origin;
  Vertex* tmp = nullptr;
  while (v_puller != v_pusher_origin) {  // avoid loop
    auto n = v_puller->neighbor.size();  // count #(possible locations) to pull
    for (auto u : v_puller->neighbor) {
      auto a = occupied_now[u->id];
      if (u == v_pusher ||
          (u->neighbor.size() == 1 && a != nullptr && ins->goals[a->id] == u)) {
        --n;      // pull-impossible with u
      } else {
        tmp = u;  // pull-possible with u
      }
    }
    if (n >= 2) return true;  // able to swap
    if (n <= 0) return false;
    v_pusher = v_puller;
    v_puller = tmp;
  }
  return false;
}




// void Planner::saveTree(const std::string &fileName) const // write the CT to a file
// {
// 	std::ofstream output;
// 	output.open(fileName, std::ios::out);
// 	output << "digraph G {" << std::endl;
// 	output << "size = \"5,5\";" << std::endl;
// 	output << "center = true;" << std::endl;
// 	for (auto node : High_node_info)
// 	{
// 		output << node.second.node_id << " [label=\"#" << node.second.node_id
// 					<< "\ng+h="<< node.second.g << "+" << node.second.h
// 					<< "\nd=" << node.second.node_expanded << "\"]" << std::endl;
// 		// if (node == dummy_start)
// 		// 	continue;
// 		output << node.second.parent_id << " -> " << node.second.node_id 
//     // << " [label=\"";
// 		// for (const auto &constraint : node->constraints)
// 		// 	output << constraint;
// 		// output << "\nAgents ";
//     //     for (const auto &path : node->paths)
//     //         output << path.first << "(+" << path.second.size() - paths_found_initially[path.first].size() << ") ";
//     //     output << "\"]"
//     << std::endl;
// 	}
//   for (auto node : High_node_info){
//     if(node.second.contain_goal_nodes){
//       output << node.second.node_id << " [color=red]" << std::endl;
//     }
//   }
// 	// auto node = goal_node;
// 	// while (node != nullptr)
// 	// {
// 	// 	output << node->time_generated << " [color=red]" << std::endl;
// 	// 	node = node->parent;
// 	// }
// 	output << "}" << std::endl;
// 	output.close();
// }


// void Planner::get_node_generated() const // write the CT to a file
// {
//   std::vector<Search_Node_Info> node_info = std::vector<Search_Node_Info>(0); 
//   for (auto n : High_node_info){
//     node_info.push_back(n.second);
//   }

//   std::sort(node_info.begin(), node_info.end(), [](const Search_Node_Info& n1, const Search_Node_Info& n2) {
//       return n1.node_id < n2.node_id;
//   });

//   uint node_generate = 0; 
//   uint target_found = 0; 
//   for (auto node : node_info)
// 	{
//     node_generate ++ ; 
//     if(node.contain_goal_nodes){
//       target_found ++; 
//     }
// 	}
//   // std::cout<< "#Node Expanded: "<< node_generate << " #Target Found: " << target_found<< std::endl;
// 	output.close();
// }



void Planner::saveTree(const std::string &fileName) const // write the CT to a file
{
	std::ofstream output;
	output.open(fileName, std::ios::out);
  output<< "version: 1.4.0\n";
  output<< "events:\n";

  std::vector<Search_Node_Info> node_info = std::vector<Search_Node_Info>(0); 
  for (auto n : High_node_info){
    node_info.push_back(n.second);
  }

  std::sort(node_info.begin(), node_info.end(), [](const Search_Node_Info& n1, const Search_Node_Info& n2) {
      return n1.node_id < n2.node_id;
  });

  uint node_generate = 0; 
  uint target_found = 0; 
  for (auto node : node_info)
	{
    std::string parentId = (node.parent_id == -1) ? "1" : std::to_string(node.parent_id+1);
    output << "  - { type: decision, id: " 
              << std::to_string(node.node_id+1) << ", pId: " 
              << parentId << ", f: "<< node.f << ", g: "<< node.g 
              << ", target: " << (node.contain_goal_nodes ? "true" : "false") << " }\n";
    node_generate ++ ; 
    if(node.contain_goal_nodes){
      target_found ++; 
    }
	}
  // std::cout<< "#Node Expanded: "<< node_generate << " #Target Found: " << target_found<< std::endl;
	output.close();
}


std::ostream& operator<<(std::ostream& os, const Objective obj)
{
  if (obj == OBJ_NONE) {
    os << "none";
  } else if (obj == OBJ_MAKESPAN) {
    os << "makespan";
  } else if (obj == OBJ_SUM_OF_LOSS) {
    os << "sum_of_loss";
  }
  return os;
}


std::ostream& operator<<(std::ostream& os, const HNode& node) {
    os << "HNode { ";
    os << "Node ID: "<< node.node_id << ", ";
    os << "Config: [";
    for (size_t i = 0; i < 20; ++i) {
        if (node.C[i]) {
            os << "(ID: " << node.C[i]->id << ", Index: " << node.C[i]->index << ")";
        } else {
            os << "Null";
        }
        if (i < node.C.size() - 1) os << ", ";
    }
    os << "], ";
    os << "g: " << node.g << ", ";
    os << "h: " << node.h << ", ";
    os << "f: " << node.f << ", ";
    os << "Parent: " << (node.parent ? node.parent->node_id : -1) << ", ";
    os << "CT-tree: " << node.search_tree.size()<< ", ";
    // os << "Neighbors: { ";
    // for (const auto& n : node.neighbor) {
    //     os << n->C << " "; 
    // }
    // os << "}, ";
    // os << "Priorities: [ ";
    // for (const auto& p : node.priorities) {
    //     os << p << " ";
    // }
    // os << "], ";
    // os << "Order: [ ";
    // for (const auto& o : node.order) {
    //     os << o << " ";
    // }
    // os << "]";
    os << "}";
    os << " ";
    os << " ";
    os << " ";
    return os;
}

std::ostream& operator<<(std::ostream& os, const LNode& node) {
    os << "LNode Details:\n";
    os << "Depth: " << node.depth << "\n";

    os << "Who: [";
    for (size_t i = 0; i < node.who.size(); ++i) {
        os << node.who[i];
        if (i < node.who.size() - 1) os << ", ";
    }
    os << "]\n";

    os << "Where: [";
    for (size_t i = 0; i < node.where.size(); ++i) {
        os << node.where[i]; // Assuming Vertex has an `id` or similar property
        if (i < node.where.size() - 1) os << ", ";
    }
    os << "]\n";

    return os;
}