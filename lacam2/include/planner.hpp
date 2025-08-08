/*
 * lacam-star
 */

 #pragma once

 #include "dist_table.hpp"
 #include "AstarDist_table.hpp"
 #include "graph.hpp"
 #include "instance.hpp"
 #include "utils.hpp"
 #include "min_max_stats.hpp"
 #include <unordered_set>
 #include "q_table.hpp"
 // objective function
 enum Objective { OBJ_NONE, OBJ_MAKESPAN, OBJ_SUM_OF_LOSS };
 std::ostream& operator<<(std::ostream& os, const Objective objective);
 
 // PIBT agent
 struct Agent {
   const uint id;
   Vertex* v_now;   // current location
   Vertex* v_next;  // next location
   Agent(uint _id) : id(_id), v_now(nullptr), v_next(nullptr) {}
 };
 using Agents = std::vector<Agent*>;
 
 // low-level node
 struct LNode {
   std::vector<uint> who;
   Vertices where;
   uint depth;
   LNode(LNode* parent = nullptr, uint i = 0,
         Vertex* v = nullptr);  // who and where
 };
 std::ostream& operator<<(std::ostream& os, const LNode& node);
 
 // high-level node
 struct HNode {
   static uint HNODE_CNT;  // count #(high-level node)
   const Config C;
 
   // tree
   HNode* parent;
   std::set<HNode*> neighbor;
 
   // costs
   uint g;        // g-value (might be updated)
   const uint h;  // h-value
   uint f;        // g + h (might be updated)
   uint current_make_span = 0;
 
   uint cached_g = 0;
   // for low-level search
   std::vector<float> priorities;
   std::vector<uint> order;
   std::queue<LNode*> search_tree;
 
   uint node_id; 
   uint node_width = 0;
   uint simulation_cost = 0;
   uint order_updated = 0;
   uint visit_times = 0; 
   HNode* simulation_parent;
   HNode(const Config& _C, DistTable& D, HNode* _parent, const uint _g,
         const uint _h);
   ~HNode();
   
  void set_simulation_parent( HNode* s_parent ){
    simulation_parent = s_parent;
  }

  void set_visit_times(int v_times){
    visit_times = v_times;
  }

  void setNodeID(uint id) {
    node_id = id;
  }


  void set_priority_and_order(const std::vector<uint>& input_order) {
    // Set the order directly from the input
    order = input_order;

    // Update priorities based on the new order
    for (size_t i = 0; i < order.size(); ++i) {
        priorities[order[i]] = (float) (order.size() - i)/ order.size(); // Higher priority for earlier indices
    }
  }

  void reordering(size_t N, DistTable& D) {
    for (size_t i = 0; i < N; ++i) {
      if (D.get(i, C[i])== 0) {
        priorities[i] = priorities[i] - (int)priorities[i];
      }else{
        priorities[i] = priorities[i] + 1;
      }
    }
    // set order
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](uint i, uint j) { return priorities[i] > priorities[j]; });
  }

   void setMakeSpan(uint _current_make_span) {
     current_make_span = _current_make_span;
   }
 
 };
 using HNodes = std::vector<HNode*>;
 std::ostream& operator<<(std::ostream& os, const HNode& node);
 
 
 
 const double C = 1.41;
 
 struct MCTNode {
     MCTNode* parent;
     HNode* hNode;
     int visits = 0;
     int failuare_times = 0;
     double reward = 0; 
     bool completed_node = false;
     uint curr_g_value = 0;
     bool first_branch = true;
     uint node_id; 
     std::vector<double> agent_ratio;
     std::unordered_map<Config, HNode*, ConfigHasher> generated_configs; // Track generated configurations
     uint best_simulation_cost = std::numeric_limits<uint>::max();

     MCTNode(MCTNode* _parent, HNode* _HNode, uint _node_id);
     ~MCTNode();

     void setCurrentGValue(uint g_value) {
       curr_g_value = g_value;
     }
   
     void initialize_agent_ratio(size_t N){
        agent_ratio.resize(N); 
     }
     // double MCTNode::compute_uct_value(HNode* goal_node) {
     //   double exploitation = reward / (visits + 1e-6);
     //   // std::cout<< " "<< std::endl; 
     //   // std::cout<< std::to_string(parent == nullptr ? 0 : parent->visits) << std::endl;
     //   // std::cout<< visits << std::endl; 
     //   double exploration = C * std::sqrt(std::log( (parent == nullptr ? 0 : parent->visits) + 1) / (visits + 1e-6) );
     //   // std::cout<< "Exploitation Score: " << exploitation << "; Exploration Score: "<< exploration << "; UTC value: "<<exploitation + exploration <<std::endl;
     //   return exploitation + exploration;
     // }
     // double MCTNode::print_uct_value(HNode* goal_node) {
     //   double exploitation = reward / (visits + 1e-6);
     
     //   double exploration = C * std::sqrt(std::log( (parent == nullptr ? 0 : parent->visits) + 1) / (visits + 1e-6) );
     //   // std::cout<< "Exploitation Score: " << exploitation << "; Exploration Score: "<< exploration << "; UTC value: "<<exploitation + exploration <<std::endl;
     //   std::cout<< "      "<< "reward value: "<< reward <<std::endl;
     //   std::cout<< "      "<< "exploitation: "<< exploitation <<std::endl;
     //   std::cout<< "      "<< "exploration: "<< exploration <<std::endl;
     //   std::cout<< "      "<< "utc value: "<< exploitation + exploration <<std::endl;
     //   return exploitation + exploration;
     // }
     
 };
 
 
 
 struct Search_Node_Info{
   uint g = 0;
   uint h = 0;
   uint f = 0;
   int node_id = 0;
   int parent_id = -1;
   bool node_expanded = false;
   bool contain_goal_nodes = false;
 
     // Default constructor
   Search_Node_Info() = default;
   Search_Node_Info(const HNode& hNode, bool expanded, bool is_goal)
       : g(hNode.g),
         h(hNode.h),
         f(hNode.f),
         node_id(int(hNode.node_id)),
         parent_id(hNode.parent == nullptr ? -1 : int(hNode.parent->node_id)),
         node_expanded(expanded),
         contain_goal_nodes(is_goal)
   {
     // if(hNode.parent != nullptr) std::cout<< hNode.parent->node_id << std::endl;
     // std::cout<< hNode.node_id << std::endl;
   }
 
   friend std::ostream& operator<<(std::ostream& os, const Search_Node_Info& info) {
     os << "Search_Node_Info:\n"
        << "  g: " << info.g << "\n"
        << "  h: " << info.h << "\n"
        << "  f: " << info.f << "\n"
        << "  node_id: " << info.node_id << "\n"
        << "  parent_id: " << info.parent_id << "\n"
        << "  node_expanded: " << (info.node_expanded ? "true" : "false") << "\n"
        << "  contain_goal_nodes: " << (info.contain_goal_nodes ? "true" : "false") << "\n";
     return os;
   }
 };
 
 
  struct PairHash {
    template <typename T1, typename T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ (std::hash<T2>()(pair.second) << 1);
    }
  };

 struct Planner {
   const Instance* ins;
   const Deadline* deadline;
   std::mt19937* MT;
   const int verbose;
 
   // hyper parameters
   const Objective objective;
   const float RESTART_RATE;  // random restart
 
   // solver utils
   const uint N;       // number of agents
   const uint V_size;  // number o vertices
   DistTable D;
   AstarDistTable PIBT_D;

   std::vector<QTable> Q_tables; // Q-tables for each agent

   uint loop_cnt;      // auxiliary
   uint rewrite_calls;   // auxiliary;
   uint sampleing_times;   // auxiliary;
   // used in PIBT
   std::vector<std::array<Vertex*, 5> > C_next;  // next locations, used in PIBT
   std::vector<float> tie_breakers;              // random values, used in PIBT
   Agents A;
   Agents occupied_now;                          // for quick collision checking
   Agents occupied_next;                         // for quick collision checking
 
   MinMaxStats minMaxStats;
   int rollout_times = 1; 
   std::unordered_map<uint, Search_Node_Info> High_node_info;
   const std::string tree_file; 
   uint global_node_id = 0;
   uint global_MCT_node_id = 1 ;
   uint walk_back_times = 0 ;
   uint simulation_times = 0 ;
   uint node_generation_times = 0;

   std::vector<uint> global_order;
   std::vector<unsigned int> tentative_g_cost = std::vector<unsigned int>(10000, 0);
   uint min_f_cost = std::numeric_limits<uint>::max();
   MCTNode* best_mct_node = nullptr;
   uint order_updated_times = 0;
   uint node_visit_times = 0; 
   double decay_factor = 0.9;
   double decay_counter = 1; 
   std::vector<std::vector<std::pair<std::pair<int, int>, double>>> punishment_vector;
   std::vector<std::unordered_map<std::pair<int, int>, double, PairHash>> individual_transition_frequency;
   std::vector<HNode*> SOLUTION_NODES;
   
   bool reset_congestion_map = false;
   Planner(const Instance* _ins, const Deadline* _deadline, std::mt19937* _MT,
           const int _verbose = 0,
           // other parameters
           const Objective _objective = OBJ_NONE,
           const float _restart_rate = 0.001, 
           const std::string save_tree_file = "none"
         );
 
  ~Planner();
  Solution solve(std::string& additional_info);
  Solution MCT_solve(std::string& additional_info);
  void rewrite_rollout(HNode* H_from, HNode* H_to, HNode* H_goal);
  void expand_lowlevel_tree_avoid_agent(HNode* H, LNode* L,std::unordered_set<int>& fixed_agents);
  std::pair<double,int> downward_rollout_policy(HNode* H, HNode*& H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, int num_sample, int max_depth);
  bool get_next_configuration_rollout(const Config& current_config, const std::vector<uint>& order, LNode* L);
  MCTNode* MCT_selection(std::vector<MCTNode*>& node_pool,HNode* goal_node);
  void MCT_backpropagate(MCTNode* mct_node, int sample_times,  double reward, int failuare_times);
  int get_rollout_increase();
  int get_makespan_lower_bound(const Config& C);
  double compute_uct_value(MCTNode* mcts_node);
  void print_utc_value(std::vector<MCTNode*>& node_pool);
  int Lacam_simulator(HNode* H_init, HNode*& H_goal, 
  std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, 
  bool first_run, int max_depth,std::vector<double>& agent_ratio);
  int rewrite_find_goal(HNode* H_from, HNode* H_to, HNode* H_goal, std::stack<HNode*>& OPEN);
  uint Config_A_Star_Search(std::unordered_map<Config, HNode*, ConfigHasher>& Connection_config) ;
  uint multiple_simulations(int num_simulations, HNode*& H_goal, MCTNode* mcts_node, 
    std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, std::vector<Config>& solution); 
  
  Solution MCT_lacam(std::string& additional_info);
  Solution MCT_vanilla_lacam(std::string& additional_info);
  Solution MCT_multiple_lacam(std::string& additional_info);
  Solution backpropagate_solve(std::string& additional_info);
  
  void export_interacted_agents_graph(const std::vector<std::set<int>>& interacted_agents, const std::string& filename);
  
  void export_cluster_solution(
    const std::vector<std::vector<int>>& solution_nodes,
    const std::vector<std::vector<int>>& depth_clustered_agents,
    const std::string& filename);

  std::vector<std::set<int>> transitiveClosureAll(const std::vector<std::set<int>>& graph);

  std::vector<int> depth_cluster(const std::vector<std::set<int>>& graph, int start, int input_depth, std::vector<bool>& visited);
  void build_dependence_graph(HNode* input_H_goal);
  void backpropagate_order(HNode* H_goal);
  void decay_punishment();

  HNode* copy_node(HNode* original_node); 
  MCTNode* MCT_select_successors_from_pool(MCTNode* mcts_node,  
    std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, 
    std::vector<Config>& solution, HNode*& H_goal);

  MCTNode* MCT_random_successor_generator(Config& C_new, MCTNode* mcts_node, 
  std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, HNode* H_goal);
  MCTNode* MCT_Learn_order_to_branch(Config& C_new, MCTNode* mcts_node, 
    std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, HNode* H_goal);
  MCTNode* MCT_Learn_order_to_branch(Config& C_new, MCTNode* mcts_node);
  void propagate_order_to_neighbors(HNode* current_node);
  MCTNode* MCT_selection(std::vector<MCTNode*>& node_pool);
  uint run_completed_lacam(HNode* H_init,std::vector<double>& agent_ratio);
  uint run_completed_lacam(HNode* H_init,std::vector<double>& agent_ratio, 
    std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, HNode*& H_goal);
  
  void propogate_q_value_k_steps(int agent_id, Vertex* const& v_curr, Vertex* const& v_next, double cost_to_go, int k_steps);
  double get_q_value(int agent_id, Vertex* const& v_curr, Vertex* const& v_next);
  void update_q_value(int agent_id, Vertex* const& v_curr, Vertex* const& v_next, double cost_to_go);
  int get_action(Vertex* const& v_curr, Vertex* const& v_next, int width);

  void retrieve_solution(HNode* H_init,std::vector<Config>& solution);
  void add_punishment(HNode* input_H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED); 
  void pick_restart_nodes(std::stack<HNode*>& OPEN);



  void rewrite_backpropagate(HNode* H_from, HNode* H_to, HNode* H_goal,
    std::stack<HNode*>& OPEN, HNode* H_init,std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED);
    
  void rewrite_no_push(HNode* H_from, HNode* H_to, HNode* H_goal);
  void update_ordering(HNode* h_from, HNode* h_to, size_t N, DistTable& D);

  uint run_constrainted_completed_lacam(HNode* H_init_node,std::vector<double>& agent_ratio,
    std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, 
    HNode*& H_goal, std::vector<Config>& solution, std::unordered_set<int>& fixed_agents);
  
  void record_each_agent_frequency(HNode* input_H_goal);
  void increase_visited_node(std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED);
  void increase_solution_weight(HNode* input_H_goal);
  void expand_lowlevel_tree(HNode* H, LNode* L);
  void rewrite(HNode* H_from, HNode* T, HNode* H_goal,
              std::stack<HNode*>& OPEN);
  uint get_edge_cost(const Config& C1, const Config& C2);
  void get_edge_cost_per_agent(std::vector<double>& agent_cost, const Config& C1, const Config& C2);
  void compute_agent_increase_ratio(std::vector<double>& agent_ratio, 
  const Config& C1, const Config& C2);
  uint get_edge_cost(HNode* H_from, HNode* H_to);
  uint get_h_value(const Config& C);
  uint get_sum_of_distance(const Config& C);
  void record_highlevel_node(HNode* H_node, bool expanded, bool contain_target);
  bool get_new_config(HNode* H, LNode* L);
  bool funcPIBT(Agent* ai);
  void increase_each_agent_cost(std::vector<double>& agent_cost, const Config& C1, const Config& C2);
  void increase_all_nodes_visited(std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED);
  void increase_solution_congestion_cost(HNode* input_H_goal);
  void increase_each_agent_cost(HNode* input_H_goal);
  void increase_weight_map(HNode* input_H_goal, bool is_goal);
  void learning_Q_value(HNode* input_H_goal,double is_goal);
  void set_individual_congestion_map(HNode* H_init);

  


  // swap operation
  Agent* swap_possible_and_required(Agent* ai);
  bool is_swap_required(const uint pusher, const uint puller,
                        Vertex* v_pusher_origin, Vertex* v_puller_origin);
  bool is_swap_possible(Vertex* v_pusher_origin, Vertex* v_puller_origin);

  void saveTree(const std::string &fileName) const;

  // utilities
  template <typename... Body>
  void solver_info(const int level, Body&&... body)
  {
    if (verbose < level) return;
    std::cout << "elapsed:" << std::setw(6) << elapsed_ms(deadline) << "ms"
              << "  loop_cnt:" << std::setw(8) << loop_cnt
              << "  node_cnt:" << std::setw(8) << HNode::HNODE_CNT << "\t";
    info(level, verbose, (body)...);
  }
 };
 