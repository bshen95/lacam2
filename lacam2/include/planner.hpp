/*
 * lacam-star
 */

 #pragma once

 #include "dist_table.hpp"
 #include "graph.hpp"
 #include "instance.hpp"
 #include "utils.hpp"
 #include "min_max_stats.hpp"
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
   const uint depth;
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
 
 
   // for low-level search
   std::vector<float> priorities;
   std::vector<uint> order;
   std::queue<LNode*> search_tree;
 
   uint node_id; 
 
   HNode(const Config& _C, DistTable& D, HNode* _parent, const uint _g,
         const uint _h);
   ~HNode();
   
   void setNodeID(uint id) {
     node_id = id;
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
     uint node_id; 
     
     MCTNode(MCTNode* _parent, HNode* _HNode, uint _node_id);
     ~MCTNode();
   
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
   std::pair<double,int> downward_rollout_policy(HNode* H, HNode*& H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED, int num_sample, int max_depth);
   bool get_next_configuration_rollout(const Config& current_config, const std::vector<uint>& order, LNode* L);
   MCTNode* MCT_selection(std::vector<MCTNode*>& node_pool,HNode* goal_node);
   void MCT_backpropagate(MCTNode* mct_node, int sample_times,  double reward, int failuare_times);
   int get_rollout_increase();
   int get_makespan_lower_bound(const Config& C);
   double compute_uct_value(MCTNode* mcts_node);
   void print_utc_value(std::vector<MCTNode*>& node_pool);
 
 
 
 
 
 
   void expand_lowlevel_tree(HNode* H, LNode* L);
   void rewrite(HNode* H_from, HNode* T, HNode* H_goal,
                std::stack<HNode*>& OPEN);
   uint get_edge_cost(const Config& C1, const Config& C2);
   uint get_edge_cost(HNode* H_from, HNode* H_to);
   uint get_h_value(const Config& C);
   uint get_sum_of_distance(const Config& C);
   void record_highlevel_node(HNode* H_node, bool expanded, bool contain_target);
   bool get_new_config(HNode* H, LNode* L);
   bool funcPIBT(Agent* ai);
 
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
 