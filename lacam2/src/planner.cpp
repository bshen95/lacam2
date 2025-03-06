#include "../include/planner.hpp"

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

  // set order
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
      loop_cnt(0),
      C_next(N),
      tie_breakers(V_size, 0),
      A(N, nullptr),
      occupied_now(V_size, nullptr),
      occupied_next(V_size, nullptr),
      tree_file(save_tree_file)
{
}

Planner::~Planner() {}

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



void Planner::MCT_backpropagate(MCTNode* mct_node, double reward){
  mct_node->visits +=1 ;
  mct_node->reward += reward;
  if(mct_node->parent != nullptr )
  {
    MCT_backpropagate(mct_node->parent, reward );
  }
}


MCTNode* Planner::MCT_selection(std::vector<MCTNode*>& node_pool, HNode* goal_node){
  MCTNode* selected_node = nullptr; 
  double max_utc_value = -std::numeric_limits<double>::infinity();
  // if(goal_node == nullptr){
  //   return  node_pool[node_pool.size()-1];
  // }
  for(auto node : node_pool){
    if(!node->completed_node){
      // std::cout<<"Node ID: "<< node->node_id <<" ";
      double UCT_value = node->compute_uct_value(goal_node);
      if(max_utc_value <= UCT_value){
        selected_node = node; 
        max_utc_value = UCT_value;
      }
    }
  }
  // std::cout<<" "<< std::endl; 
  return selected_node;
}

void Planner::print_utc_value(std::vector<MCTNode*>& node_pool, HNode* goal_node){
  MCTNode* selected_node = nullptr; 
  double max_utc_value = -std::numeric_limits<double>::infinity();
  // if(goal_node == nullptr){
  //   return  node_pool[node_pool.size()-1];
  // }
    std::cout<< "   node_selection:" <<std::endl;
  
  for(auto node : node_pool){
    if(!node->completed_node){
      // std::cout<<"Node ID: "<< node->node_id <<" ";
      double UCT_value = node->compute_uct_value(goal_node);
      std::cout<< "    - "<< "node_id: "<< node->node_id<<std::endl;
      node->print_uct_value(goal_node);
      if(max_utc_value <= UCT_value){
        selected_node = node; 
        max_utc_value = UCT_value;
      }
    }
  }
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

  while (!is_expired(deadline)) {
    if(H_goal != nullptr){
      roll_out_cnt = H_goal->current_make_span;
      // std::cout<<  H_goal->current_make_span<< std::endl;
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
      print_utc_value(OPEN,H_goal);
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
        MCT_backpropagate(MCT_NODE,0);
        k --; 
        continue;
      } 

      // create new configuration
      for (auto a : A) C_new[a->id] = a->v_next;

      // for (size_t i = 0; i < N; ++i) {
      //   auto v_i_from = H->C[i];
      //   auto v_i_to = C_new[i];
      //   // check connectivity
      //   if (v_i_from != v_i_to &&
      //       std::find(v_i_to->neighbor.begin(), v_i_to->neighbor.end(),
      //                 v_i_from) == v_i_to->neighbor.end()) {
      //     std::cout<< "wrong generated !!!" <<std::endl;
      //   }
      // }

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
          double reward = downward_rollout_policy(H_new, H_goal, EXPLORED, 10, 1.2*(roll_out_cnt - H_new->current_make_span));
          sampleing_times += roll_out_cnt - H_new->current_make_span;
          MCT_backpropagate(MCT_new,reward);
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
          MCT_backpropagate(MCT_NODE,0);
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
          double reward = downward_rollout_policy(H_new, H_goal, EXPLORED, 10, 1.2*(roll_out_cnt - H_new->current_make_span));
          // std::cout<< roll_out_cnt << std::endl;
          sampleing_times += roll_out_cnt - H_new->current_make_span;
          MCT_backpropagate(MCT_new,reward);
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
          MCT_backpropagate(MCT_NODE,0);
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


double Planner::downward_rollout_policy(HNode* H, HNode*& H_goal, std::unordered_map<Config, HNode*, ConfigHasher>& EXPLORED,  int num_sample, int max_depth){
  // std::cout<< "Miximual depth: "<< max_depth << std::endl; 
  // max_depth =150;
  int original_sample = num_sample;
  int original_max_depth = max_depth;
  double total_sum_of_f_value = 0;
  uint num_of_goal_hit = 0;
  // std::cout<<"Start sampling: "<< std::endl; 
  double success_times = 0;
  std::vector<std::vector<std::unordered_map<uint, uint>>> sample_array; 
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
    std::vector<std::unordered_map<uint, uint>> agent_to_vertex;
    // std::cout<< "Start sampling: "<< std::endl;
    while (max_depth > 0){
      agent_to_vertex.push_back(computeCellFrequency(current_sample_node->C));
      // check goal condition
      if (is_same_config(current_sample_node->C, ins->goals)) {
        if( H_goal == nullptr || current_sample_node->f < H_goal->f){
          solver_info(1, "rollout-found solution, cost: ", current_sample_node->g);
          H_goal = current_sample_node;
        }
        num_of_goal_hit++;
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

      // if(is_same_config(current_sample_node->C, C_new)){
      //   std::shuffle(order.begin(), order.end(), *MT);
      //   num_of_depth ++; 
      //   max_depth --;
      //   continue;
      // }
      // check explored list
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
        current_sample_node = H_new;
        simulated_g_value = simulated_g_value + get_edge_cost(current_sample_node->C, C_new);
      }
      max_depth --; 
      num_of_depth ++; 
    }
    sample_array.push_back(agent_to_vertex);
    // std::cout<< "Finish sampling: " << num_of_depth<< std::endl; 
    // std::cout<< " "<< std::endl; 
    // std::cout<< " "<< std::endl; 
    // std::cout<< " "<< std::endl; 
    // std::cout<< " "<< std::endl; 

    if(sample_success){
      success_times ++;
      total_sum_of_f_value += simulated_g_value + current_sample_node->h;
    }
    num_sample --;
  }

  for ( int i = 0 ; i < 9 ; i ++){
    for ( int j = i+1 ; j < 9 ; j ++){
    for( int t = 1; t < 2 ; t ++){
            const auto& map1 = sample_array[i][t];
            const auto& map2 = sample_array[j][t];
            int num_of_diff = 0;
            for (const auto& entry : map1) {
                uint cell_index = entry.first;
                uint freq1 = entry.second;
                uint freq2 = map2.count(cell_index) ? map2.at(cell_index) : 0;
                if (freq1 != freq2) {
                    num_of_diff ++;
                } 
                // if (freq1 != freq2) {
                //     std::cout << "Cell Index: " << cell_index << ", Frequency in map " << i << ": " << freq1 << ", Frequency in map " << i+1 << ": " << freq2 << "\n";
                // }
            }
            std::cout << "Number of different cells between map " << i << " and map " << j << ": " << num_of_diff << "\n";
            std::cout << "time step"<< t  << "\n";
          }
    }
  }

  // std::cout<< "Finish all sampling........... " <<std::endl; 
  // std::cout<< " "<< std::endl; 
  // std::cout<< " "<< std::endl; 
  // std::cout<< " "<< std::endl; 
  // std::cout<< " "<< std::endl; 
  // std::cout<< " "<< std::endl; 
  // std::cout<< " "<< std::endl; 
  // std::cout<< " "<< std::endl; 
  // std::cout<< " "<< std::endl; 
  // double f_value_ratio = 1 - total_sum_of_f_value / original_sample / max_f_value ;
  // double f_value_ratio = 1 - total_sum_of_f_value / original_sample /  max_improvement ;
  // std::cout<<(double)(min_f_value) <<std::endl;
  double f_value_ratio = 1 - total_sum_of_f_value / success_times / (N * get_makespan_lower_bound(ins->starts));
  if(H_goal != nullptr){
    if(success_times == 0){
      f_value_ratio = 0 ;
    }else{
      f_value_ratio = 1 - total_sum_of_f_value / success_times / H_goal->f;
    }
  }
  double success_ratio = success_times / original_sample;

  // double ratio = (H_goal->f - H->g) / (total_sum_of_f_value / original_sample - H->g);
  // std::cout<< ratio << std::endl;
  // std::cout<< (H_goal->f - H->g) << std::endl;
  // std::cout<< total_sum_of_f_value / original_sample << std::endl;

  // std::cout <<" success_ratio: " << success_ratio << " f_value_ratio: "<< f_value_ratio << std::endl;
  // double reward = std::max((closeness_ratio + congestion_ratio + f_value_ratio)/3, goal_hit_ratio);
  // double reward = ( (closeness_ratio + congestion_ratio + f_value_ratio)/3 + goal_hit_ratio ) / 2 ;
  double reward =  success_ratio + f_value_ratio;
  //  + g_increase_ratio;
  // std::cout<<"reward value: " << 2 * success_ratio + 2 * goal_ratio + g_ratio <<"  "<< std::endl;

  // std::cout<<"Finishing sampling.       "<< std::endl; 
  // std::cout<<"Reward: "<< reward<< std::endl; 
  // std::cout<<"                          "<< std::endl; 
  // std::cout<<"                          "<< std::endl; 

  return reward;
}


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
        if (n_to == H_goal)
          solver_info(1, "cost update: ", n_to->g, " -> ", g_val);
        n_to->g = g_val;
        n_to->f = n_to->g + n_to->h;
        n_to->parent = n_from;

        for (size_t i = 0; i < N; ++i) {
          auto v_i_from = n_from->C[i];
          auto v_i_to = n_to->C[i];
          // check connectivity
          if (v_i_from != v_i_to &&
              std::find(v_i_to->neighbor.begin(), v_i_to->neighbor.end(),
                        v_i_from) == v_i_to->neighbor.end()) {
            std::cout<< "here !!!" <<std::endl;
          }
        }
        Q.push(n_to);
        if (H_goal != nullptr && n_to->f < H_goal->f) {
          OPEN.push(n_to);
          record_highlevel_node(n_to,false,false);
        }
      }
    }
  }
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
  for (uint k = 0; k < L->depth; ++k) {
    const auto i = L->who[k];        // agent
    const auto l = L->where[k]->id;  // loc
    // check vertex collision
    if (occupied_next[l] != nullptr) return false;
    // check swap collision
    auto l_pre = H->C[i]->id;
    if (occupied_next[l_pre] != nullptr && occupied_now[l] != nullptr &&
        occupied_next[l_pre]->id == occupied_now[l]->id)
      return false;

    // set occupied_next
    A[i]->v_next = L->where[k];
    occupied_next[l] = A[i];
  }

  // perform PIBT
  for (auto k : H->order) {
    auto a = A[k];
    if (a->v_next == nullptr && !funcPIBT(a)) return false;  // planning failure
  }
  return true;
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
  std::sort(C_next[i].begin(), C_next[i].begin() + K + 1,
            [&](Vertex* const v, Vertex* const u) {
              return D.get(i, v) + tie_breakers[v->id] <
                     D.get(i, u) + tie_breakers[u->id];
            });

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