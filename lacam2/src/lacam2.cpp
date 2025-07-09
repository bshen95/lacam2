#include "../include/lacam2.hpp"

Solution solve(const Instance& ins, std::string& additional_info,
               const int verbose, const Deadline* deadline, std::mt19937* MT,
               const Objective objective, const float restart_rate, const std::string save_tree_file)
{
  auto planner = Planner(&ins, deadline, MT, verbose, objective, restart_rate, save_tree_file);
  // auto plan = planner.solve(additional_info);
  // return planner.MCT_vanilla_lacam(additional_info);
  // return planner.MCT_multiple_lacam(additional_info);
  return planner.backpropagate_solve(additional_info);
  // return planner.MCT_solve(additional_info);
  // return planner.MCT_lacam(additional_info);
  // return  planner.solve(additional_info);
}



// void print_plan_state(const Planner& planner, const Solution& solution){
//   std::cout<< get_makespan(solution)<<  ","<< get_sum_of_costs(solution) << "," << get_sum_of_loss(solution)<<std::endl;
// }