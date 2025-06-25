#ifndef STATUS_HPP
#define STATUS_HPP

#include "gurobi_c++.h"

enum Status{
    Not_Solved = 0,
    Unknown = 1,
    Optimal = 2,
    Optimal_Integer = 3,   // Only used for B&B lagrangian and cg
    Time_Limit = 4,
    Node_Limit = 5,
    Iteration_Limit = 6,
    Unbounded = 7,
    Infeasible = 8,
    Infeasible_or_Unbounded = 9
};

Status statusFromGurobi(int gurobi_status);
std::ostream& operator<<(std::ostream& lhs, const Status & status);

#endif