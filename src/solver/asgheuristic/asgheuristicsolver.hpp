#ifndef ASG_HEURISTIC_SOLVER_HPP
#define ASG_HEURISTIC_SOLVER_HPP

#include "gurobi_c++.h"

#include <unordered_map>
#include <map>
#include <tuple>
#include <list>
#include <ctime>
#include "../abstractsolver.hpp"

class ASGheuristicSolver : public Solver {
    protected: 
        std::vector<std::vector<std::list<int>>> matching;              /* matching[s][c]: matching for school c in scenario s (list of feasible pairs) */
        std::vector<std::vector<std::list<int>>> feas_pairs;            /* feas_pairs[s][st]: list of admissible pairs of student st */
        std::vector<std::unordered_map<int,double>> rank_feas_pairs;    /* rank_feas_pairs[s][i]: rank of feasible pair i */
        
        void upd_solution();
        void warm_start();
    
    public:
        ASGheuristicSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, bool warmstart, bool linear_relax): 
                            Solver(input,instance,pr,name,env,warmstart,linear_relax) {
            params();
            create_1st();
            create_2st();  

            matching.resize(instance.getNbScenarios(pr));
            feas_pairs.resize(instance.getNbScenarios(pr));
            rank_feas_pairs.resize(instance.getNbScenarios(pr));
            for(int s = 0; s < instance.getNbScenarios(pr); ++s){
                feas_pairs[s].resize(instance.getNbStudents());
                matching[s].resize(instance.getNbSchools()+1);
            }
        }
        ~ASGheuristicSolver(){}
        void callback(){} 
};

#endif