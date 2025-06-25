#ifndef LOCAL_SEARCH_SOLVER_HPP
#define LOCAL_SEARCH_SOLVER_HPP

#include <ctime>
#include <random>
#include "../abstractsolver.hpp"

class LocalSearchSolver : public Solver {
    protected: 
        bool used_in_ws;
        int nb_max_it = 1000;
        int it;

        double start_time;
        
        double ub_best_it;
        double * tested_t_ = NULL;
        double * best_it_t_ = NULL;

        void warm_start(){}
        void upd_iteration();
        bool time_limit_reached() { return ((used_in_ws == false && (time(NULL) - start_time >= input.getTimeLimit())) ||
                                           (used_in_ws == true && (time(NULL) - start_time >= 600.0))); }
    
    public:
        LocalSearchSolver(const Input &, const Instance &, int, std::string, GRBEnv *, bool, bool);
        ~LocalSearchSolver(){}
        void callback(){} 
        void solve();
};

#endif