#ifndef SIMUL_ANNEA_SOLVER_HPP
#define SIMUL_ANNEA_SOLVER_HPP

#include <ctime>
#include <random>
#include "../abstractsolver.hpp"

class SimulAnnealingSolver : public Solver {
    protected: 

        int nb_max_it = 300;
        int it;

        double start_time;
        
        double ub_current;
        double * tested_t_ = NULL;
        double * current_t_ = NULL;

        // Simulated Annealing parameters
        double T = 100.0;         // Initial temperature
        double alpha = 0.95;      // Cooling rate
        int iter_per_temp = 25;   // Number of iterations per temperature

        void warm_start(){}
        bool time_limit_reached() { return ((time(NULL) - start_time >= input.getTimeLimit())); }
    
    public:
        SimulAnnealingSolver(const Input &, const Instance &, int, std::string, GRBEnv *, bool);
        ~SimulAnnealingSolver(){}
        void callback(){} 
        void solve();
};

#endif



