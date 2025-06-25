#ifndef LAGRANGIAN_HPP
#define LAGRANGIAN_HPP

#include "../abstractsolver.hpp"
#include "../lconstraints/lconstraintssolver.hpp"
#include "gurobi_c++.h"
#include <iterator>
#include <iostream>

class Lagrangian : public Solver {
    private: 
        int MAX_NB_IT;
        int MAX_NB_IT_WITHOUT_IMPROV;
        double INIT_SCALE;
        double MIN_STEP;
        bool MOD_SUB;

        double start_time;
        
        int nb_it;
        int nb_it_without_improv;
        int nb_total_it_without_improv;
        int nb_infeasible;
        double step;
        double norm;
        double scale;
        bool sub_feasible;
        bool sub_respect_slackness; 

        double curr_lb;       /* Obj. value current subproblem. */
        double curr_ub;       /* Real obj. value current subproblem. */
        int curr_nb_sol;      /* Nb. solutions found solving current subproblem.*/
        Status status_sub;    /* Status current subproblem. */

        std::vector<double> stability_mult;
        std::vector<double> stability_grad;
        std::vector<double> stability_dire;

        double * init_coeff_y = NULL;    /* Initial cost variables y. */
        double * curr_coeff_y = NULL;    /* Current cost variables y (current subproblem). */

        double * curr_t_ = NULL;  /* Current solution subproblem of a node. */
        double * zero_t_ = NULL;  /* Zero extra capacity. */
        double * last_t_ = NULL;  /* Last tested extra capacity at heuristic. */

        double *tooth_ = NULL;

        double time_lp;

        int count_stab;                                                 /* Count number of stability constraints. */             
        std::vector<std::tuple<int,int,int>> inv_index_stab;

        std::vector<std::vector<std::list<int>>> matching;              /* matching[s][c]: matching for school c in scenario s (list of feasible pairs) */
        std::vector<std::vector<std::list<int>>> feas_pairs;            /* feas_pairs[s][st]: list of admissible pairs of student st */
        std::vector<std::unordered_map<int,double>> rank_feas_pairs;    /* rank_feas_pairs[s][i]: rank of feasible pair i */
    
    protected:
        void warm_start();
        void init();
        void upd_mult();
        void upd_grad();
        void upd_norm();
        void upd_prod();
        void upd_dire();

        void upd_RL_lb();
        void upd_RL_ub();
        void upd_RL_gap();
        void upd_scale();
        void upd_step();

        void tooth(int);
        void init_sub_coeff();
        void upd_sub_coeff();
        void upd_curr_bounds();
        void solve_sub();
        
        void init_heuristic();
        void init_true_ranks_heuristic();
        void init_feas_pairs_heuristic();
        void solve_heuristic(const double *);

        void params();
        void variables_2st();
        void index_2st(int);
        
        void upd_solution();
        void write_sol();
        void print_specific(std::string &);

        double get_timelimit() { return std::max(0.0, input.getTimeLimit() - (time(NULL) - start_time)); }

    public:
        Lagrangian(const Input &, const Instance &, int, std::string, GRBEnv *, bool);
        ~Lagrangian();

        void solve();
        void callback(){}  // No callback for lagrangian.
};

#endif
