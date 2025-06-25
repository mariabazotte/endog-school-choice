#ifndef LSHAPED_HPP
#define LSHAPED_HPP

#include "../abstractsolver.hpp"
#include "gurobi_c++.h"
#include <iterator>
#include <iostream>

class IntegerLShaped : public Solver {
    private:
        GRBVar theta;
        double theta_;

        double start_time;
        int nb_it;
        int nb_int_opt_cuts;
        int nb_cont_opt_cuts;
        double intfeastol;
        double precision;      // Precision for MIP gap
        bool lshaped_callback; // 0-> Do not use callback, 1-> Use callback
        bool alternate_method; // 0-> Do not use alternate cuts, 1 -> Use alternate cuts
        double curr_2st_value; // Second stage problem value for current first stage solution
        double lower_bound;
        Status substatus;

        void params();
        void variables_2st();
        void separate(GRBLinExpr &);
        void index_1st();
        void variables_1st();
        void constraints_1st();
        
        void compute_lower_bound();
        void varx_from_vart();
        void warm_start_capazero();
        void warm_start_ls();
        void warm_start_asg();
        void get_variables();

        void init();
        void verify_vars();
        double get_timelimit(){ return std::max(0.0, input.getTimeLimit() - (time(NULL) - start_time)); }

        void integer_loop();
        void solve_main();
        void solve_integer_sub();
        void add_integer_cut();
        void separate_integer_cut(GRBLinExpr &);

        /////////////////////////////////
        bool alternating_method;

        GRBModel *subproblem = NULL;
        GRBVar *sub_x = NULL;
        GRBVar *sub_y = NULL;
        GRBVar *sub_z = NULL;
        GRBVar *sub_f = NULL;
        GRBConstr *constr = NULL;
        double * duals = NULL;
        
        int count_zf;
        std::map<std::tuple<int,int>,int> index_zf;
        GRBVar & getX(int c, int j) { return sub_x[index_x[c][j]]; }
        GRBVar & getY(int s, int i) { return sub_y[index_y[s][i]]; }
        GRBVar & getZ(int s, int c) { return sub_z[index_zf[{s,c}]]; }
        GRBVar & getF(int s, int c) { return sub_f[index_zf[{s,c}]]; }
        GRBModel * getModel() { return subproblem; }

        void create_sub_2st();
        void index_sub_2st(int);
        void variables_sub_2st();
        void constraints_sub_2st(int);
        void solve_continuous_sub();
        void separate_continuous_cut(GRBLinExpr &);
        bool add_continuous_cut();

    public:
        IntegerLShaped(const Input&,const Instance&,int,std::string,GRBEnv*,bool); 
        ~IntegerLShaped();

        void callback();
        void solve();
};
        
#endif
