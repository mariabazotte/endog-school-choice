#ifndef CUTOFF_SCORE_SOLVER_HPP
#define CUTOFF_SCORE_SOLVER_HPP

#include "../abstractsolver.hpp"
#include "gurobi_c++.h"

#include <map>
#include <tuple>
#include <ctime>

class CutoffScoreSolver : public Solver  {
    protected:  
        GRBVar *z = NULL;
        GRBVar *f = NULL;

        int count_zf;

        std::map<std::tuple<int,int>,int> index_zf;

        GRBVar & getZ(int s, int c) { return z[index_zf[{s,c}]]; }
        GRBVar & getF(int s, int c) { return f[index_zf[{s,c}]]; }

        virtual void index_2st(int);
        virtual void variables_2st();
        virtual void constraints_2st(int);

    public:
        CutoffScoreSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, bool warmstart, bool linear_relax): 
                                     Solver(input,instance,pr,name,env,warmstart,linear_relax) { 
            count_zf = 0; 
            params();
            create_1st();
            create_2st();
        }
        ~CutoffScoreSolver(){
            delete[] z;
            delete[] f;
        } 
};

#endif