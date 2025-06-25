#ifndef LCONSTRAINTS_SOLVER_HPP
#define LCONSTRAINTS_SOLVER_HPP

#include "../abstractsolver.hpp"
#include "gurobi_c++.h"

#include <set>
#include <list>
#include <map>
#include <tuple>
#include <ctime>

class LConstraintsSolver : public Solver  {
    protected: 
        void constraints_2st(int);

    public:
        LConstraintsSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, bool warmstart, bool linear_relax): 
                            Solver(input,instance,pr,name,env,warmstart,linear_relax) {
            params();
            create_1st();         
            create_2st();  
        }
        ~LConstraintsSolver(){}
};

class ModLConstraintsSolver : public Solver  {
    protected: 
        std::vector<GRBConstr> constrs;
        double *dual_ = NULL;

        void constraints_2st(int);
        void upd_solution();

    public:
        ModLConstraintsSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, bool warmstart, bool linear_relax): 
                            Solver(input,instance,pr,name,env,warmstart,linear_relax) {
            params();
            create_1st();         
            create_2st();  
        }
        ~ModLConstraintsSolver(){
            if(dual_) delete[] dual_;
        }

        const double * getDual() const { return dual_; }
};

#endif