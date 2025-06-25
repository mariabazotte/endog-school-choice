#ifndef ABSTRACT_TRANSFORM_HPP
#define ABSTRACT_TRANSFORM_HPP

#include "../input/input.hpp"
#include "../instance/instance.hpp"
#include "gurobi_c++.h"
#include <limits>
#include <list>

#define UNUSED(x) (void)(x)

class Solver;
class Test1stSolver;

class AbstractTransform {
    protected:
        const Input & input;
        const Instance & instance;

        int pr;
        Solver *solver = NULL;

    public:
        AbstractTransform(const Input & input, const Instance & instance, int pr, Solver *solver): input(input),instance(instance),pr(pr),solver(solver){}
        virtual ~AbstractTransform(){}
                   
        virtual void transf_2st() = 0;                                     // Modify the model (variables + constraints) according to the selected transformation
        virtual void get_variables() = 0;                                  // Get the value of the transformation variables after optimizing the model (optimal/best solution)

        virtual void set_inisol_warm_start(AbstractTransform*) = 0;        // Set initial solution for transformation when doing warm start from heuristic based on model gurobi
        virtual void set_inisol_warm_start(Test1stSolver*) = 0;            // Set initial solution for transformation when doing warm start from heuristic based on da 
        virtual void set_constr_warm_start(Test1stSolver*) = 0;            // Include constraints according to transformation when doing warm start

        virtual void term_xi(int, int, double, double, GRBLinExpr &) = 0;  // Include a term with xi variable in a gurobi expression

        virtual void compute_feas_pairs(int,std::vector<std::list<int>> &,
                                        std::unordered_map<int,double> &) = 0; // Compute feasible pairs following the solution obtained after optimizing the model
};

#endif