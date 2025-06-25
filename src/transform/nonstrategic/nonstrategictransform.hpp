#ifndef NONSTRAT_TRANSFORM_HPP
#define NONSTRAT_TRANSFORM_HPP

#include "../abstracttransform.hpp"

class NonStrategicTransform : public AbstractTransform {       
    public:
        NonStrategicTransform(const Input & input, const Instance & instance, int pr, Solver *solver) : AbstractTransform(input,instance,pr,solver){}
        ~NonStrategicTransform(){}
        
        void transf_2st();
        void get_variables(){}     // No extra variables needed in this case

        void set_inisol_warm_start(AbstractTransform* transf){ UNUSED(transf); } // No extra variables needed in this case 
        void set_inisol_warm_start(Test1stSolver* test_1st){ UNUSED(test_1st); } // No extra variables needed in this case 
        void set_constr_warm_start(Test1stSolver*);                              // We can preprocess many variables/constraints with solution from zero extra capa for nonstrategic version

        void term_xi(int,int,double,double,GRBLinExpr &);               

        void compute_feas_pairs(int,std::vector<std::list<int>> &,std::unordered_map<int,double> &); 
};

#endif