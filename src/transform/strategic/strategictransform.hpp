#ifndef STRAT_TRANSFORM_HPP
#define STRAT_TRANSFORM_HPP

#include "../abstracttransform.hpp"

/* For compact and enumeration versions. 
If some school has probability equals to 1.0, then optimizing the expected utility has multiple optimal solutions.
-> safety schools = false: 
    Student st can list less than K (max rank schools) schools even if the nb of schools with positivy utility for st is >= K.
    If a school c with pref^st_c has probability p^st_c = 1.0, then student st does not include safety schools, i.e., he does not list any school c': pref^st_c' >= pref^st_c.
-> safety schools = true: Student st always list K schools if the nb of schools with positivy utility for st is >= K.
    If a school c with pref^st_c <= K has probability p^st_c = 1.0, then student st includes safety schools with best utilities in his list until we have K schools (pref^st_c <= pref^st_c' <= K). */

/* We also have reported or true preferences as objective. */

class StrategicTransform : public AbstractTransform {
    protected: 
        GRBVar *xi = NULL;
        GRBVar *phi = NULL;
        
        int count_xi;
        int count_phi;
        
        std::vector<std::vector<int>> index_xi;
        std::vector<std::vector<int>> index_phi;

        double *xi_ = NULL;
        double *phi_ = NULL;

        GRBVar & getXi(int s, int i) { return xi[index_xi[s][i]]; }
        GRBVar & getPhi(int s, int st) { return phi[index_phi[s][st]]; }

        GRBLinExpr getX(int, int);
        double getX_(int,int);

        void params();
        void objective_2st();
        virtual void index_2st(int,int);
        virtual void variables_2st();
        virtual void constraints_2st(int,int);
        virtual void safety_strategy_2st(int,int);
        void constraints_preprocess_2st(int,int);

    public:
        StrategicTransform(const Input & input, const Instance & instance, int pr, Solver *solver) : 
                                                        AbstractTransform(input,instance,pr,solver) { 
            count_xi = 0;
            count_phi = 0; 
        }
        
        virtual ~StrategicTransform() { 
            if(xi) delete[] xi;
            if(phi) delete[] phi; 

            if(xi_) delete[] xi_;
            if(phi_) delete[] phi_;
        }
        
        void transf_2st();
        void get_variables();
        
        void set_inisol_warm_start(AbstractTransform*);
        void set_inisol_warm_start(Test1stSolver*);
        void set_constr_warm_start(Test1stSolver* test_1st){ UNUSED(test_1st); } // We cannot preprocess when solving capacity zero
        
        void term_xi(int, int, double, double, GRBLinExpr &);
        void compute_feas_pairs(int,std::vector<std::list<int>> &,std::unordered_map<int,double> &);

        // Getters
        const double * getXi_() const { return xi_; }
        const double * getPhi_() const { return phi_; }

        GRBVar * getXi() { return xi; }
        GRBVar * getPhi() { return phi; }

        void print_xi();
};

#endif