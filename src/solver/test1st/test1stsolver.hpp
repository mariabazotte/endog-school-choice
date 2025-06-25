#ifndef TEST_1ST_SOLVER_HPP
#define TEST_1ST_SOLVER_HPP

#include "../abstractsolver.hpp"

#define UNUSED(x) (void)(x)

class Test1stSolver : public AbstractSolver {
    protected:
        double * obj_ = NULL;                            // obj_[s]: objective value scenario s
        int ** matching_students = NULL;                 // matching_students[s][st]: feasible pair i student st is matched in scenario s
        int ** nb_assign = NULL;                         // nb_assign[s][c]: nb of students school c is assigned to in scenario s
        bool * unsatured = NULL;                         // unsatured[c]: true if school c is unsatures (considering all scenarios)
        std::vector<std::vector<double>> nb_per_rank;    // Nb of assigned students per student rank and scenario
        std::vector<double> avg_per_rank;                // Average number of assigned students per student rank
        std::vector<std::vector<double>> rank_st;        // Rank student is assigned to in scenario s. rank_st[s][st]

        double *y_ = NULL;
        double *xi_ = NULL;

        int count_y;
        int count_xi;

        std::vector<std::vector<int>> index_y;
        std::vector<std::vector<int>> index_xi;

        int curr_s;
        
        virtual void index_2st(int);
        virtual void get_variables(int) = 0;
        virtual void transf_nonstrategic() = 0;
        virtual void transf_strategic(const double *,Input::Transform) = 0;
        
        virtual void solve(const double *) = 0;
        virtual void transf(const double *,Input::Transform) = 0;

        void reset();
        void write_sol(); 
        virtual void write_inf();
        void verify_stat();
    
    public:
        Test1stSolver(const Input &, const Instance &, int, std::string);
        ~Test1stSolver();

        double getObjScenario(int s) const { return obj_[s]; }
        int getMatch(int s, int st) const { return matching_students[s][st]; } // return feasible pair i
        bool getUnsatured(int c) const { return unsatured[c]; }

        const double * getY_() { return y_; }
        const double * getXi_() { return xi_; }

        void solve(const double *,Input::Transform);

        double getAvgRank(int r) const { return avg_per_rank[r]; }
        double getRankStudent(int s, int st) const { return rank_st[s][st]; }
};

class Test1stSolverModel : public Test1stSolver {
    protected:
        GRBModel **model = NULL;
        GRBVar **y = NULL;
        GRBVar **z = NULL;
        GRBVar **f = NULL;

        GRBConstr **constr_capa = NULL;
        GRBConstr **constr_waste = NULL;
        GRBConstr **constr_cutoff = NULL;

        double ** xis_ =  NULL;
        double ** rho_ = NULL;

        int count_zf;
        std::map<std::tuple<int,int>,int> index_zf;

        void params();
        void create();
        void index_2st(int);
        void variables_2st(int);
        void constraints_2st(int);

        void write_inf();

        void get_variables(int);
        void transf_nonstrategic();
        void transf_strategic(const double *,Input::Transform);

        void solve(const double *);
        void transf(const double *, Input::Transform);
    
    public: 
        Test1stSolverModel(const Input &, const Instance &, int, std::string, GRBEnv *);
        ~Test1stSolverModel();

        void solve();
};

class Test1stSolverDA : public Test1stSolver {
    protected:
        std::vector<std::vector<std::list<int>>> matching_schools;      /* matching_schools[s][c]: matching for school c in scenario s (list of feasible pairs) */
        std::vector<std::vector<std::list<int>>> feas_pairs;            /* feas_pairs[s][st]: list of admissible pairs of student st */
        std::vector<std::unordered_map<int,double>> rank_feas_pairs;    /* rank_feas_pairs[s][i]: rank of feasible pair i */

        void params(){}
        void get_variables(int);
        void transf_nonstrategic();
        void transf_strategic(const double *,Input::Transform);

        void solve(const double *);
        void transf(const double *, Input::Transform);
        
    public: 
        Test1stSolverDA(const Input &, const Instance &, int, std::string);
        ~Test1stSolverDA() {}

        void solve();

        double getRankFeasPair(int s, int i) const { return rank_feas_pairs[s].at(i); }                 // Rank feasible pair i at scenario s
        int getFeasPairByStudentAndPos(int s, int st, int k) const {                                    // School reported by student s in scenario st at position k
            auto it = feas_pairs[s][st].begin();
            std::advance(it, k);
            return *it; 
        }   
};



#endif