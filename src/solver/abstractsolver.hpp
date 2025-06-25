#ifndef ABSTRACT_SOLVER_HPP
#define ABSTRACT_SOLVER_HPP

#include "../input/input.hpp"
#include "../instance/instance.hpp"
#include "../transform/transformfactory.hpp"
#include "status/status.hpp"
#include "gurobi_c++.h"

#include <limits>
#include <cmath>
#include <set>
#include <optional>

class AbstractTransform;
class ASGheuristicSolver;
class LocalSearchSolver;

#define UNUSED(x) (void)(x)

class AbstractSolver {
    protected:
        const Input & input;
        const Instance & instance;

        int pr;
        std::string name;
        int verbose;

        double lb;
        double ub;
        double gap;
        double time_;
        double nb_bnb;
        int nb_sol;
        Status status;

        double ub_ws;

        virtual void params() = 0;
        virtual void write_inf() = 0;
        virtual void write_sol() = 0;
        virtual void print();
        virtual void verify_stat();

    public:
        AbstractSolver(const Input & input,const Instance & instance,int pr,std::string name) : input(input), instance(instance), pr(pr), name(name), verbose(input.getVerbose()) {
            lb = -std::numeric_limits<double>::infinity();
            ub = std::numeric_limits<double>::infinity();
            gap = 1.0;
            time_ = 0.0;
            nb_bnb = 0.0;
            nb_sol = 0;
            status = Status::Not_Solved;
        } 

        virtual ~AbstractSolver() {}

        double getLB() const { return lb; }
        double getUB() const { return ub; }
        double getGAP() const { return gap; }
        double getTime() const { return time_; }
        double getNbBnB() const { return nb_bnb; }
        double getNbSol() const { return nb_sol; }
        Status getStatus() const { return status; }  
        void setVerbose(int v) { verbose = v; }
};

class Solver : public AbstractSolver, public GRBCallback {
    protected:
        AbstractTransform *transf = NULL;
        
        ASGheuristicSolver *asg = NULL;
        Test1stSolver *test_1st = NULL;
        LocalSearchSolver *ls = NULL;

        bool warmstart;
        bool linear_relax;
        bool combs;

        double ub_callback;
        double time_best_ub;
        
        bool include_frac_cuts;
        bool y_continuous;

        GRBModel *model = NULL;
        GRBVar *t = NULL;
        GRBVar *y = NULL;
        GRBVar *x = NULL;
        GRBVar *xa = NULL;
        GRBVar * vtooth = NULL;

        int count_tt;
        int count_y;  
        int count_x;
        
        std::vector<std::vector<int>> index_y;
        std::vector<std::vector<int>> index_x;

        double *t_ = NULL;
        double *y_ = NULL;
        double *x_ = NULL;
        double *xa_ = NULL;
        double *tooth_minus_ = NULL;

        int count_callback;

        virtual void params();
        virtual void write_inf();
        virtual void write_sol();

        void create_1st();
        virtual void index_1st();
        virtual void variables_1st();
        virtual void constraints_1st();
        
        void create_2st();
        virtual void index_2st(int);
        virtual void variables_2st();
        virtual void constraints_2st(int);

        virtual void get_variables();
        virtual void warm_start();
        virtual void upd_solution();

        virtual void warm_start_asg();
        virtual void warm_start_capazero();
        virtual void warm_start_ls();

        void tooth_minus(int);
        void comb_baiou(int &, int &, bool &, GRBLinExpr &, double *);
        void generate_comb_expr(int &, int &, GRBLinExpr &, std::set<int> &);

    public:
        Solver(const Input&,const Instance&,int,std::string,GRBEnv*,bool,bool); 
        virtual ~Solver();

        AbstractTransform * getTransform() { return transf; }
        Test1stSolver * getTest1st() { return test_1st; }

        const double * getT_() const { return t_; }
        const double * getY_() const { return y_; }

        const GRBVar * getT() const { return t; }
        const GRBVar * getY() const { return y; }

        int getCountY() const { return count_y; }
        int getCountT() const { return count_tt; }
        
        GRBVar & getT(int c) { return t[c]; }
        virtual GRBVar & getX(int c, int j) { return x[index_x[c][j]]; }
        GRBVar & getXA(int c, int j) { return xa[index_x[c][j]]; }
        GRBVar & getVTooth(int s, int i) { return vtooth[index_y[s][i]]; }
        virtual GRBModel * getModel() { return model; }
        virtual GRBVar & getY(int s, int i) { return y[index_y[s][i]]; }
        virtual double getY_(int s, int i) const { return y_[index_y[s][i]]; }

        bool getLinearRelax() const { return linear_relax; }

        virtual void callback();  
        virtual void solve();
};

#endif





