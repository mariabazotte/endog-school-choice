#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include "gurobi_c++.h"

#include "../input/input.hpp"
#include "../instance/instance.hpp"

#include "../solver/status/status.hpp"
#include "../solver/solverfactory.hpp"
#include "../solver/abstractsolver.hpp"
#include "../solver/test1st/test1stfactory.hpp"
#include "../solver/test1st/test1stsolver.hpp"

#include "../transform/nonstrategic/nonstrategictransform.hpp"
#include "../transform/strategic/strategictransform.hpp"

class Algorithm {
    protected:
        const Input & input;
        const Instance & instance;
        
        GRBEnv *env = NULL;
        GRBEnv *env_test = NULL;

        double lb;
        double ub;
        double gap;
        double time_;
        double nb_bnb;
        Status status;

        std::vector<Solver*> solvers;
        Test1stSolver * test_stochastic = NULL;

        double ub_test;
        double var_test;
        double time_test;
        std::vector<double> avg_stu_per_rank;          /*Average nb of students assigned to 1st, 2nd, 3rd ... options */
        double nb_improv;
        double nb_enter;

        double ub_test_strat;                          /* When strategic. */
        double var_test_strat;
        double time_test_strat;
        std::vector<double> avg_stu_per_rank_strat;
        double nb_improv_strat;
        double nb_enter_strat;

        double ub_test_strat_simple;                   /* When simple strategic. */
        double var_test_strat_simple;
        double time_test_strat_simple;
        std::vector<double> avg_stu_per_rank_strat_simple;
        double nb_improv_strat_simple;
        double nb_enter_strat_simple;

        double getLB() const { return lb; }
        double getUB() const { return ub; }
        double getGap() const { return gap; }
        double getTime() const { return time_; }
        double getNbBranchNodes() const { return nb_bnb; }
        Status getStatus() const { return status; } 

        void create();
        void computeVariance(double &,double &);
        void computeTest(int,double &, double &, double &, std::vector<double> &, double &, double &, Input::Transform);
        
    public:
        Algorithm(const Input & input, const Instance & instance) : input(input), instance(instance), env(new GRBEnv()), env_test(new GRBEnv()) {
            lb = -GRB_INFINITY;
            ub = GRB_INFINITY;
            gap = 1.0;
            time_ = 0.0;
            nb_bnb = 0.0;
            status = Status::Not_Solved;
            
            ub_test = 0.0;
            var_test = 0.0;
            time_test = 0.0;
            avg_stu_per_rank.resize(instance.getNbMaxRankSchools()+1,0.0);
            nb_improv = 0.0;
            nb_enter = 0.0;

            ub_test_strat = 0.0;
            var_test_strat = 0.0;
            time_test_strat = 0.0;
            avg_stu_per_rank_strat.resize(instance.getNbMaxRankSchools()+1,0.0);
            nb_improv_strat = 0.0;
            nb_enter_strat = 0.0;

            ub_test_strat_simple = 0.0;
            var_test_strat_simple = 0.0;
            time_test_strat_simple = 0.0;
            avg_stu_per_rank_strat_simple.resize(instance.getNbMaxRankSchools()+1,0.0);
            nb_improv_strat_simple = 0.0;
            nb_enter_strat_simple = 0.0;

            create();
        }

        virtual ~Algorithm(){
            delete test_stochastic;
            for(int i = 0; i < instance.getNbProblems(); ++i) delete solvers[i];
            delete env;
            delete env_test;
        }

        virtual void solve();

        virtual std::string write() const {
            std::string output = std::to_string(status) + std::string(";") + 
                                std::to_string(lb) + std::string(";") + 
                                std::to_string(ub) + std::string(";") + 
                                std::to_string(gap) + std::string(";") + 
                                std::to_string(time_) + std::string(";") +
                                std::to_string(nb_bnb) + std::string(";");

            if(input.getTransform() == Input::Transform::NonStrategic){
                output += std::to_string(ub_test) + std::string(";") + std::to_string(var_test) + std::string(";") + std::to_string(time_test) + std::string(";");
                for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k) 
                    output += std::to_string(avg_stu_per_rank[k]) + std::string(";");
                output += std::to_string(nb_enter) + std::string(";") + std::to_string(nb_improv) + std::string(";");
            }
            if(input.getTransform() == Input::Transform::NonStrategic || input.getTransform() == Input::Transform::Simple){
                output += std::to_string(ub_test_strat_simple) + std::string(";") + std::to_string(var_test_strat_simple) + std::string(";") + std::to_string(time_test_strat_simple) + std::string(";");
                for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k) 
                    output += std::to_string(avg_stu_per_rank_strat_simple[k]) + std::string(";");
                output += std::to_string(nb_enter_strat_simple) + std::string(";") + std::to_string(nb_improv_strat_simple) + std::string(";");
            }
            output += std::to_string(ub_test_strat) + std::string(";") + std::to_string(var_test_strat) + std::string(";") + std::to_string(time_test_strat) + std::string(";");
            for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k) 
                output += std::to_string(avg_stu_per_rank_strat[k]) + std::string(";");
            output += std::to_string(nb_enter_strat) + std::string(";") + std::to_string(nb_improv_strat) + std::string(";");
            return output;
        }

        virtual std::string writeHead() const {
            std::string output = "STATUS;LB;UB;GAP;TIME;NB_BNB;";
            
            if(input.getTransform() == Input::Transform::NonStrategic){
                output += "VAL_STOCH_NON_STRAT;VAR_STOCH_NON_STRAT;TIME_STOCH_NON_STRAT;";
                for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k)
                    output += "AVG_RANK_" + std::to_string(k+1) +"_NON_STRAT;";
                output += "NB_ENTER_NON_STRAT;NB_IMPROV_NON_STRAT;";
            }
            if(input.getTransform() == Input::Transform::NonStrategic || input.getTransform() == Input::Transform::Simple){
                output += "VAL_STOCH_SIMPLE;VAR_STOCH_SIMPLE;TIME_STOCH_SIMPLE;";
                for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k)
                    output += "AVG_RANK_" + std::to_string(k+1) +"_SIMPLE;";
                output += "NB_ENTER_STRAT_SIMPLE;NB_IMPROV_STRAT_SIMPLE;";
            }
            output += "VAL_STOCH_STRAT;VAR_STOCH_STRAT;TIME_STOCH_STRAT;";
            for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k)
                output += "AVG_RANK_" + std::to_string(k+1) +"_STRAT;";  
            output += "NB_ENTER_STRAT;NB_IMPROV_STRAT;";
            return output;
        }
};


#endif