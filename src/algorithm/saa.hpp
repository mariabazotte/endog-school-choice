#ifndef SAA_HPP
#define SAA_HPP

#include "gurobi_c++.h"
#include "../input/input.hpp"
#include "../instance/instance.hpp"
#include "algorithm.hpp"

class SAA : public Algorithm {
    private:
        double var_lb;
        double var_ub;
        double var_gap;
        double var_nb_bnb;

        double stat_lb;
        double stat_ub;
        double stat_gap;
        int nb_notopt;
        int best_problem;

        double time_lb;
        double time_eval;

        double nb_enter_aux;
        double nb_improv_aux;
        
        std::vector<double> lb_vector;
        std::vector<double> ub_vector;
        std::vector<double> gap_vector;
        std::vector<double> nb_bnb_vector;

        std::vector<double> avg_rank;

        double start_time;

        Test1stSolver *saatest = NULL;
        
        void updIt(Status,int);
        void print(Status,int);
        void estimators();
        void testProblem(int);

    public:

        SAA(const Input & input, const Instance & instance) :  Algorithm(input,instance) {
            var_lb = 0.0;
            var_ub = 0.0;
            var_gap = 0.0;
            stat_lb = 0.0;
            stat_ub = 0.0;
            stat_gap = 0.0;
            var_nb_bnb = 0.0;
            nb_notopt = 0;
            best_problem = -1;

            time_lb = 0.0;
            time_eval = 0.0;

            nb_enter_aux = 0.0;
            nb_improv_aux = 0.0;

            lb_vector.resize(instance.getNbProblems());
            ub_vector.resize(instance.getNbProblems());
            gap_vector.resize(instance.getNbProblems());
            nb_bnb_vector.resize(instance.getNbProblems());
            saatest = test_stochastic;

            avg_rank.resize(instance.getNbMaxRankSchools()+1);
        }

        ~SAA(){}

        void solve();

        virtual std::string write() const {
            std::string output = "";
            for(int i = 0; i < instance.getNbProblems(); ++i){
                output += std::to_string(lb_vector[i]) + std::string(";") + 
                          std::to_string(ub_vector[i]) + std::string(";") + 
                          std::to_string(gap_vector[i]) + std::string(";") +
                          std::to_string(nb_bnb_vector[i]) + std::string(";");
            }
            output += std::to_string(status) + std::string(";") + 
                      std::to_string(lb) + std::string(";") + 
                      std::to_string(ub) + std::string(";") + 
                      std::to_string(gap) + std::string(";") + 
                      std::to_string(time_) + std::string(";") +
                      std::to_string(time_lb) + std::string(";") +
                      std::to_string(time_eval) + std::string(";") +
                      std::to_string(nb_bnb) + std::string(";") +
                      std::to_string(var_nb_bnb) + std::string(";") + 
                      std::to_string(stat_lb) + std::string(";") + 
                      std::to_string(stat_ub) + std::string(";") + 
                      std::to_string(stat_gap) + std::string(";") + 
                      std::to_string(var_lb) + std::string(";") + 
                      std::to_string(var_ub) + std::string(";") + 
                      std::to_string(var_gap) + std::string(";") + 
                      std::to_string(nb_notopt) + std::string(";");

            if(input.getTransform() == Input::Transform::NonStrategic){
                // Nonstrategic evaluation
                output += std::to_string(ub_test) + std::string(";") + std::to_string(var_test) + std::string(";") + std::to_string(time_test) + std::string(";");
                for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k) 
                    output += std::to_string(avg_stu_per_rank[k]) + std::string(";");
                output += std::to_string(nb_enter) + std::string(";") + std::to_string(nb_improv) + std::string(";");
            }
            if(input.getTransform() == Input::Transform::NonStrategic || input.getTransform() == Input::Transform::Simple){
                // Simple reported evaluation
                output += std::to_string(ub_test_strat_simple) + std::string(";") + std::to_string(var_test_strat_simple) + std::string(";") + std::to_string(time_test_strat_simple) + std::string(";");
                for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k) 
                    output += std::to_string(avg_stu_per_rank_strat_simple[k]) + std::string(";");
                output += std::to_string(nb_enter_strat_simple) + std::string(";") + std::to_string(nb_improv_strat_simple) + std::string(";");
            }

            // Strategic reported evaluation
            output += std::to_string(ub_test_strat) + std::string(";") + std::to_string(var_test_strat) + std::string(";") + std::to_string(time_test_strat) + std::string(";");
            for(int k = 0; k <= instance.getNbMaxRankSchools(); ++k) 
                output += std::to_string(avg_stu_per_rank_strat[k]) + std::string(";");
            output += std::to_string(nb_enter_strat) + std::string(";") + std::to_string(nb_improv_strat) + std::string(";");
            return output;
        }

        virtual std::string writeHead() const {
            std::string output = "";
            for(int i = 0; i < instance.getNbProblems(); ++i){
                output += "LB_" + std::to_string(i+1) + ";" +
                          "UB_" + std::to_string(i+1) + ";" +
                          "GAP_" + std::to_string(i+1) + ";" +
                          "NB_BNB_" + std::to_string(i+1) + ";";
            }
            output += "STATUS;LB;UB;GAP;TIME;TIME_LB;TIME_EVAL;NB_BNB;VAR_NB_BNB;STAT_LB;STAT_UB;STAT_GAP;VAR_LB;VAR_UB;VAR_GAP;NB_NOT_OPT;";

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