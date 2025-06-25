#ifndef INSTANCE_HPP
#define INSTANCE_HPP

#include <map>
#include <vector>
#include <list>
#include <unordered_map>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "../input/input.hpp"

typedef std::pair<int,int> feas_pair;

class Instance {
    private:
        const Input & input;                       /* Input. */
        
        std::string instance_file;                 /* File with instance information. */
        int nb_max_rank_schools;                   /* Number maximum of schools each student can rank. */
        int nb_saa_problems;                       /* Number of SAA problems. */
        int nb_saa_scenarios;                      /* Number of scenarios for each SAA problems. */
        int nb_saa_val_scenarios;                  /* Number of scenarios for the SAA validation test. */

        double seed;                               /* Seed for random data generation (scenarios). */
        std::default_random_engine rd_engine;      /* Random engine for data generation (scenarios).*/ 

        std::string name;                          /* Name of the instance. */ 

        int nb_students;                           /* Number of students. */
        int nb_schools;                            /* Number of schools. */
        int nb_max_extra_capacity;                 /* Number maximum of extra capacity to be installed. */

        double critical_tstudent;                  /* Critical value for t-student distribution. */
        double critical_normal;                    /* Critical value for normal distribution. */

        // TODO : Take location and value from instance file
        const double euler_cte = 0.57721;           /* Euler-mascheroni's constant. */
        double location_ = 0.0;                     /* Location parameter of utility distribution.*/
        double scale_ = 4.0;                        /* Scale parameter of utility distribution. */
        std::extreme_value_distribution<> ext_val;  /* Distribution. */
        std::uniform_real_distribution<> uniform;   /* Distribution. */
        std::normal_distribution<> normal;          /* Distribution. */

        std::map<std::string,int> id_schools;                /* Map name school and its id. */
        std::map<std::string,int> id_students;               /* Map name student and its id. */
 
        std::vector<int> schools_capacity;                   /* schools_capacity[c]: Initial capacity of school c. */    
        std::vector<feas_pair> all_feas_pairs;               /* List with all possible pairs of schools and students. */
        std::map<feas_pair,double> scores;                   /* scores[{st,c}]: Score of feasible pair {st,c}. */    
        std::map<feas_pair,double> utilities;                /* utilities[{st,c}]: Utility of feasible pair {st,c}. */  
        std::map<feas_pair,double> probab;                   /* probab[{st,c}]: Probability of feasible pair {st,c} (with initial capacity). */
        std::vector<std::map<feas_pair,double>> inc_probab;  /* inc_probab[j][{st,c}]: Increase in the probability of feasible pair {st,c} when j extra capacity is installed. */
        
        std::vector<std::vector<int>> nb_feasible_pairs;                              /* nb_feasible_pairs[pr][s]: Number of feasible pairs {st,c} in scenario s of problem pr. */
        std::vector<std::vector<std::vector<feas_pair>>> feasible_pairs;              /* feasible_pairs[pr][s][i]: Feasible pair {st,c} (student,school) in scenario s of problem pr. (obs: it includes the unassigned option) */             
        std::vector<std::vector<std::vector<double>>> random_utilities;               /* random_utilities[pr][s][i]: Utility of feasible pair i in scenario s of problem pr. */
        std::vector<std::vector<std::vector<std::vector<int>>>> feas_pair_by_school;  /* feas_pair_by_school[pr][s][c]: List of feasible pairs of school c in scenario s of problem pr. (obs: in order) */
        std::vector<std::vector<std::vector<std::vector<int>>>> feas_pair_by_student; /* feas_pair_by_student[pr][s][st][i]: List of feasible pairs of student st in scenario s of problem pr. (obs: in order) (obs: it does not include unassigned option) */
        std::vector<std::vector<std::vector<int>>> not_assigned_pair;                 /* not_assigned_pair[pr][s][st]: Unassigned option. */
        std::vector<std::vector<std::vector<std::map<int,int>>>> students_preference; /* students_preference[pr][s][st][i]: Rank (index) of feasible pair i by student st in scenario s of problem pr. (obj: starts at 0) (obs: it includes unassigned option) */
        std::vector<std::vector<std::vector<std::map<int,int>>>> schools_preference;  /* schools_preference[pr][s][c][i]: Rank (index) of feasible pair i by school c in scenario s of problem pr. (obj: starts at 0) */

        // Only available when using enumeration as strategic modelling
        std::vector<std::vector<std::vector<int>>> number_ranks;                                /* number_ranks[pr][s][st]: number of possible ranks. */
        std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> enumerated_ranks;  /* enumerated_ranks[pr][s][st][i]: i^th possible ranking (preference order) set. */
        
        // Define and read instance
        void defineName();
        void readInstance();
        void defineCriticalValues();
        
        // Generate scenarios
        void generateScenarios();
        void generateScenario(int,int);
        bool zeroProbability(const feas_pair &);
        void generateUtility(int, feas_pair &, double &);
        void preprocessing_strategic(std::vector<feas_pair> &,std::map<feas_pair,double> &);
        
        // Enumerate all possible ranks
        void enumAllRanks();
        void enumRank(int,int,int,std::vector<std::vector<int>>&,std::vector<int>&,int,int,int);

        void transf_strategic(int,const double *t_, double **xi_, double **rho_) const;
        void transf_nonstrategic(int,double **xi_, double **rho_) const;

    public:
        Instance(const Input &);
        ~Instance(){}

        std::string getName() const { return name; }
        int getNbProblems() const { return nb_saa_problems; }
        int getNbScenarios() const { return nb_saa_scenarios; }
        int getValidateProblem() const { return nb_saa_problems; }
        int getNbValidateScenarios() const { return nb_saa_val_scenarios; }
        int getNbScenarios(int pr) const { if(pr==nb_saa_problems) return nb_saa_val_scenarios; return nb_saa_scenarios; }

        int getNbStudents() const { return nb_students; }
        int getNbSchools() const { return nb_schools; }
        int getNbMaxExtraSeats() const { return nb_max_extra_capacity; }
        int getNbMaxExtraSeatsWithZero() const { return nb_max_extra_capacity + 1; }
        int getNbMaxRankSchools() const { return nb_max_rank_schools; }

        double getCriticalNormal() const { return critical_normal;}
        double getCriticalTStudent() const { return critical_tstudent;}
        double getEpsilon() const { return 0.01; }

        int getNbFeasPairs(int pr, int s) const { return nb_feasible_pairs[pr][s]; } // With unassigned option
        int getNbAssignFeasPairs(int pr, int s) const { return nb_feasible_pairs[pr][s] - nb_students; } // Without unassigned option
        
        const feas_pair & getFeasiblePair(int pr, int s, int i) const { return feasible_pairs[pr][s][i]; }
        int getSchoolFeasiblePair(int pr, int s, int i) const { return feasible_pairs[pr][s][i].second; }
        int getStudentFeasiblePair(int pr, int s,int i) const { return feasible_pairs[pr][s][i].first; }

        int getCapacity(int c) const { return schools_capacity[c]; }
        double getScore(int pr, int s, int i) const { return scores.at(feasible_pairs[pr][s][i]); }
        double getProbability(int pr, int s, int i, int j) const {  double p = probab.at(feasible_pairs[pr][s][i]) + inc_probab[j].at(feasible_pairs[pr][s][i]);
                                                                    return std::max(1e-3,p); }
        double getUtility(int pr, int s, int i) const { return random_utilities[pr][s][i]; } 
        
        int getUnassignedPair(int pr, int s, int st) const { return not_assigned_pair[pr][s][st]; }
        bool isUnassignedPair(int pr, int s, int st, int i) const { return (i == not_assigned_pair[pr][s][st]); }
        const std::vector<int> & getAssignFeasPairsBySchool(int pr, int s, int c) const { return feas_pair_by_school[pr][s][c]; }
        const std::vector<int> & getAssignFeasPairsByStudent(int pr, int s, int st) const { return feas_pair_by_student[pr][s][st]; } // Without unassigned option
        std::vector<int> getFeasPairsByStudent(int pr, int s, int st) const { std::vector<int> vec(feas_pair_by_student[pr][s][st]); vec.push_back(not_assigned_pair[pr][s][st]); return vec; } // With unassigned option
        
        int getFeasPairStudentRank(int pr, int s, int i) const { return students_preference[pr][s][feasible_pairs[pr][s][i].first].at(i); }
        int getNbAssignFeasPairsByStudent(int pr, int s, int st) const { return feas_pair_by_student[pr][s][st].size(); }
        int getAssignFeasPairByStudentPos(int pr, int s, int st, int pos) const { return feas_pair_by_student[pr][s][st][pos]; }

        int getFeasPairSchoolRank(int pr, int s, int i) const { return schools_preference[pr][s][feasible_pairs[pr][s][i].second].at(i); }
        int getLeastPrefFeasPairBySchool(int pr, int s, int c) const { int p = feas_pair_by_school[pr][s][c].size()-1; return feas_pair_by_school[pr][s][c][p]; }
        int getNbAssignFeasPairsBySchool(int pr, int s, int c) const { return feas_pair_by_school[pr][s][c].size(); }
        int getAssignFeasPairBySchoolPos(int pr, int s, int c, int pos) const { return feas_pair_by_school[pr][s][c][pos]; }
        int getNbAssignFeasPairsLessPrefBySchool(int pr, int s, int i) const { return feas_pair_by_school[pr][s][feasible_pairs[pr][s][i].second].size() - schools_preference[pr][s][feasible_pairs[pr][s][i].second].at(i) - 1; }
        
        /* Returns a list of the feasible pairs' indexes that student st prefers over feasible pair i. (It includes feasible pair i) */
        std::vector<int> getStudentPreference(int pr, int s, int st, int i) const { 
            if(i == not_assigned_pair[pr][s][st]) return getFeasPairsByStudent(pr, s, st); // unassigned option
            auto it = std::find(feas_pair_by_student[pr][s][st].begin(), feas_pair_by_student[pr][s][st].end(), i);
            if(it == feas_pair_by_student[pr][s][st].end()) throw std::invalid_argument("Incorrect argument for feasible pair.");
            std::vector<int> preference(feas_pair_by_student[pr][s][st].begin(),it+1);
            return preference; 
        }

        /* Returns a list of the feasible pairs' indexes that school c prefers over feasible pair i. (It does not include feasible pair i) */
        std::vector<int> getSchoolPreference(int pr, int s, int c, int i) const { 
            auto it = std::find(feas_pair_by_school[pr][s][c].begin(), feas_pair_by_school[pr][s][c].end(), i);
            if(it == feas_pair_by_school[pr][s][c].end()) throw std::invalid_argument("Incorrect argument for feasible pair.");
            std::vector<int> preference(feas_pair_by_school[pr][s][c].begin(),it);
            return preference; 
        }

        /* Returns a list of the feasible pairs' indexes that school c doe not prefer over feasible pair i. (It does not include feasible pair i) */
        std::vector<int> getSchoolNotPreference(int pr, int s, int c, int i) const { 
            auto it = std::find(feas_pair_by_school[pr][s][c].begin(), feas_pair_by_school[pr][s][c].end(), i);
            if(it == feas_pair_by_school[pr][s][c].end()) throw std::invalid_argument("Incorrect argument for feasible pair.");
            std::vector<int> preference(it+1,feas_pair_by_school[pr][s][c].end());
            return preference; 
        }

        /* Given a fixed t_, return the optimal rank (ranks) of a student st in scenario s of problem pr. MI algorithm*/
        double mi_algorithm(int pr, int s, int st, const double * t_, std::map<int,int> & ranks) const;

        /* Given a fixed t_, return the optimal rank (ranks) of a student st in scenario s of problem pr. Simple case : ordering p*u */
        double order_algorithm(int pr, int s, int st, const double * t_, std::map<int,int> & ranks) const;

        /* Given a fixed t_ and a fixed set of feasible pairs, return the matching and its value. */
        double da_algorithm(int pr, int s, const double* t_, std::vector<std::list<int>> feas_pairs, std::unordered_map<int,double> rank_feas_pairs,std::vector<std::list<int>> & matching) const;

        // Available methods when all ranks are enumerated.
        int getNbRanks(int pr, int s, int st) const { return number_ranks[pr][s][st]; }
        int getAssignFeasPairRank(int pr, int s, int st, int r, int p) const { return enumerated_ranks[pr][s][st][r][p]; }
        int getSizeRank(int pr, int s, int st, int r) const { return enumerated_ranks[pr][s][st][r].size(); }
        const std::vector<int> & getRank(int pr, int s, int st, int r) const { return enumerated_ranks[pr][s][st][r]; }
        
        // Return the maximum utility of a rank.
        double getMaxUtilityRank(int pr, int s, int st, int r) const { 
            double sum = 0.0;
            for(int i: enumerated_ranks[pr][s][st][r]) sum += getUtility(pr,s,i);
            return sum;
        }

        // Return true if assign feasible pair i is in the rank r.
        bool isAssignFeasPairInRank(int pr, int s, int st, int r, int i) const {
            auto it = std::find(enumerated_ranks[pr][s][st][r].begin(), enumerated_ranks[pr][s][st][r].end(), i);
            return it != enumerated_ranks[pr][s][st][r].end();
        } 

        void display();
};  

#endif