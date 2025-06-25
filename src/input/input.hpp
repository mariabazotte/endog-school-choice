#ifndef INPUT_HPP
#define INPUT_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <limits>

class Input {
    public:
        enum Solver {
            LConstraints = 0,         /* L-Constraints formulation. */
            ModLConstraints = 1,      /* Modified L-Constraints formulation. */
            CutoffScore = 2,          /* Cutoff Score formulation. */
            ASGHeuristic = 3,         /* Assignment heuristic. */
            Lagrangian = 4,           /* Lagrangian. */
            LocalSearch = 5,          /* Local Search. */
            IntegerLShaped = 6,       /* Integer LShaped. */
            SimulAnnealing = 7        /* Simulated Annealing. */
        };

        enum Transform {
            NonStrategic = 0,         /* Students give their true preferences. */
            Compact = 1,              /* Students strategizy, using compact formulation. */
            Enumeration = 2,          /* Students strategizy, using compact formulation. */
            Simple = 3                /* Students strategizy, approximating it with simple formulation. */
        };

        enum Type {
            Deterministic = 0,        /* Solving deterministic problem. (Expected Value problem.) */
            Stochastic = 1            /* Solving stochastic problem with SAA. */
        };

        enum Test1stSolver {
            Gurobi = 0,               /* Solving problem with fixed capacity using gurobi. */
            DA = 1,                   /* Solving problem with fixed capacity using deferred algorithm. */
        };

    private:
        std::string instance_file;     /* File with instance dataset. */
        int nb_max_rank_schools;       /* Number maximum of schools each student can rank. */
        int nb_max_extra_seats;        /* Number maximum of extra seats. */
        bool obj_access_improv;        /* Objective : 0: access. (unassigned option cost |C|+1 -> total number of schools)
                                                      1: improv. (unassigned option cost |C_s^w|+1 -> size reported preference student)*/

        std::string solution_file;     /* File to write results. */
        std::ofstream solFile;         /* Ofstream to write the results in solution_file. */

        double seed;                   /* Seed for scenario generation. */
        double time_limit;             /* Time limit for optimization. */
        double nb_threads;             /* Number of threads for gurobi. */
        double verbose;                /* Verbose (print configuration). */

        int nbproblemsSAA;             /* Number of problems used in the SAA. */
        int nbscenariosSAA;            /* Number of scenarios for each SAA problem. */
        int nbvalidatescenariosSAA;    /* Number of scenarios for the validation problem. */

        bool warm_start;               /* Use warm start procedure. */

        // For Cutoff Score, L-Constraints, Mod-L-Constraints
        bool lazy_stability;           /* Define stability constraints as lazy. */
        bool branch_tooth;             /* Increase priority on branching on tooth. */
        bool combs;                    /* Use combs. */
        bool linear_relax;             /* 0: Solve mixed integer problem, 1: Solve linear relaxation. */

        // For Lagrangian
        int max_nb_it_lag;                  /* Maximum number of lagrangian iterations. */
        int max_nb_it_without_improv_lag;   /* Maximum number of lag iterations without improvement.*/
        double init_scale_lag;              /* Initial value for scale parameter in subgradient. */
        double min_step_lag;                /* Minumum step size for lagrangian. */
        bool modified_sub;                  /* 1: Solve modified version of the subgradient. */             
        bool init_mult_lp;                  /* 1: Solve linear relaxation before the subgradient to initiaize the duals. */

        // For Integer LShaped
        bool lshaped_callback;
        bool lshaped_alternating;

        // Strategic version
        bool barrier;
        bool resprep; 
        bool largenumstab;

        int int_solver;                /* Type of solver. */
        Solver solver;

        int int_transform;             /* Type of transformation. */
        Transform transform;

        int int_type;                  /* Type of problem. */
        Type type;

        void defaultParams();
        void testParameters();
        void defineSolutionFile();

    public:
        Input() { defaultParams(); }
        Input(int argc, char* argv[]);
        ~Input(){ solFile.close(); }

        // Instance
        std::string getInstanceFile() const { return instance_file; }
        int getNbMaxRankSchools() const { return nb_max_rank_schools; }
        int getNbMaxExtraSeats() const { return nb_max_extra_seats; }
        bool getObjAccessImprov() const { return obj_access_improv; }

        // Scenarios 
        int getNbProblemsSAA() const { return nbproblemsSAA; }
        int getNbScenariosSAA() const { return nbscenariosSAA; }
        int getNbValidateScenariosSAA() const { return nbvalidatescenariosSAA; }

        // General parameters
        double getTimeLimit() const { return time_limit; }
        int getNbThreads() const { return nb_threads; }
        double getSeed() const { return seed; }
        int getVerbose() const { return verbose; }

        // Parameters for all models
        bool useWarmStart() const { return warm_start; }
        bool useLazyStability() const { return lazy_stability; }
        bool useBranchTooth() const { return branch_tooth; }
        bool useCombs() const { return combs; }
        bool getLinearRelax() const { return linear_relax; }

        // Lagrangian
        int getNbMaxItLag() const { return max_nb_it_lag; }
        int getNbMaxItWithoutImprovLag() const { return max_nb_it_without_improv_lag; }
        double getInitScale() const { return init_scale_lag; }
        double getMinStepLag() const { return min_step_lag; }
        bool getModSubUpd() const { return modified_sub; }
        bool getIniMultLP() const { return init_mult_lp; }

        // Integer LShaped
        bool useLShapedCallback() const { return lshaped_callback; }
        bool useAlternating() const { return lshaped_alternating; }

        // Strategy
        bool useBarrier() const { return barrier; }
        bool restPrep() const { return resprep; }
        bool useLargeNumStab() const { return largenumstab; }

        // Defining the problem
        Input::Solver getSolver() const { return solver; }
        Input::Transform getTransform() const { return transform; }
        Input::Type getType() const { return type; }

        // Setters
        void setSolver(Input::Solver a) { solver = a; }
        void setTransform(Input::Transform a) { transform = a; }
        void setType(Input::Type a) { type = a; if(type == Input::Type::Deterministic){ nbproblemsSAA = 1; nbscenariosSAA = 1; }  }

        void display();
        void writeHead(std::string);
        void write(std::string);
};

std::ostream& operator<<(std::ostream&, const Input::Solver &);
std::ostream& operator<<(std::ostream&, const Input::Transform &);
std::ostream& operator<<(std::ostream&, const Input::Type &);

#endif