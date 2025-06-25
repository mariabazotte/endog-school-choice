#include "localsearch.hpp"
#include "../test1st/test1stsolver.hpp"

LocalSearchSolver::LocalSearchSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, bool warmstart, bool used_in_ws): 
                            Solver(input,instance,pr,name,env,warmstart,false), used_in_ws(used_in_ws) {
    count_tt = instance.getNbSchools();

    if(t_ == NULL) t_ = new double[count_tt];
    tested_t_ = new double[count_tt];
    best_it_t_ = new double[count_tt];

    if(test_1st == NULL) {
        test_1st = new Test1stSolverDA(input,instance,pr,"Local Search Test Capacity -- Problem " + std::to_string(pr)); 
        test_1st->setVerbose(0);
    }
    
}

void LocalSearchSolver::solve(){
    start_time = time(NULL); 

    // Init information
    it = 0;
    lb = 0.0;

    // Initial feasible solution: no extra capacity
    std::fill(t_,t_+count_tt,0.0);
    test_1st->solve(t_,input.getTransform());
    ub = test_1st->getUB();
    
    // Best solution from current iteration
    ub_best_it = std::numeric_limits<double>::infinity();
    std::fill(best_it_t_,best_it_t_+count_tt,0.0);

    // Tested solution at current operation at current iteration
    std::fill(tested_t_,tested_t_+count_tt,0.0);

    // Number of seats distributed to current t_
    double nb_total_capa = 0.0;

    bool improv = true;
    while(improv == true){
        ++it;
        if(verbose >= 1) std::cout << "Iteration: " << it << "\t UB:" << ub << std::endl;
        ub_best_it = std::numeric_limits<double>::infinity();

        // Operation of including a capacity at each school
        if(nb_total_capa < instance.getNbMaxExtraSeats()){
            for(int c = 0; c < instance.getNbSchools(); ++c){
                // Exchange from current solution
                std::copy(t_,t_+count_tt,tested_t_);
                tested_t_[c] = t_[c]+1;

                // Update
                upd_iteration();
                if(time_limit_reached()) break;
            }
        }

        // Operation of exchanging capacities of two schools
        for(int c = 0; c < instance.getNbSchools()-1; ++c){
            for(int cc = c+1; cc < instance.getNbSchools(); ++cc){
                if(t_[c] > t_[cc] || t_[c] < t_[cc]){
                    // Exchange from current solution
                    std::copy(t_,t_+count_tt,tested_t_);
                    tested_t_[c] = t_[cc];
                    tested_t_[cc] = t_[c];

                    // Update
                    upd_iteration();
                    
                    if(time_limit_reached()) break;
                }
            }
        }

        // Operation of removing a capacity at each school
        // for(int c = 0; c < instance.getNbSchools(); ++c){
        //     if(t_[c]-1.0 >= 0.0){
        //         // Exchange from current solution
        //         std::copy(t_,t_+count_tt,tested_t_);
        //         tested_t_[c] = t_[c]-1;

        //         // Update
        //         upd_iteration();
        //     }
        //     if(time_limit_reached()) break;
        // }

        // Operation of removing a capacity from one school and given it to another
        // for(int c = 0; c < instance.getNbSchools(); ++c){
        //     for(int cc = 0; cc < instance.getNbSchools(); ++cc){
        //         if(c != cc && (t_[c]+1.0 <= instance.getNbMaxExtraSeats()) && (t_[cc]-1.0 >= 0.0)){
        //             // Exchange from current solution
        //             std::copy(t_,t_+count_tt,tested_t_);
        //             tested_t_[c] = t_[c]+1;
        //             tested_t_[cc] = t_[cc]-1;
        //             // Update
        //             upd_iteration();
        //         }
        //     }
        // }

        if(ub_best_it <= ub - 1e-6){
            ub = ub_best_it;
            std::copy(best_it_t_,best_it_t_+count_tt,t_);
            nb_total_capa = 0.0;
            for(int c = 0; c < instance.getNbSchools(); ++c) nb_total_capa += t_[c];
        }else{
            improv = false;
            status = Status::Optimal;
        }
        if(it >= nb_max_it) {
            improv = false;
            status = Status::Iteration_Limit;
        }
        if(time_limit_reached()){
            improv = false;
            status = Status::Time_Limit;
        }
    }
    // We always get a feasible solution
    nb_sol = 1;

    // Getting again best solution
    test_1st->solve(t_,input.getTransform());   

    time_ = time(NULL) - start_time;
    if(verbose >= 1) print();
}

void LocalSearchSolver::upd_iteration(){
    // Solve mi or ordered algorithm + da algorithm
    test_1st->solve(tested_t_,input.getTransform());

    // Save if best solution from iteration
    if(test_1st->getUB() <= ub_best_it - 1e-6){
        ub_best_it = test_1st->getUB();
        std::copy(tested_t_,tested_t_+count_tt,best_it_t_);
    }
}