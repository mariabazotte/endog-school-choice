
#include "simulannealing.hpp"
#include "../test1st/test1stsolver.hpp"

SimulAnnealingSolver::SimulAnnealingSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, bool warmstart): 
                            Solver(input,instance,pr,name,env,warmstart,false){
    count_tt = instance.getNbSchools();

    if(t_ == NULL) t_ = new double[count_tt];
    tested_t_ = new double[count_tt];
    current_t_ = new double[count_tt];

    if(test_1st == NULL) {
        test_1st = new Test1stSolverDA(input,instance,pr,"Simul Annealing Test Capacity -- Problem " + std::to_string(pr)); 
        test_1st->setVerbose(0);
    }
    
}

void SimulAnnealingSolver::solve(){
    start_time = time(NULL); 

    // Init information
    it = 0;
    lb = 0.0;

    // Initial feasible solution: no extra capacity
    std::fill(t_,t_+count_tt,0.0);
    test_1st->solve(t_,input.getTransform());
    ub = test_1st->getUB();

    // Current solution 
    std::copy(t_, t_ + count_tt, current_t_);
    ub_current = ub;

    // Tested solution at current operation at current iteration
    std::fill(tested_t_,tested_t_+count_tt,0.0);

    std::default_random_engine rng(0.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::uniform_int_distribution<int> school_dist(0, instance.getNbSchools() - 1);

    bool stop = false;
    while(!stop){
        ++it;
        for(int inner = 0; inner < iter_per_temp; ++inner){
            // Generate new tested solution
            std::copy(current_t_, current_t_ + count_tt, tested_t_);
            double nb_total_capa = 0.0;
            for(int c = 0; c < instance.getNbSchools(); ++c) nb_total_capa += current_t_[c];
            
            // Ramdomly choose operation
            int op;
            if(nb_total_capa == 0) op = 0;
            else if(nb_total_capa >= instance.getNbMaxExtraSeats()) op = 1 + (it % 2); // Total capacity - use only operations 2 and 3
            else op = it % 3; // Use any operation

            if(op == 0) { // Operation of including a capacity
                int c = school_dist(rng);
                if(tested_t_[c] + 1.0 <= instance.getNbMaxExtraSeats()) {
                    tested_t_[c] += 1.0;
                }
            } else if(op == 1) { // Operation of exchanging capacities of two schools
                int c = school_dist(rng);
                int cc = school_dist(rng);
                while(c == cc) cc = school_dist(rng);
                if(c != cc && (tested_t_[c] != tested_t_[cc])) {
                    int val_cc = tested_t_[cc];
                    int val_c = tested_t_[c];
                    tested_t_[c] = val_cc;
                    tested_t_[cc] = val_c;
                }
            } else if(op == 2) { // Operation of removing a capacity from one school and given it to another
                int c = school_dist(rng);
                int cc = school_dist(rng);
                while(c == cc) cc = school_dist(rng);
                if(c != cc && tested_t_[cc] - 1.0 >= 0.0 && tested_t_[c] + 1.0 <= instance.getNbMaxExtraSeats()) {
                    tested_t_[c] += 1.0;
                    tested_t_[cc] -= 1.0;
                }
            }

            // Evaluate tested solution
            test_1st->solve(tested_t_, input.getTransform());
            double ub_tested = test_1st->getUB();
            double delta = ub_tested - ub_current;

            if(delta <= 1e-6 || uniform(rng) <= std::exp(-delta / T)) {
                // Accept new solution
                std::copy(tested_t_, tested_t_ + count_tt, current_t_);
                ub_current = ub_tested;

                if(ub_current <= ub - 1e-6){
                    ub = ub_current;
                    std::copy(current_t_, current_t_ + count_tt, t_);
                }
            }

            if(time_limit_reached()) break;
        }

        if(verbose >= 1) std::cout << "Iteration: " << it << "\t UB:" << ub << std::endl;

        if(time_limit_reached()){
            stop = true;
            status = Status::Time_Limit;
        }
        if(it >= nb_max_it){
            stop = true;
            status = Status::Iteration_Limit;
        }

        T *= alpha;
        if(T < 0.1) {
            stop = true;
            status = Status::Optimal;
        }
    }
    // We always get a feasible solution
    nb_sol = 1;

    // Getting again best solution
    test_1st->solve(t_,input.getTransform());   

    time_ = time(NULL) - start_time;
    if(verbose >= 1) print();
}
