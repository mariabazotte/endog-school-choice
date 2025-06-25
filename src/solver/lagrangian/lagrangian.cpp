#include "lagrangian.hpp"
#include "../test1st/test1stsolver.hpp"

Lagrangian::Lagrangian(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, bool warmstart): 
                            Solver(input,instance,pr,name,env,warmstart,false), MAX_NB_IT(input.getNbMaxItLag()),
                            MAX_NB_IT_WITHOUT_IMPROV(input.getNbMaxItWithoutImprovLag()), INIT_SCALE(input.getInitScale()), 
                            MIN_STEP(input.getMinStepLag()), MOD_SUB(input.getModSubUpd()){ 
    
    count_stab = 0;

    // Initialize heuristic
    init_heuristic();
    
    // Create lagrangian subproblem
    params();
    create_1st();
    create_2st();  

    // Initialize info lagrangian
    init();
    
    // Initialize multipliers with dual lp
    time_lp = 0.0;
    if(input.getIniMultLP() == true){
        time_lp = time(NULL);
        ModLConstraintsSolver* lp_solver = new ModLConstraintsSolver(input,instance,pr,name,env,true,true);
        // solve linear relaxation of modified l-constraints model
        lp_solver->solve();

        // take dual values of the constraints
        const double* dual_ = lp_solver->getDual();
        std::copy(dual_,dual_+count_stab,stability_mult.begin());

        delete lp_solver;
        time_lp = time(NULL) - time_lp;
    }
}

Lagrangian::~Lagrangian(){
    if(curr_t_) delete[] curr_t_;
    if(zero_t_) delete[] zero_t_;
    if(last_t_) delete[] last_t_;
    if(init_coeff_y) delete[] init_coeff_y;
    if(curr_coeff_y) delete[] curr_coeff_y;
    if(tooth_) delete[] tooth_;
}

// First iteration is already the lph, so there is no need for the lph warm-start
void Lagrangian::warm_start(){
    Solver::warm_start_capazero();
    // Taking solution from warm-start
    ub = test_1st->getUB();
    std::fill(t_,t_+count_tt,0.0);
    ++nb_sol;
}

void Lagrangian::params(){
    Solver::params();
    model->set(GRB_IntParam_OutputFlag,0);
}

void Lagrangian::index_2st(int s){
    Solver::index_2st(s);

    for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        
        int pref_c = instance.getFeasPairSchoolRank(pr,s,i);
        int nb_pairs = instance.getNbAssignFeasPairsBySchool(pr,s,c);

        if((pref_c+1) < nb_pairs){
            for(int nn = pref_c+1; nn < nb_pairs; ++nn){
                int ii = instance.getAssignFeasPairBySchoolPos(pr,s,c,nn);
                inv_index_stab.push_back({s,i,ii});
                ++count_stab;
            }
        }
    }

    // Auxiliar parameter to keep value of tooth
    tooth_ = new double[count_y];
    std::fill(tooth_,tooth_+count_y,0.0);
}

void Lagrangian::variables_2st(){
    Solver::variables_2st();

    // Init auxiliars for objective (subproblem)
    init_sub_coeff();
}

void Lagrangian::init_sub_coeff(){
    init_coeff_y = new double[count_y];
    curr_coeff_y = new double[count_y];

    std::fill(init_coeff_y,init_coeff_y+count_y,0.0);
    std::fill(curr_coeff_y,curr_coeff_y+count_y,0.0);
    
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i){
            init_coeff_y[index_y[s][i]] = ((double)(instance.getFeasPairStudentRank(pr,s,i) + 1.0))/((double)instance.getNbScenarios(pr));
        }
        if(input.getObjAccessImprov() == false){
            for(int st = 0; st < instance.getNbStudents(); ++st){
                init_coeff_y[index_y[s][instance.getUnassignedPair(pr,s,st)]] = ((double)(instance.getNbSchools() + 1.0))/((double)instance.getNbScenarios(pr));
            }
        }
    }
}

void Lagrangian::write_sol(){}

////////////////////////////////////////////////////////////////////////////////////////////

void Lagrangian::init(){
    stability_mult.resize(count_stab,0.0);
    stability_grad.resize(count_stab,0.0);
    if(MOD_SUB == true) stability_dire.resize(count_stab,0.0);

    lb = 0.0;                                                          // initial lower bound
    ub = instance.getNbStudents()*(instance.getNbMaxRankSchools()+1);  // solution where all students are unassigned
    gap = (ub-lb)/ub;

    curr_lb = -std::numeric_limits<double>::infinity(); // For the current value of the subproblem
    curr_ub =  std::numeric_limits<double>::infinity(); // For the current real objective of the subproblem
    curr_nb_sol = 0;
    curr_t_ = new double[instance.getNbSchools()]; // Current solution of subproblem in node

    status = Status::Not_Solved;
    status_sub = Status::Not_Solved;

    // Initialize Lagrangian info
    nb_it = 0;
    nb_it_without_improv = 0;
    nb_total_it_without_improv = 0;
    step = 0.0;
    scale = INIT_SCALE;
    norm = 0.0;

    sub_feasible = true;
    sub_respect_slackness = true;
}

void Lagrangian::upd_mult(){
    for(int i = 0; i < count_stab; ++i){
        stability_mult[i] = (std::max(0.0, stability_mult[i] - step*stability_grad[i]));
    }
}

void Lagrangian::upd_grad(){
    sub_feasible = true;
    sub_respect_slackness = true;
    nb_infeasible = 0;

    for(int s = 0; s < instance.getNbScenarios(pr); ++s) tooth(s);

    for(int l = 0; l < count_stab; ++l){
        const auto& [s, i, ii] = inv_index_stab[l];

        stability_grad[l] = tooth_[index_y[s][i]] - getY_(s,ii);

        if(stability_grad[l] <= -0.0001) {
            sub_feasible = false;
            ++nb_infeasible; 
        }
        if((stability_grad[l] <= -0.0001 || stability_grad[l] >= 0.0001) 
        && (stability_mult[l] <= -0.0001 || stability_mult[l] >= 0.0001))
            sub_respect_slackness = false;
    } 
}

void Lagrangian::tooth(int s){
    for(int st = 0; st < instance.getNbStudents(); ++st){
        double cumulated = 0.0;
        
        // assign pairs
        int nb_assignpairs = instance.getNbAssignFeasPairsByStudent(pr,s,st);
        if(nb_assignpairs > 0){
            for(int n = 0; n < nb_assignpairs; ++n){
                int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
                cumulated += getY_(s,i);
                tooth_[index_y[s][i]] = cumulated;
            }
        }

        // unassign pair
        int i = instance.getUnassignedPair(pr,s,st);
        cumulated += getY_(s,i);
        tooth_[index_y[s][i]] = cumulated;
    }
}

void Lagrangian::upd_norm(){
    norm = 0.0;
    for(int i = 0; i < count_stab; ++i)
        norm += std::pow(stability_grad[i],2);
}

void Lagrangian::upd_dire(){  
    // grad*last_dire 
    double grad_last_dire = 0.0;
    for(int i = 0; i < count_stab; ++i)
        grad_last_dire += stability_grad[i]*stability_dire[i];

    double theta = 0.0;
    if(grad_last_dire <= -1e-5){
        
        // last_dire*last_dire
        double norm_last_dire = 0.0;
        for(int i = 0; i < count_stab; ++i)
            norm_last_dire += std::pow(stability_dire[i],2);
        
        // Until now, norm is grad*grad
        theta = norm/norm_last_dire;
    }

    // upd dire
    for(int i = 0; i < count_stab; ++i)
        stability_dire[i] = stability_grad[i] + theta*stability_dire[i];
    
    // Norm will be grad*dire
    norm = 0.0;
    for(int i = 0; i < count_stab; ++i)
        norm += stability_grad[i]*stability_dire[i];
    
    // Use dire instead of grad to update the multipliers
    std::copy(stability_dire.begin(),stability_dire.end(),stability_grad.begin());
}

void Lagrangian::upd_RL_lb(){
    if(curr_lb >= lb + 1e-6){
        lb = curr_lb;
        nb_it_without_improv = 0;
        nb_total_it_without_improv = 0;
    }else{
        ++nb_it_without_improv;
        ++nb_total_it_without_improv;
    } 
}

void Lagrangian::upd_RL_ub(){
    if(sub_feasible == true && sub_respect_slackness == true){
        lb = curr_ub;
        ub = curr_ub;
        ++nb_sol;
        for(int c = 0; c < count_tt; ++c) t_[c] = curr_t_[c];
        if(curr_ub < ub) time_best_ub = time(NULL) - start_time;
    }else if(sub_feasible == true){ 
        if(curr_ub < ub){
            ++nb_sol;
            ub = curr_ub;
            for(int c = 0; c < count_tt; ++c) t_[c] = curr_t_[c];
            if(curr_ub < ub) time_best_ub = time(NULL) - start_time;
        }
    }
    solve_heuristic(curr_t_);
    // else if(lb <= curr_lb + 1e-6) solve_heuristic(curr_t_); 
}

void Lagrangian::upd_scale(){
    if(scale >= 1e-3 && nb_it_without_improv >= MAX_NB_IT_WITHOUT_IMPROV){
        nb_it_without_improv = 0;
        scale = scale/2.0;
    }
}

void Lagrangian::upd_step(){
    step = scale*(ub - curr_lb)/norm;
}

void Lagrangian::upd_RL_gap(){
    gap = (ub - lb)/ub;
}

void Lagrangian::solve_sub(){
    model->set(GRB_DoubleParam_TimeLimit, get_timelimit());

    // Update coefficients according to current dual
    upd_sub_coeff();

    // Solve model
    model->update();
    model->optimize();

    // Get the solution
    status_sub = statusFromGurobi(model->get(GRB_IntAttr_Status));
    if(status_sub == Status::Optimal) { 
        upd_solution();
    }
}

void Lagrangian::upd_sub_coeff(){
    // Initialize
    std::copy(init_coeff_y, init_coeff_y + count_y, curr_coeff_y);

    // Dual part
    for(int l = 0; l < count_stab; ++l){
        const auto& [s, i, ii] = inv_index_stab[l];

        int st = instance.getStudentFeasiblePair(pr,s,i);
        int pref_st = instance.getFeasPairStudentRank(pr,s,i);
        for(int nnn = 0; nnn <= pref_st; ++nnn){
            int iii = instance.getAssignFeasPairByStudentPos(pr,s,st,nnn);
            curr_coeff_y[index_y[s][iii]] -= stability_mult[l];
        }
        curr_coeff_y[index_y[s][ii]] += stability_mult[l];
    }
    
    model->set(GRB_DoubleAttr_Obj,y,curr_coeff_y,count_y);
}

void Lagrangian::upd_solution(){
    curr_lb = model->get(GRB_DoubleAttr_ObjVal);
    curr_nb_sol = model->get(GRB_IntAttr_SolCount);
    
    if(curr_nb_sol > 0){
        if(curr_t_) delete[] curr_t_;
        curr_t_ = model->get(GRB_DoubleAttr_X,t,count_tt);
        for(int c = 0; c < count_tt; ++c) curr_t_[c] = std::round(curr_t_[c]); // avoid integer feasibility errors
        
        if(y_) delete[] y_;
        y_ = model->get(GRB_DoubleAttr_X,y,count_y);
        for(int l = 0; l < count_y; ++l) y_[l] = std::round(y_[l]); // avoid integer feasibility errors

        transf->get_variables();
    }

    upd_curr_bounds();
}

void Lagrangian::upd_curr_bounds(){
    curr_ub = curr_lb;
    for(int l = 0; l < count_y; ++l)
        curr_ub -= y_[l]*(curr_coeff_y[l] - init_coeff_y[l]);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Lagrangian::solve(){
    start_time = time(NULL) - time_lp;
    if(warmstart == true) warm_start();
    
    model->update();

    std::string output;
    bool stop = false;
    while(stop == false){
        ++nb_it;

        solve_sub();
        if(status_sub == Status::Infeasible || status_sub == Status::Infeasible_or_Unbounded){
            status = Status::Infeasible;
            output = "Infeasible subproblem.\n";
            break;
        }

        // Update info iteration
        upd_RL_lb();
        upd_grad(); 
        upd_norm();
        if(MOD_SUB == true) upd_dire();
        upd_scale();
        upd_step();
        
        upd_RL_ub();  
        upd_RL_gap();
        
        if(get_timelimit() <= 1e-3){
            stop = true;
            status = Status::Time_Limit;
            output = "Time limit was reached.\n";
        }
        else if(sub_feasible == true && sub_respect_slackness == true){
            stop = true;
            status = Status::Optimal;
            output = "Optimal by slackness.\n"; 
        }
        else if(gap <= 1e-4){
            stop = true;
            status = Status::Optimal; 
            output = "Optimal by gap.\n"; 
        }
        else if(lb >= (ub - 1.0/(double)instance.getNbScenarios(pr))){
            lb = ub;
            gap = (ub-lb)/ub;
            stop = true;
            status = Status::Optimal;
            output = "Optimal by gap.\n";
        }
        else if(nb_it >= MAX_NB_IT){ 
            stop = true;
            status = Status::Iteration_Limit;
            output = "Maximum number of iterations was reached.\n"; 
        }
        else if(nb_total_it_without_improv >= 2*MAX_NB_IT_WITHOUT_IMPROV){
            stop = true;
            status = Status::Iteration_Limit;
            output = "Maximum number of iterations was reached.\n"; 
        }
        else if(step < MIN_STEP/instance.getNbScenarios(pr)){
            stop = true; 
            status = Status::Iteration_Limit;
            output = "Minimum step size was reached.\n";
        }
        if(verbose >= 3) std::cout << "it: " << nb_it << " lb: " << lb << " ub: " << ub << std::endl; 
        if(stop == true) break;
        
        // Update multiplier if stop = false
        upd_mult();
    }

    time_ = time(NULL) - start_time;

    if(verbose >= 1) print_specific(output);
    if(verbose >= 1) print();
}

void Lagrangian::print_specific(std::string & output){
    std::cout << "Nb it:" << nb_it <<  " RL lb: " << lb << " RL ub: " << ub << " Time: " << time_ << ". " <<  output;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void Lagrangian::init_heuristic(){
    // For initial UB (zero extra capacity)
    zero_t_ = new double[instance.getNbSchools()];
    std::fill(zero_t_,zero_t_+instance.getNbSchools(),0.0);

    // Last tested extra capacity
    last_t_ = new double[instance.getNbSchools()];
    std::fill(last_t_,last_t_+instance.getNbSchools(),0.0);

    // For heuristic
    matching.resize(instance.getNbScenarios(pr));
    feas_pairs.resize(instance.getNbScenarios(pr));
    rank_feas_pairs.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        feas_pairs[s].resize(instance.getNbStudents());
        matching[s].resize(instance.getNbSchools()+1);    
    }

    init_true_ranks_heuristic();
    init_feas_pairs_heuristic();
}

void Lagrangian::init_true_ranks_heuristic(){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            int r = 1;
            for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){ 
                rank_feas_pairs[s][i] = r;
                ++r;
            } 
            int i = instance.getUnassignedPair(pr,s,st);
            if(input.getObjAccessImprov() == false) rank_feas_pairs[s][i] = (instance.getNbSchools()+1);
            else rank_feas_pairs[s][i] = r;
        }
    }
}

void Lagrangian::init_feas_pairs_heuristic(){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            for(int i: instance.getAssignFeasPairsByStudent(pr,s,st))
                feas_pairs[s][st].push_back(i);
        }
    }
}

void Lagrangian::solve_heuristic(const double * tested_t_){
    bool equal_zero = std::equal(tested_t_,tested_t_+instance.getNbSchools(),zero_t_);
    bool equal_last = std::equal(tested_t_,tested_t_+instance.getNbSchools(),last_t_);

    if(equal_last == false || equal_zero == true ){
        double heur_ub = 0.0;
        bool feasible = true;
        for(int s = 0; s < instance.getNbScenarios(pr); ++s){
            for(int c = 0; c < instance.getNbSchools() + 1; ++c) matching[s][c].clear(); // clear previous matching
            
            heur_ub += instance.da_algorithm(pr,s,tested_t_,feas_pairs[s],rank_feas_pairs[s],matching[s])/(double)instance.getNbScenarios(pr);
            
            for(int i: matching[s][instance.getNbSchools()]){
                if(getY(s,i).get(GRB_DoubleAttr_UB) <= 0.0001) feasible = false;
            }

            if(feasible == false) break;
        }
        if(heur_ub < ub && feasible == true){
            ++nb_sol;
            ub = heur_ub;
            for(int c = 0; c < count_tt; ++c) t_[c] = tested_t_[c];
        }
        for(int c = 0; c < instance.getNbSchools(); ++c) last_t_[c] = tested_t_[c];
        
        if(heur_ub < ub && feasible == true){
            time_best_ub = time(NULL) - start_time;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
