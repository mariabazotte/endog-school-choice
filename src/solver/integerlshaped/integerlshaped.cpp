#include "integerlshaped.hpp"
#include "../test1st/test1stsolver.hpp"
#include "../asgheuristic/asgheuristicsolver.hpp"
#include "../localsearch/localsearch.hpp"

IntegerLShaped::IntegerLShaped(const Input& input,const Instance& instance,int pr,std::string name,GRBEnv* env,bool warmstart): 
                                                        Solver(input,instance,pr,name,env,warmstart,false),
                                                        lshaped_callback(input.useLShapedCallback()),
                                                        alternating_method(input.useAlternating()) {
    // Initialize info
    start_time = 0.0;
    nb_it = 0;
    nb_int_opt_cuts = 0;
    nb_cont_opt_cuts = 0;
    intfeastol = 0.0;
    precision = 1e-4;
    ub_callback = (instance.getNbMaxRankSchools()+1)*instance.getNbStudents();

    if(test_1st == NULL) {
        test_1st = new Test1stSolverDA(input,instance,pr,"Test Capacity -- Problem " + std::to_string(pr));
    }
    test_1st->setVerbose(0);

    params();
    create_1st();    // Create main problem.
    variables_2st();   
    init();

    // Alternating method
    if(alternating_method == true){
        count_zf = 0;
        subproblem = new GRBModel(*env);
        create_sub_2st();
    }
}

IntegerLShaped::~IntegerLShaped(){
    if(sub_x) delete[] sub_x;
    if(sub_y) delete[] sub_y;
    if(sub_z) delete[] sub_z;
    if(sub_f) delete[] sub_f;
    if(constr) delete[] constr;
    if(duals) delete[] duals;
    if(subproblem) delete subproblem;
}

void IntegerLShaped::params(){
    // Params for Main problem.
    model->set(GRB_IntParam_Threads, input.getNbThreads());
    model->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
    model->set(GRB_DoubleParam_MIPGap, 1e-4);
    model->set(GRB_DoubleParam_IntFeasTol, 1e-5);
    model->set(GRB_DoubleParam_FeasibilityTol, 1e-6);
    if(verbose >= 1 && lshaped_callback == false) model->set(GRB_IntParam_OutputFlag, 0);
    else if(verbose >= 1) model->set(GRB_IntParam_OutputFlag, 1);
    else model->set(GRB_IntParam_OutputFlag, 0);
    intfeastol =  model->get(GRB_DoubleParam_IntFeasTol);  
    if(lshaped_callback == true) model->set(GRB_IntParam_LazyConstraints, 1); 
    if(lshaped_callback == true) model->setCallback(this);
}

///////////////////////////////////////////////////////////////////////////////////

void IntegerLShaped::index_1st(){
    index_x.resize(instance.getNbSchools());
    for(int c = 0; c < instance.getNbSchools(); ++c){
        index_x[c].resize(instance.getNbMaxExtraSeats());
        for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
            index_x[c][j] = count_x;
            ++count_x;
        }
    }
    count_tt = instance.getNbSchools();
}

void IntegerLShaped::variables_1st(){
    t = model->addVars(count_tt,GRB_INTEGER);
    x  = model->addVars(count_x,GRB_BINARY);
}

void IntegerLShaped::constraints_1st(){
    int nb_schools = instance.getNbSchools();
    int nb_seats = instance.getNbMaxExtraSeats();

    GRBLinExpr budget = 0;
    for(int c = 0; c < nb_schools; ++c) budget += getT(c); 
    model->addConstr(budget <= nb_seats);

    // Linking x with t
    for(int c = 0; c < instance.getNbSchools(); ++c){
        GRBLinExpr exp = getT(c);
        for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j) exp -= x[index_x[c][j]];
        model->addConstr(exp == 0.0);
    }

    // Defining variables x
    for(int c = 0; c < instance.getNbSchools(); ++c){
        for(int j = 0; j < instance.getNbMaxExtraSeats() - 1; ++j){
            model->addConstr(x[index_x[c][j]] >= x[index_x[c][j+1]]);
        }
    }
}

////////////////////////////////////////////////////////////

// Get value variable t_, and obtain value variable x_
void IntegerLShaped::varx_from_vart(){
    std::fill(x_,x_+count_x,0.0);
    for(int c = 0; c < count_tt; ++c){
        if(t_[c] > 0.001){
            for(int j = 0; j < (int)t_[c]; ++j)
                x_[index_x[c][j]] = 1.0;
        }
    }
}

void IntegerLShaped::warm_start_capazero(){
    if(verbose >= 1) std::cout << "\n### Starting warm-start zero-capacity ###\n\n";
    
    std::fill(t_,t_+count_tt,0.0);
    std::fill(x_,x_+count_x,0.0);
    test_1st->solve(t_,input.getTransform());

    // Setting initial solution
    model->set(GRB_DoubleAttr_Start,t,t_,count_tt);
    model->set(GRB_DoubleAttr_Start,x,x_,count_x);
    
    // Add initial cuts
    curr_2st_value = test_1st->getUB();
    add_integer_cut(); 

    if(verbose >= 1) std::cout << "\n#### Ending warm-start zero-capacity ###\n\n";
}

void IntegerLShaped::warm_start_asg(){
    if(verbose >= 1) std::cout << "\n### Starting warm-start ASG ###\n\n";
    
    asg->solve(); 
    if(asg->getNbSol() > 0){ 
        // Get first stage solutions
        std::copy(asg->getT_(),asg->getT_()+count_tt,t_);
        varx_from_vart();
        
        // Set initial solution
        model->set(GRB_DoubleAttr_Start,t,asg->getT_(),count_tt);
        model->set(GRB_DoubleAttr_Start,x,x_,count_x);

        // Set initial upper bound
        ub = asg->getUB();
        ub_callback = asg->getUB();

        // Add initial cuts
        curr_2st_value = asg->getUB();
        add_integer_cut(); 
    }

    if(verbose >= 1) std::cout << "\n#### Ending warm-start ASG ###\n\n";
}

void IntegerLShaped::warm_start_ls(){
    if(verbose >= 1) std::cout << "\n### Starting warm-start Local Search ###\n\n";

    ls->solve();
    
    // Get first stage solutions
    std::copy(ls->getT_(),ls->getT_()+count_tt,t_);
    varx_from_vart();
    
    // Set initial upper bound
    ub = ls->getUB();
    ub_callback = ls->getUB();

    // Set initial solution
    model->set(GRB_DoubleAttr_Start,t,ls->getT_(),count_tt);
    model->set(GRB_DoubleAttr_Start,x,x_,count_x);

    // Add initial cuts
    curr_2st_value = ls->getUB();
    add_integer_cut(); 

    if(verbose >= 1) std::cout << "\n#### Ending warm-start Local Search ###\n\n";
}

///////////////////////////////////////////////////////////////////////////

void IntegerLShaped::variables_2st(){
    theta = model->addVar(0.0,GRB_INFINITY,1.0,GRB_CONTINUOUS);
}

// Obtain lower bound by assuming the maximum number of 
// extra seats is installed to each school.
void IntegerLShaped::compute_lower_bound(){
    if(input.getTransform() == Input::Transform::NonStrategic){
        double *all_capa_t_ = new double[count_tt];
        std::fill(all_capa_t_,all_capa_t_+count_tt,instance.getNbMaxExtraSeats());
        test_1st->solve(all_capa_t_,input.getTransform());
        lower_bound = test_1st->getUB();
    }else{
        lower_bound = instance.getNbStudents();
    }
    model->addConstr(theta >= lower_bound);
}

/////////////////////////////

void IntegerLShaped::init(){
    if(t_) delete[] t_;
    if(x_) delete[] x_;
    t_ = new double[count_tt];
    x_ = new double[count_x]; 
    theta_ = 0.0;
}

void IntegerLShaped::verify_vars(){
    // Avoid problems with Integer Feasibility tolerance
    for(int i = 0; i < count_x; ++i){
        if(x_[i] <= 0.0 + intfeastol && x_[i] >= 0.0 - intfeastol)
            x_[i] = 0;
        else if(x_[i] <= 1.0 + intfeastol && x_[i] >= 1.0 - intfeastol)
            x_[i] = 1;
    }
    for(int c = 0; c < count_tt; ++c) t_[c] = std::round(t_[c]);
}

////////////////////////////////////////////////////

void IntegerLShaped::solve_main(){
    model->set(GRB_DoubleParam_TimeLimit, get_timelimit());
    model->update();
    model->optimize();

    status = statusFromGurobi(model->get(GRB_IntAttr_Status));
    verify_stat();

    if(status == Status::Time_Limit){
        lb = std::max(lb, model->get(GRB_DoubleAttr_ObjBound));
        return;
    }else if(status == Status::Optimal) {
        lb = model->get(GRB_DoubleAttr_ObjBound);
        nb_sol += 1;

        if(t_) delete[] t_;
        if(x_) delete[] x_;
        
        t_ = model->get(GRB_DoubleAttr_X,t,count_tt);
        x_ = model->get(GRB_DoubleAttr_X,x,count_x);
        theta_ = theta.get(GRB_DoubleAttr_X);
        verify_vars();

        curr_2st_value = 0.0;  
    }
}

void IntegerLShaped::solve_integer_sub(){
    // Solve second stage problem for current first stage decisions
    // MIA or OIA + DAA
    test_1st->solve(t_,input.getTransform());
    curr_2st_value = test_1st->getUB();
}

void IntegerLShaped::separate_integer_cut(GRBLinExpr & cut){
    double curr_sum_x = std::accumulate(x_,x_+count_x,0.0);
    for(int i = 0; i < count_x; ++i){
        if(x_[i] >= 0.999) cut += (curr_2st_value - lower_bound)*x[i];
        else if(x_[i] <= 0.001) cut -= (curr_2st_value - lower_bound)*x[i];
    }
    cut -= (curr_2st_value - lower_bound)*curr_sum_x;
    cut += curr_2st_value;
    cut -= theta;
}

void IntegerLShaped::add_integer_cut(){
    if(curr_2st_value >= theta_ + 1e-6){
        GRBLinExpr cut = 0;
        separate_integer_cut(cut);
        model->addConstr(cut <= 0, "Int_opt_cut_" + std::to_string(nb_int_opt_cuts)); ++nb_int_opt_cuts;
    }
}

////////////////////////////////////////////////////

void IntegerLShaped::get_variables(){
    if(nb_sol > 0){
        if(ub <= ub_ws || warmstart == false){
            if(t_) delete[] t_;
            t_ = model->get(GRB_DoubleAttr_X,t,count_tt);
            for(int c = 0; c < count_tt; ++c) t_[c] = std::round(t_[c]);
        }
    }
}

void IntegerLShaped::solve(){
    start_time = time(NULL); 

    compute_lower_bound();
    if(warmstart == true) warm_start();
    
    model->set(GRB_DoubleParam_TimeLimit, get_timelimit());
    if(lshaped_callback == true){
        model->update();
        model->optimize();

        status = statusFromGurobi(model->get(GRB_IntAttr_Status));
        if(status == Status::Optimal || status == Status::Time_Limit) upd_solution();
        else verify_stat(); 
    } 
    else integer_loop();

    time_ = time(NULL) - start_time;
    if(verbose >= 1) print();
    if(verbose >= 1) std::cout << "Nb int cuts: " << nb_int_opt_cuts << ", Nb cont cuts: " <<  nb_cont_opt_cuts << std::endl;
}

void IntegerLShaped::integer_loop(){
    bool stop = false;
    curr_2st_value = 0.0;
    
    ++nb_it;
    while(stop == false){
        solve_main(); 
        if(status == Status::Infeasible) return;

        bool execute_integer = true;
        if(alternating_method == true){ 
            curr_2st_value = 0.0;
            solve_continuous_sub();
            bool included = add_continuous_cut();
            if(included == true) execute_integer = false;
        }
        if(execute_integer == true){
            curr_2st_value = 0.0;
            solve_integer_sub(); 
            add_integer_cut(); 
        }
        
        if(get_timelimit() <= 0.001){
            status = Status::Time_Limit;
            stop = true;
        }else{
            ub = std::min(ub,curr_2st_value);
            if(ub <= 0.001 && ub >= -0.001){ gap = (ub-lb); } else{ gap = (ub-lb)/ub; } 
        }
        if((gap < precision) && (curr_2st_value < GRB_INFINITY)){
            status = Status::Optimal;
            stop = true;
        }

        if(verbose >= 1){
            if(ub >= GRB_INFINITY && lb <= -GRB_INFINITY){
                printf("i %3d | -infinity | infinity | %12.4f | %3d | %12.4f \n", nb_it, gap, nb_int_opt_cuts, time(NULL)-start_time);
            }else if(ub >= GRB_INFINITY ){
                printf("i %3d | %12.4f | infinity | %12.4f | %3d | %12.4f \n", nb_it, lb, gap, nb_int_opt_cuts, time(NULL)-start_time);
            }else if(lb <= -GRB_INFINITY){
                printf("i %3d | -infinity | %12.4f | %12.4f | %3d | %12.4f \n", nb_it, ub, gap, nb_int_opt_cuts, time(NULL)-start_time);
            } else{
                printf("i %3d | %12.6f | %12.6f | %12.4f | %3d | %12.4f \n", nb_it, lb, ub, gap, nb_int_opt_cuts, time(NULL)-start_time);
            }
        }
        ++nb_it;
    }
    if(verbose >= 1) print();
}

void IntegerLShaped::callback(){
    try{
        if(where == GRB_CB_MIPSOL){
            
            // Get main problem's solution
            if(t_) delete[] t_;
            if(x_) delete[] x_;
            
            t_ = getSolution(t,count_tt);
            x_ = getSolution(x,count_x);
            theta_ = getSolution(theta);
            verify_vars();

            // If objective value is better than the best current solution
            if(theta_ <= ub_callback + 1e-8){
                bool execute_integer = true;
                
                if(alternating_method == true){
                    // Solve continuous subproblem
                    substatus = Status::Not_Solved;
                    curr_2st_value = 0.0;
                    solve_continuous_sub();
                    // std::cout << "cont:" << curr_2st_value << " " << theta_ << std::endl;

                    // Include continuous cut if it is violated
                    if(curr_2st_value >= (theta_ + 1e-3) && substatus == Status::Optimal){
                        GRBLinExpr cut = 0;
                        separate_continuous_cut(cut);
                        addLazy(cut <= 0);
                        ++nb_cont_opt_cuts;
                        execute_integer = false;
                    }
                }

                if(execute_integer == true){                
                    // Solve integer subproblem
                    curr_2st_value = 0.0;
                    solve_integer_sub();
                    // std::cout << "int:" << curr_2st_value << " " << theta_ << std::endl;

                    // Verify if adding new iteger lshaped cut
                    if(curr_2st_value >= (theta_ + 1e-6)){
                        GRBLinExpr cut = 0;
                        separate_integer_cut(cut);
                        addLazy(cut <= 0);
                        ++nb_int_opt_cuts;
                    }

                    // Update upper bound
                    if(ub_callback >= (curr_2st_value + 1e-6)) {
                        ub_callback = curr_2st_value;
                        setSolution(x,x_,count_x);
                        setSolution(theta,curr_2st_value);
                        // std::cout << "sol:" << curr_2st_value  << std::endl;
                    }
                }
            }
        }
    } catch (GRBException e) {
        std::cout << "Error number: " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch (...) {
        std::cout << "Error during callback" << std::endl;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

void IntegerLShaped::create_sub_2st(){
    // Params
    subproblem->set(GRB_IntParam_OutputFlag, 0);
    subproblem->set(GRB_IntParam_Threads, input.getNbThreads());
    subproblem->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());

    // Definition
    index_y.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) index_sub_2st(s);
    variables_sub_2st();
    linear_relax = true; // Command to define strategic part as continuous (for continuous subproblem)
    transf->transf_2st();
    linear_relax = false;
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) constraints_sub_2st(s);
    y_ = new double[count_y];
}

void IntegerLShaped::index_sub_2st(int s){
    Solver::index_2st(s);
    for(int c = 0; c < instance.getNbSchools(); ++c){
        index_zf.insert({{s,c},count_zf});
        ++count_zf;
    }
}

void IntegerLShaped::variables_sub_2st(){
    sub_y = subproblem->addVars(count_y,GRB_CONTINUOUS);
    sub_x = subproblem->addVars(count_x,GRB_CONTINUOUS);
    if(input.getTransform() != Input::Transform::NonStrategic){
        sub_z = subproblem->addVars(count_zf,GRB_CONTINUOUS);
        sub_f = subproblem->addVars(count_zf,GRB_CONTINUOUS);
    }

    // Constraints for continuous cut
    constr = new GRBConstr[count_x];
    duals = new double[count_x];
    for(int c = 0; c < instance.getNbSchools(); ++c){
        for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
            constr[index_x[c][j]] = subproblem->addConstr(getX(c,j) == 0,"cut["+std::to_string(c)+"]["+std::to_string(j)+"]");
        }
    }
}

void IntegerLShaped::constraints_sub_2st(int s){
    int nb_schools = instance.getNbSchools();
    int nb_students = instance.getNbStudents();
    int nb_seats = instance.getNbMaxExtraSeats();

    // Assignment constraints
    for(int st = 0; st < nb_students; ++st){
        GRBLinExpr exp = 0;
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) exp += getY(s,i);
        exp += getY(s,instance.getUnassignedPair(pr,s,st));
        subproblem->addConstr(exp == 1);
    }

    // Capacity constraints
    for(int c = 0; c < nb_schools; ++c){
        GRBLinExpr exp = 0;
        for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j) exp -= getX(c,j);
        for(int i : instance.getAssignFeasPairsBySchool(pr,s,c)) exp += getY(s,i);    
        subproblem->addConstr(exp <= instance.getCapacity(c));
    }

    if(input.getTransform() != Input::Transform::NonStrategic){

        // Cutoff constraints for stability
        for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
            int v = nb_students + 1;
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            int r = instance.getScore(pr,s,i) + nb_students + 1;
            subproblem->addConstr(getZ(s,c) + v*getY(s,i) <= r);
        }

        for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            int st = instance.getStudentFeasiblePair(pr,s,i);

            int v1 = instance.getScore(pr,s,i) + instance.getEpsilon();
            int v2 = nb_students + 1;

            GRBLinExpr exp = -getZ(s,c);
            transf->term_xi(s,i,v1,0.0,exp);
            for(int l: instance.getStudentPreference(pr,s,st,i)) exp -= v2*getY(s,l);
            subproblem->addConstr(exp <= 0.0);
        }
            
        // Waste constraints
        for(int c = 0; c < nb_schools; ++c){
            GRBLinExpr exp = instance.getCapacity(c) - (1-getF(s,c))*(instance.getCapacity(c)+nb_seats);
            for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j) exp += getX(c,j);
            for(int i : instance.getAssignFeasPairsBySchool(pr,s,c)) exp -= getY(s,i);
            subproblem->addConstr(exp <= 0.0);
        }

        for(int c = 0; c < nb_schools; ++c){
            subproblem->addConstr(getZ(s,c) -(nb_students+1)*getF(s,c) <= 0.0);
        }
    }else{
        // Modified l-constraints - stability constraints
        for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            int st = instance.getStudentFeasiblePair(pr,s,i);
    
            for(int ii: instance.getSchoolNotPreference(pr,s,c,i)){
                GRBLinExpr exp = -getY(s,ii);
                for(int l: instance.getStudentPreference(pr,s,st,i)) exp += getY(s,l);
                transf->term_xi(s,i,-1,1,exp);
                subproblem->addConstr(exp >= 0.0);
            } 
        }
    }
}

void IntegerLShaped::solve_continuous_sub(){
    // Modify constraints for cut
    for(int c = 0; c < instance.getNbSchools(); ++c){
        for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
            constr[index_x[c][j]].set(GRB_DoubleAttr_RHS,x_[index_x[c][j]]);
        }
    }
    
    // Solve updated subproblem
    subproblem->set(GRB_DoubleParam_TimeLimit, get_timelimit());
    subproblem->update();
    subproblem->optimize();
    
    // Get solution
    substatus = statusFromGurobi(subproblem->get(GRB_IntAttr_Status));
    if(substatus == Status::Optimal || substatus == Status::Time_Limit){
        // optimal value
        curr_2st_value = subproblem->get(GRB_DoubleAttr_ObjVal);
        // dual variables
        if(duals) delete[] duals;
        duals = subproblem->get(GRB_DoubleAttr_Pi,constr,count_x);
        // ((StrategicTransform*)transf)->print_xi();
        // subproblem->write(name + ".lp");
    }else{
        for(int c = 0; c < instance.getNbSchools(); ++c){
            std::cout << "c:" << c << std::endl;
            for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
                std::cout << x_[index_x[c][j]] << " "; 
            }
            std::cout << std::endl;
        }  
        subproblem->set(GRB_IntParam_OutputFlag, 1);
        subproblem->write(name + ".lp");
        std::cout << "Problem " << substatus << std::endl;
        subproblem->computeIIS();
        subproblem->write(name + ".ilp");
        exit(0);
    }
}

void IntegerLShaped::separate_continuous_cut(GRBLinExpr & cut){
    for(int i = 0; i < count_x; ++i){
        cut += (x[i] - x_[i])*duals[i];
    }
    cut += curr_2st_value;
    cut -= theta;
}

bool IntegerLShaped::add_continuous_cut(){
    if(curr_2st_value >= theta_ + 1e-6){
        GRBLinExpr cut = 0;
        separate_continuous_cut(cut);
        model->addConstr(cut <= 0, "Cont_opt_cut_" + std::to_string(nb_cont_opt_cuts)); ++nb_cont_opt_cuts;
        return true;
    }
    return false;
}
