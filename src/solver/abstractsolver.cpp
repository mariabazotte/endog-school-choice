#include "abstractsolver.hpp"
#include "../transform/abstracttransform.hpp"
#include "test1st/test1stsolver.hpp"
#include "asgheuristic/asgheuristicsolver.hpp"
#include "localsearch/localsearch.hpp"

void AbstractSolver::verify_stat(){
    if(status == Status::Unbounded){
        std::cerr << "The " << name << " is unbounded, this cannot happen according to the problem definition." << std::endl;
        throw std::runtime_error( name + " is unbounded.\n");
        exit(0);
    }
    if(status == Status::Infeasible_or_Unbounded){
        std::cerr << "The " << name << " is unbounded or infeasible, this cannot happen according to the problem definition." << std::endl;
        write_inf();
        throw std::runtime_error(name + " is unbounded or infeasible.\n");
        exit(0); 
    } 
    if(status == Status::Infeasible){
        std::cerr << "The " << name << " is infeasible, this cannot happen according to the problem definition." << std::endl;
        write_inf();
        throw std::runtime_error(name + " is infeasible.\n");
        exit(0); 
    }
}

void AbstractSolver::print(){
    std::string output;
    if(status == Status::Optimal){
        output = "OPTIMAL: The optimal solution of the " + name + " was " + std::to_string(ub) + ".\n";
    }else if(status == Status::Time_Limit){
        output = "TIME: The best bound found for " + name  + " was " + std::to_string(lb) + ".\n";
        output = "TIME: The best feasible solution found for " + name + " was " + std::to_string(ub) + ".\n";
    }
    else if(status == Status::Iteration_Limit){
        output = "ITERATION: The best bound found for " + name  + " was " + std::to_string(lb) + ".\n";
        output = "ITERATION: The best feasible solution found for " + name + " was " + std::to_string(ub) + ".\n";
    }
    std::cout << output; 
    if(verbose == 2) write_sol(); 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Solver::Solver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env,bool warmstart, bool linear_relax) : 
                            AbstractSolver(input,instance,pr,name),warmstart(warmstart),linear_relax(linear_relax),combs(input.useCombs()),
                            model(new GRBModel(*env)) {
    ub_callback = std::numeric_limits<double>::infinity();
    time_best_ub = 0.0;

    include_frac_cuts = true;
    y_continuous = false;
    
    // Create transformation
    TransformFactory factory;
    transf = factory.createTransform(input,instance,pr,this);

    // Initialize auxiliary heuristics to find initial feasible solutions on warm start
    if(warmstart == true && linear_relax == false && input.getTransform() == Input::Transform::NonStrategic && input.getSolver() != Input::Solver::ASGHeuristic) 
        asg = new ASGheuristicSolver(input,instance,pr,"Warm Start -- ASG -- Problem " + std::to_string(pr),env,false,false);
    if(warmstart == true && linear_relax == false && input.getTransform() != Input::Transform::NonStrategic && input.getSolver() != Input::Solver::LocalSearch)
        ls = new LocalSearchSolver(input,instance,pr,"Warm Start -- Local Search -- Problem " + std::to_string(pr),env,false,true);
    // Initialize test 1st to test the value of first-stage solutions on warm start (using DA algorithm)
    if(warmstart == true) 
        test_1st = new Test1stSolverDA(input,instance,pr,"Test Capacity -- Problem " + std::to_string(pr));
    ub_ws = GRB_INFINITY;

    // Initialize count for variables
    count_tt = 0;
    count_y = 0;
    count_x = 0;
}

Solver::~Solver(){
    if(t_) delete[] t_;
    if(y_) delete[] y_;
    if(x_) delete[] x_;
    if(xa_) delete[] xa_;
    if(tooth_minus_) delete[] tooth_minus_;

    if(t) delete[] t;
    if(y) delete[] y;
    if(x) delete[] x;
    if(xa) delete[] xa;
    if(vtooth) delete[] vtooth;
    delete model;    

    delete transf;

    if(asg) delete asg;
    if(test_1st) delete test_1st;
}

void Solver::params(){
    model->set(GRB_IntParam_Threads, input.getNbThreads());
    model->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
    model->set(GRB_DoubleParam_MIPGap, 1e-4);
    model->set(GRB_DoubleParam_IntFeasTol, 1e-5);
    model->set(GRB_DoubleParam_FeasibilityTol, 1e-6);
    if(verbose < 1) model->set(GRB_IntParam_OutputFlag, 0);
    else model->set(GRB_IntParam_OutputFlag, 1);
    if(input.useLazyStability() == true || combs == true) model->set(GRB_IntParam_LazyConstraints,1);
    if(combs == true ) model->setCallback(this);
}

void Solver::write_inf(){
    model->computeIIS();
    model->write(name + ".ilp");
    model->write(name + ".lp");
}

void Solver::write_sol(){
    model->write(name + ".lp");
    model->write(name + ".sol");
}

////////////////////////////////////////////////////////////////////////////////////////

void Solver::create_1st(){
    index_1st();
    variables_1st();
    constraints_1st();  
}

void Solver::index_1st(){
    if((combs == true) || (input.getTransform() != Input::Transform::NonStrategic)){
        index_x.resize(instance.getNbSchools());
        for(int c = 0; c < instance.getNbSchools(); ++c){
            index_x[c].resize(instance.getNbMaxExtraSeats());
            for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
                index_x[c][j] = count_x;
                ++count_x;
            }
        }
    }
    count_tt = instance.getNbSchools();
    t_ = new double[count_tt];
}

void Solver::variables_1st(){
    if(linear_relax == false) t = model->addVars(count_tt,GRB_INTEGER);
    else t = model->addVars(count_tt,GRB_CONTINUOUS);

    for(int c = 0; c < count_tt; ++c) {
        t[c].set(GRB_DoubleAttr_UB, instance.getNbMaxExtraSeats());
        t[c].set(GRB_StringAttr_VarName, "t[" + std::to_string(c) + "]");
    }

    if((combs == true) || (input.getTransform() != Input::Transform::NonStrategic)){
        x  = model->addVars(count_x,GRB_BINARY);
        
        for(int i = 0; i < count_x; ++i) {
            if(linear_relax == true) x[i].set(GRB_CharAttr_VType,GRB_CONTINUOUS);
            x[i].set(GRB_DoubleAttr_LB,0.0);
            x[i].set(GRB_DoubleAttr_UB,1.0);
        }

        for(int c = 0; c < instance.getNbSchools(); ++c){
            for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
                x[index_x[c][j]].set(GRB_StringAttr_VarName,
                        "x[" + std::to_string(c) + 
                         "," + std::to_string(j) + "]");
            }
        }

        // The t variables can be defined as continuous (sum of binary variables)
        for(int c = 0; c < count_tt; ++c) t[c].set(GRB_CharAttr_VType,GRB_CONTINUOUS);
    }

    if(combs == true){
        xa = model->addVars(count_x,GRB_CONTINUOUS);

        for(int i = 0; i < count_x; ++i) {
            xa[i].set(GRB_DoubleAttr_LB,0.0);
            xa[i].set(GRB_DoubleAttr_UB,1.0);
        }

        for(int c = 0; c < instance.getNbSchools(); ++c){
            for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
                xa[index_x[c][j]].set(GRB_StringAttr_VarName,
                        "xa[" + std::to_string(c) + 
                          "," + std::to_string(j) + "]");
            }
        }
    }

    // Trying to find better bounds for t
    for(int c = 0; c < instance.getNbSchools(); ++c){
        int max_value = instance.getNbAssignFeasPairsBySchool(pr,0,c);
        for(int s = 0; s < instance.getNbScenarios(pr); ++s){
            if(instance.getNbAssignFeasPairsBySchool(pr,s,c) > max_value)
                max_value = instance.getNbAssignFeasPairsBySchool(pr,s,c);
        }
        if(max_value < instance.getNbMaxExtraSeats()){
            t[c].set(GRB_DoubleAttr_UB,max_value);
        }
    } 
}

void Solver::constraints_1st(){
    int nb_schools = instance.getNbSchools();
    int nb_seats = instance.getNbMaxExtraSeats();

    GRBLinExpr budget = 0;
    for(int c = 0; c < nb_schools; ++c) budget += getT(c); 
    model->addConstr(budget <= nb_seats);

    if((combs == true) || (input.getTransform() != Input::Transform::NonStrategic)){
        // Linking x with t
        for(int c = 0; c < instance.getNbSchools(); ++c){
            GRBLinExpr exp = getT(c);
            for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j) exp -= getX(c,j);
            model->addConstr(exp == 0.0);
        }

        // Defining variables x
        for(int c = 0; c < instance.getNbSchools(); ++c){
            for(int j = 0; j < instance.getNbMaxExtraSeats() - 1; ++j){
                model->addConstr(getX(c,j) >= getX(c,j+1));
            }
        }
    }

    if(combs == true){
        for(int c = 0; c < instance.getNbSchools(); ++c){
            for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j){
                model->addConstr(getX(c,j) + getXA(c,j) == 1.0);
            }
        }
    }
}

void Solver::create_2st(){
    index_y.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) index_2st(s);
    variables_2st();
    transf->transf_2st();
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) constraints_2st(s);
    y_ = new double[count_y];
    if(combs == true) tooth_minus_ = new double[count_y];
}

void Solver::index_2st(int s){
    index_y[s].resize(instance.getNbFeasPairs(pr,s));
    for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i){
        index_y[s][i] = count_y;
        ++count_y;
    }
}

void Solver::variables_2st(){
    if(linear_relax == false) y = model->addVars(count_y,GRB_BINARY);
    else y = model->addVars(count_y,GRB_CONTINUOUS);
    
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i){
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            int st = instance.getStudentFeasiblePair(pr,s,i);
            y[index_y[s][i]].set(GRB_StringAttr_VarName,
                        "y[" + std::to_string(s) + "," + std::to_string(st) + "," + std::to_string(c) + "]");
            if(linear_relax == true) y[index_y[s][i]].set(GRB_DoubleAttr_LB,0.0);
            if(linear_relax == true) y[index_y[s][i]].set(GRB_DoubleAttr_UB,1.0);
        }
    }

    // The y variables can be defined as continuous when combs == true
    // Only used if y_continuous == true
    // (only implemented for non-strategic version)
    if(combs == true && y_continuous == true){
        for(int i = 0; i < count_y; ++i){
            y[i].set(GRB_CharAttr_VType,GRB_CONTINUOUS);
            y[i].set(GRB_DoubleAttr_LB,0.0);
            y[i].set(GRB_DoubleAttr_UB,1.0);
        }
    }

    // Priority to branching on tooth: define new auxiliary variables for the tooth
    if(input.useBranchTooth() == true){
        vtooth = model->addVars(count_y,GRB_CONTINUOUS); 

        for(int s = 0; s < instance.getNbScenarios(pr); ++s){
            for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i){
                int c = instance.getSchoolFeasiblePair(pr,s,i);
                int st = instance.getStudentFeasiblePair(pr,s,i);
                vtooth[index_y[s][i]].set(GRB_IntAttr_BranchPriority,1); // branch on tooth 
                vtooth[index_y[s][i]].set(GRB_StringAttr_VarName,
                                "tooth[" + std::to_string(s) + "," + std::to_string(st) + "," + std::to_string(c) + "]");
            }
        } 
    }
}

void Solver::constraints_2st(int s){
    int nb_schools = instance.getNbSchools();
    int nb_students = instance.getNbStudents();

    // Assignment constraints
    for(int st = 0; st < nb_students; ++st){
        GRBLinExpr exp = 0;
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) exp += getY(s,i);
        exp += getY(s,instance.getUnassignedPair(pr,s,st));
        model->addConstr(exp == 1);
    }

    // Capacity constraints
    for(int c = 0; c < nb_schools; ++c){
        GRBLinExpr exp = -getT(c);
        for(int i : instance.getAssignFeasPairsBySchool(pr,s,c)) exp += getY(s,i);    
        model->addConstr(exp <= instance.getCapacity(c));
    }

    // Defining tooth to branch
    if(input.useBranchTooth() == true){
        for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
            int st = instance.getStudentFeasiblePair(pr,s,i);

            GRBLinExpr exp = getVTooth(s,i);
            for(int ii: instance.getStudentPreference(pr,s,st,i)) exp -= getY(s,ii);
            model->addConstr(exp == 0.0);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////

void Solver::solve(){
    double start_time = time(NULL); 
    if(warmstart == true) warm_start();

    model->update();
    model->optimize();

    status = statusFromGurobi(model->get(GRB_IntAttr_Status));
    if(status == Status::Optimal || status == Status::Time_Limit) upd_solution();
    else verify_stat(); 

    time_ = time(NULL) - start_time;
    if(verbose >= 1) print();
}

void Solver::get_variables(){
    if(nb_sol > 0){
        if(ub <= ub_ws || warmstart == false){
            if(t_) delete[] t_;
            t_ = model->get(GRB_DoubleAttr_X,t,count_tt);
            for(int c = 0; c < count_tt; ++c) t_[c] = std::round(t_[c]);
        }

        if(y_) delete[] y_;
        y_ = model->get(GRB_DoubleAttr_X,y,count_y);
        transf->get_variables();
    }
    if(nb_sol == 0 && t_ != NULL) nb_sol = 1;
}

void Solver::upd_solution(){
    if(linear_relax == true) {
        ub = model->get(GRB_DoubleAttr_ObjVal);
    }else{
        lb = std::max(0.0,model->get(GRB_DoubleAttr_ObjBound)); // gap between 0 and 1 (0.0 is a valid lower bound)
        nb_bnb = model->get(GRB_DoubleAttr_NodeCount);
        nb_sol = model->get(GRB_IntAttr_SolCount);
        if((model->get(GRB_DoubleAttr_ObjVal) <= ub_ws && nb_sol > 0) || warmstart == false)
            ub = model->get(GRB_DoubleAttr_ObjVal);
    }
    gap = (ub-lb)/ub;
    get_variables();
}

void Solver::warm_start(){
    model->set(GRB_IntAttr_NumStart,2); 

    // First warm start initial solution
    model->set(GRB_IntParam_StartNumber,0);
    warm_start_capazero();

    // Second warm start initial solution
    model->set(GRB_IntParam_StartNumber,1);
    if(input.getTransform() == Input::Transform::NonStrategic && linear_relax == false) warm_start_asg();
    if(input.getTransform() != Input::Transform::NonStrategic && linear_relax == false) warm_start_ls();
}

void Solver::warm_start_asg(){
    if(verbose >= 1) std::cout << "\n### Starting warm-start ASG ###\n\n";
    
    asg->solve(); 
    if(asg->getNbSol() > 0){ 
        model->set(GRB_DoubleAttr_Start,t,asg->getT_(),count_tt);
        model->set(GRB_DoubleAttr_Start,y,asg->getY_(),count_y);
        transf->set_inisol_warm_start(asg->getTransform());
        ub = asg->getUB();
        ub_ws = ub;
        std::copy(asg->getT_(),asg->getT_()+count_tt,t_);
    }

    if(verbose >= 1) std::cout << "\n#### Ending warm-start ASG ###\n\n";
}

void Solver::warm_start_ls(){
    if(verbose >= 1) std::cout << "\n### Starting warm-start Local Search ###\n\n";

    ls->solve();
    model->set(GRB_DoubleAttr_Start,t,ls->getT_(),count_tt);
    model->set(GRB_DoubleAttr_Start,y,ls->getTest1st()->getY_(),count_y);  
    transf->set_inisol_warm_start(ls->getTest1st());
    ub = ls->getUB();
    ub_ws = ub;
    std::copy(ls->getT_(),ls->getT_()+count_tt,t_);

    if(verbose >= 1) std::cout << "\n#### Ending warm-start Local Search ###\n\n";
}

void Solver::warm_start_capazero(){
    if(verbose >= 1) std::cout << "\n### Starting warm-start zero-capacity ###\n\n";
    
    double *zero_t_ = new double[count_tt];
    std::fill(zero_t_,zero_t_+count_tt,0.0);

    test_1st->solve(zero_t_,input.getTransform());

    // Setting initial solution
    model->set(GRB_DoubleAttr_Start,t,zero_t_,count_tt);
    model->set(GRB_DoubleAttr_Start,y,test_1st->getY_(),count_y);
    transf->set_inisol_warm_start(test_1st);

    // Setting auxiliar constraints
    transf->set_constr_warm_start(test_1st);

    delete[] zero_t_;
    if(verbose >= 1) std::cout << "\n#### Ending warm-start zero-capacity ###\n\n";
}

void Solver::callback(){
    try{
        // If combs == true and y_continuous == true, 
        // only x variables are binary and the y variables are continuous,
        // thus MIPSOL finds a solution with binary x and possibly fractional y
        if(where == GRB_CB_MIPSOL && combs == true){

            if(t_) delete[] t_;
            if(y_) delete[] y_;
            if(x_) delete[] x_;
            if(xa_) delete[] xa_;

            t_ = getSolution(t,count_tt);
            y_ = getSolution(y,count_y);

            x_ = getSolution(x,count_x);
            xa_ = getSolution(xa,count_x);
            for(int s = 0; s < instance.getNbScenarios(pr); ++s)
                tooth_minus(s);

            // Avoid numerical issues
            for(int c = 0; c < count_tt; ++c)  
                t_[c] = std::round(t_[c]);

            // Solve DA algorithm to obtain actual (non-fractional) matching
            test_1st->solve(t_,input.getTransform()); 
            
            // If a better feasible solution was found, update best upper bound
            if(test_1st->getUB() < ub_callback){ 
                ub_callback = test_1st->getUB();

                setSolution(t,t_,count_tt);
                setSolution(y,test_1st->getY_(),count_y);

                if(x) setSolution(x,x_,count_x);
                if(xa) setSolution(xa,xa_,count_x);

                // TODO: set value auxiliary variables if strategic version

                useSolution();
            }

            // If current objective is worst than best upper bound, then we don't add cuts
            if(ub_callback > getDoubleInfo(GRB_CB_MIPSOL_OBJ) - 10e-4){
                // Add combs
                for(int c = 0; c < instance.getNbSchools(); ++c){
                    int capacity = instance.getCapacity(c) + instance.getNbMaxExtraSeats();

                    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
                        bool found = false;
                        GRBLinExpr exp_comb = 0;
                        comb_baiou(c,s,found,exp_comb,x_);
                        if(found == true) addLazy(exp_comb >= capacity);
                    }
                } 
            }
        }
        // If combs == true and y_continuous == true, 
        // MIPNODE always find fractional capacities (x is fractional)
        // If y_continuous == false, this is not true
        if(where == GRB_CB_MIPNODE && combs == true && include_frac_cuts == true) {
            if(getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL){

                if(y_) delete[] y_;
                if(x_) delete[] x_;
                if(xa_) delete[] xa_;
                
                y_ = getNodeRel(y,count_y);
                x_ = getNodeRel(x,count_x);
                xa_ = getNodeRel(xa,count_x);
                for(int s = 0; s < instance.getNbScenarios(pr); ++s)
                    tooth_minus(s);

                // Add combs
                for(int c = 0; c < instance.getNbSchools(); ++c){
                    int capacity = instance.getCapacity(c) + instance.getNbMaxExtraSeats();

                    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
                        bool found = false;
                        GRBLinExpr exp_comb = 0;
                        comb_baiou(c,s,found,exp_comb,x_);
                        if(found == true) addLazy(exp_comb >= capacity);
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

///////////////////////////////////////////////////////////////////////////////////

// Auxiliary function to compute tooth minus
void Solver::tooth_minus(int s){
    for(int st = 0; st < instance.getNbStudents(); ++st){
        double cumulated = 0.0;
        
        // assign pairs
        int nb_assignpairs = instance.getNbAssignFeasPairsByStudent(pr,s,st);
        if(nb_assignpairs > 0){
            for(int n = 0; n < nb_assignpairs; ++n){
                int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
                tooth_minus_[index_y[s][i]] = cumulated;
                cumulated += getY_(s,i);
            }
        }

        // unassign pair
        int i = instance.getUnassignedPair(pr,s,st);
        tooth_minus_[index_y[s][i]] = cumulated;
    }
}

//////////////////////////////////////////////////////////////////////////////////

// Function to separate combs 
void Solver::comb_baiou(int &c, int &s, bool &found, GRBLinExpr &exp_comb, double * used_x_){
    found = false;
    exp_comb = 0;

    if(instance.getNbAssignFeasPairsBySchool(pr,s,c) < instance.getCapacity(c)) return;

    int capacity = instance.getCapacity(c);
    int nb_seats = instance.getNbMaxExtraSeats();

    // Function to order feas_pairs + artificial_pairs according to tooth-minus
    auto cmp = [&c,&s,&used_x_,this](const int &a, const int &b){ if(a >= 0 && b >= 0) return tooth_minus_[index_y[s][a]] < tooth_minus_[index_y[s][b]];
                                                         else if(a >= 0) return tooth_minus_[index_y[s][a]] < used_x_[index_x[c][-(b)-1]];
                                                         else if(b >= 0) return used_x_[index_x[c][-(a)-1]] < tooth_minus_[index_y[s][b]];
                                                         else return used_x_[index_x[c][-(a)-1]] < used_x_[index_x[c][-(b)-1]]; };
    
    std::vector<int> feas_pairs(instance.getAssignFeasPairsBySchool(pr,s,c)); // Defining vector with feas pairs student c 
    for(int j = 1; j <= nb_seats; ++j) feas_pairs.push_back(-j); // Including artificial pairs
    std::sort(feas_pairs.begin(),feas_pairs.end(),cmp); // Sort according to tooth-minus

    std::set<int> feas_pairs_Q;  // Defining sets Q and U
    std::list<int> feas_pairs_U;
    for(int n = 0; n < (int)feas_pairs.size(); ++n){
        if(n < (capacity+nb_seats-1)) feas_pairs_Q.insert(feas_pairs[n]); 
        else feas_pairs_U.push_back(feas_pairs[n]);
    }

    int i = instance.getLeastPrefFeasPairBySchool(pr,s,c);
    int pref_school_i = instance.getFeasPairSchoolRank(pr,s,i); // start with 0

    while(true){
        if(pref_school_i < (capacity - 1)) break;

        if(feas_pairs_Q.find(i) != feas_pairs_Q.end()){
            int i_bar =  feas_pairs_U.front(); // Feas pair in U with minimum val of tooth minus
            feas_pairs_U.pop_front();
        
            feas_pairs_Q.erase(i);
            feas_pairs_Q.insert(i_bar); 
        }
        else{
            feas_pairs_U.remove(i);
        }

        /* Computing possible violated comb. */
        double value_comb = 0.0;
        std::set<int> comb;

        // Shaft of feas pair (i)
        for(int j = 0; j < nb_seats; ++j) {     // artificial pairs - always in shaft (in every comb)
            value_comb += (1.0 - used_x_[index_x[c][j]]);
        }     
        for(int n = 0; n <= pref_school_i; ++n) { // feas pairs
            int ii = instance.getAssignFeasPairBySchoolPos(pr,s,c,n);
            comb.insert(ii);
            value_comb += getY_(s,ii);
        } 

        // Tooth minus of feas pair (i)  - It cannot be artificial feas pairs (pref aux feas pairs < capacity + nb_seats -1 )
        int st_i = instance.getStudentFeasiblePair(pr,s,i);
        int pref_student_i = instance.getFeasPairStudentRank(pr,s,i);
        if(pref_student_i > 0) { 
            for(int n = 0; n < pref_student_i; ++n){
                int ii = instance.getAssignFeasPairByStudentPos(pr,s,st_i,n);
                comb.insert(ii);
                value_comb += getY_(s,ii); 
            }
        }

        // Tooth of every feasible pair in Q - They can be the artificial feas pairs 
        for(auto it = feas_pairs_Q.begin(); it != feas_pairs_Q.end(); it++){
            if(*it < 0) { // artificial pair
                comb.insert(*it);
                int j = -(*it)-1;
                value_comb += used_x_[index_x[c][j]];
            }
            else{
                int st_it = instance.getStudentFeasiblePair(pr,s,*it);
                int pref_student_it = instance.getFeasPairStudentRank(pr,s,*it); 
                for(int n = 0; n <= pref_student_it; ++n) {
                    int ii = instance.getAssignFeasPairByStudentPos(pr,s,st_it,n);
                    if(comb.find(ii) == comb.end()){
                        comb.insert(ii); 
                        value_comb += getY_(s,ii);
                    }
                }
            }
        }

        if(value_comb <= (capacity + nb_seats - 0.00001)) {
            found = true; 
            exp_comb = 0;
            generate_comb_expr(c,s,exp_comb,comb);
            break;
        }
        else{ 
            --pref_school_i;
            i = instance.getAssignFeasPairBySchoolPos(pr,s,c,pref_school_i);
        }
    }
}

void Solver::generate_comb_expr(int &c, int &s, GRBLinExpr &exp_comb, std::set<int> & comb){
    for(int j = 0; j < instance.getNbMaxExtraSeats(); ++j) exp_comb += getXA(c,j); // in every comb
    for(int i: comb) { if(i < 0) { exp_comb += getX(c,-(i)-1);} else { exp_comb += getY(s,i); }}
}