#include "test1stsolver.hpp"

void Test1stSolver::verify_stat(){
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

Test1stSolver::Test1stSolver(const Input & input, const Instance & instance, int pr, std::string name): 
                            AbstractSolver(input,instance,pr,name) {
    
    obj_ = new double[instance.getNbScenarios(pr)];
    matching_students = new int*[instance.getNbScenarios(pr)];
    nb_assign = new int*[instance.getNbScenarios(pr)];
    unsatured = new bool[instance.getNbSchools()]; 

    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        matching_students[s] = new int[instance.getNbStudents()];
        nb_assign[s] = new int[instance.getNbSchools()];
    }

    curr_s = -1;
    count_y = 0;
    count_xi = 0;

    nb_per_rank.resize(instance.getNbScenarios(pr));
    rank_st.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        nb_per_rank[s].resize((instance.getNbMaxRankSchools()+1),0.0);
        rank_st[s].resize(instance.getNbStudents(),0.0);
    }
    avg_per_rank.resize((instance.getNbMaxRankSchools()+1),0.0);
}

Test1stSolver::~Test1stSolver(){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        delete[] matching_students[s];
        delete[] nb_assign[s];
    }
    delete[] obj_; 
    delete[] matching_students;
    delete[] nb_assign;
    delete[] unsatured;
    delete[] y_;  
    delete[] xi_;
}

void Test1stSolver::index_2st(int s){
    index_y[s].resize(instance.getNbFeasPairs(pr,s));
    index_xi[s].resize(instance.getNbFeasPairs(pr,s));

    // Same way it is computed in Solver
    for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i){
        index_y[s][i] = count_y;
        ++count_y;
    }

    // Same way it is computed in Transform
    for(int st = 0; st < instance.getNbStudents(); ++st){
        for(int i: instance.getFeasPairsByStudent(pr,s,st)){
            index_xi[s][i] = count_xi;
            ++count_xi;
        }
    }
}

void Test1stSolver::solve(const double *t_,Input::Transform type_transform){
    transf(t_,type_transform);
    solve(t_);
}

void Test1stSolver::reset(){
    lb = 0.0;
    ub = 0.0;
    std::fill(obj_,obj_+instance.getNbScenarios(pr),0.0);
    std::fill(unsatured,unsatured+instance.getNbSchools(),true);
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        std::fill(nb_assign[s],nb_assign[s]+instance.getNbSchools(),0);
        std::fill(matching_students[s],matching_students[s]+instance.getNbStudents(),-1);
    }
    std::fill(y_,y_+count_y,0.0);
    std::fill(xi_,xi_+count_xi,0.0);
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        std::fill(nb_per_rank[s].begin(),nb_per_rank[s].end(),0.0);
        std::fill(rank_st[s].begin(),rank_st[s].end(),0.0);
    }
    std::fill(avg_per_rank.begin(),avg_per_rank.end(),0.0);
}

void Test1stSolver::write_sol(){}

void Test1stSolver::write_inf(){}

//////////////////////////////////////////////////////////////////////////////////////

Test1stSolverModel::Test1stSolverModel(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env): 
                            Test1stSolver(input,instance,pr,name) {
    y = new GRBVar*[instance.getNbScenarios(pr)];
    z = new GRBVar*[instance.getNbScenarios(pr)]; 
    f = new GRBVar*[instance.getNbScenarios(pr)];

    constr_capa = new GRBConstr*[instance.getNbScenarios(pr)];
    constr_waste = new GRBConstr*[instance.getNbScenarios(pr)];
    constr_cutoff = new GRBConstr*[instance.getNbScenarios(pr)];

    xis_ = new double*[instance.getNbScenarios(pr)];
    rho_ = new double*[instance.getNbScenarios(pr)];

    model = new GRBModel*[instance.getNbScenarios(pr)];
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        model[s] = new GRBModel(*env);

        xis_[s] = new double[instance.getNbFeasPairs(pr,s)];
        rho_[s] = new double[instance.getNbFeasPairs(pr,s)];

        std::fill(xis_[s],xis_[s]+instance.getNbFeasPairs(pr,s),0.0);
        std::fill(rho_[s],rho_[s]+instance.getNbFeasPairs(pr,s),0.0);
    }

    count_zf = 0;

    params();
    create();
}

Test1stSolverModel::~Test1stSolverModel(){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        if(y) delete[] y[s];
        if(z) delete[] z[s];
        if(f) delete[] f[s];
        if(constr_capa) delete[] constr_capa[s];
        if(constr_waste) delete[] constr_waste[s];
        if(constr_cutoff) delete[] constr_cutoff[s];
        if(xis_) delete[] xis_[s];
        if(rho_) delete[] rho_[s];
        if(model)delete model[s];
    } 
    if(y) delete[] y;
    if(z) delete[] z;
    if(f) delete[] f;
    if(constr_capa) delete[] constr_capa;
    if(constr_waste) delete[] constr_waste;
    if(constr_cutoff) delete[] constr_cutoff;
    if(xis_) delete[] xis_;
    if(rho_) delete[] rho_;
    if(model) delete[] model;  
}

void Test1stSolverModel::params(){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        model[s]->set(GRB_IntParam_Threads, input.getNbThreads());
        model[s]->set(GRB_DoubleParam_TimeLimit, input.getTimeLimit());
        model[s]->set(GRB_DoubleParam_MIPGap, 1e-4);
        model[s]->set(GRB_DoubleParam_IntFeasTol, 1e-5);
        model[s]->set(GRB_DoubleParam_FeasibilityTol, 1e-6);
        model[s]->set(GRB_IntParam_OutputFlag, 0);
    }
}

void Test1stSolverModel::create(){
    index_y.resize(instance.getNbScenarios(pr));
    index_xi.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        index_2st(s);
        variables_2st(s);
        constraints_2st(s);
    }
    y_ = new double[count_y];
    xi_ = new double[count_xi];
}

void Test1stSolverModel::index_2st(int s){
    Test1stSolver::index_2st(s);

    for(int c = 0; c < instance.getNbSchools(); ++c){
        index_zf.insert({{s,c},count_zf});
        ++count_zf;
    }
}

void Test1stSolverModel::variables_2st(int s){
    y[s] = model[s]->addVars(instance.getNbFeasPairs(pr,s),GRB_BINARY);
    z[s] = model[s]->addVars(instance.getNbSchools(),GRB_CONTINUOUS);
    f[s] = model[s]->addVars(instance.getNbSchools(),GRB_BINARY);

    for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i) {
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int st = instance.getStudentFeasiblePair(pr,s,i);
        y[s][i].set(GRB_StringAttr_VarName,"y[" + std::to_string(s) + "," + std::to_string(st) + "," + std::to_string(c) + "]");
    }

    for(int c = 0; c < instance.getNbSchools(); ++c) 
        z[s][c].set(GRB_StringAttr_VarName,"z[" + std::to_string(s) + "," + std::to_string(c) + "]");

    for(int c = 0; c < instance.getNbSchools(); ++c) 
        f[s][c].set(GRB_StringAttr_VarName,"f[" + std::to_string(s) + "," + std::to_string(c) + "]");
}

void Test1stSolverModel::constraints_2st(int s){
    int nb_schools = instance.getNbSchools();
    int nb_students = instance.getNbStudents();

    constr_capa[s] = new GRBConstr[instance.getNbSchools()];
    constr_waste[s] = new GRBConstr[instance.getNbSchools()];
    constr_cutoff[s] = new GRBConstr[instance.getNbAssignFeasPairs(pr,s)];

    // Assignment constraints
    for(int st = 0; st < nb_students; ++st){
        GRBLinExpr exp = 0;
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) exp += y[s][i];
        exp += y[s][instance.getUnassignedPair(pr,s,st)];
        model[s]->addConstr(exp == 1);
    }

    for(int c = 0; c < nb_schools; ++c){
        GRBLinExpr exp = 0;
        for(int i : instance.getAssignFeasPairsBySchool(pr,s,c)) exp += y[s][i];    
        constr_capa[s][c] = model[s]->addConstr(exp <= 0.0);
    }

    // Cutoff constraints for stability
    for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
        int v = nb_students + 1;
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int r = instance.getScore(pr,s,i) + nb_students + 1;
        model[s]->addConstr(z[s][c] + v*y[s][i] <= r);
    }

    for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int st = instance.getStudentFeasiblePair(pr,s,i);

        int v = nb_students + 1;

        GRBLinExpr exp = -z[s][c];
        for(int l: instance.getStudentPreference(pr,s,st,i)) exp -= v*y[s][l];
        constr_cutoff[s][i] = model[s]->addConstr(exp <= 0.0);
    }
        
    // Waste constraints
    for(int c = 0; c < nb_schools; ++c){
        GRBLinExpr exp = f[s][c];
        for(int i : instance.getAssignFeasPairsBySchool(pr,s,c)) exp -= y[s][i];
        constr_waste[s][c] = model[s]->addConstr(exp <= 0.0);
    }

    for(int c = 0; c < nb_schools; ++c) model[s]->addConstr(z[s][c] -(nb_students+1)*f[s][c] <= 0.0);
}

void Test1stSolverModel::write_inf(){
    model[curr_s]->computeIIS();
    model[curr_s]->write(name + "_sce_" + std::to_string(curr_s) + ".ilp");
}

void Test1stSolverModel::solve(const double * t_){
    reset(); // reset information
    // Fixing capacity
    for(int c = 0; c < instance.getNbSchools(); ++c){
        double capacity = instance.getCapacity(c)+t_[c];
        for(int s = 0; s < instance.getNbScenarios(pr); ++s){
            constr_capa[s][c].set(GRB_DoubleAttr_RHS,capacity);
            model[s]->chgCoeff(constr_waste[s][c],f[s][c],capacity);
        }
    }

    double start_time = time(NULL); 
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        model[s]->update();
        model[s]->optimize();
        int gurobi_status = model[s]->get(GRB_IntAttr_Status);
        status = std::max(status,statusFromGurobi(gurobi_status));
        if(status == Status::Optimal || status == Status::Time_Limit) { 
            lb += model[s]->get(GRB_DoubleAttr_ObjBound)/instance.getNbScenarios(pr); 
            ub += model[s]->get(GRB_DoubleAttr_ObjVal)/instance.getNbScenarios(pr);
            gap = (ub-lb)/ub;
            obj_[s] = model[s]->get(GRB_DoubleAttr_ObjVal);
            get_variables(s);
        }
        else{ 
            curr_s = s;
            verify_stat();
        } 
    }
    if(verbose >= 1) print();
    time_ = time(NULL) - start_time;
}

void Test1stSolverModel::get_variables(int s){
    double *ys_ = model[s]->get(GRB_DoubleAttr_X,y[s],instance.getNbFeasPairs(pr,s));

    for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i){
        y_[index_y[s][i]] = ys_[i];
        xi_[index_xi[s][i]] = xis_[s][i];
        
        if(ys_[i] >= 0.999){
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            int st = instance.getStudentFeasiblePair(pr,s,i);
            
            matching_students[s][st] = i;
            if(c != instance.getNbSchools()) 
                ++nb_assign[s][c];
        }
    }
    
    for(int c = 0; c < instance.getNbSchools(); ++c) { 
        if(nb_assign[s][c] >= instance.getCapacity(c)) 
            unsatured[c] = false; 
    }

    delete[] ys_;

    for(int st = 0; st < instance.getNbStudents(); ++st){
        int r = 1;
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) {
            if(xi_[index_xi[s][i]] >= 0.9999){
                if(y_[index_y[s][i]] >= 0.9999){
                    ++nb_per_rank[s][r-1];
                    rank_st[s][st] = r;
                }
                ++r;
            } 
        }
        int i = instance.getUnassignedPair(pr,s,st); // xi unassigned always equal to one
        if(y_[index_y[s][i]] >= 0.9999){
            ++nb_per_rank[s][r-1];
            rank_st[s][st] = r;
        }
    }

    for(int r = 0; r <= instance.getNbMaxRankSchools(); ++r)
        avg_per_rank[r] += (1.0/instance.getNbScenarios(pr))*((double) nb_per_rank[s][r]);
}

void Test1stSolverModel::transf(const double *t_, Input::Transform transf){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        std::fill(xis_[s],xis_[s]+instance.getNbFeasPairs(pr,s),0.0);
        std::fill(rho_[s],rho_[s]+instance.getNbFeasPairs(pr,s),0.0);
    }
    if(transf == Input::Transform::NonStrategic) transf_nonstrategic();
    else transf_strategic(t_,transf);

    // Fixing endogenous uncertainty
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
            double score = -(instance.getScore(pr,s,i) + instance.getEpsilon())*xis_[s][i];
            constr_cutoff[s][i].set(GRB_DoubleAttr_RHS,score);
        }
        model[s]->set(GRB_DoubleAttr_UB,y[s],xis_[s],instance.getNbFeasPairs(pr,s));
        model[s]->set(GRB_DoubleAttr_Obj,y[s],rho_[s],instance.getNbFeasPairs(pr,s));
    }
}

void Test1stSolverModel::transf_strategic(const double *t_,Input::Transform transf){
    std::map<int,int> ranks;
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            if(transf == Input::Transform::Compact || transf == Input::Transform::Enumeration)
                instance.mi_algorithm(pr,s,st,t_,ranks);
            else if(transf == Input::Transform::Simple)
                instance.order_algorithm(pr,s,st,t_,ranks);
            int r = 1;
            if(ranks.empty() == false){
                for(auto it = ranks.begin(); it != ranks.end(); it++){
                    xis_[s][it->second] = 1;
                    rho_[s][it->second] = (double)r;
                    ++r;
                }
            }
            int i = instance.getUnassignedPair(pr,s,st);
            xis_[s][i] = 1;
            if(input.getObjAccessImprov() == false) rho_[s][i] = (double)(instance.getNbSchools()+1);
            else rho_[s][i] = (double)r;
            ranks.clear();
        }
    } 
}

void Test1stSolverModel::transf_nonstrategic(){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            int num = std::min(instance.getNbAssignFeasPairsByStudent(pr,s,st),
                               instance.getNbMaxRankSchools());
            for(int r = 0; r < num; ++r){
                int i = instance.getAssignFeasPairByStudentPos(pr,s,st,r);
                xis_[s][i] = 1;
                rho_[s][i] = (double)(r+1);
            }
            int i = instance.getUnassignedPair(pr,s,st);
            xis_[s][i] = 1;
            rho_[s][i] = (double)(num+1);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////

Test1stSolverDA::Test1stSolverDA(const Input & input, const Instance & instance, int pr, std::string name): 
                            Test1stSolver(input,instance,pr,name) {

    matching_schools.resize(instance.getNbScenarios(pr));
    feas_pairs.resize(instance.getNbScenarios(pr));
    rank_feas_pairs.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        matching_schools[s].resize(instance.getNbSchools()+1);
        feas_pairs[s].resize(instance.getNbStudents());
    }

    index_y.resize(instance.getNbScenarios(pr));
    index_xi.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) index_2st(s);
    y_ = new double[count_y];
    xi_ = new double[count_xi];
}

void Test1stSolverDA::solve(const double * t_){
    reset();
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        for(int c = 0; c < instance.getNbSchools()+1; ++c){
            matching_schools[s][c].clear();
        }
    }
    double start_time = time(NULL); 
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        double val = instance.da_algorithm(pr,s,t_,feas_pairs[s],rank_feas_pairs[s],matching_schools[s]);
        lb += val/instance.getNbScenarios(pr);
        ub += val/instance.getNbScenarios(pr);
        obj_[s] = val;
        get_variables(s);
    }
    status = Status::Optimal;
    if(verbose >= 1) print();
    time_ = time(NULL) - start_time;
}

void Test1stSolverDA::get_variables(int s){
    for(int c = 0; c < instance.getNbSchools()+1; ++c){
        for(int i: matching_schools[s][c]) {
            int st = instance.getStudentFeasiblePair(pr,s,i);
            y_[index_y[s][i]] = 1.0;
            matching_students[s][st] = i;
        }
        if(c != instance.getNbSchools()) nb_assign[s][c] = matching_schools[s][c].size();
    }

    for(int c = 0; c < instance.getNbSchools(); ++c) { 
        if(nb_assign[s][c] >= instance.getCapacity(c)) 
            unsatured[c] = false; 
    }

    for(int st = 0; st < instance.getNbStudents(); ++st){
        if(matching_students[s][st] == -1) std::cout << "Unassigned student: s:" << s << " st:" << st << std::endl;
    }

    for(int st = 0; st < instance.getNbStudents(); ++st){
        int r = 1;
        for(int i: feas_pairs[s][st]) {
            xi_[index_xi[s][i]] = 1.0;
            if(y_[index_y[s][i]] >= 0.999){
                ++nb_per_rank[s][r-1];
                rank_st[s][st] = r;
            }
            ++r;
        }
        int i = instance.getUnassignedPair(pr,s,st);
        xi_[index_xi[s][i]] = 1.0;
        if(y_[index_y[s][i]] >= 0.999){ 
            ++nb_per_rank[s][r-1];
            rank_st[s][st] = r;
        }
    }
    for(int r = 0; r <= instance.getNbMaxRankSchools(); ++r)
        avg_per_rank[r] += (1.0/instance.getNbScenarios(pr))*((double) nb_per_rank[s][r]);
}

void Test1stSolverDA::transf(const double *t_, Input::Transform transf){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s) {
        for(int st = 0; st < instance.getNbStudents(); ++st){
            feas_pairs[s][st].clear();
        }
        rank_feas_pairs[s].clear();
    }
    if(transf == Input::Transform::NonStrategic) transf_nonstrategic();
    else transf_strategic(t_,transf);
}

void Test1stSolverDA::transf_strategic(const double *t_,Input::Transform transf){
    std::map<int,int> ranks;
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            if(transf == Input::Transform::Compact || transf == Input::Transform::Enumeration)
                instance.mi_algorithm(pr,s,st,t_,ranks);
            else if(transf == Input::Transform::Simple)
                instance.order_algorithm(pr,s,st,t_,ranks);
            int r = 1;
            if(ranks.empty() == false){
                for(auto it = ranks.begin(); it != ranks.end(); it++){
                    feas_pairs[s][st].push_back(it->second);
                    rank_feas_pairs[s][it->second] = (double)r;
                    ++r;
                }
            }
            int i = instance.getUnassignedPair(pr,s,st);
            if(input.getObjAccessImprov() == false) rank_feas_pairs[s][i] = (double)(instance.getNbSchools()+1);
            else rank_feas_pairs[s][i] = (double)r;
            ranks.clear();
        }
    } 
}

void Test1stSolverDA::transf_nonstrategic(){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            int num = std::min(instance.getNbAssignFeasPairsByStudent(pr,s,st),
                               instance.getNbMaxRankSchools());
            for(int r = 0; r < num; ++r){
                int i = instance.getAssignFeasPairByStudentPos(pr,s,st,r);
                feas_pairs[s][st].push_back(i);
                rank_feas_pairs[s][i] = (double)(r+1);
            }
            int i = instance.getUnassignedPair(pr,s,st);
            rank_feas_pairs[s][i] = (double)(num+1);
        }
    }
}
