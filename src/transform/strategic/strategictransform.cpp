#include "strategictransform.hpp"
#include "../../solver/abstractsolver.hpp"
#include "../../solver/test1st/test1stsolver.hpp"

void StrategicTransform::print_xi(){
    if(xi_) delete[] xi_;
    xi_ = solver->getModel()->get(GRB_DoubleAttr_X,xi,count_xi);
    for(int s = 0; s < instance.getNbScenarios(pr);++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            std::cout << "s:" << s << "\tst:" << st << std::endl;
            for(int i: instance.getFeasPairsByStudent(pr,s,st)) std::cout << xi_[index_xi[s][i]] << " ";
            std::cout << std::endl;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////

// Functions to simplify the definition of variables x
GRBLinExpr StrategicTransform::getX(int c, int j) { 
    GRBLinExpr exp = 0.0;
    if(j == 0) exp = (1.0 - solver->getX(c,j));
    else if(j == instance.getNbMaxExtraSeats()) exp = solver->getX(c,j-1);
    else exp = (solver->getX(c,j-1) - solver->getX(c,j)); 
    return exp;
}

double StrategicTransform::getX_(int c, int j) { 
    if(j == 0)  return (1.0 - solver->getX(c,0).get(GRB_DoubleAttr_X));
    else if(j == instance.getNbMaxExtraSeats()) return solver->getX(c,j-1).get(GRB_DoubleAttr_X);
    else return (solver->getX(c,j-1).get(GRB_DoubleAttr_X) - solver->getX(c,j).get(GRB_DoubleAttr_X)); 
}

/////////////////////////////////////////////////////////////////////////////////////////////

// Functions to define transformation part of the model
void StrategicTransform::params(){
    // Avoid problems with numerical instability
    // solver->getModel()->set(GRB_DoubleParam_FeasibilityTol, 1e-9); 
    if(input.useLargeNumStab() == true) 
        solver->getModel()->set(GRB_IntParam_NumericFocus, 3);
    // else solver->getModel()->set(GRB_IntParam_NumericFocus, 1);
    
    // Scaleflag for parameters 
    solver->getModel()->set(GRB_IntParam_ScaleFlag, 2);

    // Disable aggregate to avoid numerical issues
    solver->getModel()->set(GRB_IntParam_Aggregate, 0);

    // Disable preprocessing
    if(input.restPrep() == true){
        solver->getModel()->set(GRB_IntParam_Presolve, 0);
        // solver->getModel()->set(GRB_IntParam_PrePasses, 1)
    }

    // Use barrier method for linear relaxation
    if(input.useBarrier() == true) {
        solver->getModel()->set(GRB_IntParam_Method, 2); 
        solver->getModel()->set(GRB_IntParam_Crossover, 0);
    }
}

void StrategicTransform::transf_2st(){
    params();
    
    index_xi.resize(instance.getNbScenarios(pr));
    index_phi.resize(instance.getNbScenarios(pr));
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){ 
        index_xi[s].resize(instance.getNbFeasPairs(pr,s));
        index_phi[s].resize(instance.getNbStudents());
        for(int st = 0; st < instance.getNbStudents(); ++st) 
            index_2st(s,st);
    }

    variables_2st();

    for(int s = 0; s < instance.getNbScenarios(pr); ++s){ 
        for(int st = 0; st < instance.getNbStudents(); ++st) {
            int nb_ranks = instance.getNbMaxRankSchools();
            int nb_pairs = instance.getNbAssignFeasPairsByStudent(pr,s,st);
            if(nb_pairs > 0 && nb_pairs <= nb_ranks) constraints_preprocess_2st(s,st);
            else if(nb_pairs > 0) constraints_2st(s,st);
            else {
                // Unassigned pair always selected
                int i = instance.getUnassignedPair(pr,s,st);
                solver->getModel()->addConstr(getXi(s,i) == 1.0);
                // getXi(s,i).set(GRB_DoubleAttr_LB,1.0);
                // getXi(s,i).set(GRB_DoubleAttr_UB,1.0);
            }
        }      
    }
    objective_2st();

    solver->getModel()->update();
}

void StrategicTransform::index_2st(int s,int st){
    for(int i: instance.getFeasPairsByStudent(pr,s,st)){
        index_xi[s][i] = count_xi;
        ++count_xi;
    }

    index_phi[s][st] = count_phi;
    ++count_phi;
}

void StrategicTransform::variables_2st(){
    if(solver->getLinearRelax() == false) xi = solver->getModel()->addVars(count_xi,GRB_BINARY);
    else xi = solver->getModel()->addVars(count_xi,GRB_CONTINUOUS);
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int i = 0; i < instance.getNbFeasPairs(pr,s); ++i){
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            int st = instance.getStudentFeasiblePair(pr,s,i);
            xi[index_xi[s][i]].set(GRB_StringAttr_VarName,
                        "xi[" + std::to_string(s) + "," + std::to_string(st) + "," + std::to_string(c) + "]");
            if(solver->getLinearRelax() == true) solver->getModel()->addConstr(xi[index_xi[s][i]] <= 1.0);
        }
    }

    phi = solver->getModel()->addVars(count_phi,GRB_CONTINUOUS);
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            phi[index_phi[s][st]].set(GRB_StringAttr_VarName,
                    "phi[" + std::to_string(s) + "," + std::to_string(st) + "]");
            // phi[index_phi[s][st]].set(GRB_DoubleAttr_LB,1.0);
            phi[index_phi[s][st]].set(GRB_DoubleAttr_Obj,(1.0/(double)instance.getNbScenarios(pr)));
        }
    }
}

void StrategicTransform::safety_strategy_2st(int s, int st){
    int nb_pairs = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    for(int n = 0; n < nb_pairs-1; ++n){
        int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
        int c = instance.getSchoolFeasiblePair(pr,s,i);

        GRBLinExpr exp = 0.0;
        bool nconstr = false;
        for(int j = 0; j < nb_seats; ++j){
            if(instance.getProbability(pr,s,i,j) >= 1.0){
                if(j == 0) exp += 1.0;
                else exp += solver->getX(c,j-1);
                nconstr = true; break;
            }
        }
        
        if(nconstr == true){
            // These constraints are sufficient
            if((n+1) < instance.getNbMaxRankSchools()){
                // for(int nn = n+1; nn < instance.getNbMaxRankSchools(); ++nn){
                for(int nn = 0; nn < instance.getNbMaxRankSchools(); ++nn){
                    int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,nn);
                    solver->getModel()->addConstr(getXi(s,ii) >= exp);
                }   
            }
        }
    }  
}

void StrategicTransform::constraints_2st(int s, int st){  
    // Linking variables xi with y
    for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) {
        solver->getModel()->addConstr(solver->getY(s,i) <= getXi(s,i)); 
    }

    // Unassigned pair always selected
    int i = instance.getUnassignedPair(pr,s,st);
    solver->getModel()->addConstr(getXi(s,i) == 1.0);
    // getXi(s,i).set(GRB_DoubleAttr_LB,1.0);
    // getXi(s,i).set(GRB_DoubleAttr_UB,1.0);
}

void StrategicTransform::constraints_preprocess_2st(int s, int st){
    // Preprocess, all pairs are selected
    for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) {
        solver->getModel()->addConstr(getXi(s,i) == 1.0);
        // getXi(s,i).set(GRB_DoubleAttr_LB,1.0);
        // getXi(s,i).set(GRB_DoubleAttr_UB,1.0);
    } 
    int i = instance.getUnassignedPair(pr,s,st);
    solver->getModel()->addConstr(getXi(s,i) == 1.0);
    // getXi(s,i).set(GRB_DoubleAttr_LB,1.0);
    // getXi(s,i).set(GRB_DoubleAttr_UB,1.0);
}  

void StrategicTransform::objective_2st(){
    GRBLinExpr obj = 0.0;
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            int nb_ranks = instance.getNbMaxRankSchools();
            int nb_pairs = instance.getNbAssignFeasPairsByStudent(pr,s,st);

            solver->getModel()->addConstr(getPhi(s,st) >= 1.0);
            obj += getPhi(s,st);

            if(nb_pairs <= nb_ranks && nb_pairs > 0){
                // Preprocess ranking in this case
                GRBLinExpr exp = getPhi(s,st);
                for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) {
                    exp -= (instance.getFeasPairStudentRank(pr,s,i)+1)*solver->getY(s,i);
                }
                int i = instance.getUnassignedPair(pr,s,st);
                double r = (input.getObjAccessImprov() == true) ? 
                        (std::min(instance.getFeasPairStudentRank(pr,s,i),instance.getNbMaxRankSchools())+1) : 
                        (instance.getNbSchools()+1);
                exp -= r*solver->getY(s,i);
                solver->getModel()->addConstr(exp == 0.0);
            }
            else if(nb_pairs > 0){
                // Defining variables phi used in the objective function
                for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){
                    int pref_i = instance.getFeasPairStudentRank(pr,s,i);
                    int bigM = std::min((pref_i+1),nb_ranks);
                    
                    GRBLinExpr exp = getPhi(s,st);
                    for(int nn = 0; nn <= pref_i; ++nn){
                        int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,nn);
                        exp -= getXi(s,ii);
                    }
                    // solver->getModel()->addGenConstrIndicator(solver->getY(s,i),1,exp >= 0.0);
                    for(int nn = pref_i; nn < nb_pairs; ++nn){
                        int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,nn);
                        exp -= bigM*solver->getY(s,ii);
                    }
                    int ii = instance.getUnassignedPair(pr,s,st);
                    exp -= bigM*solver->getY(s,ii);
                    solver->getModel()->addConstr(exp >= (-bigM));
                }

                // Defining objective for unassigned pair
                int i = instance.getUnassignedPair(pr,s,st);

                if(input.getObjAccessImprov() == false){
                    // Objective is access
                    solver->getModel()->addConstr(getPhi(s,st) >= (instance.getNbSchools()+1)*solver->getY(s,i));
                }
                else{
                    // If safety strategy, we cannot select less than K schools when (Rank unassigned is fixed)
                    solver->getModel()->addConstr(getPhi(s,st) >= (nb_ranks+1)*solver->getY(s,i));
                }
            }else{
                int i = instance.getUnassignedPair(pr,s,st);
                double r = (input.getObjAccessImprov() == true) ? (1) : (instance.getNbSchools()+1);
                solver->getModel()->addConstr(getPhi(s,st) == r*solver->getY(s,i));
            }
        }
    }
    solver->getModel()->addConstr(obj >= instance.getNbStudents());
}

/////////////////////////////////////////////////////////////////////////////////////////////

void StrategicTransform::set_inisol_warm_start(AbstractTransform* transf){
    StrategicTransform * strat_transf = (StrategicTransform *) transf;
    solver->getModel()->set(GRB_DoubleAttr_Start,xi,strat_transf->getXi_(),count_xi);
    solver->getModel()->set(GRB_DoubleAttr_Start,phi,strat_transf->getPhi_(),count_phi);
}

void StrategicTransform::set_inisol_warm_start(Test1stSolver* test_1st){
    solver->getModel()->set(GRB_DoubleAttr_Start,xi,test_1st->getXi_(),count_xi);
    // solver->getModel()->set(GRB_DoubleAttr_LB,xi,test_1st->getXi_(),count_xi);
    // solver->getModel()->set(GRB_DoubleAttr_UB,xi,test_1st->getXi_(),count_xi);
}

void StrategicTransform::get_variables(){
    if(xi_) delete[] xi_;
    if(phi_) delete[] phi_;
    xi_ = solver->getModel()->get(GRB_DoubleAttr_X,xi,count_xi);
    phi_ = solver->getModel()->get(GRB_DoubleAttr_X,phi,count_phi);
}

// v*xi[s,i] + c
void StrategicTransform::term_xi(int s, int i, double v, double c, GRBLinExpr & exp){
    int st = instance.getStudentFeasiblePair(pr,s,i);
    if(instance.getNbAssignFeasPairsByStudent(pr,s,st) <= instance.getNbMaxRankSchools())
        exp += v*1.0 + c;
    else if(getXi(s,i).get(GRB_DoubleAttr_LB) >= 0.99 && getXi(s,i).get(GRB_DoubleAttr_UB) >= 0.99)
        exp += v*1.0 + c;
    else
        exp += v*getXi(s,i) + c;
}

// compute feas_pairs and rank_feas_pairs according to optimized solution -> lph
void StrategicTransform::compute_feas_pairs(int s,std::vector<std::list<int>> & feas_pairs,
                                            std::unordered_map<int,double> & rank_feas_pairs){
    for(int st = 0; st < instance.getNbStudents(); ++st){
        int r = 1;
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){ 
            if(xi_[index_xi[s][i]] > 0.9999){ 
                feas_pairs[st].push_back(i); 
                rank_feas_pairs[i] = (double)r;
                ++r;
            }
        } 
        int i = instance.getUnassignedPair(pr,s,st);
        if(input.getObjAccessImprov() == false) rank_feas_pairs[i] = (instance.getNbSchools()+1);
        else rank_feas_pairs[i] = r;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////