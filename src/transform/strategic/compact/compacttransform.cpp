#include "compacttransform.hpp"
#include "../../../solver/abstractsolver.hpp"

CompactTransform::~CompactTransform(){
    if(xia) delete[] xia;
    if(q) delete[] q;
    if(qa) delete[] qa;
    if(L) delete[] L;
    if(v) delete[] v;
}

void CompactTransform::index_2st(int s, int st){
    StrategicTransform::index_2st(s,st);

    std::vector<int> assign_feas_pairs = instance.getAssignFeasPairsByStudent(pr,s,st);

    if(assign_feas_pairs.empty() == false && (int)assign_feas_pairs.size() > instance.getNbMaxRankSchools()){
        for(int i: assign_feas_pairs){
            for(int k = 0; k < instance.getNbMaxRankSchools()-1; ++k){
                index_xia.insert({{s,i,k},count_xia});
                ++count_xia;
            }
        }

        for(int i: assign_feas_pairs){
            for(int j = 0; j < instance.getNbMaxExtraSeatsWithZero(); ++j){
                for(int k = 0; k < instance.getNbMaxRankSchools(); ++k){
                    index_q.insert({{s,i,j,k},count_q});
                    ++count_q;
                }
            }
        }

        for(int i: assign_feas_pairs){
            for(int k = 0; k < instance.getNbMaxRankSchools(); ++k){
                index_qa.insert({{s,i,k},count_qa});
                ++count_qa;

                index_v.insert({{s,i,k},count_v});
                ++count_v;
            }
        }

        for(int k = 0; k < instance.getNbMaxRankSchools(); ++k){
            index_L.insert({{s,st,k},count_L});
            ++count_L;
        }
    }
}

void CompactTransform::variables_2st(){
    StrategicTransform::variables_2st();

    if(solver->getLinearRelax() == false) xia = solver->getModel()->addVars(count_xia,GRB_BINARY);
    else xia = solver->getModel()->addVars(count_xia,GRB_CONTINUOUS);
    q = solver->getModel()->addVars(count_q,GRB_CONTINUOUS);
    qa = solver->getModel()->addVars(count_qa,GRB_CONTINUOUS);
    L = solver->getModel()->addVars(count_L,GRB_CONTINUOUS);
    v = solver->getModel()->addVars(count_v,GRB_CONTINUOUS);

    for(auto it: index_xia){
        int s = std::get<0>(it.first);
        int i = std::get<1>(it.first);
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int st = instance.getStudentFeasiblePair(pr,s,i);
        xia[it.second].set(GRB_StringAttr_VarName,
                    "xia[" + std::to_string(std::get<0>(it.first)) + 
                    "," + std::to_string(st) + "," + std::to_string(c) +
                    "," + std::to_string(std::get<2>(it.first)) + "]");
        xia[it.second].set(GRB_DoubleAttr_LB,0.0);
        // xia[it.second].set(GRB_DoubleAttr_UB,1.0);
        if(solver->getLinearRelax() == true) solver->getModel()->addConstr(xia[it.second] <= 1.0);
    }

    for(auto it: index_q){
        int c = instance.getSchoolFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        int st = instance.getStudentFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        q[it.second].set(GRB_StringAttr_VarName,
                    "q[" + std::to_string(std::get<0>(it.first)) + 
                    "," + std::to_string(st) + "," + std::to_string(c) +
                    "," + std::to_string(std::get<2>(it.first)) +
                    "," + std::to_string(std::get<3>(it.first)) + "]");
        q[it.second].set(GRB_DoubleAttr_LB,0.0);
        // q[it.second].set(GRB_DoubleAttr_UB,100.0);
    }

    for(auto it: index_qa){
        int c = instance.getSchoolFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        int st = instance.getStudentFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        qa[it.second].set(GRB_StringAttr_VarName,
                    "qa[" + std::to_string(std::get<0>(it.first)) + 
                    "," + std::to_string(st) + "," + std::to_string(c) +
                    "," + std::to_string(std::get<2>(it.first)) + "]");
        qa[it.second].set(GRB_DoubleAttr_LB,0.0);
        // qa[it.second].set(GRB_DoubleAttr_UB,100.0);
    }

    for(auto it: index_L){
        L[it.second].set(GRB_StringAttr_VarName,
                    "L[" + std::to_string(std::get<0>(it.first)) + 
                    "," + std::to_string(std::get<1>(it.first)) +
                    "," + std::to_string(std::get<2>(it.first)) + "]");
        L[it.second].set(GRB_DoubleAttr_LB,0.0);
        L[it.second].set(GRB_DoubleAttr_UB,GRB_INFINITY);
    }

    for(auto it: index_v){
        int c = instance.getSchoolFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        int st = instance.getStudentFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        v[it.second].set(GRB_StringAttr_VarName,
            "v[" + std::to_string(std::get<0>(it.first)) + 
            "," + std::to_string(st) + "," + std::to_string(c) +
            "," + std::to_string(std::get<2>(it.first)) + "]");
    }
}

///////////////////////////////////////////////////////////////////////////

void CompactTransform::constraints_2st(int s, int st){
    StrategicTransform::constraints_2st(s,st);

    assignment_2st(s,st);
    probability_chain_2st(s,st);
    safety_strategy_2st(s,st);

    recurr_def_utility_2st(s,st);
    recurr_cut_max_utility_2st(s,st);
    
    preprocess_2st(s,st);
}

///////////////////////////////////////////////////////////////////////////

void CompactTransform::assignment_2st(int s, int st){
    int nb_ranks = instance.getNbMaxRankSchools();

    // Assignment constraints
    for(int k = 0; k < nb_ranks; ++k){
        GRBLinExpr exp = 0;
        if(k == nb_ranks-1) for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) exp += getXi(s,i);
        else for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) exp += getXia(s,i,k);
        if(k > 0) for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)) exp -= getXia(s,i,k-1);
        solver->getModel()->addConstr(exp == 1.0);
    }

    for(int k = 0; k < nb_ranks; ++k){
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){
            GRBLinExpr exp = 0.0;
            if(k == nb_ranks-1) exp -= getXi(s,i);
            else exp -= getXia(s,i,k);
            if(k > 0) exp += getXia(s,i,k-1);
            solver->getModel()->addConstr(exp <= 0.0);
        }
    }
}

void CompactTransform::probability_chain_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    // Probability chain constraints
    for(int k = 0; k < nb_ranks; ++k) {
        if(nb_assign_fp >= 2){
            for(int n = 1; n < nb_assign_fp; ++n){
                int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
                int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,n-1);
                
                GRBLinExpr exp = getQa(s,i,k) - getQa(s,ii,k);
                for(int j = 0; j < nb_seats; ++j){
                    double v = (1.0-instance.getProbability(pr,s,ii,j));
                    exp += getQ(s,i,j,k);
                    if(v >= 1e-6) exp -= v*getQ(s,ii,j,k);
                }
                solver->getModel()->addConstr(exp == 0.0);
            }  
        }           

        for(int n = 0; n < nb_assign_fp; ++n){
            int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
            
            GRBLinExpr exp = 0.0;
            if(k == nb_ranks-1) exp -= 100*getXi(s,i);
            else exp -= 100*getXia(s,i,k);
            for(int j = 0; j < nb_seats; ++j) exp += getQ(s,i,j,k);
            if(n == 0){solver->getModel()->addConstr(exp == 0.0);}else{solver->getModel()->addConstr(exp <= 0.0);}

            GRBLinExpr exp_qa = getQa(s,i,k);
            if(k == nb_ranks-1) exp_qa += 100*getXi(s,i);
            else exp_qa += 100*getXia(s,i,k);
            if(n == 0){solver->getModel()->addConstr(exp_qa == 100.0);}else{solver->getModel()->addConstr(exp_qa <= 100.0);}

            int c = instance.getSchoolFeasiblePair(pr,s,i);
            for(int j = 0; j < nb_seats; ++j) solver->getModel()->addConstr(getQ(s,i,j,k) <= 100*getX(c,j));
        } 
    }
}

void CompactTransform::preprocess_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    for(int n = 0; n < nb_assign_fp; ++n){
        int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
        double v = instance.getUtility(pr,s,i)*instance.getProbability(pr,s,i,0);

        int count = 0;
        if(n+1 < nb_assign_fp){
            for(int nn = n+1; nn < nb_assign_fp; ++nn){
                int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,nn);
                double vv = instance.getUtility(pr,s,ii)*instance.getProbability(pr,s,ii,nb_seats-1);
                
                if(v >= vv + 1e-6) {
                    ++count;
                    // for(int k = 0; k < nb_ranks-1; ++k){
                    //     solver->getModel()->addConstr(getXia(s,i,k) >= getXia(s,ii,k));
                    // }
                    // solver->getModel()->addConstr(getXi(s,i) >= getXi(s,ii));
                }
            }
        }
        if(count >= (nb_assign_fp-nb_ranks)) {
            solver->getModel()->addConstr(getXi(s,i) == 1.0);
            // getXi(s,i).set(GRB_DoubleAttr_LB,1.0);
            // getXi(s,i).set(GRB_DoubleAttr_UB,1.0);
        }
    }
}

/////////////////////////////////////////////////////////////////////////

void CompactTransform::recurr_def_utility_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    // Define variables v
    for(int k = 0; k < nb_ranks; ++k) {
        for(int n = 0; n < nb_assign_fp; ++n){
            int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
            GRBLinExpr exp = getV(s,i,k);
            for(int j = 0; j < nb_seats; ++j){
                double v = instance.getProbability(pr,s,i,j)*instance.getUtility(pr,s,i);
                if(v >= 1e-6) exp -= v*getQ(s,i,j,k);
            }
            solver->getModel()->addConstr(exp == 0.0);
        }
    }

    // Defining variables L
    for(int k = 0; k < nb_ranks; ++k) {
        GRBLinExpr exp = getL(s,st,k);
        if(k > 0) exp -= getL(s,st,k-1);

        for(int n = 0; n < nb_assign_fp; ++n){
            int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
            exp -= getV(s,i,k);
            if(k > 0) exp += getV(s,i,k-1);
        }
        solver->getModel()->addConstr(exp == 0.0);
    }
}

void CompactTransform::recurr_cut_max_utility_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    for(int k = 1; k < nb_ranks; ++k){
        solver->getModel()->addConstr(getL(s,st,k) >= getL(s,st,k-1)); 
    }

    for(int k = 0; k < nb_ranks; ++k){
        for(int n = 0; n < nb_assign_fp; ++n){ 
            int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
            int c = instance.getSchoolFeasiblePair(pr,s,i);

            for(int j = 0; j < nb_seats; ++j){
                double p_i = instance.getProbability(pr,s,i,j);
                double v_i = p_i*instance.getUtility(pr,s,i);
                if(k == 0){
                    if(v_i >= 1e-6){
                        if(j == 0) solver->getModel()->addConstr(getL(s,st,k) >= v_i*100); 
                        else solver->getModel()->addConstr(getL(s,st,k) >= v_i*100*solver->getX(c,j-1));  
                    } 
                }
                else{
                    if(v_i >= 1e-6) {
                        GRBLinExpr exp = getL(s,st,k) - getL(s,st,k-1);
                        exp -= v_i*getQa(s,i,k-1);
                        if(j > 0) exp -= v_i*(100*solver->getX(c,j-1) - 100.0);  // if j == 0 then (x_c^0) = 1.0 and (x_c^0)-1.0 = 0.0

                        if(n+1 < nb_assign_fp){
                            for(int nn = n+1; nn < nb_assign_fp; ++nn){
                                int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,nn);
                                exp += p_i*getV(s,ii,k-1);
                            }
                        } 
                        solver->getModel()->addConstr(exp >= 0.0);    
                    }               
                }
            }
        }
    }

    if(solver->getLinearRelax() == false){
        // Constraint as implication
        bool use_imp = false;
        if(nb_ranks >= 2 && use_imp == true){
            for(int k = 1; k < nb_ranks; ++k){
                for(int n = 0; n < nb_assign_fp; ++n){ 
                    int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
                    int c = instance.getSchoolFeasiblePair(pr,s,i);
                    for(int j = 0; j < nb_seats; ++j){
                        double p_i = instance.getProbability(pr,s,i,j);
                        double v_i = p_i*instance.getUtility(pr,s,i);
                        if(v_i >= 1e-6) {
                            GRBLinExpr exp = getL(s,st,k) - getL(s,st,k-1);
                            exp -= v_i*getQa(s,i,k-1);
                            if(n+1 < nb_assign_fp){
                                for(int nn = n+1; nn < nb_assign_fp; ++nn){
                                    int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,nn);
                                    exp += p_i*getV(s,ii,k-1);
                                }
                            }
                            if(j > 0) solver->getModel()->addGenConstrIndicator(solver->getX(c,j-1),1,exp >= 0.0);
                        }
                    }
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////