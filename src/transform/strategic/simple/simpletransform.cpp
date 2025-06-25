#include "simpletransform.hpp"
#include "../../../solver/abstractsolver.hpp"

SimpleTransform::~SimpleTransform(){
    delete[] xia;
    delete[] L;
}

void SimpleTransform::index_2st(int s, int st){
    StrategicTransform::index_2st(s,st);

    std::vector<int> assign_feas_pairs = instance.getAssignFeasPairsByStudent(pr,s,st);

    if(assign_feas_pairs.empty() == false && (int)assign_feas_pairs.size() > instance.getNbMaxRankSchools()){
        
        for(int i: assign_feas_pairs){
            for(int j = 0; j < instance.getNbMaxExtraSeatsWithZero(); ++j){
                for(int k = 0; k < instance.getNbMaxRankSchools(); ++k){
                    index_xia.insert({{s,i,j,k},count_xia});
                    ++count_xia;
                }
            }
        }

        for(int k = 0; k < instance.getNbMaxRankSchools(); ++k){
            index_L.insert({{s,st,k},count_L});
            ++count_L;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void SimpleTransform::variables_2st(){
    StrategicTransform::variables_2st();

    if(solver->getLinearRelax() == false) xia = solver->getModel()->addVars(count_xia,GRB_BINARY);
    else xia = solver->getModel()->addVars(count_xia,GRB_CONTINUOUS);
    L = solver->getModel()->addVars(count_L,GRB_CONTINUOUS);

    for(auto it: index_xia){
        int c = instance.getSchoolFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        int st = instance.getStudentFeasiblePair(pr,std::get<0>(it.first),std::get<1>(it.first));
        xia[it.second].set(GRB_StringAttr_VarName,
                    "xia[" + std::to_string(std::get<0>(it.first)) + 
                    "," + std::to_string(st) + "," + std::to_string(c) +
                    "," + std::to_string(std::get<2>(it.first)) +
                    "," + std::to_string(std::get<3>(it.first)) + "]");
        xia[it.second].set(GRB_DoubleAttr_LB,0.0);
        if(solver->getLinearRelax() == true) solver->getModel()->addConstr(xia[it.second] <= 1.0);
    }

    for(auto it: index_L){
        L[it.second].set(GRB_StringAttr_VarName,
                    "L[" + std::to_string(std::get<0>(it.first)) + 
                    "," + std::to_string(std::get<1>(it.first)) +
                    "," + std::to_string(std::get<2>(it.first)) + "]");
    }
}

void SimpleTransform::assignment_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    // Assignment constraints
    for(int k = 0; k < nb_ranks; ++k){
        GRBLinExpr exp;
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){ 
            for(int j = 0; j < nb_seats; ++j) exp += getXia(s,i,j,k);
        }
        solver->getModel()->addConstr(exp == 1.0);
    }

    for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){
        GRBLinExpr exp;
        for(int k = 0; k < nb_ranks; ++k) { 
            for(int j = 0; j < nb_seats; ++j) exp += getXia(s,i,j,k);
        }
        solver->getModel()->addConstr(exp <= 1.0);
    }

    // Linking with variables x
    for(int n = 0; n < nb_assign_fp; ++n){
        int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        for(int j = 0; j < nb_seats; ++j){
            GRBLinExpr exp = -getX(c,j);
            for(int k = 0; k < nb_ranks; ++k) exp += getXia(s,i,j,k);
            solver->getModel()->addConstr(exp <= 0.0);
        }
    }
}

void SimpleTransform::cut_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    // Cut constraints: definition of L
    for(int k = 0; k < nb_ranks; ++k){
        GRBLinExpr exp = getL(s,st,k);
        for(int n = 0; n < nb_assign_fp; ++n){
            int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
            for(int j = 0; j < nb_seats; ++j){
                double v = instance.getProbability(pr,s,i,j)*instance.getUtility(pr,s,i);
                if(v >= 1e-6) exp -= v*100*getXia(s,i,j,k);
            }
        }
        solver->getModel()->addConstr(exp == 0.0);
    }
    
    for(int k = 0; k < nb_ranks; ++k){
        for(int n = 0; n < nb_assign_fp; ++n){
            int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            
            GRBLinExpr exp = getL(s,st,k);
            for(int j = 0; j < nb_seats; ++j){
                double v = instance.getProbability(pr,s,i,j)*instance.getUtility(pr,s,i);
                if(v >= 1e-6){
                    exp -= v*100*getX(c,j);
                    if(k > 0){ for(int kk = 0; kk < k; ++kk){ exp += v*100*getXia(s,i,j,kk); }}
                }
            }
            solver->getModel()->addConstr(exp >= 0.0);   
        }
    }
}

void SimpleTransform::linking_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    // Linking with variables xi
    for(int n = 0; n < nb_assign_fp; ++n){
        int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
        GRBLinExpr exp = getXi(s,i);
        for(int j = 0; j < nb_seats; ++j){
            for(int k = 0; k < nb_ranks; ++k) exp -= getXia(s,i,j,k);
        }
        solver->getModel()->addConstr(exp == 0.0);
    }
}

void SimpleTransform::preprocess_2st(int s, int st){
    int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);
    int nb_ranks = instance.getNbMaxRankSchools(); 
    int nb_seats = instance.getNbMaxExtraSeatsWithZero();

    // See if we can fix some variables
    for(int n = 0; n < nb_assign_fp; ++n){
        int i = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
        double v = instance.getUtility(pr,s,i)*instance.getProbability(pr,s,i,0);
        
        int nb = 0;
        for(int nn = 0; nn < nb_assign_fp; ++nn){
            int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,nn);
            double vv = instance.getUtility(pr,s,ii)*instance.getProbability(pr,s,ii,nb_seats-1);
            if(nn != n && v >= vv + 1e-6) {
                ++nb;
                // See if constraint can help
                // for(int k = 0; k < nb_ranks; ++k){
                //     GRBLinExpr exp;
                //     for(int kk = 0; kk <= k; ++kk){
                //         for(int j = 0; j < nb_seats; ++j) 
                //             exp += getXia(s,i,j,kk) - getXia(s,ii,j,kk);
                //     }
                //     solver->getModel()->addConstr(exp >= 0.0);
                // }
            }
        }
        
        if(nb >= (nb_assign_fp-nb_ranks)) {
            solver->getModel()->addConstr(getXi(s,i) == 1.0);
            // getXi(s,i).set(GRB_DoubleAttr_LB,1.0);
            // getXi(s,i).set(GRB_DoubleAttr_UB,1.0);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void SimpleTransform::constraints_2st(int s, int st){
    StrategicTransform::constraints_2st(s,st);

    assignment_2st(s,st);
    cut_2st(s,st);
    linking_2st(s,st);
    preprocess_2st(s,st);
}

///////////////////////////////////////////////////////////////////