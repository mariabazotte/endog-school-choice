#include "enumtransform.hpp"
#include "../../../solver/abstractsolver.hpp"

EnumTransform::~EnumTransform(){
    if(L) delete[] L;
    if(theta) delete[] theta;
    if(pi) delete[] pi;
}

void EnumTransform::index_2st(int s, int st){
    StrategicTransform::index_2st(s,st);

    for(int r = 0; r < instance.getNbRanks(pr,s,st); ++r){
        index_L.insert({{s,st,r},count_L});
        ++count_L;

        index_theta.insert({{s,st,r},count_theta});
        ++count_theta;
    }

    for(int r = 0; r < instance.getNbRanks(pr,s,st); ++r){
        for(int i: instance.getRank(pr,s,st,r)){
            for(int j = 0; j < instance.getNbMaxExtraSeatsWithZero(); ++j){
                index_pi.insert({{s,st,r,i,j},count_pi});
                ++count_pi;
            }
        }
    }
}

void EnumTransform::variables_2st(){
    StrategicTransform::variables_2st();

    L = solver->getModel()->addVars(count_pi, GRB_CONTINUOUS);
    if(solver->getLinearRelax() == false) theta = solver->getModel()->addVars(count_theta, GRB_BINARY);
    else theta = solver->getModel()->addVars(count_theta, GRB_CONTINUOUS);
    pi = solver->getModel()->addVars(count_pi, GRB_CONTINUOUS);

    for(auto it: index_L) L[it.second].set(GRB_StringAttr_VarName,
                        "L[" + std::to_string(std::get<0>(it.first)) + 
                        "," + std::to_string(std::get<1>(it.first)) + 
                        "," + std::to_string(std::get<2>(it.first)) + "]");

    for(auto it: index_theta) theta[it.second].set(GRB_StringAttr_VarName,
                            "theta[" + std::to_string(std::get<0>(it.first)) + 
                            "," + std::to_string(std::get<1>(it.first)) + 
                            "," + std::to_string(std::get<2>(it.first)) + "]");
    
    for(auto it: index_pi) pi[it.second].set(GRB_StringAttr_VarName,
                            "pi[" + std::to_string(std::get<0>(it.first)) + 
                            "," + std::to_string(std::get<1>(it.first)) + 
                            "," + std::to_string(std::get<2>(it.first)) + 
                            "," + std::to_string(std::get<3>(it.first)) +
                            "," + std::to_string(std::get<4>(it.first)) + "]");
}

void EnumTransform::constraints_2st(int s, int st){
    StrategicTransform::constraints_2st(s,st);
    StrategicTransform::safety_strategy_2st(s,st);

    for(int r = 0; r < instance.getNbRanks(pr,s,st); ++r) {
        
        // Probability chain constraints
        if(instance.getSizeRank(pr,s,st,r) > 1) {
            for(int k = 1; k < instance.getSizeRank(pr,s,st,r); ++k){
                int i = instance.getAssignFeasPairRank(pr,s,st,r,k);
                int ii = instance.getAssignFeasPairRank(pr,s,st,r,k-1);
                
                GRBLinExpr exp = 0;
                for(int j = 0; j < instance.getNbMaxExtraSeatsWithZero(); ++j){
                    double v = 1.0-instance.getProbability(pr,s,ii,j);
                    exp += getPi(s,st,r,i,j);
                    if(v > 1e-4) exp -= v*getPi(s,st,r,ii,j);
                }
                solver->getModel()->addConstr(exp == 0.0);
            }
        }

        for(int k = 0; k < instance.getSizeRank(pr,s,st,r); ++k) {
            int i = instance.getAssignFeasPairRank(pr,s,st,r,k);
            int c = instance.getSchoolFeasiblePair(pr,s,i);
            
            for(int j = 0; j < instance.getNbMaxExtraSeatsWithZero(); ++j){
                GRBLinExpr exp = getPi(s,st,r,i,j) - getX(c,j);
                if(k == 0) solver->getModel()->addConstr(exp == 0.0); else solver->getModel()->addConstr(exp <= 0.0);
            }
        }

        // Upper and lower bounds of the utility
        GRBLinExpr lb_utility; for(int rr = 0; rr < instance.getNbRanks(pr,s,st); ++rr) lb_utility += getL(s,st,rr);
        GRBLinExpr ub_utility = getL(s,st,r);
        for(int i: instance.getRank(pr,s,st,r)) {
            for(int j = 0; j < instance.getNbMaxExtraSeatsWithZero(); ++j){ 
                double v = instance.getProbability(pr,s,i,j)*instance.getUtility(pr,s,i);
                lb_utility -= v*getPi(s,st,r,i,j);
                ub_utility -= v*getPi(s,st,r,i,j);
            } 
        }
        solver->getModel()->addConstr(lb_utility >= 0.0);
        solver->getModel()->addConstr(ub_utility <= 0.0);

        // Big-M constraint for utility and variables theta
        solver->getModel()->addConstr(getL(s,st,r) <= instance.getMaxUtilityRank(pr,s,st,r)*getTheta(s,st,r));
    }

    // Select only one rank
    GRBLinExpr exp = 0;
    for(int r = 0; r < instance.getNbRanks(pr,s,st); ++r) exp += getTheta(s,st,r);
    solver->getModel()->addConstr(exp == 1);

    // Defining variables xi
    for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){
        GRBLinExpr exp = getXi(s,i);
        for(int r = 0; r < instance.getNbRanks(pr,s,st); ++r){
            if(instance.isAssignFeasPairInRank(pr,s,st,r,i)) exp -= getTheta(s,st,r);
        } 
        solver->getModel()->addConstr(exp == 0.0);
    }
}

/////////////////////////////////////////////////////////////////////
