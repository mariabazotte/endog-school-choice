#include "cutoffscoresolver.hpp"

void CutoffScoreSolver::index_2st(int s){
    Solver::index_2st(s);
    for(int c = 0; c < instance.getNbSchools(); ++c){
        index_zf.insert({{s,c},count_zf});
        ++count_zf;
    }
}

void CutoffScoreSolver::variables_2st(){
    Solver::variables_2st();

    z = model->addVars(count_zf,GRB_CONTINUOUS);
    f = model->addVars(count_zf,GRB_BINARY);
    
    for(auto it: index_zf) z[it.second].set(GRB_StringAttr_VarName,
                                        "z[" + std::to_string(std::get<0>(it.first)) + 
                                        "," + std::to_string(std::get<1>(it.first)) + "]");

    for(auto it: index_zf) f[it.second].set(GRB_StringAttr_VarName,
                                        "f[" + std::to_string(std::get<0>(it.first)) + 
                                        "," + std::to_string(std::get<1>(it.first)) + "]");
}

void CutoffScoreSolver::constraints_2st(int s){
    Solver::constraints_2st(s);

    int nb_schools = instance.getNbSchools();
    int nb_seats = instance.getNbMaxExtraSeats();
    int nb_students = instance.getNbStudents();

    int lazy = (input.useLazyStability() == true) ? 1 : 0;

    // Cutoff constraints for stability
    for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
        int v = nb_students + 1;
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int r = instance.getScore(pr,s,i) + nb_students + 1;

        GRBConstr constr = model->addConstr(getZ(s,c) + v*getY(s,i) <= r);
        if(lazy == 1) constr.set(GRB_IntAttr_Lazy,lazy);
    }

    for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int st = instance.getStudentFeasiblePair(pr,s,i);

        int v1 = instance.getScore(pr,s,i) + instance.getEpsilon();
        int v2 = nb_students + 1;

        GRBLinExpr exp = -getZ(s,c);
        transf->term_xi(s,i,v1,0.0,exp);
        for(int l: instance.getStudentPreference(pr,s,st,i)) exp -= v2*getY(s,l);

        GRBConstr constr = model->addConstr(exp <= 0.0);
        if(lazy == 1) constr.set(GRB_IntAttr_Lazy,lazy);
    }
        
    // Waste constraints
    for(int c = 0; c < nb_schools; ++c){
        GRBLinExpr exp = instance.getCapacity(c) + getT(c) - (1-getF(s,c))*(instance.getCapacity(c)+nb_seats);
        for(int i : instance.getAssignFeasPairsBySchool(pr,s,c)) exp -= getY(s,i);
        
        GRBConstr constr = model->addConstr(exp <= 0.0);
        if(lazy == 1) constr.set(GRB_IntAttr_Lazy,lazy);
    }

    for(int c = 0; c < nb_schools; ++c){
        GRBConstr constr = model->addConstr(getZ(s,c) -(nb_students+1)*getF(s,c) <= 0.0);
        if(lazy == 1) constr.set(GRB_IntAttr_Lazy,lazy);
    }
}
