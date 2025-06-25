#include "lconstraintssolver.hpp"

void LConstraintsSolver::constraints_2st(int s){
    Solver::constraints_2st(s);

    int lazy = (input.useLazyStability() == true) ? 1 : 0;

    for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int st = instance.getStudentFeasiblePair(pr,s,i);

        GRBLinExpr exp = getT(c) + instance.getCapacity(c);
        double val = instance.getNbMaxExtraSeats() + instance.getCapacity(c);

        for(int l: instance.getStudentPreference(pr,s,st,i)) exp -= val*getY(s,l);
        for(int l: instance.getSchoolPreference(pr,s,c,i)) exp -= getY(s,l);
        transf->term_xi(s,i,val,-val,exp);

        GRBConstr constr = model->addConstr(exp <= 0.0);
        if(lazy == 1) constr.set(GRB_IntAttr_Lazy,lazy);
    }
}

void ModLConstraintsSolver::constraints_2st(int s){
    Solver::constraints_2st(s);

    int lazy = (input.useLazyStability() == true) ? 1 : 0;

    for(int i = 0; i < instance.getNbAssignFeasPairs(pr,s); ++i){
        int c = instance.getSchoolFeasiblePair(pr,s,i);
        int st = instance.getStudentFeasiblePair(pr,s,i);

        for(int ii: instance.getSchoolNotPreference(pr,s,c,i)){
            GRBLinExpr exp = -getY(s,ii);
            for(int l: instance.getStudentPreference(pr,s,st,i)) exp += getY(s,l);
            transf->term_xi(s,i,-1,1,exp);
            
            GRBConstr constr = model->addConstr(exp >= 0.0);
            if(lazy == 1) constr.set(GRB_IntAttr_Lazy,lazy);
            constrs.push_back(constr);
        } 
    }
}

void ModLConstraintsSolver::upd_solution(){
    Solver::upd_solution();
    if(linear_relax == true) dual_ = model->get(GRB_DoubleAttr_Pi,constrs.data(),constrs.size());
}

