#include "nonstrategictransform.hpp"
#include "../../solver/abstractsolver.hpp"
#include "../../solver/test1st/test1stsolver.hpp"

void NonStrategicTransform::transf_2st(){ 
    // We only need to define the objective function
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            int num = std::min(instance.getNbMaxRankSchools(),instance.getNbAssignFeasPairsByStudent(pr,s,st));
            for(int r = 0; r < num; ++r){
                int i = instance.getAssignFeasPairByStudentPos(pr,s,st,r);
                solver->getY(s,i).set(GRB_DoubleAttr_Obj,(double)(r+1)/(double)instance.getNbScenarios(pr));
            }
            int i = instance.getUnassignedPair(pr,s,st);
            if(input.getObjAccessImprov() == false) num = instance.getNbSchools();
            solver->getY(s,i).set(GRB_DoubleAttr_Obj,(double)(num+1)/(double)instance.getNbScenarios(pr));
        }
    }
}

void NonStrategicTransform::set_constr_warm_start(Test1stSolver * test_1st){
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            int i = test_1st->getMatch(s,st);
            int pref_i = instance.getFeasPairStudentRank(pr,s,i);
            int nb_assign_fp = instance.getNbAssignFeasPairsByStudent(pr,s,st);

            // The students matching can only be improved 
            if(instance.isUnassignedPair(pr,s,st,i) == false){
                GRBLinExpr exp;
                for(int n = 0; n <= pref_i; ++n){
                    int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
                    exp += solver->getY(s,ii);
                }
                solver->getModel()->addConstr(exp >= 1.0);

                if(pref_i + 1 < nb_assign_fp){ // Setting worst matches to zero
                    for(int n = pref_i + 1; n < nb_assign_fp; ++n){
                        int ii = instance.getAssignFeasPairByStudentPos(pr,s,st,n);
                        solver->getY(s,ii).set(GRB_DoubleAttr_UB,0.0);
                        solver->getModel()->addConstr(solver->getY(s,ii) == 0.0);
                    }
                }
                int ii = instance.getUnassignedPair(pr,s,st);
                solver->getY(s,ii).set(GRB_DoubleAttr_UB,0.0);
                solver->getModel()->addConstr(solver->getY(s,ii) == 0.0);
            } 

            // If student was matched with his most preferred school, we can set this matching to 1.0
            if(pref_i == 0){
                solver->getY(s,i).set(GRB_DoubleAttr_LB,1.0);
                solver->getY(s,i).set(GRB_DoubleAttr_UB,1.0);
                solver->getModel()->addConstr(solver->getY(s,i) == 1.0);
            }
        }

        // If a school is unsatured, we can set its extra seats to zero
        for(int c = 0; c < instance.getNbSchools(); ++c){ 
            if(test_1st->getUnsatured(c) == true){
                solver->getT(c).set(GRB_DoubleAttr_LB,0.0);
                solver->getT(c).set(GRB_DoubleAttr_UB,0.0);
                solver->getModel()->addConstr(solver->getT(c) == 0.0);
            }
        }
    }
}

// v*xi[s,i] + c
void NonStrategicTransform::term_xi(int s, int i, double v, double c, GRBLinExpr & exp){
    UNUSED(s);
    UNUSED(i);
    // xi always equal to one in the non-strategic version
    exp += v*1.0 + c;
}

void NonStrategicTransform::compute_feas_pairs(int s,std::vector<std::list<int>> & feas_pairs,
                                            std::unordered_map<int,double> & rank_feas_pairs){
    for(int st = 0; st < instance.getNbStudents(); ++st){
        int r = 1;
        for(int i: instance.getAssignFeasPairsByStudent(pr,s,st)){ 
            feas_pairs[st].push_back(i); 
            rank_feas_pairs[i] = r;
            ++r;
        } 
        rank_feas_pairs[instance.getUnassignedPair(pr,s,st)] = r;
    }
}
