#include "asgheuristicsolver.hpp"

void ASGheuristicSolver::warm_start(){
    // if(input.getTransform() != Input::Transform::NonStrategic && linear_relax == false) warm_start_prefrd();
    if(input.getTransform() != Input::Transform::NonStrategic && linear_relax == false) warm_start_ls();
    if(input.getTransform() == Input::Transform::NonStrategic) warm_start_capazero();
}

void ASGheuristicSolver::upd_solution(){
    // Take solution without stability constraints
    Solver::upd_solution();
    if(y_) delete[] y_;  

    // Solve deferred acceptance algorithm (DA)  
    // lb = 0.0; // A lower bound of the problem is the result of the first part of the algorithm
    ub = 0.0;
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        transf->compute_feas_pairs(s,feas_pairs[s],rank_feas_pairs[s]);
        ub += instance.da_algorithm(pr,s,t_,feas_pairs[s],rank_feas_pairs[s],matching[s])/instance.getNbScenarios(pr);
    }
    // lb = ub;
    gap = (ub-lb)/ub;

    // Taking solution y_
    y_ = new double[count_y];
    std::fill(y_,y_+count_y,0.0);
    for(int s = 0; s < instance.getNbScenarios(pr); ++s){
        for(int c = 0; c < instance.getNbSchools()+1; ++c){
            for(int i: matching[s][c]) y_[index_y[s][i]] = 1.0;
        }
    }
}