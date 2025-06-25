#ifndef SOLVER_FACTORY_HPP
#define SOLVER_FACTORY_HPP

#include "gurobi_c++.h"
#include "../instance/instance.hpp"
#include "../input/input.hpp"

#include "abstractsolver.hpp"
#include "lconstraints/lconstraintssolver.hpp"
#include "cutoffscore/cutoffscoresolver.hpp"
#include "asgheuristic/asgheuristicsolver.hpp"
#include "lagrangian/lagrangian.hpp"
#include "localsearch/localsearch.hpp"
#include "integerlshaped/integerlshaped.hpp"
#include "simulannealing/simulannealing.hpp"

class SolverFactory{
    public:
        inline Solver* createSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env){
            Input::Solver chosenSolver = input.getSolver();
            switch (chosenSolver)
            {   
                case Input::Solver::LConstraints :{
                    return new LConstraintsSolver(input,instance,pr,name,env,input.useWarmStart(),input.getLinearRelax());
                    break;
                }
                case Input::Solver::ModLConstraints :{
                    return new ModLConstraintsSolver(input,instance,pr,name,env,input.useWarmStart(),input.getLinearRelax());
                    break;
                }
                case Input::Solver::CutoffScore :{
                    return new CutoffScoreSolver(input,instance,pr,name,env,input.useWarmStart(),input.getLinearRelax());
                    break;
                }
                case Input::Solver::ASGHeuristic :{
                    return new ASGheuristicSolver(input,instance,pr,name,env,input.useWarmStart(),input.getLinearRelax());
                    break;         
                }
                case Input::Solver::Lagrangian :{
                    return new Lagrangian(input,instance,pr,name,env,input.useWarmStart());
                    break;         
                }
                case Input::Solver::LocalSearch:{
                    return new LocalSearchSolver(input,instance,pr,name,env,input.useWarmStart(),false);
                    break;  
                }
                case Input::Solver::IntegerLShaped:{
                    return new IntegerLShaped(input,instance,pr,name,env,input.useWarmStart());
                    break;  
                }
                case Input::Solver::SimulAnnealing:{
                    return new SimulAnnealingSolver(input,instance,pr,name,env,input.useWarmStart());
                    break;  
                }
                default:{
                    throw std::string("ERROR: Invalid Solver option in command line.");
                    exit(0);
                    break;
                }
            }
            return NULL;
        }
};

#endif