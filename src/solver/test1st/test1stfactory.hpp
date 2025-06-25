#ifndef TEST_FACTORY_HPP
#define TEST_FACTORY_HPP

#include "gurobi_c++.h"
#include "../../instance/instance.hpp"
#include "../../input/input.hpp"

#include "test1stsolver.hpp"

class Test1stSolverFactory{
    public:
        inline Test1stSolver* createTest1stSolver(const Input & input, const Instance & instance, int pr, std::string name, GRBEnv *env, Input::Test1stSolver chosenSolver){
            switch (chosenSolver)
            {   
                case Input::Test1stSolver::Gurobi :{
                    return new Test1stSolverModel(input,instance,pr,name,env);
                    break;
                }
                case Input::Test1stSolver::DA :{
                    return new Test1stSolverDA(input,instance,pr,name);
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