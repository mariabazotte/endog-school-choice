#ifndef ALGO_FACTORY_HPP
#define ALGO_FACTORY_HPP

#include "../instance/instance.hpp"
#include "../input/input.hpp"
#include "algorithm.hpp"
#include "saa.hpp"

class AlgoFactory {
    public:
        inline Algorithm* createAlgorithm(const Input & input, const Instance & instance){
            Input::Type chosenType = input.getType();
            switch (chosenType)
            {
                case Input::Type::Deterministic :{
                    return new Algorithm(input,instance);
                    break;
                }
                
                case Input::Type::Stochastic :{
                    return new SAA(input,instance);
                    break;         
                }

                default:{
                    throw std::string("ERROR: Invalid Type option in command line.");
                    exit(0);
                    break;
                }
            }
            return NULL;
        }
};

#endif