#ifndef ENUM_TRANSFORM_HPP
#define ENUM_TRANSFORM_HPP

#include "gurobi_c++.h"
#include "../strategictransform.hpp"

#include <map>
#include <tuple>
#include <ctime>

class EnumTransform : public StrategicTransform {
    protected:
        GRBVar *L = NULL;
        GRBVar *theta = NULL; 
        GRBVar *pi = NULL;

        int count_L;
        int count_theta;
        int count_pi;

        std::map<std::tuple<int,int,int>,int> index_L;
        std::map<std::tuple<int,int,int>,int> index_theta;
        std::map<std::tuple<int,int,int,int,int>,int> index_pi;

        GRBVar & getL(int s, int st, int r) { return L[index_L[{s,st,r}]]; }
        GRBVar & getTheta(int s, int st, int r) { return theta[index_theta[{s,st,r}]]; }
        GRBVar & getPi(int s, int st, int r, int i, int j) { return pi[index_pi[{s,st,r,i,j}]]; }

        void index_2st(int,int);
        void variables_2st();
        void constraints_2st(int,int);

    public:
        EnumTransform(const Input & input, const Instance & instance, int pr, Solver *solver) : 
                            StrategicTransform(input,instance,pr,solver) {
            count_L = 0;
            count_theta = 0;
            count_pi = 0;
        }
        ~EnumTransform();  
};

#endif