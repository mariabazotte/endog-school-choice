#ifndef SIMPLE_TRANSFORM_HPP
#define SIMPLE_TRANSFORM_HPP

#include "gurobi_c++.h"
#include "../strategictransform.hpp"


class SimpleTransform : public StrategicTransform {
    protected:
        GRBVar* xia = NULL; 
        GRBVar* L = NULL;

        int count_xia;
        int count_L;

        std::map<std::tuple<int,int,int,int>,int> index_xia;
        std::map<std::tuple<int,int,int>,int> index_L;

        GRBVar & getXia(int s, int i, int j, int k) { return xia[index_xia[{s,i,j,k}]]; }
        GRBVar & getL(int s, int st, int k) { return L[index_L[{s,st,k}]]; }

        void index_2st(int,int);
        void variables_2st();
        void constraints_2st(int,int);

        void assignment_2st(int,int);
        void cut_2st(int,int);
        void linking_2st(int,int);
        void preprocess_2st(int,int);

    public:
        SimpleTransform(const Input & input, const Instance & instance, int pr, Solver *solver) : 
                                            StrategicTransform(input,instance,pr,solver) {
            count_xia = 0;
            count_L = 0;
        }
        ~SimpleTransform();
};

#endif