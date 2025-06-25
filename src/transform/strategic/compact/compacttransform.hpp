#ifndef COMPACT_TRANSFORM_HPP
#define COMPACT_TRANSFORM_HPP

#include "gurobi_c++.h"
#include "../strategictransform.hpp"
#include <iomanip> // for std::fixed and std::setprecision

#include <map>
#include <tuple>
#include <ctime>

class CompactTransform : public StrategicTransform {
    protected: 
        GRBVar *xia = NULL;
        GRBVar *q = NULL;
        GRBVar *qa = NULL;
        GRBVar *L = NULL;
        GRBVar *v = NULL;

        int count_xia;
        int count_q;
        int count_qa;
        int count_L;
        int count_v;

        std::map<std::tuple<int,int,int>,int> index_xia;
        std::map<std::tuple<int,int,int,int>,int> index_q;
        std::map<std::tuple<int,int,int>,int> index_qa;
        std::map<std::tuple<int,int,int>,int> index_L;
        std::map<std::tuple<int,int,int>,int> index_v;

        GRBVar & getXia(int s, int i, int k) { return xia[index_xia[{s,i,k}]]; }
        GRBVar & getQ(int s, int i, int j, int k) { return q[index_q[{s,i,j,k}]]; }
        GRBVar & getQa(int s, int i, int k) { return qa[index_qa[{s,i,k}]]; }
        GRBVar & getL(int s, int st, int k) { return L[index_L[{s,st,k}]]; }
        GRBVar & getV(int s, int i, int k) { return v[index_v[{s,i,k}]]; }

        void index_2st(int,int);
        void variables_2st();
        void constraints_2st(int,int);

        void assignment_2st(int,int);
        void probability_chain_2st(int,int);

        void recurr_cut_max_utility_2st(int,int);
        void recurr_def_utility_2st(int,int);

        void preprocess_2st(int,int);
        
    public:
        CompactTransform(const Input & input, const Instance & instance, int pr, Solver *solver) : 
                                            StrategicTransform(input,instance,pr,solver) {
            count_xia = 0;
            count_q = 0;
            count_qa = 0;
            count_L = 0;
            count_v = 0;
        }
        ~CompactTransform();
};

#endif