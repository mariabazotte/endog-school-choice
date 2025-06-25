#include "algorithm.hpp"

void Algorithm::create(){
    Test1stSolverFactory testFactory;
    test_stochastic = testFactory.createTest1stSolver(input,instance,instance.getValidateProblem(),"Test Problem",env_test,Input::Test1stSolver::DA);

    SolverFactory factory;
    solvers.resize(instance.getNbProblems(),NULL);
    for(int i = 0; i < instance.getNbProblems(); ++i){
        solvers[i] = factory.createSolver(input,instance,i,"Problem " + std::to_string(i+1),env);
    }
}

void Algorithm::solve(){
    solvers[0]->solve();
    status = solvers[0]->getStatus();
    if(status == Status::Optimal || status == Status::Time_Limit || status == Status::Iteration_Limit){
        lb = solvers[0]->getLB();
        ub = solvers[0]->getUB();
        gap = solvers[0]->getGAP();
        time_ = solvers[0]->getTime();
        nb_bnb = solvers[0]->getNbBnB();
        
        if(input.getTransform() == Input::Transform::NonStrategic){
            computeTest(0,ub_test,var_test,time_test,avg_stu_per_rank,nb_enter,nb_improv,Input::Transform::NonStrategic); // nonstrategic
            computeTest(0,ub_test_strat_simple,var_test_strat_simple,time_test_strat_simple,avg_stu_per_rank_strat_simple,nb_enter_strat_simple,nb_improv_strat_simple,Input::Transform::Simple); // simple reported
            computeTest(0,ub_test_strat,var_test_strat,time_test_strat,avg_stu_per_rank_strat,nb_enter_strat,nb_improv_strat,Input::Transform::Compact); // complex reported
        }else if(input.getTransform() == Input::Transform::Simple){
            computeTest(0,ub_test_strat_simple,var_test_strat_simple,time_test_strat_simple,avg_stu_per_rank_strat_simple,nb_enter_strat_simple,nb_improv_strat_simple,Input::Transform::Simple); // simple reported
            computeTest(0,ub_test_strat,var_test_strat,time_test_strat,avg_stu_per_rank_strat,nb_enter_strat,nb_improv_strat,Input::Transform::Compact); // complex reported
        }else{
            computeTest(0,ub_test_strat,var_test_strat,time_test_strat,avg_stu_per_rank_strat,nb_enter_strat,nb_improv_strat,Input::Transform::Compact); // complex reported
        }
    }
}

void Algorithm::computeTest(int pr, double & ub_, double & var_, double & time_test_, std::vector<double> & avg_rank, double & nb_enter_test, double & nb_improv_test, Input::Transform type_transform){
    double start_time = time(NULL);
    test_stochastic->solve(solvers[pr]->getT_(),type_transform);
    if(test_stochastic->getStatus() == Status::Optimal) 
        ub_ = test_stochastic->getUB();
    computeVariance(ub_,var_);
    for(int r = 0; r <= instance.getNbMaxRankSchools(); ++r) 
        avg_rank[r] = test_stochastic->getAvgRank(r);
    time_test_ = time(NULL) - start_time;

    // Solve test_stochastic with no extra capacity to compute nb entering and nb improving
    nb_enter_test = 0.0;
    nb_improv_test = 0.0;
    std::vector<std::vector<double>> rank_t_;
    rank_t_.resize(instance.getNbValidateScenarios());
    for(int s = 0; s < instance.getNbValidateScenarios(); ++s) {
        for(int st = 0; st < instance.getNbStudents(); ++st){
            rank_t_[s].push_back(test_stochastic->getRankStudent(s,st));
        }
    }
    double * zero_capa = new double[instance.getNbSchools()];
    std::fill(zero_capa,zero_capa+instance.getNbSchools(),0.0);
    test_stochastic->setVerbose(0);
    test_stochastic->solve(zero_capa,type_transform);
    test_stochastic->setVerbose(1);
    for(int s = 0; s < instance.getNbValidateScenarios(); ++s){
        for(int st = 0; st < instance.getNbStudents(); ++st){
            int unassigned = std::min(instance.getNbAssignFeasPairsByStudent(instance.getValidateProblem(),s,st),instance.getNbMaxRankSchools())+1;
            if((int)test_stochastic->getRankStudent(s,st) >= unassigned){
                if((int)rank_t_[s][st] < unassigned)
                    nb_enter_test += 1.0/instance.getNbValidateScenarios();

            }else{
                if(test_stochastic->getRankStudent(s,st) > rank_t_[s][st]) 
                    nb_improv_test += 1.0/instance.getNbValidateScenarios();
            }
        }
    }
}

void Algorithm::computeVariance(double & ub_,double & var_){
    var_ = 0.0;
    for(int s = 0; s < instance.getNbValidateScenarios(); ++s)
        var_ += std::pow(test_stochastic->getObjScenario(s)-ub_,2);
    long long N = instance.getNbValidateScenarios();
    if(N > 1) var_ = var_/(N*(N-1));
}