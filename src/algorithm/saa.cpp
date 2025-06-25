#include "saa.hpp"

void SAA::solve(){
   double start_time = time(NULL); 
   for(int i = 0; i < instance.getNbProblems(); ++i){
      solvers[i]->solve();
      Status curr_status = solvers[i]->getStatus(); 
      status = std::max(status,curr_status);
      if(curr_status == Status::Optimal || curr_status == Status::Time_Limit || status == Status::Iteration_Limit){
         lb_vector[i] = solvers[i]->getLB();
         ub_vector[i] = solvers[i]->getUB();
         gap_vector[i] = solvers[i]->getGAP();
         nb_bnb_vector[i] = solvers[i]->getNbBnB();
         time_lb += solvers[i]->getTime();
         if(curr_status == Status::Time_Limit) ++nb_notopt;
      }
      if(solvers[i]->getNbSol() > 0) testProblem(i);
   }
   estimators();
   time_ = time(NULL) - start_time;

   if(input.getTransform() == Input::Transform::NonStrategic && best_problem > -1){
      ub_test = ub;
      var_test = var_ub;
      time_test = time_eval;
      nb_enter = nb_enter_aux;
      nb_improv = nb_improv_aux;
      std::copy(avg_rank.begin(),avg_rank.end(),avg_stu_per_rank.begin());
      
      computeTest(best_problem,ub_test_strat_simple,var_test_strat_simple,time_test_strat_simple,avg_stu_per_rank_strat_simple,nb_enter_strat_simple,nb_improv_strat_simple,Input::Transform::Simple);
      computeTest(best_problem,ub_test_strat,var_test_strat,time_test_strat,avg_stu_per_rank_strat,nb_enter_strat,nb_improv_strat,Input::Transform::Compact);
   }
   if(input.getTransform() == Input::Transform::Simple && best_problem > -1){
      computeTest(best_problem,ub_test_strat,var_test_strat,time_test_strat,avg_stu_per_rank_strat,nb_enter_strat,nb_improv_strat,Input::Transform::Compact);
      ub_test_strat_simple = ub;
      var_test_strat_simple = var_ub;
      time_test_strat_simple = time_eval;
      nb_enter_strat_simple = nb_enter_aux;
      nb_improv_strat_simple = nb_improv_aux;
      std::copy(avg_rank.begin(),avg_rank.end(),avg_stu_per_rank_strat_simple.begin());
   }
   if((input.getTransform() == Input::Transform::Compact || input.getTransform() == Input::Transform::Enumeration) && best_problem > -1){
      ub_test_strat = ub;
      var_test_strat = var_ub;
      time_test_strat = time_eval;
      nb_enter_strat = nb_enter_aux;
      nb_improv_strat = nb_improv_aux;
      std::copy(avg_rank.begin(),avg_rank.end(),avg_stu_per_rank_strat.begin());
   }
}

void SAA::estimators(){
   var_lb = 0.0;
   var_nb_bnb = 0.0;

   long long M = instance.getNbProblems();
   lb = std::accumulate(lb_vector.begin(),lb_vector.end(),0.0)/M;
   nb_bnb = std::accumulate(nb_bnb_vector.begin(),nb_bnb_vector.end(),0.0)/M;

   for(int i = 0; i < instance.getNbProblems(); ++i) {
      var_lb += std::pow((lb_vector[i]-lb),2);
      var_nb_bnb += std::pow((nb_bnb_vector[i]-nb_bnb),2);
   }
   if(M > 1) var_lb = var_lb/(M*(M-1));
   if(M > 1) var_nb_bnb = var_nb_bnb/(M*(M-1));
   var_gap = var_lb + var_ub;
   gap = (ub - lb)/ub;  

   // Statistical estimators considering 95% accuracy (alpha = 5 %)
   stat_gap = (ub - lb) + instance.getCriticalNormal()*std::sqrt(var_gap);             
   stat_lb = lb - instance.getCriticalTStudent()*std::sqrt(var_lb);     
   stat_ub = ub + instance.getCriticalNormal()*std::sqrt(var_ub);
}

void SAA::testProblem(int i){
   // Solving SAA test problem
   saatest->solve(solvers[i]->getT_(),input.getTransform());
   time_eval += saatest->getTime();

   if(saatest->getStatus() == Status::Optimal && saatest->getUB() < ub){
      ub = std::min(ub, saatest->getUB());
      for(int r = 0; r <= instance.getNbMaxRankSchools(); ++r) avg_rank[r] = saatest->getAvgRank(r);
      var_ub = 0.0;
      for(int s = 0; s < instance.getNbValidateScenarios(); ++s)
         var_ub += std::pow(saatest->getObjScenario(s)-ub,2);
      long long N = instance.getNbValidateScenarios();
      if(N > 1) var_ub = var_ub/(N*(N-1));
      best_problem = i;

      // Solve saatest with no extra capacity to compute nb entering and nb improving
      nb_enter_aux = 0.0;
      nb_improv_aux = 0.0;
      std::vector<std::vector<double>> rank_t_;
      rank_t_.resize(instance.getNbValidateScenarios());
      for(int s = 0; s < instance.getNbValidateScenarios(); ++s) {
         for(int st = 0; st < instance.getNbStudents(); ++st){
            rank_t_[s].push_back(saatest->getRankStudent(s,st));
         }
      }
      double * zero_capa = new double[instance.getNbSchools()];
      std::fill(zero_capa,zero_capa+instance.getNbSchools(),0.0);
      saatest->setVerbose(0);
      saatest->solve(zero_capa,input.getTransform());
      saatest->setVerbose(1);
      for(int s = 0; s < instance.getNbValidateScenarios(); ++s){
         for(int st = 0; st < instance.getNbStudents(); ++st){
            int unassigned = std::min(instance.getNbAssignFeasPairsByStudent(instance.getValidateProblem(),s,st),instance.getNbMaxRankSchools())+1;
            if((int)saatest->getRankStudent(s,st) >= unassigned){
               if((int)rank_t_[s][st] < unassigned)
                  nb_enter_aux += 1.0/instance.getNbValidateScenarios();

            }else{
               if(saatest->getRankStudent(s,st) > rank_t_[s][st]) 
                  nb_improv_aux += 1.0/instance.getNbValidateScenarios();
            }
         }
      }
   }
}
