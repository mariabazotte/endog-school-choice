#include "instance.hpp"

Instance::Instance(const Input & input): input(input), instance_file(input.getInstanceFile()), nb_max_rank_schools(input.getNbMaxRankSchools()),
                                        nb_saa_problems(input.getNbProblemsSAA()),nb_saa_scenarios(input.getNbScenariosSAA()),
                                        nb_saa_val_scenarios(input.getNbValidateScenariosSAA()),seed(input.getSeed()), 
                                        rd_engine(input.getSeed()),ext_val(location_,scale_),uniform(-5.0,5.0),normal(location_,scale_){
    defineName();
    defineCriticalValues(); 
    readInstance(); 
    std::cout << "Read instances finished..." << std::endl;
    double time_ = time(NULL);
    generateScenarios();
    std::cout << "Time to generate scenarios: " << time(NULL) - time_ << std::endl;
    std::cout << "Generate scenarios finished..." << std::endl;
    if(input.getTransform() == Input::Transform::Enumeration) {
        enumAllRanks(); 
        std::cout << "Rank enumeration finished..." << std::endl;
    }
}

void Instance::defineName(){
    name = instance_file;
    while(name.find("/") != std::string::npos)
        name = name.substr(name.find("/")+1);
    name.substr(0,name.size()-4);
} 

void Instance::defineCriticalValues(){
    critical_normal = 1.64;
    if(nb_saa_problems == 1) critical_tstudent = 6.31;
    if(nb_saa_problems == 5) critical_tstudent = 2.02;
    if(nb_saa_problems == 10) critical_tstudent = 1.81;
    else if(nb_saa_problems == 11) critical_tstudent = 1.80;
    else if(nb_saa_problems == 12) critical_tstudent = 1.78;
    else if(nb_saa_problems == 13) critical_tstudent = 1.77;
    else if(nb_saa_problems == 14) critical_tstudent = 1.76;
    else if(nb_saa_problems == 15 || nb_saa_problems == 16) critical_tstudent = 1.75;
    else if(nb_saa_problems == 17) critical_tstudent = 1.74;
    else if(nb_saa_problems == 18 || nb_saa_problems == 19) critical_tstudent = 1.73;
    else if(nb_saa_problems >= 20 && nb_saa_problems <= 22) critical_tstudent = 1.72;
    else if(nb_saa_problems >= 23 && nb_saa_problems <= 26) critical_tstudent = 1.71;
    else if(nb_saa_problems >= 27 && nb_saa_problems <= 39) critical_tstudent = 1.70;
    else if(nb_saa_problems >= 40 && nb_saa_problems <= 59) critical_tstudent = 1.68;
    else if(nb_saa_problems >= 60 && nb_saa_problems <= 79) critical_tstudent = 1.67;
    else if(nb_saa_problems >= 80 && nb_saa_problems <= 199) critical_tstudent = 1.66;
}

void Instance::readInstance(){
    std::ifstream file(instance_file.c_str());
    if(!file.fail()) { 
        file >> nb_schools;
        file >> nb_students;
        file >> nb_max_extra_capacity;

        int capa;
        std::string name;
        for(int c = 0; c < nb_schools; ++c){
            file >> name;
            file >> capa;
            id_schools.insert({name,c});
            schools_capacity.push_back(capa);
        }

        double val;
        int ist = 0;
        std::string school;
        std::string student;
        for(int i = 0; i < nb_schools*nb_students; ++i){
            file >> student;
            file >> school;

            if(id_students.count(student)==0){
                id_students.insert({student,ist});
                ++ist;
            }

            int c = id_schools[school];
            int st = id_students[student];
            all_feas_pairs.push_back({st,c});
            
            file >> val;
            scores.insert({{st,c},val});

            file >> val;
            utilities.insert({{st,c},val});

            file >> val;
            probab.insert({{st,c},val});

            inc_probab.resize(nb_max_extra_capacity+1);
            inc_probab[0].insert({{st,c},0.0});
            for(int j = 1; j <= nb_max_extra_capacity; ++j){
                file >> val;
                inc_probab[j].insert({{st,c},val});
            }
        }
        for(int st = 0; st < nb_students; ++st){
            scores.insert({{st,nb_schools},st+1.0}); 
        }
        if(ist != nb_students){
            std::cerr << "The number of students in the file does not correspond to the remaining information.\n";
        }
        file.close();
    }else{
        std::string errorMessage = std::string("Error: Could not open file ")+ instance_file + std::string("\n");
        throw std::runtime_error(errorMessage);
        exit(0);
    }
    if(input.getNbMaxExtraSeats() < nb_max_extra_capacity){
        nb_max_extra_capacity = input.getNbMaxExtraSeats();
    }
}

bool Instance::zeroProbability(const feas_pair & pair){
    bool always_zero = true;
    if(probab[pair] >= 1e-3) always_zero = false;
    for(int j = 0; j < nb_max_extra_capacity+1; ++j){
        if(inc_probab[j][pair] >= 1e-3) always_zero = false;
    }
    return always_zero;
}

void Instance::generateScenarios(){
    if(input.getType() == Input::Type::Deterministic){
        if(nb_saa_problems != 1) throw std::invalid_argument("Nb. SAA problems should be equal to 1 when solving the deterministic problem.");
        if(nb_saa_scenarios != 1)  throw std::invalid_argument("Nb. SAA scenarios should be equal to 1 when solving the deterministic problem.");
    }

    nb_feasible_pairs.resize(nb_saa_problems+1);
    feasible_pairs.resize(nb_saa_problems+1);
    random_utilities.resize(nb_saa_problems+1);
    feas_pair_by_school.resize(nb_saa_problems+1);
    feas_pair_by_student.resize(nb_saa_problems+1);
    not_assigned_pair.resize(nb_saa_problems+1);
    students_preference.resize(nb_saa_problems+1);
    schools_preference.resize(nb_saa_problems+1);

    int num_scenarios = nb_saa_scenarios;
    for(int pr = 0; pr < nb_saa_problems+1; ++pr){
        if(pr == nb_saa_problems) num_scenarios = nb_saa_val_scenarios;
        nb_feasible_pairs[pr].resize(num_scenarios,0);
        feasible_pairs[pr].resize(num_scenarios);
        random_utilities[pr].resize(num_scenarios);
        feas_pair_by_school[pr].resize(num_scenarios);
        feas_pair_by_student[pr].resize(num_scenarios);
        not_assigned_pair[pr].resize(num_scenarios);
        students_preference[pr].resize(num_scenarios);
        schools_preference[pr].resize(num_scenarios);

        for(int s = 0; s < num_scenarios; ++s) generateScenario(pr,s);
    }
}

void Instance::generateScenario(int pr, int s){
    std::vector<std::vector<feas_pair>> aux_feas_pairs;
    aux_feas_pairs.resize(nb_students);
    std::map<feas_pair,double> aux_rd_utilities;

    // Generate utilities of all possible pairs
    for(feas_pair pair: all_feas_pairs){ 
        double u = 0.0;
        generateUtility(pr,pair,u);
        // Keep only pairs with positive utility and probability
        if((u >= 1e-2)){ // && (zeroProbability(pair) == false)){
            aux_feas_pairs[pair.first].push_back(pair);
            aux_rd_utilities[pair] = u;
        }
    }

    feas_pair_by_school[pr][s].resize(nb_schools);
    feas_pair_by_student[pr][s].resize(nb_students);
    students_preference[pr][s].resize(nb_students);
    for(int st = 0; st < nb_students; ++st){
        
        // Preprocessing for strategic version
        if(input.getTransform() != Input::Transform::NonStrategic) preprocessing_strategic(aux_feas_pairs[st],aux_rd_utilities);
        
        // Sorting feasible pairs according to utility
        std::sort(aux_feas_pairs[st].begin(),aux_feas_pairs[st].end(), 
                    [&pr,&s,&aux_rd_utilities](feas_pair &a, feas_pair &b){ return aux_rd_utilities.at(a)>aux_rd_utilities.at(b); });

        // If this is not the evaluation problem (which we use for both strategic and non-strategic versions)
        // and if it is non-strategic, we must select only the k best schools
        int nb_pairs = aux_feas_pairs[st].size();
        if(pr != nb_saa_problems && input.getTransform() == Input::Transform::NonStrategic) 
            nb_pairs = std::min(nb_pairs,nb_max_rank_schools);

        for(int k = 0; k < nb_pairs; ++k){
            feas_pair pair = aux_feas_pairs[st][k];
            double u = aux_rd_utilities[pair];

            if(input.getTransform() != Input::Transform::NonStrategic){
                if(zeroProbability(pair) == true && ((int)feas_pair_by_student[pr][s][pair.first].size()) >= nb_max_rank_schools) {
                    continue; // will not choose schools with rank greater than k if their probability is zero (only choose schools with zero prob if rank < nb schools allowed)      
                }
            }

            feasible_pairs[pr][s].push_back(pair);
            random_utilities[pr][s].push_back(u);
            feas_pair_by_school[pr][s][pair.second].push_back(nb_feasible_pairs[pr][s]);
            feas_pair_by_student[pr][s][pair.first].push_back(nb_feasible_pairs[pr][s]);
            students_preference[pr][s][pair.first][nb_feasible_pairs[pr][s]] = feas_pair_by_student[pr][s][pair.first].size()-1;
            ++nb_feasible_pairs[pr][s];

            // Preprocessing according to probability
            if((input.getTransform() == Input::Transform::Compact || 
                input.getTransform() == Input::Transform::Enumeration) 
                && (probab[pair] >= 1.0)){
                if(((int)feas_pair_by_student[pr][s][pair.first].size()) >= nb_max_rank_schools) { break; } // will not choose schools with pref > K (max rank)
                else { nb_pairs = std::min((int)aux_feas_pairs[st].size(),nb_max_rank_schools); } // will select safety schools until size list = K (max rank)
            }
        }
    }
    
    // Unassigned option in the end of the list 
    for(int st = 0; st < nb_students; ++st){
        feasible_pairs[pr][s].push_back({st,nb_schools});
        random_utilities[pr][s].push_back(0.0);
        not_assigned_pair[pr][s].push_back(nb_feasible_pairs[pr][s]);
        students_preference[pr][s][st][nb_feasible_pairs[pr][s]] = feas_pair_by_student[pr][s][st].size();
        ++nb_feasible_pairs[pr][s];
    }

    schools_preference[pr][s].resize(nb_schools);
    for(int c = 0; c < nb_schools; ++c){
        std::sort(feas_pair_by_school[pr][s][c].begin(),feas_pair_by_school[pr][s][c].end(), 
                    [pr,s,this](int &a, int &b){ return getScore(pr,s,a) > getScore(pr,s,b); });
        if((int)feas_pair_by_school[pr][s][c].size() > 0){
            for(int n = 0; n < (int)feas_pair_by_school[pr][s][c].size(); ++n){
                int i = feas_pair_by_school[pr][s][c][n];
                schools_preference[pr][s][c][i] = n;
            }
        }
    }

    for(int st = 0; st < nb_students; ++st) aux_feas_pairs[st].clear();
    aux_feas_pairs.clear();
    aux_rd_utilities.clear();
}

void Instance::preprocessing_strategic(std::vector<feas_pair> & aux_feas_pairs,std::map<feas_pair,double> & aux_rd_utilities){
    std::vector<feas_pair> pairs_to_remove;

    for (const feas_pair& pair : aux_feas_pairs) {
        int nb_better = 0; 
        double best_pu = aux_rd_utilities[pair]*std::max(0.001,(probab[pair]+inc_probab[nb_max_extra_capacity][pair]));

        for (const feas_pair& comp_pair : aux_feas_pairs) {
            if(comp_pair.second != pair.second) {
                if((aux_rd_utilities[comp_pair]*std::max(0.001,probab[comp_pair])) >= best_pu + 1e-6) {   
                    if((input.getTransform() == Input::Transform::Simple) || 
                        ((input.getTransform() == Input::Transform::Compact || 
                        input.getTransform() == Input::Transform::Enumeration) 
                        && aux_rd_utilities[comp_pair] >= aux_rd_utilities[pair] + 1e-6)) {
                        ++nb_better;
                    }
                }
            }
        }
        if (nb_better >= nb_max_rank_schools) pairs_to_remove.push_back(pair);
    }

    for (const feas_pair& pair : pairs_to_remove) {
        auto it = std::find(aux_feas_pairs.begin(), aux_feas_pairs.end(), pair);
        if (it != aux_feas_pairs.end()) aux_feas_pairs.erase(it);
    }
    pairs_to_remove.clear();

    std::vector<feas_pair> pairs_to_choose;
    for(const feas_pair& pair : aux_feas_pairs){
        int nb_worst  = 0;
        double worst_pu = aux_rd_utilities[pair]*std::max(0.001,probab[pair]);

        for(const feas_pair& comp_pair : aux_feas_pairs) {
            if(comp_pair.second != pair.second) {
                if(worst_pu >= aux_rd_utilities[comp_pair]*std::max(0.001,(probab[comp_pair]+inc_probab[nb_max_extra_capacity][comp_pair])) + 1e-6){
                    if((input.getTransform() == Input::Transform::Simple) || 
                        ((input.getTransform() == Input::Transform::Compact || 
                        input.getTransform() == Input::Transform::Enumeration) 
                        && aux_rd_utilities[pair] >= aux_rd_utilities[comp_pair] + 1e-6)) {
                        ++nb_worst;
                    }
                }
            }
        }

        if(nb_worst >= ((int)aux_feas_pairs.size()-nb_max_rank_schools)) {
            pairs_to_choose.push_back(pair);
        }
    }

    if((int)pairs_to_choose.size() >= nb_max_rank_schools){
        for(const feas_pair& pair : aux_feas_pairs) {
            auto it = std::find(pairs_to_choose.begin(), pairs_to_choose.end(), pair);
            if(it == pairs_to_choose.end()) pairs_to_remove.push_back(pair);
        }
    }

    for (const feas_pair& pair : pairs_to_remove) {
        auto it = std::find(aux_feas_pairs.begin(), aux_feas_pairs.end(), pair);
        if (it != aux_feas_pairs.end()) aux_feas_pairs.erase(it);
    }
}

void Instance::generateUtility(int pr,feas_pair & pair, double & u){
    if(input.getType() == Input::Type::Stochastic || pr == nb_saa_problems){ 
        u = (utilities[pair] + (double)ext_val(rd_engine));
        // u = (utilities[pair] + (double)normal(rd_engine));
        // u = utilities[pair] + (double)uniform(rd_engine);
    }
    else if(input.getType() == Input::Type::Deterministic){
        u = (utilities[pair] + location_ + euler_cte*scale_);
        // u = utilities[pair] + location_;
        // u = utilities[pair];
    } 
    else u = 0.0;
}

/////////////////////////////////

void Instance::enumRank(int pr, int s, int st, std::vector<std::vector<int> >& ans, 
                        std::vector<int>& tmp, int n, int left, int k){
    if(k == 0){
        ans.push_back(tmp);
        return;
    }
    for (int i = left; i < n; ++i){
        tmp.push_back(feas_pair_by_student[pr][s][st][i]);
        enumRank(pr, s, st, ans, tmp, n, i + 1, k - 1);
        tmp.pop_back();
    }
}

void Instance::enumAllRanks(){
    enumerated_ranks.resize(nb_saa_problems+1);
    number_ranks.resize(nb_saa_problems+1);

    int num_scenarios = nb_saa_scenarios;
    for(int pr = 0; pr < nb_saa_problems+1; ++pr){
        if(pr == nb_saa_problems) num_scenarios = nb_saa_val_scenarios;
        enumerated_ranks[pr].resize(num_scenarios);
        number_ranks[pr].resize(num_scenarios);
        for(int s = 0; s < num_scenarios; ++s){
            enumerated_ranks[pr][s].resize(nb_students);
            number_ranks[pr][s].resize(nb_students);
            for(int st = 0; st < nb_students; ++st){
                int n = feas_pair_by_student[pr][s][st].size();
                int k = nb_max_rank_schools;
                if(n > k){
                    std::vector<std::vector<int>> ranks;
                    std::vector<int> tmp;
                    enumRank(pr, s, st, ranks, tmp, n, 0, k);
                    for(std::vector<int> r : ranks) enumerated_ranks[pr][s][st].push_back(r);
                    number_ranks[pr][s][st] = ranks.size();
                }
                else {
                    enumerated_ranks[pr][s][st].push_back(feas_pair_by_student[pr][s][st]); 
                    number_ranks[pr][s][st] = 1;
                }
            }
        }
    }
}

/////////////////////////////////

double Instance::mi_algorithm(int pr, int s, int st, const double * t_,
                                    std::map<int,int> & ranks) const{  
    std::vector<int> feas_pairs(feas_pair_by_student[pr][s][st]); 
    
    double val = 0.0;
    if((int)feas_pairs.size() <= nb_max_rank_schools && (int)feas_pairs.size() > 0){
        double p = 1.0;
        for(int pos_i = 0; pos_i < (int)feas_pairs.size(); ++pos_i){
            int i = feas_pairs[pos_i];
            int c = getSchoolFeasiblePair(pr,s,i);
            val += p*getProbability(pr,s,i,t_[c])*getUtility(pr,s,i);
            p *= (1.0-getProbability(pr,s,i,t_[c]));

            ranks.insert({pos_i,i});
        }
    }
    else if ((int)feas_pairs.size() > 0){
        while((int)ranks.size() < nb_max_rank_schools){
            val = 0.0;
            int fea = -1;
            int ran = -1;
            for(int pos_i = feas_pairs.size()-1; pos_i >= 0; --pos_i){
                double v = 0.0;
                double p = 1.0;
                
                int i = feas_pairs[pos_i];
                int r = getFeasPairStudentRank(pr,s,i);
                ranks.insert({r,i});
            
                for(auto const & [key, pair] : ranks){
                    int c = getSchoolFeasiblePair(pr,s,pair);
                    v += p*getProbability(pr,s,pair,t_[c])*getUtility(pr,s,pair);
                    p *= (1.0-getProbability(pr,s,pair,t_[c]));
                }

                if(v >= val){
                    fea = i;
                    ran = r;
                    val = v;
                }

                ranks.erase(r);
            }
            if(fea == -1) break; 

            ranks.insert({ran,fea});

            auto it = find(feas_pairs.begin(), feas_pairs.end(), fea);
            if(it == feas_pairs.end()) throw std::invalid_argument("The feasible pair should exist.");
            feas_pairs.erase(it);  
        }
    }
    return val;
}

double Instance::order_algorithm(int pr, int s, int st, const double * t_, 
                                    std::map<int,int> & ranks) const{
    
    std::vector<int> feas_pairs(feas_pair_by_student[pr][s][st]);
    std::map<int,double> exp_utilities;

    double u = 0.0;
    if((int)feas_pairs.size() <= nb_max_rank_schools && (int)feas_pairs.size() > 0){
        for(int pos_i = 0; pos_i < (int)feas_pairs.size(); ++pos_i){
            int i = feas_pairs[pos_i];
            int c = getSchoolFeasiblePair(pr,s,i);
            u += getProbability(pr,s,i,t_[c])*getUtility(pr,s,i);
            ranks.insert({pos_i,i});
        }
    }else if((int)feas_pairs.size() > 0){
        for(int i: feas_pairs) {
            int c = getSchoolFeasiblePair(pr,s,i);
            exp_utilities[i] = getProbability(pr,s,i,t_[c])*getUtility(pr,s,i);
        }

        std::sort(feas_pairs.begin(),feas_pairs.end(), 
            [&pr,&s,&exp_utilities](int &a, int &b){ 
                return exp_utilities.at(a) > exp_utilities.at(b); 
        });

        for(int k = 0; k < nb_max_rank_schools; ++k){
            int i = feas_pairs[k];
            u += exp_utilities[i];
            ranks.insert({getFeasPairStudentRank(pr,s,i),i});
        }
    }
    return u;
}

double Instance::da_algorithm(int pr, int s, const double* t_,
                                    std::vector<std::list<int>> feas_pairs, 
                                    std::unordered_map<int,double> rank_feas_pairs,
                                    std::vector<std::list<int>> & matching) const{

    // List of students not in the matching (we can still try to match them)
    std::vector<int> not_matched;

    for(int st = 0; st < nb_students; ++st){
        if(feas_pairs[st].empty() == false) not_matched.push_back(st);
        else matching[nb_schools].push_back(getUnassignedPair(pr,s,st));  
    }

    // Deffered acceptance algorithm
    while(!not_matched.empty()){

        for(auto it = not_matched.begin(); it != not_matched.end(); it++){
            int st = *it;

            int i = feas_pairs[st].front();
            feas_pairs[st].pop_front();

            int c = getSchoolFeasiblePair(pr,s,i);
            matching[c].push_back(i);
        }
        not_matched.clear();

        for(int c = 0; c < nb_schools; ++c){
            if((int) matching[c].size() > (schools_capacity[c] + t_[c])) 
                matching[c].sort([s,pr,this](int &a, int &b){ 
                                return getScore(pr,s,a) > getScore(pr,s,b); });
            
            while((int)matching[c].size() > (schools_capacity[c] + t_[c])){
                int i = matching[c].back();
                matching[c].pop_back();

                int st = getStudentFeasiblePair(pr,s,i);
                if(feas_pairs[st].empty() == false) not_matched.push_back(st); 
                else matching[nb_schools].push_back(getUnassignedPair(pr,s,st));   
            }
        }
    }

    double val = 0.0;
    for(int c = 0; c < nb_schools + 1; ++c){
        for(int i: matching[c]) val += (double)rank_feas_pairs[i];
    }
    return val;
}

/////////////////////////////////

void Instance::display(){
    std::cout << "----------- INSTANCE PARAMS -------------" << std::endl;
    std::cout << "NB STUDENTS          : " << nb_students << std::endl;
    std::cout << "NB SCHOOLS           : " << nb_schools << std::endl;
    std::cout << "NB MAX EXTRA CAPACITY: " << nb_max_extra_capacity << std::endl;
    std::cout << "NB MAX RANK          : " << nb_max_rank_schools << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
}