#include "input.hpp"

void Input::defaultParams(){
    // General parameters
    time_limit = 3600;
    nb_threads = 1;
    seed = 0;
    verbose = 0;

    // Instance parameters
    nb_max_rank_schools = 4; 
    nb_max_extra_seats = 5;   
    obj_access_improv = true;                    

    // Parameters for SAA or evaluation analysis
    nbproblemsSAA = 5;
    nbscenariosSAA = 100;
    nbvalidatescenariosSAA = 5000;

    // For all formulations
    warm_start = true;

    // For exact formulations: Cutoff Score and L-Constraints
    lazy_stability = false;
    combs = false;
    branch_tooth = false;
    linear_relax = false;

    // For Lagrangian
    max_nb_it_lag = 300;
    max_nb_it_without_improv_lag = 15;
    init_scale_lag = 2.0;
    min_step_lag = 1e-4;
    modified_sub = false;
    init_mult_lp = false;

    // For Integer LShaped
    lshaped_callback = true;
    lshaped_alternating = false;

    // For strategic version
    barrier = false;
    resprep = false;
    largenumstab = false;
}

Input::Input(int argc, char* argv[]){
    defaultParams();
    int mandatory = 0;
    for (int i = 1; i < argc; i += 2){
        if(std::string(argv[i]) == "-instancefile"){ // Mandatory parameters
            instance_file = argv[i+1];
            mandatory += 1;
        }else if(std::string(argv[i]) == "-solver"){
            int_solver = std::stoi(argv[i+1]);
            solver = Solver(int_solver);
            mandatory += 1;
        }else if(std::string(argv[i]) == "-transform"){
            int_transform = std::stoi(argv[i+1]);
            transform = Transform(int_transform);
            mandatory += 1;
        }else if(std::string(argv[i]) == "-type"){
            int_type = std::stoi(argv[i+1]);
            type = Type(int_type);
            mandatory += 1;
        }
        else if(std::string(argv[i]) == "-nbranks"){
            nb_max_rank_schools = std::stoi(argv[i+1]);  
            mandatory += 1;
        }
        else if(std::string(argv[i]) == "-nbseats"){
            nb_max_extra_seats = std::stoi(argv[i+1]);  
            mandatory += 1;
        }
        else if(std::string(argv[i]) == "-objtype")  // Optional :: Instance
            obj_access_improv = std::stoi(argv[i+1]);  
        else if(std::string(argv[i]) == "-timelimit") // Optional :: General parameters
            time_limit = std::stod(argv[i+1]);   
        else if(std::string(argv[i]) == "-nbthreads") 
            nb_threads = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-seed") 
            seed = std::stod(argv[i+1]);
        else if(std::string(argv[i]) == "-verbose")
            verbose = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbproblemsSAA") // Optional :: SAA problem and validation
            nbproblemsSAA = std::stoi(argv[i+1]); 
        else if(std::string(argv[i]) == "-nbscenariosSAA")
            nbscenariosSAA = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-nbvalidatescenariosSAA")
            nbvalidatescenariosSAA = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-warmstart") // Optional :: improving formulations 
            warm_start = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-lazy")  
            lazy_stability = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-branchtooth")
            branch_tooth = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-combs")
            combs = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-linrelax")
            linear_relax = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-maxitlag") // Optional :: Lagrangian
            max_nb_it_lag = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-maxitimplag")
            max_nb_it_without_improv_lag = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-modsub")
            modified_sub = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-inilp")
            init_mult_lp = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-lshapedcall") // Optional :: Integer LShaped
            lshaped_callback = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-lshapedalt") // Optional :: Integer LShaped
            lshaped_alternating = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-barrier") // Optional :: Strategy
            barrier = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-resprep") 
            resprep = std::stoi(argv[i+1]);
        else if(std::string(argv[i]) == "-largenumstab") 
            largenumstab = std::stoi(argv[i+1]);
        else{
            std::cerr << "ERROR: Argument '" << argv[i] << "' not defined." << std::endl;
            throw std::runtime_error(std::string("Incorrect line of command"));
        }
    }
    testParameters();
    if(mandatory < 6){
        std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
        if(transform == Transform::NonStrategic) 
            std::cerr << "You need to define: -instancefile -solver -transform -type -nbranks -nbseats" << std::endl;
        else 
            std::cerr << "You need to define: -instancefile -solver -transform -type -nbranks -nbseats -objstrat" << std::endl;
        throw std::runtime_error(std::string("Incorrect line of command"));
    }
    defineSolutionFile();
}

void Input::defineSolutionFile(){
    std::string instfile = instance_file;

    while(instfile.find("/") != std::string::npos)
        instfile = instfile.substr(instfile.find("/")+1);
    
    solution_file = "../results/";
    
    if(type == Input::Type::Deterministic) 
        solution_file += "deterministic/";
    if(type == Input::Type::Stochastic) 
        solution_file += "stochastic/";
    
    if(transform == Input::Transform::NonStrategic)
        solution_file += "nonstrategic/";
    if(transform == Input::Transform::Compact)
        solution_file += "compact/";
    if(transform == Input::Transform::Enumeration)
        solution_file += "enumeration/";
    if(transform == Input::Transform::Simple)
        solution_file += "simple/";
    
    if(solver == Input::Solver::LConstraints)
        solution_file += "lconstraints/";
    if(solver == Input::Solver::ModLConstraints)
        solution_file += "modlconstraints/";
    if(solver == Input::Solver::CutoffScore)
        solution_file += "cutoffscore/";
    if(solver == Input::Solver::ASGHeuristic)
        solution_file += "asgheuristic/";
    if(solver == Input::Solver::Lagrangian)
        solution_file += "lagrangian/";
    if(solver == Input::Solver::LocalSearch)
        solution_file += "localsearch/";
    if(solver == Input::Solver::IntegerLShaped)
        solution_file += "integerlshaped/";
    if(solver == Input::Solver::SimulAnnealing)
        solution_file += "simulannealing/";

    // Configuration corresponding to instance
    int seed_aux = seed;
    solution_file += instfile.substr(0, instfile.size()-4) + "_maxseats" + std::to_string(nb_max_extra_seats) + "_maxrank" + std::to_string(nb_max_rank_schools) + "_obj" + std::to_string(obj_access_improv) + "_seed" + std::to_string(seed_aux);
    
    // Configuration corresponding to type
    if(type == Input::Type::Deterministic)
        solution_file += "_nbval" + std::to_string(nbvalidatescenariosSAA);
    if(type == Input::Type::Stochastic)
        solution_file += "_nbpr" + std::to_string(nbproblemsSAA) + "_nbsc" + std::to_string(nbscenariosSAA) + "_nbval" + std::to_string(nbvalidatescenariosSAA); 
    
    // Configuration corresponding to solver
    if(solver != Input::Solver::LocalSearch && solver != Input::Solver::SimulAnnealing)
        solution_file += "_ws" + std::to_string(warm_start);
    if(solver == Input::Solver::LConstraints || solver == Input::Solver::ModLConstraints || solver == Input::Solver::CutoffScore)
        solution_file += "_ls" + std::to_string(lazy_stability) + "_bt" + std::to_string(branch_tooth) + "_comb" + std::to_string(combs) + "_rl" + std::to_string(linear_relax);
    
    if(solver == Input::Solver::Lagrangian)
        solution_file += "_it" + std::to_string(max_nb_it_lag) + "_itimp" + std::to_string(max_nb_it_without_improv_lag) + "_md" + std::to_string(modified_sub) + "_inilp" + std::to_string(init_mult_lp);

    if(solver == Input::Solver::IntegerLShaped){
        solution_file += "_call" + std::to_string(lshaped_callback);
        solution_file += "_alt" + std::to_string(lshaped_alternating);
    }
    
    if(solver != Input::Solver::LocalSearch && solver != Input::Solver::IntegerLShaped && solver != Input::Solver::Lagrangian && solver != Input::Solver::SimulAnnealing){
        if(transform != Input::Transform::NonStrategic && barrier == true)
            solution_file += "_barr";
        if(transform != Input::Transform::NonStrategic && resprep == true)
            solution_file += "_consprep";
        if(transform != Input::Transform::NonStrategic && largenumstab == true)
            solution_file += "_largenumstab";
    }

    int aux_time = time_limit;
    solution_file +=  "_time" + std::to_string(aux_time);
    solution_file += ".csv";
        
    solFile.open(solution_file, std::ios_base::out);
}

void Input::testParameters(){
    if(int_solver >= 8){
        throw std::runtime_error(std::string("Wrong solver value."));
        exit(0);
    }
    if(int_transform >= 4){
        throw std::runtime_error(std::string("Wrong transform value."));
        exit(0);
    }
    if(int_type >= 2){
        throw std::runtime_error(std::string("Wrong type value."));
        exit(0);
    }

    // If deterministic we generate only one problem, with one scenario. 
    // For testing, we have possibly more scenarios.
    if(type == Input::Type::Deterministic){ 
        nbproblemsSAA = 1;
        nbscenariosSAA = 1;
    }

    // Combs are implemented only for non-strategic version.
    if(transform != Input::Transform::NonStrategic){
        if(combs == true){
            std::cout << "Combs are not defined for Strategic version. Setting combs to false." << std::endl;
            combs = false;
        }
    }
    
    // ASGHeuristic, Lagrangian, Local Search and Integer L-Shaped have no lazy stability constraints, priority to branch-on-tooth, combs or linear relaxation options.
    if(solver == Input::Solver::ASGHeuristic || solver == Input::Solver::Lagrangian || solver == Input::Solver::LocalSearch || solver == Input::Solver::IntegerLShaped || solver == Input::Solver::SimulAnnealing){
        if(lazy_stability == true){
            std::cout << "Lazy stability constraints are not defined for " << solver << ". Setting lazy stability constraints to false." << std::endl;
            lazy_stability = false;
        }
        if(branch_tooth == true){
            std::cout << "Priority to branch on \"tooth\" is not defined for " << solver << ". Setting branch on tooth to false." << std::endl;
            branch_tooth = false;
        }
        if(combs == true){
            std::cout << "Combs are not defined for " << solver << ". Setting combs to false." << std::endl;
            combs = false;
        }
        if(linear_relax == true){
            std::cout << "Linear relaxation option is not implemented for " << solver << ". Setting linear relaxation to false." << std::endl;
            linear_relax = false;
        }
    }

    // Warm start not defined for local search.
    if(solver == Input::Solver::LocalSearch || solver == Input::Solver::SimulAnnealing){
        if(warm_start == true){
            std::cout << "Warm start option is not implemented for " << solver << ". Setting warm start to false." << std::endl;
            warm_start = false;
        }
    }

    // Lagrangian only for non-strategic version.
    if(transform != Input::Transform::NonStrategic && solver == Input::Solver::Lagrangian){
        std::cout << solver << " is not defined for Strategic versions. Setting transform to nonstrategic (option 0)." << std::endl;
        transform = Input::Transform::NonStrategic;
        int_transform = 0;
    }

}

void Input::display(){
    std::cout << "----------- INPUT PARAMS -------------" << std::endl;
	std::cout << "INSTANCE  FILE: '" << instance_file << "'" << std::endl;
    std::cout << "NB MAX RANK   :  " << nb_max_rank_schools << std::endl;
    std::cout << "NB MAX SEAT   :  " << nb_max_extra_seats << std::endl;
    std::cout << "ACCESS/IMPROV :  " << obj_access_improv << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "SOLUTION  FILE: '" << solution_file << "'" << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "TYPE          :  " << type << std::endl;
    std::cout << "TRANSFORM     :  " << transform << std::endl;
    std::cout << "SOLVER        :  " << solver << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "SEED          :  " << seed << std::endl;
    std::cout << "TIME LIMIT    :  " << time_limit << std::endl;
    std::cout << "NUMBER THREADS:  " << nb_threads << std::endl;
    std::cout << "VERBOSE       :  " << verbose << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "NB PROBLEMS   :  " << nbproblemsSAA << std::endl;
    std::cout << "NB SCENARIOS  :  " << nbscenariosSAA << std::endl;
    std::cout << "NB V SCENARIOS:  " << nbvalidatescenariosSAA << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "WARM START    :  " << warm_start << std::endl;
    std::cout << "LAZY          :  " << lazy_stability << std::endl;
    std::cout << "BRANCH TOOTH  :  " << branch_tooth << std::endl;
    std::cout << "COMBS         :  " << combs << std::endl;
    std::cout << "LINEAR RELAX  :  " << linear_relax << std::endl;
    std::cout << "--------------------------------------" << std::endl;
    if(solver == Input::Solver::Lagrangian){
        std::cout << "NB MAX IT     :  " << max_nb_it_lag << std::endl;
        std::cout << "NB MAX IT IMP :  " << max_nb_it_without_improv_lag << std::endl;
        std::cout << "INIT MULT LP  :  " << init_mult_lp << std::endl;
        std::cout << "MOD SUBGRAD   :  " << modified_sub << std::endl;
        std::cout << "--------------------------------------" << std::endl;
    }
    if(solver == Input::Solver::IntegerLShaped){
        std::cout << "CALLBACK      :  " << lshaped_callback << std::endl;
        std::cout << "ALTERNATING   :  " << lshaped_alternating << std::endl;
        std::cout << "--------------------------------------" << std::endl;
    }
    if(transform != Input::Transform::NonStrategic){
        std::cout << "BARRIER       :  " << barrier << std::endl;
        std::cout << "CONSERV. PREP.:  " << resprep << std::endl;
        std::cout << "LARGE NUM STAB:  " << largenumstab << std::endl;
        std::cout << "--------------------------------------" << std::endl;
    }
}

void Input::write(std::string output){
    if(!solFile.fail()){
        solFile << type << ";";
        solFile << transform << ";";
        solFile << solver << ";";
        solFile << instance_file << ";";
        solFile << nb_max_rank_schools << ";";
        solFile << nb_max_extra_seats << ";";
        if(obj_access_improv == false) solFile << "access;";
        else solFile << "improv;";
        solFile << seed << ";";
        solFile << nbproblemsSAA << ";";
        solFile << nbscenariosSAA << ";";
        solFile << nbvalidatescenariosSAA << ";";
        solFile << output;
        solFile << std::endl;
    }
}

void Input::writeHead(std::string output){
    if(!solFile.fail()){
        solFile << "TYPE;TRANSFORM;SOLVER;INSTANCE_FILE;MAX_RANK;MAX_SEAT;OBJ_TYPE;SEED;";
        solFile << "NB_PROBLEMS;NB_SCENARIOS;NB_TEST_SCENARIOS;";
        solFile << output;
        solFile << std::endl;
    }
}

std::ostream& operator<<(std::ostream& lhs, const Input::Solver & solver) {
    switch(solver) {
        case Input::Solver::LConstraints: {
            lhs << "L-Constraints";
            break;
        }
        case Input::Solver::ModLConstraints: {
            lhs << "Modified L-Constraints";
            break;
        }
        case Input::Solver::CutoffScore: {
            lhs << "Cutoff Score";
            break;
        }
        case Input::Solver::ASGHeuristic: {
            lhs << "ASG heuristic";
            break;
        }
        case Input::Solver::Lagrangian: {
            lhs << "Lagrangian";
            break;
        }
        case Input::Solver::LocalSearch: {
            lhs << "Local Search";
            break;
        }
        case Input::Solver::IntegerLShaped: {
            lhs << "Integer L-Shaped";
            break;
        }
        case Input::Solver::SimulAnnealing: {
            lhs << "Simulated Annealing";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::Transform & transform) {
    switch(transform) {
        case Input::Transform::NonStrategic: {
            lhs << "Non Strategic";
            break;
        }
        case Input::Transform::Compact: {
            lhs << "Strategic Compact";
            break;
        }
        case Input::Transform::Enumeration: {
            lhs << "Strategic Enumeration";
            break;
        }
        case Input::Transform::Simple: {
            lhs << "Strategic Simple";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}

std::ostream& operator<<(std::ostream& lhs, const Input::Type & type) {
    switch(type) {
        case Input::Type::Deterministic: {
            lhs << "Deterministic";
            break;
        }
        case Input::Type::Stochastic: {
            lhs << "Stochastic";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}