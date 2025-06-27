# endog-school-choice

## Description

This repository contains all instances, results, and code of the paper named "Capacity Planning in Stable Matching with Truthful
or Strategic Preference Uncertainty", from Maria Bazotte, Margarida Carvalho, and Thibaut Vidal.

## Building

### Instances

All instances used in the paper are available in the folder `data/instances`. Moreover, the code to generate new instances is available in the folder `data/generation`, as well as an auxiliary script `instances.sh`. To generate a random instance, execute:
```
cd data/generator/
python rd_instance_generator.py -parameters
```

The following parameters are available:
- -size_boots  'T'  # Number of bootstrapping samples used to generate the instance
- -nbschools   'C'  # Number of schools
- -nbstudents  'S'  # Number of students
- -nbextracapa 'B'  # Number of scholarships (extra capacity)
- -seed        'E'  # Seed for instance generation

The resulting instance file, named `'T'bts_'C'nbschools'S'nbstudents'B'extracapa'E'seed.txt`, is saved in the folder `data/instances`.

### Main Code

The main code is available in the folder `src/`. To run the main code, execute:
```
cd src/
make
./exe -parameters
``` 
The following parameters are available:

#### Configuration Parameters (Mandatory)
- -solver            # Solution method (programmed usng Gurobi): 0-> L-Constraints, 1-> Modified L-Constraints, 2-> Cutoff Score, 3-> Assignment (ASG) Heuristic, 4-> Lagrangian (LR) Heuristic, 5-> Local Search (LS) Heuristic, 6-> Integer L-Shaped, 7-> Simulated Annealing (SA) Heuristic
- -transform         # Student behavior model: 0-> Nonstrategic Utility Maximization (UM), 1-> Compact Strategic Conjoint Expected Utility Maximization (CEUM), 2-> Enumeration Strategic CEUM (enumerates all possible ranks for CEUM behavior), 3-> Simple Strategic Expected Utility Maximization (IEUM)
- -type              # Uncertainty modeling: 0-> Deterministic Problem (average scenario/ EV problem), 1-> Stochastic Problem (SAA method)

#### Instance Parameters (Mandatory)
- -instancefile      # Path to instance file, e.g., `"../data/instances/20bts_6schools50students100extracapa0seed.txt"`
- -nbranks           # Number maximum of schools each student can list (K)
- -nbseats           # Number of maximum extra seats or budget (B)

#### General Parameters
- -timelimit         # Time limit for each instance of the SAA program
- -nbthreads         # Number of threads for Gurobi solver 
- -seed              # Seed for scenario generation
- -verbose           # Level of print: 0 -> No print, 1-> Moderate print, 2-> Complete print
- -warmstart         # 0-> Do not use warm start procedure, 1-> Use the procedure

#### SAA and Evaluation Parameters
- -nbproblemsSAA            # Number of SAA programs (Set to one for paper)
- -nbscenariosSAA           # Number of scenarios for each SAA program (N)
- -nbvalidatescenariosSAA   # Number of scenarios for the evaluation of fixed first stage decisions (not restricted to SAA method) (N')

#### Mathematical Formulation Parameters
- -lazy               # 0-> Do not define  stability constraints as lazy constraints, 1-> Define them as lazy contraints
- -linrelax           # 0-> Do not solve the linear relaxation, 1-> Solve the linear relaxation

#### Lagrangian Parameters
- -maxitlag            # Maximum number of iterations 
- -maxitimplag         # Maximum number of iterations without improvement
- -modsub              # 0-> Do not use modified version for direction, 1-> Use modified version with direction (Set to zero for paper)
- -inilp               # 0-> Initialize multipliers with zero vector, 1-> Initialize multipliers with dual of the lp relaxation (Set to zero for paper)

#### L-Shaped Parameters
- -lshapedcall         # 0-> Do not use callback for cutting plane, 1-> Use callback
- -lshapedalt          # 0-> Do not include alternating cuts, 1-> Include alternating cuts

#### Folders Main Code

The main code contains the following folders:

- -`src/input`     # This folder contains the class Input, which stores information on input parameters provided by the user to the code during execution.
- -`src/instance`  # This folder contains the class Instance, which has general information about an instance, including data from schools and students. It also generates all scenarios corresponding to students' preferences. 
- -`src/algorithm` # This folder contains the classes Algorithm and SAA, which define the solution procedures to solve the expected/average scenario (deterministic) and stochastic versions of the problem.
- -`src/solver`    # This folder has the coding of the solution methods proposed in the paper. Namely, the L-constraints, modified L-constraints, and cutoff score formulation programmed using Gurobi, the ASG, LR, LS, and SA heuristics, and the Integer L-Shaped. It also contains the Test1stSolver class, which uses the implementation of the Marginal Improvement (MI), Individual Ordering (IO), and Deferred Acceptance (DA) algorithms to compute distinct metric evaluations for fixed first-stage capacity allocation decisions, considering a large number of scenarios.
- -`src/transform` # This folder has the mathematical formulations of the transformation functions for each student behavior. Namely, the nonstrategic (UM), compact (CEUM), and simple (IEUM) behaviors.

## Results

The results of the paper are saved in a repository such as: 

- results
    - compiled 
        - graphics
            - det_stoch
            - perf
            - strat
    - deterministic 
        - compact
            - simulannealing
        - nonstrategic
            - modlconstraints
        - simple 
            - simulannealing
    - stochastic 
        - compact
            - cutoffscore
            - localsearch
            - simulannealing
        - nonstrategic
            - asgheuristic
            - cutoffscore
            - lagrangian
            - localsearch
            - modlconstraints
            - simulannealing
        - simple 
            - cutoffscore
            - localsearch
            - simulannealing

To create these folders before running the main code, execute: 
```
python folders.py
```

Each execution of the main code generates a different `.csv` file in the folder corresponding to the modeling of uncertainty (deterministic for the EV problem, and stochastic for the SAA method), the modeling of student behavior (nonstrategic for UM behavior, compact for strategic CEUM behavior, and simple for strategic IEUM behavior), and the solution procedure (cutoffscore, modlconstraints, asgheuristic, lagrangian, localsearch, and simulannealing). The folder `results/compiled` contains two Excel files that compile the experimental results for the paper. The file `results_paper_final.xlsx` corresponds to the results for large instances in the paper, and the file `results_paper_small.xlsx` for small instances.

## Replicating

The folder `scripts/` contains the sbatch files for replicating the results of the paper. Specifically, the files `jobs_det.sh`, `job_nonstrat_comp_perf.sh`, `jobs_strat_comp_perf_large.sh`, and `jobs_strat_comp_perf_small.sh` contain the commands for the results in the paper. After executing these sbatch files, the resulting `.csv` files are written in the folder `results`, according to the configuration explained before.
        
### Unifying Paper Results

The file `unify_results.py` unifies the results contained in the `.csv` files, creating the Excel files with the compiled results of the paper. At line 24 of this file, the parameter `small_instances` is defined. Set small_instances = False for compiling results on large instances (to generate file `results_paper_final.xlsx`), and set small_instances = True for compiling results on small instances (to generate file `results_paper_small.xlsx`). After running the sbatch files and defining this parameter, execute:
```
python unify_results.py
```
The instance and configuration parameters in this file can also be modified to test different experiments. As previously mentioned, the unified results are written in the folder `results/compiled`.

### Generating Graphics

The file `graphics.py` generates all the graphics in the paper according to the unified results on the excel files. They are written in the folder `results/compiled/graphics`. Folder `results/compiled/graphics/det_stoch` contains the results related to comparing the deterministic and stochastic versions (e.g., VSS), folder `results/compiled/graphics/perf` contains the results related to the computational performance of the proposed solution procedures, and folder `results/compiled/graphics/strat` contains the results related to comparing different student behaviors,
