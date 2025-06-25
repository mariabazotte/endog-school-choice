import pandas as pd
import openpyxl
import argparse
import numpy as np
import math
from scipy import stats

############################################################################################################
#                                          Folders with files
############################################################################################################
general_folder = "../results/"
folder_type    = ["deterministic/","stochastic/"]
folder_transf  = ["nonstrategic/","compact/","enumeration/","simple/"]
folder_solver  = ["lconstraints/","modlconstraints/","cutoffscore/","asgheuristic/",\
                  "lagrangian/","localsearch/","integerlshaped/","simulannealing/"]
name_solver    = ["lconstraints","modlconstraints","cutoffscore","asgheuristic",\
                  "lagrangian","localsearch","integerlshaped","simulannealing"]

############################################################################################################
#                                              Instances
############################################################################################################
# You can modify these parameters if trying to test a different set of instances.

small_instances = False # If compile results for small instances or large instances

path = "../results/compiled/results_paper_final.xlsx"
seeds      = [0,1000,2000,3000,4000]
nbschools  = [20,40,60]
nbstudents = [500,1000]
nbseats    = [1,5,15,30,45,60]
nbranks    = [2,3,4]
objectives = [1]
maxnbseats = 100

if small_instances == True:
    path = "../results/compiled/results_paper_small.xlsx"
    nbschools  = [6]
    nbstudents = [50,100]
    nbseats    = [1,3,5,10,15]

instances = [[school,student,seat,rank,obj] for school in nbschools \
                                                for student in nbstudents \
                                                for seat in nbseats \
                                                for rank in nbranks \
                                                for obj in objectives]

############################################################################################################
#                                              Approches 
############################################################################################################
# types:      0-> deterministic, 1-> stochastic
# transforms: 0-> nonstrategic (UM), 1-> complex (CEUM), 3-> simple (IEUM)
# solvers:    0-> l-constraints, 1-> modified l-constraints, 2-> cutoff score, 3-> ASG heuristic, 
#             4-> lagrangian, 5-> local search, 6-> integer lshaped, 7-> simulated annealing
############################################################################################################
# These are the types and transforms tested in the paper.

types      = [0,1]                                
transforms = [0,1,3]   

if small_instances == True:
    types      = [1]
    transforms = [1,3]

############################################################################################################
#                            Parameters for comparing computational performance
############################################################################################################
# We compare the computational performance of the different methods to solve the problem

############## Solvers to compare performance ##############  
# You can modify these parameters if trying to test a different set of approaches.
# We evaluate different approaches considering type and transform:                   
solvers = {
    0:{ 0:[1], 1:[5,7], 3:[5,7] },
    1:{ 0:[1,3,4,5,7], 1:[5,7], 3:[5,7] } 
}   

if small_instances == True:
    solvers = {
        0: { 0:[], 1:[], 3:[]},
        1:{ 1:[2,5,7], 3:[2,5,7] } 
    }  

models = [(tp,tr,sol) for tp in types for tr in transforms for sol in solvers[tp][tr]]

########## Lagrangian ##########
# We test different parameters for the lagrangian: (Each configuration is one different approach to compare results)
maxnbit  = [300]      # Maximum number of iterations
maxitimp = [15]       # Maximum number of iterations without improvement
modified = [0]        # 0-> do not use modified version with direction. 1-> use modified version with direction
inilp    = [0]        # 0-> initialize multipliers with zero vector. 2-> initialize multipliers with dual of the lp relaxation.

conf_lag = [(it,itimp,mod,lp) for it in maxnbit for itimp in maxitimp for mod in modified for lp in inilp]
name_lag = {(it,itimp,mod,lp): "_maxit" + str(it) + "_maxitimprov" + str(itimp) + "_mod" + str(mod) + "_inilp" + str(lp)\
                              for it in maxnbit for itimp in maxitimp for mod in modified for lp in inilp}

########## Config. Scenarios According to the Paper ##########
# 0-> deterministic 1-> stochastic
# 0-> nonstrategic, 1-> compact, 3-> simple
conf_scenarios = {
    0:{tr: "nbval5000" for tr in transforms},\
    1:{tr: "nbpr1_nbsc100_nbval5000" if tr == 0 else "nbpr1_nbsc100_nbval5000" for tr in transforms}
}

if small_instances == True:
    conf_scenarios = {
        1:{tr: "nbpr1_nbsc50_nbval5000" for tr in transforms}
    }

########## Config. Time According the Paper ##########
# 0-> deterministic 1-> stochastic
# 0-> nonstrategic, 1-> compact, 3-> simple
time_perf = {
    0:{0: 18000, 1: 18000, 3: 18000}, 
    1:{0: 7200, 1: 28800, 3: 28800}
} 

########## Name file for general solvers config. used in the paper ##########
def get_name_solfile_perf(instance_file,tp,tr,sol,cf_lag):
    name_file = general_folder + folder_type[tp] + folder_transf[tr] + folder_solver[sol] + instance_file + conf_scenarios[tp][tr]
    if sol != 5 and sol != 7:
        name_file += "_ws1"
    if sol == 0 or sol == 1 or sol == 2:
        name_file += "_ls0_bt0_comb0_rl0" 
    if sol == 4:
        name_file += "_it" + str(cf_lag[0]) + "_itimp" + str(cf_lag[1]) + "_md" + str(cf_lag[2]) + "_inilp" + str(cf_lag[3]) 
    if sol == 6:
        name_file += "_call1_alt1"   
    name_file += "_time" + str(time_perf[tp][tr]) + ".csv" 
    return name_file                

############################################################################################################
#                         Parameters for VSS and Comparison between strategies/transforms
############################################################################################################
# Use cutoff score for all versions, with longer time (one day)
# Choose a solver for each type and transform to compute vss and comparison between strategies
vss_solver = {
    0:{tr:1 if tr == 0 else 7 for tr in transforms},
    1:{tr:4 if tr == 0 else 7 for tr in transforms}
}
vss_conf_sce = {
    0:{tr: "nbval5000" for tr in transforms},\
    1:{tr: "nbpr1_nbsc100_nbval5000" if tr == 0 else "nbpr1_nbsc100_nbval5000" for tr in transforms}
}
vss_conf_sol = {
    0:{tr:"_ws1_ls0_bt0_comb0_rl0_time18000" if tr == 0 else "_time18000" for tr in transforms},
    1:{tr:"_ws1_it300_itimp15_md0_inilp0_time7200" if tr == 0 else "_time28800" for tr in transforms}
}

# _ws1_ls0_bt0_comb0_rl0_time86400   _ws1_ls0_bt0_comb0_rl0_time86400
def get_name_solfile_vss(instance_file,tp,tr):
    return general_folder + folder_type[tp] + folder_transf[tr] + folder_solver[vss_solver[tp][tr]] + instance_file + vss_conf_sce[tp][tr] + vss_conf_sol[tp][tr] + ".csv"
       
############################################################################################################
#                                               Dictionaries 
############################################################################################################

####################################### Initialize dict with instances #####################################

def init_seed_dict(seed,instance):
    for tr in transforms:
        dict_det[tr]["seed"].append(seed)
        dict_det[tr]["schools"].append(instance[0])
        dict_det[tr]["students"].append(instance[1])
        dict_det[tr]["maxseats"].append(instance[2])
        dict_det[tr]["maxranks"].append(instance[3])
        
        dict_stoch[tr]["seed"].append(seed)
        dict_stoch[tr]["schools"].append(instance[0])
        dict_stoch[tr]["students"].append(instance[1])
        dict_stoch[tr]["maxseats"].append(instance[2])
        dict_stoch[tr]["maxranks"].append(instance[3])
        
        dict_vss[tr]["seed"].append(seed)
        dict_vss[tr]["schools"].append(instance[0])
        dict_vss[tr]["students"].append(instance[1])
        dict_vss[tr]["maxseats"].append(instance[2])
        dict_vss[tr]["maxranks"].append(instance[3])
        
def init_mean_dict(instance):
    for tr in transforms:
        dict_mean_det[tr]["schools"].append(instance[0])
        dict_mean_det[tr]["students"].append(instance[1])
        dict_mean_det[tr]["maxseats"].append(instance[2])
        dict_mean_det[tr]["maxranks"].append(instance[3])
        
        dict_mean_stoch[tr]["schools"].append(instance[0])
        dict_mean_stoch[tr]["students"].append(instance[1])
        dict_mean_stoch[tr]["maxseats"].append(instance[2])
        dict_mean_stoch[tr]["maxranks"].append(instance[3])
        
        dict_mean_vss[tr]["schools"].append(instance[0])
        dict_mean_vss[tr]["students"].append(instance[1])
        dict_mean_vss[tr]["maxseats"].append(instance[2])
        dict_mean_vss[tr]["maxranks"].append(instance[3])

def init_dict_perf(tp,tr,dict,info):  
    for sol in solvers[tp][tr]:
        if sol == 4:
            for c in conf_lag:
                for i in info:
                    key = f"{i} {name_solver[sol]} {name_lag[c]}"
                    dict[key] = []
        else:
            for i in info:
                key = f"{i} {name_solver[sol]}"
                dict[key] = []

############################### Dictionary for deterministic problem ###############################
dict_det = {}
dict_mean_det = {}
for tr in transforms:
    dict_det[tr] = {"seed":[],"schools":[],"students":[],"maxseats":[],"maxranks":[]}
    dict_mean_det[tr] = {"schools":[],"students":[],"maxseats":[],"maxranks":[]}
    
    init_dict_perf(0,tr,dict_det[tr],["status","lb","ub","gap","nb_bnb","time","gap_ub"])
    init_dict_perf(0,tr,dict_mean_det[tr],["status","lb","ub","gap","nb_bnb","time","gap_ub"])

def set_det(tr,sol,status,lb,ub,gap,nb_bnb,time,name,best_ub,list_ub):
    dict_det[tr]["status "+name].append(status)
    dict_det[tr]["lb "+name].append(lb)
    dict_det[tr]["ub "+name].append(ub)
    dict_det[tr]["gap "+name].append(100*gap)
    dict_det[tr]["nb_bnb "+name].append(nb_bnb)
    dict_det[tr]["time "+name].append(time)

    if ub != "-":
        if ub <= best_ub[0,tr] - 1e-8:
            best_ub[0,tr] = ub
    list_ub[name] = ub

def add_det(tr,sol,df,best_ub,list_ub):
    set_det(tr,sol,int(df["STATUS"][0]),float(df["LB"][0]),float(df["UB"][0]),\
            float(df["GAP"][0]),int(df["NB_BNB"][0]),float(df["TIME"][0]),name_solver[sol],best_ub,list_ub)

def add_det_lag(tr,sol,lcf,df,best_ub,list_ub):
    set_det(tr,sol,int(df["STATUS"][0]),float(df["LB"][0]),float(df["UB"][0]),\
            float(df["GAP"][0]),int(df["NB_BNB"][0]),float(df["TIME"][0]),name_solver[sol]+" "+name_lag[lcf],best_ub,list_ub)

def set_det_gapub(tr,best_ub,list_ub):
    for (key,value) in list_ub.items():
        if value != "-" and best_ub != "-":
            dict_det[tr]["gap_ub "+key].append(100*(value - best_ub)/best_ub)
        else:
            dict_stoch[tr]["gap_ub "+key].append("-")

def add_empty_det(tr,sol,best_ub,list_ub):
    set_det(tr,sol,"-","-","-","-","-","-",name_solver[sol],best_ub,list_ub) 

def add_empty_det_lag(tr,sol,lcf,best_ub,list_ub):
    set_det(tr,sol,"-","-","-","-","-","-",name_solver[sol]+" "+name_lag[lcf],best_ub,list_ub) 
    
def set_mean_dict_det(tr):
    options = ["status","lb","ub","gap","nb_bnb","time","gap_ub"]
    
    for sol in solvers[0][tr]:
        if sol == 4:
            for c in conf_lag:
                for opt in options:
                    key = f"{opt} {name_solver[sol]} {name_lag[c]}"
                    dict_mean_det[tr][key].append(str(round(np.mean(dict_det[tr][key][-len(seeds):]),2)) + "+/-" + str(round((gettstudent()*np.std(dict_det[tr][key][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))
        else:
            for opt in options:
                key = f"{opt} {name_solver[sol]}"
                dict_mean_det[tr][key].append(str(round(np.mean(dict_det[tr][key][-len(seeds):]),2)) + "+/-" + str(round((gettstudent()*np.std(dict_det[tr][key][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))   
   
############################### Dictionary for stochastic problem ###############################
dict_stoch = {}
dict_mean_stoch = {}
for tr in transforms:
    dict_stoch[tr] = {"seed":[],"schools":[],"students":[],"maxseats":[],"maxranks":[]}
    dict_mean_stoch[tr] = {"schools":[],"students":[],"maxseats":[],"maxranks":[]}
    info = ["lb1","ub1","gap1","status","nb_bnb","time_step1","gap_ub"]
    init_dict_perf(1,tr,dict_stoch[tr],info)
    init_dict_perf(1,tr,dict_mean_stoch[tr],info)

def set_stoch(tr,sol,lb,ub,gap,status,nb_bnb,time_lb,name,best_ub,list_ub):
    dict_stoch[tr]["lb1 "+name].append(lb)
    dict_stoch[tr]["ub1 "+name].append(ub)
    dict_stoch[tr]["gap1 "+name].append(100*gap)
    dict_stoch[tr]["status "+name].append(status)
    dict_stoch[tr]["nb_bnb "+name].append(nb_bnb)
    dict_stoch[tr]["time_step1 "+name].append(time_lb)
    
    if ub != "-":
        if ub <= best_ub[1,tr] - 1e-8:
            best_ub[1,tr] = ub
    list_ub[name] = ub
 
def add_stoch(tr,sol,df,best_ub,list_ub):
    if float(df["UB_1"][0]) <= 0.00001:
        print('error')
    set_stoch(tr,sol,float(df["LB_1"][0]),float(df["UB_1"][0]),float(df["GAP_1"][0]),int(df["STATUS"][0]),\
                int(df["NB_BNB"][0]),float(df["TIME_LB"][0]),name_solver[sol],best_ub,list_ub)
    
def add_stoch_lag(tr,sol,lcf,df,best_ub,list_ub):
    set_stoch(tr,sol,float(df["LB_1"][0]),float(df["UB_1"][0]),float(df["GAP_1"][0]),int(df["STATUS"][0]),\
                int(df["NB_BNB"][0]),float(df["TIME_LB"][0]),name_solver[sol]+" "+name_lag[lcf],best_ub,list_ub)
    
def set_stoch_gapub(tr,best_ub,list_ub):
    for (key,value) in list_ub.items():
        if value != "-" and best_ub != "-":
            dict_stoch[tr]["gap_ub "+key].append(100*(value - best_ub)/max(1e-3,best_ub))
        else:
            dict_stoch[tr]["gap_ub "+key].append("-")
        
def add_empty_stoch(tr,sol,best_ub,list_ub):
    set_stoch(tr,sol,"-","-","-","-","-","-",name_solver[sol],best_ub,list_ub)

def add_empty_stoch_lag(tr,sol,lcf,best_ub,list_ub):
    set_stoch_lag(tr,sol,"-","-","-","-","-","-",name_solver[sol]+" "+name_lag[lcf],best_ub,list_ub)

def set_mean_dict_stoch(tr):
    options = ["lb1","ub1","gap1","status", "nb_bnb", "time_step1", "gap_ub"]
    
    for sol in solvers[1][tr]:
        if sol == 4:
            for c in conf_lag:
                for opt in options:
                    key = f"{opt} {name_solver[sol]} {name_lag[c]}"
                    dict_mean_stoch[tr][key].append(str(round(np.mean(dict_stoch[tr][key][-len(seeds):]),2)) + "+/-" + str(round((gettstudent()*np.std(dict_stoch[tr][key][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))
        else:
            for opt in options:
                key = f"{opt} {name_solver[sol]}"
                dict_mean_stoch[tr][key].append(str(round(np.mean(dict_stoch[tr][key][-len(seeds):]),2)) + "+/-" + str(round((gettstudent()*np.std(dict_stoch[tr][key][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))   

################################### Dictionary for VSS ####################################
dict_vss = dict()
dict_mean_vss = dict()

for tr in transforms:
    dict_vss[tr] = {"seed":[],"schools":[],"students":[],"maxseats":[],"maxranks":[],"ev":[],\
                    "time_ev":[],"eev":[],"std_eev":[],"time_eev":[],"rp":[],"std_rp":[],"time_total_rp":[],\
                    "time_lb_rp":[],"time_eval_rp":[],"vss(%)":[],"nb_enter(%)":[],"nb_improv(%)":[],"nb_total(%)":[]}
    dict_mean_vss[tr] = {"schools":[],"students":[],"maxseats":[],"maxranks":[],"vss(%)":[],"nb_enter(%)":[],"nb_improv(%)":[],"nb_total(%)":[]}
    
    if tr == 3:
        for info in ['eval','nb_enter',"nb_improv","nb_total"]:
            dict_vss[tr][f'rp x stoch nonstrat {info}(%)'] = []
            dict_vss[tr][f'rp x det nonstrat {info}(%)'] = []
            # dict_vss[tr][f'eev x stoch nonstrat {info}(%)'] = []
            # dict_vss[tr][f'eev x det nonstrat {info}(%)'] = []
            
            dict_mean_vss[tr][f'rp x stoch nonstrat {info}(%)'] = []
            dict_mean_vss[tr][f'rp x det nonstrat {info}(%)'] = []
            # dict_mean_vss[tr][f'eev x stoch nonstrat {info}(%)'] = []
            # dict_mean_vss[tr][f'eev x det nonstrat {info}(%)'] = []
        
    if tr == 1 or tr == 2:
        for info in ['eval','nb_enter',"nb_improv","nb_total"]:
            dict_vss[tr][f"rp x stoch nonstrat {info}(%)"] = []
            dict_vss[tr][f"rp x det nonstrat {info}(%)"] = []
            # dict_vss[tr][f"eev x stoch nonstrat {info}(%)"] = []
            # dict_vss[tr][f"eev x det nonstrat {info}(%)"] = []
            
            dict_mean_vss[tr][f"rp x stoch nonstrat {info}(%)"] = []
            dict_mean_vss[tr][f"rp x det nonstrat {info}(%)"] = []
            # dict_mean_vss[tr][f"eev x stoch nonstrat {info}(%)"] = []
            # dict_mean_vss[tr][f"eev x det nonstrat {info}(%)"] = []
            
            dict_vss[tr][f"rp x stoch simple {info}(%)"] = []
            dict_vss[tr][f"rp x det simple {info}(%)"] = []
            # dict_vss[tr][f"eev x stoch simple {info}(%)"] = []
            # dict_vss[tr][f"eev x det simple {info}(%)"] = []
            
            dict_mean_vss[tr][f"rp x stoch simple {info}(%)"] = []
            dict_mean_vss[tr][f"rp x det simple {info}(%)"] = []
            # dict_mean_vss[tr][f"eev x stoch simple {info}(%)"] = []
            # dict_mean_vss[tr][f"eev x det simple {info}(%)"] = []
        
    for i in range(max(nbranks)+1):
        dict_vss[tr]["diff_rank"+str(i+1)] = []
        dict_mean_vss[tr]["diff_rank"+str(i+1)] = []

def set_seed_dict_vss(tr,vss_det,vss_stoch,avg_rank_det,avg_rank_stoch,max_rank):
    dict_vss[tr]["ev"].append(vss_det[tr]["ev"])
    dict_vss[tr]["time_ev"].append(vss_det[tr]["time_ev"])
    dict_vss[tr]["eev"].append(vss_det[tr]["eev"])
    dict_vss[tr]["std_eev"].append(vss_det[tr]["std_eev"])
    dict_vss[tr]["time_eev"].append(vss_det[tr]["time_eev"])
    
    dict_vss[tr]["rp"].append(vss_stoch[tr]["rp"])
    dict_vss[tr]["std_rp"].append(vss_stoch[tr]["std_rp"])
    dict_vss[tr]["time_total_rp"].append(vss_stoch[tr]["time_total_rp"])
    dict_vss[tr]["time_lb_rp"].append(vss_stoch[tr]["time_lb_rp"])
    dict_vss[tr]["time_eval_rp"].append(vss_stoch[tr]["time_eval_rp"])
    
    dict_vss[tr]["vss(%)"].append(100*(vss_det[tr]["eev"] - vss_stoch[tr]["rp"])/vss_stoch[tr]["rp"])
    dict_vss[tr]["nb_enter(%)"].append(100*(-vss_det[tr]["nb_enter"] + vss_stoch[tr]["nb_enter"])/vss_stoch[tr]["nb_enter"])
    dict_vss[tr]["nb_improv(%)"].append(100*(-vss_det[tr]["nb_improv"] + vss_stoch[tr]["nb_improv"])/vss_stoch[tr]["nb_improv"])
    dict_vss[tr]["nb_total(%)"].append(100*(-vss_det[tr]["nb_enter"] -vss_det[tr]["nb_improv"] + vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"])/(vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"]))
    
    if tr == 3: # simple strategy
        dict_vss[tr]["rp x stoch nonstrat eval(%)"].append(100*(vss_stoch[0]["eval_simple"] - vss_stoch[tr]["rp"])/vss_stoch[tr]["rp"])
        dict_vss[tr]["rp x det nonstrat eval(%)"].append(100*(vss_det[0]["eval_simple"] - vss_stoch[tr]["rp"])/vss_stoch[tr]["rp"])
            
        # dict_vss[tr]["eev x stoch nonstrat eval(%)"].append(100*(vss_stoch[0]["eval_simple"] - vss_det[tr]["eev"])/vss_det[tr]["eev"])
        # dict_vss[tr]["eev x det nonstrat eval(%)"].append(100*(vss_det[0]["eval_simple"] - vss_det[tr]["eev"])/vss_det[tr]["eev"])
        
        for info in ['nb_enter',"nb_improv"]:
            dict_vss[tr][f"rp x stoch nonstrat {info}(%)"].append(100*(-vss_stoch[0][f"{info}_simple"] + vss_stoch[tr][info])/vss_stoch[tr][info])
            dict_vss[tr][f"rp x det nonstrat {info}(%)"].append(100*(-vss_det[0][f"{info}_simple"] + vss_stoch[tr][info])/vss_stoch[tr][info])
        
        dict_vss[tr]["rp x stoch nonstrat nb_total(%)"].append(100*(-vss_stoch[0]["nb_enter_simple"] -vss_stoch[0]["nb_improv_simple"] + vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"])/(vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"]))
        dict_vss[tr]["rp x det nonstrat nb_total(%)"].append(100*(-vss_det[0]["nb_enter_simple"] -vss_det[0]["nb_improv_simple"] + vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"])/(vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"]))
            
    if tr == 1 or tr == 2: # complex strategy
        dict_vss[tr]["rp x stoch nonstrat eval(%)"].append(100*(vss_stoch[0]["eval_complex"] - vss_stoch[tr]["rp"])/vss_stoch[tr]["rp"])
        dict_vss[tr]["rp x det nonstrat eval(%)"].append(100*(vss_det[0]["eval_complex"] - vss_stoch[tr]["rp"])/vss_stoch[tr]["rp"])
        # dict_vss[tr]["eev x stoch nonstrat eval(%)"].append(100*(vss_stoch[0]["eval_complex"] - vss_det[tr]["eev"])/vss_det[tr]["eev"])
        # dict_vss[tr]["eev x det nonstrat eval(%)"].append(100*(vss_det[0]["eval_complex"] - vss_det[tr]["eev"])/vss_det[tr]["eev"])
        
        for info in ['nb_enter',"nb_improv"]:
            dict_vss[tr][f"rp x stoch nonstrat {info}(%)"].append(100*(-vss_stoch[0][f"{info}_complex"] + vss_stoch[tr][info])/vss_stoch[tr][info])
            dict_vss[tr][f"rp x det nonstrat {info}(%)"].append(100*(-vss_det[0][f"{info}_complex"] + vss_stoch[tr][info])/vss_stoch[tr][info])
        
        dict_vss[tr]["rp x stoch nonstrat nb_total(%)"].append(100*(-vss_stoch[0]["nb_enter_complex"] -vss_stoch[0]["nb_improv_complex"] + vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"])/(vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"]))
        dict_vss[tr]["rp x det nonstrat nb_total(%)"].append(100*(-vss_det[0]["nb_enter_complex"] -vss_det[0]["nb_improv_complex"] + vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"])/(vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"]))
        
        dict_vss[tr]["rp x stoch simple eval(%)"].append(100*(vss_stoch[3]["eval_complex"] - vss_stoch[tr]["rp"])/vss_stoch[tr]["rp"])
        dict_vss[tr]["rp x det simple eval(%)"].append(100*(vss_det[3]["eval_complex"] - vss_stoch[tr]["rp"])/vss_stoch[tr]["rp"])
        # dict_vss[tr]["eev x stoch simple eval(%)"].append(100*(vss_stoch[3]["eval_complex"] - vss_det[tr]["eev"])/vss_det[tr]["eev"])
        # dict_vss[tr]["eev x det simple eval(%)"].append(100*(vss_det[3]["eval_complex"] - vss_det[tr]["eev"])/vss_det[tr]["eev"])
        
        for info in ['nb_enter',"nb_improv"]:
            dict_vss[tr][f"rp x stoch simple {info}(%)"].append(100*(-vss_stoch[3][f"{info}_complex"] + vss_stoch[tr][info])/vss_stoch[tr][info])
            dict_vss[tr][f"rp x det simple {info}(%)"].append(100*(-vss_det[3][f"{info}_complex"] + vss_stoch[tr][info])/vss_stoch[tr][info])
        
        dict_vss[tr]["rp x stoch simple nb_total(%)"].append(100*(-vss_stoch[3]["nb_enter_complex"] -vss_stoch[3]["nb_improv_complex"] + vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"])/(vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"]))
        dict_vss[tr]["rp x det simple nb_total(%)"].append(100*(-vss_det[3]["nb_enter_complex"] -vss_det[3]["nb_improv_complex"] + vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"])/(vss_stoch[tr]["nb_enter"] + vss_stoch[tr]["nb_improv"]))
        
    for i in range(max_rank+1):
        val = avg_rank_stoch[tr][0][i] - avg_rank_det[tr][0][i]
        dict_vss[tr]["diff_rank"+str(i+1)].append(val)
    for i in range(max_rank+1,max(nbranks)+1):
        dict_vss[tr]["diff_rank"+str(i+1)].append("-")

def set_mean_dict_vss(tr,max_rank):
    options = ["vss(%)","nb_enter(%)","nb_improv(%)","nb_total(%)"]
    if tr == 3:
        for info in ['eval','nb_enter','nb_improv','nb_total']:
            options.extend([f"rp x stoch nonstrat {info}(%)",f"rp x det nonstrat {info}(%)"])
            # options.extend([f"eev x stoch nonstrat {info}(%)",f"eev x det nonstrat {info}(%)"])
        
    if tr == 1 or tr == 2:
        for info in ['eval','nb_enter','nb_improv','nb_total']:   
            options.extend([f"rp x stoch nonstrat {info}(%)",f"rp x det nonstrat {info}(%)"])
            options.extend([f"rp x stoch simple {info}(%)",f"rp x det simple {info}(%)"])
            # options.extend([f"eev x stoch nonstrat {info}(%)",f"eev x det nonstrat {info}(%)"])
            # options.extend([f"eev x stoch simple {info}(%)",f"eev x det simple {info}(%)"])
    
    for i in range(max_rank+1):
        options.append("diff_rank"+str(i+1))

    for opt in options:
        dict_mean_vss[tr][opt].append(str(round(np.mean(dict_vss[tr][opt][-len(seeds):]),2)) + "+/-" + str(round((gettstudent()*np.std(dict_vss[tr][opt][-len(seeds):],ddof=1))/math.sqrt(len(seeds)),2)))

    for i in range(max_rank+1,max(nbranks)+1):
        dict_mean_vss[tr]["diff_rank"+str(i+1)].append("-")

def set_info_vss(df,tp,tr,max_rank,vss_det,vss_stoch,avg_rank_det,avg_rank_stoch):
    if tp == 0: # deterministc
        vss_det[tr]["ev"] = float(df["UB"][0])
        vss_det[tr]["time_ev"] = float(df["TIME"][0])
        
        if tr == 0: # non strategic
            vss_det[tr]["eev"] = float(df["VAL_STOCH_NON_STRAT"][0])
            vss_det[tr]["std_eev"] = math.sqrt(float(df["VAR_STOCH_NON_STRAT"][0]))
            vss_det[tr]["time_eev"] = float(df["TIME_STOCH_NON_STRAT"][0])
            vss_det[tr]["nb_enter"] = float(df["NB_ENTER_NON_STRAT"][0])
            vss_det[tr]["nb_improv"] = float(df["NB_IMPROV_NON_STRAT"][0])
            
            vss_det[tr]["eval_simple"] = float(df["VAL_STOCH_SIMPLE"][0])
            vss_det[tr]["time_eval_simple"] = float(df["TIME_STOCH_SIMPLE"][0])
            vss_det[tr]["nb_enter_simple"] = float(df["NB_ENTER_STRAT_SIMPLE"][0])
            vss_det[tr]["nb_improv_simple"] = float(df["NB_IMPROV_STRAT_SIMPLE"][0])
            
            vss_det[tr]["eval_complex"] = float(df["VAL_STOCH_STRAT"][0])
            vss_det[tr]["time_eval_complex"] = float(df["TIME_STOCH_STRAT"][0])
            vss_det[tr]["nb_enter_complex"] = float(df["NB_ENTER_STRAT"][0])
            vss_det[tr]["nb_improv_complex"] = float(df["NB_IMPROV_STRAT"][0])
            
            avg_rank_det[tr].append([float(df["AVG_RANK_"+str(i+1)+"_NON_STRAT"][0]) for i in range(max_rank+1)])
            avg_rank_det[tr].append([float(df["AVG_RANK_"+str(i+1)+"_SIMPLE"][0]) for i in range(max_rank+1)])
            avg_rank_det[tr].append([float(df["AVG_RANK_"+str(i+1)+"_STRAT"][0]) for i in range(max_rank+1)])
        
        elif tr == 3: # simple strategy
            vss_det[tr]["eev"] = float(df["VAL_STOCH_SIMPLE"][0])
            vss_det[tr]["std_eev"] = math.sqrt(float(df["VAR_STOCH_SIMPLE"][0]))
            vss_det[tr]["time_eev"] = float(df["TIME_STOCH_SIMPLE"][0])
            vss_det[tr]["nb_enter"] = float(df["NB_ENTER_STRAT_SIMPLE"][0])
            vss_det[tr]["nb_improv"] = float(df["NB_IMPROV_STRAT_SIMPLE"][0])
            
            vss_det[tr]["eval_complex"] = float(df["VAL_STOCH_STRAT"][0])
            vss_det[tr]["time_eval_complex"] = float(df["TIME_STOCH_STRAT"][0])
            vss_det[tr]["nb_enter_complex"] = float(df["NB_ENTER_STRAT"][0])
            vss_det[tr]["nb_improv_complex"] = float(df["NB_IMPROV_STRAT"][0])
            
            avg_rank_det[tr].append([float(df["AVG_RANK_"+str(i+1)+"_SIMPLE"][0]) for i in range(max_rank+1)])
            avg_rank_det[tr].append([float(df["AVG_RANK_"+str(i+1)+"_STRAT"][0]) for i in range(max_rank+1)])
        
        elif tr == 1 or tr == 2: # complex strategy
            vss_det[tr]["eev"] = float(df["VAL_STOCH_STRAT"][0])
            vss_det[tr]["std_eev"] = math.sqrt(float(df["VAR_STOCH_STRAT"][0]))
            vss_det[tr]["time_eev"] = float(df["TIME_STOCH_STRAT"][0])
            vss_det[tr]["nb_enter"] = float(df["NB_ENTER_STRAT"][0])
            vss_det[tr]["nb_improv"] = float(df["NB_IMPROV_STRAT"][0])
            
            avg_rank_det[tr].append([float(df["AVG_RANK_"+str(i+1)+"_STRAT"][0]) for i in range(max_rank+1)])
            
    elif tp == 1: # stochastic
        vss_stoch[tr]["rp"] = float(df["UB"][0])
        vss_stoch[tr]["std_rp"] = math.sqrt(float(df["VAR_UB"][0]))
        vss_stoch[tr]["time_total_rp"] = float(df["TIME"][0])
        vss_stoch[tr]["time_lb_rp"] = float(df["TIME_LB"][0])
        vss_stoch[tr]["time_eval_rp"] = float(df["TIME_EVAL"][0])
        
        if tr == 0: # non strategic
            vss_stoch[tr]["nb_enter"] = float(df["NB_ENTER_NON_STRAT"][0])
            vss_stoch[tr]["nb_improv"] = float(df["NB_IMPROV_NON_STRAT"][0])
            
            vss_stoch[tr]["eval_simple"] = float(df["VAL_STOCH_SIMPLE"][0])
            vss_stoch[tr]["time_eval_simple"] = float(df["TIME_STOCH_SIMPLE"][0])
            vss_stoch[tr]["nb_enter_simple"] = float(df["NB_ENTER_STRAT_SIMPLE"][0])
            vss_stoch[tr]["nb_improv_simple"] = float(df["NB_IMPROV_STRAT_SIMPLE"][0])
            
            vss_stoch[tr]["eval_complex"] = float(df["VAL_STOCH_STRAT"][0])
            vss_stoch[tr]["time_eval_complex"] = float(df["TIME_STOCH_STRAT"][0])
            vss_stoch[tr]["nb_enter_complex"] = float(df["NB_ENTER_STRAT"][0])
            vss_stoch[tr]["nb_improv_complex"] = float(df["NB_IMPROV_STRAT"][0])
            
            avg_rank_stoch[tr].append([float(df["AVG_RANK_"+str(i+1)+"_NON_STRAT"][0]) for i in range(max_rank+1)])
            avg_rank_stoch[tr].append([float(df["AVG_RANK_"+str(i+1)+"_SIMPLE"][0]) for i in range(max_rank+1)])
            avg_rank_stoch[tr].append([float(df["AVG_RANK_"+str(i+1)+"_STRAT"][0]) for i in range(max_rank+1)])
            
        if tr == 3: # simple strategy
            vss_stoch[tr]["nb_enter"] = float(df["NB_ENTER_STRAT_SIMPLE"][0])
            vss_stoch[tr]["nb_improv"] = float(df["NB_IMPROV_STRAT_SIMPLE"][0])
            
            vss_stoch[tr]["eval_complex"] = float(df["VAL_STOCH_STRAT"][0])
            vss_stoch[tr]["time_eval_complex"] = float(df["TIME_STOCH_STRAT"][0])
            vss_stoch[tr]["nb_enter_complex"] = float(df["NB_ENTER_STRAT"][0])
            vss_stoch[tr]["nb_improv_complex"] = float(df["NB_IMPROV_STRAT"][0])
            
            avg_rank_stoch[tr].append([float(df["AVG_RANK_"+str(i+1)+"_SIMPLE"][0]) for i in range(max_rank+1)])
            avg_rank_stoch[tr].append([float(df["AVG_RANK_"+str(i+1)+"_STRAT"][0]) for i in range(max_rank+1)])
        
        if tr == 1 or tr == 2:
            vss_stoch[tr]["nb_enter"] = float(df["NB_ENTER_STRAT"][0])
            vss_stoch[tr]["nb_improv"] = float(df["NB_IMPROV_STRAT"][0])
            
            avg_rank_stoch[tr].append([float(df["AVG_RANK_"+str(i+1)+"_STRAT"][0]) for i in range(max_rank+1)])              

def set_empty_info_vss(tp,tr,max_rank,vss_det,vss_stoch,avg_rank_det,avg_rank_stoch):
    if tp == 0: # deterministc
        vss_det[tr]["ev"] = "-"
        vss_det[tr]["time_ev"] = "-"
        
        if tr == 0: # non strategic
            vss_det[tr]["eev"] = "-"
            vss_det[tr]["std_eev"] = "-"
            vss_det[tr]["time_eev"] = "-"
            vss_det[tr]["eval_simple"] = "-"
            vss_det[tr]["time_eval_simple"] = "-"
            vss_det[tr]["eval_complex"] = "-"
            vss_det[tr]["time_eval_complex"] = "-"
            
            avg_rank_det[tr].append(["-" for i in range(max_rank+1)])
            avg_rank_det[tr].append(["-" for i in range(max_rank+1)])
            avg_rank_det[tr].append(["-" for i in range(max_rank+1)])
        
        elif tr == 3: # simple strategy
            vss_det[tr]["eev"] = "-"
            vss_det[tr]["std_eev"] = "-"
            vss_det[tr]["time_eev"] = "-"
            vss_det[tr]["eval_complex"] = "-"
            vss_det[tr]["time_eval_complex"] = "-"
            
            avg_rank_det[tr].append(["-" for i in range(max_rank+1)])
            avg_rank_det[tr].append(["-" for i in range(max_rank+1)])
        
        elif tr == 1 or tr == 2: # complex strategy
            vss_det[tr]["eev"] = "-"
            vss_det[tr]["std_eev"] = "-"
            vss_det[tr]["time_eev"] = "-"
            
            avg_rank_det[tr].append(["-" for i in range(max_rank+1)])
            
    elif tp == 1: # stochastic
        vss_stoch[tr]["rp"] = "-" 
        vss_stoch[tr]["std rp"] = "-" 
        vss_stoch[tr]["time_total_rp"] = "-" 
        vss_stoch[tr]["time_lb_rp"] = "-" 
        vss_stoch[tr]["time_eval_rp"] = "-" 
        
        if tr == 0: # non strategic
            vss_stoch[tr]["eval_simple"] = "-" 
            vss_stoch[tr]["time_eval_simple"] = "-" 
            vss_stoch[tr]["eval_complex"] = "-" 
            vss_stoch[tr]["time_eval_complex"] = "-" 
            
            avg_rank_stoch[tr].append(["-" for i in range(max_rank+1)])
            avg_rank_stoch[tr].append(["-" for i in range(max_rank+1)])
            avg_rank_stoch[tr].append(["-" for i in range(max_rank+1)])
            
        if tr == 3: # simple strategy
            vss_stoch[tr]["eval_complex"] = "-" 
            vss_stoch[tr]["time_eval_complex"] = "-" 
            
            avg_rank_stoch[tr].append(["-" for i in range(max_rank+1)])
            avg_rank_stoch[tr].append(["-" for i in range(max_rank+1)])
        
        if tr == 1:
            avg_rank_stoch[tr].append(["-" for i in range(max_rank+1)]) 

#################################### Final tables ###################################

def get_interval(df):
    mean = df.mean()
    margin = stats.sem(df)*stats.t.ppf(0.975, df=df.shape[0]-1)
    return str(round(mean,3)) + '+/-' + str(round(margin,3))

def get_final_table_det_stoch():
    ### Final table for comparison between deterministic and stochastic versions ###
    dict_det_stoch = dict()
    
    # Creating dictionary
    dict_det_stoch['schools'] = []
    dict_det_stoch['students'] = []
    dict_det_stoch['value'] = []
    for seat in nbseats:
        for rank in nbranks:
            dict_det_stoch['seat'+str(seat)+'_rank'+str(rank)] = []
            
    # Convert to DataFrame
    df = {tr:pd.DataFrame(dict_vss[tr]) for tr in transforms}
    
    ### Including information ###
    # Each school and student and transform
    for c in nbschools:
        for st in nbstudents:
            for tr in transforms:
                for info in ['vss','nb_enter','nb_improv']:
                    dict_det_stoch['schools'].append(c)
                    dict_det_stoch['students'].append(st)
                    dict_det_stoch['value'].append(f'tr{tr}_{info}')
                
                # For each column
                for seat in nbseats:
                    for rank in nbranks:
                        filt_df = df[tr][(df[tr]['schools'] == c) & (df[tr]['students'] == st) & (df[tr]['maxseats'] == seat) & (df[tr]['maxranks'] == rank)]
                        for info in ['vss','nb_enter','nb_improv']:
                            dict_det_stoch[f'seat{seat}_rank{rank}'].append(get_interval(filt_df[f'{info}(%)']))
                    
    # Mean
    for tr in transforms:
        for info in ['vss','nb_enter','nb_improv']:
            dict_det_stoch['schools'].append('mean')
            dict_det_stoch['students'].append('mean')
            dict_det_stoch['value'].append(f'tr{tr}_{info}')
        
        # For each column
        for seat in nbseats:
            for rank in nbranks:
                filt_df = df[tr][(df[tr]['maxseats'] == seat) & (df[tr]['maxranks'] == rank)]
                for info in ['vss','nb_enter','nb_improv']:
                    dict_det_stoch[f'seat{seat}_rank{rank}'].append(get_interval(filt_df[f'{info}(%)']))
    
    return dict_det_stoch

def get_final_table_strat(): 
    ### Final table for comparison between srategic versions ###
    dict_strat = dict()  
    
    # Initializing dictionary
    dict_strat['value'] = []
    for seat in nbseats:
        for rank in nbranks:
            dict_strat['seat'+str(seat)+'_rank'+str(rank)] = []   
    
    # Convert to DataFrame
    df = {tr:pd.DataFrame(dict_vss[tr]) for tr in transforms}
    
    # Mean
    for tp in ['det','stoch']:
        for info in ['eval','nb_enter','nb_improv','nb_total']:
            dict_strat['value'].append(f'{tp}_complex_nonstrat_{info}')
        for info in ['eval','nb_enter','nb_improv','nb_total']:
            dict_strat['value'].append(f'{tp}_complex_simple_{info}')
        for info in ['eval','nb_enter','nb_improv','nb_total']:
            dict_strat['value'].append(f'{tp}_simple_nonstrat_{info}')
            
        # For each column
        for seat in nbseats:
            for rank in nbranks:
                filt_tr1_df = df[1][(df[1]['maxseats'] == seat) & (df[1]['maxranks'] == rank)]
                filt_tr3_df = df[3][(df[3]['maxseats'] == seat) & (df[3]['maxranks'] == rank)]
                
                for info in ['eval','nb_enter','nb_improv','nb_total']:
                    dict_strat['seat'+str(seat)+'_rank'+str(rank)].append(get_interval(filt_tr1_df[f'rp x {tp} nonstrat {info}(%)']))
                for info in ['eval','nb_enter','nb_improv','nb_total']:
                    dict_strat['seat'+str(seat)+'_rank'+str(rank)].append(get_interval(filt_tr1_df[f'rp x {tp} simple {info}(%)']))
                for info in ['eval','nb_enter','nb_improv','nb_total']:
                    dict_strat['seat'+str(seat)+'_rank'+str(rank)].append(get_interval(filt_tr3_df[f'rp x {tp} nonstrat {info}(%)']))
    
    return dict_strat

def get_final_table_perf(tr):
    ### Final table for comparison between solvers for tansform tr ###
    dict_perf = dict()
    
    # Convert to datafram
    df = pd.DataFrame(dict_stoch[tr])
    
    # Creating dictionary
    dict_perf['maxseas'] = []
    dict_perf['maxranks'] = []
    for sol in solvers[1][tr]:
        dict_perf[name_solver[sol] + '_opt_gap'] = []
        dict_perf[name_solver[sol] + '_time'] = []
        dict_perf[name_solver[sol] + '_ub_gap'] = []
        
    
    # Mean per budget and maxrank
    # Lines
    for seat in nbseats:
        for rank in nbranks:
            # Columns
            dict_perf['maxseas'].append(seat)
            dict_perf['maxranks'].append(rank)
            
            for sol in solvers[1][tr]:
                filt_df = df[(df['maxseats'] == seat) & (df['maxranks'] == rank)]
                
                # Name solver at dictionary
                name = name_solver[sol]
                if sol == 4:
                    name += " " + name_lag[conf_lag[0]]
                
                # Include info
                dict_perf[name_solver[sol] + '_opt_gap'].append(get_interval(filt_df[f"gap1 {name}"]))
                dict_perf[name_solver[sol] + '_time'].append(get_interval(filt_df[f"time_step1 {name}"]))
                dict_perf[name_solver[sol] + '_ub_gap'].append(get_interval(filt_df[f"gap_ub {name}"]))
       
    # Total mean
    dict_perf['maxseas'].append('mean')
    dict_perf['maxranks'].append('mean')
    for sol in solvers[1][tr]:        
        # Name solver at dictionary
        name = name_solver[sol]
        if sol == 4:
            name += " " + name_lag[conf_lag[0]]
        # Include info
        dict_perf[name_solver[sol] + '_opt_gap'].append(get_interval(df[f"gap1 {name}"]))
        dict_perf[name_solver[sol] + '_time'].append(get_interval(df[f"time_step1 {name}"]))
        dict_perf[name_solver[sol] + '_ub_gap'].append(get_interval(df[f"gap_ub {name}"]))
    
    # Mininum
    dict_perf['maxseas'].append('min')
    dict_perf['maxranks'].append('min')
    for sol in solvers[1][tr]:        
        # Name solver at dictionary
        name = name_solver[sol]
        if sol == 4:
            name += " " + name_lag[conf_lag[0]]
        # Include info
        dict_perf[name_solver[sol] + '_opt_gap'].append(df[f"gap1 {name}"].min())
        dict_perf[name_solver[sol] + '_time'].append(df[f"time_step1 {name}"].min())
        dict_perf[name_solver[sol] + '_ub_gap'].append(df[f"gap_ub {name}"].min())
    
    # Maximum
    dict_perf['maxseas'].append('max')
    dict_perf['maxranks'].append('max')
    for sol in solvers[1][tr]:        
        # Name solver at dictionary
        name = name_solver[sol]
        if sol == 4:
            name += " " + name_lag[conf_lag[0]]
        # Include info
        dict_perf[name_solver[sol] + '_opt_gap'].append(df[f"gap1 {name}"].max())
        dict_perf[name_solver[sol] + '_time'].append(df[f"time_step1 {name}"].max())
        dict_perf[name_solver[sol] + '_ub_gap'].append(df[f"gap_ub {name}"].max())
    
    return dict_perf

###############################################################################################

# Considering 95% confidence
def gettstudent():
    if(len(seeds)==1):
        return 0
    if(len(seeds)==2):
        return 12.71
    if(len(seeds)==3):
        return 4.30
    if(len(seeds)==4):
        return 3.18
    if(len(seeds)==5):
        return 2.78
    if(len(seeds)==6):
        return 2.57
    if(len(seeds)==7):
        return 2.45
    if(len(seeds)==8):
        return 2.36
    if(len(seeds)==7):
        return 2.31
    if(len(seeds)==10):
        return 2.26
    return 0

###############################################################################################

def main():
    with pd.ExcelWriter(path) as writer:
        for instance in instances:
            # Name instance file
            common_instance_file = "20bts_" + str(instance[0]) + "schools" + str(instance[1]) + "students" + str(maxnbseats) + "extracapa"
                
            problem_reading_instance = 0
            for seed in seeds:
                # Name instance file
                instance_file = common_instance_file + str(seed) + "seed" + "_maxseats" + str(instance[2]) + "_maxrank" + str(instance[3]) + "_obj" + str(instance[4]) + "_seed0_" 
            
                # Setting instances into the seed dictionaires
                init_seed_dict(seed,instance)
                
                ################################################################################################################
                # Read files to compare computational performance 
                
                # Auxiliar data for comparing upper bounds
                best_ub = {(tp,tr): float("inf") for tp in types for tr in transforms}
                list_ub = {(tp,tr): {} for tp in types for tr in transforms}
                
                # Reading info/analyzing all models. Keep best upper bound
                for (tp,tr,sol) in models:
                    if sol != 4:
                        name_file = get_name_solfile_perf(instance_file,tp,tr,sol,"")
                        # print(name_file)    
                        try:
                            df = pd.read_csv(name_file,delimiter=';')
                            
                            # Add information to specific sheets (comparison models for tables)
                            if tp == 0:
                                add_det(tr,sol,df,best_ub,list_ub[tp,tr])
                            elif tp == 1:
                                add_stoch(tr,sol,df,best_ub,list_ub[tp,tr])

                        except pd.errors.EmptyDataError:
                            if tp == 0:
                                add_empty_det(tr,sol,best_ub,list_ub[tp,tr])
                            elif tp == 1:
                                add_empty_stoch(tr,sol,best_ub,list_ub[tp,tr])
                            problem_reading_instance = 1
                    else:
                        for cf in conf_lag:
                            name_file = get_name_solfile_perf(instance_file,tp,tr,sol,cf)
                            # print(name_file)
                            try:
                                df = pd.read_csv(name_file,delimiter=';')
                                
                                # Add information to specific sheets (comparison models for tables)
                                if tp == 0:
                                    add_det_lag(tr,sol,cf,df,best_ub,list_ub[tp,tr])
                                elif tp == 1:
                                    add_stoch_lag(tr,sol,cf,df,best_ub,list_ub[tp,tr])
                                    
                            except pd.errors.EmptyDataError:
                                if tp == 0:
                                    add_empty_det_lag(tr,sol,cf,best_ub,list_ub[tp,tr])
                                elif tp == 1:
                                    add_empty_stoch_lag(tr,sol,cf,best_ub,list_ub[tp,tr]) 
                                problem_reading_instance = 1  
                
                # Including info about gap best ub 
                for tr in transforms:
                    if small_instances == False:
                        set_det_gapub(tr,best_ub[0,tr],list_ub[0,tr])
                    set_stoch_gapub(tr,best_ub[1,tr],list_ub[1,tr])           
                
                ################################################################################################################
                # Read files to to 1) Compare stoch. vs det. 2) Compare behaviors
                # These analysis are defined only for large instances
                
                if small_instances == False:
                    # Auxiliar data for taking information 
                    vss_det = {tr:{} for tr in transforms}
                    vss_stoch = {tr:{} for tr in transforms}
                    avg_rank_det = {tr:[] for tr in transforms}
                    avg_rank_stoch = {tr:[] for tr in transforms}
                    
                    # Read information to define different vss according to type (deterministic vs stochastic) 
                    # and transforms (nonstrat vs simple vs complex)
                    for tp in types:
                        for tr in transforms:
                            # Name solution file
                            name_file = get_name_solfile_vss(instance_file,tp,tr)
                            # print(name_file)   
                            try:
                                df = pd.read_csv(name_file,delimiter=';')
                                set_info_vss(df,tp,tr,instance[3],vss_det,vss_stoch,avg_rank_det,avg_rank_stoch)
                            
                            except pd.errors.EmptyDataError:
                                set_empty_info_vss(tp,tr,instance[3],vss_det,vss_stoch,avg_rank_det,avg_rank_stoch)
                                problem_reading_instance = 1

                    for tr in transforms:
                        # Adding info to VSS sheet          
                        set_seed_dict_vss(tr,vss_det,vss_stoch,avg_rank_det,avg_rank_stoch,instance[3])
                
                ###########################################################################################################

            # Setting information for mean dictionaries
            if(len(seeds) > 1):
                if problem_reading_instance == 0:
                    # Setting instance into the mean dictionaires
                    init_mean_dict(instance)
                    
                    # Setting info
                    for tr in transforms:
                        if small_instances == False:
                            set_mean_dict_det(tr)
                            set_mean_dict_vss(tr,instance[3]) 
                        set_mean_dict_stoch(tr)        

        #  Creating sheets
        for tr in transforms:
            if small_instances == False:
                det_df = pd.DataFrame(dict_det[tr])
                det_df.to_excel(writer,sheet_name="det_tr"+str(tr))
                
                det_mean_df = pd.DataFrame(dict_mean_det[tr])
                det_mean_df.to_excel(writer,sheet_name="det_mean_tr"+str(tr))
            
            stoch_df = pd.DataFrame(dict_stoch[tr])
            stoch_df.to_excel(writer,sheet_name="stoch_tr"+str(tr))
            
            stoch_mean_df = pd.DataFrame(dict_mean_stoch[tr])
            stoch_mean_df.to_excel(writer,sheet_name="stoch_mean_tr"+str(tr))
            
            if small_instances == False:
                vss_df = pd.DataFrame(dict_vss[tr])
                vss_df.to_excel(writer,sheet_name="vss_tr"+str(tr))
                
                vss_mean_df = pd.DataFrame(dict_mean_vss[tr])
                vss_mean_df.to_excel(writer,sheet_name="vss_mean_tr"+str(tr))
        
        # Creating final tables comp perf
        for tr in transforms:
            final_stoch_perf = pd.DataFrame(get_final_table_perf(tr))
            final_stoch_perf.to_excel(writer,sheet_name="perf_final_tr"+str(tr))
        
        if small_instances == False:
            # Creating final tables vss and strat
            final_vss_df = pd.DataFrame(get_final_table_det_stoch())
            final_vss_df.to_excel(writer,sheet_name="vss_final")
                
            final_strat_df = pd.DataFrame(get_final_table_strat())
            final_strat_df.to_excel(writer,sheet_name="strat_final")
 
if __name__ == "__main__":
    main()
            