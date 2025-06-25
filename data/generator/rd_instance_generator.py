import random
import argparse
import math
import numpy as np
import heapq
from matching.games import HospitalResident
import warnings
from collections import Counter
warnings.filterwarnings("ignore")

def parse_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('-nbschools', metavar='[nbschools]', \
                        help = 'number of schools',\
                        type = int, default=5)
    parser.add_argument('-nbstudents', metavar='[nbstudents]', \
                        help = 'number of students',\
                        type = int, default=20)
    parser.add_argument('-nbextracapa', metavar='[nbextracapa]', \
                        help = 'number of maximum extra capacity',\
                        type = int, default=5)
    parser.add_argument('-seed', metavar='[seed]', \
                        help = 'seed for data generation',\
                        type = int, default=0)
    parser.add_argument('-size_boots', metavar='[size bootstrap]', \
                        help = 'size for bootstrapping',\
                        type = int, default=100)
    arg = parser.parse_args()
    return arg

class Data():
    def __init__(self,arg):
        self.nb_boots_tries = arg.size_boots
        self.seed = arg.seed
        self.nb_schools = arg.nbschools
        self.nb_students = arg.nbstudents
        self.nb_extra_capa = arg.nbextracapa
        
        self.name_schools = ['c'+str(i) for i in range(self.nb_schools)]
        self.name_students = ['st'+str(j) for j in range(self.nb_students)]
        
        # utilities without correlation
        # self.utilities = {(st,c):random.uniform(-5,10) for c in range(self.nb_schools) for st in range(self.nb_students)}
        
        # with correlation by using distance
        loc_schools = [(random.uniform(0.0,10.0),random.uniform(0.0,10.0)) for c in range(self.nb_schools)]
        loc_students = [(random.uniform(0.0,10.0),random.uniform(0.0,10.0)) for st in range(self.nb_students)]
        
        max_dist = 10.0*math.pow(2,1/2) # Maximum possible distance between school and student
        max_accep_dist = 0.5*max_dist # 50% of maximum distance (positive utilities for schools with distance inferior to 50% of maximum distance)
        
        # Students have complete information, but do not consider schools that are too far
        self.utilities = {(st,c): max_accep_dist - math.pow((loc_schools[c][0] - loc_students[st][0])**2 + (loc_schools[c][1] - loc_students[st][1])**2, 1/2)\
                         for c in range(self.nb_schools) for st in range(self.nb_students)}
        
        # Computing scores
        self.scores = dict()
        for c in range(self.nb_schools):
            list_st = np.random.permutation(range(self.nb_students))
            aux = {(list_st[p],c): p+1 for p in range(len(list_st))}
            self.scores.update(aux)
        
        self.capacities = np.ones(self.nb_schools) + np.random.multinomial(self.nb_students - self.nb_schools, np.ones(self.nb_schools)/self.nb_schools, size=1)[0]
        self.computeProb()
        self.writeFile()
    
    def computeProb(self):
        # Bootstrapping

        self.nb_accepted = {(st,c,j):0 for st in range(self.nb_students) for c in range(self.nb_schools) for j in range(self.nb_extra_capa+1)}
        for i in range(self.nb_boots_tries):
            print("Boostrap ",i + 1)
            boots_st = np.random.choice(self.nb_students,self.nb_students)
            
            aux = [{bst:self.scores[boots_st[bst],c] for bst in range(self.nb_students)} for c in range(self.nb_schools)] 
            pref_schools = {c:[bst for bst,score in sorted(aux[c].items(), key=lambda item: item[1],reverse=True)] for c in range(self.nb_schools)}
            pref_schools[self.nb_schools] = [bst for bst in range(self.nb_students)]
            
            aux = [{c:self.utilities[boots_st[bst],c] for c in range(self.nb_schools) if self.utilities[boots_st[bst],c] > 0.0 } for bst in range(self.nb_students)]
            pref_students = {bst:[c for c,utility in sorted(aux[bst].items(), key=lambda item: item[1],reverse=True)] for bst in range(self.nb_students)}
            for bst in range(self.nb_students):
                if len(pref_students[bst]) == 0:
                    pref_students[bst] = [self.nb_schools]
            
            capa = {c:self.capacities[c] for c in range(self.nb_schools)}
            capa[self.nb_schools] = self.nb_students
            
            # Gale and Shapley algorithm - we are ignoring the ties (they are randomly dealt with)
            game = HospitalResident.create_from_dictionaries(pref_students, pref_schools, capa)
            solution = game.solve(optimal="resident") 
            dict_sol = {int(str(key)):[int(str(l)) for l in val] for (key,val) in solution.items()}

            for c in range(self.nb_schools):
                cutoff = self.nb_students + 1
                if len(dict_sol[c]) > 0:
                    for bst in dict_sol[c]:
                        cutoff = min(cutoff,self.scores[boots_st[bst],c])
                else:
                    cutoff = 0
                
                for st in range(self.nb_students):
                    if self.scores[st,c] >= cutoff:
                        self.nb_accepted[st,c,0] += 1
            
            for c in range(self.nb_schools):
                for j in range(self.nb_extra_capa):
                    capa[c] = self.capacities[c] + j + 1
                    game = HospitalResident.create_from_dictionaries(pref_students, pref_schools, capa)
                    solution = game.solve(optimal="resident") 
                    dict_sol = {int(str(key)):[int(str(l)) for l in val] for (key,val) in solution.items()}
                    
                    cutoff = self.nb_students + 1
                    if len(dict_sol[c]) > 0:
                        for bst in dict_sol[c]:
                            cutoff = min(cutoff,self.scores[boots_st[bst],c])
                    else:
                        cutoff = 0
                    
                    for st in range(self.nb_students):
                        if self.scores[st,c] >= cutoff:
                            self.nb_accepted[st,c,j+1] += 1
        
                    capa[c] = self.capacities[c]
                    
        self.prob = {(st,c,j): self.nb_accepted[st,c,j]/self.nb_boots_tries for st in range(self.nb_students) for c in range(self.nb_schools) for j in range(self.nb_extra_capa+1)}
        for st in range(self.nb_students):
            for c in range(self.nb_schools):
                for j in range(self.nb_extra_capa):
                    self.prob[st,c,j+1] = self.prob[st,c,j+1] - self.prob[st,c,0]
    
    def writeFile(self):
        # Writing file
        name = '../instances/' + str(self.nb_boots_tries) + 'bts_' + str(self.nb_schools) + 'schools' + str(self.nb_students) + 'students' + str(self.nb_extra_capa) + 'extracapa' + str(self.seed) + 'seed.txt'
        with open(name, 'w') as f:
            f.write( "%d %d %d\n" % (self.nb_schools,self.nb_students,self.nb_extra_capa))
            for c in range(self.nb_schools):
                f.write("%s %d\n" % (self.name_schools[c],self.capacities[c])) 
            for st in range(self.nb_students):
                for c in range(self.nb_schools):
                    f.write("%s %s %d %.3f " % (self.name_students[st],self.name_schools[c],self.scores[st,c],self.utilities[st,c]))
                    for j in range(self.nb_extra_capa+1):
                        f.write("%.3f " % (self.prob[st,c,j]))
                    f.write("\n")  

def main():
    arg = parse_arg()
    random.seed(arg.seed)
    np.random.seed(arg.seed)
    data = Data(arg)

if __name__ == "__main__":
    main()
    