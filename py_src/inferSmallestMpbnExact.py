# coding: utf-8

# In[1]:



import sys, getopt, os
import bonesis
import pandas as pd
from colomoto_jupyter import tabulate
import mpbn
import itertools
import math
import numpy
import networkx as nx
import pickle
#import ginsim
import re
import subprocess
import numpy as np
from matplotlib import pyplot as plt
from math import isnan
import clingo as asp
import mpbn


# In[2]:


##For testing
#os.chdir('..')
#influenceGraph = "output/regulonAnalysis/infGraphTable50.tsv"
#outDir = "output/Inference"
#obsData = "input/obsData.csv"
#ncores = 8

#home made functions
sys.path.append('py_src/')

import constraints
import funForBonesis



def main(argv):
    influenceGraph = ''
    outDir = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:d:c:",["influenceGraph=","outDir=","obsData=","ncores="])
    except getopt.GetoptError:
        print('inferMpbn.py -i <influenceGraph> -o <outDir> -d <obsData> -c <ncores>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('inferMpbn.py -i <influenceGraph> -o <outDir> -d <obsData> -c <ncores>')
            sys.exit()
        elif opt in ("-i", "--influenceGraph"):
            influenceGraph = arg
        elif opt in ("-o", "--outDir"):
            outDir = os.getcwd() +"/"+ arg
            print('Output dir is "', outDir)
        elif opt in ("-d", "--obsData"):
            obsData = arg
            print('Output dir is "', outDir)
        elif opt in ("-c", "--ncores"):
            ncores = arg
            print('number of cores used is "', ncores)


    # Load influence graph:
    influenceGraphTable = pd.read_table(influenceGraph)
    influenceGraphTable = influenceGraphTable.drop_duplicates(subset=['tf', 'target', 'mor'])
    len(influenceGraphTable)
    
    # Loading constraints
    
    # Creating influence graph
    inf = []
    for r in influenceGraphTable.index:
        inf.append((influenceGraphTable["tf"][r],influenceGraphTable["target"][r],dict(sign= influenceGraphTable["mor"][r])))

    dom0 = bonesis.InfluenceGraph(inf, maxclause = 3,exact=False)
    allEdges = funForBonesis.multiDiGraphToEdgeList(dom0)
    
    
    dataTable = pd.read_csv(obsData,index_col = 0)
    data = dataTable.to_dict("index")
    clean_data = dict() #Bug of mpbn if nan in dict
    for o in data.keys():
        clean_data[o] = {k: data[o][k] for k in data[o] if not isnan(data[o][k])}

    data= clean_data
    
    bo0 = constraints.buildConstraints(inf = dom0,data = data,exact = False)


    # In[4]:


    view = bo0.boolean_networks()   
    filename = "optimization"
    view.standalone(output_filename=filename+".asp")
    a_file = open(filename+".asp", "r")
    list_of_lines = a_file.readlines()


    # In[5]:



    insert_at = len(list_of_lines)-4   # Index starting from which multiple elements will be inserted

    insert_elements = ['% Add an edges number minimization\n',
                       'nedges(N,X) :- clause(N,_,_,_),X = #count{L,S: clause(N,_,L,S)}.\n',
                       'totedges(T) :- T=#sum{V,O: nedges(O,V)}.\n',
                       '#minimize{N:totedges(N)}.\n']


    list_of_lines[insert_at:insert_at] = insert_elements


    # In[6]:


    a_file = open(filename+".asp", "w")
    a_file.writelines(list_of_lines)
    a_file.close()
    solving = funForBonesis.shell(['sh',filename+".asp"])
    #


    # In[7]:


    results = [s for s in solving[0].split('\\n') if "Optimum" in s] 
    minimumEdges = int([s for s in solving[0].split('\\n') if "Optimization" in s][-1].split(':')[1])
        


    # In[8]:


    subprocess.call(['rm',filename+".asp"])


    # In[9]:


    view = bo0.boolean_networks()   
    filename = "miniEdgeSol"
    view.standalone(output_filename=filename+".asp")
    a_file = open(filename+".asp", "r")
    list_of_lines = a_file.readlines()


    # In[10]:


    insert_at = len(list_of_lines)-4   # Index starting from which multiple elements will be inserted

    insert_elements = ['% Add an edges number limitation\n',
                       'nedges(N,X) :- clause(N,_,_,_),X = #count{L,S: clause(N,_,L,S)}.\n',
                       'totedges(T) :- T=#sum{V,O: nedges(O,V)}.\n',
                       ':- totedges(T) ; T>'+str(minimumEdges)+'.\n']


    list_of_lines[insert_at:insert_at] = insert_elements


    # In[11]:


    a_file = open(filename+".asp", "w")
    a_file.writelines(list_of_lines)
    a_file.close()
    solving = funForBonesis.shell(['sh',filename+".asp"])
    subprocess.call(['rm',filename+".asp"])


    # In[12]:


    solving = [s for s in solving[0].split('\\n')]
    #answerList = [l.strip("\n") for l in answerList]


    # In[13]:


    answerList = [l for l in solving if l.startswith("clause")]

    len(answerList)


    # In[14]:


    solutions = [funForBonesis.minibn_of_facts2(funForBonesis.getAnswerDict(sol)) for sol in answerList]


    # In[15]:


    pickle.dump(solutions, open(outDir+"/solutions.p", "wb" ),fix_imports=True,protocol=1 )


if __name__ == "__main__":
    main(sys.argv[1:])

