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


def multiDiGraphToEdgeList(graph):
    edges = graph.edges(data=True)
    edgesTable = pd.DataFrame(edges)
    #len(edgesMinTable)
    edgesTable[3] = [edgesTable[2][i]['sign'] for i in range(len(edgesTable))]
    edgesTable = edgesTable[[0,1,3]]
    edgesTable.rename(columns={0: 'tf',1:'target',3:'mor'},inplace=True)
    
    inf = []
    for r in edgesTable.index:
        inf.append((edgesTable["tf"][r],edgesTable["target"][r],dict(sign= edgesTable["mor"][r])))
    
    return inf

def has_cyclic(bn):
        mbn = mpbn.MPBooleanNetwork(bn)
        for a in mbn.attractors():
            if "*" in a.values():
                return True
        return False

#def checkEdge(iEdges,infGraph,data):
#    edges = multiDiGraphToEdgeList(infGraph)
#    edges.pop(iEdges)
#    dom1 = bonesis.InfluenceGraph(edges, maxclause = 3,exact=False)
#    bo1 = constraints.buildConstraints(inf = dom1,data = data,exact = False)
#    view = bo1.boolean_networks(limit = 1)   
#    filename = "test"
#    view.standalone(output_filename=filename+".asp")
#    a_file = open(filename+".asp", "r")
#    list_of_lines = a_file.readlines()
#    list_of_lines[1] = 'clingo -t 24 1 --time-limit=10 --project -c bounded_nonreach=3 "${@}" - <<EOF\n'
#    a_file = open(filename+".asp", "w")
#    a_file.writelines(list_of_lines)
#    a_file.close()
#    solving = shell(['sh',filename+".asp"])
#    subprocess.call(['rm',filename+".asp"])
#    res = [x for x in  solving[0].split('\\n') if x.startswith('Models')]
#    out = res[0].find("1") != -1
#    if out:
#        return(dom1)
        
def shell(command):
    try:
        output = subprocess.check_output(command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        output = str(e.output)
    finished = output.split('\n')
    #for line in finished:
        #print(line)
    return finished

def getAnswerDict(answer):
    answerDict = answer.split(' ')
    res = []
    for e in answerDict :
        clause= {}
        clause["name"] = e.split("(")[0]
        clause["arguments"] = '('+e.split("(")[1]
        res.append(clause)
    return(res)
        
def dnfs_of_facts2(fs):
    bn = {}
    for d in fs:
        if d["name"] == "clause":
            (i,cid,lit,sign) = eval(d["arguments"])
            if i not in bn:
                bn[i] = []
            if cid > len(bn[i]):
                bn[i] += [set() for j in range(cid-len(bn[i]))]
            bn[i][cid-1].add((sign,lit))
        elif d["name"] == "constant" and len(d["arguments"]) == 2:
            (i,v) = list(map(py_of_symbol, d["arguments"]))
            bn[i] = v == 1
    return(bn)

def minibn_of_facts2(fs):
    dnfs = dnfs_of_facts2(fs)
    bn = mpbn.MPBooleanNetwork()
    def make_lit(l):
        s,v=l
        v = bn.v(v)
        if s < 0:
            v = ~v
        return v
    def make_clause(ls):
        ls = list(map(make_lit, ls))
        if len(ls) == 1:
            return ls[0]
        return bn.ba.AND(*ls)
    def make_dnf(cs):
        if isinstance(cs, bool):
            return cs
        cs = filter(len, cs)
        cs = list(map(make_clause, cs))
        if len(cs) == 1:
            return cs[0]
        return bn.ba.OR(*cs)
    for (node, cs) in sorted(dnfs.items()):
        bn[node] = make_dnf(cs)
    return bn
