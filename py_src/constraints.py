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

#def mutant_constraint(mutation, phenotypes,bo):
#    with bo.mutant(mutation) as m:
#        # each phenotype can be reached from at least one init
#        for ph in phenotypes:
#            +m.obs("iHSC") >= m.fixed(+m.obs(ph))
#        # each possible initial configuration can only reach fixed points matching with phenotypes
#        for cfg in bonesis.matching_configurations(m.obs("iHSC")):
#            cfg >> "fixpoints" ^ {m.obs(ph) for ph in phenotypes}

#def buildConstraints(inf,data, maxclause = 3,exact =True,parallel=24):
#    dom1 = bonesis.InfluenceGraph(inf, maxclause = maxclause,exact=exact)
#    bo = bonesis.BoNesis(dom1, data)
#    bo.settings["parallel"] = 24
#    fLymph = bo.fixed(~bo.obs("pLymph"))
#    fEr = bo.fixed(~bo.obs("pEr"));
#    fMk = bo.fixed(~bo.obs("pMk"));
#    fNeuMast = bo.fixed(~bo.obs("pNeuMast"));
#    start = ~bo.obs("iHSC")
#    start >= ~bo.obs("srHSC")
#    start >= fLymph;
#    start >= fEr;
#    start >= fMk;
#    start >= fNeuMast;   
#    start >= ~bo.obs("qHSC");
#    start >= ~bo.obs("diff");
#    ~bo.obs("diff") >= fEr
#    ~bo.obs("diff") >= fMk
#    ~bo.obs("diff") >= fNeuMast
#    #~bo.obs("diff") >= fLymph
#    ~bo.obs("qHSC") >= ~bo.obs("diff")
#    ~bo.obs("srHSC") >= ~bo.obs("qHSC")
#    ~bo.obs("srHSC") >= start
#    ~bo.obs("qHSC") >= start
#    ~bo.obs("diff") / ~bo.obs("qHSC")
#    ~bo.obs("diff") / ~bo.obs("srHSC")
#    ~bo.obs("diff") / start
#    ~bo.obs('zero') / fNeuMast
#    ~bo.obs('zero') / fMk
#    ~bo.obs('zero') / fLymph
#    #~bo.obs('zero') / fpEr
#    #Universal reachable fixed points
#    #~bo.obs("iHSC") >> "fixpoints" ^ {bo.obs(obs) for obs in ["pLymph", "pNeuMast","pEr","pMk"]};
#    ~bo.obs("iHSC") >> "fixpoints" ^ {bo.obs(obs) for obs in ["pLymph", "pNeuMast","pEr","pMk"]}
#    mutant_constraint({"Spi1":0}, ["pEr","pMk"],bo)
#    mutant_constraint({"Cebpa":0}, ["pEr","pMk"],bo)
#    mutant_constraint({"Junb":1}, ["G0MkHSC"],bo)
#    mutant_constraint({"Junb":0}, ["prolifNeuMast","pEr","pMk","pNeuMast","pLymph"],bo)
#    return bo
