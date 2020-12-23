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

def buildConstraints(inf,data, maxclause = 3,exact =True,parallel=24):
    dom1 = bonesis.InfluenceGraph(inf, maxclause = maxclause,exact=exact)
    dom1
    bo = bonesis.BoNesis(dom1, data)
    bo.settings["parallel"] = parallel
    fT = bo.fixed(~bo.obs("T"))
    fEr = bo.fixed(~bo.obs("Er"));
    fMk = bo.fixed(~bo.obs("Mk"));
    fMye = bo.fixed(~bo.obs("Mye"));
    start = ~bo.obs("HSC_start")
    start >= ~bo.obs("HSC_SR")
    start >= fT;
    start >= fEr;
    start >= fMk;
    start >= fMye;   
    start >= ~bo.obs("HSC_G0");
    start >= ~bo.obs("diff");
    start >= ~bo.obs("diff2");
    ~bo.obs("diff") >= fEr
    ~bo.obs("diff") >= fMk
    ~bo.obs("diff") >= fMye
    ~bo.obs("diff2") >= fT
    ~bo.obs("HSC_G0") >= ~bo.obs("diff")
    ~bo.obs("HSC_SR") >= ~bo.obs("HSC_G0")
    ~bo.obs("HSC_SR") >= start
    ~bo.obs("HSC_G0") >= start
    ~bo.obs("HSC_G0") >= ~bo.obs("diff2")
    ~bo.obs("diff") / ~bo.obs("HSC_G0")
    ~bo.obs("diff") / ~bo.obs("HSC_SR")
    ~bo.obs("diff") / start
    ~bo.obs('zero') / fMye
    ~bo.obs('zero') / fMk
    ~bo.obs('zero') / fT
    #~bo.obs('zero') / fEr
    #Universal reachable fixed points
    ~bo.obs("HSC_start") >> "fixpoints" ^ {bo.obs(obs) for obs in ["T", "Mye","Er","Mk"]};
        #~bo.obs("HSC_start") >> "fixpoints" ^ {bo.obs(obs) for obs in ["T", "Mye","Er","Mk","HSC_G0"]};
    return bo
