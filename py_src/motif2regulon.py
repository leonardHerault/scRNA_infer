import sys, getopt
#import pandas as pd
import pickle
import json
from pyscenic.utils import load_motifs
from pyscenic.prune import df2regulons



def main(argv):
    MOTIFS_FNAME = ''
    REGULONS_FNAME = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["inputFile=","outputFile="])
    except getopt.GetoptError:
        print('motif2regulon.py -i <inputfile> -o <outputFile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('motif2regulon.py -i <inputfile> -o <outputFile>')
            sys.exit()
        elif opt in ("-i", "--inputFile"):
            MOTIFS_FNAME = arg
        elif opt in ("-o", "--outputFile"):
            REGULONS_FNAME = arg
            df = load_motifs(MOTIFS_FNAME)
            regulons = df2regulons(df)
            name2targets = {r.name: list(r.gene2weight.keys()) for r in regulons}
            with open(REGULONS_FNAME, 'w') as f:
                f.write(json.dumps(name2targets))

if __name__ == "__main__":
    main(sys.argv[1:])
