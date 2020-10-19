from pyscenic.genesig import Regulon
from pyscenic.aucell import aucell
from multiprocessing import cpu_count
import sys, getopt
import json
import pandas as pd

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hr:e:c:s:o:",["aggregatedJson=","expMtx","cores","seed","outDir="])
    except getopt.GetoptError:
        print('aggregateScenicRuns.py -r <aggregatedJson> -e <expMtx>-c <cores> -s <seed> -o <outDir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('aggregateScenicRuns.py -r <aggregatedJson> -e <expMtx> -c <cores> -s <seed> -o <outDir>')
            sys.exit()
        elif opt in ("-c", "--cores"):
            cores = int(arg)
        elif opt in ("-s", "--seed"):
            seed = int(arg)
        elif opt in ("-e", "--expMtx"):
            expMtx = arg
        elif opt in ("-r", "--aggregatedJson"):
            aggregatedJson = arg
        elif opt in ("-o", "--outDir"):
            outDir = arg
            print('Outdir is "', outDir)

# Put all weight and score to 1 to be able to use aucell function but only it is only valid if noweigth = True
    with open(aggregatedJson, encoding='utf-8-sig') as json_file:
        text = json_file.read()
        regulonDic = json.loads(text)

    exp_mtx = pd.read_table(expMtx)

    regulonList = list()
    for r in regulonDic.keys():
        if ("("+r.split("(")[1] == "(+)"):
            context = frozenset({"activating",None})
        else:
            context = frozenset({"repressing",None})
        
        tf_name = r.split("(")[0]
    #print(list(zip([gene for gene in regulonDic[r]],[None for gene in regulonDic[r]])))
        regulon = Regulon(name= r,
                          score=1,
                          context=context,
                          transcription_factor=tf_name,
                          gene2weight=list(zip([gene for gene in regulonDic[r]],[1 for gene in regulonDic[r]])),
                          gene2occurrence=[])
        regulonList.append(regulon)
    
    auc_mtx = aucell(exp_mtx, regulonList, num_workers=cores,seed = 2020,noweights = True)
    auc_mtx.to_csv(outDir+"/regulons_enrichment.csv")
    
if __name__ == "__main__":
    main(sys.argv[1:])
