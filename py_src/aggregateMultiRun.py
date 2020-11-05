import sys, getopt, os
import json
from collections import Counter
import pandas as pd
import sqlite3
import statistics as stat
import time



def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hr:a:c:o:",["multipleRunFolder=","multipleAdjaFolder","cutOffProp=","outDir="])
    except getopt.GetoptError:
        print('aggregateScenicRuns.py -r <multipleRunFolder> -c <cutOffProp> -o <outDir>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('aggregateScenicRuns.py -r <multipleRunFolder> -c <cutOffProp> -o <outDir>')
            sys.exit()
        elif opt in ("-r", "--multipleRunFolder"):
            multipleRunFolder = arg
        elif opt in ("-c", "--cutOffProp"):
            cutOffProp = float(arg)
        elif opt in ("-o", "--outDir"):
            outDir = os.getcwd() +"/"+ arg
            print('Output file is "', outDir)
        elif opt in ("-a", "--multipleAdjaFolder"):
            multipleAdjaFolder = arg
            print('Adjacencies directory is "', multipleAdjaFolder)
            
    regulonsOfInterest = ["Gata2","Runx1","Klf1","Cebpa","Gata1","Fli1","Tal1","Spi1","Ikzf1","Myc","Junb","Erg","Egr1","Myb","Hoxb4","Gfi1b","Zfpm1"]
    regulonsOfInterestPos = [s + "(+)" for s in regulonsOfInterest]
    regulonsOfInterestNeg = [s + "(-)" for s in regulonsOfInterest]
    
    #multipleRunFolder = "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/"
    #multipleAdjaFolder = "output/ScenicRNA_multipleRuns/GRNboost/"
    #cutOffProp = 0.8
    #outDir = "/shared/projects/scRNA_HSPC_Aging/scRNA_infer/output/ScenicRNA_multipleRuns/cis_target_maskDropouts/"

    jsonFiles = [each for each in os.listdir(multipleRunFolder) if (each.endswith('.json')&each.startswith("run"))]
    cutOff = cutOffProp*len(jsonFiles)
    regulons = []
    aggregatedJson = dict()
    for f in jsonFiles:
        with open(multipleRunFolder+"/"+f, encoding='utf-8-sig') as json_file:
            text = json_file.read()
            json_data = json.loads(text)
            regulons += json_data.keys()
            for r in json_data.keys():
                if (r in aggregatedJson):
                    aggregatedJson[r] += json_data[r]
                else:
                    aggregatedJson[r] = json_data[r]
            

    regulonCount = Counter(regulons)
    regulonKept = {k: v for k, v in regulonCount.items() if v >= cutOff}
            
    aggregatedJsonKept = { key: aggregatedJson[key] for key in regulonKept }
            
    for r in regulonKept:
        aggregatedJsonKept[r] = {k: v for k, v in Counter(aggregatedJsonKept[r]).items() if v >= cutOff} 
            
            
    aggregatedJsonKeptNonEmpty = {k: v for k, v in aggregatedJsonKept.items() if len(v) > 0}

    regulonKept = aggregatedJsonKeptNonEmpty.keys()
            
            
            
    t = [r in regulonKept for r in regulonsOfInterestPos]
            
    t2 = [r in regulonKept for r in regulonsOfInterestNeg]
            
    print("Activating regulons of interest present:")
    print([regulonsOfInterestPos[i] for i, x in enumerate(t) if x])
            
    print("Activating regulons of interest absent:")
    print([regulonsOfInterestPos[i] for i, x in enumerate(t) if not x])
            
    print("Repressing regulons of interest present:")
    print([regulonsOfInterestNeg[i] for i, x in enumerate(t2) if x])
            
    print("Repressing regulons of interest absent:")
    print([regulonsOfInterestNeg[i] for i, x in enumerate(t2) if not x])
            
            
            #Writing Json file, scenic type
    aggregatedJsonScenic = {key : list(aggregatedJsonKeptNonEmpty[key]) for key in aggregatedJsonKeptNonEmpty.keys()}
    with open(outDir+'/aggregatedRegulons.json', 'w') as json_file:
        json.dump(aggregatedJsonScenic, json_file)
            
    #writing a Json with target recovered time number
    #adding regulons scores, mean and sd and regulon recovered time number
    
    adjaFiles = [adjaFile for adjaFile in os.listdir(multipleAdjaFolder) if (adjaFile.endswith('.tsv')&adjaFile.startswith("run"))]


    #Create a new database file:
    from sqlalchemy import create_engine

    print("creating database")
    db = sqlite3.connect(outDir+"adja.sqlite")
    print("creating engine")
    db1 = create_engine('sqlite:////'+outDir+'adja.sqlite')
    print("connecting")
    conn = db1.connect()

    regulonsRunDic = dict()
    
    for n in range(len(adjaFiles)):
        print("start loop")
        start_time = time.time()   
        adjacencies = pd.read_csv(multipleAdjaFolder+"/run"+str(n)+"_GRNboost.tsv",sep = "\t")
        adjacencies["run"] = n
        adjacencies['TF_target'] = adjacencies['TF'] + "_" + adjacencies['target']
        #adjaTable = adjaTable.append(adjacencies)
        adjacencies.to_sql("adja", conn, if_exists="append",chunksize = 1000, method = "multi")
        with open(multipleRunFolder+"/"+"run"+str(n)+"_regulons.json", encoding='utf-8-sig') as json_file:
            text = json_file.read()
            json_data = json.loads(text)
        
        regulonsRunDic[n] = json_data
        print("--- %s seconds ---" % (time.time() - start_time))
        print(n)
    
    conn.execute("CREATE INDEX TF_target ON adja(TF_target)") 
    #conn.close()
   
    def write_query(runList):
        sql = "select * from adja where TF_target=? and run in  ({1});"
        sql = sql.format('?', ','.join('?' * len(runList)))
        return sql

    def flatten(l):
        for el in l:
            try:
                yield from flatten(el)
            except TypeError:
                yield el
            


    def search(myDict, search1):
        search.a=[]
        for key, value in myDict.items():
            if search1 in value:
                search.a.append(key)
        return search.a


    
    

#adjacencies.to_sql("adja", db1, if_exists="append",chunksize = 100,method = "multi")
#print("--- %s seconds ---" % (time.time() - start_time))


# import time
# start_time = time.time()
# adjacencies = pd.read_table(multipleAdjaFolder+"/run"+str(n)+"_GRNboost.tsv")
# adjacencies["run"] = n
# adjacencies['TF_target'] = adjacencies['TF'] + "_" + adjacencies['ls Scentarget']
# #adjaTable = adjaTable.append(adjacencies)
# adjacencies.to_sql("adja", db, if_exists="append",chunksize = 1000)
# print("--- %s seconds ---" % (time.time() - start_time))
        

    #conn = sqlite3.connect("adja.sqlite")
    #for r in aggregatedJsonKeptNonEmpty.keys()[1:2]:
    for r in aggregatedJsonKeptNonEmpty.keys():
        print(r)
        tfName = r.split("(")[0]
        runRecovReg = search(regulonsRunDic, r)
        targets = aggregatedJsonKeptNonEmpty[r]
        for t in targets:
            # if(t != r.split("(")[0]): # No GRNBoost importance for TF itself
            #     aggregatedJsonKeptNonEmpty[r][t] = [aggregatedJsonKeptNonEmpty[r][t],adjaTable[(adjaTable["TF"] == r.split("(")[0])&(adjaTable["target"] == t)&(adjaTable["run"].isin(runRecovReg)) ]["importance"].tolist()]
            # else:
            #     aggregatedJsonKeptNonEmpty[r][t] = [aggregatedJsonKeptNonEmpty[r][t],[None]]
            runRecovRegTarget0 = [t in regulonsRunDic[run][r] for run in runRecovReg]
            runRecovRegTarget = [i for i, x in enumerate(runRecovRegTarget0) if x]
            if(t != r.split("(")[0]): # No GRNBoost importance for TF itself
                tfTarget =  tfName+"_"+t
                aggregatedJsonKeptNonEmpty[r][t] = [aggregatedJsonKeptNonEmpty[r][t],pd.read_sql_query(write_query(runRecovRegTarget), conn,params = tuple( [tfTarget]+runRecovRegTarget))["importance"].tolist()]
                aggregatedJsonKeptNonEmpty[r][t] = [aggregatedJsonKeptNonEmpty[r][t],stat.mean(aggregatedJsonKeptNonEmpty[r][t][1]),stat.stdev(aggregatedJsonKeptNonEmpty[r][t][1])]
            else:
                aggregatedJsonKeptNonEmpty[r][t] = [aggregatedJsonKeptNonEmpty[r][t],[None],None,None]    


    
    with open(outDir+'/aggregatedRegulonsMeta.json', 'w') as json_file:
        json.dump(aggregatedJsonKeptNonEmpty, json_file)
                
           

if __name__ == "__main__":
    main(sys.argv[1:])







