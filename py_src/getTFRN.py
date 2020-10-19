import sys, getopt
import json
import pandas as pd


def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hr:t:o:",["regulonsFile=", "tfFile=","outFile="])
    except getopt.GetoptError:
        print('getTFRN.py -r <regulonsfile> -t <tfFile -o <outputFile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('getTFRN.py -r <regulonsfile> -t <tfFile -o <outputFile>')
            sys.exit()
        elif opt in ("-r", "--regulonsFile"):
            regulonsFile = arg
        elif opt in ("-t", "--tfFile"):
            tfFile = arg
        elif opt in ("-o", "--outFile"):
            outputFile = arg
            print('Output file is "', outputfile)

            with open(regulonsFile, encoding='utf-8-sig') as json_file:
                text = json_file.read()
                json_data = json.loads(text)
                print(json_data)
                                
            json_data.keys()
            len(json_data.keys())
                                        
            with open(tfFile) as TF_file:
                tf_list = TF_file.read()
                print(tf_list)
                                                    
            tf_set = tf_list.split('\n')
            len(tf_set)
                                                        
            tf_regulons = dict()
                                                        
            for i in json_data.keys():
                print(i)
                j = json_data.get(i)
                tf_regulons[i] = set(j).intersection(tf_set)
                            
            #writing csv for GRN
            grn = pd.DataFrame(columns=['TF','interaction','Target'])
                                
            k=0
            for i in range(0,len(tf_regulons.keys())):
                TF = list(tf_regulons.keys())[i]
                TFname = TF.split('(')[0]
                interaction = TF.split('(')[1].split(')')[0]
                for j in tf_regulons.get(TF):
                    dfAdd = pd.DataFrame([[TFname,interaction,j]],columns = ['TF','interaction','Target'])
                    grn = grn.append(dfAdd,ignore_index=True)
                                                        
            grn.to_csv(outputFile)


if __name__ == "__main__":
    main(sys.argv[1:])







