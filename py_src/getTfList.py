import sys, getopt
import pandas as pd

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["inputFile=","outputFile="])
    except getopt.GetoptError:
        print('getTfList.py -i <inputfile> -o <outputFile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('getTfList.py -i <inputfile> -o <outputFile>')
            sys.exit()
        elif opt in ("-i", "--inputFile"):
            inputFile = arg
        elif opt in ("-o", "--outputFile"):
            outputFile = arg
            df_motifs_mgi = pd.read_csv(inputFile, sep='\t')
            mm_tfs = df_motifs_mgi.gene_name.unique()
            with open(outputFile, 'wt') as f:
                f.write('\n'.join(mm_tfs) + '\n')
                len(mm_tfs)

if __name__ == "__main__":
    main(sys.argv[1:])
