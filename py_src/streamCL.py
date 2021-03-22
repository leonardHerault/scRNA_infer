import sys, getopt
import stream as st
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os


def main(argv):
    matrixFile = ''
    outDir = ''
    try:
        opts, args = getopt.getopt(argv,"m:p:c:n:j:o:r:",["matrixFile=","annoFiles=","colorFiles=","annoNames=","jobNumber=","outDir=","root="])
    except getopt.GetoptError:
        print('streamCL.py -m <matrixFile> -p <annoFile1,annoFile2,..> -c <colorFile1,colorFile2,..> -n <annoName1,annoName2> -j <jobNumber> -o <outDir> -r <root>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('streamCL.py -m <matrixFile> -p <annoFile1,annoFile2,..> -c <colorFile1,colorFile2,..> -n <annoName1,annoName2> -j <jobNumber> -o <outDir> -r <root>')
            sys.exit()
        elif opt in ("-m", "--matrixFile"):
            matrixFile = arg
        elif opt in ("-p", "--annoFiles"):
            annoFiles = arg
        elif opt in ("-c", "--colorFiles"):
            colorFiles = arg
        elif opt in ("-n", "--annoNames"):
            annoNames = arg
        elif opt in ("-j", "--jobNumber"):
            jobNumber = int(arg)
        elif opt in ("-o", "--outDir"):
            print(arg)
            outDir = arg
        elif opt in ("-r", "--root"):
            root = arg
            print('Output dir is "', outDir)
            print('Matrix file is "', matrixFile)
            print('job number is "', jobNumber)

            
            fileName = matrixFile
            cellLabels = annoFiles.split(",")
            cellColors = colorFiles.split(",")
            labels = annoNames.split(",")
            
            n_components = 3
            flag_savefig = True
            
            cellLabel= cellLabels[0]
            cellColor = cellColors[0]
            label = annoNames[0]
            
            adata=st.read(file_name= fileName,workdir=outDir)
            st.add_cell_labels(adata,file_name= cellLabel)
            st.add_cell_colors(adata,file_name= cellColor)
            adata.obs[label] = adata.obs['label']
            
            #Dimensionnal reduction using PCA
            st.select_top_principal_components(adata,n_pc=15,first_pc=True)
            
            #Dimensionnal reduction using MLLE on 15 first pcs

            st.dimension_reduction(adata=adata,method ='mlle',feature='top_pcs',n_jobs=jobNumber,n_components=n_components,n_neighbors=30)

            # construct elastic graph
            st.seed_elastic_principal_graph(adata,n_clusters=10)
            
            st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=False,save_fig=flag_savefig,fig_name = "graphDimRedSeedTree.pdf")

            st.elastic_principal_graph(adata)
            
            st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=False,save_fig=flag_savefig,fig_name = "graphDimRedFirstTree.pdf")
            
            ###Optimize branching
            st.optimize_branching(adata)
            
            st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=False,save_fig=flag_savefig,fig_name = "graphDimRedOptimzedBranching.pdf")
            
            
            ###Extend leaf branch to reach further cells 
            st.extend_elastic_principal_graph(adata)
            #st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=False,save_fig=flag_savefig)#st.plot_branches_with_cells(adata)
            
            ## plot flat Tree
            st.plot_flat_tree(adata,show_graph=True, show_text=True,save_fig=flag_savefig)
            
            ## plot subway map
            st.plot_stream_sc(adata,root=root,fig_legend_ncol=1,save_fig=flag_savefig) 
            
            ## plot stream
            #st.plot_stream(adata,root=root,fig_legend_ncol=1,fig_size=(8,8),factor_min_win=1.2,log_scale=True,factor_zoomin=300,save_fig=flag_savefig)
            
            ## plot stream for other colors
            
            for i in range(len(labels)):
                st.add_cell_labels(adata,file_name= cellLabels[i])
                st.add_cell_colors(adata,file_name= cellColors[i])
                adata.obs[labels[i]] = adata.obs['label']
                st.plot_stream(adata,root=root,
                    fig_legend_ncol=1,
                    fig_size=(8,8),
                    factor_min_win=1.2,
                    log_scale=True,
                    factor_zoomin=300,
                    save_fig=True,
                    fig_path=outDir+"/"+labels[i])
                    
            cols = [root+"_pseudotime","branch_id_alias"] 
            table = adata.obs[cols]
            table.to_csv(outDir+'/phenoDataStream.csv')
            
            st.write(adata)


if __name__ == "__main__":
    main(sys.argv[1:])
