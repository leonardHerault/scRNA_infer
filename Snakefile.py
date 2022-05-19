import os
import re
import pandas as pd

baseDir = config["dir"]["base"]
inputDataDir = config["dir"]["inputData"]

seurat_monocle_env = config["seurat_monocle_env"]
pyscenic_env = config["pyscenic_env"]
bonesis_env = config["bonesis_env"]
beta_env = config["beta_env"]
stream_env = config["stream_env"]
report_env = config["reportInf_env"]



ages = ["young","old"]

scenicRuns = list(range(50))


########################### Preprocessing of cistrome database############################

infofile =  config["cistrome"]['info']
f = open(infofile,"r",encoding="utf-8")
info = f.read()
f.close()
celltypes =  config["cistrome"]['celltypes']
tissue = config["cistrome"]['tissue']
c = open(celltypes,"r",encoding="utf-8")
ct = c.read()
c.close()

ct = ct.split("\n")
  
exp = info.split("\n")
exp = [e.split("\t") for e in exp[:-1]]
  
expFiltered = [e for e in exp if e[5] in ct]
peakFiles  = []
targetFiles = []
for e in expFiltered:
  peakFiles.append(e[0] + "_sort_peaks.narrowPeak.bed")
  targetFiles.append(e[2]+"_"+e[3])
  
expToFile = dict()
for e in expFiltered:
   expToFile[e[2] + "_" + e[3]] = e[0] 

expFilteredLines = ["\t".join(exp[0])+"\n"]
for n in expFiltered:
  expFilteredLines.append("\t".join(n)+"\n")

outexpFiltered = open("input/mouse_factor_full_QC_filtered.txt","w")
outexpFiltered.writelines(expFilteredLines) 
outexpFiltered.close()
   
## Same with tissue selection
expFilteredTissue = [e for e in exp if e[6] == tissue]
peakFilesTissue  = []
targetFilesTissue = []
for e in expFilteredTissue:
  peakFilesTissue.append(e[0] + "_sort_peaks.narrowPeak.bed")
  targetFilesTissue.append(e[2]+"_"+e[3])
  
expToFileTissue = dict()
for e in expFilteredTissue:
   expToFileTissue[e[2] + "_" + e[3]] = e[0] 


print(peakFilesTissue[1:5])
print(targetFilesTissue[1:5])
#expToFile.keys[1:5]

expFilteredLinesTissue = ["\t".join(exp[0])+"\n"]
for n in expFilteredTissue:
  expFilteredLinesTissue.append("\t".join(n)+"\n")

outexpFilteredTissue = open("input/mouse_factor_full_QC_filtered_tissue.txt","w")
outexpFilteredTissue.writelines(expFilteredLinesTissue) 
outexpFilteredTissue.close()


###########################################################################

rule all:
    input:#"reportProjet2/report_final.html",
          #"output/stream/stream_result.pkl",
          #"output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json",
          #"output/ScenicRNA_multipleRuns/AUCell_maskDropouts/regulons_enrichment.csv",
          #"output/regulonAnalysis/mainRegulonTable.tsv",
          #"output/publicData/trrust.tsv",
          #"output/Inference/Net/solutions.p",
          #"output/Cistrome_BM/cistromeReg.tsv",
          #"report/solutionFinal.zginml",
          "report/figures_part2.html",
          #"output/Inference/bonesis/solution_space_first_exploration.html"
          

          #"output/Inference/bonesis/solutions.p"
    
    
rule prepare_data_for_stream_all:
  input: seurat = inputDataDir+"report/seurat_report.rds",
         monocle = inputDataDir+"report/monocle_report.rds",
  output: data = "output/stream/all_scaleDataForStream.tsv.gz",
          monocleStates = "output/stream/all_monocleStates.tsv",
          seuratClusters = "output/stream/all_seuratClusters.tsv",
          ages = "output/stream/all_ages.tsv",
          colorAges = "output/stream/all_colorAges.tsv",
          colorCluster = "output/stream/all_colorClusters.tsv",
          colorStates = "output/stream/all_colorMonocleStates.tsv"
  threads: 12
  conda: seurat_monocle_env
  shell: "Rscript R_src/prepareStreamCL.R -i {input.seurat} -m {input.monocle} \
  -o output/stream/ -p all -r pL2;cd output/stream/;gzip -f all_scaleDataForStream.tsv"

rule stream_all:
  input: data = "output/stream/all_scaleDataForStream.tsv.gz",
         monocleStates = "output/stream/all_monocleStates.tsv",
         seuratClusters = "output/stream/all_seuratClusters.tsv",
         ages = "output/stream/all_ages.tsv",
         colorAges = "output/stream/all_colorAges.tsv",
         colorCluster = "output/stream/all_colorClusters.tsv",
         colorStates = "output/stream/all_colorMonocleStates.tsv"
  output: "output/stream/stream_result.pkl"
  params: annoFiles = "output/stream/all_monocleStates.tsv,output/stream/all_seuratClusters.tsv,output/stream/all_ages.tsv",
          colorNames = "output/stream/all_colorMonocleStates.tsv,output/stream/all_colorClusters.tsv,output/stream/all_colorAges.tsv",
          annoNames = "Monocle_state,Seurat_cluster,Age"
  threads: 12
  conda: stream_env
  shell: "python py_src/streamCL.py -m {input.data} \
  -p {params.annoFiles} -c {params.colorNames} -n {params.annoNames} \
  -j 12 -o output/stream/ -r S1"

rule prepare_tf_list:
     input: inputDataDir+"/publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
     output: "output/publicData/mm_mgi_tfs.txt"
     conda: pyscenic_env
     shell: "python py_src/getTfList.py -i {input} -o {output}"
     
rule split_by_age:
    input: inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
    output: "output/ScenicRNA_multipleRuns_{age}/expressionRawCountFilteredForScenic.tsv"
    conda: pyscenic_env
    shell:
        "head -1 {input} > {output}; cat {input} | grep '^{wildcards.age}_[A-Z]_' >> {output}"

rule multiple_grnboost_rna:
    input: inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
            "output/publicData/mm_mgi_tfs.txt"
    output: "output/ScenicRNA_multipleRuns/GRNboost/run{run}_GRNboost.tsv"
    threads: 20
    conda: pyscenic_env
    shell:
        "pyscenic grn --seed {wildcards.run} --num_workers {threads} {input[0]} {input[1]} -o {output}"

rule multiple_grnboost_rna_by_age:
    input: "output/ScenicRNA_multipleRuns_{age}/expressionRawCountFilteredForScenic.tsv",
           "output/publicData/mm_mgi_tfs.txt"
    output: "output/ScenicRNA_multipleRuns_{age}/GRNboost/run{run}_GRNboost.tsv"
    threads: 20
    conda: pyscenic_env
    shell:
        "pyscenic grn --seed {wildcards.run} --num_workers {threads} {input[0]} {input[1]} -o {output}"


rule multiple_prune_modules_rna_maskDropouts_with_cis_target_csv:
    input:
      exp=inputDataDir+"output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
      GRN_boost= "output/ScenicRNA_multipleRuns/GRNboost/run{run}_GRNboost.tsv",
      TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
      cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicRNA_multipleRuns/cis_target_maskDropouts/run{run}_regulons.csv"
    params: wdir="output/ScenicRNA/cis_target_maskDropouts/"
    threads:20
    conda:pyscenic_env
    shell:
        "cd {params.wdir};pyscenic ctx --annotations_fname {input.TF_motif} -a --mask_dropouts\
 --expression_mtx_fname {input.exp} -o ../../../{output} \
 ../../../{input.GRN_boost} {input.cis_motif}"
 
rule multiple_prune_modules_rna_maskDropouts_with_cis_target_csv_by_age:
    input:
      exp="output/ScenicRNA_multipleRuns_{age}/expressionRawCountFilteredForScenic.tsv",
      GRN_boost="output/ScenicRNA_multipleRuns_{age}/GRNboost/run{run}_GRNboost.tsv",
      TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
      cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/run{run}_regulons.csv"
    params: wdir="output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/"
    threads:20
    conda:pyscenic_env
    shell:
        "pyscenic ctx --annotations_fname {input.TF_motif} -a --mask_dropouts\
 --expression_mtx_fname {input.exp} -o {output} \
 {input.GRN_boost} {input.cis_motif}"

rule multiple_motifEnriched2regulons_rna_maskDropouts:
    input: "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/run{run}_regulons.csv"
    output: "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/run{run}_regulons.json"
    conda:pyscenic_env
    shell: "python py_src/motif2regulon.py -i {input} -o {output}"
    
rule multiple_motifEnriched2regulons_rna_maskDropouts_by_age:
    input: "output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/run{run}_regulons.csv"
    output: "output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/run{run}_regulons.json"
    conda:pyscenic_env
    shell: "python py_src/motif2regulon.py -i {input} -o {output}"

rule merging_multiple_runs:
    input: files = expand("output/ScenicRNA_multipleRuns/cis_target_maskDropouts/run{run}_regulons.json",run = scenicRuns),
           adjacencies =  expand("output/ScenicRNA_multipleRuns/GRNboost/run{run}_GRNboost.tsv",run = scenicRuns) 
    output: "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json",
            "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulonsMeta.json",
            "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/adja.sqlite"
    shell: "python py_src/aggregateMultiRun.py -r output/ScenicRNA_multipleRuns/cis_target_maskDropouts/ -c 0.8 -a output/ScenicRNA_multipleRuns/GRNboost/ -o output/ScenicRNA_multipleRuns/cis_target_maskDropouts/"
    
rule merging_multiple_runs_by_age:
    input: files = expand("output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/run{run}_regulons.json",run = scenicRuns,age =ages),
           adjacencies =  expand("output/ScenicRNA_multipleRuns_{age}/GRNboost/run{run}_GRNboost.tsv",run = scenicRuns,age =ages), 
    output: "output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/aggregatedRegulons.json",
            "output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/aggregatedRegulonsMeta.json",
            "output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/adja.sqlite"
    shell: "python py_src/aggregateMultiRun.py -r output/ScenicRNA_multipleRuns_{wildcards.age}/cis_target_maskDropouts/ -a output/ScenicRNA_multipleRuns_{wildcards.age}/GRNboost/ -c 0.8 -o output/ScenicRNA_multipleRuns_{wildcards.age}/cis_target_maskDropouts/"
    
rule multiple_aucell_TF_seurat_regulon_rna_maskDropouts:
    input:inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv", 
          "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json"
    output:"output/ScenicRNA_multipleRuns/AUCell_maskDropouts/regulons_enrichment.csv"

    threads:12
    conda:pyscenic_env
    shell:
        "python3.6 py_src/aucellJson.py -c {threads} -s 2020 -e {input[0]} -r {input[1]} -o output/ScenicRNA_multipleRuns/AUCell_maskDropouts/"
        
rule creating_regulonTable:
    input: main = "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulonsMeta.json",
           supp = expand("output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/aggregatedRegulonsMeta.json",age=ages)
    output: "output/regulonAnalysis/mainRegulonTable.tsv"
    params: condName = "young+old",
            tf_database = "output/publicData/mm_mgi_tfs.txt",
            supp = "+".join(expand("output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/aggregatedRegulonsMeta.json",age=ages))
    conda: seurat_monocle_env
    shell:
        "Rscript R_src/makeRegulonTableCL.R -o output/regulonAnalysis/ -t {params.tf_database} -r {input.main} -s {params.supp} -n {params.condName}"


rule makeInfluenceGraph:
    input: tfList = "input/selectedTF.txt",
           regulonTable = "output/regulonAnalysis/mainRegulonTable.tsv",
           cistromeBeta = "output/Cistrome_BM/cistromeReg.tsv",
           interactionBiblio = "input/interactionsReferences.txt"
    output: "output/Inference/influenceGraph/infGraphTable45.tsv"
    params: recovTimesThreshold = config["regulonAnalysis"]["recovTimesThresholdforInfGraph"],
    conda: seurat_monocle_env
    shell: "Rscript R_src/makeInfluenceGraphCL.R -t {input.tfList} -r {input.regulonTable} \
    -o output/Inference/influenceGraph/ -c -j {input.interactionBiblio} -a {input.cistromeBeta} -v {params.recovTimesThreshold}"
    
    
###################################################################################################################################
#### Report first part
rule get_report_influence_graph_disctretization:
    input:"report/figures_part1.Rmd",
          "output/regulonAnalysis/mainRegulonTable.tsv",
          "output/publicData/mm_mgi_tfs.txt",
          "input/selectedTF.txt",
          "output/ScenicRNA_multipleRuns/AUCell_maskDropouts/regulons_enrichment.csv",
          "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json",
          "output/ScenicRNA_multipleRuns_young/cis_target_maskDropouts/aggregatedRegulons.json",
          "output/ScenicRNA_multipleRuns_old/cis_target_maskDropouts/aggregatedRegulons.json",
          "output/regulonAnalysis/youngRegulonTable.tsv",
          "output/regulonAnalysis/oldRegulonTable.tsv",
          "output/Cistrome_BM/cistromeReg.tsv",
          "output/Inference/influenceGraph/infGraphTable45.tsv",
          inputDataDir+"/report/seurat_report.rds",
          inputDataDir+"/report/monocle_report.rds"
    output:"report/seuratObs.rds",
           "output/Inference/obsDataDis.csv",
           "report/tables/regulatorNet45.tsv",
           "report/tables/nodeTable45.tsv",
           "report/tables/regulonPopulationMarker.tsv",
           "report/tables/interactionTable.tsv",
           "report/figures_part1.html"
    threads: 12
    conda: report_env  
    shell: "echo $CONDA_DEFAULT_ENV;Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
###################################################################################################################################
   



rule installBonesis:
    output: "config/bonesis/install_done"
    params: gitUrl=config["bonesis"]["gitUrl"]
    conda: bonesis_env
    shell: 
      "cd config/;git clone {params.gitUrl};cd bonesis;pip install --user -e .;touch install_done"


rule bonesis_explore:
    input: graph =  "output/Inference/influenceGraph/infGraphTable45.tsv",
           obsData = "output/Inference/obsDataDis.csv",
           install = "config/bonesis/install_done"
    output: "output/Inference/bonesis/solution_space_first_exploration.html"
    conda: bonesis_env
    threads: 24
    shell:
      "jupyter nbconvert --ExecutePreprocessor.timeout=1000000 --to HTML --execute output/Inference/bonesis/solution_space_first_exploration.ipynb"

      
rule bonesis_optimize:
    input: graph =  "output/Inference/influenceGraph/infGraphTable45.tsv",
           obsData = "output/Inference/obsDataDis.csv",
           install = "config/bonesis/install_done"
    output: "output/Inference/bonesis/regulatory_graph_optimization.html",
            "output/Inference/bonesis/possible_final_solutions.p"
    conda: bonesis_env
    threads: 24
    shell:
      "jupyter nbconvert --ExecutePreprocessor.timeout=1000000 --to HTML --execute output/Inference/bonesis/regulatory_graph_optimization.ipynb" 

rule bonesis_final_solution:
    input: graph =  "output/Inference/influenceGraph/infGraphTable45.tsv",
           obsData = "output/Inference/obsDataDis.csv",
           install = "config/bonesis/install_done",
           sol = "output/Inference/bonesis/possible_final_solutions.p"
    output: "report/reportBonesis.html",
            "report/tables/interactionTableFinalSol.csv",
            "report/solutionFinal.zginml"
    conda: bonesis_env
    threads: 24
    shell:
      "jupyter nbconvert --ExecutePreprocessor.timeout=1000000 --to HTML --execute report/reportBonesis.ipynb" 

###################################################################################################################################
#### Report second part
rule get_final_report:
    input:"report/figures_part2.Rmd",
          "output/publicData/mm_mgi_tfs.txt",
          "input/selectedTF.txt",
          inputDataDir+"/report/seurat_report.rds",
          "report/tables/interactionTableFinalSol.csv",
          "report/seuratObs.rds"
    output:"report/tables/interactionRegulonTableFinalSol.tsv",
           "report/figures_part2.html",
           "report/tables/regulonPopulationAgingMarker.tsv"
    threads: 12
    conda: report_env  
    shell: "Rscript -e 'rmarkdown::render(\"{input[0]}\")'"
###################################################################################################################################
      
rule downloadRegDatabases:
    output: trrust = "output/publicData/trrust.tsv"
    params: trrust = config["databases"]["trrust"],
            cistrome = config["databases"]["cistrome"]
    shell:
      "wget -cO - {params.trrust} > {output.trrust}\
       wget -cO - {params.cistrome} > {output.cistrome}"

############################## Cistrome analysis  #########################################

rule filter_bed:
   input: "input/mouse_factor/{peakFile}_sort_peaks.narrowPeak.bed"
   output: "output/Cistrome/peaksFiltered/{peakFile}_sort_peaks.narrowPeak_filtered.bed"
   shell: "LC_ALL=C awk '{{ if($7 >= 5) {{ print }} }}' {input} > {output}"

rule BETA_score:
   input:  lambda wildcards: "output/Cistrome/peaksFiltered/"+expToFile[wildcards.gsm_factor]+"_sort_peaks.narrowPeak_filtered.bed"
   output: "output/Cistrome/BETA/{gsm_factor}_targets.txt"
   params: outdir = " output/Cistrome/BETA"
   conda: beta_env
   shell: "BETA minus -p {input} -g mm10 -d 10000 -o {params.outdir} -n {wildcards.gsm_factor}"

rule filter_bed_bone_marrow:
   input: "input/mouse_factor/{peakFile}_sort_peaks.narrowPeak.bed"
   output: "output/Cistrome_BM/peaksFiltered/{peakFile}_sort_peaks.narrowPeak_filtered.bed"
   shell: "LC_ALL=C awk '{{ if($7 >= 5) {{ print }} }}' {input} > {output}"

rule BETA_score_bone_marrow:
   input:  lambda wildcards: "output/Cistrome_BM/peaksFiltered/"+expToFileTissue[wildcards.gsm_factor]+"_sort_peaks.narrowPeak_filtered.bed"
   output: "output/Cistrome_BM/BETA/{gsm_factor}_targets.txt"
   params: outdir = " output/Cistrome_BM/BETA"
   conda: beta_env
   shell: "BETA minus -p {input} -g mm10 -d 10000 -o {params.outdir} -n {wildcards.gsm_factor}"

rule BETA_score_aggregation:
   input: expand("output/Cistrome_BM/BETA/{gsm_factor}_targets.txt",gsm_factor = expToFileTissue.keys())
   output: "output/Cistrome_BM/cistromeReg.tsv"
   conda: seurat_monocle_env
   shell: "Rscript R_src/cistromeAnalysisCL.R -b output/Cistrome_BM/BETA/ -o output/Cistrome_BM/"

