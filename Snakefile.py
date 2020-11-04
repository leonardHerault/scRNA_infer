import os
import re
import pandas as pd

baseDir = config["dir"]["base"]
inputDataDir = config["dir"]["inputData"]

seurat_monocle_env = config["seurat_monocle_env"]
pyscenic_env = config["pyscenic_env"]
bonesis_env = config["bonesis_env"]

dca_env = config["dca_env"]
progeny_env = config["progeny_env"]
dorothea_env = config["dorothea_env"]
ages = ["young","old"]

scenicRuns = list(range(50))


rule all:
    input:"reportProjet2/report_final.html",
          "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json",
          "output/ScenicRNA_multipleRuns/AUCell_maskDropouts/regulons_enrichment.csv",
          "output/regulonAnalysis/mainRegulonTable.tsv",
          "output/publicData/trrust.tsv",
          "output/regulonAnalysis/clusterMarkerRegulonTable.tsv",
          "output/regulonAnalysis/infGraphTable.tsv",
          "output/Inference/bonesis/install_done",
          "output/Inference/bonesis/solutions.p"
    
rule getDataMatrixFromPreviousWork:
    input: inputDataDir+"/report/seurat_report.rds",
    output: "input/dataMatrix.csv"
    conda: seurat_monocle_env
    params: slot = "data"
    threads: 1
    shell: "Rscript R_src/getDataMatrixCL.R -i {input} -s {params.slot} -o input/"

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

rule analyseRegulonScore:
    input: seurat =  inputDataDir+"/report/seurat_report.rds",
           tfList = "input/selectedTF.txt",
           regulonScore =  "output/ScenicRNA_multipleRuns/AUCell_maskDropouts/regulons_enrichment.csv"
    output: "output/regulonAnalysis/clusterMarkerRegulonTable.tsv"
    params: scoreDiff = config["regulonAnalysis"]["scoreDiff"]
    conda: seurat_monocle_env
    shell:
      "Rscript R_src/regulonMarkerCL.R -i {input.seurat} -t {input.tfList} \
      -r {input.regulonScore} -o output/regulonAnalysis/ -c -a AGE"

rule makeInfluenceGraph:
    input: tfList = "input/selectedTF.txt",
           regulonTable = "output/regulonAnalysis/mainRegulonTable.tsv",
           biblioNet = "input/KrumsiekAdapted.reggraph",
    output: "output/regulonAnalysis/infGraphTable.tsv"
    params: recovTimesThreshold = config["regulonAnalysis"]["recovTimesThresholdforInfGraph"]
    conda: seurat_monocle_env
    shell: "Rscript R_src/makeInfluenceGraphCL.R -t {input.tfList} -r {input.regulonTable} \
    -o output/regulonAnalysis/ -c -b {input.biblioNet} -v {params.recovTimesThreshold}"

rule installBonesis:
    output: "output/Inference/bonesis/install_done"
    params: gitUrl=config["bonesis"]["gitUrl"]
    conda: bonesis_env
    shell: 
      "cd output/Inference/;git clone {params.gitUrl};cd bonesis;pip install --user -e .;touch install_done"

rule bonesis:
    input: "output/regulonAnalysis/infGraphTable.tsv"
    output: "output/Inference/bonesis/solutions.p"
    conda: bonesis_env
    threads: 20
    shell: 
      "python py_src/inferMpbn.py -o output/Inference/bonesis/ -c {threads} -i {input}" 

      
rule downloadRegDatabases:
    output: trrust = "output/publicData/trrust.tsv"
    params: trrust = config["databases"]["trrust"]
    shell:
      "wget -cO - {params.trrust} > {output.trrust}"
  

rule progeny:
    input: normData = "input/dataMatrix.csv"
    output: "output/progeny/rna/progeny_scores.tsv"
    conda: progeny_env
    threads:1
    shell: "Rscript R_src/progenyCL.R -i {input.normData} -o output/progeny/rna"
  
rule get_reportProjet2:
    input:"output/ScenicRNA/AUCell/regulons_enrichment.csv",
          "output/ScenicRNA/AUCell_maskDropouts/regulons_enrichment.csv",
          "output/progeny/rna/progeny_scores.tsv",
          "output/ScenicRNA/tf_grn_regulons.csv",
          "output/ScenicRNA/tf_grn_regulons_maskDropouts.csv",
    output:"reportProjet2/report_final.html"
    shell: "touch {output}"
