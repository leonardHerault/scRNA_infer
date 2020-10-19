import os
import re
import pandas as pd

baseDir = config["dir"]["base"]
inputDataDir = config["dir"]["inputData"]

seurat_monocle_env = config["seurat_monocle_env"]
pyscenic_env = config["pyscenic_env"]
dca_env = config["dca_env"]
progeny_env = config["progeny_env"]
dorothea_env = config["dorothea_env"]
ages = ["young","old"]

scenicRuns = list(range(50))


rule all:
    input:"reportProjet2/report_final.html",
          "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json",
          expand("output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/aggregatedRegulons.json",age = ages),
          "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/strict/aggregatedRegulons.json",
          expand("output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/strict/aggregatedRegulons.json",age = ages),
          "output/ScenicRNA_multipleRuns/AUCell_maskDropouts/regulons_enrichment.csv"
    
rule getDataMatrixFromPreviousWork:
    input: inputDataDir+"/report/seurat_report.rds",
    output: "input/dataMatrix.csv"
    conda: seurat_monocle_env
    params: slot = "data"
    threads: 1
    shell: "Rscript R_src/getDataMatrixCL.R -i {input} -s {params.slot} -o input/"

    

# rule dca:
#     input: "input/dataMatrix.csv"
#     output: "output/dca/results/mean.tsv"
#     threads: 24 
#     conda: dca_env
#     params : config["dca"]
#     shell: "dca {params} {input} output/dca/results"
    
# rule dca_hyper:
#     input: inputDataDir+"report/dca/dataMatrix.csv"
#     output: "output/dca/hyper/search_done"
#     threads: 24 
#     conda: dca_env
#     shell: "dca --hyper {input} output/dca/hyper --hyper  &> output/dca/hyper;touch {output}"
    
    
# rule filter_genes_for_scenic_dca:
#     input: "output/dca/results/mean.tsv"
#     output: "output/ScenicDCA/expressionRawCountFilteredForScenic.tsv"
#     conda: seurat_monocle_env
#     shell: "Rscript R_src/filterForDcaScenicCL.R -i {input} -o output/ScenicDCA/"

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

rule grnboost_rna:
    input: inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
            "output/publicData/mm_mgi_tfs.txt"
    output: "output/ScenicRNA/GRNboost/GRNboost.tsv"
    threads: 20
    conda: pyscenic_env
    shell:
        "pyscenic grn --seed 2020 --num_workers {threads} {input[0]} {input[1]} -o {output}"
        
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
    output: "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json"
    shell: "python py_src/aggregateMultiRun.py -r output/ScenicRNA_multipleRuns/cis_target_maskDropouts/ -c 0.8 -a output/ScenicRNA_multipleRuns/GRNboost/ -o output/ScenicRNA_multipleRuns/cis_target_maskDropouts/"
    
rule merging_multiple_runs_by_age:
    input: files = expand("output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/run{run}_regulons.json",run = scenicRuns,age =ages),
           adjacencies =  expand("output/ScenicRNA_multipleRuns_{age}/GRNboost/run{run}_GRNboost.tsv",run = scenicRuns,age =ages), 
    output: "output/ScenicRNA_multipleRuns_{age}/cis_target_maskDropouts/aggregatedRegulons.json"
    shell: "python py_src/aggregateMultiRun.py -r output/ScenicRNA_multipleRuns_{wildcards.age}/cis_target_maskDropouts/ -a output/ScenicRNA_multipleRuns_{wildcards.age}/GRNboost/ -c 0.8 -o output/ScenicRNA_multipleRuns_{wildcards.age}/cis_target_maskDropouts/"
    
rule multiple_aucell_TF_seurat_regulon_rna_maskDropouts:
    input:inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv", 
          "output/ScenicRNA_multipleRuns/cis_target_maskDropouts/aggregatedRegulons.json"
    output:"output/ScenicRNA_multipleRuns/AUCell_maskDropouts/regulons_enrichment.csv"

    threads:12
    conda:pyscenic_env
    shell:
        "python3.6 py_src/aucellJson.py -c {threads} -s 2020 -e {input[0]} -r {input[1]} -o output/ScenicRNA_multipleRuns/AUCell_maskDropouts/"

        
# rule grnboost_dca:
#     input: "output/ScenicDCA/expressionRawCountFilteredForScenic.tsv",
#             "output/publicData/mm_mgi_tfs.txt"
#     output: "output/ScenicDCA/GRNboost/GRNboost.tsv"
#     threads: 26
#     conda: pyscenic_env
#     shell:
#         "pyscenic grn --num_workers {threads} {input[0]} {input[1]} -o {output}"
        
# rule prune_modules_dca_with_cis_target_csv:
#     input:
#         exp="output/ScenicDCA/expressionRawCountFilteredForScenic.tsv",
#         GRN_boost="output/ScenicDCA/GRNboost/GRNboost.tsv",
#         TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
#         cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
#     output:"output/ScenicDCA/cis_target/regulons.csv"
#     params: wdir="cenicDCA/cis_target"
#     threads:15
#     conda:pyscenic_env
#     shell:
#         "cd {params.wdir};pyscenic ctx  --num_workers {threads}  --annotations_fname {input.TF_motif} -a\
# --expression_mtx_fname ../../{input.exp} -o ../../{output} \
#  ../../{input.GRN_boost} {input.cis_motif}"
 
rule prune_modules_rna_with_cis_target_csv:
    input:
      exp=inputDataDir+"output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
      GRN_boost="output/ScenicRNA/GRNboost/GRNboost.tsv",
      TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
      cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicRNA/cis_target/regulons.csv"
    params: wdir="output/ScenicRNA/cis_target/"
    threads:20
    conda:pyscenic_env
    shell:
        "cd {params.wdir};pyscenic ctx --annotations_fname {input.TF_motif} -a \
--expression_mtx_fname {input.exp} -o ../../../{output} \
 ../../../{input.GRN_boost} {input.cis_motif}"
 
rule prune_modules_rna_maskDropouts_with_cis_target_csv:
    input:
      exp=inputDataDir+"output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
      GRN_boost="output/ScenicRNA/GRNboost/GRNboost.tsv",
      TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
      cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicRNA/cis_target_maskDropouts/regulons.csv"
    params: wdir="output/ScenicRNA/cis_target_maskDropouts/"
    threads:20
    conda:pyscenic_env
    shell:
        "cd {params.wdir};pyscenic ctx --annotations_fname {input.TF_motif} -a --mask_dropouts\
 --expression_mtx_fname {input.exp} -o ../../../{output} \
 ../../../{input.GRN_boost} {input.cis_motif}"


# rule prune_modules_dca_with_cis_target_json:
#     input:
#         exp="output/ScenicDCA/expressionRawCountFilteredForScenic.tsv",
#         GRN_boost="output/ScenicDCA/GRNboost/GRNboost.tsv",
#         TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
#         cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
#     output:"output/ScenicDCA/cis_target/regulons.json"
#     threads:20
#     conda:pyscenic_env
#     shell:
#         "pyscenic ctx --annotations_fname {input.TF_motif} --min_genes 15 \
# --expression_mtx_fname {input.exp} -a -o {output} \
#  {input.GRN_boost} {input.cis_motif}"
#  
 
rule motifEnriched2regulons_rna:
    input: "output/ScenicRNA/cis_target/regulons.csv"
    output: "output/ScenicRNA/cis_target/regulons.json"
    conda:pyscenic_env
    shell: "python py_src/motif2regulon.py -i {input} -o {output}"
    
rule motifEnriched2regulons_rna_maskDropouts:
    input: "output/ScenicRNA/cis_target_maskDropouts/regulons.csv"
    output: "output/ScenicRNA/cis_target_maskDropouts/regulons.json"
    conda:pyscenic_env
    shell: "python py_src/motif2regulon.py -i {input} -o {output}"
    
# rule motifEnriched2regulons_dca:
#     input: "output/ScenicDCA/cis_target/regulons.csv"
#     output: "output/ScenicDCA/cis_target/regulons.json"
#     conda:pyscenic_env
#     shell: "python py_src/motif2regulon.py -i {input} -o {output}"
 


# rule aucell_TF_seurat_regulon_dca:
#     input:"output/ScenicDCA/expressionRawCountFilteredForScenic.tsv", 
#           "output/ScenicDCA/cis_target/regulons.csv"
#     output:"output/ScenicDCA/AUCell/regulons_enrichment.csv"
#     threads:20
#     conda:pyscenic_env
#     shell:
#         "pyscenic aucell --seed 2020 {input[0]} {input[1]} -o {output}"
        
        
rule aucell_TF_seurat_regulon_rna:
    input:inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv", 
          "output/ScenicRNA/cis_target/regulons.csv"
    output:"output/ScenicRNA/AUCell/regulons_enrichment.csv"
    params: wdir="output/ScenicRNA/AUCell/"
    threads:6
    conda:pyscenic_env
    shell:
        "cd {params.wdir};pyscenic aucell --num_workers {threads} --seed 2020 {input[0]} ../../../{input[1]} -o ../../../{output}"
        
rule aucell_TF_seurat_regulon_rna_maskDropouts:
    input:inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv", 
          "output/ScenicRNA/cis_target_maskDropouts/aggregatedRegulons.json"
    output:"output/ScenicRNA/AUCell_maskDropouts/regulons_enrichment.csv"
    params: wdir="output/ScenicRNA/AUCell_maskDropouts/"
    threads:6
    conda:pyscenic_env
    shell:
        "cd {params.wdir};pyscenic aucell --num_workers {threads} --seed 2020 {input[0]} ../../../{input[1]} -o ../../../{output}"



# rule regulons_in_traj_TF_seurat_regulon_dca:
#     input: "report/monocle_report.rds",
#             "output/ScenicDCA/AUCell/regulons_enrichment.csv"
#     output:"output/ScenicDCA/regulons_in_trajectory/distrib_AUC_sum.png",
#            "output/ScenicDCA/regulons_in_trajectory/monocleWithRegulons.rds"
#     #params:mart = config["biomaRt"]["mmusculus"]
#     threads:12
#     conda:seurat_monocle_env
#     shell:"Rscript R_src/regulonsCL.R -i {input[0]} -o output/ScenicDCA/regulons_in_trajectory/ -r {input[1]}"

# rule get_tf_regulatory_networks_TF_seurat_regulon_dca:
#     input:  "output/publicData/mm_mgi_tfs.txt",
#             "output/ScenicDCA/cis_target/regulons.json"
#     output: "output/ScenicDCA/tf_grn_regulons.csv"
#     threads:1
#     conda: pyscenic_env
#     shell: "python py_src/getTFRN.py -t {input[0]} -r {input[1]} -o {output}"
     
rule get_tf_regulatory_networks_TF_seurat_regulon_rna:
    input:  "output/publicData/mm_mgi_tfs.txt",
            "output/ScenicRNA/cis_target/regulons.json"
    output: "output/ScenicRNA/tf_grn_regulons.csv"
    threads:1
    conda: pyscenic_env
    shell: "python py_src/getTFRN.py -t {input[0]} -r {input[1]} -o {output}"

rule get_tf_regulatory_networks_TF_seurat_regulon_rna_maskDropouts:
    input:  "output/publicData/mm_mgi_tfs.txt",
            "output/ScenicRNA/cis_target_maskDropouts/regulons.json"
    output: "output/ScenicRNA/tf_grn_regulons_maskDropouts.csv"
    threads:1
    conda: pyscenic_env
    shell: "python py_src/getTFRN.py -t {input[0]} -r {input[1]} -o {output}"
    

# rule add_dca_to_Seurat:
#     input: seurat =inputDataDir+"/report/seurat_report.rds",
#            sig = inputDataDir+"/output/signatures/publicSignatures.rds",
#            rodriguez = inputDataDir+"/publicData/rodriguez_cluster_signatures.xlsx",
#            dcaMean = "output/dca/results/mean.tsv"
#     output: "output/dca/Seurat3/seuratDca.rds",
#             "output/dca/Seurat3/dca_norm.csv"
#     conda: seurat_monocle_env
#     params: ident = "numclust"
#     threads: 1
#     shell: "Rscript R_src/addDcaCL.R -i {input.seurat} -j {input.dcaMean} -o output/dca/Seurat3 -a {params.ident} -s {input.sig} -k {input.rodriguez}"

rule progeny:
    input: normData = "input/dataMatrix.csv"
    output: "output/progeny/rna/progeny_scores.tsv"
    conda: progeny_env
    threads:1
    shell: "Rscript R_src/progenyCL.R -i {input.normData} -o output/progeny/rna"
    
# rule progeny_dca:
#     input: normData = "output/dca/Seurat3/dca_norm.csv"
#     output: "output/progeny/dca/progeny_scores.tsv"
#     conda: progeny_env
#     threads:1
#     shell: "Rscript R_src/progenyCL.R -i {input.normData} -o output/progeny/dca"
    
    
# rule dorothea:
#     input: seurat = inputDataDir+"report/seurat_report.rds"
#     output:  "output/dorothea/rna/seuratDorothea.rds"
#     conda: dorothea_env
#     threads:2
#     shell: "Rscript R_src/DorotheaCL.R -c {threads} -i {input.seurat} -o output/dorothea/rna"

# rule dorothea_all_levels:
#     input: seurat = inputDataDir+"report/seurat_report.rds",
#            dorotheaCorrectlyInstalled ="output/dorothea/rna/seuratDorothea.rds"
#     output:  "output/dorothea_all_levels/rna/seuratDorothea.rds"
#     conda: dorothea_env
#     threads:2
#     shell: "Rscript R_src/DorotheaCL.R -l A_B_C_D_E -c {threads} -i {input.seurat} -o output/dorothea_all_levels/rna"

# rule dorothea_dca:
#     input: seurat = "output/dca/Seurat3/seuratDca.rds",
#            dorotheaCorrectlyInstalled ="output/dorothea/rna/seuratDorothea.rds"
#     output: "output/dorothea/dca/seuratDorothea.rds"
#     conda: dorothea_env
#     threads:2
#     shell: "Rscript R_src/DorotheaCL.R -a DCA -c {threads} -i {input.seurat} -o output/dorothea/dca"
  

# rule dorothea_dca_all_levels:
#     input: seurat = "output/dca/Seurat3/seuratDca.rds",
#            dorotheaCorrectlyInstalled ="output/dorothea/rna/seuratDorothea.rds"
#     output: "output/dorothea_all_levels/dca/seuratDorothea.rds"
#     conda: dorothea_env
#     threads:2
#     shell: "Rscript R_src/DorotheaCL.R -a DCA -l A_B_C_D_E -c {threads} -i {input.seurat} -o output/dorothea_all_levels/dca"
    
rule get_reportProjet2:
    input:"output/ScenicRNA/AUCell/regulons_enrichment.csv",
          "output/ScenicRNA/AUCell_maskDropouts/regulons_enrichment.csv",
          "output/progeny/rna/progeny_scores.tsv",
          "output/ScenicRNA/tf_grn_regulons.csv",
          "output/ScenicRNA/tf_grn_regulons_maskDropouts.csv",
    output:"reportProjet2/report_final.html"
    shell: "touch {output}"
