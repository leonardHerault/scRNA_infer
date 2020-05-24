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

os.chdir(baseDir)


rule all:
    input:"reportProjet2/report_final.html"
    
rule getDataMatrixFromPreviousWork:
    input: inputDataDir+"/report/seurat_report.rds",
    output: "input/dataMatrix.csv"
    conda: seurat_monocle_env
    params: slot = "data"
    threads: 1
    shell: "Rscript R_src/getDataMatrixCL.R -i {input} -s {params.slot} -o input/"

    

rule dca:
    input: "input/dataMatrix.csv"
    output: "output/dca/results/mean.tsv"
    threads: 24 
    conda: dca_env
    params : config["dca"]
    shell: "dca {params} {input} output/dca/results"
    
# rule dca_hyper:
#     input: inputDataDir+"report/dca/dataMatrix.csv"
#     output: "output/dca/hyper/search_done"
#     threads: 24 
#     conda: dca_env
#     shell: "dca --hyper {input} output/dca/hyper --hyper  &> output/dca/hyper;touch {output}"
    
    
rule filter_genes_for_scenic_dca:
    input: "output/dca/results/mean.tsv"
    output: "output/ScenicDCA/expressionRawCountFilteredForScenic.tsv"
    conda: seurat_monocle_env
    shell: "Rscript R_src/filterForDcaScenicCL.R -i {input} -o output/ScenicDCA/"

rule prepare_tf_list:
     input: inputDataDir+"/publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
     output: "output/publicData/mm_mgi_tfs.txt"
     conda: pyscenic_env
     shell: "python py_src/getTfList.py -i {input} -o {output}"

        
rule grnboost_rna:
    input: inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
            "output/publicData/mm_mgi_tfs.txt"
    output: "output/ScenicRNA/GRNboost/GRNboost.tsv"
    threads: 20
    conda: pyscenic_env
    shell:
        "pyscenic grn --num_workers {threads} {input[0]} {input[1]} -o {output}"

        
rule grnboost_dca:
    input: "output/ScenicDCA/expressionRawCountFilteredForScenic.tsv",
            "output/publicData/mm_mgi_tfs.txt"
    output: "output/ScenicDCA/GRNboost/GRNboost.tsv"
    threads: 26
    conda: pyscenic_env
    shell:
        "pyscenic grn --num_workers {threads} {input[0]} {input[1]} -o {output}"
        
rule prune_modules_dca_with_cis_target_csv:
    input:
        exp="output/ScenicDCA/expressionRawCountFilteredForScenic.tsv",
        GRN_boost="output/ScenicDCA/GRNboost/GRNboost.tsv",
        TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
        cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicDCA/cis_target/regulons.csv"
    threads:20
    conda:pyscenic_env
    shell:
        "pyscenic ctx --annotations_fname {input.TF_motif} -a --min_genes 15 \
--expression_mtx_fname {input.exp} -o {output} \
 {input.GRN_boost} {input.cis_motif}"
 
rule prune_modules_rna_with_cis_target_csv:
    input:
      exp=inputDataDir+"output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
      GRN_boost="output/ScenicRNA/GRNboost/GRNboost.tsv",
      TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
      cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicRNA/cis_target/regulons.csv"
    threads:20
    conda:pyscenic_env
    shell:
        "pyscenic ctx --annotations_fname {input.TF_motif} -a --min_genes 15 \
--expression_mtx_fname {input.exp} -o {output} \
 {input.GRN_boost} {input.cis_motif}"


rule prune_modules_dca_with_cis_target_json:
    input:
        exp="output/ScenicDCA/expressionRawCountFilteredForScenic.tsv",
        GRN_boost="output/ScenicDCA/GRNboost/GRNboost.tsv",
        TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
        cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicDCA/cis_target/regulons.json"
    threads:20
    conda:pyscenic_env
    shell:
        "pyscenic ctx --annotations_fname {input.TF_motif} --min_genes 15 \
--expression_mtx_fname {input.exp} -a -o {output} \
 {input.GRN_boost} {input.cis_motif}"
 
 
rule prune_modules_rna_with_cis_target_json:
    input:
        exp="output/ScenicAll/expressionRawCountFilteredForScenic.tsv",
        GRN_boost="output/ScenicRNA/GRNboost/GRNboost.tsv",
        TF_motif=inputDataDir+"publicData/database/motifs-v9-nr.mgi-m0.001-o0.0.tbl",
        cis_motif=inputDataDir+"publicData/database/mm9-tss-centered-10kb-7species.mc9nr.feather"
    output:"output/ScenicRNA/cis_target/regulons.json"
    threads:20
    conda:pyscenic_env
    shell:
        "pyscenic ctx --annotations_fname {input.TF_motif} --min_genes 15 \
--expression_mtx_fname {input.exp} -a -o {output} \
 {input.GRN_boost} {input.cis_motif}"
 
 
        
rule aucell_TF_seurat_regulon_dca:
    input:"output/ScenicDCA/expressionRawCountFilteredForScenic.tsv", 
          "output/ScenicDCA/cis_target/regulons.csv"
    output:"output/ScenicDCA/AUCell/regulons_enrichment.csv"
    threads:20
    conda:pyscenic_env
    shell:
        "pyscenic aucell {input[0]} {input[1]} -o {output}"
        
        
rule aucell_TF_seurat_regulon_rna:
    input:inputDataDir+"/output/ScenicAll/expressionRawCountFilteredForScenic.tsv", 
          "output/ScenicRNA/cis_target/regulons.csv"
    output:"output/ScenicRNA/AUCell/regulons_enrichment.csv"
    threads:20
    conda:pyscenic_env
    shell:
        "pyscenic aucell {input[0]} {input[1]} -o {output}"



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
#            "output/ScenicDCA/cis_target/regulons.json"
#     output: "output/ScenicDCA/tf_grn_regulons.csv"
#     threads:1
#     conda: pyscenic_env
#     shell: "python py_src/getTFRN.py -t {input[0]} -r {input[1]} -o {output}"
#     
# rule get_tf_regulatory_networks_TF_seurat_regulon_rna:
#     input:  "output/publicData/mm_mgi_tfs.txt",
#            "output/ScenicRNA/cis_target/regulons.json"
#     output: "output/ScenicRNA/tf_grn_regulons.csv"
#     threads:1
#     conda: pyscenic_env
#     shell: "python py_src/getTFRN.py -t {input[0]} -r {input[1]} -o {output}"

    

rule add_dca_to_Seurat:
    input: seurat =inputDataDir+"/report/seurat_report.rds",
           sig = inputDataDir+"/output/signatures/publicSignatures.rds",
           rodriguez = inputDataDir+"/publicData/rodriguez_cluster_signatures.xlsx",
           dcaMean = "output/dca/results/mean.tsv"
    output: "output/dca/Seurat3/seuratDca.rds",
            "output/dca/Seurat3/dca_norm.csv"
    conda: seurat_monocle_env
    params: ident = "numclust"
    threads: 1
    shell: "Rscript R_src/addDcaCL.R -i {input.seurat} -j {input.dcaMean} -o output/dca/Seurat3 -a {params.ident} -s {input.sig} -k {input.rodriguez}"

rule progeny:
    input: normData = "input/dataMatrix.csv"
    output: "output/progeny/rna/progeny_scores.tsv"
    conda: progeny_env
    threads:1
    shell: "Rscript R_src/progenyCL.R -i {input.normData} -o output/progeny/rna"
    
rule progeny_dca:
    input: normData = "output/dca/Seurat3/dca_norm.csv"
    output: "output/progeny/dca/progeny_scores.tsv"
    conda: progeny_env
    threads:1
    shell: "Rscript R_src/progenyCL.R -i {input.normData} -o output/progeny/dca"
    
    
rule dorothea:
    input: seurat = inputDataDir+"report/seurat_report.rds",
    output:  "output/dorothea/rna/seuratDorothea.rds"
    conda: dorothea_env
    threads:10
    shell: "Rscript R_src/DorotheaCL.R -c {threads} -i {input.seurat} -o output/dorothea/rna"
    

rule dorothea_dca:
    input: seurat = "output/dca/Seurat3/seuratDca.rds",
           dorotheaCorrectlyInstalled ="output/dorothea/rna/seuratDorothea.rds"
    output: "output/dorothea/dca/seuratDorothea.rds"
    conda: dorothea_env
    threads:20
    shell: "Rscript R_src/DorotheaCL.R -c {threads} -i {input.seurat} -o output/dorothea/dca"
    
rule get_reportProjet2:
    input:#"output/dorothea/dca/seuratDorothea.rds",
          #"output/dorothea/rna/seuratDorothea.rds",
          "output/ScenicDCA/AUCell/regulons_enrichment.csv",
          "output/ScenicRNA/AUCell/regulons_enrichment.csv",
          "output/progeny/dca/progeny_scores.tsv",
          "output/progeny/rna/progeny_scores.tsv",
          #"output/ScenicRNA/tf_grn_regulons.csv",
          #"output/ScenicDCA/tf_grn_regulons.csv",
          #"output/dca/hyper/search_done"
    output:"reportProjet2/report_final.html"
    shell: "touch {output}"
