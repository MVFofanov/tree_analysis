configfile: "/home/zo49sog/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary/tree_analysis_test/config/config.yaml"

import os
from glob import glob
import yaml

# Load the configuration file
config_file = "/home/zo49sog/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary/tree_analysis_test/config/config.yaml"

with open(config_file) as file:
    config = yaml.safe_load(file)

# Directories and paths from the config
working_dir = config["input"]["deni_data"]
tree_leaves = config["input"]["tree_leaves"].format(deni_data=working_dir)
wd = config["input"]["wd"].format(tree_leaves=tree_leaves)
phylogenetic_trees_dir = config["input"]["phylogenetic_trees_dir"].format(deni_data=working_dir)
annotation_file_id = config["input"]["annotation_file_id"].format(tree_leaves=tree_leaves)
config_dir = config["input"]["config_dir"].format(wd=wd)
clusters_file = config["input"]["clusters_file"].format(config_dir=config_dir)

base_output_dir = config["output"]["base_output_dir"].format(wd=wd)
output_dir = config["output"]["base_output_dir"].format(wd=wd)
logs_dir = config["output"]["logs_dir"].format(base_output_dir=output_dir)

# Read cluster names from the clusters.txt file
with open(clusters_file) as f:
    cluster_names = [line.strip() for line in f.readlines() if line.strip()]

# Main rule to request all outputs
rule all:
    input:
        expand(f"{base_output_dir}/{{cluster}}/rooted/{{cluster}}_log_tree_analysis.log", cluster=cluster_names),
        f"{base_output_dir}/cluster_analysis/rooted/comparison_complete.log"

# Rule for processing individual clusters
rule process_cluster:
    input:
        config="/home/zo49sog/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary/tree_analysis_test/config/config.yaml"
    output:
        log_file=f"{base_output_dir}/{{cluster}}/rooted/{{cluster}}_log_tree_analysis.log",
        biggest_clades=f"{base_output_dir}/{{cluster}}/rooted/biggest_non_intersecting_clades_all.tsv",
        tree=f"{base_output_dir}/{{cluster}}/rooted/annotated_tree.pdf"
    threads: 1
    shell:
        """
        mkdir -p {output_dir}/{wildcards.cluster}
        python3 /home/zo49sog/crassvirales/phylomes/scripts/tree_analysis/main.py --cluster {wildcards.cluster} --config {input.config}
        """

# Rule for comparing clusters after processing all clusters
rule compare_clusters:
    input:
        expand(f"{base_output_dir}/{{cluster}}/rooted/biggest_non_intersecting_clades_all.tsv", cluster=cluster_names),
        config="/home/zo49sog/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/phylome_summary/tree_analysis_test/config/config.yaml",
        clusters_file=clusters_file
    output:
        final_log=f"{base_output_dir}/cluster_analysis/rooted/comparison_complete.log",
        concatenated_clusters=f"{base_output_dir}/cluster_analysis/rooted/concatenated_clusters_data.tsv"
    threads: 1
    shell:
        """
        python3 /home/zo49sog/crassvirales/phylomes/scripts/tree_analysis/cluster_comparison.py --config {input.config} --clusters_file "{input.clusters_file}" > {output.final_log}
        """
