configfile: "/home/zo49sog/crassvirales/phylomes/tree_analysis/config/config.yaml"

import os
from glob import glob
import yaml

# Load the configuration file
config_file = "/home/zo49sog/crassvirales/phylomes/tree_analysis/config/config.yaml"

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

# Define tree types
tree_types = ["rooted", "unrooted", "midpoint"]

# Main rule to request all outputs for both rooted and unrooted trees
rule all:
    input:
        expand(f"{base_output_dir}/{{cluster}}/{{tree_type}}/{{cluster}}_log_tree_analysis.log",
               cluster=cluster_names, tree_type=tree_types),
        expand(f"{base_output_dir}/cluster_analysis/{{tree_type}}/comparison_complete.log", tree_type=tree_types)

# Rule for processing individual clusters
rule process_cluster:
    input:
        config=config_file,
        clusters_file=clusters_file
    output:
        log_file=f"{base_output_dir}/{{cluster}}/{{tree_type}}/{{cluster}}_log_tree_analysis.log",
        biggest_clades=f"{base_output_dir}/{{cluster}}/{{tree_type}}/biggest_non_intersecting_clades_all.tsv",
        tree=f"{base_output_dir}/{{cluster}}/{{tree_type}}/annotated_tree.pdf"
    params:
        tree_type=lambda wildcards: wildcards.tree_type  # Handle both rooted and unrooted
    threads: 1
    shell:
        """
        mkdir -p {output_dir}/{wildcards.cluster}/{wildcards.tree_type}
        source /home/zo49sog/mambaforge/etc/profile.d/conda.sh && conda activate tree_analysis
        python3 /home/zo49sog/crassvirales/phylomes/tree_analysis/scripts/main.py --cluster {wildcards.cluster} --config {input.config}
        """

# Rule for comparing clusters after processing all clusters (rooted and unrooted)
rule compare_clusters:
    input:
        expand(f"{base_output_dir}/{{cluster}}/{{tree_type}}/biggest_non_intersecting_clades_all.tsv",
               cluster=cluster_names, tree_type=tree_types),
        config=config_file,
        clusters_file=clusters_file
    output:
        final_log=f"{base_output_dir}/cluster_analysis/{{tree_type}}/comparison_complete.log",
        concatenated_clusters=f"{base_output_dir}/cluster_analysis/{{tree_type}}/concatenated_clusters_data.tsv"
    threads: 1
    params:
        tree_type=lambda wildcards: wildcards.tree_type
    shell:
        """
        source /home/zo49sog/mambaforge/etc/profile.d/conda.sh && conda activate tree_analysis
        python3 /home/zo49sog/crassvirales/phylomes/tree_analysis/scripts/cluster_comparison.py --config {input.config} --clusters_file "{input.clusters_file}" > {output.final_log}
        """
