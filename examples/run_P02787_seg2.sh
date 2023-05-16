#!/bin/sh

mkdir -p ~/EMBL-EBI/funclan_work/P02787_2/ca_distances
mkdir -p ~/EMBL-EBI/funclan_work/P02787_2/distance_differences
mkdir -p ~/EMBL-EBI/funclan_work/P02787_2/distance_difference_maps
mkdir -p ~/EMBL-EBI/funclan_work/P02787_2/cluster_results
# Remove files from previous executions
rm ~/EMBL-EBI/funclan_work/P02787_2/ca_distances/*
rm ~/EMBL-EBI/funclan_work/P02787_2/distance_differences/*
rm ~/EMBL-EBI/funclan_work/P02787_2/distance_difference_maps/*
rm ~/EMBL-EBI/funclan_work/P02787_2/cluster_results/*

python3 run_find_conformers.py -u "P02787" \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1a8e_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1a8f_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1b3e_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1bp5_updated.cif A B C D \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1btj_updated.cif A B \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1d3k_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1d4n_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1dtg_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1fqe_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1fqf_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1jqf_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1n7w_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1n7x_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1n84_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1oqg_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1oqh_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1ryo_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/1suv_updated.cif A B C D \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/2o7u_updated.cif B A C D E F G H I \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/2o84_updated.cif X \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3fgs_updated.cif A \
    -c ~/EMBL-EBI/funclan_work/P02787_2/ca_distances/ \
    -d ~/EMBL-EBI/funclan_work/P02787_2/distance_differences/ \
    -s ~/EMBL-EBI/funclan_work/P02787_2/cluster_results/ \
    -g ~/EMBL-EBI/funclan_work/P02787_2/cluster_results/ png svg

# -o ~/EMBL-EBI/funclan_work/P02787/distance_difference_maps/ \
