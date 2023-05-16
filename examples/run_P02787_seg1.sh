#!/bin/sh

mkdir -p ~/EMBL-EBI/funclan_work/P02787_1/ca_distances
mkdir -p ~/EMBL-EBI/funclan_work/P02787_1/distance_differences
mkdir -p ~/EMBL-EBI/funclan_work/P02787_1/distance_difference_maps
mkdir -p ~/EMBL-EBI/funclan_work/P02787_1/cluster_results
# Remove files from previous executions
rm ~/EMBL-EBI/funclan_work/P02787_1/ca_distances/*
rm ~/EMBL-EBI/funclan_work/P02787_1/distance_differences/*
rm ~/EMBL-EBI/funclan_work/P02787_1/distance_difference_maps/*
rm ~/EMBL-EBI/funclan_work/P02787_1/cluster_results/*

python3 run_find_conformers.py -u "P02787" \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/2hau_updated.cif A B \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/2hav_updated.cif A B \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3qyt_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3s9l_updated.cif A B C D \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3s9m_updated.cif A B C D \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3s9n_updated.cif A B C D \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3skp_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3v83_updated.cif A B C D E F \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3v89_updated.cif A B \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3v8x_updated.cif A B \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/3ve1_updated.cif A B C D \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/4h0w_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/4x1b_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/4x1d_updated.cif A B \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/5dyh_updated.cif A B \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/5h52_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/5wtd_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/5x5p_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/5y6k_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6ctc_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6d03_updated.cif A B C D E \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6d04_updated.cif A B C D E F \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6d05_updated.cif A B C D E F \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6jas_updated.cif A \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6soy_updated.cif A B C \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6soz_updated.cif A B C \
    -m /home/jellaway/EMBL-EBI/funclan_work/P02787_updated/6uj6_updated.cif A \
    -c ~/EMBL-EBI/funclan_work/P02787_1/ca_distances/ \
    -d ~/EMBL-EBI/funclan_work/P02787_1/distance_differences/ \
    -s ~/EMBL-EBI/funclan_work/P02787_1/cluster_results/ \
    -g ~/EMBL-EBI/funclan_work/P02787_1/cluster_results/ png svg

# -o ~/EMBL-EBI/funclan_work/P02787/distance_difference_maps/ \
