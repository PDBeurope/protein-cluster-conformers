#!/bin/sh

path_base="/media/jellaway/FlashData/EMBL-EBI/funclan_work/clustering_w_uniref"
path_cluster_results="$path_base/cluster_results_uniref50"


path_mmcifs="$path_base/reindexed_mmcifs_combined/uniref50_mmcifs"
path_ca_distances="$path_cluster_results/ca_distances"
path_distance_differences="$path_cluster_results/distance_differences"
path_cluster_results="$path_cluster_results/cluster_results"

rm $path_ca_distances/*
rm $path_ca_distances/unp_residue_ids/*
rm $path_distance_differences/*
rm $path_cluster_results/*

python3 find_conformers.py \
    -u "uniref50" \
    -m $path_mmcifs/6mkh_updated.cif A \
    -m $path_mmcifs/6bsr_updated.cif A \
    -m $path_mmcifs/8f3v_updated.cif A \
    -m $path_mmcifs/6g88_updated.cif A B C \
    -m $path_mmcifs/6mki_updated.cif A \
    -m $path_mmcifs/5dvy_updated.cif A \
    -m $path_mmcifs/8f3t_updated.cif A \
    -m $path_mmcifs/8f3x_updated.cif A \
    -m $path_mmcifs/8f3f_updated.cif A \
    -m $path_mmcifs/8u55_updated.cif A \
    -m $path_mmcifs/8f3l_updated.cif A \
    -m $path_mmcifs/8f67_updated.cif A B C \
    -m $path_mmcifs/8f3s_updated.cif A \
    -m $path_mmcifs/6g0k_updated.cif A B C \
    -m $path_mmcifs/8f3n_updated.cif A \
    -m $path_mmcifs/8f3i_updated.cif A \
    -m $path_mmcifs/8f3o_updated.cif A \
    -m $path_mmcifs/8f3y_updated.cif A \
    -m $path_mmcifs/8f3p_updated.cif A \
    -m $path_mmcifs/8f3q_updated.cif A \
    -m $path_mmcifs/8f3m_updated.cif A \
    -m $path_mmcifs/6mkf_updated.cif A \
    -m $path_mmcifs/8f3z_updated.cif A \
    -m $path_mmcifs/6bsq_updated.cif A \
    -m $path_mmcifs/6mkj_updated.cif A \
    -m $path_mmcifs/6mkg_updated.cif A \
    -m $path_mmcifs/8f3g_updated.cif A \
    -m $path_mmcifs/8f3h_updated.cif A \
    -m $path_mmcifs/8f3r_updated.cif A \
    -m $path_mmcifs/8f3w_updated.cif A \
    -m $path_mmcifs/6mka_updated.cif A \
    -m $path_mmcifs/8f3j_updated.cif A \
    -m $path_mmcifs/8f3u_updated.cif A \
    -m $path_mmcifs/5e31_updated.cif A \
    -c $path_ca_distances \
    -d $path_distance_differences \
    -s $path_cluster_results \
    -g $path_cluster_results png svg \
    -f \
    -n 8 \
    # -v \
    # -a benchmark_data/examples/O34926/O34926_alpha_fold_mmcifs \
    # -o benchmark_data/examples/O34926/O34926_distance_difference_maps/ \
    # -w benchmark_data/examples/O34926/O34926_cluster_results/ png svg \
