#!/bin/sh

path_mmcifs="benchmark_data/examples/O34926/O34926_updated_mmcif"
path_ca_distances="benchmark_data/examples/O34926/O34926_ca_distances"
path_distance_differences="benchmark_data/examples/O34926/O34926_distance_differences"
path_cluster_results="benchmark_data/examples/O34926/O34926_cluster_results"

rm $path_ca_distances/*
rm $path_ca_distances/unp_residue_ids/*
rm $path_distance_differences/*
rm $path_cluster_results/*

python3 find_conformers.py \
    -u "O34926" \
    -m $path_mmcifs/3nc7_updated.cif A B \
    -m $path_mmcifs/3nc5_updated.cif A B \
    -m $path_mmcifs/3nc3_updated.cif A B \
    -m $path_mmcifs/3nc6_updated.cif A B \
    -c $path_ca_distances \
    -d $path_distance_differences \
    -s $path_cluster_results \
    -g $path_cluster_results png svg \
    -0 1 \
    -1 405 \
    -f \
    # -v \
    # -a benchmark_data/examples/O34926/O34926_alpha_fold_mmcifs \
    # -o benchmark_data/examples/O34926/O34926_distance_difference_maps/ \
    # -w benchmark_data/examples/O34926/O34926_cluster_results/ png svg \
