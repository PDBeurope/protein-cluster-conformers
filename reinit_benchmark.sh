#!/bin/sh

# Remvoe benchmark results files
rm benchmark_data/all_uniprots/ca_distances/*
rm benchmark_data/all_uniprots/clustering_results/dendrograms/*
rm benchmark_data/all_uniprots/distance_difference_maps/*
rm benchmark_data/all_uniprots/distance_differences/all_conformers/*
rm -r benchmark_data/alphafold_mmcifs/*

# Reinitialising select examples from benchmark results
for SUBDIR in "O34926_ca_distances" "O34926_distance_differences" "O34926_cluster_results" "O34926_distance_difference_maps"
do
    rm ./benchmark_data/examples/O34926/$SUBDIR/*
done
rm -r benchmark_data/examples/O34926/O34926_alpha_fold_mmcifs/*

# Reinitialising select examples from benchmark results
for SUBDIR in "P15291_ca_distances" "P15291_distance_differences" "P15291_cluster_results" "P15291_distance_difference_maps"
do
    rm ./benchmark_data/examples/P15291/$SUBDIR/*
done
rm -r benchmark_data/examples/P15291/P15291_alpha_fold_mmcifs/*
