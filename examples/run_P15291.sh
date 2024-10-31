#!/bin/sh

# Create directories if they don't exist
mkdir -p benchmark_data/examples/P15291/P15291_updated_mmcif
mkdir -p benchmark_data/examples/P15291/P15291_ca_distances/unp_residue_ids
mkdir -p benchmark_data/examples/P15291/P15291_distance_differences
mkdir -p benchmark_data/examples/P15291/P15291_cluster_results
mkdir -p benchmark_data/examples/P15291/P15291_alpha_fold_mmcifs

# Remove files if they exist from previous runs
rm benchmark_data/examples/P15291/P15291_ca_distances/*
rm benchmark_data/examples/P15291/P15291_ca_distances/unp_residue_ids/*
rm benchmark_data/examples/P15291/P15291_distance_differences/*
rm benchmark_data/examples/P15291/P15291_cluster_results/*

# mprof run --python
python3 find_conformers.py -u "P15291" \
    -m benchmark_data/examples/P15291/P15291_updated_mmcif/2fyb_updated.cif A \
    -m benchmark_data/examples/P15291/P15291_updated_mmcif/2fyd_updated.cif B D \
    -m benchmark_data/examples/P15291/P15291_updated_mmcif/6fwu_updated.cif A B \
    -m benchmark_data/examples/P15291/P15291_updated_mmcif/2fya_updated.cif A \
    -m benchmark_data/examples/P15291/P15291_updated_mmcif/2fy7_updated.cif A \
    -m benchmark_data/examples/P15291/P15291_updated_mmcif/2fyc_updated.cif B D \
    -c benchmark_data/examples/P15291/P15291_ca_distances/ \
    -d benchmark_data/examples/P15291/P15291_distance_differences/ \
    -s benchmark_data/examples/P15291/P15291_cluster_results/ \
    -o benchmark_data/examples/P15291/P15291_distance_difference_maps/ \
    -g benchmark_data/examples/P15291/P15291_cluster_results/ png svg \
    -f \
    -0 122 \
    -1 398 \
    # -w benchmark_data/examples/P15291/P15291_cluster_results/ png svg \
    # -a benchmark_data/examples/P15291/P15291_alpha_fold_mmcifs \
    # -v \
