#!/bin/sh

# Create directories if they don't exist
mkdir -p benchmark_data/examples/O34926/O34926_updated_mmcif
mkdir -p benchmark_data/examples/O34926/O34926_ca_distances/unp_residue_ids
mkdir -p benchmark_data/examples/O34926/O34926_distance_differences
mkdir -p benchmark_data/examples/O34926/O34926_cluster_results
mkdir -p benchmark_data/examples/O34926/O34926_alpha_fold_mmcifs

# Remove files if they exist from previous runs
rm benchmark_data/examples/O34926/O34926_ca_distances/*
rm benchmark_data/examples/O34926/O34926_ca_distances/unp_residue_ids/*
rm benchmark_data/examples/O34926/O34926_distance_differences/*
rm benchmark_data/examples/O34926/O34926_cluster_results/*

# mprof run --python
python3 find_conformers.py -u "O34926" \
    -m benchmark_data/examples/O34926/O34926_updated_mmcif/3nc7_updated.cif A B \
    -m benchmark_data/examples/O34926/O34926_updated_mmcif/3nc5_updated.cif A B \
    -m benchmark_data/examples/O34926/O34926_updated_mmcif/3nc3_updated.cif A B \
    -c benchmark_data/examples/O34926/O34926_ca_distances/ \
    -s benchmark_data/examples/O34926/O34926_cluster_results/ \
    -d benchmark_data/examples/O34926/O34926_distance_differences/ \
    -m benchmark_data/examples/O34926/O34926_updated_mmcif/3nc6_updated.cif A B \
    -i 3nc6 \
    -o benchmark_data/examples/O34926/O34926_distance_difference_maps/ \
    -f \
    -g benchmark_data/examples/O34926/O34926_cluster_results/ png svg \
    -0 1 \
    -1 405 \
    # -w benchmark_data/examples/O34926/O34926_cluster_results/ png svg \
    # -a benchmark_data/examples/O34926/O34926_alpha_fold_mmcifs \
    # -v \
