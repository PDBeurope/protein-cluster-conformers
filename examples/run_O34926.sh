#!/bin/sh

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
    -f \
    -g benchmark_data/examples/O34926/O34926_cluster_results/ png svg \
    # -0 1 \
    # -1 100
    # -v \
    # -a benchmark_data/examples/O34926/O34926_alpha_fold_mmcifs

# -o benchmark_data/examples/O34926/O34926_distance_difference_maps/ \
# -w benchmark_data/examples/O34926/O34926_cluster_results/ png svg \
