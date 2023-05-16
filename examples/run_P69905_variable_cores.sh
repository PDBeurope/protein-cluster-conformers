#!/bin/sh

# rm benchmark_data/examples/O34926/O34926_ca_distances/*
# rm benchmark_data/examples/O34926/O34926_distance_differences/*
# rm benchmark_data/examples/O34926/O34926_cluster_results/*




# mprof run --
for nproc in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
do
    for repeat in 1 2 3
    do
        echo "repeat = $repeat, nproc = $nproc, nchains = 310, remaking matxs each time"
        # rm /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/ca_distances/*
        # rm /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/ca_distances/unp_residue_ids/*
        # rm -r /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/dd_matrices/
        # mkdir /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/dd_matrices/
        rm /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/cluster_results/*

        { time python /home/jellaway/EMBL-EBI/funclan_work/contact_map_difference/find_conformers.py -u "P69905" \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1a00_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1a01_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1a0u_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1a0z_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1a3n_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1a3o_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1a9w_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1abw_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1aby_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1aj9_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1b86_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1bab_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1bbb_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1bij_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1buw_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1bz0_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1bz1_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1bzz_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1c7b_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1c7c_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1c7d_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1cls_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1cmy_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1coh_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1dke_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1dxt_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1dxu_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1dxv_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1fdh_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1fn3_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1g9v_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1gbu_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1gbv_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1gli_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1gzx_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hab_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hac_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hba_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hbb_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hbs_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hco_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hdb_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hga_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hgb_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hgc_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1hho_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ird_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1j3y_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1j3z_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1j40_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1j41_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1j7s_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1j7w_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1j7y_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1jy7_updated.cif.gz A C E G I K \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1k0y_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1k1k_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1kd2_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1lfl_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1lfq_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1lft_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1lfv_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1lfy_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1lfz_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ljw_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1m9p_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1mko_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1nej_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1nih_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1nqp_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1i_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1j_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1k_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1l_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1m_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1n_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1o_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1o1p_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1qi8_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1qsh_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1qsi_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1qxd_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1qxe_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1r1x_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1r1y_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1rps_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1rq3_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1rq4_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1rqa_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1rvw_updated.cif.gz A \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1sdk_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1sdl_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1shr_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1si4_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1thb_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1uiw_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1vwt_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xxt_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xy0_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xye_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xz2_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xz4_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xz5_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xz7_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xzu_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1xzv_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y01_updated.cif.gz B \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y09_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y0a_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y0c_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y0d_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y0t_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y0w_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y22_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y2z_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y31_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y35_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y45_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y46_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y4b_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y4f_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y4g_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y4p_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y4q_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y4r_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y4v_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y5f_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y5j_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y5k_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y7c_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y7d_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y7g_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y7z_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y83_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y85_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1y8w_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ydz_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ye0_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ye1_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ye2_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yen_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yeo_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yeq_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yeu_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yev_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yff_updated.cif.gz A C E G \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yg5_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ygd_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1ygf_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yh9_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yhe_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yhr_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yie_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yih_updated.cif.gz A C \
            -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yvq_updated.cif.gz A C \
            -c /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/ca_distances/ \
            -d /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/dd_matrices/ \
            -s /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/cluster_results/ \
            -n $nproc ; } 2>> variable_nproc_results_reusing_matxs_repeated_230214_half_chains.log
    done
done

            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yvt_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1yzi_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/1z8u_updated.cif.gz B D \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2d5z_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2d60_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2dn1_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2dn2_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2dn3_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2dxm_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2h35_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hbc_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hbd_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hbe_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hbf_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hbs_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hco_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hhb_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hhd_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2hhe_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2m6z_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2w6v_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2w72_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/2yrs_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3b75_updated.cif.gz A C E G I \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3d17_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3d7o_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3dut_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3hhb_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3hxn_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3ia3_updated.cif.gz B D \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3ic0_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3ic2_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3kmf_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3nl7_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3nmm_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3odq_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3onz_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3oo4_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3oo5_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3ovu_updated.cif.gz C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3p5q_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3qjb_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3qjc_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3qjd_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3qje_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3r5i_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3s48_updated.cif.gz B D \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3s65_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3s66_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3szk_updated.cif.gz A E \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3wcp_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/3whm_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4fc3_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4hhb_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4ij2_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4l7y_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4m4a_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4m4b_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4mqc_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4mqg_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4mqh_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4mqi_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4mqj_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4mqk_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4n7n_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4n7o_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4n7p_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4n8t_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4ni0_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4ni1_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4rol_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4rom_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4wjg_updated.cif.gz A F K P U Z \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4x0l_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/4xs0_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5e29_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5e6e_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5e83_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5ee4_updated.cif.gz C E \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5hu6_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5hy8_updated.cif.gz A C E G I \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5jdo_updated.cif.gz C E \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5kdq_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5ksi_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5ksj_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5ni1_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5sw7_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5u3i_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5ucu_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5ufj_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5urc_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5vmm_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5wog_updated.cif.gz A B \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5woh_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5x2r_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5x2s_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5x2t_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/5x2u_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6bb5_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6bnr_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6bwp_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6bwu_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6di4_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6hal_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6hbw_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6hk2_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6ka9_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kae_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kah_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kai_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kao_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kap_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kaq_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kar_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kas_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kat_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kau_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kav_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6kye_updated.cif.gz A C E G I K \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6l5v_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6l5w_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6l5x_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6l5y_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6lcw_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6lcx_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6nbc_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6nbd_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6nq5_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6tb2_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6xd9_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6xdt_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/6xe7_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7aet_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7aeu_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7aev_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7cue_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7dy3_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7dy4_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7jjq_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7jxz_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7jy0_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7jy1_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7jy3_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7k4m_updated.cif.gz A C E G I \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7pcf_updated.cif.gz A \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7pch_updated.cif.gz A D \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7pcq_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7qu4_updated.cif.gz C D \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7ud7_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7ud8_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7uf6_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7uf7_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7vde_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/7xgy_updated.cif.gz A C \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/8dov_updated.cif.gz A C E G \
            # -m /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/updated_mmcifs/8egi_updated.cif.gz A C \

# -g /home/jellaway/EMBL-EBI/funclan_work/many-chain-unps/P69905/cluster_results/ png svg \
# -a benchmark_data/examples/O34926/O34926_alpha_fold_mmcifs
# -o benchmark_data/examples/O34926/O34926_distance_difference_maps/ \
# -w benchmark_data/examples/O34926/O34926_cluster_results/ png svg \
