## Overview over functions and Utils I have collected so far
# Layout
- BC-128-pass-lay.bvef
- BC-128-pass-lay.gif
- BC-128-pass-lay.mat
- Oder_electrodes.txt (changed the order for the connectivity matrices)
# Preprocessing
- Custom made function for separating the eyes open from the eyes closed: mytrialfun.m

- Code from Mike X Cohen for the Laplacian: laplacian_perrinX.m
https://github.com/mikexcohen/ANTS_youtube_videos/tree/main/ANTS6_preprocessing

# Power
FOOOF wrapper for MATLAB (needed if you want to fit a knee)
- fooof_check_settings
- fooof_get_model
- fooof_group
- fooof_plot
- fooof_unpack_results
- fooof_version
- fooof

# Connectivity
DISCOVER-EEG Automatic preprocessing pipeline for data in BIDS Format (Gil Ávila et al., 2023)
- Small worldness (calculation: compute_graph_measures.m was used and modified)
- isconnected.m: have a look if a matrix is connected
- correlation_illustration_SummerSchool: courtesy of Janus Rønn Lind Kobbersmed


# Scripts to save space
- generate_big_table
- generate_thresholds
