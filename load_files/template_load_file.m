% Module Selection
grp_proc_info.compute_phase_locked_eeg = 0;
grp_proc_info.plot_phase_locked_eeg = 1;

% User inputs needed to feed into the function generate_PL_non_PL_ITPC.m

grp_proc_info.phase_lock_run_tag = '_ERP_happe_v3_raw_seg_blcorr_2_200_notch'; %beapp run tag to pull for the phsae locking analyses
grp_proc_info.eeg_lab_dir = 'C:\Users\ch220650\beapp\Packages\eeglab14_1_2b'; %path to your local eeglab package (usually in the beapp/Packages\eeglab14_1_2b folder)
grp_proc_info.help_fcn = 'Z:\Groups\SPA\02_Data_Processing\02_Post-processing\analysis_phase_locked_eeg\helper_functions'; %path to your helper functions YAEL to IDENTIFY WHICH YOU NEED
grp_proc_info.group_sub_dirs = {'TD','ASD','SPC'}; %group names (need to be the names of your subfolders of segment directory)

% Newtimef specifications
grp_proc_info.phase_lock_base_norm = 1; 
grp_proc_info.chans_egi = [37 42 36 41 42 47 28 30 18 38]; %EGI channels to compute phase locked EEG values for 
grp_proc_info.chans_ccf = []; % Brain vision channels to compute phase locked EEG values for (leave blank/goes unused for SPA)

%Bootstrapping specifications
grp_proc_info.min_trials = 90; %Code will only take participants that have at least X trials + 9 (ignores first 9 to avoid habituation effects) 
grp_proc_info.trial_subset = 50; %will take random sampling of X trials for N iterations specified below 
grp_proc_info.phase_locked_cond_names = {'F20','F40','F75','F125','F200'}; %conditions to compute phase locked power for, must be a subset of file_proc_info.grp_wide_possible_cond_names_at_segmentation
grp_proc_info.phase_locked_cond_values = [20 40 75 125 200]; %The numberical representation of the conditions, needs to be same number of elements as phase_locked_cond_names
grp_proc_info.iterations = 100;

%Split half reliability
grp_proc_info.run_shr = 0;
grp_proc_info.conditions_to_estimate = [20 40 75];% 125 200]; %Must contain all or a subset of the condition values run for your batch

%% Required if your scripts is not already loading these specifications from the user inputs in your segmented run
grp_proc_info.src_dir={'Z:\Groups\SPA\01_Data_Raw_Summary_Processed\EEG\Participant_Data\03_Processed_Data\01_In_Progress\09_Aurora'}; %the directory containing your source filesopen 
grp_proc_info.evt_seg_win_start = -.3;
grp_proc_info.evt_seg_win_end = .4;
grp_proc_info.evt_trial_baseline_win_start = -.1;
grp_proc_info.evt_trial_baseline_win_end = 0;

%% User Inputs needed to run characterize_PL_outputs.m
%if not loading previous ui will need
grp_proc_info.src_dir={'Z:\Groups\SPA\01_Data_Raw_Summary_Processed\EEG\Participant_Data\03_Processed_Data\01_In_Progress\09_Aurora'}; %the directory containing your source filesopen 
grp_proc_info.group_sub_dirs = {'TD','ASD','SPC'}; %group names (need to be the names of your subfolders of segment directory)

grp_proc_info.freq_band_lims = [12 31];
grp_proc_info.early_late_windows = [50 150; 150 300]; %2X2 matrix first row 

% Identify which batch to plot  % which batch do you want to plot, will be of the form beapp_run_tag '_', trialsubset 't_',iterations,'i'
grp_proc_info.pl_run_tag_toplot = '_ERP_happe_v3_raw_seg_blcorr_2_200_notch';
grp_proc_info.pl_trial_subset_toplot = 60;
grp_proc_info.pl_iterations_toplot = 100;
grp_proc_info.pl_chans_to_plot = [36];%[18 37 42 36]; %EGI channels to plot phase locked heat maps and input/output curves for
grp_proc_info.pl_grp_colors = {[0.00,0.45,0.74],[0.47,0.67,0.19],[0.49,0.18,0.56]}; %colors for each diagnostic group
grp_proc_info.pl_erp_win_to_highlight = [128 172]; %start and end in ms of time window to draw lines at over heat maps
grp_proc_info.conditions_to_estimate = [20 40 75 125]; %Must contain all or a subset of the condition values run for your batch
grp_proc_info.plot_grp_av = 0; %Toggle on/off plotting group averages
grp_proc_info.plot_indi_in_out = 1; %Toggle on/off plotting individual input output curves
grp_proc_info.plot_indi_lines_of_fit = 0; %Toggle on/off individual lines of best fit