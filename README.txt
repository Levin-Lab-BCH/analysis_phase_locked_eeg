READ ME

Author: Yael Braverman - Code for Generating Phase-Locked and Non-Phase-Locked EEG Responses


% You will need :
%   - BEAPP segment files seperated into subfolders based on Diagnostic
%   Group
%   - Helper functions folder that contains modified functions needed to
%   run this script
% - User inputs to specify parameters to run the script (see template_load_file.m in helper_functions subfolder)
%   - EEGLAB dev version released 08/10/2023
%       - dip_fit v5.3
%       - firfilt v2.7.1
%       - cleanrawdata 2.0

% Configure load file in load_files
% Configure paths in run_phase_locked_eeg_anlaysis
% Run run_phase_locked_eeg_analysis

%Details on code executed to genrate the phase locked data:
% Adapted from KVS INSAR 2024 Code 
% Changes from KVS code:
% - Uses architecture of the newtimef and diagnostic group sorting, and
% modified beapp2eeglab
% - Does not plot output in code, uses a struct for all data outputs
% tracked by subject id rather than cell arrays to prevent overwriting
% across multiple runs
% - Added bootstrapping iterations, Added computation of phase locked vs
% not phase locked, add metadata and added

% Built to interface with the SPA "run_beapp_spa_[TASK]" as part of the postbeapp
% step. Parent function will:
% 1. First, reference the user input script for the segmentation script to
% load base user inputs (e.g source directory, event tags, baseline and
% segment windows etc)
% 2. Reference the postbeapp user input script to load inputs specific to
% this phase locked EEG module (selects which conditions, channels to
% analyze outputs for, how many trials to pull for each group)
% Main function generates estimates of itpc and phase locked and non-phase
% locked power through a bootstrapping procedure; it will select a random
% subset of N trials (specified by user), generate estimates, and repeate
% iti times (specified by user), the code will then average over these
% iterations and save the averaged estimates.
% Main function will only compute new estimates by:
% 1. Checking and loading previous dataset with this beapp run tag and trial subset/ iteration combination
% 2. Checking that files have not previously been computed for the given
% beapp run tag for this channel/trial subset/iteration combination, if so
% will skip to next
% 3. Checking that files have enough trials for the trial subset in all
% conditions presented to them