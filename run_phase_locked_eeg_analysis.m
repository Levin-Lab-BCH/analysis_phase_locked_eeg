%% Load file 

%edit path to your load file
load_file_path = 'Z:\Groups\SPA\02_Data_Processing\02_Post-processing\analysis_phase_locked_eeg\load_files\template_load_file.m';
code_dir =    'Z:\Groups\SPA\02_Data_Processing\02_Post-processing\analysis_phase_locked_eeg'; % path to where analysis_phase_locked_eeg is
%% DO NOT EDIT

% Load in inputs
run(load_file_path)
addpath(fullfile(code_dir,'helper_functions'))
addpath(code_dir)
if grp_proc_info.compute_phase_locked_eeg
    % Run Phase Locked EEG analysis

    generate_PL_non_PL_ITPC(grp_proc_info)
end
if grp_proc_info.plot_phase_locked_eeg
    % Characterize and Plot Outputs

    characterize_PL_outputs(grp_proc_info)
end


