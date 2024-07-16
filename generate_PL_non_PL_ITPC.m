function [] = generate_PL_non_PL_ITPC(grp_proc_info_in)
% potential future improvements: making iterating piece its own function, check about flowing into plotting code
% completed stuff: took out excess code, started mapping inputs to user input script, tested if code breaks, tested if the channel thing works in updating/saving - yes works!
%% Author: Yael Braverman - Code for Generating Phase-Locked and Non-Phase-Locked EEG Responses
% Adapted from KVS INSAR 2024 Code 
% Changes from KVS code:
% - Uses architecture of the newtimef and diagnostic group sorting, and
% modified beapp2eeglab
% - Does not plot output in code, uses a struct for all data outputs
% tracked by subject id rather than cell arrays to prevent overwriting
% across multiple runs
% - Added bootstrapping iterations, Added computation of phase locked vs
% not phase locked, add metadata and added

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

%% Pull In Inputs Inputs

% reliability inputs may not use later
run_shr = grp_proc_info_in.run_shr;
run_shr = 0; %setting to 0 because it is not fully supported



% directory inputs
iterations = grp_proc_info_in.iterations; %how many shuffled estimates to generate/collapse over
eeglab_dir = grp_proc_info_in.eeg_lab_dir;
help_fcn = grp_proc_info_in.help_fcn; % 
grp_dir = fullfile(grp_proc_info_in.src_dir{1,1},['segment',grp_proc_info_in.phase_lock_run_tag])  ;
gpi_dir = fullfile(grp_proc_info_in.src_dir{1,1},['out',grp_proc_info_in.phase_lock_run_tag],['Run_Report_Variables_and_Settings',grp_proc_info_in.phase_lock_run_tag,'.mat']);
group_sub_dir_name = grp_proc_info_in.group_sub_dirs;
% output dirs and names
out_dir = fullfile(grp_proc_info_in.src_dir{1,1},'Phase_Locked_Power');
save_name = strcat('PhaseLocking',grp_proc_info_in.phase_lock_run_tag,'_',num2str(grp_proc_info_in.trial_subset),'t_',num2str(iterations),'i.mat'); % saves the phase locking mat by beapp run tag, number of trials subsampled and number of bootstrapping iterations
shr_name = strcat('PhaseLocking',grp_proc_info_in.phase_lock_run_tag,'_shr.mat'); %file to store split half reliability metrics
% ITPC ERSP Parameters to pull in from inputs
baseline_start = grp_proc_info_in.evt_trial_baseline_win_start*1000; %Pulled from segmentation user input script
baseline_end = grp_proc_info_in.evt_trial_baseline_win_end*1000; % Pulled from segmentation user input script
base_norm = grp_proc_info_in.phase_lock_base_norm;
chans_egi = grp_proc_info_in.chans_egi;% EGI channels desired for analysis
chans_ccf = grp_proc_info_in.chans_ccf; % Brain Vision channesl desired for analysis

% Switches (0 is no/off, 1 is yes/on)
add_paths =1; % if you need to add your directory paths
min_trials = grp_proc_info_in.min_trials; 
trials_subset =grp_proc_info_in.trial_subset; %35;%75;%35;%70;%84; % how many trials should be taken from each participant
if run_shr; if mod(trials_subset,2); trials_subset = trials_subset - 1; end ;end; %If running split half reliability set trials to be even
cond_names = grp_proc_info_in.phase_locked_cond_names; %conditions to compute phase locked power for, must be a subset of file_proc_info.grp_wide_possible_cond_names_at_segmentation
%iterations = 10;
%% FLAGGED PARAMETERS
start_ms =160;%100; %25; %60 74 FLAG IF NEEDED - ends up not getting used b/c not actually plotting in this code
end_ms = 240;%150;%65; %130 138 FLAG IF NEEDED- ends up not getting used b/c not actually plotting in this code
%Plotting parameters
%% HARD CODED ITPC params currently hard coded (goes into newtimef)
% see newtimef from more details for min/max cycles
alp = 0.05; %alpha to determine significance NOT CURRENTLY USED
min_cyc =1 ;%1; %1    Typical value: 'cycles', [3 (min_cyc), 0.5 (max_cyc)]
max_cyc =.5; %.65 %maximun # of cycles to use. >1, fixed cycles. <1, percentage of cycles to use at max frequceny
min_freq = 10; %2; %3 %minfreq being determined by the number of data points, cycles and sampling frequency. {default [minfreq 50]}
max_freq = 80;%30; %100 % max freq desired
pad_ratio = 2;%16;

%% Add necessary folders to path/create necessary output folders
if ~isfolder(out_dir); mkdir(out_dir); end
s       = pathsep;
pathStr = [s, path, s];
onPath_eeglab  = contains(pathStr, [s, eeglab_dir, s]);
onPath_helper_fcn = contains(pathStr,[s,help_fcn,s]);
if ~onPath_eeglab; cd(eeglab_dir); eeglab nogui ;end
if ~onPath_helper_fcn; addpath(help_fcn); end

%% If batch file currently exists load it 
if isfile(fullfile(out_dir,save_name))
    load(fullfile(out_dir,save_name))
else
    batch_info.collapsed_metrics_computed = table();
    batch_info.cond_values = grp_proc_info_in.phase_locked_cond_values;
end

conditions_to_estimate = grp_proc_info_in.conditions_to_estimate;
[~,~,conditions_to_plot] = intersect(conditions_to_estimate , batch_info.cond_values);
conditions_to_plot = conditions_to_plot';

% If split half reliability metrics exists load it
if run_shr
    if isfile(fullfile(out_dir,shr_name))
        load(fullfile(out_dir,shr_name))
    else
        shr_struct.shr_metrics = table();
    end
end
%% Build out variables to be used in call to newtimef
P1 = [start_ms end_ms]; %sort window values is center_ms_amp = INF
freqs = [min_freq max_freq];
cycles = [min_cyc max_cyc];
toggle_plot = 'off';

for ng = 1:length(group_sub_dir_name) %group loop
    %initialize file count
    file_count = 0;
    %go into diagnostic group sub-folder
    grp = group_sub_dir_name{ng};
    dir_str = fullfile(grp_dir,grp);
    mydir = dir(strcat(dir_str,'\*.mat'));
    cd(dir_str);
    n_part = length(mydir);

    for file = 1:n_part %individual loop
        name = erase(mydir(file).name,'.mat');
        name = strrep(name,'_','-');
            %Set channels based on file type
            if contains(name,'7903-01') %find in using Brain vision ch
                channels = chans_ccf;
            else
                channels = chans_egi;
            end %end find brain vision ch

        %Clean up file name
        n =strrep(name,'-','_');
        n = strcat('s',n);
        % Check if metrics have already been computed for this file
        if any(strcmp(batch_info.collapsed_metrics_computed.Properties.RowNames,n)) %Check if this subject/file is a row in the metrics table
            %Check which channels you want to analyze are present in old
            %batch metrics table
            log_chans_overlap = cell2mat(cellfun(@(x) any(strcmp(batch_info.collapsed_metrics_computed.Properties.VariableNames,x)),strsplit(num2str(channels)),'UniformOutput',false));
            chans_present_in_batch_table = channels(log_chans_overlap);
            % Of channels present in metrics table, check which are completed
            completed_chans = chans_present_in_batch_table(batch_info.collapsed_metrics_computed{n,strsplit(num2str(chans_present_in_batch_table))} == 1);
            channels = setdiff(channels,completed_chans); %identifies which channels have yet to be completed for this participant
            if isempty(channels) %if all channels have been completed, skip to next participant
                continue
            end
        end      
        %% load file_proc_info to check for trial count   
        load(mydir(file).name,'file_proc_info'); 

        %find out how many trials there are per conditions you want to analyze

        %find the rows of your evt_conditions_being_analyzed that contain your conditions of interst
        cond_row_idx = (cellfun(@(x) find(strcmp(x,file_proc_info.evt_conditions_being_analyzed.Evt_Codes)),cond_names,'UniformOutput',false));
        conds_not_seen = cellfun(@(x) isempty(x),cond_row_idx);
        cond_row_idx = cell2mat(cond_row_idx(~conds_not_seen)); %ignore conditions participant didn't see at all (due to old version of paradigm etc)
        
        % Get the number of trials post rejection for conditions
        % seen/selected for analysis
        tris = file_proc_info.evt_conditions_being_analyzed{cond_row_idx,'Num_Segs_Post_Rej'};

        %% If participant has enough trials, load them in
        if ~any(tris < min_trials | tris < trials_subset +9) %check trials loop - if participant has sufficient trials to take subset for all relevant conditions and minimum trial threshold
            %load file
            curr_file = load(mydir(file).name); 
            curr_file.file_proc_info_orig = curr_file.file_proc_info;
            file_count = file_count + 1; 

            for i_chan = channels % iterate over channels that haven't been completed yet
                %% Repeat generating the phase locked, non-phase locked, and itpc response for iti bootstrapped iterations
                for iti = 1:iterations
                    for i_c = 1:size(cond_names,2)%condition loop
                        %skip condition if participant didn't have any trials for this condition
                        if ismember(i_c,find(conds_not_seen))
                            continue;
                        end
                        curr_condition = find(strcmp(cond_names{i_c},file_proc_info.grp_wide_possible_cond_names_at_segmentation)); %gets the actual condition that it's in
                        
                        
                        tr = sort(randperm(tris(i_c,1)-9,trials_subset)); %get a random sampling of trials for the amount of trials desired
                        tr = tr+9; %offset by 10 so you dont get the first 10 trials to avoid habituation effects
                       
                        %If run_shr
                        if run_shr
                            shuffled_order = randperm(length(tr));
                            curr_subset(1,:) = tr(shuffled_order(1:floor(length(tr)/2))); % first trial set for split half reliability
                            curr_subset(2,:) = tr(shuffled_order(floor(length(tr)/2)+1:end)); % second trial set for split half reliability
                        else
                            curr_subset(1,:) = tr; %take all positions of tr
                        end                        
                        % Add loop to generate split-half reliability if
                        % desired - will just loop once if split half
                        % reliability is not turned on, otherwise loops
                        % twice
                        for shr = 1:size(curr_subset,1)
                        
                        curr_file.temp_eeg_w{curr_condition} = curr_file.eeg_w{curr_condition}(:,:,curr_subset(shr,:));
                        data.(grp).tr.(n)(shr,:,iti) = curr_subset(shr,:);
                        if isfield(curr_file.file_proc_info,'epoch')
                            curr_file.file_proc_info.epoch{curr_condition} = curr_file.file_proc_info_orig.epoch{curr_condition}(:,curr_subset(shr,:)); %get the event info for the specific trails chosen at random
                        end

                        EEG_orig = beapp2eeglab_phase_locking(curr_file.file_proc_info,curr_file.temp_eeg_w{curr_condition},1,1,grp_proc_info_in,curr_condition); %turn BEAPP files into EEGLAB files
                                               
                        data.(grp).indi.(n){1,i_c,shr} = squeeze(mean(EEG_orig.data(i_chan,:,:),1,"omitnan"));
                        data.(grp).indi_avg.(n){1,i_c,shr} = mean(data.(grp).indi.(n){1,i_c,shr},2); %average over trials
                        %Prepare "Total Power" data
                        EEG = EEG_orig;
                        %update the data to be based on channels
                        EEG.data = data.(grp).indi.(n){1,i_c,shr}; %force the EEGLAB structure to have the averaged over channels data ( timepoints x trials)

                        %Prepare "Non-phase_locked Power" data - subtract indiv_avg                        
                        data.(grp).indi_non_phase_locked.(n){1,i_c,shr} = data.(grp).indi.(n){1,i_c,shr} - data.(grp).indi_avg.(n){1,i_c,shr};
                        
                        EEG_non_phase_locked = EEG;
                        EEG_non_phase_locked.data =  data.(grp).indi_non_phase_locked.(n){1,i_c,shr};

                        %Prepare newtimef variables
                        frames = length(EEG.times);
                        epochlim = [EEG.xmin EEG.xmax]*1000 ; %turn min and max times into seconds
                        srate = EEG.srate;

                        if ~base_norm
                            %compute variables for total power
                            [ data.(grp).ersp.(n){i_chan,i_c,iti,shr},data.(grp).itc.(n){i_chan,i_c,iti,shr},data.(grp).powbase.(n){1,i_c,shr},data.(grp).ITCtimes,data.(grp).freqs,data.(grp).erspboot.(n){1,i_c,shr},data.(grp).itcboot.(n){1,i_c,shr},data.(grp).alltfx.(n){1,i_c,shr}] = ...
                                newtimef(EEG.data, frames, epochlim, srate,cycles,'freqs',freqs,...
                                'wletmethod','dftfilt3','vert',P1,'plotphasesign','off','plotphaseonly','off','plotitc',toggle_plot,'plotersp','off','padratio',pad_ratio,'alpha',alp);%'alpha',.05,'verbose','on');%,'alpha',alp,'verbose','on');
                            %compute variables for non phase locked signal
                            [ data.(grp).ersp_nonPL.(n){i_chan,i_c,iti,shr},data.(grp).itc_nonPL.(n){i_chan,i_c,iti,shr},data.(grp).powbase_nonPL.(n){1,i_c,shr},data.(grp).ITCtimes_nonPL,data.(grp).freqs_nonPL,~,data.(grp).itcboot_nonPL.(n){1,i_c,shr},data.(grp).alltfx_nonPL.(n){1,i_c,shr}] = ...
                                newtimef(EEG_non_phase_locked.data, frames, epochlim, srate,cycles,'freqs',freqs,...
                                'wletmethod','dftfilt3','vert',P1,'plotphasesign','off','plotphaseonly','off','plotitc','off','plotersp','off','padratio',pad_ratio, 'alpha',alp);%,'alpha',alp,'verbose','on');
                        else
                            %compute variables for total power
                            [ data.(grp).ersp.(n){i_chan,i_c,iti,shr},data.(grp).itc.(n){i_chan,i_c,iti,shr},data.(grp).powbase.(n){1,i_c,shr},data.(grp).ITCtimes,data.(grp).freqs,data.(grp).erspboot.(n){1,i_c,shr},data.(grp).itcboot.(n){1,i_c,shr},data.(grp).alltfx.(n){1,i_c,shr}] = ...
                                newtimef(EEG.data, frames, epochlim, srate,cycles,'freqs',freqs,...
                                'wletmethod','dftfilt3','vert',P1,'plotphasesign','off','plotphaseonly','off','basenorm','on','itcmax',[0 .4],'baseline',[baseline_start baseline_end],'plotitc',toggle_plot,'plotersp','off','verbose','on','padratio',pad_ratio,'alpha',alp);%,'alpha',alp,'verbose','on');

                            %compute variables for non phase locked signal
                            [ data.(grp).ersp_nonPL.(n){i_chan,i_c,iti,shr},data.(grp).itc_nonPL.(n){i_chan,i_c,iti,shr},data.(grp).powbase_nonPL.(n){1,i_c,shr},data.(grp).ITCtimes_nonPL,data.(grp).freqs_nonPL,data.(grp).erspboot_nonPL.(n){1,i_c,shr},data.(grp).itcboot_nonPL.(n){1,i_c,shr},data.(grp).alltfx_nonPL.(n){1,i_c,shr}] = ...
                                newtimef(EEG_non_phase_locked.data, frames, epochlim, srate,cycles,'freqs',freqs,...
                                'wletmethod','dftfilt3','vert',P1,'plotphasesign','off','plotphaseonly','off','basenorm','on','baseline',[baseline_start baseline_end],'plotitc','off','plotersp','off','padratio',pad_ratio,'alpha',alp,'verbose','on');
                        end
                        %Compute Phase Locked Response
                        data.(grp).phase_locked.(n){i_chan,i_c,iti,shr} = data.(grp).ersp.(n){i_chan,i_c,iti,shr}-data.(grp).ersp_nonPL.(n){i_chan,i_c,iti,shr};
                        
                      

                        end % end split half reliability loop



                        if run_shr
                            orig_fields = {'phase_locked','ersp_nonPL','itc'};
                            for i_field = 1:length(orig_fields)
                            %correlate measures with one another for this
                            %iteration for whole matrix
                            data.(grp).(['reliability_all' orig_fields{i_field}]).(n){i_chan,i_c,iti} = corr(real(data.(grp).(orig_fields{i_field}).(n){i_chan,i_c,iti,1}(:)),real(data.(grp).(orig_fields{i_field}).(n){i_chan,i_c,iti,2}(:)));

                            %get SHR for just P1 window
                            a = data.(grp).(orig_fields{i_field}).(n){i_chan,i_c,iti,1}((data.(grp).freqs >= 13 & data.(grp).freqs <= 31),data.(grp).ITCtimes >= 50 & data.(grp).ITCtimes <= 150);
                            b = data.(grp).(orig_fields{i_field}).(n){i_chan,i_c,iti,2}((data.(grp).freqs >= 13 & data.(grp).freqs <= 31),data.(grp).ITCtimes >= 50 & data.(grp).ITCtimes <= 150);
                            data.(grp).(['reliability_beta_50_150' orig_fields{i_field}]).(n){i_chan,i_c,iti} = corr(real(a(:)),real(b(:)));   

                            
                            end %end fields loop
                        end
                    end %end condition loop

                    if run_shr
                        for shr = 1:2
                      % Get slope for 
                        if size(data.(grp).ersp_nonPL.(n),2)<max(conditions_to_plot)
                            continue
                        else
                            tic
                            data.(grp).(['phase_locked' strrep(num2str(conditions_to_estimate),' ','_')]).(n){i_chan,iti,shr} = get_tf_slopes(conditions_to_estimate, data.(grp).phase_locked.(n)(i_chan,conditions_to_plot,iti,shr));
                            toc
                            B = get_max_slope_table(data.(grp).(['phase_locked' strrep(num2str(conditions_to_estimate),' ','_')]).(n){i_chan,iti,shr},data.(grp).ITCtimes,data.(grp).freqs,n,i_chan);
                            shr_struct.(grp).shr_phase_locked_max_slope.(n){iti,shr,i_chan} = B.("Max Slope_Early");
                            shr_struct.(grp).shr_phase_locked_Freq_Early.(n){iti,shr,i_chan} = B.("Freq (Hz) Early");
                            shr_struct.(grp).shr_phase_locked_Time_Early.(n){iti,shr,i_chan} = B.("Time (ms) Early");
                        end
                        end

                        shr_struct.(grp).(['reliability_slope_all' 'phase_locked']).(n){i_chan,iti} = corr(real(data.(grp).(['phase_locked' strrep(num2str(conditions_to_estimate),' ','_')]).(n){i_chan,iti,1}(:)),real(data.(grp).(['phase_locked' strrep(num2str(conditions_to_estimate),' ','_')]).(n){i_chan,iti,2}(:)));

                    end
                end %end iteration loop

                % Compute collapsed / averaged phase locked, itpc, non
                % phase locked responses from iterations and toss the
                % temporary variables to save memory
                orig_fields = {'phase_locked','ersp_nonPL','itc','itc_nonPL'}; %corresponding fields to pull dat afrom at indexes in the max slope tables
                for i_field = 1:length(orig_fields)
                    y_axis_start = squeeze(data.(grp).(orig_fields{i_field}).(n)(i_chan,:,:,1)); % defaulting to taking first split of data
                    ii = 0;
                    for i_c = 1:size(cond_names,2)
                        if ismember(i_c,find(conds_not_seen))
                            continue;
                        end


                        y_axis_temp = reshape(y_axis_start(i_c,:),1,1,iterations);
                        data.(grp).([orig_fields{i_field},'_collapsed']).(n){i_chan,i_c} =  mean(cell2mat(y_axis_temp),3);
                        data.(grp).([orig_fields{i_field},'_collapsed_std']).(n){i_chan,i_c} = std( cell2mat(y_axis_temp),0,3);
                        data.(grp).([orig_fields{i_field},'_collapsed_std_error']).(n){i_chan,i_c} =  data.(grp).([orig_fields{i_field},'_collapsed_std']).(n){i_chan,i_c}/sqrt(iterations);
                        % fill in shr metrics
                        if run_shr
                            shr_struct.shr_metrics{n,['beta_50_150_cond_',cond_names{i_c}, orig_fields{i_field},'_chan',num2str(i_chan),'_',num2str(trials_subset/2),'_Trials_SHR']} = mean([data.(grp).(['reliability_beta_50_150' orig_fields{i_field}]).(n){i_chan,i_c,:}]);
                            % after test fill in overall
                            shr_struct.shr_metrics{n,['all_cond_',cond_names{i_c}, orig_fields{i_field},'_chan',num2str(i_chan),'_',num2str(trials_subset/2),'_Trials_SHR']} = mean([data.(grp).(['reliability_all' orig_fields{i_field}]).(n){i_chan,i_c,:}]);
                        end
                    end
                    % Once collapsed is set for this field, clear the
                    % iteration variables for memory/speed
                    data.(grp).(orig_fields{i_field}) = rmfield(data.(grp).(orig_fields{i_field}),(n));
                end
                % Iteratively save by file / channell
                batch_info.collapsed_metrics_computed{n,num2str(i_chan)} = 1; % mark in batch info that this channel has been computed
                save(fullfile(out_dir,save_name),'data','batch_info','-v7.3') %save file and batch table at this point
                if run_shr
                save(fullfile(out_dir,shr_name),'shr_struct');
                end
            end %end channel loop       
        else
            disp(strcat("Participant ",convertCharsToStrings(name)," does not have enough trials ( has ",num2str(tris'),", needs ",num2str(min_trials)," ) to be included in group avg."));
            continue;
        end %end check to see if participants have enough trials
    end %end individual loop
    if ~exist('data','var') || ~isfield(data,grp)
        continue
    end
end % end group  loop
