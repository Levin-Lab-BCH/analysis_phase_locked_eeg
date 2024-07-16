%% beapp2eeglab (file_proc_info,eeg_curr_rec_period,curr_rec_period,event_data)

% convert file (one recording period) from beapp format to EEGLAB format
% only applied to data pre-segmentation
% Inputs: 
% file_proc_info: BEAPP file_proc_info structure
%
% eeg_curr_rec_period : eeg data for desired rec_period (in eeg) ex eeg{1} 
% 
% curr_rec_period : rec_period to pull file information from (bad chans,
% etc). should be integer ex. 1 

% event_data : 1 if you'd like to add events, 0 if not
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% The Batch Electroencephalography Automated Processing Platform (BEAPP)
% Copyright (C) 2015, 2016, 2017
% Authors: AR Levin, AS Méndez Leal, LJ Gabard-Durnam, HM O'Leary
%
% This software is being distributed with the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See GNU General
% Public License for more details.
%
% In no event shall Boston Children’s Hospital (BCH), the BCH Department of
% Neurology, the Laboratories of Cognitive Neuroscience (LCN), or software
% contributors to BEAPP be liable to any party for direct, indirect,
% special, incidental, or consequential damages, including lost profits,
% arising out of the use of this software and its documentation, even if
% Boston Children’s Hospital,the Laboratories of Cognitive Neuroscience,
% and software contributors have been advised of the possibility of such
% damage. Software and documentation is provided “as is.” Boston Children’s
% Hospital, the Laboratories of Cognitive Neuroscience, and software
% contributors are under no obligation to provide maintenance, support,
% updates, enhancements, or modifications.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License (version 3) as
% published by the Free Software Foundation.
%
% You should receive a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function EEG = beapp2eeglab_phase_locking(file_proc_info,eeg_curr_rec_period,curr_rec_period,event_data,gpi,curr_condition)
% initialize EEGLab structure using EEGLab script
EEG=eeg_emptyset;
EEG.data = eeg_curr_rec_period;
EEG.srate=file_proc_info.beapp_srate;
EEG.setname = file_proc_info.beapp_fname{1};
EEG.nbchan=file_proc_info.beapp_nchans_used(curr_rec_period);
EEG.chanlocs=file_proc_info.net_vstruct;
EEG.history=sprintf(['EEG = pop_importdata(''dataformat'',''matlab'',''nbchan'',0,''data'',',file_proc_info.beapp_fname{1},'''',',''setname'',''BEAPP_Dataset'',''srate'',0,''pnts'',0,''xmin'',0);\n','EEG = eeg_checkset( EEG );']);
EEG.noiseDetection.status = [];
EEG.noiseDetection.errors.status = '';
if isfield(file_proc_info,'evt_info') && event_data
    EEG =add_events_eeglab_struct(EEG,file_proc_info.evt_info{curr_rec_period});
end
%% KVS addition 9/1/23
for cevt = 1:length(file_proc_info.evt_info{curr_rec_period})
    EEG.event(cevt).epoch = cevt;
end
    if isfield(file_proc_info,'epoch')
        EEG.epoch = file_proc_info.epoch{curr_condition};
        %EEG.xmin = file_proc_info.epoch_xmin;

        EEG.xmin = gpi.evt_seg_win_start;
        EEG.xmax = gpi.evt_seg_win_end;
        %EEG.xmax = file_proc_info.epoch_xmax;
        EEG.times = [EEG.xmin*EEG.srate:EEG.srate/1000:EEG.xmax*EEG.srate-1];%file_proc_info.epoch_times;
        EEG.pnts = length(EEG.times);% file_proc_info.epoch_pnts;
        EEG.trials = size(EEG.data,3);
        EEG.urevent = EEG.event;
        %EEG.event = EEG.event([EEG.epoch.event].');
    end
%     EEG.icaact =
%     EEG.icawinv =
   % EEG.icachansind = file_proc_info.post_mara_icachansind{curr_rec_period};
    %EEG.icasphere = file_proc_info.post_mara_icasphere{curr_rec_period};
    %EEG.icaweights = file_proc_info.post_mara_icaweights{curr_rec_period};
     
%%
EEG=orderfields(EEG);
[~,EEG]=evalc('eeg_checkset(EEG)');
end
