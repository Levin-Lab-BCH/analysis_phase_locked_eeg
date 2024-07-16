%Input Slopes is a Nfreq by Ntimes matrix where each index represent the
%slope of a linear fit of the Nfreq by Ntimes original data across X
%conditions

%returns a table giving the maximum slope and which frequency/time position
%it occurred in post stimulus

function [output] = get_max_slope_table(input_slopes,times,freqs,n,i_chan,freq_window,time_window)

%freq_window : 1 X 2 matrix indicating the min and max frequency to look
%for
%time_window : 2x2 matrix indication min and max time windows to look at for early and late windows number splitting the early and late windows

i = 1;
temp = input_slopes;
%pad pre-stimulus window with NaNs
temp(:,times<time_window(1,1)) = NaN; 
temp((freqs<freq_window(1) | freqs >=freq_window(2)),:) = NaN; %only look at beta range

% do early window
temp_early = temp;
temp_early(:,times>=time_window(1,2)) = NaN;
[max_across_time, freq_ind] = max(temp_early);
[max_across_both, time_ind] = max(max_across_time);
max_slope_freq_ind = freq_ind(time_ind);

%do later window
temp_late = temp;
temp_late(:,times<time_window(2,1)) = NaN;
temp_late(:,times>=time_window(2,2)) = NaN;
[max_across_time_late, freq_ind_late] = max(temp_late);
[max_across_both_late, time_ind_late] = max(max_across_time_late);
max_slope_freq_ind_late = freq_ind_late(time_ind);

output = table(i_chan,max_across_both,max_slope_freq_ind,freqs(max_slope_freq_ind), time_ind, times(time_ind),...
                max_across_both_late,max_slope_freq_ind_late,freqs(max_slope_freq_ind_late),time_ind_late,times(time_ind_late),...
    'VariableNames',{'Channel','Max Slope_Early', 'Freq_Ind_Early',        'Freq (Hz) Early',         'Time_Ind Early','Time (ms) Early', ...
    'Max Slope_Late', 'Freq_Ind_Late',        'Freq (Hz) Late',         'Time_Ind Late','Time (ms) Late'},'RowNames',{[n(1:5),'_E',num2str(i_chan)]});
end