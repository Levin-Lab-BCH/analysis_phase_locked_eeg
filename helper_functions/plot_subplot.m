function [] = plot_subplot(curr_title,data,times,P,grp,color_axis_limits,P1,indi_erp,freq_window,time_window)
title(curr_title)
fig = gca;
colormap(jet(256));
pos = get(gca,'position');
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)];

axis off;
h(1) =  axes('Position',[.1 .1 .9 .9].*s+q);
set(h(1), 'tag', 'ersp');
imagesc(data.(grp).ITCtimes,data.(grp).freqs,real(P));  set(gca,'ydir','normal'); colorbar;caxis(color_axis_limits)


pos = get(h(1), 'position');
dim =  pos.*[1.5 1.5 0.1 0.1];
%annotation(fX, 'textbox', dim ,'String', 'test 1');

%xline(P1,'-',{'st','end'},'Color','m','LineWidth',2)
ylim([min(freq_window) max(freq_window)]);
xlim([min(time_window(:)) max(time_window(:))])
if ~isempty(indi_erp)
    h(4) = axes('Position',[.1 .1-0.2 .8 .1].*s+q);
    plot(times,indi_erp); xlim([min(data.(grp).ITCtimes),max(data.(grp).ITCtimes)]); title({'ERP'})
end
end