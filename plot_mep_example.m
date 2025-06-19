%gen_mep_curve_fig.m

proc_dir = '/Users/ERCOLE/Documents/Research/Repos/mEP_processing';
load('EP_040224.mat')

pt_ind = 8;
tr_nos = [1 5 7];

%%
pt_ind = 4;
settings_temp = ep_struct(pt_ind).stim_settings.param_strings;
%settings_temp = settings_temp(contains(settings_temp,'60'));
[settings_sorted, order] = sort(settings_temp);
for kk = 1:length(settings_temp)
    if ~contains(settings_sorted{kk},'2')
       continue 
    end
    fprintf('%d: %s\n',order(kk),settings_sorted{kk})
end

%% 
%ephys087:
%tr_nos = [26 7 15 24 27]; %2,3,4,5,6 mA

%ephys089
%tr_nos = [2,18,16,36]; %2,3,4,5 mA
%emg_ind = 2;

%ephys017
%tr_nos = [9,28,32]; %[3,21,24]; %3B  %[38,30,28]; %2B %1,3,5 mA
tr_nos = [13,22,23];
emg_ind = 2;

col = [1,0,0];
scale = linspace(0,1,length(tr_nos));
t_plot = linspace(0,100,2200);

figure; hold on
for kk = 1:length(tr_nos)
    temptrace = ep_struct(pt_ind).mep.raw.mep_means{tr_nos(kk),emg_ind};
    temptrace = temptrace- median(temptrace);
    plot(t_plot,temptrace,'Color',scale(kk)*col)    
end
ylim([-80,80])
xlim([0,50])
xlabel('Time (ms)'); ylabel('Voltage (\muV)')
%legend({'2A','2B', '2C'},'Location','Northeast')
legend({'1 mA','3 mA', '5 mA'},'Location','Northeast')
title(sprintf('%s: %s',ep_struct(pt_ind).patient_ID, ep_struct(pt_ind).mep.emg_labels{emg_ind}))
set(gca,'FontSize',16)

%% plot vs. segment
emg_ind = 4;
sp_array = [11,7,8,9,4,5,6,2];
%tr_nos = [27,30,8,34,21,37,16,10];  %3 mA
tr_nos = [2,28,9,32,24,25,18,31];

for kk = 1:length(sp_array)
    subplot(4,3,sp_array(kk))
    temptrace = ep_struct(pt_ind).mep.raw.mep_means{tr_nos(kk),emg_ind};
    temptrace = temptrace- median(temptrace);
    plot(t_plot,temptrace,'Color','black')   
    
    ylim([-40,40])
    xlim([0,50])
    xlabel('Time (ms)'); ylabel('Voltage (\muV)')
    title(settings_temp{tr_nos(kk)})
end


%% 
for kk = 1:length(ep_struct)
    fprintf('%s: %d\n',ep_struct(kk).patient_ID,length(unique(ep_struct(kk).stim_settings.amplitudes)))
end