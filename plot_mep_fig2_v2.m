%%% plot_mep_fig2.m
%
%generates plot for diagram showing mEP detection overview

%% 
load('EP_112624_final_all.mat');


%%

pt = 9;
trial = 18;

disp(ep_struct(pt).stim_settings.param_strings{trial})

y_offset = 300;

mep_traces = ep_struct(pt).mep.raw.mep_means(trial,:);
emg_labels = ep_struct(pt).mep.emg_labels;

mep_traces{4} = ep_struct(18).mep.mep_means{1,1};

tt = linspace(0,100,2200);

figure
for kk = 1:length(mep_traces)
    subplot(4,2,kk)
    temp = mep_traces{kk} -median(mep_traces{kk}) + 20*(rand - 0.5);
    %temp = mep_traces{kk};
    %plot(tt, mep_traces{kk} + y_offset*(kk-1), 'k');
    plot(tt, temp, 'k');
    hold on
    title(emg_labels{kk})
    
    %xlim([0,100])
    ylim([-25, 25])
    if kk == length(mep_traces)
        ylim([-75, 75])
    end
    if kk == 7
        xlabel('Time (ms)')
    end
end

%ylim([-100, y_offset*9])

%%

pt = 9;
trial = 18;
chans = [1,2,5,8];

disp(ep_struct(pt).stim_settings.param_strings{trial})

y_offset = 300;

mep_traces = ep_struct(pt).mep.raw.mep_means(trial,chans);
emg_labels = ep_struct(pt).mep.emg_labels;

tt = linspace(0,100,2200);

figure
for kk = 1:length(mep_traces)
    subplot(4,1,kk)
    temp = mep_traces{kk} -median(mep_traces{kk}) + 20*(rand - 0.5);
    %temp = mep_traces{kk};
    %plot(tt, mep_traces{kk} + y_offset*(kk-1), 'k');
    plot(tt, temp, 'k');
    hold on
    title(emg_labels{chans(kk)})
    
    %xlim([0,100])
    ylim([-25, 25])
    if kk == length(mep_traces)
        ylim([-50, 50])
    end
    if kk == 4
        xlabel('Time (ms)')
    end
end

%%

% pt = 9;
% trial = 19;
% chans = [1,2,5,8];

save_figs = true;
overwrite_examples = false;

% pt_list = [9,9,9,11];
% trial_list = [1,19,19,3];
% ch_list = [1,2,5,7];
pt = 7;
ch = 4;
trial = 32;

pt_list = [9,9,9,11];
trial_list = [1,19,19,3];
ch_list = [1,2,5,7];

ltcs = {[14,26],[14,32],[5 16], [5,16]};
thrs = [7,7,4,4];
lims = [30,75,60,75];

cols = {'k','b','b','r'};

disp(ep_struct(pt).stim_settings.param_strings{trial})

y_offset = 300;

if overwrite_examples
mep_traces = cell(4,1);
emg_labels = cell(4,1);
for kk = 1:length(pt_list)
    pt = pt_list(kk); trial = trial_list(kk); 
    mep_traces{kk} = ep_struct(pt).mep.raw.mep_means{trial,ch_list(kk)};
    emg_labels{kk} = ep_struct(pt).mep.emg_labels{ch_list(kk)};
end
mep_traces_corr = mep_traces;
end

emg_labels =  {'L Bicep', 'L ECR', 'L Nasalis', 'R Orb. Oris'};

tt = linspace(0,100,2200);

base_length = 220;

figure('Position',[100 100 450 600])
for kk = 1:length(mep_traces)
    subplot(4,1,kk)
    if kk < 4
    temp = mep_traces{kk} -median(mep_traces{kk}) + 20*(rand - 0.5);
    else
    temp = mep_traces{kk};
    end
    %temp = mep_traces{kk};
    %plot(tt, mep_traces{kk} + y_offset*(kk-1), 'k');
    plot(tt, temp, 'k');
    hold on
    set(gca,'FontSize',12)
    title(emg_labels{kk},'FontSize',18)
    
    %xlim([0,100])
    ylim([-50, 50])
    if kk == length(mep_traces)
        ylim([-50, 50])
    end
    if kk == 4
        xlabel('Time (ms)','FontSize',16)
        ylim([-100,100])
    end


end
if save_figs
    exportgraphics(gcf,'/Users/ercole/Documents/mep_paper/figs/fig2a_p1.png','Resolution',300)
end

figure('Position',[700 100 450 600])
for kk = 1:length(mep_traces)
    subplot(4,1,kk)
    temptrace = mep_traces_corr{kk};
    if kk < 4
    baseline = median((temptrace(end-base_length:end)));
    sd = std(temptrace(end-base_length:end));
    temptrace = (temptrace -baseline)/sd;
    end
    % if kk == 4
    %     temptrace = temp;
    % end
    %temp = mep_traces{kk};
    %plot(tt, mep_traces{kk} + y_offset*(kk-1), 'k');
    plot(tt, temptrace, cols{kk});
    hold on
    set(gca,'FontSize',12)
    title(emg_labels{kk},'FontSize',18)
    
    if kk < 4
        xline(ltcs{kk});
        yline([thrs(kk), -thrs(kk)], 'k--')
    end


    %xlim([0,100])
    %ylim([-30, 30])
    % if kk < length(mep_traces)
    %     ylim([-75, 75])
    % else
    %     ylim([-100, 100])
    % end
    ylim([-lims(kk), lims(kk)])
    
    if kk == 4
        xlabel('Time (ms)','FontSize',16)
        ylim([-100,100])
    end
    
end
if save_figs
    exportgraphics(gcf,'/Users/ercole/Documents/mep_paper/figs/fig2a_p2.png','Resolution',300)
end

%% plot individual examples
tr = 4;
ch = 3;
pt = 30;

tr = 3; 
ch = 2; 
pt = 10;

% tr = 3;
% ch = 6;
% pt = 12;

% tr = 3;
% ch = 7;
% pt = 11;

% tr = 3;
% ch = 3;
% pt = 29;

% tr = 4;
% ch = 3;
% pt = 30;



tt = linspace(0,100,2200);

disp(size(ep_struct(pt).mep.raw.mep_means))

figure
plot(tt,ep_struct(pt).mep.raw.mep_means{tr, ch})
%ylim([-100, 100])

%% plot channels en masse to look for good examples (notch filt + clear stim artifact)
y_offset = 40;
x_offset = 110;

pt = 9;
chans = [1,2,3,4,5,6,7,8];
trials = 1:39;

if max(trials) > size(ep_struct(pt).mep.raw.mep_means,1)
    trials = trials(1):size(ep_struct(pt).mep.raw.mep_means,1);
end

figure('Position', [100 100 1000 600])
for k1 = 1:length(chans)
    for k2 = 1:length(trials)
        temptrace = ep_struct(pt).mep.raw.mep_means{trials(k2), chans(k1)};
        temptrace = (temptrace - median(temptrace(end-220:end)))/std(temptrace(end-220:end));
        plot(tt+ (k1-1)*x_offset,(k2-1)*y_offset+temptrace, 'k')
        hold on

    end
end
ylim([(min(trials)-1)*y_offset, (max(trials)+1)*y_offset])
yticks((trials-1)*y_offset)
yticklabels(trials)

disp(ep_struct(pt).mep.emg_labels(chans))

%% plot example for fig with/without notch filt
save_fig = false;

addpath('./functions')
pt = 7;
ch = 4;
trial = 32;

% temptrace = ep_struct(7).mep.raw.mep_means{26, 4};
% temptrace = (temptrace - median(temptrace))/std(temptrace(end-220:end));

figure
subplot(2,2,1)
temptrace = ep_struct(pt).mep.raw.mep_means{trial, ch};
temptrace = (temptrace - median(temptrace));%/std(temptrace(end-220:end));

plot(tt,temptrace, 'k')
hold on
ylim([-35,35])
set(gca,'FontSize',12)

title('Pulse-Averaged','FontSize',24)
ylabel('Voltage (\muV)', 'FontSize',18)

[yfilt, yest, A, phi] = td_notch(temptrace,tt/1000,60,22000);

plot(tt,yest,'b--')

subplot(2,2,2)
plot(tt,yfilt,'k')
ylim([-35,35])

set(gca,'FontSize',12)
title('Filtered','FontSize',24)

%xlabel('Time (ms)', 'FontSize',14)

pt = 7;
ch = 4;
trial = 26;

subplot(2,2,3)
temptrace = ep_struct(pt).mep.raw.mep_means{trial, ch};
temptrace = (temptrace - median(temptrace));%/std(temptrace(end-220:end));

plot(tt,temptrace, 'k')
set(gca,'FontSize',12)
hold on
ylim([-25,25])
xlabel('Time (ms)', 'FontSize',18)
ylabel('Voltage (\muV)', 'FontSize',18)

[yfilt, yest, A, phi] = td_notch(temptrace,tt/1000,60,22000);
%yfilt = yfilt/std(yfilt(end-220:end));

plot(tt,yest,'b--')
legend({'EMG','Sine fit'},'Location','Northeast')

subplot(2,2,4)
plot(tt,yfilt,'k')
ylim([-25,25])
set(gca,'FontSize',12)
xlabel('Time (ms)', 'FontSize',18)

if save_fig
    exportgraphics(gcf,'/Users/ercole/Documents/mep_paper/figs/fig2b_v1.png','Resolution',300)
    savefig('/Users/ercole/Documents/mep_paper/figs/fig2b_v1.fig')
end

%% plot example with/without notch filt

save_fig = true;

addpath('./functions')
pt = 7;
ch = 4;
trial = 32;

figure
subplot(1,2,1)
temptrace = ep_struct(pt).mep.raw.mep_means{trial, ch};
temptrace = (temptrace - median(temptrace))/std(temptrace(end-220:end));

plot(tt,temptrace, 'k')
hold on
ylim([-10,10])

title('Original')

[yfilt, yest, A, phi] = td_notch(temptrace,tt/1000,60,22000);

plot(tt,yest,'b--')

subplot(1,2,2)
plot(tt,yfilt,'k')
ylim([-10,10])
title('Filtered')

%% plot artifact rejection example

save_fig = false;
% 
% pt = 9;
% ch = 5;
% trial = 2;

% good mep + art
pt = 9;
ch = 5;
trial = 24;

% pt = 9;
% ch = 5;
% trial = 10;

% pt2 = 9;
% ch2 = 8;
% trial2 =25;

% example with small mep
% pt2 = 9;
% ch2 = 5;
% trial2 = 10;

%artifact example
pt2 = 9;
ch2 = 5;
trial2 = 27;

%artifact example 
% pt2 = 9;
% ch2 = 5;
% trial2 = 22;

load('MEP_ampstats_opt.mat')
%load('MEP_latencies_opt.mat')

facial_inds = [5,6,7,8,9]; limb_inds = [1 2 3 4 8 10];

sd_thresh_facial = 4; sd_thresh_limb = 7;

latency_labels = {'bicep','ECR','FCR','FDI','nasalis','orb oris','geniogl','trap','ocul','tib'};

use_saved_latencies = false;
if use_saved_latencies
    latencies = [];
    latency_tolerance = 1;
else
    %original latencies;
    %latencies = [16.92,24.56;16.59,29.24;18.31,28.69;24,26;8,13;8.17,14.21; 8,13; 14.71,14.71; 8,13; 16.92,24.56];
    %latency_tolerance = 3;
    
    %initial wide latency ranges to then determine final latencies in a
    %data-driven way
    %latencies = [10,40;10,40;10,40;10,40;4,20;4,20;4,20; 4,20; 10,40; 4,20; 10,40];

    % optimized data-driven latencies after running on whole data set with
    % wide range LTC

    latencies = [11,30.5; 16,34; 14,32; 19.5,30.5; 5,18.5; 5,18.5; 5,17; 4,20; 14,32; 5,18.5; 14,32];
    latency_tolerance = 1;
    latencies(:,1) = latencies(:,1) - latency_tolerance;
    latencies(:,2) = latencies(:,2) + latency_tolerance;
end

mep_settings = struct();
mep_settings.art_window = [3,25];
mep_settings.art_detect = [3,6]; %window to check whether to apply art correction
mep_settings.art_thresh = 3;
mep_settings.sd_thresh = sd_thresh_facial;
mep_settings.fs = 22000;
mep_settings.base_length = 10;   %round(10*mep_settings.fs/1000);
mep_settings.pctl_thresh = 50;
mep_settings.use_dbs = false;
mep_settings.notch_filt_emg = true;
mep_settings.art_type = 'template';   %whether to fit artifact with PCA on EMG, LFP, or both chanels, or use 'template'
mep_settings.art_chans = 'facial'; %which chans to apply correction to
mep_settings.fit_type = 'median'; %how to fit/remove artifact template for individual channels

base_length = 220;
tt = linspace(0,100,2200);

if strcmp(mep_settings.art_type, 'template') || strcmp(mep_settings.art_type, 'template_pca')
    load('EMG_artifacts_templib.mat')
    templib = cell2mat({EMG_save(:).data}');
    for kl = 1:size(templib,1)
        temptrace = templib(kl,:);

        baseline = median((temptrace(end-base_length:end)));
        sd = std(temptrace(end-base_length:end));
        temptrace = (temptrace -baseline)/sd;

        templib(kl,:) = temptrace;
    end
    %freq=300;
    order=2;

    %framelen=round(1/freq*fs);
    framelen = 51;
    if mod(framelen,2)==0
        framelen=framelen+1;
    end
    templib = transpose(sgolayfilt(double(templib'),order,framelen));
    templib = [templib; -templib];  %flip traces to accommodate both negative and positive shapes
    mep_settings.templib = templib;
end

exc_baseline = [1100,2200];
mean_thresh = 10; var_thresh = 50;

key_ind = 5;

ltc_to_plot = latencies(key_ind,:);
latencies_pca = repmat(ltc_to_plot,[length(ep_struct(pt).mep.emg_labels) 1]);

mep_settings.is_facial = true(size(ep_struct(pt).mep.emg_labels));

mep_means = ep_struct(pt).mep.raw.mep_means(trial,:)';
%dlep_means =  ep_struct(pt_ind).dlep.raw_means(ko,:)';
dlep_means = [];
[detected, amps, delays, amps_z, corr_traces, templates_temp] = detect_mep_rejpca(mep_means, dlep_means, latencies_pca, mep_settings);
pc_template = templates_temp{ch};
corr_trace = corr_traces{ch};

disp('Done')

t_art = linspace(mep_settings.art_window(1),mep_settings.art_window(2),length(pc_template));

figure
subplot(2,2,1)
temptrace = mep_means{ch};
temptrace = (temptrace - median(temptrace(70:440)))/std(temptrace(end-220:end));

mep_traces{3} = temptrace;

plot(tt,temptrace,'k')
hold on
set(gca,'FontSize',12)
title('Pulse-Averaged','FontSize',24)
ylabel('Voltage (Z-scored)', 'FontSize',20)

plot(t_art,pc_template-median(pc_template),'b--')

%yline(-4, 'k--'); yline(4, 'k--');

ylim([-25 25])
xlim([-1 50])

subplot(2,2,2)
plot(tt,corr_traces{ch},'k') 
mep_traces_corr{3} = corr_traces{ch};
ylim([-25 25])
xlim([-1 50])
yline(-4, 'k--'); yline(4, 'k--'); 
xline(4, 'k--'); xline(19, 'k--'); 
set(gca,'FontSize',12)
title('Corrected','FontSize',24)
%ylabel('Voltage (Z-scored)', 'FontSize',14)
mep_means = ep_struct(pt2).mep.raw.mep_means(trial2,:)';
%dlep_means =  ep_struct(pt_ind).dlep.raw_means(ko,:)';
dlep_means = [];
[detected, amps, delays, amps_z, corr_traces, templates_temp] = detect_mep_rejpca(mep_means, dlep_means, latencies_pca, mep_settings);
pc_template = templates_temp{ch2};
corr_trace = corr_traces{ch2};

subplot(2,2,3)

temptrace = mep_means{ch2};
baseline_sd = 4*std(temptrace(end-220:end));
%temptrace = temptrace - median(temptrace(70:440));
temptrace = (temptrace - median(temptrace(70:440)))/std(temptrace(end-220:end));
plot(tt,temptrace,'k')
hold on
plot(t_art,pc_template-median(pc_template),'b--')
%yline(-4, 'k--'); yline(4, 'k--');
set(gca,'FontSize',12)
legend({'EMG','Artifact fit','',''},'Location','Northeast')
xlabel('Time (ms)', 'FontSize',18)
ylabel('Voltage (Z-scored)', 'FontSize',18)


ylim([-25 25])
xlim([-1 50])


subplot(2,2,4)
plot(tt,corr_traces{ch2},'k')
set(gca,'FontSize',12)
xlabel('Time (ms)', 'FontSize',18)
%ylabel('Voltage (Z-scored)', 'FontSize',14)
hold on
ylim([-25 25])
xlim([-1 50])
yline(-4, 'k--'); yline(4, 'k--'); 
xline(4, 'k--'); xline(19, 'k--'); 

if save_fig
    exportgraphics(gcf,'/Users/ercole/Documents/mep_paper/figs/fig2d_v1.png','Resolution',300)
    savefig('/Users/ercole/Documents/mep_paper/figs/fig2d_v1.fig')
end