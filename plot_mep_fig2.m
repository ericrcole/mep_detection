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

pt = 9;
trial = 18;
chans = [1,2,5,8];

ltcs = {[14,26],[14,32],[5 16], [5,16]};
thrs = [7,7,4,4];

cols = {'k','b','b','r'};

disp(ep_struct(pt).stim_settings.param_strings{trial})

y_offset = 300;

mep_traces = ep_struct(pt).mep.raw.mep_means(trial,chans);
emg_labels = ep_struct(pt).mep.emg_labels;

tt = linspace(0,100,2200);

base_length = 220;

figure
for kk = 1:length(mep_traces)
    subplot(4,1,kk)
    temptrace = mep_traces{kk};
    baseline = median((temptrace(end-base_length:end)));
    sd = std(temptrace(end-base_length:end));
    temptrace = (temptrace -baseline)/sd;
    if kk == 4
        temptrace = temp;
    end
    %temp = mep_traces{kk};
    %plot(tt, mep_traces{kk} + y_offset*(kk-1), 'k');
    plot(tt, temptrace, cols{kk});
    hold on
    title(emg_labels{chans(kk)})
    
    if kk < 4
        xline(ltcs{kk});
        yline([thrs(kk), -thrs(kk)], 'k--')
    end


    %xlim([0,100])
    ylim([-30, 30])
    if kk == length(mep_traces)
        ylim([-50, 50])
    end
    if kk == 4
        xlabel('Time (ms)')
    end
end

%% plot individual examples
tr = 2;
ch = 7;
pt = 11;

tt = linspace(0,100,2200);

disp(size(ep_struct(pt).mep.raw.mep_means))

figure
plot(tt,ep_struct(pt).mep.raw.mep_means{tr, ch})
ylim([-150, 150])