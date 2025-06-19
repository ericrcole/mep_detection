%mep_amp_rescale_fig.m
%
% extension of tune_mep_ampvals to plot figure showing mep amplitude
% differences between muscles, and amplitude distributions after muscle
% re-scaling

load('EP_102924_final_STN.mat');

%% get all unique muscle names
emg_labels_all = {};

for kk = 1:length(ep_struct)
    emg_labels_all = [emg_labels_all, ep_struct(kk).mep.emg_labels];
end

disp(unique(emg_labels_all)')

%%
emg_keys = {'bicep','ecr','fcr','fdi','nasal','ocul','oris','genio'};
emg_names = {'Biceps','ECR','FCR','FDR','Nasalis','Orb. Oculi','Orb. Oris','Genioglossus'};
facial_keys = {'nasal','ocul','oris','genio'};
limb_keys = {'bicep','ecr','fcr','fdi','tib','trap'};

keys_facial = [0 0 0 0 1 1 1 1];

mode = 'opt';
exc_outliers = false; 

base_length = 220;

mep_settings = struct();
mep_settings.art_window = [3,20];
mep_settings.art_detect = [3,6]; %window to check whether to apply art correction
mep_settings.art_thresh = 3;
mep_settings.sd_thresh = 4;
mep_settings.fs = 22000;
mep_settings.base_length = 10;   %round(10*mep_settings.fs/1000);
mep_settings.pctl_thresh = 50;
mep_settings.use_dbs = true;
mep_settings.notch_filt_emg = true;
mep_settings.art_type = 'template';   %whether to fit artifact with PCA on EMG, LFP, or both chanels, or use 'template'
mep_settings.art_chans = 'facial'; %which chans to apply correction to
mep_settings.fit_type = 'median'; %how to fit/remove artifact template for individual channels

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

emg_peak_thresh = 50;

outlier_trials = [];
outlier_muscles = [];
outlier_pts = [];
outlier_amps = [];

all_amps_by_muscle = cell(size(emg_keys));
all_ltcs_by_muscle = cell(size(emg_keys));

tic

for kk = 1:length(ep_struct)
    if kk == 5 || (ep_struct(kk).pt_info.impute_emg)
        continue
    end

    temp_str = struct();

    labels_detect = -1*ones(size(ep_struct(kk).mep.(mode).labels_detect));
    amps_detect = cell(size(labels_detect));
    amps_detect_z = ep_struct(kk).mep.(mode).amps_detect_z;
    delays_detect = ep_struct(kk).mep.(mode).delays_detect;
    mep_means = ep_struct(kk).mep.raw.mep_means;
    
    amps_max = nan(size(labels_detect));

    for kl = 1:length(ep_struct(kk).mep.emg_labels)
        if contains(ep_struct(kk).mep.emg_labels,'hyo')
            ep_struct(kk).mep.emg_labels{kl} = 'genioglossus';
        end

        key_ind = 0;
        for km = 1:length(emg_keys)
            if contains(lower(ep_struct(kk).mep.emg_labels{kl}), emg_keys{km})
                key_ind = km;
                break
            end
        end
        if key_ind == 0
            continue
        end

        for km = 1:size(labels_detect,1)
            
%             tempmat = ep_struct(ptind).mep.raw.mep_means(:,kl);
%             tempmat = double(cell2mat(tempmat(~cellfun(@isempty,tempmat))));
%             for km = 1:size(tempmat,1)
%                 tempmat(km,:) = tempmat(km,:) - median(tempmat(km,end-220:end));
%             end
            temptrace = ep_struct(kk).mep.raw.mep_means{km,kl};
            temptrace = temptrace - median(temptrace(end-220:end));
            
            exc_flag = false;
            if exc_outliers
                if abs(mean(temptrace(exc_baseline(1):exc_baseline(2)))) > mean_thresh
                    exc_flag = true;
                end
                if var(temptrace(exc_baseline(1):exc_baseline(2))) > var_thresh
                    exc_flag = true;
                end
            end

            if exc_flag
                labels_detect(km,kl) = nan;
                delays_detect{km,kl} = [];
                amps_detect_z{km,kl} =[];
                amps_detect{km,kl} = [];
                %mep_means{km,kl} = [];
            else

            end
            if ~isempty(amps_detect_z{km,kl})
                amps_max(km,kl) = max(amps_detect_z{km,kl});
                tempamp = max(amps_detect_z{km,kl})';
                templtc = delays_detect{km,kl}';

                all_amps_by_muscle{key_ind} = [all_amps_by_muscle{key_ind}; tempamp];
                all_ltcs_by_muscle{key_ind} = [all_ltcs_by_muscle{key_ind}; templtc];

                if max(amps_detect_z{km,kl}) > emg_peak_thresh
                    outlier_trials = [outlier_trials; km];
                    outlier_muscles = [outlier_muscles; kl];
                    outlier_pts = [outlier_pts; kk];
                    outlier_amps = [outlier_amps; tempamp];

                    fprintf('%s: trial %d, muscle %d, %.1f\n', ep_struct(kk).patient_ID, km, kl, tempamp)
                end
            end
        end
    end
    
end

toc

%% plot specific mEP example
% inds = [1;2;3;4;5;7;8;6];
% 
% tt = linspace(0,100,2200);
% 
% figure;
% for kk = 1:8
%     ind = inds(kk);
% 
%     ptind = outlier_pts(ind);
%     km = outlier_trials(ind);
%     kl = outlier_muscles(ind);
% 
%     subplot(4,2,kk)
%     temptrace = ep_struct(ptind).mep.raw.mep_means{km,kl};
%     plot(tt,temptrace - median(temptrace));
%     xlabel('time (ms)')
%     ylabel('Voltage')
%     ylim([-300 300])
% 
%     title(sprintf('Amplitude = %.2f',outlier_amps(ind)))
% end
% 
% sgtitle('Greatest outlier trials')

%% Now, get amplitude distributions for different muscles

save_amp_dist = false;

amp_dist_log = [];
amp_dist_z = [];

figure
subplot(2,1,1)
categs = emg_names;
for kk = 1:length(emg_keys)
    temp = abs(all_amps_by_muscle{kk});
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    amp_dist_z = [amp_dist_z; prctile(temp,25), prctile(temp,75), nanmin(temp)];
end
xticks(1:length(categs))
xticklabels(categs)
ylim([0 50])
ylabel('Amplitude (z-score)')
%title('Raw Z-amplitudes')

% subplot(2,2,2)
% for kk = 1:length(emg_keys)
%     temp = abs(all_amps_by_muscle{kk});
%     if keys_facial(kk)
%         temp(temp<5) = 5;
%         temp = log2(1+temp - 5)+3;
%     else
%         temp(temp<8) = 8;
%         temp = log2(1+temp - 8)+3;
%     end
%     
%     %temp = log2(1+abs(all_amps_by_muscle{kk}));
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
% 
%     fprintf('%s: median = %.2f, SD = %.2f\n',emg_keys{kk},nanmedian(temp), nanstd(temp))
%     %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
% 
%     amp_dist_log = [amp_dist_log; nanmean(temp), nanstd(temp)];
% end
% fprintf('\n')
% xticks(1:length(categs))
% xticklabels(categs)
% xlabel('Muscle')
% ylabel('Amplitude (log-scaled)')
% title('Log-scaled Z-amplitudes')


amp_dist = amp_dist_log;
mean_stats = mean(amp_dist);

stat_adj_log = zeros(size(amp_dist)); %differences needed to adjust each muscle to target mean distribution
for kl = 1:size(stat_adj_log,1)
    stat_adj_log(kl,1) = mean_stats(1) - amp_dist(kl,1);
    stat_adj_log(kl,2) = mean_stats(2)./amp_dist(kl,2);
end

mep_amp_stats = struct();
mep_amp_stats.keys = emg_keys;
mep_amp_stats.mean_stats_log = mean_stats;
mep_amp_stats.all_stats_log = amp_dist;
mep_amp_stats.stat_adj_log = stat_adj_log;

amp_dist = amp_dist_z;
mean_stats = mean(amp_dist);

stat_adj_z = zeros(size(amp_dist,1),1); %differences needed to adjust each muscle to target mean distribution
for kl = 1:size(stat_adj_z,1)
    stat_adj_z(kl) = (mean_stats(2) - mean_stats(1))./(amp_dist(kl,2) - amp_dist(kl,1));
end

mep_amp_stats.mean_stats_z = mean_stats;
mep_amp_stats.all_stats_z = amp_dist;
mep_amp_stats.stat_adj_z = stat_adj_z;

if save_amp_dist
    save('MEP_ampstats_opt.mat','mep_amp_stats')
end


subplot(2,1,2)
categs = emg_names;
for kk = 1:length(emg_keys)
    temp = abs(all_amps_by_muscle{kk});
    temp = temp - amp_dist_z(kk,3);
    
    temp = temp*stat_adj_z(kk);
    temp = temp + 8;

    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    
end
xticks(1:length(categs))
xticklabels(categs)
xlabel('Muscle')
ylabel('Amplitude (normalized)')
%title('Adjusted Z-amplitudes')
ylim([0 50])


% subplot(2,2,4)
% for kk = 1:length(emg_keys)
%     temp = abs(all_amps_by_muscle{kk});
%     if keys_facial(kk)
%         temp(temp<5) = 5;
%         temp = log2(1+temp - 5)+3;
%     else
%         temp(temp<8) = 8;
%         temp = log2(1+temp - 8)+3;
%     end
%     
%     tempmean = nanmean(temp);
%     temp = temp - tempmean;
%     temp = temp*stat_adj_log(kk,2) + tempmean+stat_adj_log(kk,1);
% 
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
% 
%     %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
% 
% end
% fprintf('\n')
% xticks(1:length(categs))
% xticklabels(categs)
% xlabel('Muscle')
% ylabel('Amplitude (adjusted)')
% title('Adjusted log-scaled Z-amplitudes')

%% Latencies

save_ltc = false;

latencies_to_save = [];

figure
for kk = 1:length(emg_keys)
    temp = abs(all_ltcs_by_muscle{kk});
    boxchart(ones(size(temp))*kk,temp);
    hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    fprintf('%s: %.2f - %.2f ms\n', emg_keys{kk}, nanmin(temp),nanmax(temp))

    latencies_to_save = [latencies_to_save; nanmin(temp),nanmax(temp)];
end
fprintf('\n')
xticks(1:length(categs))
xticklabels(categs)
xlabel('Muscle')
ylabel('Latency')

if save_ltc
    latencies_saved = struct();
    latencies_saved.latencies = latencies_to_save;
    latencies_saved.emg_labels = emg_keys;
    
    save('MEP_latencies_opt.mat','latencies_saved')
end

%%
total_count = 0;
for kk = 1:length(ep_struct)
    total_count = total_count + length(ep_struct(kk).stim_settings.param_strings);
end
disp(total_count)
