%% quantify effectiveness of different artifact rejection strategies:
% this version: evaluates frequency and time-domain based notch filtering

%%
categs = 0:0.05:0.95;   %corrcoef thresholds

addpath('C:\Users\ORStream\Documents\mEP_processing\functions')

tt = linspace(0,0.1,2200);

rem_art_pca = true;
use_tdnotch = false;
use_frnotch = false;
filt_sg = false;
rem_ringing_art = false;
filt_hp = false;
filt_lp = false;

mode = 'raw';

mep_settings = struct();
mep_settings.art_window = [3,20];
mep_settings.art_detect = [3,6]; %window to check whether to apply art correction
mep_settings.art_thresh = 3;
mep_settings.sd_thresh = 4;
mep_settings.fs = 22000;
mep_settings.base_length = 10;   %round(10*mep_settings.fs/1000);
mep_settings.pctl_thresh = 50;
mep_settings.use_dbs = true;
mep_settings.notch_filt_emg = false;
mep_settings.art_type = 'template';   %whether to fit artifact with EMG, LFP, or both, or use pre-saved templates
mep_settings.art_chans = 'facial'; %which chans to apply correction to
mep_settings.fit_type = 'median'; %how to fit/remove artifact template for individual channels

if strcmp(mep_settings.art_type, 'template')
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
    templib = [templib; -templib];
    mep_settings.templib = templib;
end


acc_ampsweep = cell(length(categs),3);
prec_ampsweep = cell(length(categs),3);
rec_ampsweep = cell(length(categs),3);
f1_ampsweep = cell(length(categs),3);

for ampind = 1:length(categs)

    corr_thresh = categs(ampind);

    acc_by_muscle = cell(size(emg_mapping_labels));
    prec_by_muscle = cell(size(emg_mapping_labels));
    rec_by_muscle = cell(size(emg_mapping_labels));
    f1_by_muscle = cell(size(emg_mapping_labels));

    for pt_ind = 1:length(ep_struct)
        [facial_log,facial_inds] = ismember(emg_mapping_all{pt_ind},is_facial);
        mep_settings.is_facial = logical(facial_log);
        for kp = 1:length(emg_mapping_all{pt_ind})
            n_trials = length(ep_struct(pt_ind).stim_settings.param_strings);
            detected_temp = zeros(n_trials,1);
            ltc_to_plot = latencies(emg_mapping_all{pt_ind}(kp),:) + [-latency_tolerance, latency_tolerance];
            
            temptrace_mat = cell(n_trials, 1);

            for ko = 1:n_trials
                temptrace = double(ep_struct(pt_ind).mep.(mode).mep_means{ko,kp});
                if isempty(temptrace)
                    temptrace_mat{ko} = zeros(size(tt));
                    continue
                end
                if filt_sg
                    freq=150;
                    order=2;

                    %framelen=round(1/freq*fs);
                    framelen = 25;
                    if mod(framelen,2)==0
                        framelen=framelen+1;
                    end
                    temptrace = sgolayfilt(temptrace,order,framelen);

                end
                if rem_ringing_art
                    [temptrace, artifact_fit] = ep_remove_ringing_artifact_erc((temptrace), [.005 .04], fs);
                end
                if filt_hp
                    temptrace = filtfilt(b_high, a_high, temptrace);
                end
                if filt_lp
                    temptrace = filtfilt(b_low, a_low, temptrace);
                end
                if use_frnotch
                    temptrace = filtfilt(b_notch, a_notch, temptrace);
                end

                if use_tdnotch
                    %temptrace = tdnotch(temptrace);
                    [temptrace, yest, A, phi] = td_notch(temptrace,tt,60,fs);
                end

                if ~rem_art_pca
                    temptrace_mat{ko} = temptrace; 
                end

                if rem_art_pca
                    mep_means = ep_struct(pt_ind).mep.raw.mep_means(ko,:)';
                    dlep_means =  ep_struct(pt_ind).dlep.raw_means(ko,:)';
                    latencies_pca = latencies(emg_mapping_all{pt_ind},:) + [-latency_tolerance, latency_tolerance];
                    [detected, amps, delays, amps_z, corr_traces, templates] = detect_mep_rejpca(mep_means, dlep_means, latencies_pca, mep_settings);
                    pc_template = templates{kp};

                    corr_trace = corr_traces{kp};
                    temptrace_mat{ko} = corr_trace;

                    t_art_samps = round(mep_settings.art_window*fs/1000);
                    detected_temp(ko) = detected(kp);

                    baseline = median((temptrace(end-base_length:end)));
                    sd = std(temptrace(end-base_length:end));
                    temptrace = (temptrace -baseline)/sd;
                    latency_pts = round(ltc_to_plot/1000*fs);
                else
%                     baseline = median((temptrace(end-base_length:end)));
%                     sd = std(temptrace(end-base_length:end));
%                     temptrace = (temptrace -baseline)/sd;
%                     latency_pts = round(ltc_to_plot/1000*fs);
                end
                %
                if ~rem_art_pca
                %[pkval, loc, widths, prom] = findpeaks(abs(temptrace(latency_pts(1):latency_pts(2))),'MinPeakHeight',sd_thresh,'MaxPeakWidth',fs*0.015,'MinPeakWidth',fs*0.0015,'MinPeakProminence',sd_thresh);
                %[detected, amps, delays, amps_z] = ep_thresh_windowed(temptrace,10,mep_settings.sd_thresh,fs,ltc_to_plot);
                %         delay_temp = ep_struct(pt_ind).mep.delays_detect{ko,kp};
                %         amps_temp = ep_struct(pt_ind).mep.amps_detect_z{ko,kp};
% %                 amps_temp = nan(size(pkval));
% %                 delays_temp = nan(size(pkval));
% %                 amps_z_temp = nan(size(pkval));
% % 
% %                 if ~isempty(pkval)
% %                     detected = true;
% %                     p_widths{ko,kp} = widths;
% %                     p_prominences{ko,kp} = prom;
% %                     for kl = 1:length(pkval)
% %                         amps_temp(kl) = pkval(kl);
% %                         delays_temp(kl) = x_plot(loc(kl))+ltc_to_plot(1);
% %                         amps_z_temp(kl) =  amps_temp(kl);
% %                     end
% %                     p_locs{ko,kp} = delays_temp;
% %                     p_amps_z{ko,kp} = amps_z_temp;
% %                 else
% %                     detected = false;
% %                 end
                %detected_temp(ko) = detected;
                end

               
            end
            
            temptrace_mat = cell2mat(temptrace_mat);
            
            tic
            [detected_temp, sel_inds, corr_vals] = detect_mep_corr(temptrace_mat, corr_thresh, ltc_to_plot, fs);
            toc

            acc_by_muscle{emg_mapping_all{pt_ind}(kp)} = [acc_by_muscle{emg_mapping_all{pt_ind}(kp)};  mean(detected_temp == mep_labels_sm{pt_ind}(:,kp))];
            acctemp = mean(detected_temp == mep_labels_sm{pt_ind}(:,kp));
            cmat = confusionmat(logical(detected_temp),mep_labels_sm{pt_ind}(:,kp));
            nanflag = false;
            if (size(cmat,1) == 1) && (mean(detected_temp)==0)
                cmat = [cmat 0; 0 0];
                nanflag = true;
            elseif (size(cmat,1) == 1) && (mean(detected_temp)==1)
                cmat = [0 0; 0 cmat];
            end
            
            if nanflag
                precision = 1;
                recall = 1;
                f1 = 1;

                prec_by_muscle{emg_mapping_all{pt_ind}(kp)} = [prec_by_muscle{emg_mapping_all{pt_ind}(kp)}; precision];
                rec_by_muscle{emg_mapping_all{pt_ind}(kp)} = [rec_by_muscle{emg_mapping_all{pt_ind}(kp)}; recall];
            else
                   
                precision = cmat(2,2) ./ (cmat(2,2)+cmat(2,1));
                recall = cmat(2,2) ./ (cmat(2,2)+cmat(1,2));
    
                prec_by_muscle{emg_mapping_all{pt_ind}(kp)} = [prec_by_muscle{emg_mapping_all{pt_ind}(kp)}; precision];
                rec_by_muscle{emg_mapping_all{pt_ind}(kp)} = [rec_by_muscle{emg_mapping_all{pt_ind}(kp)}; recall];
    
                f1 = 2*(precision.*recall)/(precision+recall);
            end

            if isnan(f1) && ((precision == 0) && (recall == 0))
                f1 = 0;
            end

            f1_by_muscle{emg_mapping_all{pt_ind}(kp)} = [f1_by_muscle{emg_mapping_all{pt_ind}(kp)}; f1];

        end

    end
    DD = acc_by_muscle;
    is_facial = [5,6,7,8,11,12,13];
    is_limb = [1,2,3,4,9,10];
    q1 =DD(is_facial); q2 = DD(is_limb);
    acc_ampsweep(ampind,:) = {cell2mat(DD(:)),cell2mat(q1(:)),cell2mat(q2(:))};

    DD = prec_by_muscle;
    is_facial = [5,6,7,8,11,12,13];
    is_limb = [1,2,3,4,9,10];
    q1 =DD(is_facial); q2 = DD(is_limb);
    prec_ampsweep(ampind,:) = {cell2mat(DD(:)),cell2mat(q1(:)),cell2mat(q2(:))};

    DD = rec_by_muscle;
    is_facial = [5,6,7,8,11,12,13];
    is_limb = [1,2,3,4,9,10];
    q1 =DD(is_facial); q2 = DD(is_limb);
    rec_ampsweep(ampind,:) = {cell2mat(DD(:)),cell2mat(q1(:)),cell2mat(q2(:))};

    DD = f1_by_muscle;
    is_facial = [5,6,7,8,11,12,13];
    is_limb = [1,2,3,4,9,10];
    q1 =DD(is_facial); q2 = DD(is_limb);
    f1_ampsweep(ampind,:) = {cell2mat(DD(:)),cell2mat(q1(:)),cell2mat(q2(:))};

end

disp('Done')

%%

acc_ampsweep_med = cellfun(@nanmean,acc_ampsweep);
rec_ampsweep_med = cellfun(@nanmean,rec_ampsweep);
prec_ampsweep_med = cellfun(@nanmean,prec_ampsweep);
f1_ampsweep_med = cellfun(@nanmean,f1_ampsweep);

% acc_ampsweep_med = cellfun(@(x) prctile(x, 10),acc_ampsweep);
% rec_ampsweep_med = cellfun(@(x) prctile(x, 10),rec_ampsweep);
% prec_ampsweep_med = cellfun(@(x) prctile(x, 10),prec_ampsweep);
% f1_ampsweep_med = cellfun(@(x) prctile(x, 10),f1_ampsweep);

sem = @(x) (std(x)./sqrt(length(x)));
acc_ampsweep_sem = cellfun(sem,acc_ampsweep);
rec_ampsweep_sem = cellfun(sem,rec_ampsweep);
prec_ampsweep_sem = cellfun(sem,prec_ampsweep);
f1_ampsweep_sem = cellfun(sem,f1_ampsweep);

figure('Position',[100 100 1000 600])
subplot(2,2,1)
plot(categs, acc_ampsweep_med(:,2))
ylabel('Accuracy')
xlabel('Corr threshold')
title('Facial')

subplot(2,2,2)
plot(categs, f1_ampsweep_med(:,2))
ylabel('Accuracy')
xlabel('Corr threshold')
ylabel('F1')
xlabel('Filter type')
title('Facial')

subplot(2,2,3)
plot(categs, acc_ampsweep_med(:,3))
ylabel('Accuracy')
xlabel('Corr threshold')
title('Limb')

subplot(2,2,4)
plot(categs, f1_ampsweep_med(:,3))
ylabel('Accuracy')
xlabel('Corr threshold')
ylabel('F1')
xlabel('Filter type')
title('Limb')

sgtitle('Correlation-based detection')

cols = {'r','b','k'};
labels = {'All','Facial','Limb'};

figure
subplot(2,2,1)
hold on
plot(categs,acc_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
plot(categs,acc_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
plot(categs,acc_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
ylim([0 1.05])
%title(sprintf('Accuracy: %s', mode))
title('Accuracy')
xlabel('Threshold')
legend(labels,'Location','Southeast')

subplot(2,2,2)
hold on
plot(categs,prec_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
plot(categs,prec_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
plot(categs,prec_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
%title(sprintf('Precision: %s', mode))
title('Precision')
xlabel('Threshold')
ylim([0 1.05])

subplot(2,2,3)
hold on
plot(categs,rec_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
plot(categs,rec_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
plot(categs,rec_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
%title(sprintf('Recall: %s',mode))
title('Recall')
xlabel('Threshold')
ylim([0 1.05])

subplot(2,2,4)
hold on
plot(categs,f1_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
plot(categs,f1_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
plot(categs,f1_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
%title(sprintf('F1: %s', mode))
title('F1')
xlabel('Threshold')
ylim([0 1.05])

sgtitle('Correlation-based detection')

% %%
% 
% figure('Position',[100 100 1000 600])
% subplot(1,2,1)
% for kk = 1:length(categs)
%     temp = acc_ampsweep{kk,2};
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
%     scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
% end
% xticks(1:length(categs))
% xticklabels(categ_names)
% ylabel('Accuracy')
% xlabel('Filter type')
% 
% subplot(1,2,2)
% for kk = 1:length(categs)
%     temp = f1_ampsweep{kk,2};
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
%     scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp,20, [0.5 0.5 0.5],'filled')
% end
% xticks(1:length(categs))
% xticklabels(categ_names)
% ylabel('F1')
% xlabel('Filter type')
% 
% sgtitle('Notch filtering')
% 
% %%
% 
% acc_ampsweep_med = cellfun(@nanmean,acc_ampsweep);
% rec_ampsweep_med = cellfun(@nanmean,rec_ampsweep);
% prec_ampsweep_med = cellfun(@nanmean,prec_ampsweep);
% f1_ampsweep_med = cellfun(@nanmean,f1_ampsweep);
% 
% % acc_ampsweep_med = cellfun(@(x) prctile(x, 10),acc_ampsweep);
% % rec_ampsweep_med = cellfun(@(x) prctile(x, 10),rec_ampsweep);
% % prec_ampsweep_med = cellfun(@(x) prctile(x, 10),prec_ampsweep);
% % f1_ampsweep_med = cellfun(@(x) prctile(x, 10),f1_ampsweep);
% 
% sem = @(x) (std(x)./sqrt(length(x)));
% acc_ampsweep_sem = cellfun(sem,acc_ampsweep);
% rec_ampsweep_sem = cellfun(sem,rec_ampsweep);
% prec_ampsweep_sem = cellfun(sem,prec_ampsweep);
% f1_ampsweep_sem = cellfun(sem,f1_ampsweep);
% 
% figure('Position',[100 100 1000 600])
% subplot(2,2,1)
% for kk = 1:length(categs)
%     temp = acc_ampsweep{kk,2};
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
%     scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
% end
% xticks(1:length(categs))
% xticklabels(categ_names)
% ylabel('Accuracy')
% xlabel('Filter type')
% title('Facial')
% 
% subplot(2,2,2)
% for kk = 1:length(categs)
%     temp = f1_ampsweep{kk,2};
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
%     scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp,20, [0.5 0.5 0.5],'filled')
% end
% xticks(1:length(categs))
% xticklabels(categ_names)
% ylabel('F1')
% xlabel('Filter type')
% title('Facial')
% 
% subplot(2,2,3)
% for kk = 1:length(categs)
%     temp = acc_ampsweep{kk,3};
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
%     scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
% end
% xticks(1:length(categs))
% xticklabels(categ_names)
% ylabel('Accuracy')
% xlabel('Filter type')
% title('Limb')
% 
% subplot(2,2,4)
% for kk = 1:length(categs)
%     temp = f1_ampsweep{kk,3};
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
%     scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp,20, [0.5 0.5 0.5],'filled')
% end
% xticks(1:length(categs))
% xticklabels(categ_names)
% ylabel('F1')
% xlabel('Filter type')
% title('Limb')
% 
% sgtitle('Correlation-based detection')


function [detected, sel_inds, corr_vals] = detect_mep_corr(dataMatrix, thresh, latency, fs)
    %dataMatrix: N_trials x time matrix of EMG traces from same EMG channel
    %threshold: minimum expected correlation between trials (0-1)
    %latency: time subset for expected muscle latency response
    %
    %finds submatrix that is above correlation threshold, to find subset of
    %mEP trials with high inter-correlation
    %
    %Algo to do this in linear time: first, sort all trials by mean
    %correlation with all other trials;
    %then find the biggest contiguous submatrix after sorting,
    %and iterate through remaining trials to see if they satisfy the
    %threshold when appended to the submatrix

    latency_pts = round(latency/1000*fs);
    
    % Compute the correlation coefficients for the selected time points
    subsetData = dataMatrix(:, latency_pts(1):latency_pts(2))'; % Extract the subset of the matrix
    corrMatrix = corrcoef(subsetData); % Calculate the correlation coefficients
    
    N = size(corrMatrix, 1);
    meanCorrelations = nanmean(corrMatrix - eye(N), 2); % Subtract self-correlation (which is 1)
    [corr_means, corr_inds] = sort(meanCorrelations,'descend');
    corrMatrix_sort = corrMatrix(corr_inds, corr_inds);

    stopflag = false;  %after finding submatrix, switch to append mode

    contig_inds = false(size(corrMatrix_sort,1),1);
    for k = 2:N
        %heuristic to find correlated submatrix in linear time
        submat_temp = corrMatrix_sort(1:k,1:k);
        if all(nanmean(submat_temp)>thresh) && ~stopflag
            %contiguous matrix is above threshold
            
            %sel_inds = corr_inds(1:k);
            contig_inds(1:k) = true;
        elseif (k == 2) && ~(all(nanmean(submat_temp)>thresh))
            %no correlated trials, return empty result (no mEP)
            sel_inds = []; corr_vals = [];
        else
            stopflag = true;  %added trial is non-correlated; switch to filling in new trials
        end

        if stopflag
            tempsel = contig_inds; tempsel(k) = true;
            submat_temp = corrMatrix_sort(tempsel,tempsel);
            if all(nanmean(submat_temp)>thresh)
                contig_inds = tempsel;
            end
        end
    end
    sel_inds = corr_inds(contig_inds);
    detected = false(N,1);
    detected(sel_inds) = true;
    corr_vals = mean(corrMatrix_sort(contig_inds,contig_inds));
end

% figure
% subplot(2,2,1)
% hold on
% plot(amp_sweep,acc_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
% plot(amp_sweep,acc_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
% plot(amp_sweep,acc_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
% ylim([0 1.05])
% %title(sprintf('Accuracy: %s', mode))
% title('Recall')
% xlabel('Threshold (z-score)')
% legend(labels,'Location','Southeast')
% 
% subplot(2,2,2)
% hold on
% plot(amp_sweep,prec_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
% plot(amp_sweep,prec_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
% plot(amp_sweep,prec_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
% %title(sprintf('Precision: %s', mode))
% title('Precision')
% xlabel('Threshold (z-score)')
% ylim([0 1.05])
% 
% subplot(2,2,3)
% hold on
% plot(amp_sweep,rec_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
% plot(amp_sweep,rec_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
% plot(amp_sweep,rec_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
% %title(sprintf('Recall: %s',mode))
% title('Recall')
% xlabel('Threshold (z-score)')
% ylim([0 1.05])
% 
% subplot(2,2,4)
% hold on
% plot(amp_sweep,f1_ampsweep_med(:,1),'Color',cols{1},'LineWidth',1)
% plot(amp_sweep,f1_ampsweep_med(:,2),'Color',cols{2},'LineWidth',1)
% plot(amp_sweep,f1_ampsweep_med(:,3),'Color',cols{3},'LineWidth',1)
% %title(sprintf('F1: %s', mode))
% title('F1')
% xlabel('Threshold (z-score)')
% ylim([0 1.05])
% 
% sgtitle('Time-domain notch filtered')
% 
% fprintf('\n')
% for kk = 1:3
%     [~, maxind] = max(f1_ampsweep_med(:,kk));
%     fprintf('%s: max F1 = %.3f, threshold = %.2f\n', labels{kk}, max(f1_ampsweep_med(:,kk)), amp_sweep(maxind));
% end
% fprintf('\n')