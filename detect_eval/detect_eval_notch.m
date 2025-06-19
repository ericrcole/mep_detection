%% quantify effectiveness of different artifact rejection strategies:
% this version: evaluates frequency and time-domain based notch filtering

categs = {'std','freq','td'};
categ_names = {'Standard','Butterworth','Time-domain'};

addpath('C:\Users\ORStream\Documents\mEP_processing\functions')

tt = linspace(0,0.1,2200);

rem_art_pca = false;
use_tdnotch = false;
use_frnotch = false;
filt_sg = false;
rem_ringing_art = false;
filt_hp = false;
filt_lp = false;

mode = 'raw';

mep_settings = struct();
mep_settings.art_window = [1,20];
mep_settings.art_detect = [3,6]; %window to check whether to apply art correction
mep_settings.art_thresh = 3;
mep_settings.sd_thresh = 4;
mep_settings.fs = 22000;
mep_settings.base_length = 10;   %round(10*mep_settings.fs/1000);
mep_settings.pctl_thresh = 50;
mep_settings.use_dbs = true;
mep_settings.notch_filt_emg = false;
mep_settings.art_type = 'both';   %whether to fit artifact with EMG, LFP, or both
mep_settings.art_chans = 'facial'; %which chans to apply correction to
mep_settings.fit_type = 'median'; %how to fit/remove artifact template for individual channels


acc_ampsweep = cell(length(categs),3);
prec_ampsweep = cell(length(categs),3);
rec_ampsweep = cell(length(categs),3);
f1_ampsweep = cell(length(categs),3);

for ampind = 1:length(categs)

    acc_by_muscle = cell(size(emg_mapping_labels));
    prec_by_muscle = cell(size(emg_mapping_labels));
    rec_by_muscle = cell(size(emg_mapping_labels));
    f1_by_muscle = cell(size(emg_mapping_labels));

    if strcmp(categs{ampind},'std')
        rem_art_pca = false;
        filt_hp = false;
    elseif strcmp(categs{ampind},'freq')
%         rem_art_pca = true;
%         mep_settings.art_type = 'emg';  
%         mep_settings.art_chans = 'facial'; 
        use_frnotch = true; use_tdnotch = false;
        [b_notch, a_notch] = butter(2, [59 61]./(fs/2), 'stop');
    elseif strcmp(categs{ampind},'td')
        use_tdnotch = true;
    end

    for pt_ind = 1:length(ep_struct)
        [facial_log,facial_inds] = ismember(emg_mapping_all{pt_ind},is_facial);
        mep_settings.is_facial = logical(facial_log);
        for kp = 1:length(emg_mapping_all{pt_ind})
            n_trials = length(ep_struct(pt_ind).stim_settings.param_strings);
            detected_temp = zeros(n_trials,1);
            ltc_to_plot = latencies(emg_mapping_all{pt_ind}(kp),:) + [-latency_tolerance, latency_tolerance];

            for ko = 1:n_trials
                temptrace = double(ep_struct(pt_ind).mep.(mode).mep_means{ko,kp});
                if isempty(temptrace)
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
                    temptrace(40:end) = filtfilt(b_notch, a_notch, temptrace(40:end));
                end

                if use_tdnotch
                    %temptrace = tdnotch(temptrace);
                    [temptrace, yest, A, phi] = td_notch(temptrace,tt,60,fs);
                end

                if rem_art_pca
                    mep_means = ep_struct(pt_ind).mep.raw.mep_means(ko,:)';
                    dlep_means =  ep_struct(pt_ind).dlep.raw_means(ko,:)';
                    latencies_pca = latencies(emg_mapping_all{pt_ind},:) + [-latency_tolerance, latency_tolerance];
                    [detected, amps, delays, amps_z, corr_traces, templates] = detect_mep_rejpca(mep_means, dlep_means, latencies_pca, mep_settings);
                    pc_template = templates{kp};
                    corr_trace = corr_traces{kp};
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
                [detected, amps, delays, amps_z] = ep_thresh_windowed(temptrace,10,mep_settings.sd_thresh,fs,ltc_to_plot);
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
                detected_temp(ko) = detected;
                end

                %             if detected && (mep_labels_sm{pt_ind}(ko,kp))
                %                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{1})
                %                 pred_type(ko) = 1;
                %             elseif detected && ~(mep_labels_sm{pt_ind}(ko,kp)) %false positive; red
                %                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{2})
                %                 pred_type(ko) = 2;
                %             elseif ~detected && (mep_labels_sm{pt_ind}(ko,kp)) %false negative; blue
                %                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{3})
                %                 pred_type(ko) = 3;
                %             else
                %                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{4})
                %                 pred_type(ko) = 4;
                %             end
                %             if ~isempty(delays_temp) && show_marks
                %                 for kq = 1:length(delays_temp)
                %                     delay_ind = find(x_plot>delays_temp(kq),1);
                %                     plot(x_plot(delay_ind) +kp*offset_x, temptrace(delay_ind) + ko*offset_y, 'k*', 'MarkerSize', 20);
                %                 end
                %             end
                %             if show_thresh
                %                 plot([x_plot(1), x_plot(end)]+kp*offset_x,[-thresh_to_plot,-thresh_to_plot ]+ko*offset_y,'k--')
                %                 plot([x_plot(1), x_plot(end)]+kp*offset_x,[thresh_to_plot,thresh_to_plot ]+ko*offset_y,'k--')
                %             end

            end
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
            %             if isnan(f1) && (mean(mep_labels_sm{pt_ind}(:,kp)) == 0)
            %                 f1 = acctemp;
            %             elseif isnan(f1) && ((precision == 0) && (recall == 0))
            %                 f1 = 0;
            %             end
            if isnan(f1) && ((precision == 0) && (recall == 0))
                f1 = 0;
            end


            f1_by_muscle{emg_mapping_all{pt_ind}(kp)} = [f1_by_muscle{emg_mapping_all{pt_ind}(kp)}; f1];
            %         if show_thresh
            %             ltc_to_plot = latencies(emg_mapping_all{pt_ind}(kp),:);
            %             %         plot([x_plot(1), x_plot(end)]+kp*offset_x,[mean(temptrace) - thresh_to_plot,mean(temptrace) - thresh_to_plot ]+ko*offset_y,'k--')
            %             %         plot([x_plot(1), x_plot(end)]+kp*offset_x,[mean(temptrace) + thresh_to_plot,mean(temptrace) + thresh_to_plot ]+ko*offset_y,'k--')
            %             xline(kp*offset_x + ltc_to_plot(1) - latency_tolerance, 'k--'); xline(kp*offset_x + ltc_to_plot(2) + latency_tolerance, 'k--');
            %             st_to_print{kp} = sprintf('Accuracy for %s: %.3f           Latency: %.3f - %.3f ms',ep_struct(pt_ind).mep.emg_labels{kp},mean(detected_temp == mep_labels_sm{pt_ind}(:,kp)), ltc_to_plot(1)-latency_tolerance, ltc_to_plot(2)+latency_tolerance);
            %         else
            %             fprintf('Accuracy for %s: %.3f\n',ep_struct(pt_ind).mep.emg_labels{kp},acc_by_pt{pt_ind}(kp))
            %         end

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

%figure('Position',[100 100 1000 600])
subplot(1,2,2)
for kk = 1:length(categs)
    temp = acc_ampsweep{kk,3};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('Accuracy')
xlabel('Filter type')

subplot(1,2,2)
for kk = 1:length(categs)
    temp = f1_ampsweep{kk,3};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp,20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('F1')
xlabel('Filter type')

%sgtitle('Notch filtering')

%%
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

%% 
figure
subplot(1,2,1)
for kk = 1:length(categs)
    temp = acc_ampsweep{kk,2};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('Accuracy')
xlabel('Filter type')
title('Facial')

subplot(1,2,2)
for kk = 1:length(categs)
    temp = acc_ampsweep{kk,3};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('Accuracy')
xlabel('Filter type')
title('Limb')


%%
figure('Position',[100 100 1000 600])
subplot(2,2,1)
for kk = 1:length(categs)
    temp = acc_ampsweep{kk,2};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('Accuracy')
xlabel('Filter type')
title('Facial')

subplot(2,2,2)
for kk = 1:length(categs)
    temp = f1_ampsweep{kk,2};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp,20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('F1')
xlabel('Filter type')
title('Facial')

subplot(2,2,3)
for kk = 1:length(categs)
    temp = acc_ampsweep{kk,3};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('Accuracy')
xlabel('Filter type')
title('Limb')

subplot(2,2,4)
for kk = 1:length(categs)
    temp = f1_ampsweep{kk,3};
    boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    hold on
    scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp,20, [0.5 0.5 0.5],'filled')
end
xticks(1:length(categs))
xticklabels(categ_names)
ylabel('F1')
xlabel('Filter type')
title('Limb')

sgtitle('Notch filtering')

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