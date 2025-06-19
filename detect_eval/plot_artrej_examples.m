%% test mEP artifact rejection on individual examples
pt_ind =6;
ko = 26;

offset_y = 100; offset_x = 110; x_plot = linspace(1,100,2200); 
%base_length = 220;
sd_thresh = 8;

use_thresh = false;
show_marks = true;
show_thresh = true; thresh_to_plot = 4;

base_length = 10;

is_facial = [5,6,7,8,11,12,13];

show_original = false; 

rem_ringing_art = false;
rem_art_mscale = false;

rem_art_pca = true;
art_plot_template = true;
art_plot_corrected = false;

filt_lp = false;
filt_hp = false;
filt_sg = false;

mep_settings = struct();
mep_settings.art_window = [2,20];
mep_settings.art_detect = [3,6]; %window to check whether to apply art correction
mep_settings.art_thresh = 3; 
mep_settings.sd_thresh = 3;
mep_settings.fs = 22000;
mep_settings.base_length = 10;   %round(10*mep_settings.fs/1000);
mep_settings.pctl_thresh = 50;
mep_settings.use_dbs = true;
mep_settings.notch_filt_emg = true;
mep_settings.art_type = 'template';   %whether to fit artifact with PCA on EMG, LFP, or both chanels, or use 'template'
mep_settings.art_chans = 'facial'; %which chans to apply correction to
mep_settings.fit_type = 'median'; %how to fit/remove artifact template for individual channels

plot_corrected = true;

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

tr_counter = 1;
figure
hold on

emg_loop = [1 2 3 4 8 6];
for ind = 1:length(emg_loop)
    kp = emg_loop(ind);
    detected_temp = zeros(size(ep_struct(pt_ind).mep.raw.labels_detect,1),1);
    ltc_to_plot = latencies(emg_mapping_all{pt_ind}(kp),:) + [-latency_tolerance, latency_tolerance];
    [facial_log,facial_inds] = ismember(emg_mapping_all{pt_ind},is_facial);

    mep_settings.is_facial = logical(facial_log);

    temptrace = ep_struct(pt_ind).mep.raw.mep_means{ko,kp};
    temptrace_og = temptrace;

    if isempty(temptrace)
        continue
    end
    if filt_sg
        freq=150;
        order=2;

        framelen=round(1/freq*fs);
        if mod(framelen,2)==0
            framelen=framelen+1;
        end
        temptrace = sgolayfilt(temptrace,order,15);

    end
    if rem_ringing_art
        [temptrace] = ep_remove_ringing_artifact_erc(temptrace, [.005 .015], fs);
    end

    if filt_hp
        temptrace = filtfilt(b_high, a_high, temptrace);
    end
    if filt_lp
        temptrace = filtfilt(b_low, a_low, temptrace);
    end

    if (kp == 5) && (ko == 5)
        qqq = 2;
    end

    if rem_art_pca
        mep_means = ep_struct(pt_ind).mep.raw.mep_means(ko,:)';
        dlep_means =  ep_struct(pt_ind).dlep.raw_means(ko,:)';
        latencies_pca = latencies(emg_mapping_all{pt_ind},:) + [-latency_tolerance, latency_tolerance];
        [detected, amps, delays, amps_z, corr_traces, templates] = detect_mep_rejpca(mep_means, dlep_means, latencies_pca, mep_settings);
        pc_template = templates{kp};
        corr_trace = corr_traces{kp};
        t_art_samps = round(mep_settings.art_window*fs/1000);
        detected = detected(kp);

        baseline = median((temptrace(end-base_length:end)));
        sd = std(temptrace(end-base_length:end));
        temptrace = (temptrace -baseline)/sd;
        latency_pts = round(ltc_to_plot/1000*fs);
    else
        baseline = median((temptrace(end-base_length:end)));
        sd = std(temptrace(end-base_length:end));
        temptrace = (temptrace -baseline)/sd;
        latency_pts = round(ltc_to_plot/1000*fs);
    end

    if rem_art_mscale
        scales = 0.75:.025:1.25;
        [temptrace, template_fit, time, scale_min, dist_min] = remove_mep_artifact_mscale(temptrace, template_smoothed, scales, 2, 22000);

    end

    latency_pts = round(ltc_to_plot*fs/1000);
    [pkval, loc, widths, prom] = findpeaks(abs(corr_traces{kp}(latency_pts(1):latency_pts(2))),'MinPeakHeight',2.5,'MaxPeakWidth',fs*0.015,'MinPeakWidth',fs*0.0015,'MinPeakProminence',2.5);

    %         delay_temp = ep_struct(pt_ind).mep.delays_detect{ko,kp};
    %         amps_temp = ep_struct(pt_ind).mep.amps_detect_z{ko,kp};
    amps_temp = transpose(nan(size(pkval)));
    delays_temp = transpose(nan(size(pkval)));
    amps_z_temp = transpose(nan(size(pkval)));
    if ~rem_art_pca
        if ~isempty(pkval)
            detected = true;
            p_widths{ko,kp} = transpose(widths)/fs*1000;
            p_prominences{ko,kp} = transpose(prom);
            for kl = 1:length(pkval)
                amps_temp(kl) = pkval(kl);
                delays_temp(kl) = x_plot(loc(kl))+ltc_to_plot(1);
                amps_z_temp(kl) =  amps_temp(kl);
            end
            p_locs{ko,kp} = delays_temp;
            p_amps_z{ko,kp} = amps_z_temp;
        else
            detected = false;
        end
    end
    detected_temp(ko) = detected;

    if show_original
        temptrace = temptrace_og;
    end

%     if detected && (mep_labels_sm{pt_ind}(ko,kp))
%         plot(x_plot,temptrace+kp*offset_y,col_code{1})
%         pred_type(ko,kp) = 1;
%     elseif detected && ~(mep_labels_sm{pt_ind}(ko,kp)) %false positive; red
%         plot(x_plot,temptrace+kp*offset_y,col_code{2})
%         pred_type(ko,kp) = 2;
%     elseif ~detected && (mep_labels_sm{pt_ind}(ko,kp)) %false negative; blue
%         plot(x_plot,temptrace+kp*offset_y,col_code{3})
%         pred_type(ko,kp) = 3;
%     else
%         plot(x_plot,temptrace+kp*offset_y,col_code{4})
%         pred_type(ko,kp) = 4;
%     end
    if ~plot_corrected
        plot(x_plot,temptrace+ind*offset_y,'k')

        if rem_art_pca && art_plot_template
            templ_offset = median(temptrace(t_art_samps(1):t_art_samps(2)) - pc_template);
            plot(x_plot(t_art_samps(1):t_art_samps(2)), pc_template+templ_offset+ind*offset_y, 'm--')
        elseif rem_art_pca && art_plot_corrected

            plot(x_plot,corr_trace+ind*offset_y,col_code{4})
        end
    else
        %templ_offset = median(temptrace(t_art_samps(1):t_art_samps(2)) - pc_template);
        temp_to_plot = temptrace;
        temp_to_plot(t_art_samps(1):t_art_samps(2)) = temp_to_plot(t_art_samps(1):t_art_samps(2)) - pc_template;
        plot(x_plot,temp_to_plot+ind*offset_y,'b')
    end

    

    if ~isempty(delays_temp) && show_marks
        for kq = 1:length(delays_temp)
            delay_ind = find(x_plot>delays_temp(kq),1);
            plot(x_plot(delay_ind), temptrace(delay_ind) + ind*offset_y, 'k*', 'MarkerSize', 20);
        end
    end
    if show_thresh
        plot([x_plot(1), x_plot(end)],[-thresh_to_plot,-thresh_to_plot ]+ind*offset_y,'k--')
        plot([x_plot(1), x_plot(end)],[thresh_to_plot,thresh_to_plot ]+ind*offset_y,'k--')
    end

    tr_counter = tr_counter + 1;

end
yticks((1:length(emg_loop))*offset_y)
yticklabels(ep_struct(pt_ind).mep.emg_labels(emg_loop))
ylim([-100,900])
