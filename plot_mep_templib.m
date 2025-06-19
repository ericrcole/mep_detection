% plot_mep_templib.m
%
% plotting stim-averaged EMG traces to select examples containing only
% stimulation artifact, to build a library that can be leveraged for
% artifact rejection

%% load data
load('EP_042023.mat')
% 
% latency_labels = {'bicep','ECR','FCR','FDI','nasalis','orb oris','geniogl','trap'};
% 
% latencies = [16.92,24.56;16.59,29.24;18.31,28.69;24,26;8,13;8.17,14.21; 8,13; 14.71,14.71];
% latency_tolerance = 3;
% 
% latencies(:,1) = latencies(:,1) - latency_tolerance;
% latencies(:,2) = latencies(:,2) + latency_tolerance;

ltc = [5,16];

%% plotting block
pt_ind = 12;

use_thresh = false;
show_marks = true;
show_thresh = true; thresh_to_plot = 2.5;

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

fs = 22000;

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
mep_settings.art_type = 'emg';
mep_settings.art_chans = 'all';
mep_settings.fit_type = 'recon';

f_hp = 10; f_lp = 2000;
[b_high, a_high] = butter(2,[20]/22000,'high');
[b_low, a_low] = butter(2,[1000]/22000,'low');

figure
hold on
col_code = {'g','r','b','k'};
offset_y = 50; offset_x = 110; x_plot = linspace(1,100,2200); 
base_length = 220;

p_prominences = cell(size(ep_struct(pt_ind).mep.labels_detect));
p_widths = cell(size(ep_struct(pt_ind).mep.labels_detect)); 
p_locs = cell(size(ep_struct(pt_ind).mep.labels_detect));
p_amps_z = cell(size(ep_struct(pt_ind).mep.labels_detect)); 
pred_type = zeros(size(ep_struct(pt_ind).mep.labels_detect));

art_fits = cell(size(ep_struct(pt_ind).mep.labels_detect));

xtick_pos = (1:offset_x:offset_x*length(ep_struct(pt_ind).mep.emg_labels))+offset_x;
ytick_pos = (1:offset_y:offset_y*size(ep_struct(pt_ind).mep.labels_detect,1))+offset_y;

st_to_print = cell(length(ep_struct(pt_ind).mep.emg_labels),1);

for kp = 1:length(ep_struct(pt_ind).mep.emg_labels)
    detected_temp = zeros(size(ep_struct(pt_ind).mep.labels_detect,1),1);
    ltc_to_plot = [5,16];

    facial_log = zeros(size(ep_struct(pt_ind).mep.emg_labels));
    facial_log(5:end) = 1;

    mep_settings.is_facial = logical(facial_log);

    tr_counter = 1;

    for ko = 1:size(ep_struct(pt_ind).mep.labels_detect,1)
        temptrace = ep_struct(pt_ind).mep.mep_means{ko,kp};
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
            mep_means = ep_struct(pt_ind).mep.mep_means(ko,:)';
            %dlep_means =  ep_struct(pt_ind).dlep.raw_means(ko,:)';
            latencies_pca = repmat([5,16],[8 1]);
            [detected, amps, delays, amps_z, corr_traces, templates] = detect_mep_rejpca(mep_means, [], latencies_pca, mep_settings);
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

        if detected
            plot(x_plot+kp*offset_x,temptrace+ko*offset_y,'b')
        else
            plot(x_plot+kp*offset_x,temptrace+ko*offset_y,'k')
        end

        if rem_art_pca && art_plot_template
            templ_offset = median(temptrace(t_art_samps(1):t_art_samps(2)) - pc_template);
            plot(x_plot(t_art_samps(1):t_art_samps(2))+kp*offset_x, pc_template+templ_offset+ko*offset_y, 'm--')
        elseif rem_art_pca && art_plot_corrected

            plot(x_plot+kp*offset_x,corr_trace+ko*offset_y,col_code{4})
        end

        if ~isempty(delays_temp) && show_marks
            for kq = 1:length(delays_temp)
                delay_ind = find(x_plot>delays_temp(kq),1);
                plot(x_plot(delay_ind) +kp*offset_x, temptrace(delay_ind) + ko*offset_y, 'k*', 'MarkerSize', 20);
            end
        end
        if show_thresh
            plot([x_plot(1), x_plot(end)]+kp*offset_x,[-thresh_to_plot,-thresh_to_plot ]+ko*offset_y,'k--')
            plot([x_plot(1), x_plot(end)]+kp*offset_x,[thresh_to_plot,thresh_to_plot ]+ko*offset_y,'k--')
        end

        tr_counter = tr_counter + 1;
    end
    
    if show_thresh
        ltc_to_plot = [5 16];
%         plot([x_plot(1), x_plot(end)]+kp*offset_x,[mean(temptrace) - thresh_to_plot,mean(temptrace) - thresh_to_plot ]+ko*offset_y,'k--')
%         plot([x_plot(1), x_plot(end)]+kp*offset_x,[mean(temptrace) + thresh_to_plot,mean(temptrace) + thresh_to_plot ]+ko*offset_y,'k--')
        xline(kp*offset_x + ltc_to_plot(1), 'k--'); xline(kp*offset_x + ltc_to_plot(2), 'k--'); 
        %%
        %st_to_print{kp} = sprintf('Accuracy for %s: %.3f           Latency: %.3f - %.3f ms',ep_struct(pt_ind).mep.emg_labels{kp},mean(detected_temp == mep_labels_sm{pt_ind}(:,kp)), ltc_to_plot(1)-latency_tolerance, ltc_to_plot(2)+latency_tolerance);

    end
    
end

set(gca,'xTick',xtick_pos)
set(gca,'XtickLabel',ep_struct(pt_ind).mep.emg_labels)
ylim([- 150 2500])

set(gca,'yTick',ytick_pos)
set(gca,'YtickLabel',1:size(ep_struct(pt_ind).mep.labels_detect,1))

title(sprintf('patient %d: %s',pt_ind,ep_struct(pt_ind).patient_ID))


function [detected, amps, delays, amps_z, corr_traces, templates] = detect_mep_rejpca(mep_means, dlep_means, latencies, mep_settings)
    %detect meps in facial muscles. this version derives a template of the
    %artifact by performing PCA across recording channels
    %
    %example settings:
%     mep_settings = struct();
%     mep_settings.art_window = [1,20];
%     mep_settings.art_detect = [3,6]; %window to check whether to apply art correction
%     mep_settings.art_thresh = 3;
%     mep_settings.sd_thresh = 3;
%     mep_settings.fs = 22000;
%     mep_settings.base_length = 10;   %round(10*mep_settings.fs/1000);
%     mep_settings.pctl_thresh = 50;
%     mep_settings.use_dbs = true;
%     mep_settings.notch_filt_emg = true;
%     mep_settings.art_type = 'both';
    art_type = mep_settings.art_type;
    art_window = mep_settings.art_window;
    z_thresh = mep_settings.sd_thresh;
    base_length = mep_settings.base_length;
    fs = mep_settings.fs;
    fit_type = mep_settings.fit_type;
    is_facial = mep_settings.is_facial;
    
    ortho_proj = @(A,B) (sum(A.*B)/(norm(B)^2))*B;
    
    n_emg = length(mep_means);
    detected = nan(1,n_emg);
    amps = cell(1,n_emg);
    delays = cell(1,n_emg);
    amps_z = cell(1,n_emg);
    corr_traces = cell(1,n_emg);
    templates = cell(1,n_emg);

    t_art_samps = round(art_window*fs/1000);
    t_art_samps_dlep = t_art_samps;

    z_traces = nan(n_emg,length(mep_means{1}));
    for kl = 1:n_emg
        temptrace = mep_means{kl};

        baseline = median(abs(temptrace(end-base_length:end)));
        sd = std(temptrace(end-base_length:end));
        temptrace = (temptrace -baseline)/sd;

        z_traces(kl,:) = temptrace;
    end
    %dlep_traces = cellfun(@minus,ep_struct(pt_ind).dlep.dlep_means(kk,:),ep_struct(pt_ind).dlep.dlep_means_filt(kk,:),'Un',0);
    if ~(strcmp(art_type,'emg'))
        dlep_inds = cellfun(@(x) length(x)>1,dlep_means,'UniformOutput',true);
        dlep_length  = max(cellfun(@length,dlep_means));
        z_traces_dlep = zeros(sum(dlep_inds),dlep_length);
    else
        dlep_inds = [];
    end

    counter = 1;
    for kl = find(dlep_inds)
        temptrace = dlep_means{kl};

        baseline = median(abs(temptrace(end-base_length:end)));
        sd = std(temptrace(end-base_length:end));
        temptrace = (temptrace -baseline)/sd;

        z_traces_dlep(counter,:) = temptrace;

        counter = counter + 1;
    end
    %z_traces_dlep(~dlep_inds,:) = [];
    
    if strcmp(art_type,'dlep')
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(z_traces_dlep(:,t_art_samps_dlep(1):t_art_samps_dlep(2)),'Centered',false);
    elseif strcmp(art_type,'both')
        dlep_temp = z_traces_dlep(:,t_art_samps_dlep(1):t_art_samps_dlep(2));
        dlep_temp = (dlep_temp - median(dlep_temp(:)))./(max(dlep_temp(:))-min(dlep_temp(:)));
        if strcmp(mep_settings.art_chans, 'facial')
            mep_temp = z_traces(mep_settings.is_facial,t_art_samps(1):t_art_samps(2));
        else
            mep_temp = [z_traces(mep_settings.is_facial,t_art_samps(1):t_art_samps(2)); z_traces(~mep_settings.is_facial,t_art_samps(1):t_art_samps(2))];
        end
        %mep_temp = (mep_temp - median(mep_temp(:)))./(max(mep_temp(:))-min(mep_temp(:)));
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca([mep_temp; dlep_temp],'Centered',false);
    elseif strcmp(art_type,'emg')
        if strcmp(mep_settings.art_chans, 'facial')
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(z_traces(is_facial,t_art_samps(1):t_art_samps(2)),'Centered',false);
        else
            mep_temp = [z_traces(mep_settings.is_facial,t_art_samps(1):t_art_samps(2)); z_traces(~mep_settings.is_facial,t_art_samps(1):t_art_samps(2))];
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(mep_temp,'Centered',false);
        end
    end
    art_template = COEFF(:,1)';
    %art_template = mean(z_traces(:,t_art_samps(1):t_art_samps(2)));
    %art_template = median(z_traces(:,t_art_samps(1):t_art_samps(2)));

    ch_counter = 1;
    for kl = 1:n_emg
        %%% ADD: reconstruction mode
        if (strcmp(mep_settings.art_chans,'facial') && mep_settings.is_facial(kl)) && strcmp(fit_type,'recon')
            ch_sel = ch_counter; %if only using some channels for PCA, make sure we use the right one
            ch_counter = ch_counter + 1;
        else
            ch_sel = kl;
        end

        if strcmp(fit_type,'proj')
            template_fit = ortho_proj(z_traces(kl,t_art_samps(1):t_art_samps(2)),art_template);
        elseif strcmp(fit_type,'median')
            temp1 = z_traces(kl,t_art_samps(1):t_art_samps(2)) - median(z_traces(kl,t_art_samps(1):t_art_samps(2)));
            temp2 = art_template - median(art_template);
            template_fit = art_template*median(temp1./temp2);
        elseif strcmp(fit_type,'recon')
            temp1 = SCORE(ch_sel,:); temp1(2:end) = 0;
            template_fit = temp1*COEFF';
        end


        if strcmp(fit_type,'recon')
            temp1 = SCORE(ch_sel,:); temp1(1) = 0;
            temptrace = z_traces(kl,:);
            temptrace(t_art_samps(1):t_art_samps(2)) = temp1*COEFF';
        else
            temptrace = z_traces(kl,:);
            temptrace(t_art_samps(1):t_art_samps(2)) = temptrace(t_art_samps(1):t_art_samps(2)) - template_fit;
        end
%         else
%             temptrace = z_traces(kl,:); 
%             template_fit = zeros(size(temptrace));
%         end
        if (strcmp(mep_settings.art_chans,'facial') && ~mep_settings.is_facial(kl))
            temptrace = z_traces(kl,:);
            template_fit = zeros(size(temptrace(t_art_samps(1):t_art_samps(2))));
        end
        
        [detected(kl), amps{kl}, delays{kl}, amps_z{kl}] = ep_thresh_windowed(temptrace,base_length,z_thresh,fs,latencies(kl,:));
        corr_traces{kl} = temptrace;
        templates{kl} = template_fit;
        
    end

end

function [detected, amps, delays, amps_z] = ep_thresh_windowed(ep_means,window_base,sd_thresh,fs,latency)
    detected = zeros(size(ep_means,1),1);
    amps = cell(size(detected));
    delays = cell(size(detected));
    amps_z = cell(size(detected));

    t = linspace(0,100,size(ep_means,2));
    base_length = window_base*fs/1000;
    latency_pts = round(latency/1000*fs);
    for kk = 1:size(ep_means,1)
%         baseline = median(abs(ep_means(kk,end-base_length:end)));
%         thresh = sd_thresh*std(ep_means(kk,end-base_length:end));
%         [pkval,loc] =  findpeaks(abs(ep_means(kk,latency_pts(1):latency_pts(2))),'MinPeakProminence',baseline+thresh,'MaxPeakWidth',fs*0.010,'MinPeakWidth',fs*0.0025);
%         
        temptrace = ep_means(kk,:);
        baseline = median(ep_means(kk,end-base_length:end));
        thresh = sd_thresh*std(ep_means(kk,end-base_length:end));
        temptrace = temptrace - baseline;
            %old settings from 1/17
        %[pkval,loc] =  findpeaks(abs(temptrace(latency_pts(1):latency_pts(2))),'MinPeakProminence',thresh,'MaxPeakWidth',fs*0.010,'MinPeakWidth',fs*0.0025);
        
            %medium version
        %[pkval,loc] =  findpeaks(abs(temptrace(latency_pts(1):latency_pts(2))),'MinPeakHeight',thresh,'MaxPeakWidth',fs*0.010,'MinPeakWidth',fs*0.001);
        
            %after tuning:
        [pkval,loc] =  findpeaks(abs(temptrace(latency_pts(1):latency_pts(2))),'MinPeakHeight',2.5,'MaxPeakWidth',fs*0.015,'MinPeakWidth',fs*0.0015,'MinPeakProminence',2.5);


        amps_temp = nan(size(pkval));
        delays_temp = nan(size(pkval));
        amps_z_temp = nan(size(pkval));

        if ~isempty(pkval)
            detected(kk) = 1;
            for kl = 1:length(pkval)
                amps_temp(kl) = pkval(kl);
                delays_temp(kl) = t(loc(kl))+latency(1);
                amps_z_temp(kl) = (amps_temp(kl)-baseline)/(thresh/sd_thresh);
            end
        end
        amps{kk} = amps_temp;
        delays{kk} = delays_temp;
        amps_z{kk} = amps_z_temp;
    end
    detected = logical(detected);
end

