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

        baseline = median((temptrace(end-round(base_length*fs/1000):end)));
        sd = std(temptrace(end-round(base_length*fs/1000):end));
        temptrace = (temptrace -baseline)/sd;

        z_traces(kl,:) = temptrace;
    end
    %dlep_traces = cellfun(@minus,ep_struct(pt_ind).dlep.dlep_means(kk,:),ep_struct(pt_ind).dlep.dlep_means_filt(kk,:),'Un',0);
    if ~(strcmp(art_type,'emg')) && ~isempty(dlep_means)
        dlep_inds = cellfun(@(x) length(x)>1,dlep_means,'UniformOutput',true);
        dlep_length  = max(cellfun(@length,dlep_means));
        z_traces_dlep = zeros(sum(dlep_inds),dlep_length);
    else
        dlep_inds = [];
    end

    counter = 1;
    if ~isempty(dlep_means)
        for kl = find(dlep_inds)
            temptrace = dlep_means{kl};
    
            baseline = median((temptrace(end-round(base_length*fs/1000):end)));
            sd = std(temptrace(end-round(base_length*fs/1000):end));
            temptrace = (temptrace -baseline)/sd;
    
            z_traces_dlep(counter,:) = temptrace;
    
            counter = counter + 1;
        end
    end
    %z_traces_dlep(~dlep_inds,:) = [];
    
    if strcmp(art_type,'dlep')
        %only use LFP data to estimate artifact with PCA
        [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(z_traces_dlep(:,t_art_samps_dlep(1):t_art_samps_dlep(2)),'Centered',false);
    elseif strcmp(art_type,'both')
        %stack EMG and LFP data to estimate artifact with PCA
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
        %EMG data only for PCA estimation
        if strcmp(mep_settings.art_chans, 'facial')
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(z_traces(is_facial,t_art_samps(1):t_art_samps(2)),'Centered',false);
        else
            mep_temp = [z_traces(mep_settings.is_facial,t_art_samps(1):t_art_samps(2)); z_traces(~mep_settings.is_facial,t_art_samps(1):t_art_samps(2))];
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(mep_temp,'Centered',false);
        end
    end
    if ~strcmp(art_type,'template') && ~strcmp(art_type,'template_pca')
        art_template = COEFF(:,1)';
    end
    %art_template = mean(z_traces(:,t_art_samps(1):t_art_samps(2)));
    %art_template = median(z_traces(:,t_art_samps(1):t_art_samps(2)));

    ch_counter = 1;
    for kl = 1:n_emg
        if strcmp(art_type,'template')
            %find most correlated template to fit and subtract from given channel
            corrmat = corrcoef([z_traces(kl,t_art_samps(1):t_art_samps(2)); mep_settings.templib(:,t_art_samps(1):t_art_samps(2))]');
            [~,maxind] = max(corrmat(1,2:end));
            art_template = mep_settings.templib(maxind,t_art_samps(1):t_art_samps(2));
        end

        if strcmp(art_type,'template_pca')
            %find top N most correlated pre-saved templates, stack them
            %with facial channel data when performing PCA
            k_emg = sum(mep_settings.is_facial);
            corrmat = corrcoef([z_traces(kl,t_art_samps(1):t_art_samps(2)); mep_settings.templib(:,t_art_samps(1):t_art_samps(2))]');
            [~,maxinds] = maxk(corrmat(1,2:end), k_emg);

            mep_temp = [z_traces(mep_settings.is_facial,t_art_samps(1):t_art_samps(2)); mep_settings.templib(maxinds,t_art_samps(1):t_art_samps(2))];
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(mep_temp,'Centered',false);
            art_template = COEFF(:,1)';
        end
        
        if (strcmp(mep_settings.art_chans,'facial') && mep_settings.is_facial(kl)) && strcmp(fit_type,'recon')
            ch_sel = ch_counter; %if only using some channels for PCA, make sure we use the right one
            ch_counter = ch_counter + 1;
        elseif (strcmp(mep_settings.art_chans,'facial') && ~mep_settings.is_facial(kl)) && strcmp(fit_type,'recon')
            ch_sel = 1; %limb channel, but only using facial so don't need to do anything
        else
            ch_sel = kl;
        end

        if strcmp(fit_type,'proj')
            template_fit = ortho_proj(z_traces(kl,t_art_samps(1):t_art_samps(2)),art_template);
        elseif strcmp(fit_type,'median')
            temp1 = z_traces(kl,t_art_samps(1):t_art_samps(2)) - median(z_traces(kl,t_art_samps(1):t_art_samps(2)));
            temp2 = art_template - median(art_template);
            template_fit = art_template*median(temp1./temp2);
        elseif strcmp(fit_type,'pctl')
            temp1 = z_traces(kl,t_art_samps(1):t_art_samps(2)) - median(z_traces(kl,t_art_samps(1):t_art_samps(2)));
            temp2 = art_template - median(art_template);
            template_fit = art_template*prctile(temp1./temp2, 20);
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
