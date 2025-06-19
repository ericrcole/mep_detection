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
        [pkval,loc] =  findpeaks(abs(temptrace(latency_pts(1):latency_pts(2))),'MinPeakHeight',thresh,'MaxPeakWidth',fs*0.015,'MinPeakWidth',fs*0.0015,'MinPeakProminence',thresh);


        amps_temp = nan(size(pkval));
        delays_temp = nan(size(pkval));
        amps_z_temp = nan(size(pkval));

        if ~isempty(pkval)
            detected(kk) = 1;
            for kl = 1:length(pkval)
                amps_temp(kl) = pkval(kl);
                delays_temp(kl) = t(loc(kl))+latency(1);
                amps_z_temp(kl) = (amps_temp(kl))/(thresh/sd_thresh);
            end
        end
        amps{kk} = amps_temp;
        delays{kk} = delays_temp;
        amps_z{kk} = amps_z_temp;
    end
    detected = logical(detected);
end