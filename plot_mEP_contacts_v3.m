%plot_EP_contacts_v3.m
%
%visualize DLEP/mEP amplitudes, biomarkers, etc from processed data
%vs. stimulation parameter space for given patient

plot_predicted = false;
log_transform = false;
montage = 'monopolar_segmented';  
    %'monopolar' = only 4 ring contacts
    %'monopolar_segmented' = 8 contacts
    
lead = 'None';  %keep as 'None', or set to 'ABB' to rotate abbott leads

patient_ID = 'ephys030';
[~,struct_ind]=find(ismember({ep_struct(:).patient_ID},patient_ID),1);

pt_struct = ep_struct(struct_ind);

% exc_freq = contains(pt_struct.stim_settings.param_strings,'130Hz');

if plot_predicted
    mep_vals = transpose(data.mep{struct_ind});
    dlep_vals = transpose(data.dlep{struct_ind});
else
    mep_vals = pt_struct.mep.amps_avg;
    dlep_vals = pt_struct.dlep.rms_max;
end


mep_raw = mep_vals;
dlep_raw = dlep_vals;

if log_transform
   mep_vals = log10(mep_vals+0.01);
   dlep_vals = log10(dlep_vals+0.01);
end

% mep_vals(exc_freq) = []; 
% pt_struct.stim_settings.param_strings(exc_freq) = [];
% pt_struct.stim_settings.amplitudes(exc_freq) = [];
% pt_struct.stim_settings.cathodes(exc_freq) = [];
% pt_struct.stim_settings.anodes(exc_freq) = [];
% pt_struct.stim_settings.contacts(exc_freq) = [];

norm_gs = @(x) (x-min(x))/(max(x)-min(x));

mep_vals = norm_gs((mep_vals));
dlep_vals = norm_gs((dlep_vals));

ep_trials = 1:length(mep_vals);

amps = unique(pt_struct.stim_settings.amplitudes);

%trial_inds = parse_params('monopolar_segmented', pt_struct.stim_settings,ep_trials, 'all');
%trial_inds = parse_params('monopolar_segmented', pt_struct.stim_settings,ep_trials, amps, true, lead);

trial_inds = parse_params(montage, pt_struct.stim_settings,ep_trials, amps, true, lead);

figure
amps = unique(pt_struct.stim_settings.amplitudes);

grid_cols = length(amps);

for kk = 1:length(amps)
    subplot(2,grid_cols,kk)
    trials_amp = trial_inds(:,kk);
    plot_contacts(dlep_vals,trials_amp,montage)
    if kk == 1
        ylabel('DLEP')
    end
    title(sprintf('%d mA',amps(kk)))
end
for kk = 1:length(amps)
    subplot(2,grid_cols,kk+grid_cols)
    trials_amp = trial_inds(:,kk);
    plot_contacts(mep_vals,trials_amp,montage)
    if kk == 1
        ylabel('mEP')
    end
    title(sprintf('%d mA',amps(kk)))
end
suptitle(sprintf('%s; %s; %s',patient_ID,pt_struct.pt_info.lead_model_name,pt_struct.pt_info.lfp_target))

display(max(mep_raw));

figure
subplot(1,2,1)
histogram(dlep_raw)
title('DLEP vals')
subplot(1,2,2)
histogram(mep_raw)
title('mEP vals')

