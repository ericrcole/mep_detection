%split saved file into GPi and STN-specific parts, and save

load('EP_102924_final_all.mat')

for kk = length(ep_struct):-1:1
    if ~contains(ep_struct(kk).pt_info.dbs_target, 'STN')
        ep_struct(kk) = [];
    end
end

save('EP_102924_final_STN.mat', 'ep_struct')

%% 

load('EP_102924_final_all.mat')

for kk = length(ep_struct):-1:1
    if ~contains(ep_struct(kk).pt_info.dbs_target, 'GPi')
        ep_struct(kk) = [];
    end
end

save('EP_102924_final_GPi.mat', 'ep_struct')

%% check GPI data

for kk = 1:length(ep_struct)
    if ep_struct(kk).pt_info.impute_emg_thresh > 1
        disp(ep_struct(kk).patient_ID)
    end
end

%% check that biomarker values are reasonable
pt_ind = 12;

fprintf('\n%s:\n', ep_struct(pt_ind).patient_ID)

pstrings = ep_struct(pt_ind).stim_settings.param_strings;
dlep_temp = ep_struct(pt_ind).dlep.rms_max;
mep_temp = ep_struct(pt_ind).mep.opt.mepscore_log;

for kk = 1:length(ep_struct(pt_ind).stim_settings.param_strings)
    fprintf('%s: DLEP = %.2f, MEP = %.2f\n', pstrings{kk}, dlep_temp(kk), mep_temp(kk))
end


%% check number of NANs in mep data

for kk = 1:length(ep_struct)
    mep_temp = ep_struct(kk).mep.opt.mepscore_log;
    fprintf('%s: %.2f\n', ep_struct(kk).patient_ID, mean(isnan(mep_temp)))
end