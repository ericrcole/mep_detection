% %plot_mEP_characteristics.m
% %
% %

%plot_mEP_characteristics.m
%
%
load('M:/cEP_mEP_analysis/EPdata_5-17-2020_ANALYZED_postopMRI.mat')
    %saved data structure of svjetlana's labels for intraoperative EMG data
load('EP_mEP_011722.mat')
    %my processed data

%addpath(genpath('D:\ERC\EP_MetaOpt\Violinplot-Matlab'));

%%grab info necessary to match trials and muscles between my processed data
%%and svjetlana's manually labeled data
EP_reordered = readtable('EP_reordered.xlsx','Format','auto');

emg_mapping_pts = {'ephys010','ephys011','ephys012','ephys014','ephys015',...
        'ephys017','ephys018','ephys025','ephys026'}; %manual channel assignment because names don't match
emg_mapping_sm = {[1:8], [1:8], [1:8], [1:8], [1:8], [1:8],[1 2 3 4 6 5 8 7], [1 2 4 5 6 7], [1, 2, 4, 5]};

emg_mapping_all = {[1 2 3 4 10 5 7 9], [1 2 3 4 10 5 11 8], [1 2 3 4 10 5 7 11], [1 2 3 4 7 5 12 11], [1 2 3 4 7 5 12 11], [1 2 3 4 5 8 11 13],...
            [1 2 3 4 5 7 11 12], [1 2 4 5 8 11], [1 2 4 6]};

emg_mapping_labels = {'Biceps', 'ECR', 'FCR', 'FDI', 'Nasalis', 'Hypoglossus', 'Genioglossus',...
    'Orb oris','Trap', 'Tibia','Ipsi nasalis','Ipsi genioglossus', 'Ipsi orb oris'};
%labels for muscles

%% goes through SVJ's data and organizes them into structure for analysis
mep_labels_sm = cell(0,1);
mep_amps_sm = cell(0,1);
mep_ampsp2p_sm = cell(0,1); 
mep_delays_sm = cell(0,1); %one delay value of biggest peak per mEP
mep_delays_sm_all = cell(0,1); %all delay values

amplitudes = cell(0,1);
cathodes = cell(0,1);
anodes = cell(0,1);

ep_struct_inds_sm = [];
ep_struct_inds_co = [];

reordered_inds = cell(0,1);

for kk = 1:length(ep_struct)
    patient_ID = ep_struct(kk).patient_ID;
    
    amplitudes = [amplitudes; {ep_struct(kk).stim_settings.amplitudes}];
    cathodes = [cathodes; {ep_struct(kk).stim_settings.cathodes}];
    anodes = [anodes; {ep_struct(kk).stim_settings.cathodes}];

    %mep_reordered = EP_reordered.(patient_ID);

    pts_labeled_sm = {EP(:).code};
    if ~any(strcmp(pts_labeled_sm,patient_ID))
        continue
    end

    ep_struct_inds_sm = [ep_struct_inds_sm; kk];

    sm_ind = find(strcmp({EP.code},patient_ID));
    amps_full_sm = EP(sm_ind).amplitudesEMG;
    delays_full_sm = EP(sm_ind).latenciesEMG;

    labels_sm_temp = zeros(size(EP(sm_ind).latenciesEMG,2),length(emg_mapping_sm{kk}));
    amps_sm_temp = nan(size(EP(sm_ind).latenciesEMG,2),length(ep_struct(kk).mep.emg_labels));
    ampsp2p_sm_temp = nan(size(EP(sm_ind).latenciesEMG,2),length(ep_struct(kk).mep.emg_labels));
    delays_sm_temp = nan(size(EP(sm_ind).latenciesEMG,2),length(ep_struct(kk).mep.emg_labels));
    delays_sm_temp_all = cell(size(EP(sm_ind).latenciesEMG,2),length(ep_struct(kk).mep.emg_labels));

    for kl = 1:length(emg_mapping_sm{kk})
        emg_ind = kl;
        amps_full_temp = amps_full_sm(:,:,emg_ind);
        delays_full_temp = delays_full_sm(:,:,emg_ind);

        for km = 1:size(amps_full_temp,2)
            if all(isnan(amps_full_temp(:,km)))
                labels_sm_temp(km,kl) = 1;
            else
                [~,amp_ind] = max(abs(amps_full_temp(:,km)-amps_full_temp(1,km)));
                amps_sm_temp(km,kl) = amps_full_temp(amp_ind,km);
                delays_sm_temp(km,kl) = delays_full_temp(amp_ind,km);
                temp = ~isnan(delays_full_temp(:,km)); temp(1) = false;
                delays_sm_temp_all{km,kl} = delays_full_temp(temp,km);

                p2p_temp = diff(amps_full_temp(:,kl)-amps_full_temp(1,km));
                [~,p2p_ind] = max(p2p_temp);
                ampsp2p_sm_temp (km,kl) = p2p_temp(p2p_ind);
            end
        end
        
    end
    if strcmp(patient_ID, 'ephys025')
        labels_sm_temp(55:56,:) = [];
        amps_sm_temp(55:56,:) = [];
        ampsp2p_sm_temp(55:56,:) = [];
        delays_sm_temp(55:56,:) = [];
        delays_sm_temp_all(55:56,:) = [];
        if kk == 1
        amplitudes{kk}(55:56) = [];
        cathodes{kk}(55:56) = [];
        anodes{kk}(55:56) = [];
        end
    end
    mep_labels_sm = [mep_labels_sm; {~logical(labels_sm_temp)}];
    mep_amps_sm = [mep_amps_sm; {abs(amps_sm_temp)}];
    mep_ampsp2p_sm = [mep_ampsp2p_sm; {abs(ampsp2p_sm_temp)}];
    mep_delays_sm = [mep_delays_sm; {delays_sm_temp}];
    mep_delays_sm_all = [mep_delays_sm_all; {delays_sm_temp_all}];
end

% for kk = 1:length(ep_struct)
%     stim_reordered = EP_reordered.( ep_struct(kk).patient_ID);
%     stim_reordered = stim_reordered(~isnan(stim_reordered));
%     ep_struct(kk).mep.labels_detect = ep_struct(kk).mep.labels_detect(stim_reordered,:);
%     ep_struct(kk).mep.delays_detect = ep_struct(kk).mep.delays_detect(stim_reordered,:);
%     ep_struct(kk).mep.amps_detect = ep_struct(kk).mep.amps_detect(stim_reordered,:);
%     ep_struct(kk).stim_settings.param_strings = ep_struct(kk).stim_settings.param_strings(stim_reordered);
%     ep_struct(kk).mep.mep_means = ep_struct(kk).mep.mep_means(stim_reordered,:);
%     ep_struct(kk).mep.amps_detect_z = ep_struct(kk).mep.amps_detect_z(stim_reordered,:);
% end

for kk = 1:length(ep_struct)
    stim_reordered = EP_reordered.( ep_struct(kk).patient_ID);
    stim_reordered = stim_reordered(~isnan(stim_reordered));
    emg_reordered = emg_mapping_sm{kk};
    
    ep_struct(kk).mep.emg_labels = ep_struct(kk).mep.emg_labels(emg_reordered);
    ep_struct(kk).mep.labels_detect = ep_struct(kk).mep.labels_detect(stim_reordered,emg_reordered);
    ep_struct(kk).mep.delays_detect = ep_struct(kk).mep.delays_detect(stim_reordered,emg_reordered);
    ep_struct(kk).mep.amps_detect = ep_struct(kk).mep.amps_detect(stim_reordered,emg_reordered);
    ep_struct(kk).stim_settings.param_strings = ep_struct(kk).stim_settings.param_strings(stim_reordered);
    ep_struct(kk).mep.mep_means = ep_struct(kk).mep.mep_means(stim_reordered,emg_reordered);
    ep_struct(kk).mep.amps_detect_z = ep_struct(kk).mep.amps_detect_z(stim_reordered,emg_reordered);
end

%% start plotting the data
amps = [1 3 5];

mep_delays_by_muscle = cell(size(emg_mapping_labels));
stim_total_by_muscle = cell(size(emg_mapping_labels));
stim_mep_by_muscle = cell(size(emg_mapping_labels));
stim_trials_by_pt = cell(size(emg_mapping_sm));

for kl = 1:length(emg_mapping_labels)
    stim_mep_by_muscle{kl} = zeros(8,3);
    stim_total_by_muscle{kl} = zeros(8,3);
end

param_structs = {ep_struct(:).stim_settings};

pts_loop = [1 3 4 5 6 7 8 9];

%for kk = 1:length(mep_delays_sm_all)
for kk = pts_loop
    mono_trials =  parse_params('monopolar_segmented', ep_struct(kk).stim_settings,[1:size(mep_delays_sm{kk},1)], amps, true, ep_struct(kk).pt_info.lead_model_name);
    stim_trials_by_pt{kk} = mono_trials;
    for kl = 1:length(emg_mapping_sm{kk})
        mono_mep_labels = parse_params('monopolar_segmented', ep_struct(kk).stim_settings,mep_labels_sm{kk}(:,kl), amps, true, ep_struct(kk).pt_info.lead_model_name);
    
        mep_delays_by_muscle{emg_mapping_all{kk}(kl)} = [mep_delays_by_muscle{emg_mapping_all{kk}(kl)}; cell2mat(transpose({mep_delays_sm_all{kk}{:,kl}}))];
        
        mono_mep_labels(isnan(mono_mep_labels)) = 0; 
        stim_total_by_muscle{emg_mapping_all{kk}(kl)} = stim_total_by_muscle{emg_mapping_all{kk}(kl)} + (mono_trials > 0);
        stim_mep_by_muscle{emg_mapping_all{kk}(kl)} = stim_mep_by_muscle{emg_mapping_all{kk}(kl)} + mono_mep_labels;

    end
end

mep_delays_vplot = padcat(mep_delays_by_muscle{1},mep_delays_by_muscle{2},mep_delays_by_muscle{3},mep_delays_by_muscle{4}, ...
                    mep_delays_by_muscle{5},mep_delays_by_muscle{6},mep_delays_by_muscle{7},mep_delays_by_muscle{8},...
                    mep_delays_by_muscle{11},mep_delays_by_muscle{12});

muscle_loop = [1 2 3 4 5 6 7 8 11 12]; %non-empty muscles

% for kk = 2:length(mep_delays_by_muscle)
%     mep_delays_vplot = padcat(mep_delays_vplot,mep_delays_by_muscle{kk});
% end
for kl = 1:length(stim_mep_by_muscle)
    stim_mep_by_muscle{kl} = stim_mep_by_muscle{kl} ./ stim_total_by_muscle{kl};
end
%violin plot of mEP delays for diff muscles
figure
violinplot(mep_delays_vplot,emg_mapping_labels([1,2,3,4,5,6,7,8,11,12]),'Bandwidth',0.5)
xlabel('Muscle')
ylabel('Delay (ms)')

for kl = muscle_loop
    figure
    for kk = 1:3
    subplot(1,3,kk)
    plot_contacts(stim_mep_by_muscle{kl},'monopolar_segmented');
    if kk == 1
        ylabel('mEP')
    end
    title(sprintf('%d mA',amps(kk)))
    end
    suptitle(emg_mapping_labels{kl})
end
%suptitle(sprintf('%s; %s; %s',patient_ID,pt_struct.pt_info.lead_model_name,pt_struct.pt_info.lfp_target))


%trial_inds = parse_params('monopolar_segmented', ep_struct(kk).stim_settings,[1:size(mep_delays_sm{kk},1)], 'all', true, ep_struct(kk).pt_info.lead_model_name);

% mep_delays_by_muscle = cell(size(emg_mapping_labels));
% stim_total_by_muscle = cell(size(emg_mapping_labels));
% stim_mep_by_muscle = cell(size(emg_mapping_labels));
% stim_trials_by_pt = cell(size(emg_mapping_sm));
% 
% for kl = 1:length(emg_mapping_labels)
%     stim_mep_by_muscle{kl} = zeros(8,3);
%     stim_total_by_muscle{kl} = zeros(8,3);
% end
% 
% param_structs = {ep_struct(:).stim_settings};
% 
% for kk = 1:length(mep_delays_sm_all)
%     mono_trials =  parse_params('monopolar_segmented', ep_struct(kk).stim_settings,[1:size(mep_delays_sm{kk},1)], 'all', true, ep_struct(kk).pt_info.lead_model_name);
%     stim_trials_by_pt{kk} = mono_trials;
%     for kl = 1:length(emg_mapping_sm{kk})
%         mono_mep_labels = parse_params('monopolar_segmented', ep_struct(kk).stim_settings,mep_labels_sm{kk}(:,kl), 'all', true, ep_struct(kk).pt_info.lead_model_name);
%     
%         mep_delays_by_muscle{emg_mapping_all{kk}(kl)} = [mep_delays_by_muscle{emg_mapping_all{kk}(kl)}; cell2mat(transpose({mep_delays_sm_all{kk}{:,kl}}))];
%         
%         stim_total_by_muscle{emg_mapping_all{kk}(kl)} = stim_total_by_muscle{kl} + (mono_trials > 0);
%         stim_mep_by_muscle{emg_mapping_all{kk}(kl)} = stim_total_by_muscle{kl} + mono_mep_labels;
% 
%     end
% end
% 
% mep_delays_vplot = padcat(mep_delays_by_muscle{1},mep_delays_by_muscle{2},mep_delays_by_muscle{3},mep_delays_by_muscle{4}, ...
%                     mep_delays_by_muscle{5},mep_delays_by_muscle{6},mep_delays_by_muscle{7},mep_delays_by_muscle{8},...
%                     mep_delays_by_muscle{11},mep_delays_by_muscle{12});
% 
% % for kk = 2:length(mep_delays_by_muscle)
% %     mep_delays_vplot = padcat(mep_delays_vplot,mep_delays_by_muscle{kk});
% % end
% 
% %violin plot of mEP delays for diff muscles
% figure
% violinplot(mep_delays_vplot,emg_mapping_labels([1,2,3,4,5,6,7,8,11,12]),'Bandwidth',0.5)
% xlabel('Muscle')
% ylabel('Delay (ms)')
% 
% ep_trials = 1:length(mep_vals);
% 
% trial_inds = parse_params('monopolar_segmented', ep_struct(kk).stim_settings,[1:size(mep_delays_sm{kk},1)], 'all', true, ep_struct(kk).pt_info.lead_model_name);