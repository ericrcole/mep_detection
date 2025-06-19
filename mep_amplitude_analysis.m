%%% MEP analysis for designing biomarker encoding:
% analyzing variability in amplitudes

load('EP_042023.mat');
load('EPdata_5-17-2020_ANALYZED_postopMRI.mat');

%% how much do muscles vary in mep amplitude?

%first, get unique muscles in mep test group:

all_muscles = {};
for kk = 1:length(ep_struct)
    all_muscles = [all_muscles, ep_struct(kk).mep.emg_labels];
end

display(unique(all_muscles));

%% define muscles for analysis:

use_labeled = true;
ampval = []; 

m_groups = {'ECR','FCR','FDI','bicep','genio','nasal','oris'};
pt_names = {ep_struct(:).patient_ID};
ep_inds = [10,11,12,14,15,17,18,25,26];

amps_all = cell(size(m_groups));
amps_by_pt = cell(length(m_groups),length(ep_struct));

amps_by_amp = cell(length(m_groups),3);
pct_by_amp = nan(length(m_groups),3);
amp_arr = [1,3,5];

if use_labeled
    loopend = length(ep_struct);
else
    loopend = 9;
end

for kk = 1:loopend
    if use_labeled
        tempmat = squeeze(nanmax(abs(EP(ep_inds(kk)).amplitudesEMG),[],1));
    else
        tempmat = ep_struct(kk).mep.amps_detect;
        tempmat(cellfun(@isempty,tempmat)) = {nan};
        tempmat = cellfun(@nanmax, tempmat);
    end
       
    if kk == 8
        emg_lim = 6;
    elseif kk == 9
        emg_lim = 4;
    else
        emg_lim = 8;
    end
    
    for kl = 1:emg_lim
        if use_labeled
            if ~isempty(ampval)
                param_inds = EP(ep_inds(kk)).StimAmplitude == ampval;
                %param_inds = and(param_inds', ep_struct(kk).mep.labels_detect(:,kl));
            else
                param_inds = logical(ones(size(EP(ep_inds(kk)).StimAmplitude)));
                %param_inds = and(param_inds, ep_struct(kk).mep.labels_detect(:,kl));
            end
        else
            if ~isempty(ampval)
                param_inds = ep_struct(kk).stim_settings.amplitudes == ampval;
                param_inds = and(param_inds, ep_struct(kk).mep.labels_detect(:,kl));
            else
                param_inds = logical(ones(size(EP(ep_inds(kk)).StimAmplitude)));
                param_inds = and(param_inds, ep_struct(kk).mep.labels_detect(:,kl));
            end
        end
        
        for km = 1:length(m_groups)
            if use_labeled
                if contains(EP(ep_inds(kk)).channelEMG{kl},m_groups{km})
                    m_ind = km;
                    break
                end
            else
                if contains(ep_struct(kk).mep.emg_labels{kl},m_groups{km})
                    m_ind = km;
                    break
                end
            end
        end
        
        amps_all{m_ind} = [amps_all{m_ind}; tempmat(param_inds,kl)];
        amps_by_pt{m_ind,kk} = [amps_by_pt{m_ind,kk}; tempmat(param_inds,kl)];
        
    end
    
    for kn = 1:length(amp_arr)
        for kl = 1:emg_lim
            if use_labeled
                if ~isempty(ampval)
                    param_inds = EP(ep_inds(kk)).StimAmplitude == amp_arr(kn);
                    %param_inds = and(param_inds, ~isnan(tempmat(param_inds,kl)));
                else
                    %param_inds = and(param_inds, ep_struct(kk).mep.labels_detect(:,kl));
                end
            else
                if ~isempty(ampval)
                    param_inds = ep_struct(kk).stim_settings.amplitudes == amp_arr(kn);
                    param_inds = and(param_inds, ep_struct(kk).mep.labels_detect(:,kl));
                else
                    param_inds = and(param_inds, ep_struct(kk).mep.labels_detect(:,kl));
                end
            end
            
            for km = 1:length(m_groups)
                if use_labeled
                    if contains(EP(ep_inds(kk)).channelEMG{kl},m_groups{km})
                        m_ind = km;
                        break
                    end
                else
                    if contains(ep_struct(kk).mep.emg_labels{kl},m_groups{km})
                        m_ind = km;
                        break
                    end
                end
            end
            
            %amps_all{m_ind} = [amps_all{m_ind}; tempmat(param_inds,kl)];
            amps_by_amp{m_ind,kn} = [amps_by_amp{m_ind,kn}; tempmat(param_inds,kl)];
            
        end
    end
end

%% analysis: amps vs. muscle
addpath('/Users/ERCOLE/Documents/Research/Repos/param_spaces_analysis/Violinplot-Matlab')

figure
subplot(1,2,1)
for kk = 1:length(amps_all)
    Violin(amps_all{kk},kk)
end
xticks(1:length(amps_all))
xticklabels(m_groups)
xtickangle(30)
ylabel('Peak amplitude (\muV)')

subplot(1,2,2)
for kk = 1:length(amps_all)
    scatter(kk, nanmedian(amps_all{kk}),75, 'k','filled')
    hold on
    plot([kk kk], [prctile(amps_all{kk}, 5), prctile(amps_all{kk}, 95)], 'k','LineWidth',2)
end
xticks(1:length(amps_all))
xticklabels(m_groups)
xlim([0,8])
ylabel('Peak amplitude (\muV)')
xtickangle(30)

title('Amplitudes by muscle')

%% analysis by patient/muscle
figure
for kk = 1:length(m_groups)
    subplot(2,4,kk)
    for kl = 1:size(amps_by_pt, 2)
        if ~isempty(amps_by_pt{kk,kl})
            Violin(amps_by_pt{kk,kl},kl);
        end
    end
    title(m_groups{kk})
    xticks(1:size(amps_by_pt, 2))
    xlabel('Patient no.')
    ylabel('Peak amplitude (\muV)')
end

suptitle('Amplitudes by patient')

med_amp_by_pt = cellfun(@nanmedian,amps_by_pt);
amps_by_pt(cellfun(@isempty,amps_by_pt)) = {nan};
max_amp_by_pt = cellfun(@nanmax,amps_by_pt);

figure
subplot(1,2,1)
for kk = 1:length(m_groups)
    scatter(kk, nanmedian(med_amp_by_pt(kk,:)),75, 'k','filled')
    hold on
    plot([kk kk], [prctile(med_amp_by_pt(kk,:), 5), prctile(med_amp_by_pt(kk,:), 95)], 'k','LineWidth',2)
end
ylabel('Median amplitude across patients')
xticks(1:length(amps_all))
xticklabels(m_groups)
xlim([0,8])
xtickangle(30)
set(gca,'FontSize',14)

subplot(1,2,2)
for kk = 1:size(med_amp_by_pt,2)
    scatter(kk, nanmedian(med_amp_by_pt(:,kk)),75, 'k','filled')
    hold on
    plot([kk kk], [prctile(med_amp_by_pt(:,kk), 5), prctile(med_amp_by_pt(:,kk), 95)], 'k','LineWidth',2)
end

ylabel('Median amplitude across muscles')
xlabel('Patient number')
xlim([0,10])
set(gca,'FontSize',14)

% figure
% subplot(1,2,1)
% for kk = 1:length(m_groups)
%     scatter(kk, nanmedian(max_amp_by_pt(kk,:)),75, 'k','filled')
%     hold on
%     plot([kk kk], [prctile(max_amp_by_pt(kk,:), 5), prctile(max_amp_by_pt(kk,:), 95)], 'k','LineWidth',2)
% end
% ylabel('Max amplitude across patients')
% xticks(1:length(amps_all))
% xticklabels(m_groups)
% xlim([0,8])
% xtickangle(30)
% set(gca,'FontSize',14)
% 
% subplot(1,2,2)
% for kk = 1:size(med_amp_by_pt,2)
%     scatter(kk, nanmedian(max_amp_by_pt(:,kk)),75, 'k','filled')
%     hold on
%     plot([kk kk], [prctile(max_amp_by_pt(:,kk), 5), prctile(max_amp_by_pt(:,kk), 95)], 'k','LineWidth',2)
% end
% 
% ylabel('Max amplitude across muscles')
% xlabel('Patient number')
% xlim([0,10])
% set(gca,'FontSize',14)

%%
med_amp_by_pt = cellfun(@nanmedian,cellfun(@log,amps_by_pt,'UniformOutput',false));
amps_by_pt(cellfun(@isempty,amps_by_pt)) = {nan};
max_amp_by_pt = cellfun(@nanmax,cellfun(@log,amps_by_pt,'UniformOutput',false));

med_diff = [nan nan];
var_diff = [nan nan];
skew_diff = [nan nan];

figure
subplot(1,2,1)
for kk = 1:length(m_groups)
    scatter(kk, nanmedian(med_amp_by_pt(kk,:)),75, 'k','filled')
    hold on
    plot([kk kk], [prctile(med_amp_by_pt(kk,:), 5), prctile(med_amp_by_pt(kk,:), 95)], 'k','LineWidth',2)
end
ylabel('Median amplitude across patients')
xticks(1:length(amps_all))
xticklabels(m_groups)
xlim([0,8])
xtickangle(30)
set(gca,'FontSize',14)

subplot(1,2,2)
for kk = 1:size(med_amp_by_pt,2)
    scatter(kk, nanmedian(med_amp_by_pt(:,kk)),75, 'k','filled')
    hold on
    plot([kk kk], [prctile(med_amp_by_pt(:,kk), 5), prctile(med_amp_by_pt(:,kk), 95)], 'k','LineWidth',2)
end

ylabel('Median amplitude across muscles')
xlabel('Patient number')
xlim([0,10])
set(gca,'FontSize',14)
suptitle('Log-scaled amplitudes')

%% Analysis vs. stim amplitude

figure
for kk = 1:length(m_groups)
    subplot(2,4,kk)
    for kl = 1:length(amp_arr)
        if ~isempty(amps_by_amp{kk,kl})
            Violin(amps_by_amp{kk,kl},amp_arr(kl))
        end
    end
    xlim([0.5,5.5])
    xlabel('Current (mA)')
    title(m_groups{kk})
    ylabel('Peak amplitude (\muV)')
end

suptitle('mEP amplitudes by stim current')


figure
for kk = 1:length(m_groups)
    subplot(2,4,kk)
%     for kl = 1:length(amp_arr)
%         if ~isempty(amps_by_amp{kk,kl})
%             Violin(amps_by_amp{kk,kl},kl)
%         end
%     end
    numfun = @(x) sum(~isnan(x));
    temp = cellfun(numfun, amps_by_amp(kk,:));
    plot(amp_arr,temp)
    xlim([0.5,5.5])
    xlabel('Current (mA)')
    title(m_groups{kk})
    ylabel('Setting count')
end

suptitle('Num. responding settings by stim current')