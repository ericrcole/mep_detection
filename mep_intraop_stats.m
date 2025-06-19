% mep_intraop_stats.m
%
% first outline of analysis for the mep paper:
% creates mep data table for stats and tests several hypotheses w/ mixed
% effects models

%% load data

%proc_dir = '/Users/ERCOLE/Documents/Research/Repos/mEP_processing';
% load('EP_040224.mat')
% 
% for kk = length(ep_struct):-1:1
%     %remove ~5 patients/leads with synthetic EMG data
%     if ep_struct(kk).pt_info.impute_emg
%         ep_struct(kk) = [];
%     end
% end

load('EP_112624_final_all.mat')

for kk = length(ep_struct):-1:1
    %remove ~5 patients/leads with synthetic EMG data
    if ep_struct(kk).pt_info.impute_emg || contains(ep_struct(kk).pt_info.dbs_target,'GPi')
        ep_struct(kk) = [];
    end
end

% briefly check range of all tested parameters

% create and save data table
%has following variables: subject, stim parameters (ampl, pulse width,
%contacts), mep amplitude, mep binary response, recorded muscle,
%hemisphere, mono vs. bipolar

save_table = true;
%resp_settings_only = true;  %only keep data for settings that produced an mEP in at least one muscle

mode = 'opt';

dbs_ampl = [];
contact_config = {};
cathode = {};
dbs_pw = [];
pstrings = {};

dbs_target = {};

response = [];
patient_ID = {};
stim_type = {};
muscle_og = {};
muscle_hemi = {};   %ipsi or contra
muscle_full = {};
muscle_group = {};
muscle_loc = {};

mep_ampl = [];
mep_ampl_z = [];
mep_ltc = [];

n_response = [];
mepscore = [];

m_groups = {'ECR','FCR','FDI','bicep','genio','nasal','ocul','oris','hypo','hyo','tib','trap'};
m_facial = {'genio','nasal','oris','hypo','hyo'};
m_limb = {'ECR','FCR','FDI','bicep','tib','trap'};

cathodes_all = {};
anodes_all = {};
dbs_ampl_all = [];
pw_all = [];

total_settings = 0;
total_muscles = 0;

for kk = 1:length(ep_struct)
    n_settings = length(ep_struct(kk).stim_settings.param_strings);
    total_settings = total_settings + n_settings;
    total_muscles = total_muscles + length(ep_struct(kk).mep.emg_labels);
    
    
    amp_mat = ep_struct(kk).mep.(mode).amps_detect;
    amp_fun = @(x) abs(max([x, 0]) - min([x, 0]));
    amp_mat = cellfun(amp_fun, amp_mat,'UniformOutput',false);
    %amp_mat(cellfun(@isempty, amp_mat)) = {nan};
    amp_mat = cell2mat(amp_mat);
    amp_mat(amp_mat == 0) = nan;
    
    ampz_mat = ep_struct(kk).mep.(mode).amps_detect_z;
    ampz_fun = @(x) abs(max([x, 0]) - min([x, 0]));
    ampz_mat = cellfun(amp_fun, ampz_mat,'UniformOutput',false);
    %amp_mat(cellfun(@isempty, amp_mat)) = {nan};
    ampz_mat = cell2mat(ampz_mat);
    ampz_mat(ampz_mat == 0) = nan;
    
    ltc_mat = cellfun(@min, ep_struct(kk).mep.(mode).delays_detect,'UniformOutput',false);
    ltc_mat(cellfun(@isempty, ltc_mat)) = {nan};
    ltc_mat = cell2mat(ltc_mat);
    
    cathode_temp = convert_contacts(ep_struct(kk).stim_settings.cathodes);
    anode_temp = convert_contacts(ep_struct(kk).stim_settings.anodes);
    
    cathodes_all = [cathodes_all; cathode_temp];
    anodes_all = [anodes_all; anode_temp];
    dbs_ampl_all = [dbs_ampl_all; ep_struct(kk).stim_settings.amplitudes];
    pw_all = [pw_all; ep_struct(kk).stim_settings.pulse_widths];
    
    polarity_temp = cell(size(ep_struct(kk).stim_settings.contacts));
    for kl= 1:length(polarity_temp)
        if contains(ep_struct(kk).stim_settings.contacts{kl},"C+")
            polarity_temp{kl} = 'cathodic';
        elseif ep_struct(kk).stim_settings.is_anodic(kl) == 1
            polarity_temp{kl} = 'anodic';
        else
            polarity_temp{kl} = 'bipolar';
        end
    end
    
    target_temp = ep_struct(kk).pt_info.dbs_target;
    
    for kl = 1:length(ep_struct(kk).mep.emg_labels)
        emgtemp = ep_struct(kk).mep.emg_labels{kl};
        emgsplit = strsplit(emgtemp, ' ');
        
        if contains(emgsplit{1},'R') && contains(target_temp, 'R')
            hemi_temp = 'ipsi';
        elseif contains(emgsplit{1},'L') && contains(target_temp, 'L')
            hemi_temp = 'ipsi';
        else
            hemi_temp = 'contra';
        end
        
        if any(arrayfun(@(substr) contains(emgtemp, substr, 'IgnoreCase', true), m_facial))
            mloc_temp = 'facial';
        else
            mloc_temp = 'limb';
        end
        mgroup_temp = m_groups{find(arrayfun(@(substr) contains(emgtemp, substr, 'IgnoreCase', true), m_groups),1)};
        if strcmp(mgroup_temp,'hyo')
           mgroup_temp = 'hypo'; 
%         elseif strcmp(mgroup_temp,'oris')
%            mgroup_temp = 'orb'; 
        end
        
        muscle_hemi = [muscle_hemi; repmat({hemi_temp},[n_settings 1])];
        muscle_og = [muscle_og; repmat({emgtemp},[n_settings 1])];
        muscle_full = [muscle_full; repmat({[hemi_temp, ' ', mgroup_temp]},[n_settings 1])];
        muscle_group = [muscle_group; repmat({mgroup_temp},[n_settings 1])];
        muscle_loc = [muscle_loc; repmat({mloc_temp},[n_settings 1])];

        patient_ID = [patient_ID; repmat({ep_struct(kk).patient_ID},[n_settings 1])];
        dbs_target = [dbs_target; repmat({ep_struct(kk).pt_info.dbs_target},[n_settings 1])];
        dbs_ampl = [dbs_ampl; ep_struct(kk).stim_settings.amplitudes];
        cathode = [cathode; cathode_temp];
        contact_config = [contact_config; ep_struct(kk).stim_settings.contacts];
        dbs_pw = [dbs_pw; ep_struct(kk).stim_settings.pulse_widths];
        stim_type = [stim_type; polarity_temp];
        pstrings = [pstrings; ep_struct(kk).stim_settings.param_strings];
        
        n_response = [n_response; nanmean(ep_struct(kk).mep.(mode).labels_detect,2)];
        mepscore = [mepscore; ep_struct(kk).mep.(mode).mepscore_z];


        response = [response; ep_struct(kk).mep.(mode).labels_detect(:,kl)];
        mep_ampl = [mep_ampl; amp_mat(:,kl)];
        mep_ampl_z = [mep_ampl_z; ampz_mat(:,kl)];
        mep_ltc = [mep_ltc; ltc_mat(:,kl)];
        
        
    end
end

mep_table = table();
mep_table.patient_ID = patient_ID;
mep_table.dbs_target = dbs_target;

mep_table.muscle_full = muscle_og;
mep_table.muscle_adj = muscle_full;
mep_table.muscle_hemi = muscle_hemi;
mep_table.muscle_group = muscle_group;
mep_table.muscle_loc = muscle_loc;

mep_table.n_response = n_response;
mep_table.mepscore = mepscore;

mep_table.response = response;
mep_table.mep_ampl = mep_ampl;
mep_table.mep_ampl_z = mep_ampl_z;
mep_table.mep_ltc = mep_ltc;

mep_table.dbs_ampl = dbs_ampl;
mep_table.pulse_width = dbs_pw;
mep_table.contact = contact_config;
mep_table.cathode = cathode;
mep_table.stim_type = stim_type;
mep_table.param_strings = pstrings;

mep_table(strcmp(mep_table.muscle_adj,'contra tib'),:) = [];
mep_table(strcmp(mep_table.muscle_adj,'ipsi bicep'),:) = [];
mep_table(strcmp(mep_table.muscle_adj,'ipsi FDI'),:) = [];
mep_table(strcmp(mep_table.muscle_adj,'contra trap'),:) = [];
mep_table(strcmp(mep_table.muscle_adj,'ipsi ECR'),:) = [];
mep_table(strcmp(mep_table.muscle_adj,'ipsi FCR'),:) = [];

%save('mep_table_opt.mat','mep_table')

%% clear weird muscles from table
% mep_table(strcmp(mep_table.muscle_adj,'contra tib'),:) = [];
% mep_table(strcmp(mep_table.muscle_adj,'ipsi bicep'),:) = [];
% mep_table(strcmp(mep_table.muscle_adj,'ipsi FDI'),:) = [];

%% summarize which parameters were tested (bipolar only)
bp_clist = {'1','2','2A','2B','2C','3','3A','3B','3C','4'};
invalid_pairs = {'2-2A','2A-2','2-2B','2B-2','2-2C','2C-2',...
    '3-3A','3A-3','3-3B','3B-3','3-3C','3C-3'};

n_settings = zeros(length(bp_clist));

c_map = containers.Map(bp_clist, [1 2 3 4 5 6 7 8 9 10]);

for kk = 1:10
   n_settings(kk,kk) = -1; 
end
for kk = 1:length(invalid_pairs)
    strtemp = strsplit(invalid_pairs{kk},'-');
    n_settings(c_map(strtemp{1}),c_map(strtemp{2})) = -1;
end

for kk = 1:length(cathodes_all)
    if strcmp(cathodes_all{kk},'C') || strcmp(anodes_all{kk},'C')
    
    else
        n_settings(c_map(cathodes_all{kk}),c_map(anodes_all{kk})) = n_settings(c_map(cathodes_all{kk}),c_map(anodes_all{kk}))+1;
        n_settings(c_map(anodes_all{kk}),c_map(cathodes_all{kk})) = n_settings(c_map(anodes_all{kk}),c_map(cathodes_all{kk}))+1;
        
    end
end

fprintf('Valid bipolar combs: %d\n',sum(n_settings(:) ~= -1))

%% mep latencies
addpath('/Users/ERCOLE/Documents/Research/Repos/param_spaces_analysis/Violinplot-Matlab')

m_unique = unique(mep_table.muscle_adj);
ltc_by_muscle = cell(size(m_unique));

for kk = 1:length(m_unique)
    ltc_by_muscle{kk} = mep_table.mep_ltc(strcmp(mep_table.muscle_adj, m_unique{kk}));
end

med_ltcs = cellfun(@nanmedian, ltc_by_muscle);
[~, m_order_ltc] = sort(med_ltcs);
m_sorted = m_unique(m_order_ltc);
ltc_by_muscle = ltc_by_muscle(m_order_ltc);

%[p,tbl,stats] = kruskalwallis(nresp_by_pt);

figure
for kk = 1:length(ltc_by_muscle)
    hold on
    boxchart(kk*ones(size(ltc_by_muscle{kk})),ltc_by_muscle{kk})
    %Violin(ltc_by_muscle{kk},kk)
end
xticks(1:length(m_sorted))
xticklabels(m_sorted)
xtickangle(45)
xlabel('Muscle'); ylabel('Latency (ms)')
set(gca,'FontSize',14)

%% mep amplitude by muscle
addpath('/Users/ERCOLE/Documents/Research/Repos/param_spaces_analysis/Violinplot-Matlab')

%m_unique = unique(mep_table.muscle_adj);
amp_by_muscle = cell(size(m_sorted));

for kk = 1:length(m_sorted)
    amp_by_muscle{kk} = mep_table.mep_ampl(strcmp(mep_table.muscle_adj, m_sorted{kk}));
end

% med_ltcs = cellfun(@nanmedian, ltc_by_muscle);
% [~, m_order_ltc] = sort(med_ltcs);
% m_sorted = m_unique(m_order_ltc);
%amp_by_muscle = amp_by_muscle(m_order_ltc);

amp_by_muscle_z = cell(size(m_sorted));

for kk = 1:length(m_sorted)
    amp_by_muscle_z{kk} = mep_table.mep_ampl_z(strcmp(mep_table.muscle_adj, m_sorted{kk}));
end

%amp_by_muscle_z = amp_by_muscle_z(m_order_ltc);


%figure
subplot(2,2,2)
for kk = 1:length(amp_by_muscle)
    hold on
    boxchart(kk*ones(size(amp_by_muscle{kk})),amp_by_muscle{kk},'BoxFaceColor',plotcols{kk},'MarkerColor',plotcols{kk})
    %Violin(amp_by_muscle{kk},kk)
end
xticks((1:length(m_sorted))+0.5)
xticklabels(m_sorted_labels)
xtickangle(25)
ylabel('Amplitude (\muV)')
legend({'Facial','','','','Limb','','',''},'Location','Northwest')
set(gca,'FontSize',14)
ylim([0,150])

% subplot(1,2,2)
% for kk = 1:length(amp_by_muscle_z)
%     hold on
%     boxchart(kk*ones(size(amp_by_muscle_z{kk})),amp_by_muscle_z{kk})
%     %Violin(amp_by_muscle{kk},kk)
% end
% xticks(1:length(m_sorted))
% xticklabels(m_sorted)
% xtickangle(45)
% xlabel('Muscle'); ylabel('Amplitude (z-scored)')
% set(gca,'FontSize',14)
% ylim([0 150])

amp_dist_z = [];

for kk = 1:length(emg_keys)
    temp = abs(amp_by_muscle_z{kk});
    % boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
    % hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    amp_dist_z = [amp_dist_z; prctile(temp,25), prctile(temp,75), nanmin(temp)];
end

amp_dist = amp_dist_z;
mean_stats = mean(amp_dist);

stat_adj_z = zeros(size(amp_dist,1),1); %differences needed to adjust each muscle to target mean distribution
for kl = 1:size(stat_adj_z,1)
    stat_adj_z(kl) = (mean_stats(2) - mean_stats(1))./(amp_dist(kl,2) - amp_dist(kl,1));
end

subplot(2,2,4)
for kk = 1:length(amp_by_muscle_z)
    temp = abs(amp_by_muscle_z{kk});
    temp = temp - amp_dist_z(kk,3);
    
    temp = temp*stat_adj_z(kk);
    temp = temp + 8;

    boxchart(ones(size(temp))*kk,temp,'BoxFaceColor',plotcols{kk},'MarkerColor',plotcols{kk});
    hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    
end
xticks((1:length(m_sorted))+0.5)
xticklabels(m_sorted_labels)
xtickangle(25)
ylim([0,60])
xlabel('Muscle')
ylabel('Amplitude (norm.)')
set(gca,'FontSize',14)

%% stats for fig 3 (response prob/latency/amplitude by muscle)
temp = mep_table.muscle_loc;
temp2 = zeros(size(temp));
temp2(strcmp(temp, 'limb')) = 1;

mep_table.loc_code = temp2;

disp('Latency test (facial vs limb):')
lme = fitlme(mep_table,'mep_ltc ~ loc_code + (loc_code|patient_ID)')

disp('Amplitude test (facial vs limb):')
lme = fitlme(mep_table,'mep_ampl ~ loc_code + (loc_code|patient_ID)')

disp('Response prob. test (facial vs limb):')
lme = fitglme(mep_table,'response ~ loc_code + (loc_code|patient_ID)',Distribution="Binomial")

%% How does mep amp vary by both pt/muscle?
unique_pts = unique(mep_table.patient_ID);
unique_muscles = m_sorted;

amps_z_by_ptm = cell(length(unique_pts), length(unique_muscles));
amps_by_ptm = cell(length(unique_pts), length(unique_muscles));

for kk = 1:length(unique_pts)
    for kl = 1:length(unique_muscles)
        inds = and(strcmp(mep_table.patient_ID,unique_pts{kk}),strcmp(mep_table.muscle_adj,m_sorted{kl}));
        sub_table = mep_table(inds,:);
        if ~isempty(sub_table)
            amps_by_ptm{kk,kl} = sub_table.mep_ampl;
            amps_z_by_ptm{kk,kl} = sub_table.mep_ampl_z;
        end
    end
end

amps_by_ptm = transpose(amps_by_ptm);
mvar = cellfun(@iqr, amps_by_ptm);
% mvarz =cellfun(@iqr, amps_z_by_ptm);
% mvar_across_pt = cat(1, amps_by_ptm{:});
% mvar_across_muscle = cat(1, amps_by_ptm{:});
mvar_across_muscle = nanmedian(mvar,1);
mvar_across_pt = nanmedian(mvar,2);

figure; 
subplot(1,2,1)
bar(mvar_across_muscle)
xticks(1:length(m_sorted))
xticklabels(m_sorted)
xtickangle(30)

subplot(1,2,2)
bar(mvar_across_pt)

% xticks(1:length(m_sorted))
% xticklabels(m_sorted)
% xtickangle(30)
%%

figure
subplot(1,2,1)
boxchart(mvar_across_pt)


%% probability of response for each muscle?
%load('mep_table_opt.mat')
%filter by settings that produced a response in at least one muscle

m_sorted = [{'contra oris'   }
    {'ipsi oris' }
    {'contra nasal'}
    {'ipsi nasal'  }
    {'contra bicep'}
    {'contra FCR'  }
    {'contra ECR'  }
    {'contra FDI'  }];

unique_settings = unique(mep_table(:,[1,19]),'rows');
for kk = 1:size(unique_settings,1)
    inds = and(strcmp(mep_table.patient_ID,unique_settings.patient_ID{kk}),strcmp(mep_table.param_strings,unique_settings.param_strings{kk}));
    if ~any(mep_table.response(inds))
       mep_table(inds,:) = []; 
    end
end

% unique_settings = unique(mep_table(:,[1,17]),'rows');
%
ntot_by_muscle = nan(size(m_sorted));
nresp_by_muscle = nan(size(m_sorted));

unique_pts = unique(mep_table.patient_ID);
nresp_by_pt = nan(length(unique_pts), length(m_sorted));
ntot_by_pt = nan(length(unique_pts), length(m_sorted));

for kk = 1:length(m_sorted)
    ntot_by_muscle(kk) = sum(strcmp(mep_table.muscle_adj, m_sorted{kk}));
    nresp_by_muscle(kk) = 100*nanmean(mep_table.response(strcmp(mep_table.muscle_adj, m_sorted{kk})));

    for kl = 1:length(unique_pts)
        pttable = mep_table(strcmp(mep_table.patient_ID, unique_pts{kl}),:);

        if isempty(pttable)
            continue
        end

        ntot_by_pt(kl,kk) = sum(strcmp(pttable.muscle_adj, m_sorted{kk}));
        nresp_by_pt(kl,kk) = 100*nanmean(pttable.response(strcmp(pttable.muscle_adj, m_sorted{kk})));

    end
end

%figure
% subplot(2,1,1)
% bar(nresp_by_muscle)
% ylabel('Response rate (pct.)')
% xticks(1:length(m_sorted))
% xticklabels(m_sorted)
% xtickangle(45)
% set(gca,'FontSize',14)

subplot(2,1,1)
for kk = 1:length(m_sorted)
    hold on
    temp = nresp_by_pt(:,kk);
    boxchart(kk*ones(size(temp)),temp)
    %Violin(ltc_by_muscle{kk},kk)
end
ylabel('Response rate (pct.)')
xticks(1:length(m_sorted))
xticklabels(m_sorted)
xtickangle(45)
set(gca,'FontSize',14)

subplot(2,1,2)
bar(ntot_by_muscle)
ylabel('Total no. settings')
xticks(1:length(m_sorted))
xticklabels(m_sorted)
xtickangle(45)
set(gca,'FontSize',14)

[p,tbl,stats] = kruskalwallis(nresp_by_pt);


%% fig. 4: top part (response probabilities and latency by muscle)
%load('mep_table_opt.mat')
%mep_table = mep_table_v1;
mep_table(contains(mep_table.muscle_adj,'ocul'),:) = [];

m_sorted = [{'contra oris'   }
    {'ipsi oris' }
    {'contra nasal'}
    {'ipsi nasal'  }
    {'contra bicep'}
    {'contra FCR'  }
    {'contra ECR'  }
    {'contra FDI'  }];

m_sorted_labels = [{'CL oris'   }
    {'IL oris' }
    {'CL nasalis'}
    {'IL nasalis'  }
    {'CL bicep'}
    {'CL FCR'  }
    {'CL ECR'  }
    {'CL FDI'  }];

plotcols = {
    [0.7500    0.5250    0.0980],
    [0.7500    0.5250    0.0980],
    [0.7500    0.5250    0.0980],
    [0.7500    0.5250    0.0980],
    [0    0.4470    0.7410],
    [0    0.4470    0.7410],
    [0    0.4470    0.7410],
    [0    0.4470    0.7410]};

m_unique = unique(mep_table.muscle_adj);
ltc_by_muscle = cell(size(m_sorted));

for kk = 1:length(m_sorted)
    disp(kk)
    ltc_by_muscle{kk} = mep_table.mep_ltc(strcmp(mep_table.muscle_adj, m_sorted{kk}));
end

med_ltcs = cellfun(@nanmedian, ltc_by_muscle);
[~, m_order_ltc] = sort(med_ltcs);
%m_sorted = m_unique(m_order_ltc);
ltc_by_muscle = ltc_by_muscle(m_order_ltc);

subplot(2,2,1)
for kk = 1:length(m_sorted)
    hold on
    temp = nresp_by_pt(:,kk);
    boxchart(kk*ones(size(temp)),temp,'BoxFaceColor',plotcols{kk},'MarkerColor',plotcols{kk})
    %Violin(ltc_by_muscle{kk},kk)
end
ylabel('Response rate (pct.)')
xticks((1:length(m_sorted))+0.5)
xticklabels(m_sorted_labels)
xtickangle(30)
set(gca,'FontSize',14)

subplot(2,2,3)
for kk = 1:length(ltc_by_muscle)
    hold on
    boxchart(kk*ones(size(ltc_by_muscle{kk})),ltc_by_muscle{kk},'BoxFaceColor',plotcols{kk},'MarkerColor',plotcols{kk})
    %Violin(ltc_by_muscle{kk},kk)
end
xticks((1:length(m_sorted))+0.5)
xticklabels(m_sorted_labels)
%xticks([])
xtickangle(30)
xlabel('Muscle'); ylabel('Latency (ms)')
set(gca,'FontSize',14)

unique_settings = unique(mep_table(:,[1,19]),'rows');
for kk = 1:size(unique_settings,1)
    inds = and(strcmp(mep_table.patient_ID,unique_settings.patient_ID{kk}),strcmp(mep_table.param_strings,unique_settings.param_strings{kk}));
    if ~any(mep_table.response(inds))
       mep_table(inds,:) = []; 
    end
end


% unique_settings = unique(mep_table(:,[1,17]),'rows');
%
% ntot_by_muscle = nan(size(m_unique));
% nresp_by_muscle = nan(size(m_unique));

%mep_table(isnan(mep_table.response),:) = [];

% for kk = 1:length(m_sorted)
%     ntot_by_muscle(kk) = sum(strcmp(mep_table.muscle_adj, m_sorted{kk}));
%     nresp_by_muscle(kk) = 100*nanmean(mep_table.response(strcmp(mep_table.muscle_adj, m_sorted{kk})));
% end
% 
% subplot(2,1,1)
% bar(nresp_by_muscle)
% ylabel('Response rate (pct.)')
% xticks(1:length(m_sorted))
% xticklabels(m_sorted)
% %xticks([])
% xtickangle(30)
% xlim([0,11])
%set(gca,'FontSize',14)

save_amp_dist = false;

%% old version of amplitude fig
figure
subplot(2,2,2)
categ_inds = [5,7,1,2,3,4];
categs = emg_keys(categ_inds);
categ_labels = {'Nasalis','Orb. Oris','Bicep', 'ECR', 'ECR', 'FDI'};

ind = 1;
for kk = categ_inds
    temp = abs(all_amps_by_muscle_raw{kk});
    
    boxchart(ones(size(temp))*ind,temp,'MarkerStyle','none');
    hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    amp_dist_z = [amp_dist_z; prctile(temp,25), prctile(temp,75), nanmin(temp)];
    ind = ind + 1;
end
xticks(1:length(categs)+0.5)
xticklabels(categ_labels)
ylim([0 110])
set(gca, 'FontSize',14)
%xlabel('Muscle')
xlabel([])
ylabel('Amplitude (\muV)')

%
amp_dist_log = [];
amp_dist_z = [];

%%%%%
%figure
%subplot(3,1,2)
categ_inds = [7,5,1,2,3,4];
categs = emg_keys(categ_inds);
categ_labels = {'Nasalis','Orb. Oris','Bicep', 'ECR', 'ECR', 'FDI'};

ind = 1;
for kk = categ_inds
    temp = abs(all_amps_by_muscle{kk});

    %disp(categs{ind})
    %boxchart(ones(size(temp))*ind,temp,'MarkerStyle','none');
    %hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    amp_dist_z = [amp_dist_z; prctile(temp,25), prctile(temp,75), nanmin(temp)];
    ind = ind + 1;
end
% xticks(1:length(categs)+0.5)
% xticklabels(categ_labels)
% ylim([0 50])
% set(gca, 'FontSize',12)
% xlabel('Muscle')
% ylabel('Amplitude')
%title('Raw Z-amplitudes')
%%%%%


% subplot(4,1,4)
% for kk = 1:length(emg_keys)
%     temp = abs(all_amps_by_muscle{kk});
%     if keys_facial(kk)
%         temp(temp<5) = 5;
%         temp = log2(1+temp - 5)+3;
%     else
%         temp(temp<8) = 8;
%         temp = log2(1+temp - 8)+3;
%     end
% 
%     %temp = log2(1+abs(all_amps_by_muscle{kk}));
%     boxchart(ones(size(temp))*kk,temp,'MarkerStyle','none');
%     hold on
% 
%     fprintf('%s: median = %.2f, SD = %.2f\n',emg_keys{kk},nanmedian(temp), nanstd(temp))
%     %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
% 
%     amp_dist_log = [amp_dist_log; nanmean(temp), nanstd(temp)];
% end
% fprintf('\n')
% xticks(1:length(categs))
% xticklabels(categs)
% xlabel('Muscle')
% ylabel('Amplitude (log-scaled)')
% title('Log-scaled Z-amplitudes')


amp_dist = amp_dist_log;
mean_stats = mean(amp_dist);

stat_adj_log = zeros(size(amp_dist)); %differences needed to adjust each muscle to target mean distribution
for kl = 1:size(stat_adj_log,1)
    stat_adj_log(kl,1) = mean_stats(1) - amp_dist(kl,1);
    stat_adj_log(kl,2) = mean_stats(2)./amp_dist(kl,2);
end

mep_amp_stats = struct();
mep_amp_stats.keys = emg_keys;
mep_amp_stats.mean_stats_log = mean_stats;
mep_amp_stats.all_stats_log = amp_dist;
mep_amp_stats.stat_adj_log = stat_adj_log;

amp_dist = amp_dist_z;
mean_stats = mean(amp_dist);

stat_adj_z = zeros(size(amp_dist,1),1); %differences needed to adjust each muscle to target mean distribution
for kl = 1:size(stat_adj_z,1)
    stat_adj_z(kl) = (mean_stats(2) - mean_stats(1))./(amp_dist(kl,2) - amp_dist(kl,1));
end

mep_amp_stats.mean_stats_z = mean_stats;
mep_amp_stats.all_stats_z = amp_dist;
mep_amp_stats.stat_adj_z = stat_adj_z;

if save_amp_dist
    save('MEP_ampstats_opt.mat','mep_amp_stats')
end


subplot(2,2,4)
%categs = emg_keys;
tempind = 1;
for kk = categ_inds
    temp = abs(all_amps_by_muscle{kk});
    temp = temp - amp_dist_z(tempind,3);
    
    temp = temp*stat_adj_z(tempind);
    temp = temp + 8;

    boxchart(ones(size(temp))*tempind,temp);
    hold on
    %scatter(ones(size(temp))*kk+randn(size(temp))*0.1, temp, 20, [0.5 0.5 0.5],'filled')
    tempind = tempind + 1;
end
xticks(1:length(categs)+0.5)
xticklabels(categ_labels)
set(gca, 'FontSize',14)
xlabel('Muscle')
ylabel('Amplitude (norm.)')
%title('Adjusted Z-amplitudes')
ylim([0 50])

%% create table w/ muscle statistics

load('mep_table_opt.mat')
muscle_info = table();

unique_muscles = unique(mep_table.muscle_adj);
muscle_info.muscle = unique_muscles;

n_patients = zeros(size(unique_muscles));
num_settings = zeros(size(unique_muscles));
num_excluded = zeros(size(unique_muscles));

for kk = 1:length(unique_muscles)
    temptable = mep_table(strcmp(mep_table.muscle_adj, unique_muscles{kk}),:);
    n_patients(kk) = size(unique(temptable.patient_ID),1);
    num_settings(kk) = size(temptable,1);
    num_excluded(kk) = sum(isnan(temptable.response));
end

muscle_info.num_patients = n_patients;
muscle_info.num_settings = num_settings;
muscle_info.pct_excluded = 100*num_excluded./num_settings;

writetable(muscle_info,'muscle_info.csv')

%% fig. 5ab: does stim amplitude affect mEP amplitude?
% addpath('/Users/ERCOLE/Documents/Research/Repos/densityScatterChart')
% pct variability between muscles?
%amp_model = fitglme(mep_table, 'mep_ampl ~ 1 + dbs_ampl*muscle_adj + (dbs_ampl*muscle_adj|patient_ID) + (dbs_ampl*muscle_adj|contact)');
% amp_model = fitglme(mep_table, 'mep_ampl ~ dbs_ampl + (dbs_ampl|patient_ID) + (dbs_ampl|muscle_adj)');
% 
% disp(amp_model)

save_fig =  false;
cols = [transpose(linspace(0,1,6)), [0;0;0;0;0;0], [0;0;0;0;0;0]];

load('mep_table_opt.mat')
%mep_table(mep_table.dbs_ampl > 5,:) = [];
mep_table(mep_table.mep_ampl > 400,:) = [];

figure('Position',[100 100 400 600])
subplot(2,1,1)
scatter(mep_table.dbs_ampl+ randn(size(mep_table.dbs_ampl))*.05, mep_table.mep_ampl,'filled')
% for kk = 1:6
%     temp = mep_table.mep_ampl(mep_table.dbs_ampl == kk);
%     boxchart(kk*ones(size(temp)),temp, 'BoxFaceColor',cols(kk,:),'MarkerColor',cols(kk,:));
%     hold on
% end
xlim([0,6.5])

set(gca, 'FontSize', 16)
xticks([0,1,2,3,4,5,6]);
% xlabel('DBS amplitude (mA)')
% ylabel('mEP amplitude (\muV)')
%set(gca, 'FontSize', 12)

%densityScatterChart(mep_table.dbs_ampl, mep_table.mep_ampl)

% subplot(2,1,1)
% scatter(mep_table.dbs_ampl + randn(size(mep_table.dbs_ampl))*.05, mep_table.mep_ampl,'filled')


%amp_model = fitlme(mep_table, 'mep_ampl ~ dbs_ampl + (dbs_ampl|patient_ID) + (dbs_ampl|cathode)');
amp_model = fitlme(mep_table, 'mep_ampl ~ dbs_ampl + (dbs_ampl|patient_ID) + (dbs_ampl|muscle_adj) ');
%disp(amp_model)

fprintf('Peak amp. vs. DBS amp: %.4e\n',amp_model.Coefficients.pValue(2))

lcoeffs = amp_model.Coefficients.Estimate;

hold on
plot([0, 7], lcoeffs(1)+lcoeffs(2)*[0,7], 'k--')

% xlim([0 6.5])
% ylim([0 100])

temptable = mep_table;
keyCols = temptable(:, {'patient_ID', 'param_strings'});
disp(size(temptable));



% Find unique combinations of the two columns
[~, ia] = unique(keyCols, 'rows', 'stable');

% Keep only rows with unique column combinations
temptable = temptable(ia, :);
disp(size(temptable));

setting_table = temptable;

subplot(2,1,2)
scatter(temptable.dbs_ampl + randn(size(temptable.dbs_ampl))*.1, temptable.n_response,'filled')
% for kk = 1:6
%     temp = temptable.n_response(temptable.dbs_ampl == kk);
%     boxchart(kk*ones(size(temp)),temp,'BoxFaceColor',cols(kk,:),'MarkerColor',cols(kk,:));
%     hold on
% end
xlim([0,6.5])
xticks([0,1,2,3,4,5,6]);

amp_model = fitlme(temptable, 'n_response ~ dbs_ampl + (dbs_ampl|patient_ID)');
%disp(amp_model)

fprintf('Response pct. vs. DBS amp: %.4e\n',amp_model.Coefficients.pValue(2))

lcoeffs = amp_model.Coefficients.Estimate;

hold on
plot([0, 7], lcoeffs(1)+lcoeffs(2)*[0,7], 'k--')

set(gca, 'FontSize', 16)

if save_fig
    savestr = 'fig6_ab';
    exportgraphics(gcf,sprintf('/Users/ercole/Documents/mep_paper/figs/%s.png',savestr),'Resolution',300)
end
% xlabel('DBS amplitude (mA)')
% ylabel('Proportion of muscles')

%% mep score anlysis vs dbs amplitude
figure
%scatter(mep_table.dbs_ampl + randn(size(mep_table.dbs_ampl))*.05, mep_table.n_response,'filled')
for kk = 1:6
    temp = temptable.mepscore(temptable.dbs_ampl == kk);
    boxchart(kk*ones(size(temp)),temp,'BoxFaceColor',cols(kk,:),'MarkerColor',cols(kk,:));
    hold on
end

amp_model = fitlme(temptable, 'mepscore ~ dbs_ampl + (dbs_ampl|patient_ID)');
%disp(amp_model)

fprintf('mEP score vs. DBS amp: %.4e\n',amp_model.Coefficients.pValue(2))

% plot([0, 15], lcoeffs(1)+lcoeffs(2)*[0,15], 'k--')
% xlabel('Stim amplitude (mA)')
% ylabel('Proportion of muscles')
% set(gca, 'FontSize', 12)
% xlim([0 6.5])
% ylim([0 1])

%% check/summarize which patients had which types of settings tested:
unique_pts = unique(mep_table.patient_ID);
for kk = 1:length(unique_pts)
    subtable = mep_table(strcmp(mep_table.patient_ID, unique_pts{kk}),:);

    %disp(unique_pts{kk})

    fprintf('%s: %d settings\n', unique_pts{kk}, length(ep_struct(strcmp({ep_struct.patient_ID},unique_pts{kk})).stim_settings.param_strings));
    disp(unique(subtable.cathode));
    disp(unique(subtable.dbs_ampl))
end

segment_pts = {'ephys014','ephys015','ephys017','ephys020','ephys025','ephys026','ephys034','ephys062','ephys063','ephys077','ephys085'};

monopolar_pts = {'ephys002','ephys007','ephys008','ephys010','ephys011','ephys014','ephys015','ephys017','ephys020','ephys025','ephys026','ephys030','ephys034','ephys062','ephys063','ephys085','ephys088_R','ephys098'};

%% fig. 5c: plotting mep score1+ vs stim contact
amp_groups = {[0 2], [2 4], [4 6.5]};
gnames = {'Low','Mid','High'};
cols = {'k','b',[255 74 20]/255};

type = 'monopolar';
save_fig = true;

use_seg_patients = false;
use_mono_patients = false;

if strcmp(type,'segmented')
    sps = [11,7,8,9,4,5,6,2];
    n_settings = 8;
    caths = {'1','2B','2A','2C','3B','3A','3C','4'};
    cath_labels = {'1','2M','2A','2L','3M','3A','3L','4'};
    posarr = [100 100 600 600];
    xlablist = [11,7,9];
    savestr = 'fig6c_seg_v1';
else
    sps = [4,3,2,1];
    n_settings = 4;
    caths = {'1','2','3','4'};
    cath_labels = {'1','2','3','4'};
    posarr = [100 100 225 600];
    xlablist = [4];
    savestr = 'fig6d_mono_v1';
end

if strcmp(type,'segmented')
    scores_by_all_contact = cell(0);
end
scores_all = cell(n_settings,3);
scores_by_all_contact = [];

figure('Position',posarr)
counter = 1;
for kk = 1:n_settings
    temptable = setting_table;

    if use_seg_patients
        
        rowsToKeep = ismember(temptable.patient_ID, segment_pts);

        % Filter the table
        temptable = temptable(rowsToKeep, :);
        
    elseif use_mono_patients
        rowsToKeep = ismember(temptable.patient_ID, monopolar_pts);

        % Filter the table
        temptable = temptable(rowsToKeep, :);
    end
    temptable = temptable(cellfun(@(x) contains(x,'C+'), temptable.contact),:);
    temptable = temptable(strcmp(temptable.cathode,caths{kk}),:);

    if strcmp(type,'segmented')
        subplot(4,3,sps(kk))
    else
        subplot(4,1,sps(kk))
    end

    for kl = 1:length(amp_groups)
        amptemp = amp_groups{kl};
        subtable = temptable;
        subtable = subtable(and(subtable.dbs_ampl<=amptemp(2), subtable.dbs_ampl>=amptemp(1)),:);
        
        tempscores = subtable.mepscore;

        scores_all{kk,kl} = tempscores;
        if kl == 2
            scores_by_all_contact = [scores_by_all_contact; {tempscores}];
        end
        
        boxchart(ones(size(tempscores))*kl, tempscores, 'BoxFaceColor',cols{kl},'MarkerStyle','none');
        hold on

        title(cath_labels{kk});
        yticks([0,2,4,6])
        xticks(1:3);
        if ismember(sps(kk),xlablist)
            xticklabels(gnames)
        else
            xticklabels([])
        end
        ylim([0 8])
        set(gca,'FontSize',16)

        counter = counter + 1;
    end
end

if save_fig
    exportgraphics(gcf,sprintf('/Users/ercole/Documents/mep_paper/figs/%s.png',savestr),'Resolution',300)
end

%%
%scores_by_all_contact = scores_by_all_contact([1 9 2 3 4 10 5 6 7 8]);

%ch_labels = {'1','2','2A','2B','2C','3','3A','3B','3C','4'};
ch_labels = {'1','2A','2B','2C','3A','3B','3C','4'};

figure
for kk = 1:length(scores_by_all_contact)
    tempscores = scores_by_all_contact{kk};
    hold on
    boxchart(ones(size(tempscores))*kk, tempscores, 'BoxFaceColor',cols{2},'MarkerStyle','none');
end

xticks(1:length(scores_by_all_contact))
xlim([0, length(scores_by_all_contact)+1])
xticklabels(ch_labels)
ylim([0 5.5])


% statistical comparison of mep score for different contacts separated for 3/5 mA
q = nanpad(scores_by_all_contact);

[p, tbl, stats] = kruskalwallis(q)

c = multcompare(stats)

%% box plot for mep score vs. cathodes/anodes
n_settings = 10;
caths = {'1','2','2B','2A','2C','3','3B','3A','3C','4'};
cath_labels = {'1','2','2M','2A','2L','3','3M','3A','3L','4'};
amps = [3,5];

use_seg_patients = false;

scores_by_all_contact = cell(length(amps),length(caths));

figure('Position',[100 100 1500 350])
for k1 = 1:length(amps)
    subtable = setting_table;
    if use_seg_patients
        disp(size(subtable));
        rowsToKeep = ismember(subtable.patient_ID, segment_pts);

        % Filter the table
        subtable = subtable(rowsToKeep, :);
        disp(size(subtable));
    else
        
    end
    %subtable = subtable(and(subtable.dbs_ampl<=amptemp(2), subtable.dbs_ampl>=amptemp(1)),:);
    subtable = subtable(subtable.dbs_ampl == amps(k1),:);
    for k2 = 1:length(caths)
        
        temptable = subtable(cellfun(@(x) contains(x,'C+'), subtable.contact),:);
        temptable = temptable(strcmp(temptable.cathode,caths{k2}),:);
        tempscores = temptable.mepscore;

        scores_by_all_contact{k1,k2} = tempscores;

    end
end

subplot(1,3,1)
for kk = 1:size(scores_by_all_contact,2)
    tempscores = scores_by_all_contact{1,kk};
    
    boxchart(ones(size(tempscores))*kk, tempscores, 'BoxFaceColor',cols{2},'MarkerStyle','none');
    hold on
end
xticks(1:length(caths))
xticklabels(cath_labels)
xtickangle(0)
%title('3 mA')
%ylabel('mEP Score')
set(gca,'FontSize',16)

q = nanpad(scores_by_all_contact(1,:));

% [p, tbl, stats] = kruskalwallis(q)
% 
% c = multcompare(stats)

subplot(1,3,2)
for kk = 1:size(scores_by_all_contact,2)
    tempscores = scores_by_all_contact{2,kk};
    
    boxchart(ones(size(tempscores))*kk, tempscores, 'BoxFaceColor',cols{3},'MarkerStyle','none');
    hold on
end
xticks(1:length(caths))
xticklabels(cath_labels)
%title('5 mA')
xtickangle(0)
%ylabel('mEP Score')
set(gca,'FontSize',16)

% stats for mEP by stim contact:

prox_contacts = {'1','2C','3C','4'};

load('mep_table_opt.mat');
%mep_table = mep_table(mep_table.dbs_ampl == 3, :);
mep_table = mep_table(and(mep_table.dbs_ampl < 4,mep_table.dbs_ampl >2),:);
mep_table = mep_table(cellfun(@(x) contains(x,'C+'), mep_table.contact),:);
%mep_table = mep_table(mep_table.dbs_ampl>2,:);
% temptable = setting_table(setting_table.dbs_ampl == 5, :);
temptable = setting_table;

mep_table.prox_to_ic = contains(mep_table.cathode, prox_contacts);

temptable.prox_to_ic = contains(temptable.cathode, prox_contacts);
%temptable = temptable(and(temptable.dbs_ampl < 4,temptable.dbs_ampl >2),:);
temptable = temptable(temptable.dbs_ampl == 3,:);
temptable = temptable(cellfun(@(x) contains(x,'C+'), temptable.contact),:);

%figure
subplot(1,3,3)

temp = temptable.mepscore(temptable.prox_to_ic);
boxchart(ones(size(temp))*1, temp,'BoxFaceColor',[0.8 0 0],'MarkerColor',[0.8 0 0]);
hold on
temp = temptable.mepscore(~temptable.prox_to_ic);
boxchart(ones(size(temp))*2, temp,'BoxFaceColor',[0.4 0 0],'MarkerColor',[0.8 0 0]);
hold on

xticks([1,2])
xticklabels({'Near IC', 'Far from IC'})
%ylabel('mEP Score')
set(gca,'FontSize',16)

%figure

% temp = temptable.mepscore(temptable.prox_to_ic);
% boxchart(ones(size(temp))*1, temp);
% hold on
% 
% temp = temptable.mepscore(~temptable.prox_to_ic);
% boxchart(ones(size(temp))*2, temp);
% hold on
% 
% xticks([1,2])
% xticklabels({'Near IC', 'Far from IC'})
% ylabel('mEP Score')

disp('mEP response test (near vs far from capsule):')
%lme = fitlme(temptable,'mepscore ~ prox_to_ic')
%lme = fitglme(mep_table,'response ~ prox_to_ic + (prox_to_ic|patient_ID) + (prox_to_ic|muscle_adj)','Distribution','binomial')


lme = fitglme(mep_table,'response ~ prox_to_ic + (prox_to_ic|patient_ID)+ (prox_to_ic|muscle_adj)','Distribution','binomial')
fprintf('mEP response prob. test (near vs far from capsule): %.3f\n', lme.Coefficients.pValue(2));

disp('mEP Score test (near vs far from capsule):')
lme = fitlme(temptable,'mepscore ~ prox_to_ic + (prox_to_ic|patient_ID)')

fprintf('mEP Score test (near vs far from capsule): %.3f\n', lme.Coefficients.pValue(2));

if save_fig
    savestr = 'fig6_ef';
    exportgraphics(gcf,sprintf('/Users/ercole/Documents/mep_paper/figs/%s.png',savestr),'Resolution',300)
end

%% does mep amplitude vs. stim amplitude effect change vs. muscle?

save_fig = true;

% visualizing and modeling mEP thresholds?
load('mep_table.mat')
mep_table(mep_table.dbs_ampl > 5,:) = [];
%mep_table(mep_table.mep_ampl > 1000,:) = [];

thr_model = fitglme(mep_table, 'response ~ 1 + (1|patient_ID) + (1|muscle_group)','Distribution','binomial');
%thr_model = fitglme(mep_table, 'response ~ dbs_ampl','Distribution','binomial');
disp(thr_model)

sig_coeffs = thr_model.Coefficients.Estimate;

amps_on = mep_table.dbs_ampl(mep_table.response==1);
amps_off = mep_table.dbs_ampl(mep_table.response==0);

figure
hold on
histogram(amps_on,20,'FaceColor','red','Normalization','pdf')
histogram(amps_off,20,'FaceColor','blue','Normalization','pdf')
xlim([-1.5,10])
xlabel('Stim Amplitude')
legend({'mEP', 'No mEP'},'Location','Northeast')
ylabel('Normalized frequency')
set(gca,'FontSize',14)

%histogram(amps_off,20,'FaceColor','blue','Normalization','pdf')
[f_on,x_on] = ksdensity(amps_on,'Bandwidth',1);
[f_off,x_off] = ksdensity(amps_off,'Bandwidth',1);

plotlim = max([f_on, f_off]);

figure
hold on
q1 = area(x_on,f_on,'FaceColor','red');
q1.FaceAlpha = 0.4;

q2 = area(x_off,f_off,'FaceColor','blue');
q2.FaceAlpha = 0.4;

sigmoid = @(coeffs, x) exp(coeffs(1) + coeffs(2)*x)./(1+ exp(coeffs(1) + coeffs(2)*x));
sigmoid_int = @(coeffs, x) exp(coeffs(1) + x)./(1+ exp(coeffs(1) + x));

x_sig = linspace(-2.5,15,100);
plot(x_sig, plotlim*sigmoid_int(-2.4,x_sig),'k')

xlim([-1.5,10])
xlabel('Stim Amplitude')
legend({'mEP', 'No mEP'},'Location','Northeast')
ylabel('Probability (normalized)')
set(gca,'FontSize',14)

%histogram(amps_on,20,'FaceColor','red','Normalization','pdf')


%% testing thresholds for individual contacts

unique_contacts = unique(mep_table.cathode);

labels = unique_contacts;
labels{end} = 'bipolar';

ints = nan(size(unique_contacts));
ints_ci = nan(length(unique_contacts),2);
ints_p = nan(size(unique_contacts));

for kk = 1:length(unique_contacts)
    if kk<length(unique_contacts)
        inds = and(strcmp(mep_table.stim_type,'cathodic'),strcmp(mep_table.cathode,unique_contacts{kk}));
    else
        inds = strcmp(mep_table.stim_type,'bipolar');
    end
    sub_table = mep_table(inds,:);
    
    thr_model = fitglme(sub_table, 'response ~ 1 + (1|patient_ID) + (1|muscle_group)','Distribution','binomial');

    ints(kk) = -thr_model.Coefficients.Estimate;
    ints_ci(kk,1) = -thr_model.Coefficients.Lower;
    ints_ci(kk,2) = -thr_model.Coefficients.Upper;
    ints_p(kk) =thr_model.Coefficients.pValue;
end

figure
hold on
for kk = 1:length(ints)
    plot([kk, kk],[ints_ci(kk,1), ints_ci(kk,2)],'k','LineWidth',5)
    scatter(kk,ints(kk),400,'black','filled');
end

xticks(1:length(ints))
xticklabels(labels)
xlim([0.5,length(unique_contacts)+0.5])
ylabel('Est. Amplitude Threshold (mA)')
xlabel('Contact')
set(gca,'FontSize',14)
% x_sig = linspace(-1,15,100);
% plot(x_sig, sigmoid(sig_coeffs,x_sig),'k')

%% testing thresholds for individual muscles

unique_contacts = mep_sorted;

labels = unique_contacts;

ints = nan(size(unique_contacts));
ints_ci = nan(length(unique_contacts),2);
ints_p = nan(size(unique_contacts));

for kk = 1:length(unique_contacts)
    if kk<length(unique_contacts)
        inds = and(strcmp(mep_table.stim_type,'cathodic'),strcmp(mep_table.cathode,unique_contacts{kk}));
    else
        inds = strcmp(mep_table.stim_type,'bipolar');
    end
    sub_table = mep_table(inds,:);
    
    thr_model = fitglme(sub_table, 'response ~ 1 + (1|patient_ID) + (1|muscle_group)','Distribution','binomial');

    ints(kk) = -thr_model.Coefficients.Estimate;
    ints_ci(kk,1) = -thr_model.Coefficients.Lower;
    ints_ci(kk,2) = -thr_model.Coefficients.Upper;
    ints_p(kk) =thr_model.Coefficients.pValue;
end

figure
hold on
for kk = 1:length(ints)
    plot([kk, kk],[ints_ci(kk,1), ints_ci(kk,2)],'k','LineWidth',5)
    scatter(kk,ints(kk),400,'black','filled');
end

xticks(1:length(ints))
xticklabels(labels)
xlim([0.5,length(unique_contacts)+0.5])
ylabel('Est. Amplitude Threshold (mA)')
xlabel('Contact')
set(gca,'FontSize',14)

%% test nanpad
q = nanpad(ltc_by_muscle);


function M = nanpad(cellArray)
    % padCellArrayWithNaN converts a cell array of vectors into a matrix
    % padded with NaNs to match the length of the longest vector.
    %
    % Input:
    %   cellArray - a 1-D cell array containing vectors of varying lengths
    %
    % Output:
    %   M - a matrix where each column (or row) corresponds to a vector
    %       from the cell array, padded with NaNs

    % Find the length of the longest vector
    maxLength = max(cellfun(@length, cellArray));

    % Number of vectors
    n = numel(cellArray);

    % Preallocate matrix with NaNs
    M = NaN(maxLength, n);

    % Fill in each column with the vector values
    for i = 1:n
        vec = cellArray{i};
        M(1:length(vec), i) = vec(:);  % Ensure column vector
    end
end

function [converted] = convert_contacts(contacts)
    converted = cell(size(contacts));
    for kk = 1:length(contacts)
        if length(contacts{kk}) == 1
            converted{kk} = contacts{kk}{1};
        else
            converted{kk} = contacts{kk}{1}(1);
        end
        
        if contains(converted{kk},'9')
            converted{kk} = replace(converted{kk},'9','1');
        elseif contains(converted{kk},'10')
            converted{kk} = replace(converted{kk},'10','2');
        elseif contains(converted{kk},'11')
            converted{kk} = replace(converted{kk},'11','3');
        elseif contains(converted{kk},'12')
            converted{kk} = replace(converted{kk},'12','4');
        end
        
    end
end