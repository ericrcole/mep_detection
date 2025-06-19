%analyze_mep_correlation.m
%
%run "plot_mEP_results_v3" initialization block first
%analyzes relative correlation strength between mEP waveforms within
%a specific subject/muscle, and between, to see if personalized mEP
%templates can improve detection

corr_within_muscle = cell(9,length(emg_mapping_labels));
corr_across_muscle = cell(size(emg_mapping_labels));
means_by_muscle = cell(size(emg_mapping_labels));
means_by_muscle_full = cell(size(emg_mapping_labels));
ptinds_by_muscle = cell(size(emg_mapping_labels));

t_plot = linspace(0,100,2200);

for kk = 1:length(emg_mapping_sm)
    for kl = 1:length(emg_mapping_sm{kk})
        tempmat = cell2mat(ep_struct(kk).mep.raw.mep_means(mep_labels_sm{kk}(:,kl),kl));
        m_ind = emg_mapping_all{kk}(kl);

        ltctemp = latencies(m_ind,:) + [-latency_tolerance, latency_tolerance];
        ltc1 = find(t_plot>=ltctemp(1),1); ltc2 = find(t_plot>=ltctemp(2),1); 
        
        if ~isempty(tempmat)
            means_by_muscle{m_ind} = [means_by_muscle{m_ind}; tempmat(:,ltc1:ltc2)];
            means_by_muscle_full{m_ind} = [means_by_muscle_full{m_ind}; tempmat];
            ptinds_by_muscle{m_ind} = [ptinds_by_muscle{m_ind}; kk*ones(size(tempmat,1),1)];

            corr_within_muscle{kk,m_ind} = [corr_within_muscle{kk,m_ind}; mean(mean(corrcoef(tempmat')))];
        end
%     pt_accs = [pt_accs; mean(ep_struct(kk).mep.labels_detect == mep_labels_sm{kk})];
%     pt_accs_co = [pt_accs_co; mean(ep_struct(pt_ind).mep.labels_detect == mep_labels_co(pt_inds_co(kk)).labels)];
%     
%     tempacc = zeros(length(emg_mapping_sm{kk}),length(ampz_sweep));
%     for kl = 1:length(ampz_sweep)
%         if iscell(ep_struct(kk).mep.raw.amps_detect_z)
%             tempmat = cellfun(@nanmax, ep_struct(kk).mep.raw.amps_detect_z, 'UniformOutput', false);
%             tempmat(cellfun(@isempty,tempmat)) = {0};
%             temp_detect = cell2mat(tempmat)>ampz_sweep(kl);
%         else
%             temp_detect = ep_struct(kk).mep.raw.amps_detect_z>ampz_sweep(kl);
%         end
%         tempacc(:,kl) = mean(temp_detect == mep_labels_sm{kk});
%         %disp(mean(temp_detect == mep_labels_sm{kk}))
%         
% %         if ampz_sweep(kl)==2
% %             acc_sweep_by_pt{kk} = (temp_detect == mep_labels_sm{kk});
% %         end
%     end
%     for km = 1:length(emg_mapping_sm{kk})
%         acc_by_muscle{emg_mapping_all{kk}(km)} = [acc_by_muscle{emg_mapping_all{kk}(km)};  mean(ep_struct(kk).mep.raw.labels_detect(:,km)==mep_labels_sm{kk}(:,km))];
%     end
%     %acc_all = cat(3,acc_sweep, tempmat);
%     acc_sweep_by_pt{kk} = tempacc;
%     
%     acc_by_pt{kk} = mean(ep_struct(kk).mep.raw.labels_detect == mep_labels_sm{kk});

    %pt_accs{emg_mapping_all{kk}(km)} = [pt_accs{emg_mapping_all{kk}(km)}; mean(ep_struct(kk).mep.labels_detect == mep_labels_sm{kk})];
    end
end

%% mEP waveform correlation matrices
figure
ctr = 1;
for kk = [1 2 3 4 5 6 7 8]
    subplot(2,4,ctr)
    imagesc(corrcoef(means_by_muscle{kk}'))
    temp = find(diff(ptinds_by_muscle{kk})>0)+0.5;
    hold on
    
    for ind = 1:length(temp)
        xline(temp(ind), 'k--')
        yline(temp(ind), 'k--')
    end
    title(emg_mapping_labels{kk})
    colorbar
    ctr = ctr+1;
end

% for kk = 1:8
%     
% 
% end

%%

pt_list = [6,2,9];
m_list = [2,2,2];
n_curves = 4;

tplot = linspace(0,100,2200);

figure
for kk = 1:length(pt_list)
    subplot(2,3,kk)
    tempmat = cell2mat(ep_struct(pt_list(kk)).mep.raw.mep_means(mep_labels_sm{pt_list(kk)}(:,m_list(kk)),m_list(kk)));
    %disp(size(tempmat))
    %trials = randi(size(tempmat,1),[n_curves 1]);
    trials = randperm(size(tempmat,1),n_curves);
    for kl = 1:length(trials)
        temp = tempmat(trials(kl),:);
        plot(tplot,temp - median(temp))
        hold on
    end
    ylim([-50 50])
    xlim([0 50])

    xlabel('time (ms)')
    ylabel('voltage (\muV)')
    title(sprintf('%s; %s',ep_struct(pt_list(kk)).patient_ID, ep_struct(pt_list(kk)).mep.emg_labels{m_list(kl)}))
end

pt_list = [6,7,8];
m_list = [6,6,5];
n_curves = 3;

tplot = linspace(0,100,2200);

for kk = 1:length(pt_list)
    subplot(2,3,kk+3)
    tempmat = cell2mat(ep_struct(pt_list(kk)).mep.raw.mep_means(mep_labels_sm{pt_list(kk)}(:,m_list(kk)),m_list(kk)));
    %disp(size(tempmat))
    %trials = randi(size(tempmat,1),[n_curves 1]);
    trials = randperm(size(tempmat,1),n_curves);
    for kl = 1:length(trials)
        temp = tempmat(trials(kl),:);
        plot(tplot,temp - median(temp))
        hold on
    end
    ylim([-50 50])
    xlim([0 50])

    xlabel('time (ms)')
    ylabel('voltage (\muV)')
    title(sprintf('%s; %s',ep_struct(pt_list(kk)).patient_ID, ep_struct(pt_list(kk)).mep.emg_labels{m_list(kl)}))

end

%% figure 3c: cross-patient correlation matrices
save_fig = true;

figure('Position',[100 100 800 450])

pos{1} = [0.065 0.735 0.2 0.2];
pos{2} = [0.065 0.41 0.2 0.2];

pos{3} = [0.295 0.735 0.2 0.2];
pos{4} = [0.295 0.41 0.2 0.2];

pos{5} = [0.065 0.1 0.425 0.2];
%pos4 = [0.1 0.66 0.2 0.3];

pos{6} = [0.54 0.15 0.19 0.325];
pos{7} = [0.7725 0.15 0.19 0.325];
pos{8} = [0.54 0.59 0.19 0.325];
pos{9} = [0.7725 0.59 0.19 0.325];

annotation('textbox',[.002 .98 .05 .05],'String','a','EdgeColor','none','FontSize',36)
annotation('textbox',[.002 .66 .05 .05],'String','b','EdgeColor','none','FontSize',36)
annotation('textbox',[.5 .98 .05 .05],'String','c','EdgeColor','none','FontSize',36)
annotation('textbox',[.002 .35 .05 .05],'String','d','EdgeColor','none','FontSize',36)

% subplot('Position',pos1)
% subplot('Position',pos2)
% subplot('Position',pos3)
% subplot('Position',pos4)
% subplot('Position',pos5)
% 
% subplot('Position',pos6)
% subplot('Position',pos7)
% subplot('Position',pos8)
% subplot('Position',pos9)

pt_list = [6,6];
m_list = [2,6];

n_curves = 4;

trials1 = [2,4,5,8];
trials2 = [21,5,15]; %15, 21, 5 good

cols = {'r','b','k','m'};
%yoffset = 25;

%figure
for kk = 1:length(pt_list)
    if kk == 1
        %subplot(2,3,1)
        subplot('Position',pos{1})
    else
        subplot('Position',pos{3})
    end
    tempmat = cell2mat(ep_struct(pt_list(kk)).mep.raw.mep_means(mep_labels_sm{pt_list(kk)}(:,m_list(kk)),m_list(kk)));
    %disp(size(tempmat))
    %trials = randi(size(tempmat,1),[n_curves 1]);
    disp(ep_struct(pt_list(kk)).mep.emg_labels{m_list(kk)});
    if kk == 1
    trials = trials1;
    elseif kk == 2
    trials = trials2;
    else
    trials = randperm(size(tempmat,1),n_curves);
    end
    for kl = 1:length(trials)
        temp = tempmat(trials(kl),:);
        if (kk == 2) && (kl == 1)
            temp = temp - median(temp) -2;
            plot(tplot,temp, cols{kl})
        else
            plot(tplot,temp - median(temp), cols{kl})
        end
       
        hold on
    end
    ylim([-50 50])
    %ylim([-50 150])
    xlim([0 50])

    xlabel('Time (ms)','FontSize',12)
    if kk == 1
    ylabel('Voltage (\muV)','FontSize',14)
    end
    if kk == 1
        title('ECR','FontSize',14)
    else
        title('Orb. Oris','FontSize',14)
    end
    %title(ep_struct(pt_list(kk)).mep.emg_labels{m_list(kk)})
    %title(sprintf('%s; %s',ep_struct(pt_list(kk)).patient_ID, ep_struct(pt_list(kk)).mep.emg_labels{m_list(kl)}))
end

m_list = [2,8];

n_curves = 3;

trials1 = [1,3,6,11];
trials2 = [1,7,26];

cols = {'b','r','k','m'};

for kk = 1:length(m_list)
    if kk == 1
        subplot('Position',pos{2})
        trials = trials1;
    else
        trials = trials2;
        subplot('Position',pos{4})
    end

    %tempmat = cell2mat(ep_struct(pt_list(kk)).mep.raw.mep_means(mep_labels_sm{pt_list(kk)}(:,m_list(kk)),m_list(kk)));
    tempmat = means_by_muscle_full{m_list(kk)}(trials,:);

    for kl = 1:length(trials)
        temp = tempmat(kl,:);
        plot(tplot,temp - median(temp), cols{kl})
        hold on
    end
    ylim([-50 50])
    xlim([0 50])

    xlabel('Time (ms)','FontSize',12)
    if kk == 1
    ylabel('Voltage (\muV)','FontSize',14)
    end

    if kk == 1
        title('ECR','FontSize',14)
    else
        title('Orb. Oris','FontSize',14)
    end
    %title(ep_struct(pt_list(kk)).mep.emg_labels{m_list(kk)})
    %title(sprintf('%s; %s',ep_struct(pt_list(kk)).patient_ID, ep_struct(pt_list(kk)).mep.emg_labels{m_list(kl)}))
end


%splist = [2,3,5,6];
splist = [6 7 8 9];
cmat_list = [3,4,5,8];

inds = [];

ctr=1;
for kk = cmat_list
    %subplot(2,3,splist(ctr))
    subplot('Position', pos{splist(ctr)})
    imagesc(corrcoef(means_by_muscle{kk}'))
    temp = find(diff(ptinds_by_muscle{kk})>0)+0.5;
    hold on
    
    maxind = size(corrcoef(means_by_muscle{kk}'),2);
    prevind = 0;

    inds = [];
    for ind = 1:length(temp)
        xline(temp(ind), 'k-', 'LineWidth', 2)
        yline(temp(ind), 'k-', 'LineWidth', 2)

        inds = [inds; (temp(ind) + prevind)/2];
        prevind = temp(ind);
    end
    inds = [inds; (maxind + prevind)/2];
    
    ctr = ctr+1;
%     if (kk == 3) || (kk==6)
%         colorbar
%     end
    inds = unique(inds);
    % inds = inds(1:end-1)+diff(inds);

    xticks(inds);
    labels = {};
    for kl = 1:length(inds)
        labels = [labels; {sprintf('pt. %d', kl)}];
    end
    xticklabels(labels)
    xtickangle(90)
    
    yticks(inds);
    yticklabels(labels)
    set(gca,'FontSize',10)
    title(emg_mapping_labels{kk},'FontSize',14)
    %ytickangle(45)
end

% bar graph and statistical test comparing within vs. across-muscle correlation

muscle_labels = cell(size(means_by_muscle));
within_pt_corr = cell(size(means_by_muscle));
across_pt_corr = cell(size(means_by_muscle));

ctr=1;
for kk = 1:length(means_by_muscle)

    muscle_labels{kk} = emg_mapping_labels{kk};
    corrmat = corrcoef(means_by_muscle{kk}');
    temp = find(diff(ptinds_by_muscle{kk})>0)+1;
    if isempty(temp)
        continue
    end
    %temp = [temp; size(corrmat,2)-1];
    temp(end) = temp(end)-1;
    
    for kl = 1:(length(temp)-1)
        for km = 1:(length(temp)-1)
            ind1 = temp(kl);
            ind2 = temp(km);
            
            submat = corrmat(temp(kl):temp(kl)+1, temp(km):temp(km)+1);
            if kl == km
               N = size(submat,1);
               idx = eye(N,N);
               Y = (1-idx).*submat;
               Z = submat(~idx);

               within_pt_corr{kk} = [within_pt_corr{kk}; mean(Z)];
            else
               across_pt_corr{kk} = [across_pt_corr{kk}; mean(submat(:))];
            end
        end
    end

end

h = gca;
c = colorbar(h,'Position',[0.97 0.375 0.0125 0.325]);
c.Ruler.TickLabelRotation=270;
clim([-1 1]);


across_pt_all = [];
within_pt_all = [];

subplot('Position',pos{5})
for kk = 1:length(across_pt_corr)
    temp = across_pt_corr{kk};
    across_pt_all = [across_pt_all; temp];
    scatter(ones(size(temp)) + (randn(size(temp)))/10, temp, 'filled')
    hold on

    temp = within_pt_corr{kk};
    within_pt_all = [within_pt_all; temp];
    scatter(ones(size(temp))+1 + (randn(size(temp)))/10, temp, 'filled')
    hold on
end

boxchart(ones(size(across_pt_all)), across_pt_all, 'BoxFaceColor',[0.5 0.5 0.5]);
boxchart(ones(size(within_pt_all))+1, within_pt_all, 'BoxFaceColor',[0.5 0.5 0.5]);

ylim([-1,1.08])
%set(gca,'FontSize',16)

xticks(1:2)
%xticklabels({'Across-Patient', 'Within-Patient'},'FontSize',14)
set(gca,'XTickLabel',{'Across-Patient', 'Within-Patient'},'fontsize',12)
ylabel('Correlation','FontSize',14)


%[h,p] = ttest2(across_pt_all, within_pt_all);
[p,h] = ranksum(across_pt_all, within_pt_all);
fprintf('across vs. within-patient: p = %0.4f\n', p)

if save_fig
    cd('/Users/ercole/Documents/mep_paper/figs/')
    exportgraphics(gcf,'fig4_corr_v2.png','Resolution',300)
end

%% Develop mep detection with template-based convolution

pt_ind = 6;
emg_ind = 5;
t_plot = linspace(0,100,2200);

ncorr = 3;

kk = pt_ind; kl = emg_ind;
tempmat = cell2mat(ep_struct(kk).mep.raw.mep_means(mep_labels_sm{kk}(:,kl),kl));
m_ind = emg_mapping_all{kk}(kl);

ltctemp = latencies(m_ind,:) + [-latency_tolerance, latency_tolerance];
ltc1 = find(t_plot>=ltctemp(1),1); ltc2 = find(t_plot>=ltctemp(2),1); 

if ~isempty(tempmat)
    corrmat = corrcoef(tempmat(:,ltc1:ltc2)');
    
    [vals, inds] = maxk(mean(corrmat),3);
    conv_template = mean(tempmat(inds,ltc1:ltc2));
    conv_template = conv_template - movmin(conv_template,75);
%     means_by_muscle{m_tempmatind} = [means_by_muscle{m_ind}; tempmat(:,ltc1:ltc2)];
%     ptinds_by_muscle{m_ind} = [ptinds_by_muscle{m_ind}; kk*ones(size(tempmat,1),1)];
% 
%     corr_within_muscle{kk,m_ind} = [corr_within_muscle{kk,m_ind}; mean(mean(corrcoef(tempmat')))];
end

figure
subplot(1,2,1)
histogram(corrmat(3,:),25,'FaceColor','b')
hold on
histogram(corrmat(5,:),25,'FaceColor','r')
histogram(corrmat(6,:),25,'FaceColor','y')

subplot(1,2,2)
plot(tempmat(inds,ltc1:ltc2)')
hold on
plot(conv_template, 'k--')
ylim([-100,100])
xlim([-100,500])

%%
figure; offset_y = 50;
for kk = 1:size(tempmat,1)
    temptrace = tempmat(kk,:);
    temptrace2 = conv(temptrace, conv_template, 'same');

% %     temptrace2 = temptrace; 
% %     temptrace2(ltc1:ltc2) = temptrace2(ltc1:ltc2).*conv_template;

    baseline = median((temptrace(end-base_length:end)));
    sd = std(temptrace(end-base_length:end));
    temptrace = (temptrace -baseline)/sd;

    baseline = median((temptrace2(end-base_length:end)));
    sd = std(temptrace2(end-base_length:end));
    temptrace2 = (temptrace2 -baseline)/sd;

    plot(t_plot, temptrace+offset_y*kk, 'b')
    hold on
    plot(t_plot, temptrace2+offset_y*kk, 'k--')
end

%%
kk = pt_ind; kl = emg_ind;
tempmat_all = ep_struct(kk).mep.raw.mep_means(:,kl);
logvals = cellfun(@isempty, tempmat_all);
tempmat_all(cellfun(@isempty, tempmat_all)) = [];
tempmat_all = cell2mat(tempmat_all);

corrvals = corrcoef([conv_template; tempmat_all(:,ltc1:ltc2)]');
corrvals = corrvals(1,2:end);

logtemp = mep_labels_sm{kk}(:,kl);
logtemp(logvals) = [];

figure
histogram(corrvals(logtemp == 1),25,'FaceColor','b')
hold on
histogram(corrvals(logtemp == 0),25,'FaceColor','r')
