
load('EP_files_040124.mat')
load('EP_042023.mat')

%% 
for kk = 1:length(ep_files)
    if ~isempty(ep_files(kk).emg_bad_channels)
        fprintf('%s\n', ep_files(kk).patient_ID)
        disp(ep_files(kk).emg_bad_channels)
    end
end

%eval_pts = {'ephys025','ephys026','ephys036','ephys038','ephys043','ephys072','ephys081'};
eval_pts = {'ephys025','ephys026', 'ephys036'};

dc_vals = {};
var_vals = {};
types = {};

baseline = [1100,2200];

for kk = 1:length(eval_pts)
    ptind = find(contains({ep_struct(:).patient_ID},eval_pts{kk}),1);
    if strcmp(ep_struct(kk).patient_ID,'ephys025')
        bad_chans = [3,4,8];
    elseif strcmp(ep_struct(kk).patient_ID,'ephys026')
        bad_chans = [3,5,7,8];
    else
        fileind = find(contains({ep_files(:).patient_ID},eval_pts{kk}),1);
        bad_chans = ep_files(fileind).emg_bad_channels;
    end
    
    for kl = 1:length(ep_struct(ptind).mep.emg_labels)
        
        tempmat = ep_struct(ptind).mep.mep_means(:,kl);
        tempmat = double(cell2mat(tempmat(~cellfun(@isempty,tempmat))));
        for km = 1:size(tempmat,1)
            tempmat(km,:) = tempmat(km,:) - median(tempmat(km,end-220:end));
        end
        if isempty(tempmat)
            continue
        end
        if any(bad_chans == kl)
            types = [types; {'bad'}];
        else
            types = [types; {'good'}];
        end
        dc_vals = [dc_vals; {mean(tempmat(:,baseline(1):baseline(2)),2)}];
        var_vals = [var_vals; {var(tempmat(:,baseline(1):baseline(2)),0,2)}];
    end
    
end

empty_inds = cellfun(@isempty, dc_vals);
types(empty_inds) = [];
dc_vals(empty_inds) = [];
var_vals(empty_inds) = [];

disp('Done')

%%
mean_thresh = 10; var_thresh = 50;

dc_good = cell2mat(dc_vals(strcmp(types,'good')));
var_good= cell2mat(var_vals(strcmp(types,'good')));

dc_bad = cell2mat(dc_vals(strcmp(types,'bad')));
var_bad= cell2mat(var_vals(strcmp(types,'bad')));

figure
subplot(2,2,1)
histogram(log10(abs(dc_good)), 25,'FaceColor','blue','Normalization','pdf')
hold on
histogram(log10(abs(dc_bad)), 25,'FaceColor','red','Normalization','pdf')
xline(log10(mean_thresh),'k--')
ylabel('Count (normalized)')
xlabel('Signal mean (log10)')

subplot(2,2,2)
histogram(log10(abs(var_good)), 25,'FaceColor','blue','Normalization','pdf')
hold on
histogram(log10(abs(var_bad)), 25,'FaceColor','red','Normalization','pdf')
xline(log10(var_thresh),'k--')
xlabel('Signal variance (log10)')

subplot(2,2,[3 4])
scatter(log10(abs(dc_good)), log10(var_good), 'b', 'filled')
hold on
scatter(log10(abs(dc_bad)), log10(var_bad), 'r', 'filled')
xline(log10(mean_thresh),'k--')
yline(log10(var_thresh),'k--')
xlabel('Signal mean (log10)')
ylabel('Signal variance (log10)')
legend({'Good channels','Bad channels'},'Location','Southeast')

%% manually check bad channels
tt = linspace(0,100,2200);

trials_plot = [1,1,1];
muscles_plot = [7,5,7];
pts_plot = [14,11,11];

figure
for kk = 1:3
    subplot(3,1,kk); 
    temptrace = ep_struct(pts_plot(kk)).mep.mep_means{trials_plot(kk),muscles_plot(kk)};
    temptrace = temptrace - median(temptrace);
    plot(tt,temptrace); 
end
xlabel('Time (ms)')
% for kk = 1:8
%     subplot(2,4,kk); plot(tt,ep_struct(11).mep.mep_means{2,kk}); ylim([-500,500])
% end


%% plot channels with exclusion criteria

mean_thresh = 5; var_thresh = 50;

pt_ind =10; 

use_thresh = false;
show_marks = false;
show_sm_marks = true;
show_thresh = false; thresh_to_plot = 3;

figure
hold on
col_code = {'g','r','b','k'};
offset_y = 50; offset_x = 110; x_plot = linspace(1,100,2200); 
base_length = 220;

xtick_pos = (1:offset_x:offset_x*length(ep_struct(pt_ind).mep.emg_labels))+offset_x;

n_exc = 0; n_tot = 0;
for kp = 1:length(ep_struct(pt_ind).mep.emg_labels)
    for ko = 1:size(ep_struct(pt_ind).mep.labels_detect,1)
        
        temptrace = ep_struct(pt_ind).mep.mep_means{ko,kp};
        if isempty(temptrace)
            continue
        end
        baseline = median((temptrace(end-base_length:end)));
        temptrace = temptrace - baseline;

        if (abs(mean(temptrace(1100:2200))) > mean_thresh) || (var(temptrace(1100:2200)) > var_thresh)
            exc_flag = true;
            coltemp = 'r';
            n_exc = n_exc + 1;
        else
            exc_flag = false;
            coltemp = 'k';
        end
        

        sd = std(temptrace(end-base_length:end));
        temptrace = (temptrace)/sd;
        
        delay_temp = ep_struct(pt_ind).mep.delays_detect{ko,kp};
        amps_temp = ep_struct(pt_ind).mep.amps_detect_z{ko,kp};
        
%             if (ep_struct(pt_ind).mep.raw.labels_detect(ko,kp)) && (mep_labels_sm{pt_ind}(ko,kp))
%                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{1})
%             elseif (ep_struct(pt_ind).mep.raw.labels_detect(ko,kp)) && ~(mep_labels_sm{pt_ind}(ko,kp)) %false positive; red
%                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{2})
%             elseif ~(ep_struct(pt_ind).mep.raw.labels_detect(ko,kp)) && (mep_labels_sm{pt_ind}(ko,kp)) %false negative; blue
%                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{3})
%             else
%                 plot(x_plot+kp*offset_x,temptrace+ko*offset_y,col_code{4})
%             end
            plot(x_plot+kp*offset_x,temptrace+ko*offset_y,coltemp)

            if ~isempty(delay_temp) && show_marks
                for kq = 1:length(delay_temp)
                    delay_ind = find(x_plot>delay_temp(kq),1);
                    plot(x_plot(delay_ind) +kp*offset_x, temptrace(delay_ind) + ko*offset_y, 'k*', 'MarkerSize', 20);
                end
            end

            if show_thresh
                plot([x_plot(1), x_plot(end)]+kp*offset_x,[-thresh_to_plot,-thresh_to_plot ]+ko*offset_y,'k--')
                plot([x_plot(1), x_plot(end)]+kp*offset_x,[thresh_to_plot,thresh_to_plot ]+ko*offset_y,'k--')
            end

        n_tot = n_tot + 1;
    end
%     
%     if show_thresh
%         ltc_to_plot = latencies(emg_mapping_all{pt_ind}(kp),:);
%         plot([x_plot(1), x_plot(end)]+kp*offset_x,[mean(temptrace) - thresh_to_plot,mean(temptrace) - thresh_to_plot ]+ko*offset_y,'k--')
%         plot([x_plot(1), x_plot(end)]+kp*offset_x,[mean(temptrace) + thresh_to_plot,mean(temptrace) + thresh_to_plot ]+ko*offset_y,'k--')
%         xline(kp*offset_x + ltc_to_plot(1) - latency_tolerance, 'k--'); xline(kp*offset_x + ltc_to_plot(2) + latency_tolerance, 'k--'); 
%         %fprintf('Accuracy for %s: %.3f           Latency: %.3f - %.3f ms\n',ep_struct(pt_ind).mep.emg_labels{kp},acc_by_pt{pt_ind}(kp), ltc_to_plot(1)-latency_tolerance, ltc_to_plot(2)+latency_tolerance)
%     else
%         %fprintf('Accuracy for %s: %.3f\n',ep_struct(pt_ind).mep.emg_labels{kp},acc_by_pt{pt_ind}(kp))
%     end
    
end

fprintf('%.2f settings excluded\n', 100*n_exc/n_tot)

set(gca,'xTick',xtick_pos)
set(gca,'XtickLabel',ep_struct(pt_ind).mep.emg_labels)

