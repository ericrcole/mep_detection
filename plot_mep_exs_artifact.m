
pt_struct = ep_struct;
t_plot = linspace(0,100,2200);

emg_chans = [1 6 2]; %cloth
%emg_chans = [1 6 7]; %dual;
%emg_chans = [1 2 3 4];
%emg_chans = [2,3];

figure
%sp_arr = [6,8];
x_offset = 110;

cols = {'r','b','k','m','g'};

disp(pt_struct.mep.emg_labels)

normalize = false;

for kk = 1:length(pt_struct.stim_settings.param_strings)
    
    for kl = 1:length(emg_chans)
        temptrace = pt_struct.mep.mep_means{kk,emg_chans(kl)};
        temptrace = temptrace- nanmedian(temptrace);
        if normalize
            temptrace = temptrace/nanstd(temptrace);
        end
        plot(t_plot+(x_offset*(kk-1)),temptrace,cols{kl})
        hold on
    end
    
    
    
end
xticks((0:x_offset:(x_offset*(kk-1)))+50)
xticklabels(pt_struct.stim_settings.param_strings)
xtickangle(10)
ylim([-50 50])

legend(pt_struct.mep.emg_labels(emg_chans))
