%% fig. 5ab: does stim amplitude affect mEP amplitude?
% addpath('/Users/ERCOLE/Documents/Research/Repos/densityScatterChart')
% pct variability between muscles?
%amp_model = fitglme(mep_table, 'mep_ampl ~ 1 + dbs_ampl*muscle_adj + (dbs_ampl*muscle_adj|patient_ID) + (dbs_ampl*muscle_adj|contact)');
% amp_model = fitglme(mep_table, 'mep_ampl ~ dbs_ampl + (dbs_ampl|patient_ID) + (dbs_ampl|muscle_adj)');
% 
% disp(amp_model)

cols = [transpose(linspace(0,1,6)), [0;0;0;0;0;0], [0;0;0;0;0;0]];

mep_table = readtable('mep_table_opt.csv')
%mep_table(mep_table.dbs_ampl > 5,:) = [];
mep_table(mep_table.mep_ampl > 400,:) = [];

figure
subplot(2,1,1)
%scatter(mep_table.dbs_ampl, mep_table.mep_ampl,'filled')
for kk = 1:6
    temp = mep_table.mep_ampl(mep_table.dbs_ampl == kk);
    boxchart(kk*ones(size(temp)),temp, 'BoxFaceColor',cols(kk,:),'MarkerColor',cols(kk,:));
    hold on
end

xlabel('DBS amplitude (mA)')
ylabel('mEP amplitude (\muV)')
set(gca, 'FontSize', 12)
%densityScatterChart(mep_table.dbs_ampl, mep_table.mep_ampl)

% subplot(2,1,1)
% scatter(mep_table.dbs_ampl + randn(size(mep_table.dbs_ampl))*.05, mep_table.mep_ampl,'filled')


%amp_model = fitlme(mep_table, 'mep_ampl ~ dbs_ampl + (dbs_ampl|patient_ID) + (dbs_ampl|cathode)');
amp_model = fitlme(mep_table, 'mep_ampl ~ dbs_ampl + (dbs_ampl|patient_ID) + (dbs_ampl|muscle_adj) ');
%disp(amp_model)

fprintf('Peak amp. vs. DBS amp: %.4e\n',amp_model.Coefficients.pValue(2))

lcoeffs = amp_model.Coefficients.Estimate;

% hold on
% plot([0, 15], lcoeffs(1)+lcoeffs(2)*[0,15], 'k--')

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
%scatter(mep_table.dbs_ampl + randn(size(mep_table.dbs_ampl))*.05, mep_table.n_response,'filled')
for kk = 1:6
    temp = temptable.n_response(temptable.dbs_ampl == kk);
    boxchart(kk*ones(size(temp)),temp,'BoxFaceColor',cols(kk,:),'MarkerColor',cols(kk,:));
    hold on
end

amp_model = fitlme(temptable, 'n_response ~ dbs_ampl + (dbs_ampl|patient_ID)+ (dbs_ampl|muscle_adj) ');
%disp(amp_model)

fprintf('Response pct. vs. DBS amp: %.4e\n',amp_model.Coefficients.pValue(2))

lcoeffs = amp_model.Coefficients.Estimate;

xlabel('DBS amplitude (mA)')
ylabel('Proportion of muscles')
set(gca, 'FontSize', 12)
