%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

clc
clear
close all

%% Build sensitivity table
files   = {'SA_results_Diesel.mat','SA_results_Jet fuel.mat','SA_results_Gasoline.mat'};
names   = {'Diesel', 'Jet Fuel', 'Gasoline'};
outFile = 'SensitivityTables.xlsx';

colNames = {'Threshold','RMSE_H','RMSE_DH','n','r_1','r_2','r_3','r_4','r_5','r_6', ...
'A_1','A_2','A_3','A_4','B_1','B_2','B_3','B_4','alpha','beta'};
nCols  = numel(colNames);
R6_COL = 10;
delta  = 0.1;

if exist(outFile, 'file')
    delete(outFile);
end

boldRanges = cell(1,numel(files));

for f = 1:numel(files)
    S    = load(files{f}, 'sens_results');
    sens = S.sens_results(:);
    nT   = numel(sens);

    data = nan(nT, nCols);
    data(:,1) = [sens.threshold]';
    data(:,2) = arrayfun(@(s) s.R2sysV(s.pareto_point_id), sens)';
    data(:,3) = arrayfun(@(s) s.R2fitV(s.pareto_point_id), sens)';
    data(:,4) = arrayfun(@(s) s.n_sys(1), sens);
for k = 1:nT
        R6full = reshape(sens(k).R_sys(1:6), 1, []);
        L6full = reshape(double(sens(k).L(1:6)), 1, []);
        data(k,5:10) = R6full ./ L6full;
end
    pipe6Fixed = arrayfun(@(s) s.R_fit(6) == 1 && s.n_fit(6) == 0, sens);
    data(pipe6Fixed, R6_COL) = NaN;
    data(:,11:14) = cell2mat(arrayfun(@(s) reshape(s.A_sys(1:4),1,[]), sens, 'UniformOutput', false));
    data(:,15:18) = cell2mat(arrayfun(@(s) reshape(s.B_sys(1:4),1,[]), sens, 'UniformOutput', false));
    data(:,19)    = [sens.alpha_sys]';
    data(:,20)    = [sens.beta_sys]';

    stat_cols = 2:nCols;
    mn      = min(data(:,stat_cols), [], 1, 'omitnan');
    mx      = max(data(:,stat_cols), [], 1, 'omitnan');
    avr     = mean(data(:,stat_cols), 1, 'omitnan');
    sd      = std(data(:,stat_cols), 0, 1, 'omitnan');
    cov_pct = 100 * sd ./ abs(avr);

    ref_row = data(abs(data(:,1)-0.95)<1e-9, :);
    lower   = ref_row(stat_cols) * (1-delta);
    upper   = ref_row(stat_cols) * (1+delta);
    lower(end) = ref_row(end) * (1+delta);
    upper(end) = ref_row(end) * (1-delta);

    withinTol = false(nT, numel(stat_cols));
    for c = 1:numel(stat_cols)
        v = data(:,stat_cols(c));
        withinTol(:,c) = v >= lower(c) & v <= upper(c);
    end

    nRows = 2 + nT + 5 + 7;
    sheet = cell(nRows, nCols);
    tableNum = f + 1;   % Diesel=S2, Jet Fuel=S3, Gasoline=S4
    sheet{1,1} = sprintf('Table S%d: Sensitivity analysis for filtering percentile (%s)', tableNum, names{f});
    sheet(2,:) = colNames;
    sheet(3:2+nT, :) = num2cell(data);
    nanRows = isnan(data(:,R6_COL));
    sheet(find(nanRows)+2, R6_COL) = {[]};

    stat_labels = {'Min','Max','Avr','Std','CoV (%)'};
    stat_vals   = [mn; mx; avr; sd; cov_pct];
    for k = 1:5
        r = 2 + nT + k;
        sheet{r,1} = stat_labels{k};
        sheet(r, stat_cols) = num2cell(stat_vals(k,:));
    end

    r0 = 2 + nT + 5;
    sheet{r0+1,2} = 'Notes: values in bold are within 10% of the value obtained in 95th percentile';

    writecell(sheet, outFile, 'Sheet', names{f}, 'Range', 'B1');
    fprintf('%s: wrote %d data rows + stats + tolerance block to sheet "%s" of %s.\n', names{f}, nT, names{f}, outFile);

    boldRanges{f} = withinTol;
end

excel = actxserver('Excel.Application');
cleanupObj = onCleanup(@() cleanupExcel(excel));

excel.DisplayAlerts = false;
excel.Visible = false;
wb = excel.Workbooks.Open(fullfile(pwd, outFile));

colLetter = @(c) [char(64+floor((c-1)/26)*(c>26)) char(65+mod(c-1,26))];

for f = 1:numel(files)
    fprintf('Formatting sheet: %s\n', names{f});
    ws = wb.Sheets.Item(names{f});
    withinTol = boldRanges{f};
    [nT, nC] = size(withinTol);
    for r = 1:nT
        for c = 1:nC
            if withinTol(r,c)
                excelRow = 2 + r;
                excelCol = 2 + c;
                addr = sprintf('%s%d', colLetter(excelCol), excelRow);
                try
                    ws.Range(addr).Font.Bold = true;
                catch ME
                    fprintf('FAILED at %s on sheet %s: %s\n', addr, names{f}, ME.message);
                end
            end
        end
    end
end
wb.Save;
wb.Close(false);
fprintf('Applied bold formatting to %s.\n', outFile);

function cleanupExcel(excel)
    try
        excel.Quit;
    catch
    end
    delete(excel);
end