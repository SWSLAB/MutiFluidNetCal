%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

clc
clear
close all
%% Build Sensitivity Table Cross validation

% ---- CONFIG: edit fluid files/names/output here ------------------------
filesCV   = {'CV_results_Diesel.mat', 'CV_results_Jet fuel.mat', 'CV_results_Gasoline.mat'};
namesCV   = {'Diesel', 'Jet Fuel', 'Gasoline'};
outFileCV = 'CrossValidationTables.xlsx';
delta     = 0.1;
% -------------------------------------------------------------------------
if exist(outFileCV, 'file')
    try
        delete(outFileCV);
    catch ME
        error('Could not delete %s — it may be open in Excel or another program. Close it and rerun. (%s)', outFileCV, ME.message);
    end
end

boldRangesCV = cell(1,numel(filesCV));

for f = 1:numel(filesCV)
    S  = load(filesCV{f}, 'CV_results');
    cv = S.CV_results;
    T  = cv.cv_table;

    pipe6Fixed = (cv.R_fit(6) == 1) && (cv.n_fit(6) == 0);
    if pipe6Fixed
        T.r_6 = nan(height(T),1);
    end

    colNamesCV = T.Properties.VariableNames;
    nColsCV    = numel(colNamesCV);
    nFolds     = height(T);
    data       = table2array(T);

    stat_cols = 2:nColsCV;
    mn      = min(data(:,stat_cols), [], 1, 'omitnan');
    mx      = max(data(:,stat_cols), [], 1, 'omitnan');
    avr     = mean(data(:,stat_cols), 1, 'omitnan');
    sd      = std(data(:,stat_cols), 0, 1, 'omitnan');
    cov_pct = 100 * sd ./ abs(avr);

    lower = avr * (1-delta);
    upper = avr * (1+delta);
    lower(end) = avr(end) * (1+delta);
    upper(end) = avr(end) * (1-delta);

    withinTol = false(nFolds, numel(stat_cols));
    for c = 1:numel(stat_cols)
        v = data(:,stat_cols(c));
        withinTol(:,c) = v >= lower(c) & v <= upper(c);
    end

    nRowsCV = 2 + nFolds + 5 + 1;
    sheet = cell(nRowsCV, nColsCV);
    tableNum = f + 7;   % Diesel=S8, Jet Fuel=S9, Gasoline=S10
    sheet{1,1} = sprintf('Table S%d: Cross-validation analysis (%s)', tableNum, namesCV{f});
    sheet(2,:) = colNamesCV;
    sheet(3:2+nFolds, :) = num2cell(data);
    if pipe6Fixed
        r6col = find(strcmp(colNamesCV,'r_6'));
        sheet(3:2+nFolds, r6col) = {[]};
    end

    stat_labels = {'Min','Max','Avr','Std','CoV (%)'};
    stat_vals   = [mn; mx; avr; sd; cov_pct];
    for k = 1:5
        r = 2 + nFolds + k;
        sheet{r,1} = stat_labels{k};
        sheet(r, stat_cols) = num2cell(stat_vals(k,:));
    end

    r0 = 2 + nFolds + 5;
        sheet{r0+1,1} = 'Notes: values in bold are within 10% of the average value';

    writecell(sheet, outFileCV, 'Sheet', namesCV{f}, 'Range', 'B1');
    fprintf('%s: wrote %d fold rows + stats to sheet "%s" of %s.\n', namesCV{f}, nFolds, namesCV{f}, outFileCV);

    boldRangesCV{f} = withinTol;
end

excel = actxserver('Excel.Application');
cleanupObj = onCleanup(@() cleanupExcel(excel));

excel.DisplayAlerts = false;
excel.Visible = false;
wb = excel.Workbooks.Open(fullfile(pwd, outFileCV));

colLetter = @(c) [char(64+floor((c-1)/26)*(c>26)) char(65+mod(c-1,26))];

for f = 1:numel(filesCV)
    fprintf('Formatting sheet: %s\n', namesCV{f});
    ws = wb.Sheets.Item(namesCV{f});
    withinTol = boldRangesCV{f};
    [nFolds, nC] = size(withinTol);
    for r = 1:nFolds
        for c = 1:nC
            if withinTol(r,c)
                excelRow = 2 + r;
                excelCol = 2 + c;
                addr = sprintf('%s%d', colLetter(excelCol), excelRow);
                try
                    ws.Range(addr).Font.Bold = true;
                catch ME
                    fprintf('FAILED at %s on sheet %s: %s\n', addr, namesCV{f}, ME.message);
                end
            end
        end
    end
end
wb.Save;
wb.Close(false);
fprintf('Applied bold formatting to %s.\n', outFileCV);

function cleanupExcel(excel)
    try
        excel.Quit;
    catch
    end
    delete(excel);
end