%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

clc
clear
close all

%% Build sensitivity table multistart

% ---- CONFIG: edit fluid files/names/output here ------------------------
filesMS   = {'MS_results_Diesel.mat', 'MS_results_Jet fuel.mat', 'MS_results_Gasoline.mat'};
namesMS   = {'Diesel', 'Jet Fuel', 'Gasoline'};
outFileMS = 'MultiStartTables.xlsx';
delta     = 0.1;
% -------------------------------------------------------------------------

colNamesMS = {'Run','RMSE_H','RMSE_DH','n','r_1','r_2','r_3','r_4','r_5','r_6', ...
              'A_1','A_2','A_3','A_4','B_1','B_2','B_3','B_4','alpha','beta'};
nColsMS   = numel(colNamesMS);
R6_COL_MS = 10;

if exist(outFileMS, 'file')
    try
        delete(outFileMS);
    catch ME
        error('Could not delete %s — it may be open in Excel or another program. Close it and rerun. (%s)', outFileMS, ME.message);
    end
end

boldRangesMS = cell(1,numel(filesMS));

for f = 1:numel(filesMS)
    S  = load(filesMS{f}, 'fluid_results');
    fr = S.fluid_results;

    pid   = fr.pareto_point_id;
    nRuns = size(fr.stab_xsol, 2);

    pipe6Fixed = (fr.R_fit(6) == 1) && (fr.n_fit(6) == 0);

    data = nan(nRuns, nColsMS);
    data(:,1) = (1:nRuns)';
    for k = 1:nRuns
        x = fr.stab_xsol{pid, k};
        data(k,2)     = fr.stab_R2sys(pid, k);
        data(k,3)     = fr.stab_R2fit(pid, k);
        data(k,4)     = x(24);
        data(k,5:10)  = x(1:6);
        data(k,11:14) = x(7:10);
        data(k,15:18) = x(11:14);
        data(k,19)    = x(22);
        data(k,20)    = x(23);
    end
    if pipe6Fixed
        data(:,R6_COL_MS) = NaN;
    end

    stat_cols = 2:nColsMS;
    mn      = min(data(:,stat_cols), [], 1, 'omitnan');
    mx      = max(data(:,stat_cols), [], 1, 'omitnan');
    avr     = mean(data(:,stat_cols), 1, 'omitnan');
    sd      = std(data(:,stat_cols), 0, 1, 'omitnan');
    cov_pct = 100 * sd ./ abs(avr);

    lower = avr * (1-delta);
    upper = avr * (1+delta);
    lower(end) = avr(end) * (1+delta);
    upper(end) = avr(end) * (1-delta);

    withinTol = false(nRuns, numel(stat_cols));
    for c = 1:numel(stat_cols)
        v = data(:,stat_cols(c));
        withinTol(:,c) = v >= lower(c) & v <= upper(c);
    end

    nRowsMS = 2 + nRuns + 5 + 1;
    sheet = cell(nRowsMS, nColsMS);
    tableNum = f + 4;   % Diesel=S5, Jet Fuel=S6, Gasoline=S7
    sheet{1,1} = sprintf('Table S%d: Multi-start stability analysis (%s)', tableNum, namesMS{f});
    sheet(2,:) = colNamesMS;
    sheet(3:2+nRuns, :) = num2cell(data);
    nanRows = isnan(data(:,R6_COL_MS));
    sheet(find(nanRows)+2, R6_COL_MS) = {[]};

    stat_labels = {'Min','Max','Avr','Std','CoV (%)'};
    stat_vals   = [mn; mx; avr; sd; cov_pct];
    for k = 1:5
        r = 2 + nRuns + k;
        sheet{r,1} = stat_labels{k};
        sheet(r, stat_cols) = num2cell(stat_vals(k,:));
    end

    r0 = 2 + nRuns + 5;
    sheet{r0+1,1} = 'Notes: values in bold are within 10% of the average value';

    writecell(sheet, outFileMS, 'Sheet', namesMS{f}, 'Range', 'B1');
    fprintf('%s: wrote %d run rows + stats to sheet "%s" of %s.\n', namesMS{f}, nRuns, namesMS{f}, outFileMS);

    boldRangesMS{f} = withinTol;
end

excel = actxserver('Excel.Application');
cleanupObj = onCleanup(@() cleanupExcel(excel));

excel.DisplayAlerts = false;
excel.Visible = false;
wb = excel.Workbooks.Open(fullfile(pwd, outFileMS));

colLetter = @(c) [char(64+floor((c-1)/26)*(c>26)) char(65+mod(c-1,26))];

for f = 1:numel(filesMS)
    fprintf('Formatting sheet: %s\n', namesMS{f});
    ws = wb.Sheets.Item(namesMS{f});
    withinTol = boldRangesMS{f};
    [nRuns, nC] = size(withinTol);
    for r = 1:nRuns
        for c = 1:nC
            if withinTol(r,c)
                excelRow = 2 + r;
                excelCol = 2 + c;
                addr = sprintf('%s%d', colLetter(excelCol), excelRow);
                try
                    ws.Range(addr).Font.Bold = true;
                catch ME
                    fprintf('FAILED at %s on sheet %s: %s\n', addr, namesMS{f}, ME.message);
                end
            end
        end
    end
end
wb.Save;
wb.Close(false);
fprintf('Applied bold formatting to %s.\n', outFileMS);

function cleanupExcel(excel)
    try
        excel.Quit;
    catch
    end
    delete(excel);
end