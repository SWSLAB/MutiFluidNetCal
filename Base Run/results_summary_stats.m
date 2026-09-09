%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

clc
clear
close all

%% Load all fluid results
fprintf('\n');
fprintf('================================================================================\n');
fprintf('                    COMPREHENSIVE RESULTS SUMMARY                               \n');
fprintf('================================================================================\n');

if ~isfile('fluid1.mat') || ~isfile('fluid2.mat') || ~isfile('fluid3.mat')
    error('Missing fluid files. Run main code 3 times with fluid_choice=1,2,3');
end

load('fluid1.mat'); fluid1 = fluid_results;
load('fluid2.mat'); fluid2 = fluid_results;
load('fluid3.mat'); fluid3 = fluid_results;

fluids = {fluid1, fluid2, fluid3};
fluid_names = {'Diesel', 'Jet fuel', 'Gasoline'};

g = 9.81; % m/s²

% Recover topology sizes from saved structs
Pipes     = fluid1.Pipes;
Pumps     = fluid1.Pumps;
n_pipes   = numel(Pipes);
num_pumps = numel(Pumps);

%% ========== DATA RETENTION SUMMARY ==========
fprintf('\n========================================================================\n');
fprintf('DATA RETENTION SUMMARY\n');
fprintf('========================================================================\n');

stage_names = {'Raw SCADA data'; 'Fluid homogeneity'; 'Non-operational conditions'; ...
    'Flow stability (signal peaks + mass balance)'; 'Missing data'; ...
    'Hydraulic feasibility (Optimization dataset)'};

combined_N = fluid1.N_combined;
print_retention_table(stage_names, combined_N, 1);

fprintf('\n--- FLUID-SPECIFIC OPTIMIZATION DATASETS ---\n');
for f = 1:3
    fprintf('%-12s: N = %d (%.1f%% of original)\n', ...
        fluid_names{f}, fluids{f}.N_selected_after, ...
        fluids{f}.N_selected_after / fluid1.N_combined(1) * 100);
end
fprintf('========================================================================\n\n');

%% ========== PER-FLUID DATA RETENTION BY STAGE ==========
fprintf('\n========================================================================\n');
fprintf('PER-FLUID DATA RETENTION BY STAGE\n');
fprintf('========================================================================\n');

N_fluid_all = [nan(3,1), fluid1.N_fluid_25(:,1:4), ...
    [fluid1.N_fluid_25(1,5); fluid2.N_fluid_25(2,5); fluid3.N_fluid_25(3,5)]];

for f = 1:3
    fprintf('\n%s:\n', fluid_names{f});
    print_retention_table(stage_names(2:end), N_fluid_all(f,2:end), 2);
end
fprintf('\n(Note: %% is relative to that fluid''s count after fluid homogeneity.)\n');
fprintf('========================================================================\n\n');

%% Local helper for retention tables
function print_retention_table(names, N, start_idx)
fprintf('%-46s %12s %12s %12s\n', 'Stage', 'N_rows', 'Kept(%%)', 'Discarded(%%)');
for k = 1:length(N)
    kept_pct = N(k) / N(1) * 100;
    if k == 1
        fprintf('%-46s %12d %12.0f %12s\n', names{k}, N(k), 100, '-');
    else
        fprintf('%-46s %12d %12.1f %12.1f\n', names{k}, N(k), kept_pct, 100-kept_pct);
    end
end
end

%% ========== SECTION 1: PARETO OPTIMIZATION RESULTS ==========
fprintf('\n');
fprintf('========== PARETO OPTIMIZATION RESULTS ==========\n');
fprintf('%-12s | %-10s | %-15s | %-15s | %-10s\n', 'Fluid', 'Knee W', 'RMSE_H (m)', 'RMSE_ΔH (m)', 'Flow Exp n');
fprintf('-------------|------------|-----------------|-----------------|------------\n');
for f = 1:3
    W = (fluids{f}.pareto_point_id - 1) * 0.1;
    rmse_h = fluids{f}.R2sysV(fluids{f}.pareto_point_id);
    rmse_dh = fluids{f}.R2fitV(fluids{f}.pareto_point_id);
    n_opt = mean(fluids{f}.n_sys);
    fprintf('%-12s | %10.2f | %15.2f | %15.2f | %10.3f\n', ...
        fluid_names{f}, W, rmse_h, rmse_dh, n_opt);
end

fprintf('\n%-12s | %-15s | %-15s\n', 'Fluid', 'RMSE_H (m)', 'RMSE_ΔH (m)');
fprintf('-------------|-----------------|------------------\n');
for f = 1:3
    fprintf('%-12s | %15.2f | %15.2f\n', ...
        fluid_names{f}, fluids{f}.R2sys_classical, fluids{f}.R2fit_classical);
end

fprintf('\nImprovement (Optimal vs Classical):\n');
fprintf('  H column: positive = optimized RMSE_H is LOWER (better) than classical.\n');
fprintf('  DH column: positive = optimized RMSE_DH is HIGHER (worse) than classical -- raw\n');
fprintf('  percent change, not an "improvement" sign, matching Table 3 in the manuscript.\n');
fprintf('%-12s | %-18s | %-22s\n', 'Fluid', 'H improve (%)', 'DH change (%, +=worse)');
fprintf('-------------|--------------------|------------------------\n');
for f = 1:3
    rmse_h_opt = fluids{f}.R2sysV(fluids{f}.pareto_point_id);
    rmse_h_cls = fluids{f}.R2sys_classical;
    rmse_dh_opt = fluids{f}.R2fitV(fluids{f}.pareto_point_id);
    rmse_dh_cls = fluids{f}.R2fit_classical;

    h_improvement = (rmse_h_cls - rmse_h_opt) / rmse_h_cls * 100;   % positive = optimal better
    dh_change     = (rmse_dh_opt - rmse_dh_cls) / rmse_dh_cls * 100; % positive = optimal worse (raw change)

    fprintf('%-12s | %18.1f | %24.1f\n', ...
        fluid_names{f}, h_improvement, dh_change);
end

%% ========== SECTION 1A: MATCHED ΔH COMPARISON (FOR FIGURE 11) ==========
fprintf('\n');
fprintf('========== MATCHED ΔH COMPARISON (Figure 11) ==========\n');
fprintf('Demonstrating improvement potential at constant head-loss error\n\n');

fprintf('%-12s | %-15s | %-15s | %-15s | %-15s\n', ...
    'Fluid', 'H_classical (m)', 'H_matched (m)', 'ΔH_classical (m)', 'ΔH_matched (m)');
fprintf('-------------|-----------------|-----------------|-----------------|------------------\n');

h_classical  = zeros(3,1);
h_matched    = zeros(3,1);
dh_classical = zeros(3,1);
dh_matched   = zeros(3,1);

for f = 1:3
    x_pareto = fluids{f}.R2sysV;
    y_pareto = fluids{f}.R2fitV;
    x_classical = fluids{f}.R2sys_classical;
    y_classical = fluids{f}.R2fit_classical;

    % Interpolate H value at EXACT classical ΔH (matches Fig. 9/11 method)
    h_matched(f) = interp1(y_pareto, x_pareto, y_classical, 'linear', 'extrap');

    h_classical(f)  = x_classical;
    dh_classical(f) = y_classical;
    dh_matched(f)   = y_classical;   % exact match by construction

    fprintf('%-12s | %15.2f | %15.2f | %15.2f | %15.2f\n', ...
        fluid_names{f}, h_classical(f), h_matched(f), dh_classical(f), dh_matched(f));
end

fprintf('\nImprovement at Matched ΔH:\n');
fprintf('%-12s | %-15s | %-15s | %-15s | %-15s\n', ...
    'Fluid', 'ΔH diff (m)', 'ΔH diff (%%)', 'H improvement (m)', 'H improvement (%%)');
fprintf('-------------|-----------------|-----------------|-----------------|------------------\n');

for f = 1:3
    dh_diff_m   = abs(dh_matched(f) - dh_classical(f));   % = 0 by construction
    dh_diff_pct = dh_diff_m / dh_classical(f) * 100;
    h_improvement_m   = h_classical(f) - h_matched(f);
    h_improvement_pct = h_improvement_m / h_classical(f) * 100;

    fprintf('%-12s | %15.3f | %15.1f | %15.2f | %15.1f\n', ...
        fluid_names{f}, dh_diff_m, dh_diff_pct, h_improvement_m, h_improvement_pct);
end

%% ========== SECTION 2: PIPE PARAMETERS ==========
fprintf('\n');
fprintf('========== PIPE PARAMETERS ==========\n');
num_pipes = length(fluids{1}.Pipes);

for i = 1:num_pipes
    fprintf('\n--- Pipe %d ---\n', i);
    fprintf('%-12s | %-12s | %-12s | %-12s | %-12s\n', ...
        'Fluid', 'r_opt', 'r_fitted', 'n_opt', 'n_fitted');
    fprintf('-------------|--------------|--------------|--------------|-------------\n');

    for f = 1:3
        r_opt = fluids{f}.R_sys(i) / fluids{f}.L(i);
        r_fit = fluids{f}.R_fit(i) / fluids{f}.L(i);
        n_opt = fluids{f}.n_sys(i);
        n_fit = fluids{f}.n_fit(i);

        % Check if pipe has flow
        has_flow = any(fluids{f}.I_Pipes(:,i) == 1) && ...
            any(fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,i)==1, i) > 1e-6);

        if has_flow
            fprintf('%-12s | %12.4f | %12.4f | %12.3f | %12.3f\n', ...
                fluid_names{f}, r_opt, r_fit, n_opt, n_fit);
        else
            fprintf('%-12s | %12s | %12s | %12s | %12s\n', ...
                fluid_names{f}, 'No flow', 'No flow', 'No flow', 'No flow');
        end
    end
end

%% ========== SECTION 3: PUMP PARAMETERS ==========
fprintf('\n');
fprintf('========== PUMP PARAMETERS ==========\n');
num_pumps = length(fluids{1}.Pumps);

for i = 1:num_pumps
    fprintf('\n--- Pump %d ---\n', i);
    fprintf('%-12s | %-12s | %-12s | %-12s | %-12s\n', ...
        'Fluid', 'A_opt (m)', 'A_fitted', 'B_opt', 'B_fitted');
    fprintf('-------------|--------------|--------------|--------------|-------------\n');

    for f = 1:3
        A_opt = fluids{f}.A_sys(i);
        A_fit = fluids{f}.A_fit(i);
        B_opt = fluids{f}.B_sys(i);
        B_fit = fluids{f}.B_fit(i);

        % Check if pump has flow
        has_flow = any(fluids{f}.I_Pumps(:,i) == 1) && ...
            any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1, i) > 1e-6);

        if has_flow
            fprintf('%-12s | %12.1f | %12.1f | %12.0f | %12.0f\n', ...
                fluid_names{f}, A_opt, A_fit, B_opt, B_fit);
        else
            fprintf('%-12s | %12s | %12s | %12s | %12s\n', ...
                fluid_names{f}, 'No flow', 'No flow', 'No flow', 'No flow');
        end
    end
end

%% ========== SECTION 4: VALVE PARAMETERS ==========
fprintf('\n');
fprintf('========== VALVE PARAMETERS ==========\n');
fprintf('%-12s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    'Fluid', 'α_opt', 'α_fitted', 'β_opt', 'β_fitted', 'd_eff (in)', 'Cv');
fprintf('-------------|--------------|--------------|--------------|-------------|--------------|-------------\n');

for f = 1:3
    d_eff_in = fluids{f}.d_valve5_inch;
    K_unitless_f = fluids{f}.alpha_sys * pi^2 * g * fluids{f}.d_valve5_m^4 / 8;
    Cv_f = sqrt(891 * d_eff_in^4 / K_unitless_f);

    fprintf('%-12s | %12.1f | %12.1f | %12.3f | %12.3f | %12.2f | %12.1f\n', ...
        fluid_names{f}, fluids{f}.alpha_sys, fluids{f}.alpha_fit, ...
        fluids{f}.beta_sys, fluids{f}.beta_fit, d_eff_in, Cv_f);
end

%% ========== SECTION 5: BOUNDS SUMMARY ==========
fprintf('\n');
fprintf('========== BOUNDS SUMMARY ==========\n');

% Pipe bounds
fprintf('\nPipe Resistance Bounds (r):\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Pipe', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('---------|--------------|--------------|-------------\n');
for i = 1:num_pipes
    fprintf('Pipe %d   | [%.3f,%.3f] | [%.3f,%.3f] | [%.3f,%.3f]\n', i, ...
        fluid1.r_min(i), fluid1.r_max(i), ...
        fluid2.r_min(i), fluid2.r_max(i), ...
        fluid3.r_min(i), fluid3.r_max(i));
end

% Flow exponent bounds
fprintf('\nFlow Exponent Bounds (n):\n');
fprintf('All pipes: [%.2f, %.2f]\n', fluid1.n_min(1), fluid1.n_max(1));

% Pump bounds
fprintf('\nPump A Bounds (m):\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Pump', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('---------|--------------|--------------|-------------\n');
for i = 1:num_pumps
    fprintf('Pump %d   | [%.1f,%.1f] | [%.1f,%.1f] | [%.1f,%.1f]\n', i, ...
        fluid1.A_min(i), fluid1.A_max(i), ...
        fluid2.A_min(i), fluid2.A_max(i), ...
        fluid3.A_min(i), fluid3.A_max(i));
end

fprintf('\nPump B Bounds (s²/m⁵):\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Pump', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('---------|--------------|--------------|-------------\n');
for i = 1:num_pumps
    fprintf('Pump %d   | [%.0f,%.0f] | [%.0f,%.0f] | [%.0f,%.0f]\n', i, ...
        fluid1.B_min(i), fluid1.B_max(i), ...
        fluid2.B_min(i), fluid2.B_max(i), ...
        fluid3.B_min(i), fluid3.B_max(i));
end

% Valve bounds
fprintf('\nValve Bounds:\n');
fprintf('%-12s | %-20s | %-20s\n', 'Fluid', 'α', 'β');
fprintf('-------------|----------------------|---------------------\n');
for f = 1:3
    fprintf('%-12s | [%.1f, %.1f] | [%.3f, %.3f]\n', ...
        fluid_names{f}, ...
        fluids{f}.alpha_min, fluids{f}.alpha_max, ...
        fluids{f}.beta_min, fluids{f}.beta_max);
end

%% ========== SECTION 6: CRITICAL BOUNDS ANALYSIS ==========
fprintf('\n');
fprintf('========== CRITICAL BOUNDS ANALYSIS ==========\n');
fprintf('(Checking if optimized values hit bounds)\n\n');

% Pipe n check
fprintf('Flow Exponent (n) - Binding Constraints:\n');
for f = 1:3
    n_opt = mean(fluids{f}.n_sys);
    n_min = fluids{f}.n_min(1);
    n_max = fluids{f}.n_max(1);

    if abs(n_opt - n_min) < 0.01
        status = '⚠ BINDING (at minimum)';
    elseif abs(n_opt - n_max) < 0.01
        status = '⚠ BINDING (at maximum)';
    else
        status = '✓ Interior';
    end

    fprintf('  %-12s: n=%.3f, bounds=[%.2f,%.2f] %s\n', ...
        fluid_names{f}, n_opt, n_min, n_max, status);
end

% Pump A/B check
fprintf('\nPump Parameters - Binding Constraints:\n');
for i = 1:num_pumps
    fprintf('  Pump %d:\n', i);
    for f = 1:3
        has_flow = any(fluids{f}.I_Pumps(:,i) == 1) && ...
            any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1, i) > 1e-6);

        if has_flow
            A_opt = fluids{f}.A_sys(i);
            A_min = fluids{f}.A_min(i);
            A_max = fluids{f}.A_max(i);

            if abs(A_opt - A_min)/A_min < 0.05
                A_status = '⚠ Near min';
            elseif abs(A_opt - A_max)/A_max < 0.05
                A_status = '⚠ Near max';
            else
                A_status = '✓ Interior';
            end

            fprintf('    %-12s: A=%.1f, bounds=[%.1f,%.1f] %s\n', ...
                fluid_names{f}, A_opt, A_min, A_max, A_status);
        end
    end
end

%% ========== SECTION 7: FLOW STATISTICS ==========
fprintf('\n');
fprintf('========== FLOW STATISTICS ==========\n');

fprintf('\nPipe Flow Coverage:\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Pipe', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('---------|--------------|--------------|-------------\n');
for i = 1:num_pipes
    for f = 1:3
        n_obs(f) = sum(fluids{f}.I_Pipes(:,i) == 1);
        if n_obs(f) > 0
            Q_avg(f) = mean(fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,i)==1, i)) * 3600;
        else
            Q_avg(f) = 0;
        end
    end
    fprintf('Pipe %d   | %4d (%.1f) | %4d (%.1f) | %4d (%.1f)\n', i, ...
        n_obs(1), Q_avg(1), n_obs(2), Q_avg(2), n_obs(3), Q_avg(3));
end
fprintf('(Numbers show: observation count (avg flow in m³/h))\n');

fprintf('\nPump Flow Coverage:\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Pump', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('---------|--------------|--------------|-------------\n');
for i = 1:num_pumps
    for f = 1:3
        n_obs(f) = sum(fluids{f}.I_Pumps(:,i) == 1);
        if n_obs(f) > 0
            Q_avg(f) = mean(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1, i)) * 3600;
        else
            Q_avg(f) = 0;
        end
    end
    fprintf('Pump %d   | %4d (%.1f) | %4d (%.1f) | %4d (%.1f)\n', i, ...
        n_obs(1), Q_avg(1), n_obs(2), Q_avg(2), n_obs(3), Q_avg(3));
end

%% ========== SECTION 8: PARAMETER RATIOS ==========
fprintf('\n');
fprintf('========== PARAMETER RATIOS (Optimal/Fitted) ==========\n');

fprintf('\nPipe Resistance Ratios:\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Pipe', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('---------|--------------|--------------|-------------\n');
for i = 1:num_pipes
    for f = 1:3
        has_flow = any(fluids{f}.I_Pipes(:,i) == 1) && ...
            any(fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,i)==1, i) > 1e-6);
        if has_flow
            ratio(f) = fluids{f}.R_sys(i) / fluids{f}.R_fit(i);
        else
            ratio(f) = NaN;
        end
    end

    if all(isnan(ratio))
        fprintf('Pipe %d   | %12s | %12s | %12s\n', i, 'No flow', 'No flow', 'No flow');
    else
        fprintf('Pipe %d   | %12.3f | %12.3f | %12.3f\n', i, ratio(1), ratio(2), ratio(3));
    end
end

fprintf('\nPump A Coefficient Ratios:\n');
fprintf('%-8s | %-12s | %-12s | %-12s\n', 'Pump', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('---------|--------------|--------------|-------------\n');
for i = 1:num_pumps
    for f = 1:3
        has_flow = any(fluids{f}.I_Pumps(:,i) == 1) && ...
            any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1, i) > 1e-6);
        if has_flow && fluids{f}.A_fit(i) > 0
            ratio(f) = fluids{f}.A_sys(i) / fluids{f}.A_fit(i);
        else
            ratio(f) = NaN;
        end
    end

    if all(isnan(ratio))
        fprintf('Pump %d   | %12s | %12s | %12s\n', i, 'No flow', 'No flow', 'No flow');
    else
        fprintf('Pump %d   | %12.3f | %12.3f | %12.3f\n', i, ratio(1), ratio(2), ratio(3));
    end
end

fprintf('\nValve Parameter Ratios:\n');
fprintf('%-12s | %-12s | %-12s\n', 'Fluid', 'α ratio', 'β ratio');
fprintf('-------------|--------------|-------------\n');
for f = 1:3
    alpha_ratio = fluids{f}.alpha_sys / fluids{f}.alpha_fit;
    beta_ratio = fluids{f}.beta_sys / fluids{f}.beta_fit;
    fprintf('%-12s | %12.3f | %12.3f\n', ...
        fluid_names{f}, alpha_ratio, beta_ratio);
end

%% ========== SECTION 9: VALVE COEFFICIENT ANALYSIS ==========
fprintf('\n');
fprintf('========== VALVE COEFFICIENT ANALYSIS ==========\n');

d_m = fluid1.d_valve5_m;
valve_ID_inches = fluid1.d_valve5_inch;
K_unitless_literature = fluid1.K_unitless_valve5;
K_Q_theoretical = fluid1.KQ_valve5_full_open;

fprintf('\nValve Geometry (from Main.m, K = %.1f):\n', K_unitless_literature);
fprintf('  Valve ID: %.2f inches (%.4f m)\n', valve_ID_inches, d_m);
fprintf('  K_Q theoretical: %.0f s²/m⁵\n', K_Q_theoretical);

% Step 3: Back-calculate K_unitless from optimized alpha values
fprintf('\nOptimized Alpha to K Conversion:\n');
fprintf('%-12s | %-12s | %-12s | %-15s | %-15s\n', ...
    'Fluid', 'α_opt (s²/m⁵)', 'K_unitless', 'K_literature', '% Difference');
fprintf('-------------|--------------|--------------|-----------------|----------------\n');

for f = 1:3
    K_unitless_f = fluids{f}.alpha_sys * pi^2 * g * fluids{f}.d_valve5_m^4 / 8;
    pct_diff_f = (K_unitless_f - K_unitless_literature) / K_unitless_literature * 100;

    fprintf('%-12s | %12.1f | %12.1f | %12.1f | %12.2f\n', ...
        fluid_names{f}, fluids{f}.alpha_sys, K_unitless_f, K_unitless_literature, pct_diff_f);
end

%% ========== SINGLE COMPREHENSIVE OPTIMAL PARAMETERS TABLE ==========
fprintf('\n');
fprintf('================================================================================\n');
fprintf('              OPTIMAL PARAMETERS SUMMARY - ALL DECISION VARIABLES              \n');
fprintf('================================================================================\n');
fprintf('\n');

% Create table header
fprintf('%-25s | %-15s | %-15s | %-15s\n', ...
    'Decision Variable', 'Diesel', 'Jet fuel', 'Gasoline');
fprintf('--------------------------|-----------------|-----------------|------------------\n');

% Global flow exponent
fprintf('%-25s | %15.3f | %15.3f | %15.3f\n', ...
    'n_avg (flow exponent)', ...
    mean(fluid1.n_sys), mean(fluid2.n_sys), mean(fluid3.n_sys));

fprintf('--------------------------|-----------------|-----------------|------------------\n');

% Pipe resistance coefficients
for i = 1:num_pipes
    row_label = sprintf('r_%d ', i);

    % Check flow for each fluid
    for f = 1:3
        has_flow = any(fluids{f}.I_Pipes(:,i) == 1) && ...
            any(fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,i)==1, i) > 1e-6);
        if has_flow
            r_val(f) = fluids{f}.R_sys(i) / fluids{f}.L(i);
            has_flow_flag(f) = true;
        else
            r_val(f) = NaN;
            has_flow_flag(f) = false;
        end
    end

    % Print row
    fprintf('%-25s |', row_label);
    for f = 1:3
        if has_flow_flag(f)
            fprintf(' %15.4f |', r_val(f));
        else
            fprintf(' %15s |', 'No flow');
        end
    end
    fprintf('\n');
end

fprintf('--------------------------|-----------------|-----------------|------------------\n');

% Pump A coefficients
for i = 1:num_pumps
    row_label = sprintf('A_%d (m)', i);

    % Check flow for each fluid
    for f = 1:3
        has_flow = any(fluids{f}.I_Pumps(:,i) == 1) && ...
            any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1, i) > 1e-6);
        if has_flow
            A_val(f) = fluids{f}.A_sys(i);
            has_flow_flag(f) = true;
        else
            A_val(f) = NaN;
            has_flow_flag(f) = false;
        end
    end

    % Print row
    fprintf('%-25s |', row_label);
    for f = 1:3
        if has_flow_flag(f)
            fprintf(' %15.1f |', A_val(f));
        else
            fprintf(' %15s |', 'No flow');
        end
    end
    fprintf('\n');
end

fprintf('--------------------------|-----------------|-----------------|------------------\n');

% Pump B coefficients
for i = 1:num_pumps
    row_label = sprintf('B_%d (s^2/m^5)', i);

    % Check flow for each fluid
    for f = 1:3
        has_flow = any(fluids{f}.I_Pumps(:,i) == 1) && ...
            any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1, i) > 1e-6);
        if has_flow
            B_val(f) = fluids{f}.B_sys(i);
            has_flow_flag(f) = true;
        else
            B_val(f) = NaN;
            has_flow_flag(f) = false;
        end
    end

    % Print row
    fprintf('%-25s |', row_label);
    for f = 1:3
        if has_flow_flag(f)
            fprintf(' %15.0f |', B_val(f));
        else
            fprintf(' %15s |', 'No flow');
        end
    end
    fprintf('\n');
end

fprintf('--------------------------|-----------------|-----------------|------------------\n');

% Valve alpha coefficient
fprintf('%-25s | %15.1f | %15.1f | %15.1f\n', ...
    'alpha (s^2/m^5)', fluid1.alpha_sys, fluid2.alpha_sys, fluid3.alpha_sys);

% Valve beta exponent
fprintf('%-25s | %15.3f | %15.3f | %15.3f\n', ...
    'beta', fluid1.beta_sys, fluid2.beta_sys, fluid3.beta_sys);

fprintf('--------------------------|-----------------|-----------------|------------------\n');

% Performance metrics
fprintf('\n--- PERFORMANCE METRICS ---\n');
fprintf('%-25s | %15.2f | %15.2f | %15.2f\n', ...
    'RMSE_H (m)', ...
    fluid1.R2sysV(fluid1.pareto_point_id), ...
    fluid2.R2sysV(fluid2.pareto_point_id), ...
    fluid3.R2sysV(fluid3.pareto_point_id));

fprintf('%-25s | %15.2f | %15.2f | %15.2f\n', ...
    'RMSE_dH (m)', ...
    fluid1.R2fitV(fluid1.pareto_point_id), ...
    fluid2.R2fitV(fluid2.pareto_point_id), ...
    fluid3.R2fitV(fluid3.pareto_point_id));

%% ========== END OF SUMMARY ==========
fprintf('\n');
fprintf('================================================================================\n');
fprintf('                         END OF RESULTS SUMMARY                                 \n');
fprintf('================================================================================\n');
fprintf('\n');
