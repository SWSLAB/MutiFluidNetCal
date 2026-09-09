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

fprintf('\n=== LOADING FLUID RESULTS ===\n');

if ~isfile('fluid1.mat') || ~isfile('fluid2.mat') || ~isfile('fluid3.mat')
    error('Missing fluid files. Run main code 3 times with fluid_choice=1,2,3');
end

load('fluid1.mat'); fluid1 = fluid_results;
load('fluid2.mat'); fluid2 = fluid_results;
load('fluid3.mat'); fluid3 = fluid_results;

fluids      = {fluid1, fluid2, fluid3};
fluid_names = {'Diesel','Jet fuel','Gasoline'};
colors      = [0 0 1; 0 0.8 0; 1 0 0];
markers     = {'s','o','+'};

Pipes    = fluid1.Pipes;
Pumps    = fluid1.Pumps;
PosValves = fluid1.PosValves;

fprintf('Loaded: Diesel, Jet, Gasoline\n');
fprintf('\n=== CREATING FIGURES ===\n');

%% Global Q range for pipelines and pumps (for unification)
% Only use positive, finite flows
all_Q_pipes = [];
all_Q_pumps = [];

for f = 1:3
    qP = fluids{f}.Q_Pipes;
    qP = qP(isfinite(qP) & qP > 0);
    all_Q_pipes = [all_Q_pipes; qP(:)];
    
    qPu = fluids{f}.Q_Pumps;
    qPu = qPu(isfinite(qPu) & qPu > 0);
    all_Q_pumps = [all_Q_pumps; qPu(:)];
end

if ~isempty(all_Q_pipes)
    Q_pipe_max = max(all_Q_pipes)*3600*1.05; % m^3/h
else
    Q_pipe_max = 100;
end

if ~isempty(all_Q_pumps)
    Q_pump_max = max(all_Q_pumps)*3600*1.05; % m^3/h
else
    Q_pump_max = 100;
end


%% Figure 1: Pipeline Bounds (All Fluids) 

fprintf('Creating Pipeline Bounds...\n');

figure('Position',[100,150,1000,500]);
set(gcf,'Color','w');

% Global r-bounds (same logic as original)
r_min_global = zeros(length(Pipes), 1);
r_max_global = zeros(length(Pipes), 1);
for i = 1:length(Pipes)
    valid_r_min = [];
    valid_r_max = [];
    for f = 1:3
        if fluids{f}.r_max(i) < 5.0 && fluids{f}.r_min(i) > 1e-5
            valid_r_min = [valid_r_min; fluids{f}.r_min(i)];
            valid_r_max = [valid_r_max; fluids{f}.r_max(i)];
        end
    end
    if ~isempty(valid_r_min)
        r_min_global(i) = min(valid_r_min);
        r_max_global(i) = max(valid_r_max);
    else
        r_min_global(i) = min([fluids{1}.r_min(i), fluids{2}.r_min(i), fluids{3}.r_min(i)]);
        r_max_global(i) = max([fluids{1}.r_max(i), fluids{2}.r_max(i), fluids{3}.r_max(i)]);
    end
end

n_min_global = min([fluid1.n_min; fluid2.n_min; fluid3.n_min]);
n_max_global = max([fluid1.n_max; fluid2.n_max; fluid3.n_max]);

% Left subplot – A: pipeline resistance
subplot('Position',[0.08, 0.18, 0.60, 0.72]); hold on; grid on; box on;
x_vals = 1:length(Pipes);

for i = 1:length(Pipes)
    fluids_with_flow = [];
    for f = 1:3
        has_flow = any(fluids{f}.I_Pipes(:,i)==1) && ...
                   any(fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,i)==1,i) > 1e-6);
        if has_flow
            fluids_with_flow = [fluids_with_flow, f];
        end
    end
    if length(fluids_with_flow) == 1
        f_single  = fluids_with_flow(1);
        r_center  = (fluids{f_single}.r_min(i) + fluids{f_single}.r_max(i))/2;
        e_low     = r_center - fluids{f_single}.r_min(i);
        e_high    = fluids{f_single}.r_max(i) - r_center;
        errorbar(x_vals(i), r_center, e_low, e_high, 'k','LineWidth',3);
    else
        r_center = (r_min_global(i) + r_max_global(i))/2;
        errorbar(x_vals(i), r_center, r_center-r_min_global(i), ...
                 r_max_global(i)-r_center, 'k','LineWidth',3);
    end
    for f = fluids_with_flow
        plot(x_vals(i), fluids{f}.R_sys(i)./fluids{f}.L(i), markers{f}, ...
            'Color',colors(f,:), 'MarkerSize',10, 'LineWidth',2);
    end
end

xlim([0.5 length(Pipes)+0.5]);
xticks(x_vals);
xticklabels(arrayfun(@(ii) sprintf('Pipeline %d',ii),1:length(Pipes),'UniformOutput',false));
ylabel('r_{l} (–)','FontSize',16);
title('A: Pipeline resistance','FontWeight','bold','FontSize',16);
set(gca,'FontSize',16);
yl = ylim; ylim([max(0,yl(1)-0.15*diff(yl)), yl(2)+0.15*diff(yl)]);

% Right subplot – B: flow exponent
subplot('Position',[0.75, 0.18, 0.18, 0.72]); hold on; grid on; box on;
n_center = (n_min_global+n_max_global)/2;
errorbar(1, n_center, n_center-n_min_global, n_max_global-n_center, ...
    'k','LineWidth',3);
for f = 1:3
    nVal = mean(fluids{f}.n_sys);
    plot(1, nVal, markers{f}, 'Color',colors(f,:), 'MarkerSize',10, 'LineWidth',2);
end
xlim([0.5 1.5]); set(gca,'XTick',[]);
ylabel('n (–)','FontSize',16);
title('B: Flow exponent','FontWeight','bold','FontSize',16);
set(gca,'FontSize',16);
yl = ylim; ylim([yl(1)-0.15*diff(yl), yl(2)+0.15*diff(yl)]);

% Legend below
h_bounds   = errorbar(NaN,NaN,NaN,'k','LineWidth',3);
h_diesel   = plot(NaN,NaN,markers{1},'Color',colors(1,:),'MarkerSize',10,'LineWidth',2);
h_jet      = plot(NaN,NaN,markers{2},'Color',colors(2,:),'MarkerSize',10,'LineWidth',2);
h_gasoline = plot(NaN,NaN,markers{3},'Color',colors(3,:),'MarkerSize',10,'LineWidth',2);
lgd = legend([h_bounds,h_diesel,h_jet,h_gasoline], ...
    {'Bounds','Diesel','Jet fuel','Gasoline'}, ...
    'Orientation','horizontal','FontSize',12);
lgd.Units = 'normalized';
lgd.Position = [0.35 0.01 0.30 0.05];

% Second subplot x-label "All Pipes"
subplot('Position',[0.75, 0.18, 0.18, 0.72]);
xlabel('All pipelines','FontSize',16);

saveas(gcf,'Fig1_Bounds_Pipelines.png');
fprintf('Saved pipeline bounds\n');


%% Figure 2: Pump Bounds (All Fluids) 
fprintf('Creating Pump Bounds...\n');

num_pumps = length(Pumps);
pump_labels = arrayfun(@(i) sprintf('P%d',i),1:num_pumps,'UniformOutput',false);

% Compute global bounds excluding wide fallbacks
A_min_global = zeros(num_pumps, 1);
A_max_global = zeros(num_pumps, 1);
B_min_global = zeros(num_pumps, 1);
B_max_global = zeros(num_pumps, 1);

for i = 1:num_pumps
    valid_A_min = [];
    valid_A_max = [];
    valid_B_min = [];
    valid_B_max = [];
    
    for f = 1:3
        % Check if pump has actual flow
        has_flow = any(fluids{f}.I_Pumps(:,i)==1) && any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1,i) > 1e-6);
        
        if has_flow
            % A bounds: use as computed by optimization_bounds.m for this fluid
            valid_A_min = [valid_A_min; fluids{f}.A_min(i)];
            valid_A_max = [valid_A_max; fluids{f}.A_max(i)];
    
            % B bounds: use as computed by optimization_bounds.m for this fluid
            valid_B_min = [valid_B_min; fluids{f}.B_min(i)];
            valid_B_max = [valid_B_max; fluids{f}.B_max(i)];
        end
    end
    
    if ~isempty(valid_A_min)
        A_min_global(i) = min(valid_A_min);
        A_max_global(i) = max(valid_A_max);
    else
        A_min_global(i) = min([fluids{1}.A_min(i), fluids{2}.A_min(i), fluids{3}.A_min(i)]);
        A_max_global(i) = max([fluids{1}.A_max(i), fluids{2}.A_max(i), fluids{3}.A_max(i)]);
    end
    if ~isempty(valid_B_min)
        B_min_global(i) = min(valid_B_min);
        B_max_global(i) = max(valid_B_max);
    else
        B_min_global(i) = min([fluids{1}.B_min(i), fluids{2}.B_min(i), fluids{3}.B_min(i)]);
        B_max_global(i) = max([fluids{1}.B_max(i), fluids{2}.B_max(i), fluids{3}.B_max(i)]);
    end
end

% Create figure with rectangular subplots (narrower width)
figure('Position',[300,150,1200,600]);
set(gcf,'Color','w');

for i = 1:num_pumps
    % 2 rows: A and B; columns: pumps
    % A parameter subplot
    subplot(2, num_pumps, i); hold on; grid on; box on;

    % Check if this pump has valid A data
    has_valid_A = ~isempty(valid_A_min) || any([fluids{1}.A_sys(i), fluids{2}.A_sys(i), fluids{3}.A_sys(i)] > 0.01);
    
    if has_valid_A && (A_max_global(i) - A_min_global(i)) > 1
        a_center = (A_min_global(i)+A_max_global(i))/2;
        errorbar(1, a_center, a_center-A_min_global(i), A_max_global(i)-a_center, ...
            'k','LineWidth',2,'CapSize',5);
        for f = 1:3
            has_flow = any(fluids{f}.I_Pumps(:,i)==1) && any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1,i) > 1e-6);
            if has_flow && fluids{f}.A_sys(i) > 0.01
                plot(1, fluids{f}.A_sys(i), markers{f}, 'Color',colors(f,:), 'MarkerSize',8, 'LineWidth',2);
            end
        end
    else
        % No valid data
        text(1, 0.5, 'No data', 'HorizontalAlignment','center', 'FontSize',12, 'Color',[0.5 0.5 0.5]);
    end
    xlim([0.5 1.5]); set(gca,'XTick',[], 'FontSize',16);
    
    % Individual y-axis scale for each pump's A parameter
    all_A_vals = [A_min_global(i); A_max_global(i)];
    for f = 1:3
        has_flow = any(fluids{f}.I_Pumps(:,i)==1) && ...
                   any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1,i) > 1e-6);
        if has_flow
            all_A_vals = [all_A_vals; fluids{f}.A_sys(i)];
        end
    end
    A_range_local = max(all_A_vals) - min(all_A_vals);
    A_ylim_local = [min(all_A_vals) - 0.15*A_range_local, max(all_A_vals) + 0.15*A_range_local];
    ylim(A_ylim_local);
    
    ylabel('A (m)', 'FontSize',16);
    title(sprintf('Pump %d', i), 'FontWeight','bold','FontSize',16);

    % B parameter subplot
    subplot(2, num_pumps, num_pumps + i); hold on; grid on; box on;

    % Check if this pump has valid B data
    B_range_check = (B_max_global(i) - B_min_global(i)) > 1;
    has_valid_B = ~isempty(valid_B_min) && B_range_check;
    
    if has_valid_B
        b_center = (B_min_global(i)+B_max_global(i))/2;
        errorbar(1, b_center, b_center-B_min_global(i), B_max_global(i)-b_center, ...
            'k','LineWidth',2,'CapSize',5);
        for f = 1:3
            has_flow = any(fluids{f}.I_Pumps(:,i)==1) && any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1,i) > 1e-6);
            B_range_f = fluids{f}.B_max(i) - fluids{f}.B_min(i);
            if has_flow && B_range_f > 1
                plot(1, fluids{f}.B_sys(i), markers{f}, 'Color',colors(f,:), 'MarkerSize',8, 'LineWidth',2);
            end
        end
    else
        % No valid data
        text(1, 0.5, 'No data', 'HorizontalAlignment','center', 'FontSize',12, 'Color',[0.5 0.5 0.5]);
    end
    xlim([0.5 1.5]); set(gca,'XTick',[], 'FontSize',16);
    
    % Individual y-axis scale for each pump's B parameter
    all_B_vals = [B_min_global(i); B_max_global(i)];
    for f = 1:3
        has_flow = any(fluids{f}.I_Pumps(:,i)==1) && ...
                   any(fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,i)==1,i) > 1e-6);
        if has_flow
            all_B_vals = [all_B_vals; fluids{f}.B_sys(i)];
        end
    end
    B_range_local = max(all_B_vals) - min(all_B_vals);
    B_ylim_local = [min(all_B_vals) - 0.15*B_range_local, max(all_B_vals) + 0.15*B_range_local];
    ylim(B_ylim_local);
    
    ylabel('B (s²/m⁵)', 'FontSize',16);
    title(sprintf('Pump %d', i), 'FontWeight','bold','FontSize',16);
end

% Common legend
h_bounds   = errorbar(NaN,NaN,NaN,'k','LineWidth',2);
h_diesel   = plot(NaN,NaN,markers{1},'Color',colors(1,:),'MarkerSize',8,'LineWidth',2);
h_jet      = plot(NaN,NaN,markers{2},'Color',colors(2,:),'MarkerSize',8,'LineWidth',2);
h_gasoline = plot(NaN,NaN,markers{3},'Color',colors(3,:),'MarkerSize',8,'LineWidth',2);
lgd = legend([h_bounds,h_diesel,h_jet,h_gasoline], ...
    {'Bounds','Diesel','Jet fuel','Gasoline'}, ...
    'Orientation','horizontal','Location','southoutside','FontSize',12);
lgd.Units = 'normalized';
lgd.Position = [0.30 0.01 0.40 0.05];

annotation('textbox',[0.5 0.96 0.0 0.0], ...
    'FitBoxToText','on','HorizontalAlignment','center', ...
    'LineStyle','none','FontWeight','bold','FontSize',16);

saveas(gcf,'Fig2_Bounds_Pumps.png');
fprintf('Saved pump bounds\n');

%% Figure 3: Valve Bounds (All Fluids) 

fprintf('Creating Valve Bounds...\n');

figure('Position',[200,200,600,450]);
set(gcf,'Color','w');

alpha_min_global = min([fluid1.alpha_min, fluid2.alpha_min, fluid3.alpha_min]);
alpha_max_global = max([fluid1.alpha_max, fluid2.alpha_max, fluid3.alpha_max]);
beta_min_global  = min([fluid1.beta_min,  fluid2.beta_min,  fluid3.beta_min]);
beta_max_global  = max([fluid1.beta_max,  fluid2.beta_max,  fluid3.beta_max]);

subplot(1,2,1); hold on; grid on; box on;
alpha_center = (alpha_min_global+alpha_max_global)/2;
errorbar(1,alpha_center,alpha_center-alpha_min_global,alpha_max_global-alpha_center, ...
    'k','LineWidth',3,'CapSize',10);
for f = 1:3
    plot(1, fluids{f}.alpha_sys, markers{f}, 'Color',colors(f,:), ...
        'MarkerSize',10, 'LineWidth',2);
end
xlim([0.8 1.2]); set(gca,'XTick',[]);
ylabel('\alpha (–)','FontSize',16);
title('(A)','FontWeight','bold','FontSize',16);
set(gca,'FontSize',16);
yl = ylim; ylim([yl(1)-0.15*diff(yl), yl(2)+0.15*diff(yl)]);

subplot(1,2,2); hold on; grid on; box on;
beta_center = (beta_min_global+beta_max_global)/2;
errorbar(1,beta_center,beta_center-beta_min_global,beta_max_global-beta_center, ...
    'k','LineWidth',3,'CapSize',10);
for f = 1:3
    plot(1, fluids{f}.beta_sys, markers{f}, 'Color',colors(f,:), ...
        'MarkerSize',10, 'LineWidth',2);
end
xlim([0.8 1.2]); set(gca,'XTick',[]);
ylabel('\beta (–)','FontSize',16);
title('(B)','FontWeight','bold','FontSize',16);
set(gca,'FontSize',16);
yl = ylim; ylim([yl(1)-0.15*diff(yl), yl(2)+0.15*diff(yl)]);

h_bounds   = errorbar(NaN,NaN,NaN,'k','LineWidth',3);
h_diesel   = plot(NaN,NaN,markers{1},'Color',colors(1,:),'MarkerSize',10,'LineWidth',2);
h_jet      = plot(NaN,NaN,markers{2},'Color',colors(2,:),'MarkerSize',10,'LineWidth',2);
h_gasoline = plot(NaN,NaN,markers{3},'Color',colors(3,:),'MarkerSize',10,'LineWidth',2);
lgd = legend([h_bounds,h_diesel,h_jet,h_gasoline], ...
    {'Bounds','Diesel','Jet fuel','Gasoline'}, ...
    'Orientation','horizontal','FontSize',12);
lgd.Units = 'normalized';
lgd.Position = [0.32 0.03 0.36 0.06];

saveas(gcf,'Fig3_Bounds_Valve.png');
fprintf('Saved valve bounds\n');

%% Figure 4: Pareto Fronts (All Fluids) 
fprintf('Creating Pareto Fronts...\n');
figure('Position',[100,100,1400,450]);
set(gcf,'Color','w');

all_y = [fluids{1}.R2fitV(:); fluids{2}.R2fitV(:); fluids{3}.R2fitV(:)];
all_y = all_y(isfinite(all_y));
if ~isempty(all_y)
    y_min_global = min(all_y);
    y_max_global = max(all_y);
    y_range = y_max_global - y_min_global;
    if y_range > eps
        y_limits = [y_min_global - 0.1*y_range, y_max_global + 0.1*y_range];
    else
        y_limits = [y_min_global - 1, y_max_global + 1];
    end
else
    y_limits = [0, 10];
end

for f = 1:3
    subplot(1,3,f); hold on; grid on; box on;

    % Pareto as line without circle markers
    h_pareto = plot(fluids{f}.R2sysV, fluids{f}.R2fitV, '-', ...
        'Color', colors(f,:), 'LineWidth', 2, 'DisplayName','Pareto Front');

    h_knee = plot(fluids{f}.R2sysV(fluids{f}.pareto_point_id), ...
                  fluids{f}.R2fitV(fluids{f}.pareto_point_id), ...
                  's', 'Color', colors(f,:), 'MarkerSize', 10, 'LineWidth', 2, ...
                  'DisplayName', 'Selected knee point');

    h_fit = plot(fluids{f}.R2sys_classical, fluids{f}.R2fit_classical, ...
                 'kx', 'MarkerSize', 10, 'LineWidth', 2, ...
                 'DisplayName', 'Fitting Solution');

    pareto_min   = min(fluids{f}.R2sysV);
    pareto_max   = max(fluids{f}.R2sysV);
    pareto_range = pareto_max - pareto_min;
    x_lim = [pareto_min - 0.1*pareto_range, pareto_max + 0.1*pareto_range];

    % Include classical solution in range
    if fluids{f}.R2sys_classical < x_lim(1)
        x_lim(1) = fluids{f}.R2sys_classical - 0.5;
    end
    if fluids{f}.R2sys_classical > x_lim(2)
        x_lim(2) = fluids{f}.R2sys_classical + 0.5;
    end

    % Round limits for cleaner axis
    x_lim(1) = floor(x_lim(1));
    x_lim(2) = ceil(x_lim(2));
    xlim(x_lim);
    
    % Determine appropriate tick spacing
    x_range = x_lim(2) - x_lim(1);
    if x_range <= 6
        tick_spacing = 1;  % Every 1 m
    elseif x_range <= 12
        tick_spacing = 2;  % Every 2 m
    else
        tick_spacing = 5;  % Every 5 m
    end
    
    % Generate tick positions
    x_tick_min = ceil(x_lim(1) / tick_spacing) * tick_spacing;
    x_tick_max = floor(x_lim(2) / tick_spacing) * tick_spacing;
    xticks(x_tick_min:tick_spacing:x_tick_max);

    ylim(y_limits);

    xlabel('Head domain (m)', 'FontSize', 16);
    ylabel('Head loss domain (m)', 'FontSize', 16);

    % Add A/B/C to subtitles
    if f == 1
        title(['(A) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    elseif f == 2
        title(['(B) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    else
        title(['(C) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    end

    set(gca, 'FontSize', 13);

    % Legend right corner
    legend([h_pareto, h_knee, h_fit], ...
        {'Pareto Front','Selected knee point','Fitting Solution'}, ...
        'Location','northeast','FontSize',12);
end

% Global A/B description (if you want it)
annotation('textbox',[0.5 0.96 0.0 0.0], ...
    'FitBoxToText','on','HorizontalAlignment','center', ...
    'LineStyle','none','FontWeight','bold','FontSize',16);

saveas(gcf, 'Fig9_Pareto_AllFluids.png');
fprintf('Saved Pareto fronts\n');

%% Figure 5: Flow Exponent Comparison 
fprintf('Creating Flow Exponent Comparison - Publication Version...\n');
figure('Position',[100,100,1600,500]);
set(gcf,'Color','w');

% Collect all valid n values to determine unified y-axis
all_n_sys = [mean(fluids{1}.n_sys), mean(fluids{2}.n_sys), mean(fluids{3}.n_sys)];
all_n_fit = [];
for f = 1:3
    valid_vals = fluids{f}.n_fit(fluids{f}.n_fit > 0);
    all_n_fit = [all_n_fit; valid_vals(:)];
end

% Unified y-axis range
y_min = min([all_n_sys, all_n_fit']) - 0.1;
y_max = max([all_n_sys, all_n_fit']) + 0.1;
y_lim = [y_min, y_max];

% Left subplot (A): narrower
subplot('Position',[0.08 0.15 0.32 0.75]);  % [left bottom width height]
hold on; grid on; box on;
n_opt = [mean(fluids{1}.n_sys), mean(fluids{2}.n_sys), mean(fluids{3}.n_sys)];
x = 1:3;
bar_width = 0.6;
for i = 1:3
    bar(x(i), n_opt(i), bar_width, 'FaceColor', colors(i,:));
    text(x(i), n_opt(i), sprintf('%.2f', n_opt(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
end
set(gca,'XTick',x,'XTickLabel',fluid_names,'FontSize',16);
ylabel('n (–)', 'FontSize',16);
title('(A) System Optimization', 'FontSize',16,'FontWeight','bold');
ylim(y_lim);

% Right subplot (B): wider
subplot('Position',[0.46 0.15 0.50 0.75]);
hold on; grid on; box on;

num_pipes = length(Pipes);
x_base    = 1:num_pipes;
bar_width = 0.22;

for f = 1:3
    x_pos = x_base + (f-2)*bar_width;

    % Read directly from data: only keep pipes where this fluid actually flows
    n_fitted_values = fluids{f}.n_fit;
    has_flow = any(fluids{f}.I_Pipes == 1, 1);      % logical 1×num_pipes
    n_fitted_values(~has_flow) = NaN;               % hide pipes with no flow

    b(f) = bar(x_pos, n_fitted_values, bar_width, 'FaceColor', colors(f,:));

    for i = 1:num_pipes
        if ~isnan(n_fitted_values(i)) && n_fitted_values(i) > 0
            text(x_pos(i), n_fitted_values(i), sprintf('%.2f', n_fitted_values(i)), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                'FontSize',9,'FontWeight','bold');
        end
    end
end

set(gca,'XTick',x_base, ...
    'XTickLabel',arrayfun(@(i) sprintf('Pipeline %d',i),1:num_pipes,'UniformOutput',false), ...
    'FontSize',16);
ylabel('n (–)', 'FontSize',16);
title('(B) Pipeline-Specific Fitting', 'FontSize',16,'FontWeight','bold');
legend(b, fluid_names, 'Location','northeast','FontSize',12);
ylim(y_lim);

saveas(gcf, 'Fig4_FlowExponent_Comparison_Paper.png');
fprintf('Saved flow exponent comparison for paper\n');

%% Figure 6: Pipeline 1 Curve (All Fluids) 
fprintf('Creating Pipeline 1 Curve (All Fluids) - Publication Version...\n');
pipe_idx = 1;  % Select Pipeline 1 for paper
figure('Position',[100,100,1800,500]);
set(gcf,'Color','w');
% Calculate unified Q range across all fluids for Pipeline 1
Q_min_all = inf;
Q_max_all = -inf;
for f=1:3
    Q_data = fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,pipe_idx)==1,pipe_idx);
    if ~isempty(Q_data)
        Q_min_all = min(Q_min_all, min(Q_data));
        Q_max_all = max(Q_max_all, max(Q_data));
    end
end
% Calculate unified y-axis range across all fluids for Pipeline 1
y_min_all = inf;
y_max_all = -inf;
for f=1:3
    Q_data = fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,pipe_idx)==1,pipe_idx);
    dh_data = fluids{f}.dh_Pipes(fluids{f}.I_Pipes(:,pipe_idx)==1,pipe_idx);
    if ~isempty(Q_data)
        y_min_all = min(y_min_all, min(dh_data));
        y_max_all = max(y_max_all, max(dh_data));
        % Also check curve values
        xx = linspace(Q_min_all, Q_max_all, 200);
        yy_sys = fluids{f}.R_sys(pipe_idx) * xx.^fluids{f}.n_sys(pipe_idx);
        yy_fit = fluids{f}.R_fit(pipe_idx) * xx.^fluids{f}.n_fit(pipe_idx);
        y_min_all = min([y_min_all, min(yy_sys), min(yy_fit)]);
        y_max_all = max([y_max_all, max(yy_sys), max(yy_fit)]);
    end
end
y_range = y_max_all - y_min_all;
y_lim = [y_min_all - 0.05*y_range, y_max_all + 0.15*y_range];
for f=1:3
    subplot(1,3,f); hold on; grid on; box on;
    
    Q_data = fluids{f}.Q_Pipes(fluids{f}.I_Pipes(:,pipe_idx)==1,pipe_idx);
    dh_data = fluids{f}.dh_Pipes(fluids{f}.I_Pipes(:,pipe_idx)==1,pipe_idx);
    
    if ~isempty(Q_data)
        % Plot SCADA data
        h_scada = scatter(Q_data*3600, dh_data, 40, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.5);
        
        % Plot curves using unified Q range
        xx = linspace(Q_min_all, Q_max_all, 200);
        yy_sys = fluids{f}.R_sys(pipe_idx) * xx.^fluids{f}.n_sys(pipe_idx);
        yy_fit = fluids{f}.R_fit(pipe_idx) * xx.^fluids{f}.n_fit(pipe_idx);
        
        h_opt = plot(xx*3600, yy_sys, '-', 'Color', [0 0 1], 'LineWidth', 3);
        h_fit = plot(xx*3600, yy_fit, '-', 'Color', [1 0 0], 'LineWidth', 3);
        
        % Text box with ratios - positioned in right corner
        r_opt = fluids{f}.R_sys(pipe_idx) / fluids{f}.L(pipe_idx);
        r_fitted = fluids{f}.R_fit(pipe_idx) / fluids{f}.L(pipe_idx);
        r_ratio = r_opt / r_fitted;
        
        txt = sprintf('r_{opt}/r_{fitted}=%.2f', r_ratio);
        
        text(0.98, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'right', 'BackgroundColor', 'white', ...
            'EdgeColor', 'black', 'FontSize', 12);
        
        % Set unified y-axis
        ylim(y_lim);
    end
    
    xlabel('Flow Rate (m³/h)', 'FontSize', 16);
    ylabel('Head Loss (m)', 'FontSize', 16);
    % Add A, B, C to subtitles
    if f == 1
        title(['(A) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    elseif f == 2
        title(['(B) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    else
        title(['(C) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    end
end
% Add figure-level legend
if exist('h_scada','var') && exist('h_opt','var') && exist('h_fit','var')
    lgd = legend([h_scada, h_opt, h_fit], {'SCADA', 'System Optimization', 'Pipeline-Specific Fitting'}, ...
        'Orientation', 'horizontal', 'FontSize', 12);
    lgd.Position = [0.30, -0.005, 0.40, 0.05];
end
% No main title for paper
saveas(gcf, sprintf('Fig5_Pipeline1_AllFluids_Paper.png', pipe_idx));
fprintf('Saved Pipeline 1 figure for paper\n');

%% Figure 7: Pump 1 Curve (All Fluids) 
fprintf('Creating Pump 1 Curve (All Fluids) - Publication Version...\n');
pump_idx = 1;  % Select Pump 1 for paper
figure('Position',[100,100,1800,500]);
set(gcf,'Color','w');
% Calculate unified Q range across all fluids for Pump 1
Q_min_all = 0;  % Pumps start at 0
Q_max_all = -inf;
for f=1:3
    Q_data = fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,pump_idx)==1,pump_idx);
    if ~isempty(Q_data) && max(Q_data) > 0
        Q_max_all = max(Q_max_all, max(Q_data));
    end
end
Q_max_all = Q_max_all * 1.1;  % Add 10% margin
% Calculate unified y-axis range across all fluids for Pump 1
y_min_all = inf;
y_max_all = -inf;
for f=1:3
    Q_data = fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,pump_idx)==1,pump_idx);
    dh_data = fluids{f}.dh_Pumps(fluids{f}.I_Pumps(:,pump_idx)==1,pump_idx);
    if ~isempty(Q_data) && max(Q_data) > 0
        y_min_all = min(y_min_all, min(dh_data));
        y_max_all = max(y_max_all, max(dh_data));
        % Also check curve values
        xx = linspace(0, Q_max_all, 200);
        yy_sys = fluids{f}.A_sys(pump_idx) - fluids{f}.B_sys(pump_idx) * xx.^2;
        yy_fit = fluids{f}.A_fit(pump_idx) - fluids{f}.B_fit(pump_idx) * xx.^2;
        y_min_all = min([y_min_all, min(yy_sys), min(yy_fit)]);
        y_max_all = max([y_max_all, max(yy_sys), max(yy_fit)]);
    end
end
y_range = y_max_all - y_min_all;
y_lim = [y_min_all - 0.05*y_range, y_max_all + 0.15*y_range];
for f=1:3
    subplot(1,3,f); hold on; grid on; box on;
    
    Q_data = fluids{f}.Q_Pumps(fluids{f}.I_Pumps(:,pump_idx)==1,pump_idx);
    dh_data = fluids{f}.dh_Pumps(fluids{f}.I_Pumps(:,pump_idx)==1,pump_idx);
    
    if ~isempty(Q_data) && max(Q_data) > 0
        % Plot SCADA data
        h_scada = scatter(Q_data*3600, dh_data, 40, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.5);
        
        % Plot curves using unified Q range
        xx = linspace(0, Q_max_all, 200);
        yy_sys = fluids{f}.A_sys(pump_idx) - fluids{f}.B_sys(pump_idx) * xx.^2;
        yy_fit = fluids{f}.A_fit(pump_idx) - fluids{f}.B_fit(pump_idx) * xx.^2;
        
        h_opt = plot(xx*3600, yy_sys, '-', 'Color', [0 0 1], 'LineWidth', 3);
        h_fit = plot(xx*3600, yy_fit, '-', 'Color', [1 0 0], 'LineWidth', 3);
        
        % Text box with ratios - positioned in right corner
        A_ratio = fluids{f}.A_sys(pump_idx) / fluids{f}.A_fit(pump_idx);
        B_ratio = fluids{f}.B_sys(pump_idx) / fluids{f}.B_fit(pump_idx);
        
        txt = sprintf('A_{opt}/A_{fitted}=%.2f\nB_{opt}/B_{fitted}=%.2f', A_ratio, B_ratio);
        
        text(0.98, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'right', 'BackgroundColor', 'white', ...
            'EdgeColor', 'black', 'FontSize', 12);
        
        % Set unified x and y axes
        xlim([0, Q_max_all*3600]);
        ylim(y_lim);
    else
        % No data for this pump-fluid combination
        text(0.5, 0.5, 'No data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', [0.5 0.5 0.5]);
    end
    
    xlabel('Flow Rate (m³/h)', 'FontSize', 16);
    ylabel('Head Gain (m)', 'FontSize', 16);
    % Add A, B, C to subtitles
    if f == 1
        title(['(A) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    elseif f == 2
        title(['(B) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    else
        title(['(C) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    end
end
% Add figure-level legend - moved below x-axis
if exist('h_scada','var') && exist('h_opt','var') && exist('h_fit','var')
    lgd = legend([h_scada, h_opt, h_fit], {'SCADA', 'System Optimization', 'Pump-Specific Fitting'}, ...
        'Orientation', 'horizontal', 'FontSize', 12);
    lgd.Position = [0.30, -0.005, 0.40, 0.05];
end
% No main title for paper
saveas(gcf, sprintf('Fig6_Pump1_AllFluids_Paper.png', pump_idx));
fprintf('Saved Pump 1 figure for paper\n');

%% Figure 8: Valve Curve (All Fluids) 
fprintf('Creating Valve Curve (All Fluids) - Publication Version...\n');
figure('Position',[100,100,1800,500]);
set(gcf,'Color','w');
for f=1:3
    subplot(1,3,f); hold on; grid on; box on;
    
    x_data = fluids{f}.x_v5(fluids{f}.I_Valves(:,PosValves)==1);
    Q_data = fluids{f}.Q_Valves(fluids{f}.I_Valves(:,PosValves)==1,PosValves);
    dh_data = fluids{f}.dh_Valves(fluids{f}.I_Valves(:,PosValves)==1,PosValves);
    K_data = dh_data./Q_data.^2;
    
    if ~isempty(x_data)
        % Plot SCADA data
        h_scada = scatter(x_data*100, K_data, 40, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.5);
        
        % Plot curves
        xx = linspace(min(x_data), max(x_data), 200);
        yy_sys = fluids{f}.alpha_sys * xx.^fluids{f}.beta_sys;
        yy_fit = fluids{f}.alpha_fit * xx.^fluids{f}.beta_fit;
        
        h_opt = plot(xx*100, yy_sys, '-', 'Color', [0 0 1], 'LineWidth', 3);
        h_fit = plot(xx*100, yy_fit, '-', 'Color', [1 0 0], 'LineWidth', 3);
        
        % Text box with ratios - positioned in right corner
        alpha_ratio = fluids{f}.alpha_sys / fluids{f}.alpha_fit;
        beta_ratio = fluids{f}.beta_sys / fluids{f}.beta_fit;
        
        txt = sprintf('α_{opt}/α_{fitted}=%.2f\nβ_{opt}/β_{fitted}=%.2f', alpha_ratio, beta_ratio);
        
        text(0.98, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'right', 'BackgroundColor', 'white', ...
            'EdgeColor', 'black', 'FontSize', 12);
    end
    
    xlabel('V5 Degree Opening (%)', 'FontSize', 16);
    ylabel('K-factor (s²/m⁵)', 'FontSize', 16);
    % Add A, B, C to subtitles
    if f == 1
        title(['(A) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    elseif f == 2
        title(['(B) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    else
        title(['(C) ' fluid_names{f}], 'FontWeight', 'bold', 'FontSize', 16);
    end
    xlim([30, 100]);  % Extended below 50% to show white space
end
% Add figure-level legend - moved below x-axis
if exist('h_scada','var') && exist('h_opt','var') && exist('h_fit','var')
    lgd = legend([h_scada, h_opt, h_fit], {'SCADA', 'System Optimization', 'Valve-Specific Fitting'}, ...
        'Orientation', 'horizontal', 'FontSize', 12);
    lgd.Position = [0.30, -0.005, 0.40, 0.05];
end
% No main title for paper
saveas(gcf, 'Fig7_Valve_AllFluids_Paper.png');
fprintf('Saved valve figure for paper\n');

%% Figure 9: Head Domain Improvement Bar Chart at Matched Head Loss
fprintf('Creating Head Domain Improvement Bar Chart...\n');
figure('Position',[100,100,1600,600]);
set(gcf,'Color','w');

% Use same colors as other figures
colors = [0 0 1; 0 0.8 0; 1 0 0];  % Diesel=blue, Jet=green, Gasoline=red

% For each fluid, interpolate Pareto point at EXACT classical ΔH
h_classical = zeros(3,1);
h_matched = zeros(3,1);
dh_classical = zeros(3,1);

for f = 1:3
    x_pareto = fluids{f}.R2sysV;
    y_pareto = fluids{f}.R2fitV;
    x_classical_val = fluids{f}.R2sys_classical;
    y_classical_val = fluids{f}.R2fit_classical;
    
    % **INTERPOLATE H value at EXACT classical ΔH**
    % Use linear interpolation on the Pareto curve
    h_matched(f) = interp1(y_pareto, x_pareto, y_classical_val, 'linear', 'extrap');
    
    h_classical(f) = x_classical_val;
    dh_classical(f) = y_classical_val;
end

% ΔH values for matched points are EXACTLY the classical values
dh_matched = dh_classical;  % **EXACT match by construction**

% **UNIFIED Y-AXIS LIMIT FOR BOTH PANELS**
y_max_unified = max([h_classical; h_matched; dh_classical; dh_matched]) * 1.15;

x = 1:3;
bar_width = 0.35;

% Panel A: Head domain (RMSE_H)
subplot(1,2,1); hold on; grid on; box on;
b1 = bar(x - bar_width/2, h_matched, bar_width, 'FaceColor', [0 0 1]);
b2 = bar(x + bar_width/2, h_classical, bar_width, 'FaceColor', [1 0 0]);

for i = 1:3
    text(x(i) - bar_width/2, h_matched(i), sprintf('%.2f', h_matched(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
    text(x(i) + bar_width/2, h_classical(i), sprintf('%.2f', h_classical(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
end

set(gca,'XTick',x,'XTickLabel',fluid_names,'FontSize',16);
ylabel('RMSE (m)', 'FontSize',16);
title('(A) Head domain', 'FontSize',16,'FontWeight','bold');
ylim([0, y_max_unified]);
legend([b1, b2], {'System Optimization','Element-Specific Fitting'}, ...
    'Location','northwest','FontSize',12);

% Panel B: Head loss domain (RMSE_ΔH) - bars will be IDENTICAL
subplot(1,2,2); hold on; grid on; box on;
b1 = bar(x - bar_width/2, dh_matched, bar_width, 'FaceColor', [0 0 1]);
b2 = bar(x + bar_width/2, dh_classical, bar_width, 'FaceColor', [1 0 0]);

for i = 1:3
    text(x(i) - bar_width/2, dh_matched(i), sprintf('%.2f', dh_matched(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
    text(x(i) + bar_width/2, dh_classical(i), sprintf('%.2f', dh_classical(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
end

set(gca,'XTick',x,'XTickLabel',fluid_names,'FontSize',16);
title('(B) Head loss domain', 'FontSize',16,'FontWeight','bold');
ylim([0, y_max_unified]);
legend([b1, b2], {'System Optimization','Element-Specific Fitting'}, ...
    'Location','northwest','FontSize',12);

saveas(gcf, 'Fig_Improvement_BarChart.png');
fprintf('Saved improvement bar chart\n');

%% Figure 10: RMSE Comparison 
fprintf('Creating RMSE Comparison Figure...\n');
figure('Position',[100,100,1600,600]);
set(gcf,'Color','w');

% Optimal (selected Pareto point) values
h_domain_opt = [fluids{1}.R2sysV(fluids{1}.pareto_point_id), ...
                fluids{2}.R2sysV(fluids{2}.pareto_point_id), ...
                fluids{3}.R2sysV(fluids{3}.pareto_point_id)];

dh_domain_opt = [fluids{1}.R2fitV(fluids{1}.pareto_point_id), ...
                 fluids{2}.R2fitV(fluids{2}.pareto_point_id), ...
                 fluids{3}.R2fitV(fluids{3}.pareto_point_id)];

% Classical fitting values (stored separately, not from Pareto sweep)
h_domain_fitted = [fluids{1}.R2sys_classical, ...
                   fluids{2}.R2sys_classical, ...
                   fluids{3}.R2sys_classical];

dh_domain_fitted = [fluids{1}.R2fit_classical, ...
                    fluids{2}.R2fit_classical, ...
                    fluids{3}.R2fit_classical];

% Unified y-axis limit
y_max = max([h_domain_opt, h_domain_fitted, dh_domain_opt, dh_domain_fitted]) * 1.15;

x = 1:3;
bar_width = 0.35;

% Panel A: Head domain
subplot(1,2,1); hold on; grid on; box on;
b1 = bar(x - bar_width/2, h_domain_opt,    bar_width, 'FaceColor', [0 0 1]);
b2 = bar(x + bar_width/2, h_domain_fitted, bar_width, 'FaceColor', [1 0 0]);
for i = 1:3
    text(x(i) - bar_width/2, h_domain_opt(i), sprintf('%.2f', h_domain_opt(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
    text(x(i) + bar_width/2, h_domain_fitted(i), sprintf('%.2f', h_domain_fitted(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
end
set(gca,'XTick',x,'XTickLabel',fluid_names,'FontSize',16);
ylabel('RMSE (m)', 'FontSize',16);
title('(A) Head domain', 'FontSize',16,'FontWeight','bold');
ylim([0, y_max]);
legend([b1, b2], {'System Optimization','Element-Specific Fitting'}, ...
    'Location','northwest','FontSize',12);

% Panel B: Head loss domain
subplot(1,2,2); hold on; grid on; box on;
b1 = bar(x - bar_width/2, dh_domain_opt,    bar_width, 'FaceColor', [0 0 1]);
b2 = bar(x + bar_width/2, dh_domain_fitted, bar_width, 'FaceColor', [1 0 0]);
for i = 1:3
    text(x(i) - bar_width/2, dh_domain_opt(i), sprintf('%.2f', dh_domain_opt(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
    text(x(i) + bar_width/2, dh_domain_fitted(i), sprintf('%.2f', dh_domain_fitted(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',12,'FontWeight','bold');
end
set(gca,'XTick',x,'XTickLabel',fluid_names,'FontSize',16);
title('(B) Head loss domain', 'FontSize',16,'FontWeight','bold');
ylim([0, y_max]);
legend([b1, b2], {'System Optimization','Element-Specific Fitting'}, ...
    'Location','northwest','FontSize',12);

saveas(gcf, 'Fig8_RMSE_Comparison.png');
fprintf('Saved RMSE comparison figure\n');





