%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

function [r_min, r_max, A_min, A_max, B_min, B_max, alpha_min, alpha_max, beta_min, beta_max, n_min, n_max] = ...
    optimization_bounds(Hobs, Q_Pipes, Q_Pumps, Q_Valves, x_v5, I_Pipes, I_Pumps, I_Valves, L, O_Pipes, D_Pipes, O_Pumps, D_Pumps, O_Valves, D_Valves, Pipes, Pumps, PosValves, R_fit, A_fit, B_fit)

%% Fixed n bounds (physical turbulent flow range)
n_min = 1.00 * ones(length(Pipes), 1);
n_max = 2.00 * ones(length(Pipes), 1);

%% Pipe resistance bounds
r_min = zeros(length(Pipes), 1);
r_max = zeros(length(Pipes), 1);

for i = 1:length(Pipes)
    dh_Pipes = Hobs(I_Pipes(:, i) == 1, O_Pipes(i)) - Hobs(I_Pipes(:, i) == 1, D_Pipes(i));
    Q_pipe = Q_Pipes(I_Pipes(:, i) == 1, i);

    if isempty(Q_pipe)
        continue;
    end

    % Calculate r at both extreme n values (no loop needed)
    r_at_n_min = dh_Pipes ./ (Q_pipe .^ n_min(i) * L(i));
    r_at_n_max = dh_Pipes ./ (Q_pipe .^ n_max(i) * L(i));
    r_all = [r_at_n_min; r_at_n_max];

    r_reasonable = r_all(r_all > 0 & r_all < 10 & isfinite(r_all));
    if ~isempty(r_reasonable)
        r_min(i) = min(r_reasonable);
        r_max(i) = max(r_reasonable);
    end
end

% Fallback for pipes with no data
for i = 1:length(Pipes)
    if r_min(i) == 0 && r_max(i) == 0
        % Check if pipe has any flow
        has_flow = any(I_Pipes(:,i) == 1) && any(Q_Pipes(I_Pipes(:,i)==1, i) > 1e-6);

        if has_flow
            % Pipe operates but bounds calculation failed
            r_min(i) = (R_fit(i) / L(i)) * 0.5;
            r_max(i) = (R_fit(i) / L(i)) * 2.0;
        else
            % Pipe never operates for this fluid - wide bounds
            r_min(i) = 1e-6;
            r_max(i) = 10.0;
            fprintf('Warning: Pipe %d has no flow data - using wide r bounds\n', i);
        end
    end
end

%% Pump bounds (A, B)
A_min = zeros(length(Pumps), 1);
A_max = zeros(length(Pumps), 1);
B_min = zeros(length(Pumps), 1);
B_max = zeros(length(Pumps), 1);

for i = 1:length(Pumps)
    if i ~= 1
        dh_Pumps = Hobs(I_Pumps(:, i) == 1, D_Pumps(i)) - Hobs(I_Pumps(:, i) == 1, O_Pumps(i));
    else
        dh_Pumps = Hobs(I_Pumps(:, i) == 1, D_Pumps(i) + 1) - Hobs(I_Pumps(:, i) == 1, O_Pumps(i));
    end
    Q_pump = Q_Pumps(I_Pumps(:, i) == 1, i);

    if length(Q_pump) < 2
        continue;
    end

    pairs = nchoosek(1:length(Q_pump), 2);
    Q1 = Q_pump(pairs(:, 1));
    Q2 = Q_pump(pairs(:, 2));
    H1 = dh_Pumps(pairs(:, 1));
    H2 = dh_Pumps(pairs(:, 2));

    flow_sep = abs(Q1 - Q2) ./ ((Q1 + Q2) / 2);
    sep_mask = flow_sep > 0.10;

    if ~any(sep_mask)
        continue;
    end

    Q1 = Q1(sep_mask);
    Q2 = Q2(sep_mask);
    H1 = H1(sep_mask);
    H2 = H2(sep_mask);

    denom = Q1.^2 - Q2.^2;
    valid = abs(denom) > 1e-4;

    if ~any(valid)
        continue;
    end

    B_pairs = -(H1(valid) - H2(valid)) ./ denom(valid);
    A_pairs = H1(valid) + B_pairs .* (Q1(valid).^2);

    finite_mask = isfinite(A_pairs) & isfinite(B_pairs) & (A_pairs > 0) & (B_pairs > 0);

    if any(finite_mask)
        A_min(i) = min(A_pairs(finite_mask));
        A_max(i) = max(A_pairs(finite_mask));
        B_min(i) = min(B_pairs(finite_mask));
        B_max(i) = max(B_pairs(finite_mask));
    end
end

% Fallback for pumps with no nchoosek data
for i = 1:length(Pumps)
    if A_min(i) == 0 && A_max(i) == 0
        % Check if pump has any flow
        has_flow = any(I_Pumps(:,i) == 1) && any(Q_Pumps(I_Pumps(:,i)==1, i) > 1e-6);

        if has_flow
            % Pump operates but nchoosek failed - use A_fit based bounds
            A_min(i) = max(1e-6, A_fit(i) * 0.5);
            A_max(i) = A_fit(i) * 2.0;
        else
            % Pump never operates for this fluid - VERY wide bounds
            A_min(i) = 1e-6;
            A_max(i) = 2000;
            fprintf('Warning: Pump %d has no flow data - using wide A bounds\n', i);
        end
    end

    if B_min(i) == 0 && B_max(i) == 0
        has_flow = any(I_Pumps(:,i) == 1) && any(Q_Pumps(I_Pumps(:,i)==1, i) > 1e-6);

        if has_flow
            B_min(i) = max(1e-6, B_fit(i) * 0.5);
            B_max(i) = B_fit(i) * 2.0;
        else
            B_min(i) = 1e-6;
            B_max(i) = 1e6;
            fprintf('Warning: Pump %d has no flow data - using wide B bounds\n', i);
        end
    end
end

%% Valve bounds (alpha, beta)
dh_Valve = Hobs(I_Valves(:, PosValves) == 1, O_Valves(PosValves)) - Hobs(I_Valves(:, PosValves) == 1, D_Valves(PosValves));
Q_valve = Q_Valves(I_Valves(:, PosValves) == 1, PosValves);
v_pos = x_v5(I_Valves(:, PosValves) == 1);

K_flow = dh_Valve ./ (Q_valve.^2);

valid_mask = (K_flow > 100) & (K_flow < 1e6) & (v_pos >= 0.20) & (v_pos <= 1.00) & isfinite(K_flow) & (Q_valve > 0);
K_clean = K_flow(valid_mask);
v_clean = v_pos(valid_mask);

if length(v_clean) >= 2
    pairs = nchoosek(1:length(v_clean), 2);
    v1 = v_clean(pairs(:, 1));
    v2 = v_clean(pairs(:, 2));
    K1 = K_clean(pairs(:, 1));
    K2 = K_clean(pairs(:, 2));

    v_sep = abs(v1 - v2) ./ max(v1, v2);
    sep_mask = (v_sep > 0.10) & (v1 > 0) & (v2 > 0) & (K1 > 0) & (K2 > 0);

    if any(sep_mask)
        v1 = v1(sep_mask);
        v2 = v2(sep_mask);
        K1 = K1(sep_mask);
        K2 = K2(sep_mask);

        denom = log(v2 ./ v1);
        valid = isfinite(denom) & (denom ~= 0) & isfinite(log(K2 ./ K1));

        if any(valid)
            beta_pairs = log(K2(valid) ./ K1(valid)) ./ denom(valid);
            alpha_pairs = K1(valid) ./ (v1(valid) .^ beta_pairs);

            % Pure nchoosek: only filter extreme/invalid values, no hardcoded physics limit
            keep = (beta_pairs < 0) & (beta_pairs >= -5) & ...  % Must be negative, not extreme
                (alpha_pairs >= 100) & (alpha_pairs <= 1e5) & ...
                isfinite(beta_pairs) & isfinite(alpha_pairs);

            if any(keep)
                alpha_min = min(alpha_pairs(keep));
                alpha_max = max(alpha_pairs(keep));
                beta_min = min(beta_pairs(keep));
                beta_max = max(beta_pairs(keep));
            end
        end
    end
end

%% Final bounds sanitation (prevents lb > ub infeasibility)

eps_rel = 0.02;    % 2% width if degenerate
eps_abs = 1e-8;    % absolute minimum width

% Helper inline: widen interval if too tight
widen = @(lo,hi) deal( ...
    min(lo,hi), ...
    max(lo,hi) );

% Pipes
for i=1:numel(r_min)
    if ~isfinite(r_min(i)) || r_min(i) <= 0, r_min(i) = 1e-6; end
    if ~isfinite(r_max(i)) || r_max(i) <= 0, r_max(i) = max(10*r_min(i), 1e-3); end
    [r_min(i), r_max(i)] = widen(r_min(i), r_max(i));
    if r_max(i) - r_min(i) < max(eps_abs, eps_rel*abs(r_max(i)))
        mid = 0.5*(r_min(i)+r_max(i));
        r_min(i) = max(1e-6, mid*(1-eps_rel));
        r_max(i) = mid*(1+eps_rel);
    end
end

% Pumps
for i=1:numel(A_min)
    if ~isfinite(A_min(i)) || A_min(i) <= 0, A_min(i) = 1e-6; end
    if ~isfinite(A_max(i)) || A_max(i) <= 0, A_max(i) = max(10*A_min(i), 1); end
    [A_min(i), A_max(i)] = widen(A_min(i), A_max(i));
    if A_max(i) - A_min(i) < max(eps_abs, eps_rel*abs(A_max(i)))
        mid = 0.5*(A_min(i)+A_max(i));
        A_min(i) = max(1e-6, mid*(1-eps_rel));
        A_max(i) = mid*(1+eps_rel);
    end
end

for i=1:numel(B_min)
    if ~isfinite(B_min(i)) || B_min(i) <= 0, B_min(i) = 1e-6; end
    if ~isfinite(B_max(i)) || B_max(i) <= 0, B_max(i) = max(10*B_min(i), 1e-3); end
    [B_min(i), B_max(i)] = widen(B_min(i), B_max(i));
    if B_max(i) - B_min(i) < max(eps_abs, eps_rel*abs(B_max(i)))
        mid = 0.5*(B_min(i)+B_max(i));
        B_min(i) = max(1e-6, mid*(1-eps_rel));
        B_max(i) = mid*(1+eps_rel);
    end
end

% Valve
if exist('alpha_min','var')
    if ~isfinite(alpha_min) || alpha_min <= 0, alpha_min = 1e2; end
    if ~isfinite(alpha_max) || alpha_max <= 0, alpha_max = 1e5; end
    [alpha_min, alpha_max] = widen(alpha_min, alpha_max);
    if alpha_max - alpha_min < max(eps_abs, eps_rel*abs(alpha_max))
        mid = 0.5*(alpha_min+alpha_max);
        alpha_min = max(1e-6, mid*(1-eps_rel));
        alpha_max = mid*(1+eps_rel);
    end
end

if exist('beta_min','var')
    if ~isfinite(beta_min), beta_min = -5; end
    if ~isfinite(beta_max), beta_max = -0.1; end
    [beta_min, beta_max] = widen(beta_min, beta_max);
    if beta_max - beta_min < 1e-6
        mid = 0.5*(beta_min+beta_max);
        beta_min = mid - 0.01;
        beta_max = mid + 0.01;
    end
end

end