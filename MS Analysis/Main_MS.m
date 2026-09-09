%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

clc
clear
close all

for fluid_choice = 1:3   % 1 = Diesel, 2 = Jet, 3 = Gasoline
    % Solve base run
    Main

    % Stability analysis: mvnrnd-randomized restarts
    CoV = 0.10;
    K_restarts = 2;
    options_stab = optimoptions(options,'Display','off');   % quiet copy, used only in this block

    stab_xsol  = cell(length(deltaV), K_restarts);
    stab_R2sys = nan(length(deltaV), K_restarts);
    stab_R2fit = nan(length(deltaV), K_restarts);
    stab_exit  = nan(length(deltaV), K_restarts);

    for pidx = 1:length(deltaV)
        if pidx==1
            x_center = x000;
            obj_fun_here = fun_h;
            nonlcon_here = [];
        elseif pidx==length(deltaV)
            x_center = x000;
            obj_fun_here = fun_dh;
            nonlcon_here = [];
        else
            x_center = xsolV{pidx-1};
            delta = deltaV(pidx);
            obj_fun_here = fun_dh;
            nonlcon_here = @(x) deal(myobjfun(x,N,Pobs,x_v5,I_Pumps,I_Valves,Q_Pipes,I_Pipes,Q_Pumps,Q_Valves,gamma,I_links,O_links,D_links,Ag,invAgI,Z,L...
                ,O_Pipes,O_Pumps,O_Valves,D_Pipes,D_Pumps,D_Valves,Jm,J,T,Ts,JT,Pipes,Pumps,Valves,PosValves,NonPosValves,1)-delta,[]);
        end

        sigma0 = CoV*abs(x_center);
        sigma0(sigma0==0) = 1e-6;

        r = 1;
        cnt = 1;
        stopflag = 0;
        while stopflag==0
            X0mat(r,:) = mvnrnd(x_center',diag(sigma0.^2),1);
            x0_rand = max(min(X0mat(r,:)',ub),lb);
            [xsol_r,~,exitflag_r] = fmincon(obj_fun_here,x0_rand,[],[],A_eq_n,b_eq_n,lb,ub,nonlcon_here,options_stab);
            fun_dh(xsol_r);
            if exitflag_r>0
                stab_xsol{pidx,cnt}  = xsol_r;
                stab_R2sys(pidx,cnt) = R2sys;
                stab_R2fit(pidx,cnt) = R2fit;
                stab_exit(pidx,cnt)  = exitflag_r;
                cnt = cnt+1;
            end
            fprintf('  restart %d (pt %d/%d): exitflag=%d, R2sys=%.4f [%d/%d stored]\n', ...
                r, pidx, length(deltaV), exitflag_r, R2sys, min(cnt,K_restarts), K_restarts);
            if cnt>K_restarts
                stopflag = 1;
            else
                r = r+1;
            end
        end

        fprintf('Pareto pt %d/%d: R2sys CoV=%.2f%%, R2fit CoV=%.2f%%\n', pidx,length(deltaV), ...
            100*std(stab_R2sys(pidx,:))/abs(mean(stab_R2sys(pidx,:))), ...
            100*std(stab_R2fit(pidx,:))/abs(mean(stab_R2fit(pidx,:))));
    end

    %% Statistical summary table
    stab_table = table((1:length(deltaV))', mean(stab_R2sys,2), std(stab_R2sys,0,2), ...
        mean(stab_R2fit,2), std(stab_R2fit,0,2), sum(stab_exit>0,2), ...
        'VariableNames', {'ParetoIdx','R2sys_mean','R2sys_std','R2fit_mean','R2fit_std','N_converged'});
    disp(stab_table)

    %% Store MS results
    fluid_results.stab_xsol  = stab_xsol;
    fluid_results.stab_R2sys = stab_R2sys;
    fluid_results.stab_R2fit = stab_R2fit;
    fluid_results.stab_exit  = stab_exit;
    fluid_results.stab_table = stab_table;

    %% Save results
    filename = sprintf('MS_results_%s.mat', fluid_names{fluid_choice});
    save(filename,'fluid_results');
    fprintf('\n✅ Multistart complete for %s across %d runs. Saved as %s.\n', fluid_names{fluid_choice}, K_restarts, filename);

end


