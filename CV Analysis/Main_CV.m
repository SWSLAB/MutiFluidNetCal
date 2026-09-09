%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

clc
clear
close all

%% Network Parameters
Z = [2.9; 21.5; 22.4; 69.5;2.9; 2.9; 21.5; 21.5; 13.6; 13.6; 26.9; 26.9; 26.9; 84.9; 22.4; 69.5]; % Nodes elevations in m
L = [12.7; 41.1; 38.4; 44.8; 25.9; 3.5]*1000;          % Pipeline lengths in m
O_Pipes = [6; 7; 10; 13; 14; 14];                      % J-2, J-3, J-6, J-9, J-10, J-10
D_Pipes = [7; 9; 11; 14; 15; 16];                      % J-3, J-5, J-7, J-10, J-11, J-12
O_Pumps = [1; 2; 9; 11];                               % T-1, T-1, T-2, J-5, J-5, J-7
D_Pumps = [5; 8; 10; 12];                              % J-1, J-1, J-4, J-6, J-6, J-8
O_Valves = [5; 8; 9; 11; 12; 15; 16; 12];              % J-1, J-4, J-5, J-7, J-8, J-11, J-12
D_Valves = [6; 7; 10; 13; 13; 3; 4; 13];               % J-2, J-3, J-6, J-9, J-9, T-3, T-4
Jm = [1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];   % Measured nodes
J=5:16;
T=1:4;
Ts=[1 2];
JT=1:16;
Pipes=1:6;
Pumps=1:4;
Valves=1:8;
PosValves=5;
NonPosValves=setdiff(Valves,PosValves);

%% Read and Clean Database
load('Database.mat');
SelPerc=0.95;
CleanData

%% 10-fold cross-validation
for fluid_choice = 2:3  % 1 = Diesel, 2 = Jet, 3 = Gasoline
    K = 10; rng(1);
    Nf_cv = size(data8{fluid_choice},1);
    perm_f = randperm(Nf_cv);
    edges_cv = round(linspace(0,Nf_cv,K+1));

    fit_data  = cell(K+1,1);   % 1..K = train folds, K+1 = full data
    test_data = cell(K+1,1);   % test-fold companion, empty for the full-data pass
    for k=1:K
        te_idx = perm_f(edges_cv(k)+1:edges_cv(k+1));
        tr_idx = setdiff(1:Nf_cv,te_idx);
        fit_data{k}  = data8{fluid_choice}(tr_idx,:);
        test_data{k} = data8{fluid_choice}(te_idx,:);
    end
    fit_data{K+1} = data8{fluid_choice};   % full-data production fit, no held-out test


    for fold = 1:K

        % Clear variables from previous run
        clear dh_Pipes0 dh_Pumps0 dh_Valves0 Q95 dhmin0 std_dhmin dh_Pipes dh_Pumps dh_Valves
        clear k_fit H_fit P_fit Dhf_fit Dhk_fit Dhg_fit Hobs_sources RHS_fit AgI invAgI
        clear xsolV R2sysV R2fitV objsolV exitflagV deltaV fh CV_results

        % Run and accumulate results per fold
        Main
        cv_xsol{fold} = xsol;
        cv_R2sysV{fold} = R2sysV;
        cv_R2fitV{fold} = R2fitV;
        cv_deltaV{fold} = deltaV;
        cv_xsolV{fold}  = xsolV;
        cv_RMSE_H_train(fold)  = R2sysV(pareto_point_id);
        cv_RMSE_DH_train(fold) = R2fitV(pareto_point_id);

        data_te = test_data{fold};
        Pobs_te=data_te(:,1:16)/1e6; x_v5_te=data_te(:,17); x_v5_te(x_v5_te==0)=0.001;
        I_Pumps_te=data_te(:,18:21)==1; I_Valves_te=data_te(:,22:29)==1;
        Q_Pipes_te=data_te(:,30:35); I_Pipes_te=Q_Pipes_te>0;
        Q_Pumps_te=data_te(:,36:39); Q_Valves_te=data_te(:,40:47);
        gamma_te = mean(data_te(:,48:49),2);
        I_links_te=[I_Pumps_te I_Valves_te I_Pipes_te];
        N_te = size(data_te,1);
        invAgI_te = cell(N_te,1);
        for i=1:N_te, invAgI_te{i}=pinv(Ag.*I_links_te(i,:)'); end
        fun_dh_te=@(x)myobjfun(x,N_te,Pobs_te,x_v5_te,I_Pumps_te,I_Valves_te,Q_Pipes_te,I_Pipes_te,Q_Pumps_te,Q_Valves_te,gamma_te,I_links_te,O_links,D_links,Ag,invAgI_te,Z,L...
            ,O_Pipes,O_Pumps,O_Valves,D_Pipes,D_Pumps,D_Valves,Jm,J,T,Ts,JT,Pipes,Pumps,Valves,PosValves,NonPosValves,0);
        fun_dh_te(xsol);
        cv_RMSE_H_test(fold)  = R2sys;
        cv_RMSE_DH_test(fold) = R2fit;

        fprintf('CV fold %d/%d: RMSE_H train=%.4f test=%.4f, RMSE_DH train=%.4f test=%.4f\n', ...
            fold, K, cv_RMSE_H_train(fold), cv_RMSE_H_test(fold), cv_RMSE_DH_train(fold), cv_RMSE_DH_test(fold));
    end
    
    %% Cross-validation summary table
    nP=length(Pipes); nPu=length(Pumps); nV=length(NonPosValves);
    r_i=1:nP; A_i=nP+1:nP+nPu; B_i=nP+nPu+1:nP+2*nPu;
    al_i=nP+2*nPu+nV+1; be_i=al_i+1; n_i=be_i+1;

    cv_params = zeros(K, nP+2*nPu+3);
    for kk=1:K
        x = cv_xsol{kk};
        cv_params(kk,:) = [x(n_i), x(r_i)', x(A_i)', x(B_i)', x(al_i), x(be_i)];
    end
    param_names = [{'n'}, arrayfun(@(i)sprintf('r_%d',i),1:nP,'UniformOutput',false), ...
        arrayfun(@(i)sprintf('A_%d',i),1:nPu,'UniformOutput',false), arrayfun(@(i)sprintf('B_%d',i),1:nPu,'UniformOutput',false), {'alpha','beta'}];
    cv_table = table((1:K)', cv_RMSE_H_train', cv_RMSE_DH_train', cv_RMSE_H_test', cv_RMSE_DH_test', ...
        'VariableNames', {'fold','RMSE_H_train','RMSE_DH_train','RMSE_H_test','RMSE_DH_test'});
    cv_table = [cv_table, array2table(cv_params,'VariableNames',param_names)];
    disp(cv_table)
    
    %% Store CV results
    CV_results=fluid_results;
    CV_results.cv_table = cv_table;
    CV_results.cv_xsol = cv_xsol;
    CV_results.cv_RMSE_H_train = cv_RMSE_H_train;
    CV_results.cv_RMSE_DH_train = cv_RMSE_DH_train;
    CV_results.cv_RMSE_H_test = cv_RMSE_H_test;
    CV_results.cv_RMSE_DH_test = cv_RMSE_DH_test;
    CV_results.cv_R2sysV = cv_R2sysV;
    CV_results.cv_R2fitV = cv_R2fitV;
    CV_results.cv_deltaV = cv_deltaV;
    CV_results.cv_xsolV  = cv_xsolV;
    
    %% Save results
    filename = sprintf('CV_results_%s.mat', fluid_names{fluid_choice});
    save(filename,'CV_results');
    fprintf('\n✅ Cross validation complete for %s across %d folds. Saved as %s.\n', fluid_names{fluid_choice}, K, filename);

end
