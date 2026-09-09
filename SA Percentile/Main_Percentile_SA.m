%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

clc
clear
close all

%% Sensitivity-analysis loop setup
threshold_vector = 0.90:0.01:0.99;
for fluid_choice = 1:3  % 1 = Diesel, 2 = Jet, 3 = Gasoline
    for th_id = 1:length(threshold_vector)

        % Clear variables from previous run
        clear dh_Pipes0 dh_Pumps0 dh_Valves0 Q95 dhmin0 std_dhmin dh_Pipes dh_Pumps dh_Valves
        clear k_fit H_fit P_fit Dhf_fit Dhk_fit Dhg_fit Hobs_sources RHS_fit AgI invAgI
        clear xsolV R2sysV R2fitV objsolV exitflagV deltaV fh fluid_results

        % Run SA threshold iteration
        SelPerc=threshold_vector(th_id);
        Main

        % Accumulate results per threshold
        fluid_results.threshold=SelPerc;
        sens_results(th_id)=fluid_results;

    end

    %% Save results
    filename = sprintf('SA_results_%s.mat', fluid_names{fluid_choice});
    save(filename,'sens_results','fluid_choice','threshold_vector');
    fprintf('\n✅ Sensitivity analysis complete for %s across %d thresholds. Saved as %s.\n', fluid_names{fluid_choice}, length(threshold_vector), filename);
end
