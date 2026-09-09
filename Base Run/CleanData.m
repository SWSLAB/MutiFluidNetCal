%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

% Data Cleaning
% Data retention tracking
fluid_names = {'Diesel', 'Jet fuel', 'Gasoline'};
N_combined = nan(6,1);
N_combined(1) = size(data0,1);   % Raw SCADA

% Define Parameters
Pobs0=data0(:,1:16)/1e6;         % Pressure in MPa
x_v50=data0(:,17);
I_Pumps0=data0(:,18:21)==1;
I_Valves0=data0(:,22:29)==1;
Q_Pipes0=data0(:,30:35);
I_Pipes0=Q_Pipes0>0;
Q_Pumps0=data0(:,36:39);
Q_Valves0=data0(:,40:47);
gamma0 = mean(data0(:,48:49),2);

% NaN pressures when nodes are not in flow path
Pobs0(Q_Pumps0(:,1)==0,[1 6])=nan;
Pobs0(Q_Pumps0(:,2)==0,[2 8])=nan;
Pobs0(Q_Pipes0(:,5)==0,15)=nan;
Pobs0(Q_Pipes0(:,6)==0,16)=nan;

% Time diffs in data signals
D=abs(diff([Pobs0 Q_Pipes0 x_v50 gamma0],1,1));
D01=[nan(1,size(D,2));D];
D02=[D;nan(1,size(D,2))];
D0=[D01 D02];

% Filter out non-homogeneous fluid
filter_out_idx0=find(~strcmp(f0(:,1),f0(:,2)));
pct_removed_0 = length(filter_out_idx0)/size(data0,1)*100;
data1=data0;
data1(filter_out_idx0,:)=[];
D1=D0;
D1(filter_out_idx0,:)=[];
f1=f0;
f1(filter_out_idx0,:)=[];
N_combined(2) = size(data1,1);  % After fluid homogeneity
fluid_id1=strcmp(f1(:,1),'Diesel')*1+strcmp(f1(:,1),'Jet')*2+strcmp(f1(:,1),'Gasoline')*3;

% Filter out zero sources and zero destinations
Q_Pipes1=data1(:,30:35);
Q_Pumps1=data1(:,36:39);
filter_out_idx1=unique([find(sum(Q_Pumps1(:,1:2),2)==0);find(sum(Q_Pipes1(:,5:6),2)==0)]);
pct_removed_1 = length(filter_out_idx1)/size(data1,1)*100;
data2=data1;
data2(filter_out_idx1,:)=[];
D2=D1;
D2(filter_out_idx1,:)=[];
fluid_id2=fluid_id1;
fluid_id2(filter_out_idx1,:)=[];
N_combined(3) = size(data2,1);  % After non-operational conditions
figure;hist(fluid_id2)

% Filter out peaks in data signals according to fluid
for f=1:3
    data3{f}=data2(fluid_id2==f,:);
    D3{f}=D2(fluid_id2==f,:);
    flagMat{f}=D3{f};
    for i=1:size(D3{f},2)
        P95(f,i)=quantile(D3{f}(:,i),SelPerc);
        flagMat{f}(~isnan(D3{f}(:,i)),i)=D3{f}(~isnan(D3{f}(:,i)),i)<=P95(f,i);
    end
    flag{f} = all(flagMat{f} == 1 | isnan(flagMat{f}), 2);
    filter_out_idx3{f}=find(flag{f}==0);
    pct_removed_3(f) = length(filter_out_idx3{f})/size(D3{f},1)*100;

    data4{f}= data3{f};
    data4{f}(filter_out_idx3{f},:)=[];
    D4{f}=D3{f};
    D4{f}(filter_out_idx3{f},:)=[];
end

% Filter out nonsteady flow
for f=1:3
    Q_Pipes4{f}=data4{f}(:,30:35);
    Q_Pumps4{f}=data4{f}(:,36:39);
    dQ{f}=abs(sum(Q_Pumps4{f}(:,1:2),2)-sum(Q_Pipes4{f}(:,5:6),2));
    dQ_P95(f)=quantile(dQ{f},SelPerc);
    filter_out_idx4{f}=find(dQ{f}>dQ_P95(f));
    pct_removed_4(f) = length(filter_out_idx4{f})/size(D4{f},1)*100;

    data5{f}= data4{f};
    data5{f}(filter_out_idx4{f},:)=[];
    D5{f}=D4{f};
    D5{f}(filter_out_idx4{f},:)=[];
end
N_combined(4) = sum(cellfun(@(x) size(x,1), data5));  % After fluid flow stability (peak filter + mass-balance nonsteady-flow filter)

% Set pressures to zero when node is not in flow path
for f=1:3
    Q_Pipes5{f}=data5{f}(:,30:35);
    Q_Pumps5{f}=data5{f}(:,36:39);
    Pobs5{f}=data5{f}(:,1:16)/1e6; % Pressure in MPa

    Pobs6{f}=Pobs5{f};
    Pobs6{f}(Q_Pumps5{f}(:,1)==0,[1 6])=0;
    Pobs6{f}(Q_Pumps5{f}(:,2)==0,[2 8])=0;
    Pobs6{f}(Q_Pipes5{f}(:,5)==0,15)=0;
    Pobs6{f}(Q_Pipes5{f}(:,6)==0,16)=0;

    data6{f}= data5{f};
    data6{f}(:,1:16)= Pobs6{f}*1e6;
end

% Check for rows with missing data
Jm = [1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];   % measured nodes
for f=1:3
    tmp=[data6{f}(:,[Jm 17:47])];
    filter_out_idx6{f}=find(any(isnan(tmp),2));
    pct_removed_6(f) = length(filter_out_idx6{f})/size(data6{f},1)*100;

    data7{f}= data6{f};
    data7{f}(filter_out_idx6{f},:)=[];
end
N_combined(5) = sum(cellfun(@(x) size(x,1), data7));  % After missing-data filter (isnan check only)

% Hydraulic feasibility
for f=1:3
    data0 = data7{f};
    N0 = size(data0,1);

    Pobs0=data0(:,1:16)/1e6; % Pressure in MPa
    x_v50=data0(:,17);
    I_Pumps0=data0(:,18:21)==1;
    I_Valves0=data0(:,22:29)==1;
    Q_Pipes0=data0(:,30:35);
    I_Pipes0=Q_Pipes0>0;
    Q_Pumps0=data0(:,36:39);
    Q_Valves0=data0(:,40:47);
    gamma0 = mean(data0(:,48:49),2);
    dh_Pipes0=[]; dh_Pumps0=[];dh_Valves0=[];

    Hobs0=repmat(Z',N0,1)+1e6*Pobs0./repmat(gamma0,1,length(JT));
    for i=1:length(Pipes)
        dh_Pipes0(:,i)=Hobs0(:,O_Pipes(i))-Hobs0(:,D_Pipes(i));
    end
    for i=1:length(Pumps)
        if i~=1
            dh_Pumps0(:,i)=Hobs0(:,D_Pumps(i))-Hobs0(:,O_Pumps(i));
        else
            dh_Pumps0(:,i)=Hobs0(:,D_Pumps(i)+1)-Hobs0(:,O_Pumps(i));
        end
        Q95(i)=quantile(Q_Pumps0(:,i),SelPerc);
        dhmin0(i)=min(dh_Pumps0(Q_Pumps0(:,i)>=Q95(i),i));
        std_dhmin(i)=std(dh_Pumps0(Q_Pumps0(:,i)>=Q95(i),i));
    end
    for i=1:length(Valves)
        dh_Valves0(:,i)=Hobs0(:,O_Valves(i))-Hobs0(:,D_Valves(i));
    end

    % Physical constraints for valve 5 - Literature-based fully-open K_Q
    Cv_design = 130;                    % [gpm] design Cv at full opening
    K_unitless_full_open = 2.5;         % Dimensionless minor loss coefficient at full opening [-]
    g = 9.81;                           % [m/s^2]

    % Cv–K–d relation (Crane): K = 891 * d_inch^4 / Cv^2
    d_valve_inch = (K_unitless_full_open * Cv_design^2 / 891)^0.25;  % [inch]
    d_valve_m    = d_valve_inch * 0.0254;                            % [m]

    % Q-based coefficient for h_L = K_Q_full_open * Q^2, Q in m^3/s, headloss in m
    K_Q_full_open = (8 * K_unitless_full_open) / (pi^2 * g * d_valve_m^4);

    % SCADA-based k5 for valve 5
    k5 = dh_Valves0(:,5) ./ (Q_Valves0(:,5).^2);

    % Physical filter for valve 5
    badK_valve5 = k5 < K_Q_full_open & I_Valves0(:,5)==1;

    fout = any(I_Pipes0==1 & dh_Pipes0<0, 2) | ...
        any(I_Pumps0==1 & dh_Pumps0<0, 2) | ...
        any(I_Pumps0==1 & dh_Pumps0<(dhmin0-std_dhmin), 2) | ...
        badK_valve5;

    data0(fout,:) = [];
    data8{f}=data0;
end

% Per-fluid retention breakdown, stages 2-6
N_combined(6) = sum(cellfun(@(x) size(x,1), data8));                          % After missing-data filter
N_fluid_26 = nan(3,4);                                                        % rows = Diesel/Jet/Gasoline, cols = stages 2:5
N_fluid_26(:,1) = [sum(fluid_id1==1); sum(fluid_id1==2); sum(fluid_id1==3)];  % stage 2
N_fluid_26(:,2) = [sum(fluid_id2==1); sum(fluid_id2==2); sum(fluid_id2==3)];  % stage 3
N_fluid_26(:,3) = cellfun(@(x) size(x,1), data5)';                            % stage 4
N_fluid_26(:,4) = cellfun(@(x) size(x,1), data7)';                            % stage 5
N_fluid_26(:,5) = cellfun(@(x) size(x,1), data8)';                            % stage 6

%% Operational envelope of clean data
% Distribution comparison: before vs after filtering
for f=1:3
    flow_before{f}=sum(data3{f}(:,36:37),2)*3600;
    flow_after{f}=sum(data8{f}(:,36:37),2)*3600;
    dens_before{f}=mean(data3{f}(:,48:49),2);
    dens_after{f}=mean(data8{f}(:,48:49),2);
    v5_before{f}=data3{f}(:,17);
    v5_after{f}=data8{f}(:,17);
end
figure
for f=1:3
    subplot(3,3,f)
    hist(flow_before{f},40)
    hold on
    hist(flow_after{f},40)
    title(fluid_names{f})
    xlabel('Flow (m3/h)')
    subplot(3,3,3+f)
    hist(dens_before{f},40)
    hold on
    hist(dens_after{f},40)
    xlabel('Density')
    xlim([7000 8500])
    subplot(3,3,6+f)
    hist(v5_before{f},40)
    hold on
    hist(v5_after{f},40)
    xlabel('Valve V5 position')
    xlim([0 1])
end
prc=[5 25 50 75 95 99];
edges_flow = 0:50:350;
edges_dens = 7000:250:8500;
edges_v5 = 0:0.1:1;
for f=1:3
    fprintf('\n%s:\n', fluid_names{f});
    fprintf('  flow before=%.1f after=%.1f, density before=%.1f after=%.1f, v5 before=%.3f after=%.3f\n', ...
        mean(flow_before{f},'omitnan'), mean(flow_after{f},'omitnan'), ...
        mean(dens_before{f},'omitnan'), mean(dens_after{f},'omitnan'), ...
        mean(v5_before{f},'omitnan'), mean(v5_after{f},'omitnan'));
    fprintf('  flow percentiles before: ');    fprintf('%.1f ', quantile(flow_before{f},prc/100)); fprintf('\n');
    fprintf('  flow percentiles after:  ');    fprintf('%.1f ', quantile(flow_after{f},prc/100));  fprintf('\n');
    fprintf('  dens percentiles before: ');    fprintf('%.1f ', quantile(dens_before{f},prc/100)); fprintf('\n');
    fprintf('  dens percentiles after:  ');    fprintf('%.1f ', quantile(dens_after{f},prc/100));  fprintf('\n');
    fprintf('  v5 percentiles before:   ');    fprintf('%.3f ', quantile(v5_before{f},prc/100));   fprintf('\n');
    fprintf('  v5 percentiles after:    ');    fprintf('%.3f ', quantile(v5_after{f},prc/100));    fprintf('\n');
    fprintf('  --- operational envelope (%% of calibration data by range) ---\n');
    cnt_flow = histcounts(flow_after{f}, edges_flow);
    pct_flow = cnt_flow/sum(cnt_flow)*100;
    fprintf('  Flow (m3/h):\n');
    for b=1:length(edges_flow)-1
        fprintf('    [%d, %d): %.1f%%\n', edges_flow(b), edges_flow(b+1), pct_flow(b));
    end
    cnt_dens = histcounts(dens_after{f}, edges_dens);
    pct_dens = cnt_dens/sum(cnt_dens)*100;
    fprintf('  Density:\n');
    for b=1:length(edges_dens)-1
        fprintf('    [%d, %d): %.1f%%\n', edges_dens(b), edges_dens(b+1), pct_dens(b));
    end
    cnt_v5 = histcounts(v5_after{f}, edges_v5);
    pct_v5 = cnt_v5/sum(cnt_v5)*100;
    fprintf('  Valve V5 position:\n');
    for b=1:length(edges_v5)-1
        fprintf('    [%.1f, %.1f): %.1f%%\n', edges_v5(b), edges_v5(b+1), pct_v5(b));
    end
end
