%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

%clc
%clear
%close all

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
%SelPerc=0.95;
CleanData

%% Run Optimization on clean dataset
%fluid_choice = 3;                   % 1 = Diesel, 2 = Jet, 3 = Gasoline

data = data8{fluid_choice};
N_selected_after = size(data8{fluid_choice},1);

Pobs=data(:,1:16)/1e6;              % Pressure in MPa
x_v5=data(:,17);
x_v5(x_v5==0)=0.001;                % To prevent numerical issues.
I_Pumps=data(:,18:21)==1;
I_Valves=data(:,22:29)==1;
Q_Pipes=data(:,30:35);
I_Pipes=Q_Pipes>0;
Q_Pumps=data(:,36:39);
Q_Valves=data(:,40:47);
gamma = mean(data(:,48:49),2);

%% Build graph matrix
I_links=[I_Pumps I_Valves I_Pipes];
O_links=[O_Pumps;O_Valves;O_Pipes];
D_links=[D_Pumps;D_Valves;D_Pipes];
OD=[O_links D_links];
Ag=zeros(size(OD,1),length(JT));
for i=1:size(OD,1)
    Ag(i,OD(i,:))=[-1 1];
end
Ag(:,1:2)=[];

%% Classic curve fitting
N = size(data,1);
Hobs=repmat(Z',N,1)+1e6*Pobs./repmat(gamma,1,length(JT));
for i=1:length(Pipes)
    dh_Pipes(:,i)=Hobs(:,O_Pipes(i))-Hobs(:,D_Pipes(i));
    X = log(Q_Pipes(I_Pipes(:,i)==1,i));
    Y = log(dh_Pipes(I_Pipes(:,i)==1,i));
    % Perform linear regression
    p = polyfit(X, Y, 1);
    n_fit(i,1) = p(1);
    R_fit(i,1) = exp(p(2));
end
for i=1:length(Pumps)
    if i~=1
        dh_Pumps(:,i)=Hobs(:,D_Pumps(i))-Hobs(:,O_Pumps(i));
    else
        dh_Pumps(:,i)=Hobs(:,D_Pumps(i)+1)-Hobs(:,O_Pumps(i));
    end
    X = Q_Pumps(I_Pumps(:,i)==1,i).^2;
    Y = dh_Pumps(I_Pumps(:,i)==1,i);
    % Perform linear regression
    p = polyfit(X, Y, 1);
    B_fit(i,1) = -p(1);
    A_fit(i,1) = p(2);
end
for i=PosValves
    dh_Valves(:,i)=Hobs(:,O_Valves(i))-Hobs(:,D_Valves(i));
    X = log(x_v5(I_Valves(:,i)==1));
    Y = log(dh_Valves(I_Valves(:,i)==1,i)./Q_Valves(I_Valves(:,i)==1,i).^2);
    % Perform linear regression
    p = polyfit(X, Y, 1);
    beta_fit = p(1);
    alpha_fit = exp(p(2));
end

%% Plotting
plotting_flag=0;
if plotting_flag==1
    close all
    for i=1:length(Pipes)
        figure
        try % In case of empty data for pipe
            xx = linspace(min(Q_Pipes(I_Pipes(:,i)==1,i)), max(Q_Pipes(I_Pipes(:,i)==1,i)), 200);
            yy_fit = R_fit(i) * xx.^n_fit(i);
            fh(i)=axes;
            plot(fh(i),Q_Pipes(I_Pipes(:,i)==1,i), dh_Pipes(I_Pipes(:,i)==1,i), 'o', xx, yy_fit, 'b-','LineWidth',1.5);
            hold all
            xlabel('Q'); ylabel('dH');
            legend('Data','Fit');
            title(['Pipe ' num2str(i)])
            grid on;
        catch ME
            title(['Pipe ' num2str(i)])
        end
    end

    for i=1:length(Pumps)
        figure
        try % In case of empty data for pump
            xx = linspace(min(Q_Pumps(I_Pumps(:,i)==1,i)), max(Q_Pumps(I_Pumps(:,i)==1,i)), 200);
            yy_fit = A_fit(i)-B_fit(i)* xx.^2;
            fh(i+6)=axes;
            plot(fh(i+6),Q_Pumps(I_Pumps(:,i)==1,i), dh_Pumps(I_Pumps(:,i)==1,i), 'o', xx, yy_fit, 'b-','LineWidth',1.5);
            hold all
            xlabel('Q'); ylabel('dH');
            legend('Data','Fit');
            title(['Pump Station ' num2str(i)])
            grid on;
        catch ME
            title(['Pump Station ' num2str(i)])
        end
    end

    for i=PosValves
        figure
        fh(11)=axes;
        xx = linspace(min(x_v5(I_Valves(:,i)==1)), max(x_v5(I_Valves(:,i)==1)), 200);
        yy_fit = alpha_fit*xx.^beta_fit;
        plot(fh(11),x_v5(I_Valves(:,i)==1),dh_Valves(I_Valves(:,i)==1,i)./Q_Valves(I_Valves(:,i)==1,i).^2, 'o', xx, yy_fit, 'b-','LineWidth',1.5);
        xlabel('x'); ylabel('K');
        hold all
        legend('Data','Fit');
        title(['Valve ' num2str(i)])
        grid on;
    end
end

%% Evaluate the classic fitting for nodal pressure simulation
k0_fit=zeros(length(NonPosValves),1);
k_fit(:,NonPosValves)=repmat(k0_fit',N,1);
k_fit(:,PosValves)=alpha_fit*x_v5.^beta_fit;

Dhf_fit=-repmat(R_fit',N,1).*Q_Pipes.^repmat(n_fit',N,1);
Dhk_fit=-k_fit.*Q_Valves.^2;
Dhg_fit=repmat(A_fit',N,1)-repmat(B_fit',N,1).*Q_Pumps.^2;
Hobs_sources=Z(1:2)'+1e6*Pobs(:,1:2)./gamma;
RHS_fit=([Dhg_fit';Dhk_fit';Dhf_fit']+[Hobs_sources';zeros(size(Ag,1)-2,N)]).*I_links';

for i=1:N
    AgI{i}=Ag.*I_links(i,:)';
    invAgI{i}=pinv(AgI{i});
    H_fit(i,:)=[Hobs_sources(i,:) (invAgI{i}*RHS_fit(:,i))'];
    H_fit(i,find(all(invAgI{i}==0,2))+2)=Hobs(i,find(all(invAgI{i}==0,2))+2);    % free nodes (nodes that do not appear in the dataset) are set to observed values
end
P_fit=(H_fit-repmat(Z',N,1)).*repmat(gamma,1,length(JT))/1e6;

%% Plotting
if plotting_flag==1
    figure
    fh(12)=axes;
    hold on
    plot(fh(12),JT,H_fit(:,[1 2 5:16 3:4]),'-b')
    xlabel('Node'); ylabel('Head (m)');
    title('Head Diagram from fitting')
    grid on;
    xticks(1:16);
    xticklabels(["T"+(1:2), "J"+(1:12),"T"+(3:4)])

    xlim([0 16])
    ylim([0 1000])
    figure
    fh(13)=axes;
    hold on
    plot(fh(13),JT,Hobs(:,[1 2 5:16 3:4]),'-k')
    xlabel('Node'); ylabel('Head (m)');
    title('Observed Head Diagram')
    grid on;
    xticks(1:16);
    xticklabels(["T"+(1:2), "J"+(1:12),"T"+(3:4)])

    xlim([0 16])
    ylim([0 1000])
end
drawnow

%% System Hydraulic Estimation
options=optimoptions('fmincon','Display','iter');
options.MaxFunEvals=1e6;
options.MaxIter=5000;

% Initial guess
x0{1,1}=R_fit./L;
x0{2,1}=A_fit;
x0{3,1}=B_fit;
x0{4,1}=k0_fit;
x0{5,1}=alpha_fit;
x0{6,1}=beta_fit;
x0{7,1}=n_fit;
x00=cell2mat(x0);

%% Bounds
[r_min, r_max, A_min, A_max, B_min, B_max, alpha_min, alpha_max, beta_min, beta_max, n_min, n_max] = optimization_bounds(Hobs, Q_Pipes, ...
    Q_Pumps, Q_Valves, x_v5, I_Pipes, I_Pumps, I_Valves, L, O_Pipes, D_Pipes, O_Pumps, D_Pumps, O_Valves, D_Valves, Pipes, Pumps, PosValves, R_fit, A_fit, B_fit);

k0_min = zeros(length(NonPosValves), 1);
k0_max = zeros(length(NonPosValves), 1);

lb0{1,1} = r_min;
lb0{2,1} = A_min;
lb0{3,1} = B_min;
lb0{4,1} = k0_min;
lb0{5,1} = alpha_min;
lb0{6,1} = beta_min;
lb0{7,1} = n_min;
lb = cell2mat(lb0);

ub0{1,1} = r_max;
ub0{2,1} = A_max;
ub0{3,1} = B_max;
ub0{4,1} = k0_max;
ub0{5,1} = alpha_max;
ub0{6,1} = beta_max;
ub0{7,1} = n_max;
ub = cell2mat(ub0);

x000 = max(min(x00, ub), lb);  % Clip x00 to be within bounds

% Equal n constraint
A_eq_n=zeros(length(Pipes)-1,length(x00));
b_eq_n=zeros(length(Pipes)-1,1);
A_eq_n(:,length(x00)-length(Pipes)+1)=1;
for i=length(x00)-length(Pipes)+2:length(x00)
    A_eq_n(i-(length(x00)-length(Pipes)+2)+1,i)=-1;
end

%% Epsilon-constraint Pareto sweep (f1=R2fit=RMSE_dh minimized, f2=R2sys=RMSE_h constrained)
global R2sys R2fit P H

% Solve two single-objective optimization problems to get the range of the Pareto
fun_dh=@(x)myobjfun(x,N,Pobs,x_v5,I_Pumps,I_Valves,Q_Pipes,I_Pipes,Q_Pumps,Q_Valves,gamma,I_links,O_links,D_links,Ag,invAgI,Z,L...
    ,O_Pipes,O_Pumps,O_Valves,D_Pipes,D_Pumps,D_Valves,Jm,J,T,Ts,JT,Pipes,Pumps,Valves,PosValves,NonPosValves,0);
fun_h=@(x)myobjfun(x,N,Pobs,x_v5,I_Pumps,I_Valves,Q_Pipes,I_Pipes,Q_Pumps,Q_Valves,gamma,I_links,O_links,D_links,Ag,invAgI,Z,L...
    ,O_Pipes,O_Pumps,O_Valves,D_Pipes,D_Pumps,D_Valves,Jm,J,T,Ts,JT,Pipes,Pumps,Valves,PosValves,NonPosValves,1);

[x_extreme_h,~]=fmincon(fun_h,x000,[],[],A_eq_n,b_eq_n,lb,ub,[],options);
fun_h(x_extreme_h);
f2_min=R2sys;
[x_extreme_dh,~]=fmincon(fun_dh,x000,[],[],A_eq_n,b_eq_n,lb,ub,[],options);
fun_dh(x_extreme_dh);
f2_max=R2sys;

xsolV{1}=x_extreme_h;
xsolV{10}=x_extreme_dh;
fun_dh(xsolV{1});
R2sysV(1)=R2sys;
R2fitV(1)=R2fit;
objsolV{1}=R2fit;

fun_dh(xsolV{10});
R2sysV(10)=R2sys;
R2fitV(10)=R2fit;
objsolV{10}=R2fit;
exitflagV(1)=1;
exitflagV(10)=1;

fprintf('\nEpsilon-constraint range: R2sys in [%.4f, %.4f]\n',f2_min,f2_max);

% Solve the inner Pareto points
deltaV=linspace(f2_min,f2_max,10);

for cnt=2:9
    delta=deltaV(cnt);
    nonlcon=@(x) deal(myobjfun(x,N,Pobs,x_v5,I_Pumps,I_Valves,Q_Pipes,I_Pipes,Q_Pumps,Q_Valves,gamma,I_links,O_links,D_links,Ag,invAgI,Z,L...
        ,O_Pipes,O_Pumps,O_Valves,D_Pipes,D_Pumps,D_Valves,Jm,J,T,Ts,JT,Pipes,Pumps,Valves,PosValves,NonPosValves,1)-delta,[]);
    [xsolV{cnt},objsolV{cnt},exitflagV(cnt)] = ...
        fmincon(fun_dh,xsolV{cnt-1},[],[],A_eq_n,b_eq_n,lb,ub,nonlcon,options);
    if exitflagV(cnt)<=0
        fprintf('⚠️ delta=%.4f (idx %d): exitflag=%d — check convergence\n',delta,cnt,exitflagV(cnt));
    end
    fun_dh(xsolV{cnt});
    R2sysV(cnt)=R2sys;
    R2fitV(cnt)=R2fit;
end

%% Auto-detect knee point using perpendicular distance method
% Normalize objectives to [0,1] for fair comparison
x_norm = (R2sysV - min(R2sysV)) ./ max(eps, max(R2sysV) - min(R2sysV));
y_norm = (R2fitV - min(R2fitV)) ./ max(eps, max(R2fitV) - min(R2fitV));

% Endpoints of Pareto front
p_start = [x_norm(1), y_norm(1)];
p_end = [x_norm(end), y_norm(end)];

% Perpendicular distance from each point to line connecting endpoints
dx = p_end(1) - p_start(1);
dy = p_end(2) - p_start(2);
line_length = sqrt(dx^2 + dy^2);
line_length = max(line_length, eps);  % Protect against coincident endpoints

numerator = abs(dy*x_norm - dx*y_norm + p_end(1)*p_start(2) - p_end(2)*p_start(1));
distances = numerator / line_length;
[~, pareto_point_id] = max(distances);

% Apply selection
delta_selected=deltaV(pareto_point_id);
xsol=xsolV{pareto_point_id};
r_sys=xsol(1:length(Pipes));
R_sys=L.*r_sys;
A_sys=xsol(length(Pipes)+1:length(Pumps)+length(Pipes));
B_sys=xsol(length(Pumps)+length(Pipes)+1:length(Pumps)+length(Pipes)+length(Pumps));
k0_sys=xsol(length(Pumps)+length(Pipes)+length(Pumps)+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves));
alpha_sys=xsol(length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1);
beta_sys=xsol(length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1);
n_sys=xsol(length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1+length(Pipes));
fun_dh(xsol);
P_sys=P;
H_sys=H;

% Save classical fitting performance
fun_dh(x00);
R2sys_classical = R2sys;  % Save classical fitting H domain RMSE
R2fit_classical = R2fit;  % Save classical fitting DH domain RMSE

%% Plotting
if exist('fh','var') && length(fh)>=12 && all(arrayfun(@ishandle, fh(1:12)))   % Plot on existing figures

    for i=1:6
        try
            xx = linspace(min(Q_Pipes(I_Pipes(:,i)==1,i)), max(Q_Pipes(I_Pipes(:,i)==1,i)), 200);
            yy_sys = R_sys(i) * xx.^n_sys(i);
            plot(fh(i), xx, yy_sys, 'r-','LineWidth',1.5);
        catch ME
        end
    end

    for i=1:4
        try
            xx = linspace(min(Q_Pumps(I_Pumps(:,i)==1,i)), max(Q_Pumps(I_Pumps(:,i)==1,i)), 200);
            yy_sys = A_sys(i)-B_sys(i)* xx.^2;
            plot(fh(i+6), xx, yy_sys, 'r-','LineWidth',1.5);
        catch ME
        end
    end

    for i=PosValves
        xx = linspace(min(x_v5(I_Valves(:,i)==1)), max(x_v5(I_Valves(:,i)==1)), 200);
        yy_sys = alpha_sys*xx.^beta_sys;
        plot(fh(11), xx, yy_sys, 'r-','LineWidth',1.5);
    end

    plot(fh(12),JT,H_sys(:,[1 2 5:16 3:4]),'-r')
    plot(fh(12),JT,Hobs(:,[1 2 5:16 3:4]),'-g')
end

if plotting_flag==1
    figure
    for i=3:13
        subplot(2,11,i-2)
        hist(abs(Pobs(:,Jm(i))-P_sys(:,Jm(i))))
        title(['sys-J' num2str(Jm(i)-4)])
        xlim([0 0.7])
    end

    for i=3:13
        subplot(2,11,11+i-2)
        hist(abs(Pobs(:,Jm(i))-P_fit(:,Jm(i))))
        title(['fit-J' num2str(Jm(i)-4)])
        xlim([0 0.7])
    end

    figure
    plot(R2sysV,R2fitV,'o-','LineWidth',2,'MarkerSize',8);
    hold on
    plot(R2sys_classical,R2fit_classical,'kx','MarkerSize',12,'LineWidth',2);  % Classical fitting solution
    plot(R2sysV(pareto_point_id),R2fitV(pareto_point_id),'rs','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','r');  % Selected pareto point
    xlabel('RMSE H Domain (m)'); ylabel('RMSE DH Domain (m)');
    legend('Pareto Front','Initial Guess','Selected Solution','Location','best');
    grid on;
end

%% Print optimized parameters

fprintf('\nSelected Pareto: idx=%d (delta=%.4f), RMSE_H=%.2f m, RMSE_DH=%.2f m', pareto_point_id, deltaV(pareto_point_id), R2sysV(pareto_point_id), R2fitV(pareto_point_id));

fprintf('\nPipe parameters:\n');
for i=1:length(Pipes)
    fprintf('Pipe %d: r=%.4f, n=%.3f, R=%.1f (fit: n=%.3f, R=%.1f)\n', i, r_sys(i), n_sys(i), R_sys(i), n_fit(i), R_fit(i));
end

fprintf('\nPump parameters:\n');
for i=1:length(Pumps)
    fprintf('Pump %d: A=%.0f, B=%.0f (fit: A=%.0f, B=%.0f)\n', i, A_sys(i), B_sys(i), A_fit(i), B_fit(i));
end

fprintf('\nValve parameters:\n');
fprintf('Valve V%d: alpha=%.0f, beta=%.3f (fit: alpha=%.0f, beta=%.3f)\n', PosValves, alpha_sys, beta_sys, alpha_fit, beta_fit);

for i=1:length(NonPosValves)
    fprintf('Valve V%d: k0=%.4f\n', NonPosValves(i), k0_sys(i));
end

%% Plotting
if plotting_flag==1
    for i=1:6
        try
            xx = linspace(min(Q_Pipes(I_Pipes(:,i)==1,i)), max(Q_Pipes(I_Pipes(:,i)==1,i)), 200);
            yy_sys = R_sys(i) * xx.^n_sys(i);
            plot(fh(i), xx, yy_sys, 'r-','LineWidth',1.5);
        catch ME
        end
    end

    for i=1:4
        try
            xx = linspace(min(Q_Pumps(I_Pumps(:,i)==1,i)), max(Q_Pumps(I_Pumps(:,i)==1,i)), 200);
            yy_sys = A_sys(i)-B_sys(i)* xx.^2;
            plot(fh(i+6), xx, yy_sys, 'r-','LineWidth',1.5);
        catch ME
        end
    end

    for i=PosValves
        xx = linspace(min(x_v5(I_Valves(:,i)==1)), max(x_v5(I_Valves(:,i)==1)), 200);
        yy_sys = alpha_sys*xx.^beta_sys;
        plot(fh(11), xx, yy_sys, 'r-','LineWidth',1.5);
    end

    plot(fh(12),JT,H_sys(:,[1 2 5:16 3:4]),'-r')
    plot(fh(12),JT,Hobs(:,[1 2 5:16 3:4]),'-g')

    figure
    for i=3:13
        subplot(2,11,i-2)
        hist(abs(Pobs(:,Jm(i))-P_sys(:,Jm(i))))
        title(['sys-J' num2str(Jm(i)-4)])
        xlim([0 0.7])
    end

    for i=3:13
        subplot(2,11,11+i-2)
        hist(abs(Pobs(:,Jm(i))-P_fit(:,Jm(i))))
        title(['fit-J' num2str(Jm(i)-4)])
        xlim([0 0.7])
    end
end
%% Save results for current fluid
fluid_results.R_fit = R_fit;
fluid_results.n_fit = n_fit;
fluid_results.A_fit = A_fit;
fluid_results.B_fit = B_fit;
fluid_results.alpha_fit = alpha_fit;
fluid_results.beta_fit = beta_fit;
fluid_results.R_sys = R_sys;
fluid_results.n_sys = n_sys;
fluid_results.A_sys = A_sys;
fluid_results.B_sys = B_sys;
fluid_results.alpha_sys = alpha_sys;
fluid_results.beta_sys = beta_sys;
fluid_results.Q_Pipes = Q_Pipes;
fluid_results.I_Pipes = I_Pipes;
fluid_results.dh_Pipes = dh_Pipes;
fluid_results.Q_Pumps = Q_Pumps;
fluid_results.I_Pumps = I_Pumps;
fluid_results.dh_Pumps = dh_Pumps;
fluid_results.Q_Valves = Q_Valves;
fluid_results.I_Valves = I_Valves;
fluid_results.dh_Valves = dh_Valves;
fluid_results.x_v5 = x_v5;
fluid_results.r_min = r_min;
fluid_results.r_max = r_max;
fluid_results.A_min = A_min;
fluid_results.A_max = A_max;
fluid_results.B_min = B_min;
fluid_results.B_max = B_max;
fluid_results.alpha_min = alpha_min;
fluid_results.alpha_max = alpha_max;
fluid_results.beta_min = beta_min;
fluid_results.beta_max = beta_max;
fluid_results.n_min = n_min;
fluid_results.n_max = n_max;
fluid_results.R2sysV = R2sysV;
fluid_results.R2fitV = R2fitV;
fluid_results.R2sys_classical = R2sys_classical;
fluid_results.R2fit_classical = R2fit_classical;
fluid_results.pareto_point_id = pareto_point_id;
fluid_results.Pipes = Pipes;
fluid_results.Pumps = Pumps;
fluid_results.PosValves = PosValves;
fluid_results.L = L;
fluid_results.Z = Z;
fluid_results.N_combined = N_combined;
fluid_results.N_fluid_25 = N_fluid_26;
fluid_results.N_selected_after = N_selected_after;
fluid_results.d_valve5_inch = d_valve_inch;
fluid_results.d_valve5_m = d_valve_m;
fluid_results.KQ_valve5_full_open = K_Q_full_open;
fluid_results.K_unitless_valve5 = K_unitless_full_open;

%% Save results to a file
% if fluid_choice==1
%     save('fluid1.mat','fluid_results');
%     fprintf('\n✅ Saved Diesel results to fluid1.mat\n');
% elseif fluid_choice==2
%     save('fluid2.mat','fluid_results');
%     fprintf('\n✅ Saved Jet results to fluid2.mat\n');
% elseif fluid_choice==3
%     save('fluid3.mat','fluid_results');
%     fprintf('\n✅ Saved Gasoline results to fluid3.mat\n');
% end

