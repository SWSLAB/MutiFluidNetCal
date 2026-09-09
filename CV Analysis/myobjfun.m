%% ####################################################################################################################
% Code for the paper:
% SCADA-Based Multiobjective Calibration of Hydraulic Models for Multi-Fluid Petroleum Pipelines
% By Alon Mandel and Mashor Housh
% University of Haifa
%% ####################################################################################################################

function [obj]=myobjfun(x,N,Pobs,x_v5,I_Pumps,I_Valves,Q_Pipes,I_Pipes,Q_Pumps,Q_Valves,gamma,I_links,O_links,D_links,Ag,invAgI,Z,L...
    ,O_Pipes,O_Pumps,O_Valves,D_Pipes,D_Pumps,D_Valves,Jm,J,T,Ts,JT,Pipes,Pumps,Valves,PosValves,NonPosValves,W)
global P H R2sys R2fit

%% Formulation
% Independent Decision Variables
r=x(1:length(Pipes));
A=x(length(Pipes)+1:length(Pumps)+length(Pipes));
B=x(length(Pumps)+length(Pipes)+1:length(Pumps)+length(Pipes)+length(Pumps));
k0=x(length(Pumps)+length(Pipes)+length(Pumps)+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves));
alpha=x(length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1);
beta=x(length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1);
n=x(length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1+1:length(Pumps)+length(Pipes)+length(Pumps)+length(NonPosValves)+1+1+length(Pipes));

%% H Domain
Hobs=Z'+1e6*Pobs./gamma;
Hobs_sources=Hobs(:,Ts);

Dhf=-r'.*L'.*Q_Pipes.^(n');
Dhk(:,NonPosValves)=-k0'.*Q_Valves(:,NonPosValves).^2;
Dhk(:,PosValves)=-alpha*x_v5.^beta.*Q_Valves(:,PosValves).^2;
Dhg=A'-B'.*Q_Pumps.^2;

RHS=([Dhg';Dhk';Dhf']+[Hobs_sources';zeros(size(Ag,1)-2,N)]).*I_links';

H=zeros(N,size(Ag,2)+2);
for i=1:N
    H(i,:)=[Hobs_sources(i,:) (invAgI{i}*RHS(:,i))'];
    H(i,find(all(invAgI{i}==0,2))+2)=Hobs(i,find(all(invAgI{i}==0,2))+2);    % Free nodes (nodes that do not appear in the dataset) are set to observed values
end
P=(H-Z').*gamma/1e6;

e2Hm=((H(:,Jm)-Hobs(:,Jm))./(Hobs(:,Jm)+1e-6)).^2;
e2H=e2Hm(:);

%% DH Domain
for i=1:6
    dh_Pipes_obs{i}=Hobs(I_Pipes(:,i)==1,O_Pipes(i))-Hobs(I_Pipes(:,i)==1,D_Pipes(i));
    dh_Pipes{i}=H(I_Pipes(:,i)==1,O_Pipes(i))-H(I_Pipes(:,i)==1,D_Pipes(i));
    if ~isempty(dh_Pipes_obs{i})
        R2pipes(i)= R2(dh_Pipes_obs{i},dh_Pipes{i});
    end
end

for i=1:4
    if i~=1
        dh_Pumps_obs{i}=Hobs(I_Pumps(:,i)==1,D_Pumps(i))-Hobs(I_Pumps(:,i)==1,O_Pumps(i));
        dh_Pumps{i}=H(I_Pumps(:,i)==1,D_Pumps(i))-H(I_Pumps(:,i)==1,O_Pumps(i));
    else
        dh_Pumps_obs{i}=Hobs(I_Pumps(:,i)==1,D_Pumps(i)+1)-Hobs(I_Pumps(:,i)==1,O_Pumps(i));
        dh_Pumps{i}=H(I_Pumps(:,i)==1,D_Pumps(i)+1)-H(I_Pumps(:,i)==1,O_Pumps(i));
    end
    if ~isempty(dh_Pumps_obs{i})
        R2pumps(i)= R2(dh_Pumps_obs{i},dh_Pumps{i});
    end
end

for i=PosValves
    dh_Valve5_obs=Hobs(I_Valves(:,i)==1,O_Valves(i))-Hobs(I_Valves(:,i)==1,D_Valves(i));
    dh_Valve5=H(I_Valves(:,i)==1,O_Valves(i))-H(I_Valves(:,i)==1,D_Valves(i));
end
R2valve5= R2(dh_Valve5_obs,dh_Valve5);

cnt=1;
for j=Jm(3:end)
    R2_Heads(cnt)=R2(Hobs(:,j),H(:,j));
    cnt=cnt+1;
end

%% Define objective
R2sys=sqrt(mean(R2_Heads));                      % RMSE H
R2fit=sqrt(mean([R2pipes R2pumps R2valve5]));    % RMSE dH
total_R2=W*R2sys+(1-W)*R2fit;
obj=total_R2;

if isnan(obj)
    error('nan in obj')
end
end

function [MSE] = R2(Yobs,Ypred)
%R2=1 - sum((Yobs-Ypred).^2)/(sum((Yobs-mean(Yobs)).^2)+1e-6);
MSE=(mean((Yobs-Ypred).^2));
end
