%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "Distributed Control of α-Heterogeneous String Interconnected Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

% Load parameters
load('parameters')
N=3;    % No. of subsystems
nui=1;  % No. of inputs per subsystems
nyi=2;  % No. of outputs per subsystems

% Construct generalized plant
% z_i=rho_i+0.5u_i 
[Gp, Gp_de, Gp_AG, Gp_AS]=build_plant(p,N);
%% Controller Synthesis
% Centralized Controller
K0=h2syn(Gp,N*nyi,N*nui);
% decentralized controller
K1s=h2syn(Gp_de,nyi,nui);
K1=blkdiag(K1s,K1s,K1s);
% Systems Interconnected over Arbitrary Graphs
sub=1.05;  % solve for a five percent suboptimal controller due to numerical issues
ny=ones(1,N)*nyi;
nu=ones(1,N)*nui;
[K_AG, gam1]= h2AG(Gp_AG,ny,nu,sub);
K2 = AG2MIMO(K_AG);
% alpha-heterogeneous string interconnected systems
[K_AS, gam5]=h2AS(Gp_AS,2*ones(1,3),ones(1,3),sub);
K3=AS2MIMO(K_AS);

%% Comparison
CL0 = lft(Gp, K0); 
disp(['H2 performance h2syn: ', num2str(norm(CL0,2),3)])

CL0 = lft(Gp, K1); 
disp(['H2 performance h2syn_d: ', num2str(norm(CL0,2),4)])

CL0 = lft(Gp, K2); 
disp(['H2 performance h2AG: ', num2str(norm(CL0,2),3), '         Bound: ', num2str(gam1,4)])

CL0 = lft(Gp, K3); 
disp(['H2 performance h2AS: ', num2str(norm(CL0,2),3), '        Bound: ', num2str(gam5,4)])

%% Comparison in time domain
% disturbance of [0.5 -1 0.2] Nm at time t=1s and duration of 0.5s.
K=K0;
t_stop=8;
sim('closed_loop',t_stop);
load u_sim
load y_sim
y_limit=[-20 20];
u_limit=[-5 5];
figure(1)
subplot(2,4,1)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('h2syn')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    legend('\rho_1','\rho_2','\rho_3')
    ylabel('angle [°]')
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,4,5)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('h2syn')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    legend('u_1','u_2','u_3')
    ylabel('voltage [V]')
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
K=K1;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,4,2)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('h2syn_d')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    %legend('\rho_1','\rho_2','\rho_3')
    %ylabel('angle [°]')
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,4,6)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('h2syn_d')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    %legend('u_1','u_2','u_3')
    %ylabel('voltage [V]')
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
K=K2;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,4,3)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AG')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    %legend('\rho_1','\rho_2','\rho_3')
    %ylabel('angle [°]')
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,4,7)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AG')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    %legend('u_1','u_2','u_3')
    %ylabel('voltage [V]')
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
K=K3;
sim('closed_loop',t_stop);
load u_sim
load y_sim
figure(1)
subplot(2,4,4)
    plot(y_sim(1,:),(180/pi)*y_sim(2,:),'b')
    title('AS')
    hold on
    grid on
    plot(y_sim(1,:),(180/pi)*y_sim(4,:),'r')
    plot(y_sim(1,:),(180/pi)*y_sim(6,:),'color',[0 0.7 0])
    %legend('\rho_1','\rho_2','\rho_3')
    %ylabel('angle [°]')
    ylim(y_limit)
    xlim([0 t_stop])
subplot(2,4,8)
    plot(u_sim(1,:),u_sim(2,:),'b')
    title('AS')
    hold on
    grid on
    plot(u_sim(1,:),u_sim(3,:),'r')
    plot(u_sim(1,:),u_sim(4,:),'color',[0 0.7 0])
    %legend('u_1','u_2','u_3')
    %ylabel('voltage [V]')
    xlabel('time [s]')
    ylim(u_limit)
    xlim([0 t_stop])
 