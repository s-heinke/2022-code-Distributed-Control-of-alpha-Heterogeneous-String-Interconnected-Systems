%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "Distributed Control of α-Heterogeneous String Interconnected Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

% Load parameters
load('parameters')

Ni=[3 10 20 30 40 50 75 100 200 300]; % number of subsystems
nyi=2;
nui=1;

Nc=7; % centralized controller only up to 75 subsystems
% Construct generalized plant
% z_i=rho_i+0.5u_i 
for k=1:5 % take average computaion time of 5 runs
    for i=1:length(Ni)
        N=Ni(i);
        [Gp, Gp_de, Gp_AG, Gp_AS]=build_plant(p,N);
        %% Controller Synthesis
        sub=1.05;
        % Centralized Controller
        if i<=Nc
            tic
            K0=h2syn(Gp,N*nyi,N*nui);
            t0(i,k)=toc
        end
        % decentralized controller
        tic
        K1s=h2syn(Gp_de,nyi,nui);
        t1(i,k)=toc
        % Systems Interconnected over Arbitrary Graphs
        ny=ones(1,N)*nyi;
        nu=ones(1,N)*nui;
        tic
        [K_AG, gam1]= h2AG(Gp_AG,ny,nu,sub);
        t2(i,k)=toc

        %% alpha-heterogeneous string interconnected systems
        tic
        [K_AS, gam5]=h2AS(Gp_AS,2*ones(1,3),ones(1,3),sub);
        t3(i,k)=toc
    end
end
%%

figure(1)
plot(Ni(1:Nc),sum(t0')./5,'b','LineWidth',1.5);
hold on
grid on
ylim([0,30])
plot(Ni,sum(t1')./5,'r','LineWidth',1.5);
plot(Ni,sum(t2')./5,'color',[0 0.7 0],'LineWidth',1.5);
plot(Ni,sum(t3')./5,'color',[1.0, 0.5, 0],'LineWidth',1.5);
legend('h2syn','h2syn_d','AG','AS')
xlabel('Number of subsystems')
ylabel('Computation time [s]')
