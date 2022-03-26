function [Gp, Gp_de, Gp_AG, Gp_AS] = build_plant(p,N)
    % This file constructs the generalized plants for the interconnected 
    % seesaw-cart system.
    % Author: Simon Heinke
    
    % cost function
    Dzu=0.5;
    Ctz=[1 0 0 0];

    %% 1) Interconnected System (Arbitrary Graph)
    % Subsystem 1 and 3
    Att1=[0 0 1 0; ...
        0 0 0 1; ...
          p.a*p.mA*p.g/p.J-(p.k*(p.L^2)/(4*p.J)) -p.mC*p.g/p.J -p.cappa/p.J -p.c*p.Ke*p.Kt/(p.J*p.Ra*(p.p^2)*(p.r^2)); ...
          p.c*p.a*p.mA*p.g/p.J-p.g-p.c*p.k*p.L^2/(4*p.J) -p.c*p.mC*p.g/p.J -p.c*p.cappa/p.J -(p.c^2/p.J+1/p.mC)*p.Ke*p.Kt/(p.Ra*p.p^2*p.r^2)];
    Ats1=[0; 0; -p.k*(p.L^2)/(4*p.J); -p.c*p.k*p.L^2/(4*p.J)];
    Ast1=[1 0 0 0];
    Ass1=0;
    Btd1=-[0; 0; 1/p.J; p.c/p.J];
    Bsd1=0;
    Btu1=[0; 0; p.c*p.Kt/(p.J*p.Ra*p.p*p.r); (p.c^2/p.J+1/p.mC)*p.Kt/(p.Ra*p.p*p.r)];
    Bsu1=0;
    Ctz1=Ctz;
    Csz1=0;
    Dzd1=0;
    Dzu1=Dzu;
    Cty1=[eye(2) zeros(2)];
    Csy1=zeros(2,1);
    Dyd1=zeros(2,1);
    Dyu1=zeros(2,1);
    
    % Subsystem 2
    Att2=[0 0 1 0; ...
        0 0 0 1; ...
          p.a*p.mA*p.g/p.J-(p.k*(p.L^2)/(2*p.J)) -p.mC*p.g/p.J -p.cappa/p.J -p.c*p.Ke*p.Kt/(p.J*p.Ra*(p.p^2)*(p.r^2)); ...
          p.c*p.a*p.mA*p.g/p.J-p.g-p.c*p.k*p.L^2/(2*p.J) -p.c*p.mC*p.g/p.J -p.c*p.cappa/p.J -(p.c^2/p.J+1/p.mC)*p.Ke*p.Kt/(p.Ra*p.p^2*p.r^2)];
      
    Ats2=[0 0; 0 0; -p.k*(p.L^2)/(4*p.J) -p.k*(p.L^2)/(4*p.J); -p.c*p.k*p.L^2/(4*p.J) -p.c*p.k*p.L^2/(4*p.J)];
    Ast2=[1 0 0 0; 1 0 0 0];
    Ass2=zeros(2);
    Btd2=-[0; 0; 1/p.J; p.c/p.J];
    Bsd2=zeros(2,1);
    Btu2=[0; 0; p.c*p.Kt/(p.J*p.Ra*p.p*p.r); (p.c^2/p.J+1/p.mC)*p.Kt/(p.Ra*p.p*p.r)];
    Bsu2=zeros(2,1);
    Ctz2=Ctz;
    Csz2=zeros(1,2);
    Dzu2=Dzu;
    Cty2=[eye(2) zeros(2)];        
    Csy2=zeros(2);
    Dyd2=zeros(2,1);
    Dyu2=zeros(2,1);
    Dzd2=0;

    G1.A=[Att1 Ats1; Ast1 Ass1];
    G1.B=[Btd1 Btu1; Bsd1 Bsu1];
    G1.C=[Ctz1 Csz1; Cty1 Csy1];
    G1.D=[Dzd1 Dzu1; Dyd1 Dyu1];                     

    G2.A=[Att2 Ats2; Ast2 Ass2];
    G2.B=[Btd2 Btu2; Bsd2 Bsu2];
    G2.C=[Ctz2 Csz2; Cty2 Csy2];
    G2.D=[Dzd2 Dzu2; Dyd2 Dyu2];                         

    Gp_AG.Sub={G1, G2};                
    Gp_AG.Int=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  
    Ord=2*ones(1,N); Ord(1)=1; Ord(N)=1;
    Gp_AG.Ord=Ord;                 
    Gp_AG.Ts=0;
    
    %% 2) %% alpha-heterogeneous string interconnected systems
    Gp_AS=Gp_AG;
    Gp_AS.Sub{3}=Gp_AS.Sub{1};
    Gp_AS.Ord(N)=3;
    Gp_AS.ns=1;
    Gp_AS.Nki=ones(1,3);
    Gp_AS.Nki(2)=N-2;

    %% 3) centralized
    Gp=AG2MIMO(Gp_AG,2*ones(1,N),ones(1,N));
    %% 5) Decentralized
    A=[0 0 1 0; ...
       0 0 0 1; ...
       p.a*p.mA*p.g/p.J    -p.mC*p.g/p.J    -p.cappa/p.J    -p.c*p.Ke*p.Kt/(p.J*p.Ra*(p.p^2)*(p.r^2)); ...
       p.c*p.a*p.mA*p.g/p.J-p.g    -p.c*p.mC*p.g/p.J   -p.c*p.cappa/p.J    -(p.c^2/p.J+1/p.mC)*p.Ke*p.Kt/(p.Ra*p.p^2*p.r^2)];

    Gp_de=ss(A,[Btd1 Btu1],[Ctz; Cty1],[Dzd1 Dzu; Dyd1 Dyu1]);
    
end


