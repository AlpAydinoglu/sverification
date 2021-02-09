%%%Uses the LCP paramaterization of the NN (LCP_param.mat) and creates the complementarity representation of the system (sys_params.mat) as in Section 3.2.

clear all
clc
close all

%load LCP parameterization of the NN
load('LCP_param.mat')

%construct the complementarity model

%physical system parameters
mc = 1;
mp = 1;
L = 1;
d = 1;
k_1 = 1;
k_2 = 1;
g = 9.81;
Ts = 0.1;
A2 = [0 0 1 0; 0 0 0 1; 0 g*mp/mc 0 0; 0 g*(mc+mp)/(L*mc) 0 0];
B2 = [0; 0; 1/mc; 1/(L*mc)];
D2 = [0 0; 0 0; 0 0; 1/(L*mp) -1/L*mp];
A2 = eye(4) + Ts*A2;
B2 = Ts*B2;
D2 = Ts*D2;
Ec2 = [-1 L 0 0; 1 -L 0 0];
Fc2 = [1/k_1 0; 0 1/k_2];
c2 = [d;d];

%LQR controller parameters
sys = ss(A2,B2,zeros(4), zeros(4,1),1);
K = lqr(sys, 10*eye(4), eye(1));
A2 = A2 - B2*K;

%holder for system parameters
Ecs = Ec;
Fcs = Fc;
cs = c';
Ds = D;
lam0 = pathlcp(Fcs,Ecs*zeros(4,1) + cs);
k = -Ds*lam0;

%Combine the physical system and LCP parameterization of the NN as in
%Section 3.2
c=c';
k = double(k);
Ec = [Ec; Ec2];
Fc = blkdiag(Fc, Fc2);
c = [c; c2];
D = [B2*D D2];
A = A2;
cons = B2*k;

%save the complementarity representation
save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'cons')