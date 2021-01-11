%%%Uses the LCP paramaterization of the NN (LCP_param.mat) and creates the complementarity representation of the system (sys_params.mat) as in Section 3.2.

clear all
clc
close all

%load LCP parameterization of the NN
load('LCP_param.mat')

%construct the complementarity model

%physical system parameters (double integrator dynamics)
A2 = [1 1; 0 1];
B2 = [0.5; 1];
D2 = [0; 0];
Ec2 = [0 0];
Fc2 = [1];
c2 = [0];

%LQR parameters
sys = ss(A2,B2,zeros(2), zeros(2,1),1);
K = lqr(sys, 0.1*eye(2), 1*eye(1));
A2 = A2 - B2*K;

%Combine the physical system and LCP parameterization of the NN as in
%Section 3.2
Ecs = Ec;
Fcs = Fc;
cs = c';
Ds = D;
lam0 = pathlcp(Fcs,Ecs*zeros(2,1) + cs);
k = -Ds*lam0;
cons = B2*k;
c=c';
A = A2;
H2 = zeros(20,1);
D = B2*D;

%save the complementarity representation
z = cons;
save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'z')