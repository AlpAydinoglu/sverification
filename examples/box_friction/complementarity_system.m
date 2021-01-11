%%%Uses the LCP paramaterization of the NN (LCP_param.mat) and creates the complementarity representation of the system (sys_params.mat) as in Section 3.2.

clear all
clc
close all

%load LCP parameterization of the NN
load('LCP_param.mat')

%construct the complementarity model

%physical system parameters
Ts = 0.1;
A2 = [1 Ts;0 1];
B2 = [0;1]; B2 = B2 * Ts;
D2 = [0 0 0; 1 -1 0]; D2 = D2 * Ts;
Ec2 = [0 1;0 -1;0 0];
Fc2 = [1 -1 1; -1 1 1; -1 -1 0];
c2 = [0;0;0.981];
H2 = [1;-1;0];

%LQR controller parameters
sys = ss(A2,B2,zeros(2), zeros(2,1),1);
K = lqr(sys, 10*eye(2), 1*eye(1));
A2 = A2 - B2*K;
Ec2 = Ec2 - H2*K;

%Combine the physical system and LCP parameterization of the NN as in
%Section 3.2
c=c';
k = double(k);
Ec = [Ec; Ec2];
Fc = blkdiag(Fc, Fc2);
c = [c; c2]; 
sD = [D zeros(1,3)];
D = [B2*D D2];
A = A2;
H2 = [zeros(20,1); H2];
H = H2*sD;
cons = zeros(size(A,2),1);
Fc = Fc+H;

%save the complementarity representation
save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'cons')