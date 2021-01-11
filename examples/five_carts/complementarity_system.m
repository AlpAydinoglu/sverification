%%%Uses the LCP paramaterization of the NN (LCP_param.mat) and creates the complementarity representation of the system (sys_params.mat) as in Section 3.2.

clear all
clc
close all

%load LCP parameterization of the NN
load('LCP_param.mat')

%physical system parameters 
n = 5;
m = n-1;
k = n;
alpha = 1;
A2 = eye(n);
B2 = eye(n);
D2 = [-eye(m); zeros(1,m)]+[zeros(1,m); eye(m)];
Ec2 = [-eye(m) zeros(m,1)] + [zeros(m,1) eye(m)];
Fc2 = alpha*eye(m);
c2 = zeros(m,1);
H2 = zeros(m,n);

%LQR parameters
sys = ss(A2,B2,zeros(n), zeros(n,k),1);
K = lqr(sys, 1*eye(n), 1*eye(k));
A2 = A2 - B2*K;

%Combine the physical system and LCP parameterization of the NN as in
%Section 3.2
c=c';
Ec = [Ec; Ec2];
Fc = blkdiag(Fc, Fc2);
c = [c; c2];
D = [B2*D D2];
A = A2;
cons = zeros(n,1);
z = cons;

%save the complementarity representation
save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'z')