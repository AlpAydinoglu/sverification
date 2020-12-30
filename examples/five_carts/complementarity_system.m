clear all
clc
close all

load('LCP_param.mat')

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
sys = ss(A2,B2,zeros(n), zeros(n,k),1);
K = lqr(sys, 1*eye(n), 1*eye(k));
A2 = A2 - B2*K;
c=c';
Ec = [Ec; Ec2];
Fc = blkdiag(Fc, Fc2);
c = [c; c2];
D = [B2*D D2];
A = A2;
cons = zeros(n,1);
z = cons;

%save system values
save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'z')