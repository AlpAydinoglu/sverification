clear all
clc
close all

load('LCP_param.mat')

%double integrator dynamics
A2 = [1 1; 0 1];
B2 = [0.5; 1];
D2 = [0; 0];
Ec2 = [0 0];
Fc2 = [1];
c2 = [0];

sys = ss(A2,B2,zeros(2), zeros(2,1),1);
K = lqr(sys, 0.1*eye(2), 1*eye(1));
A2 = A2 - B2*K;

%save values
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

rng(2)
%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts
num_iter = 150;
x = zeros(n,num_iter);
lam = zeros(m, num_iter);
x(:,1) = 10 * (rand(n,1) - 0.5);

%simulate for sanity check
for i = 1:num_iter
    lam(:,i) = pathlcp(Fc,Ec*x(:,i) + c); 
    x(:,i+1) = A*x(:,i) + D*lam(:,i) + cons;
end

%for sanity check
%plot([0:num_iter], x)

%save system values
z = cons;
save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'z')