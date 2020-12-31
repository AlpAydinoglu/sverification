% clear all
% clc
% close all

load('LCP_param.mat')

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
D2 = [0 0; 0 0; 0 0; -1/(L*mp) 1/L*mp];
A2 = eye(4) + Ts*A2;
B2 = Ts*B2;
D2 = Ts*D2;
Ec2 = [-1 L 0 0; 1 -L 0 0];
Fc2 = [1/k_1 0; 0 1/k_2];
c2 = [d;d];

sys = ss(A2,B2,zeros(4), zeros(4,1),1);
K = lqr(sys, 10*eye(4), eye(1));
A2 = A2 - B2*K;

%save values
Ecs = Ec;
Fcs = Fc;
cs = c';
Ds = D;

%check x=0
lam0 = pathlcp(Fcs,Ecs*zeros(4,1) + cs);
k = -Ds*lam0;

%system parameters
c=c';
k = double(k);
Ec = [Ec; Ec2];
Fc = blkdiag(Fc, Fc2);
c = [c; c2];
D = [B2*D D2];
A = A2;
cons = B2*k;

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

rng(4)
num_iter = 500;
x = zeros(n,num_iter);
lam = zeros(m, num_iter);
%initial condition
x(:,1) = 3 * (rand(n,1) - 0.5);

%check if x=0 is an equilibrium (sanity check)
for i = 1:num_iter
    lam(:,i) = pathlcp(Fc,Ec*x(:,i) + c);
    x(:,i+1) = A*x(:,i) + D*lam(:,i) + cons;
end 
% figure
% plot(0:num_iter, x, 'LineWidth', 2)

%save system values
save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'cons')