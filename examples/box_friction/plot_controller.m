clear all
clc
close all

load('LCP_param.mat')
load('system_param2.mat')
addpath C:\Users\alp\Desktop\pathlcp

%double integrator dynamics
A2 = [1 1; 0 1];
B2 = [0.5; 1];
D2 = [0; 0];
Ec2 = [0 0];
Fc2 = [1];
c2 = [0];
%H2 = [0];

sys = ss(A2,B2,zeros(2), zeros(2,1),1);
K = lqr(sys, 0.1*eye(2), 1*eye(1));
A2 = A2 - B2*K;

c=c';
A = A2;
H2 = zeros(20,1);
cons = zeros(size(A,2),1);
k = double(k);
%D = B2*D;

x1 = -4:0.1:4;
x2 = -4:0.1:4;
[X1,X2] = meshgrid(x1,x2);
u = zeros(size(X1,1), size(X1,2));

%simulate 
 for i = 1:size(X1,1)
     i
     for j = 1:size(X1,2)
         x = [X1(i,j); X2(i,j)];
         u(i,j) = D*pathlcp(Fc,Ec*x + c);
     end
 end
 

figure
surf(X1,X2,u)
xlabel('x_1', 'FontSize', 24)
ylabel('x_2', 'FontSize', 24)
zlabel('\phi(x)','FontSize', 24)
set(gca,'FontSize',24)
