clear all
clc
close all

%load('LCP_param.mat')
load('sys_params.mat')
load('system_param2.mat')
addpath C:\Users\alp\Desktop\pathlcp

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

x1=[-5:0.1:5]; x2=[-5:0.1:5];
[X1,X2]=meshgrid(x1,x2);
u = zeros(size(X1,1), size(X1,2));

 for i = 1:size(X1,1)
     if mod(i,50) == 0
         i
     end
     for j = 1:size(X1,2)
         x = [ X1(i,j); X2(i,j); 0; 0; 0];  
         %x = [0; 0; X1(i,j); X2(i,j)];
         lam = pathlcp(Fc,Ec*x + c); 
          x_basis = x; lam_basis = lam;
          inp = D(:,1:20)*lam(1:20); 
          u(i,j) = inp(1); 
     end
 end
 
[dfdx, dfdy] = gradient(u);
% 
% dfdx = dfdx*100; dfdy = dfdy*100;

%dfdx = 10000*rand(61,61);

figure
surf(X1,X2,u, dfdx.^10 + dfdy.^10)
xlabel('x_1', 'FontSize', 24)
ylabel('x_2', 'FontSize', 24)
zlabel('\phi(x)','FontSize', 24)
set(gca,'FontSize',30)