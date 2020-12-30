clear all
clc
close all

load('system_param.mat')
load('sys_params.mat')
addpath C:\Users\alp\Desktop\pathlcp

num_iter = 100;

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

x1=[-2:0.01:2]; x2=[-2:0.01:2];
[X1,X2]=meshgrid(x1,x2);
V = zeros(size(X1,1), size(X1,2));

rng(1)
for i = 1:size(X1,1)
     if mod(i,50) == 0
         i
     end
     for j = 1:size(X1,2)
         x = [ X1(i,j); X2(i,j); 0; 0];  
         %x = [0; 0; X1(i,j); X2(i,j)];
         lam = pathlcp(Fc,Ec*x + c); 
          x_basis = x; lam_basis = lam;
           V(i,j) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
                    + lam_basis'  * RR * lam_basis + cc1 * x_basis + cc2 * lam_basis + cc3;
     end
end

figure
contour(X1,X2,V, [1 10 20 30 40 50 60 100 150 200 500 1000], 'LineWidth', 3, 'Color', 'k')
set(gca,'FontSize',40)
xlabel('x_1', 'FontSize', 40)
ylabel('x_2', 'FontSize', 40)
 