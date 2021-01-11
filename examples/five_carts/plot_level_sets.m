%%%Plots the sublevel sets of the Lyapunov function (Figure 12) with few trajectories on top. (requires sys_params.mat and lyap_params.mat)

clear all
clc
close all

load('sys_params.mat')
load('lyap_params.mat')

num_iter = 100;

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

x1=[-2:0.01:2]; x2=[-2:0.1:2];
[X1,X2]=meshgrid(x1,x2);
V = zeros(size(X1,1), size(X1,2));

for i = 1:size(X1,1)
     if mod(i,50) == 0
         i
     end
     for j = 1:size(X1,2)
         a = X1(i,j); b = X2(i,j);
         x = [ 0; a; 0; b; 0];  
         lam = pathlcp(Fc,Ec*x + c); 
          x_basis = x; lam_basis = lam;
           V(i,j) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
                    + lam_basis'  * RR * lam_basis + cc1' * x_basis + cc2' * lam_basis + cc3;
     end
end

figure
subplot(1,2,1)
contour(X1,X2,V, 'LineWidth', 3, 'Color', 'k')
set(gca,'FontSize',40)
xlabel('x_2', 'FontSize', 40)
ylabel('x_4', 'FontSize', 40)

V = zeros(size(X1,1), size(X1,2));

for i = 1:size(X1,1)
     if mod(i,50) == 0
         i
     end
     for j = 1:size(X1,2)
         a = X1(i,j); b = X2(i,j);
         x = [ a; 0; b; 0; 0];  
         lam = pathlcp(Fc,Ec*x + c); 
          x_basis = x; lam_basis = lam;
           V(i,j) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
                    + lam_basis'  * RR * lam_basis + cc1' * x_basis + cc2' * lam_basis + cc3;
     end
end

subplot(1,2,2)
contour(X1,X2,V, 'LineWidth', 3, 'Color', 'k')
set(gca,'FontSize',40)
xlabel('x_3', 'FontSize', 40)
ylabel('x_5', 'FontSize', 40)
 