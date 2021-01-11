%%%Plots the sublevel sets of the Lyapunov function (Figure 4) with few trajectories on top. (requires sys_params.mat and lyap_params.mat)

clear all
clc
close all

load('lyap_params.mat')
load('sys_params.mat')

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

x1=[-4:0.1:4]; x2=[-4:0.1:4];
[X1,X2]=meshgrid(x1,x2);
V = zeros(size(X1,1), size(X1,2));
ls = zeros(size(X1,1), size(X1,2));

rng(1)
for i = 1:size(X1,1)
     if mod(i,50) == 0
         i
     end
     for j = 1:size(X1,2)
         x = [ X1(i,j); X2(i,j)];  
         lam = pathlcp(Fc,Ec*x + c); 
          x_basis = x; lam_basis = lam;
           V(i,j) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
                    + lam_basis'  * RR * lam_basis + cc1' * x_basis + cc2' * lam_basis + cc3;
     end
end

figure
contour(X1,X2,V, 20, 'LineWidth', 3, 'Color','k')
%surf(X1,X2,V)
hold on
contour(X1,X2,V,[170,170], 'LineWidth', 7, 'Color','b')  %480 480 oncesi
hold on

num_iter = 150;

    xh = zeros(n,num_iter);
    lamh = zeros(m, num_iter);
    %initial condition
    xh(:,1) = [-1;3];
         for i = 1:num_iter
             lamh(:,i) = pathlcp(Fc,Ec*xh(:,i) + c); 
             xh(:,i+1) = A*xh(:,i) + D*lamh(:,i) + z;
         end
    clrs = rand(3,1)';
    plot(xh(1,:), xh(2,:), '--', 'LineWidth',4, 'Color', 'r')
    hold on
    scatter(xh(1,:), xh(2,:),100, 'r', 'filled')
    
hold on

num_iter = 150;

    xh = zeros(n,num_iter);
    lamh = zeros(m, num_iter);
    %initial condition
    xh(:,1) = [-3;1];
         for i = 1:num_iter
             lamh(:,i) = pathlcp(Fc,Ec*xh(:,i) + c); 
             xh(:,i+1) = A*xh(:,i) + D*lamh(:,i) + z;
         end
    clrs = rand(3,1)';
    plot(xh(1,:), xh(2,:), '--', 'LineWidth',4, 'Color', 'm')
    hold on
    scatter(xh(1,:), xh(2,:),100, 'm', 'filled')

hold on

num_iter = 150;

    xh = zeros(n,num_iter);
    lamh = zeros(m, num_iter);
    %initial condition
    xh(:,1) = [0;-3];
         for i = 1:num_iter
             lamh(:,i) = pathlcp(Fc,Ec*xh(:,i) + c); 
             xh(:,i+1) = A*xh(:,i) + D*lamh(:,i) + z;
         end
    clrs = rand(3,1)';
    plot(xh(1,:), xh(2,:), '--', 'LineWidth',4, 'Color', [0 0.5 0])
    hold on
    scatter(xh(1,:), xh(2,:),100, [0 0.5 0], 'filled')
    
 num_iter = 150;

    xh = zeros(n,num_iter);
    lamh = zeros(m, num_iter);
    %initial condition
    xh(:,1) = [2.5;-1];
         for i = 1:num_iter
             lamh(:,i) = pathlcp(Fc,Ec*xh(:,i) + c); 
             xh(:,i+1) = A*xh(:,i) + D*lamh(:,i) + z;
         end
    clrs = rand(3,1)';
    plot(xh(1,:), xh(2,:), '--', 'LineWidth',4, 'Color', [.5 0 .5])
    hold on
    scatter(xh(1,:), xh(2,:),100, [.5 0 .5], 'filled')   
    
    
set(gca,'FontSize',40)
xlabel('x_1', 'FontSize', 40)
ylabel('x_2', 'FontSize', 40)



 