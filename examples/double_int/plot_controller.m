clear all
clc
close all

load('LCP_param.mat')
load('lyap_params.mat')

c=c';
k=-0.0147;

x1 = -4:0.1:4;
x2 = -4:0.1:4;
[X1,X2] = meshgrid(x1,x2);
u = zeros(size(X1,1), size(X1,2));

%simulate 
 for i = 1:size(X1,1)
     i
     for j = 1:size(X1,2)
         x = [X1(i,j); X2(i,j)];
         u(i,j) = D*pathlcp(Fc,Ec*x + c) + k;
     end
 end
 

figure
surf(X1,X2,u)
xlabel('x_1', 'FontSize', 40)
ylabel('x_2', 'FontSize', 40)
zlabel('\phi(x)','FontSize', 40)
set(gca,'FontSize',40)
