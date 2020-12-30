clear all
clc
close all

load('LCP_param.mat')
load('system_param.mat')
addpath C:\Users\alp\Desktop\pathlcp

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

c=c';
k = double(k);
Ec = [Ec; Ec2];
Fc = blkdiag(Fc, Fc2);
c = [c; c2];
hold = D;
D = [B2*D D2];
A = A2;

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
         x = [ X1(i,j); X2(i,j); 0; 0];  
         %x = [0; 0; X1(i,j); X2(i,j)];
         lam = pathlcp(Fc,Ec*x + c); 
          x_basis = x; lam_basis = lam;
          u(i,j) = hold*lam(1:20); 
     end
 end
 
[dfdx, dfdy] = gradient(u);

dfdx = dfdx*100; dfdy = dfdy*100;

%dfdx = 10000*rand(61,61);

figure
surf(X1,X2,u, dfdx.^10 + dfdy.^10)
xlabel('x_1', 'FontSize', 24)
ylabel('x_2', 'FontSize', 24)
zlabel('\phi(x)','FontSize', 24)
set(gca,'FontSize',30)