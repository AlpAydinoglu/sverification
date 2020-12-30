clear all
clc
close all

load('LCP_param.mat')
addpath 'C:\Users\alp1a\pathlcp\pathmexw64'
%addpath 'C:\Users\alp1a\OneDrive\Masaüstü\lcp_solve_trial'
%addpath C:\Users\alp\Desktop\pathlcp

%1 if you found the Lyapunov function (its parameters), zero otherwise
Lyap_func = 1;
if Lyap_func == 1
    load('system_param.mat')
else
    disp('Simulation without the Lyapunov function, saving system parameters (A,D,Ec,Fc,c,cons)')
end

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

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

rng(4)
num_iter = 150;
x = zeros(n,num_iter);
lam = zeros(m, num_iter);
%initial condition
x(:,1) = 3 * (rand(n,1) - 0.5);

%find equilibrium (check if x=0 is an equilibrium)
    for i = 1:num_iter
        lam(:,i) = pathlcp(Fc,Ec*x(:,i) + c);
        x(:,i+1) = A*x(:,i) + D*lam(:,i) + B2*k;
    end
    
 if Lyap_func == 0
     figure
     plot(0:num_iter, x, 'LineWidth', 2)
 end
    
%shift equilibrium to zero
% c = c + Ec*x(:,end);
% cons = B2*k + A*x(:,end) - x(:,end);
cons = B2*k;

if Lyap_func == 1

    
figure

trials = 1000;
num_iter = 80;
rng(1)

for i = 1:trials
    
    xh{i} = zeros(n,num_iter);
    lamh{i} = zeros(m, num_iter);
    xh{i}(:,1) = 3 * (rand(n,1) - 0.5);
    V{i} = zeros(1,num_iter);
    for j = 1:num_iter
            lamh{i}(:,j) = pathlcp(Fc,Ec*xh{i}(:,j) + c);
            xh{i}(:,j+1) = A*xh{i}(:,j) + D*lamh{i}(:,j) + cons;
            x_basis = xh{i}(:,j); lam_basis = lamh{i}(:,j);
           V{i}(j) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
                   + lam_basis'  * RR * lam_basis + cc1 * x_basis + cc2 * lam_basis + cc3;
    end
    
end

mx = zeros(1,num_iter);
st_x_min = zeros(n,num_iter+1);
st_x_max = zeros(n,num_iter+1);
for i = 1:trials
    mx = max(mx,V{i});
    st_x_min = min(st_x_min, xh{i});
    st_x_max = max(st_x_max, xh{i});
end

t = 0:num_iter;
t_lyap = 0:num_iter-1;
subplot(1,2,1)
plot(t, min(st_x_min), 'LineWidth', 4, 'Color', [.7 .7 .7])
hold on
plot(t, max(st_x_max), 'LineWidth', 4,'Color', [.7 .7 .7])
hold on
patch([t fliplr(t)], [min(st_x_min) fliplr(max(st_x_max))], [.5 .5 .5], 'FaceAlpha', 0.2)
hold on
plot(t, xh{1},'LineWidth', 2, 'Color', 'k')
xlabel('iteration(k)', 'FontSize', 40)
ylabel('\{x_k\}', 'FontSize', 40)
set(gca,'FontSize',40)
xlim([0 80])
subplot(1,2,2)
plot(0:num_iter-1, V{1}, 'LineWidth', 2, 'Color', 'k')
hold on
patch([t_lyap fliplr(t_lyap)], [zeros(1,num_iter) fliplr(mx)], [.5 .5 .5], 'FaceAlpha', 0.2)
hold on
plot(t_lyap, mx','LineWidth', 4,'Color', [.7 .7 .7])
xlabel('iteration(k)', 'FontSize', 40)
ylabel('V(x_k, \lambda_k)', 'FontSize', 40)
set(gca,'FontSize',40)
xlim([0 80])
    
    
    
    
    
    
%     V = zeros(1,num_iter);
%     Vdot = zeros(1,num_iter);
%     for i = 1:num_iter
%             lam(:,i) = pathlcp(Fc,Ec*x(:,i) + c); %lam(1:20,i) = zeros(20,1); %without the nn controller
%             x(:,i+1) = A*x(:,i) + D*lam(:,i) + cons;
%             x_basis = x(:,i); lam_basis = lam(:,i);
%             V(i) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
%                     + lam_basis'  * RR * lam_basis + cc1 * x_basis + cc2 * lam_basis + cc3;
%              x_plus = x(:,i+1);
%              lamd_basis = pathlcp(Fc,Ec*x_plus + c);
%              Vplus =  x_plus' * PP * x_plus + 2 * x_plus' * QQ  * lamd_basis ...
%                     + lamd_basis'  * RR * lamd_basis + cc1 * x_plus + cc2 * lamd_basis + cc3;
%              Vdot(i) = Vplus - V(i);
%     end

%     %plots
%     plam = lam(end-1:end,:);
%     figure
%     plot(0:num_iter, x, 'LineWidth', 2)
%     % figure
%     % plot(0:num_iter-1, plam,'LineWidth', 2)
%     figure
%     plot(0:num_iter-1, V, 'LineWidth', 2)
%     % figure
%     % semilogy(0:num_iter-1, V, 'LineWidth', 2)
%     figure
%     plot(0:num_iter-1, Vdot, 'LineWidth', 2)

% figure
% subplot(1,2,1)
% plot(0:num_iter, x, 'LineWidth', 4)
% legend({'x_1','x_2','x_3', 'x_4'},'FontSize', 24)
% xlabel('iteration(k)', 'FontSize', 24)
% ylabel('{x_k}', 'FontSize', 24)
% set(gca,'FontSize',20)
% subplot(1,2,2)
% %plot(0:num_iter-1, V, 'LineWidth', 4)
% semilogy(0:num_iter-1, V, 'LineWidth', 4)
% xlabel('iteration(k)', 'FontSize', 24)
% %ylabel('V(x_k, \lambda_k)', 'FontSize', 24)
% ylabel('log V(x_k, \lambda_k)', 'FontSize', 24)
% set(gca,'FontSize',20)

end

save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'cons')