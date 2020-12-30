clear all
clc
close all

%addpath 'C:\Users\alp1a\pathlcp\pathmexw64'
%addpath 'C:\Users\alp1a\OneDrive\Masaüstü\lcp_solve_trial'
addpath C:\Users\alp\Desktop\pathlcp

load('system_param2.mat')
load('LCP_param.mat')

n = 5;
m = n-1;
k = n;
alpha = 1;

A2 = eye(n);
B2 = eye(n);
D2 = [-eye(m); zeros(1,m)]+[zeros(1,m); eye(m)];
Ec2 = [-eye(m) zeros(m,1)] + [zeros(m,1) eye(m)];
Fc2 = alpha*eye(m);
c2 = zeros(m,1);
H2 = zeros(m,n);

sys = ss(A2,B2,zeros(n), zeros(n,k),1);
K = lqr(sys, 1*eye(n), 1*eye(k));
A2 = A2 - B2*K;
%Ec2 = Ec2 - H2*K;

c=c';
Ec = [Ec; Ec2];
Fc = blkdiag(Fc, Fc2);
c = [c; c2];
D = [B2*D D2];
A = A2;
cons = zeros(n,1);


% A=A2;
% D = D2;
% cons = zeros(n,1);
% c = c2;
% H = H2;
% Ec = Ec2;
% Fc = Fc2;
% B = B2;


%rng(1)
%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts
num_iter = 15;
x = zeros(n,num_iter);
lam = zeros(m, num_iter);
x(:,1) = 10 * (rand(n,1) - 0.5);


N = 10;
%simulate 
%         for i = 1:num_iter
%             i
%             u = mixed_integer_CS(x(:,i), N, A, B, D, cons, Ec, Fc, c, H); 
%             lam(:,i) = pathlcp(Fc,Ec*x(:,i) + c + H*u);
%             x(:,i+1) = A*x(:,i) + D*lam(:,i) + cons + B * u;
%         end
%      


% % %simulate 
        for i = 1:num_iter
            lam(:,i) = pathlcp(Fc,Ec*x(:,i) + c);
            x(:,i+1) = A*x(:,i) + D*lam(:,i) + cons;
        end
 
 x(:,end)' * x(:,end)

% c = c + Ec*x(:,end); 
% cons = B2*k + A*x(:,end) - x(:,end);

figure

trials = 1000;
num_iter = 6;
rng(1)

for i = 1:trials
    
    xh{i} = zeros(n,num_iter);  
    lamh{i} = zeros(m, num_iter); 
    xh{i}(:,1) = 5 * (rand(n,1) - 0.5);
    V{i} = zeros(1,num_iter);
    for j = 1:num_iter
            lamh{i}(:,j) = pathlcp(Fc,Ec*xh{i}(:,j) + c);
            xh{i}(:,j+1) = A*xh{i}(:,j) + D*lamh{i}(:,j) + cons;
            x_basis = xh{i}(:,j); lam_basis = lamh{i}(:,j);
           V{i}(j) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
                   + lam_basis'  * RR * lam_basis + cc1' * x_basis + cc2' * lam_basis + cc3;
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
subplot(1,2,2)
plot(0:num_iter-1, V{3}, 'LineWidth', 2, 'Color', 'k')
hold on
patch([t_lyap fliplr(t_lyap)], [zeros(1,num_iter) fliplr(mx)], [.5 .5 .5], 'FaceAlpha', 0.2)
hold on
plot(t_lyap, mx','LineWidth', 4,'Color', [.7 .7 .7])
xlabel('iteration(k)', 'FontSize', 40)
ylabel('V(x_k, \lambda_k)', 'FontSize', 40)
set(gca,'FontSize',40)
xlim([0 5])




% V = zeros(1,num_iter);
% Vdot = zeros(1,num_iter);
% for i = 1:num_iter
%         lam(:,i) = pathlcp(Fc,Ec*x(:,i) + c);
%         x(:,i+1) = A*x(:,i) + D*lam(:,i) + cons;
%         x_basis = x(:,i); lam_basis = lam(:,i);
%        V(i) = x_basis' * PP * x_basis + 2 * x_basis' * QQ  * lam_basis ...
%                + lam_basis'  * RR * lam_basis + cc1' * x_basis + cc2' * lam_basis + cc3;
% end

% plam = lam(end-1:end,:);
% figure
% plot(0:num_iter, x, 'LineWidth', 2)
% % figure
% % plot(0:num_iter-1, plam,'LineWidth', 2)
% figure
% plot(0:num_iter-1, V, 'LineWidth', 2)
% % figure
% % semilogy(0:num_iter-1, V, 'LineWidth', 2)
% % figure
% % plot(0:num_iter-1, Vdot, 'LineWidth', 2)
% % %x(:,end)

% figure
% subplot(1,2,1)
% plot(0:num_iter, x, 'LineWidth', 4)
% legend({'x_1','x_2','x_3', 'x_4'},'FontSize', 24)
% xlabel('iteration(k)', 'FontSize', 24)
% ylabel('{x_k}', 'FontSize', 24)
% set(gca,'FontSize',20)
% subplot(1,2,2)
% plot(0:num_iter-1, V, 'LineWidth', 4)
% %semilogy(0:num_iter-1, V, 'LineWidth', 4)
% xlabel('iteration(k)', 'FontSize', 24)
% ylabel('V(x_k, \lambda_k)', 'FontSize', 24)
% %ylabel('log V(x_k, \lambda_k)', 'FontSize', 24)
% set(gca,'FontSize',20)

z = cons;

save('sys_params','A', 'D', 'Ec', 'Fc', 'c', 'z')