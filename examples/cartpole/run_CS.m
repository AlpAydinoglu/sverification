%%%Plots envelopes for 1000 trajectories with corresponding Lyapunov functions (Figure 7). (requires sys_params.mat and lyap_params.mat)

clear all
clc
close all

load('sys_params.mat')
load('lyap_params.mat')

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts
 
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