%%%Creates 2000 samples for the optimal control policy and saves the samples as data_training.mat
clear
clc
close all

Ts = 0.1; %sampling time
N = 5; %optimal control horizon

%system parameters
A = [1 Ts;0 1];
D = [0 0 0; 1 -1 0]; D = D * Ts;
Ec = [0 1;0 -1;0 0];
Fc = [1 -1 1; -1 1 1; -1 -1 0];
c = [0;0;0.981];
H = [1;-1;0];
cons = [0;0];
B = [0;1]; B = B*Ts;

%LQR controller design
sys = ss(A,B,zeros(2), zeros(2,1),1);
K = lqr(sys, 10*eye(2), 1*eye(1));
A = A - B*K;
Ec = Ec - H*K;

%seed random number generator for repeatability
rng(0)

%sample_points
n_samples = 2000; %number of samples
n=2; %state dimension
k=1; %input dimension
data_state = zeros(n, n_samples); data_input = zeros(k, n_samples); %save the generated data

%sample from the optimal control policy
for i = 1:n_samples
    i
    data_state(:,i) = 10 * (rand(n,1) - 0.5);
    data_input(:,i) = mixed_integer_CS(data_state(:,i), N, A, B, D, cons, Ec, Fc, c, H);
end

%save the generated data
save('data_training','data_state', 'data_input')