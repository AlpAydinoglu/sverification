%%%Creates 2000 samples for the optimal control policy and saves the samples as data_training.mat
clear
clc
close all

N=10; %optimal control horizon

%system parameters
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

%LQR controller
sys = ss(A2,B2,zeros(n), zeros(n,k),1);
K = lqr(sys, 1*eye(n), 1*eye(k));
A2 = A2 - B2*K;
A=A2;
D = D2;
cons = zeros(n,1);
c = c2;
H = H2;
Ec = Ec2;
Fc = Fc2;
B = B2;

%seed random number generator for repeatability
rng(0)

%sample_points
n_samples = 2000; %number of samples
data_state = zeros(n, n_samples); data_input = zeros(k, n_samples); %save the generated data

%sample from the optimal control policy
for i = 1:n_samples
    i
    data_state(:,i) = 10 * (rand(n,1) - 0.5);
    data_input(:,i) = mixed_integer_CS(data_state(:,i), N, A, B, D, cons, Ec, Fc, c, H);
end

%save the generated data
save('data_training','data_state', 'data_input')