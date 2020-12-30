clear all
clc
close all

%addpath
%addpath(genpath('C:\Users\alp1a\OneDrive\Masaüstü\research\YALMIP-master\YALMIP-master'))
%addpath 'C:\Program Files\Mosek\9.2\toolbox\R2015a'
%addpath 'C:\Program Files\Mosek\9.2\toolbox\R2015a'

addpath 'C:\Users\alp1a\pathlcp\pathmexw64'
addpath(genpath('C:\Users\alp1a\OneDrive\Masaüstü\research\YALMIP-master\YALMIP-master'))
addpath 'C:\Program Files\Mosek\9.2\toolbox\R2015a'

eps = 10^-3;
load('sys_params.mat')

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

%define variables
%the state variables
for i = 1:n
   x{i} = sdpvar(1,1); 
end
%variables related to the contact force
for i = 1:m
    lam{i} = sdpvar(1,1); %lambda_{k}
    lamd{i} = sdpvar(1,1); %\lambda_{k+1}
end

%basis vectors
%the state vector
x_basis = [];
for i = 1:n
    x_basis = [x_basis; x{i}];
end
%the contact vector and related variables
lam_basis = [];
lamd_basis = [];
for i = 1:m
    lam_basis = [lam_basis; lam{i}];
    lamd_basis = [lamd_basis; lamd{i}];
end

%Lyap Function
%Define the variables
P = sdpvar(n,n); Q = sdpvar(n,m); R = sdpvar(m,m); c1 = sdpvar(1,n); c2 = sdpvar(1,m); c3 = sdpvar(1,1);
%Construct the Lyapunov function
V = x_basis' * P * x_basis + 2 * x_basis' * Q  * lam_basis ...
    + lam_basis'  * R * lam_basis + c1 * x_basis + c2 * lam_basis + c3;
%Define the derivative vectors
x_plus = A*x_basis + D*lam_basis + cons;
%Define the derivative of the Lyapunov function
Vdot = x_plus' * P * x_plus - x_basis' * P * x_basis + 2 * x_plus' * Q * lamd_basis - 2 * x_basis' * Q * lam_basis ...
    + lamd_basis' * R * lamd_basis - lam_basis' * R * lam_basis + c1 * x_plus - c1 * x_basis + c2 * lamd_basis - c2*lam_basis;

%initalize the constraint set
F = [];

%S-procedure terms
%(Ex + F \lam + c ) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ1{i} = Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i); 
    INEQ1_M{i} = sdpvar(1,1); F = [F, INEQ1_M{i} >= 0];
    INEQ1_MD{i} = sdpvar(1,1); F = [F, INEQ1_MD{i} >= 0];
    INEQ1_MDD{i} = sdpvar(1,1); F = [F, INEQ1_MDD{i} >= 0];
end
%(\lam) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ2{i} = lam_basis(i); INEQ2_M{i} = sdpvar(1,1); F = [F, INEQ2_M{i} >= 0];
    INEQ2_MD{i} = sdpvar(1,1); F = [F, INEQ2_MD{i} >= 0];    
    INEQ2_MDD{i} = sdpvar(1,1); F = [F, INEQ2_MDD{i} >= 0];
end
%(\lam)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ4{ind} = lam_basis(i)*lam_basis(j); INEQ4_M{ind} = sdpvar(1,1); F = [F, INEQ4_M{ind} >= 0];
        INEQ4_M2{ind} = sdpvar(1,1); F = [F, INEQ4_M2{ind} >= 0];
        INEQ4_MDD{ind} = sdpvar(1,1); F = [F, INEQ4_MDD{ind} >= 0];
    end
end 

%(Ex + F \lam + c)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ6{ind} = (Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i))*(Ec(j,:)*x_basis + Fc(j,:)*lam_basis + c(j)); 
        INEQ6_M{ind} = sdpvar(1,1); F = [F, INEQ6_M{ind} >= 0]; 
        INEQ6_M2{ind} = sdpvar(1,1); F = [F, INEQ6_M2{ind} >= 0];
        INEQ6_MDD{ind} = sdpvar(1,1); F = [F, INEQ6_MDD{ind} >= 0];
    end
end

%\lam^T (Ex + F\lam + c) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ1{i} = lam_basis(i) * (Ec(i,:)*x_basis + Fc(i,:)*lam_basis + c(i)); EQ1_M{i} = sdpvar(1,1); EQ1_MD{i} = sdpvar(1,1); EQ1_MDD{i} = sdpvar(1,1);
end

%(Ex_plus + F \lamd + c ) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ7{i} = Ec(i,:)*x_plus + Fc(i,:)*lamd_basis + c(i); 
    %INEQ7_M{i} = sdpvar(1,1); F = [F, INEQ7_M{i} >= 0];
    INEQ7_MD{i} = sdpvar(1,1); F = [F, INEQ7_MD{i} >= 0];
end

%(Ex_plus + F \lamd + c)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ8{ind} = (Ec(i,:)*x_plus + Fc(i,:)*lamd_basis + c(i))*(Ec(j,:)*x_plus + Fc(j,:)*lamd_basis + c(j)); 
        %INEQ8_M{ind} = sdpvar(1,1); F = [F, INEQ8_M{ind} >= 0]; 
        INEQ8_M2{ind} = sdpvar(1,1); F = [F, INEQ8_M2{ind} >= 0];
    end
end

%(\lamd) >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    INEQ9{i} = lamd_basis(i); %INEQ9_M{i} = sdpvar(1,1); F = [F, INEQ9_M{i} >= 0];
    INEQ9_MD{i} = sdpvar(1,1); F = [F, INEQ9_MD{i} >= 0];    
end

%(\lamd)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    for j = 1:m
        ind = j+(i-1)*m;
        INEQ10{ind} = lamd_basis(i)*lamd_basis(j); %INEQ10_M{ind} = sdpvar(1,1); F = [F, INEQ10_M{ind} >= 0];
        INEQ10_M2{ind} = sdpvar(1,1); F = [F, INEQ10_M2{ind} >= 0];
    end
end

%\lamd^T (Ex_plus + F\lamd + c) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    EQ2{i} = lamd_basis(i) * (Ec(i,:)*x_plus + Fc(i,:)*lamd_basis + c(i)); %EQ2_M{i} = sdpvar(1,1); 
    EQ2_MD{i} = sdpvar(1,1);
end

alpha = sdpvar(1,1); F = [F, alpha >= eps];

%inequalities
ineq1 = V - eps * (x_basis' * x_basis); 
ineq2 = - Vdot - eps * (x_basis' * x_basis);
ineq3 = -V + alpha * (x_basis' * x_basis);

%add S-procedure terms
for i = 1:m
    ineq1 = ineq1 - INEQ1{i}*INEQ1_M{i} - INEQ2{i}*INEQ2_M{i}  - EQ1{i}*EQ1_M{i};
    ineq2 = ineq2 - INEQ1{i}*INEQ1_MD{i} - INEQ2{i}*INEQ2_MD{i} ... 
        - EQ1{i}*EQ1_MD{i} - EQ2{i}*EQ2_MD{i} - INEQ7{i}*INEQ7_MD{i} - INEQ9{i}*INEQ9_MD{i};
    ineq3 = ineq3 - INEQ1{i}*INEQ1_MDD{i} - INEQ2{i}*INEQ2_MDD{i} - EQ1{i}*EQ1_MDD{i}; 
end

for i = 1:(m*m)
    ineq1 = ineq1 - INEQ4{i}*INEQ4_M{i} - INEQ6{i}*INEQ6_M{i};
    ineq2 = ineq2 - INEQ4{i}*INEQ4_M2{i} - INEQ6{i}*INEQ6_M2{i} - INEQ8{i}*INEQ8_M2{i} - INEQ10{i}*INEQ10_M2{i};
    ineq3 = ineq3 - INEQ4{i}*INEQ4_MDD{i} - INEQ6{i}*INEQ6_MDD{i};
    if mod(i,100) == 0
        i 
    end
end

%Construct the sos program
v = monolist([x_basis' lam_basis'],1);
Ke = sdpvar(length(v));
p_sos = v'*Ke*v;
F = [F, coefficients(ineq1-p_sos,[x_basis' lam_basis']) == 0, Ke>=0];

h = monolist([x_basis' lam_basis' lamd_basis'],1);
He = sdpvar(length(h));
q_sos = h'*He*h;

F = [F, coefficients(ineq2-q_sos,[x_basis' lam_basis' lamd_basis']) == 0, He>=0];

t = monolist([x_basis' lam_basis'],1);
Fe = sdpvar(length(t));
f_sos = t'*Fe*t;
F = [F, coefficients(ineq3-f_sos,[x_basis' lam_basis']) == 0, Fe>=0];

options = sdpsettings('solver','mosek','dualize',1);
optimize(F, [] , options)
PP = double(P);
QQ = double(Q);
RR = double(R);
cc1 = double(c1);
cc2 = double(c2);
cc3 = double(c3);

%save the matrices
save('system_param.mat','PP','QQ','RR', 'cc1', 'cc2', 'cc3')