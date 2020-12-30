clear all
clc
close all

eps = 10^-3;
load('sys_params.mat')

%extract dimension information
n = size(A,2); %dimension of state space
m = size(D,2); %number of contacts

%Lyap Function
%Define the variables
P = sdpvar(n,n); Q = sdpvar(n,m); R = sdpvar(m,m); c1 = sdpvar(n,1); c2 = sdpvar(m,1); c3 = sdpvar(1,1);

%initalize the constraint set
F = [];

%define V and DV
V = [P Q c1/2; Q' R c2/2; c1'/2 c2'/2 c3];
dv11 = A' * P * A - P; dv12 = A' * P * D - Q; dv13 = A' * Q; dv14 = A' * P * z - (c1/2); dv21 = D' * P * A - Q'; dv22 = D' * P * D - R; dv23 = D' * Q; dv24 = D' * P * z + (D' * c1 /2) - c2/2;
dv31 = Q' * A;  dv32 = Q' * D; dv33 = R; dv34= Q' * z + c2/2; dv41 = (z' * P * A - c1'/2); dv42 = z' * P * D + (c1' * D / 2) - (c2'/2)  ; dv43 = (z'*Q + c2'/2); dv44 = z' * P * z + c1' * z;
DV = [dv11  dv12  dv13  dv14; dv21 dv22 dv23 dv24 ; dv31  dv32 dv33 dv34; dv41 dv42 dv43 dv44];

%S-procedure terms
%(\lam)^2 >= 0 and (Ex + F \lam + c)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_quad = [Ec Fc c; zeros(m,n) eye(m) zeros(m,1); zeros(1,m+n) 1];
DV_quad = [Ec Fc zeros(m,m) c; zeros(m,n) eye(m) zeros(m,m+1); zeros(1,n+2*m) 1];
W = sdpvar(2*m+1,2*m+1); U = sdpvar(2*m+1,2*m+1); W2 = sdpvar(2*m+1,2*m+1);
F = [F, U(:) >= 0, W(:) >= 0, W2(:) >= 0];

%\lam^T (Ex + F\lam + c) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for ([x lam 1])
for i = 1:m
    Zc_nonsym = [zeros(n+i-1, n+m+1); Ec(i,:) Fc(i,:) c(i,:); zeros(m-i+1, n+m+1)]; %create nonsymmetric version
    Zc{i} = (Zc_nonsym + Zc_nonsym')/2; %symmetrize the matrix
    tau{i} = sdpvar(1,1);
    tau2{i} = sdpvar(1,1);
end

%extension to ([x lam lamd 1])
for i = 1:m
    Zc_nonsym = [zeros(n+i-1, n+2*m+1); Ec(i,:) Fc(i,:) zeros(1,m) c(i,:); zeros(m-i, n+2*m+1); zeros(m+1, n+2*m+1)]; %create nonsymmetric version
    Zc_e{i} = (Zc_nonsym + Zc_nonsym')/2; %symmetrize the matrix
    tau_e{i} = sdpvar(1,1);
end

%(Ex_plus + F \lamd + c)^2 >= 0, (\lamd)^2 >= 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DV_quad2 = [Ec*A Ec*D Fc Ec*z+c; zeros(m,n+m) eye(m) zeros(m,1)];
U_V2 = sdpvar(2*m,2*m); 
F = [F, U_V2(:) >= 0];

%\lamd^T (Ex_plus + F\lamd + c) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:m
    h1 = Ec*A; h2 = Ec*D; h3 = Ec*z + c;
    Zc_nonsym = [zeros(n+m+i-1, n+2*m+1); h1(i,:) h2(i,:) Fc(i,:) h3(i,:); zeros(m-i, n+2*m+1); zeros(1, n+2*m+1)]; %create nonsymmetric version
    Zc_ee{i} = (Zc_nonsym + Zc_nonsym')/2; %symmetrize the matrix
    tau_ee{i} = sdpvar(1,1);
end

%construct the inequalities
ineq1 = V - V_quad' * W * V_quad;
ineq2 = -DV - DV_quad' * U * DV_quad - DV_quad2' * U_V2 * DV_quad2;
ineq3 = -V - V_quad' * W2 * V_quad;

for i = 1:m
    ineq1 = ineq1 - tau{i}*Zc{i};
    ineq2 = ineq2 - tau_e{i}*Zc_e{i} - tau_ee{i}*Zc_ee{i};
    ineq3 = ineq3 - tau2{i}*Zc{i};
end

alpha = sdpvar(1,1); F = [F, alpha >= eps];

matl = zeros(n+m+1); matl(1:n,1:n) = eps*eye(n); matk = zeros(n+2*m+1); matk(1:n,1:n) = eps*eye(n); matz = zeros(n+m+1); matz(1:n,1:n) = eye(n);
F = [F, ineq1 >= matl, ineq2 >= matk, ineq3 >= -alpha*matz];

%solve the problem
options = sdpsettings('solver','mosek','dualize',1);
optimize(F, [] , options)
PP = double(P);
QQ = double(Q);
RR = double(R);
cc1 = double(c1);
cc2 = double(c2);
cc3 = double(c3);

%save the matrices
save('lyap_params.mat','PP','QQ','RR', 'cc1', 'cc2', 'cc3')