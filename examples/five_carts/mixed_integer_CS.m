function uout = mixed_integer_CS(x0, N, A, B, D, z, Ec, Fc, c, H)

M = 100;

%dimensions
n = size(A,2); m = size(D,2); k = size(B,2);

%check if H is zero matrix
if sum(H == 0) == size(H,1)*size(H,2)
    lam0 = pathlcp(Fc, Ec * x0 + c); 
end

%Define variables
for i = 1 : N+1
    x{i} = sdpvar(n,1);
    if i < N+1
        lam{i} = sdpvar(m,1);
        u{i} = sdpvar(k,1);
        x_int{i} = binvar(m,1);
    end
end

%constraints
F = [];
F = [F, x{1} == x0];
for i = 1:N
   F = [F, x{i+1} == A * x{i} + B * u{i} + D * lam{i} + z];
   F = [F, M * x_int{i} >= Ec * x{i} + Fc * lam{i} + H * u{i} + c >= 0, M * ( ones(m,1) - x_int{i}) >= lam{i} >= 0];
end

%objective
obj = 0;
P = 10*eye(n);
R = 1*eye(k);
Pn = 10*eye(n);

for i = 1:N
   obj = obj + x{i}' * P * x{i} + u{i}' * R * u{i};
end


%final constraint
obj = obj + x{N+1}' * Pn * x{N+1};

% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','mosek');

% Solve the problem
sol = optimize(F,obj,options);

xout = [];
lamout = [];
uout = [];

%output
for i = 1:N
    xx{i} = double(x{i});
    ll{i} = double(lam{i});
    uu{i} = double(u{i});
    xout = [xout xx{i}];
    lamout = [lamout ll{i}];
    uout = [uout uu{i}];
end
xx{N+1} = double(x{N+1});
xout = [xout xx{N+1}];

uout = uout(:,1);

end

