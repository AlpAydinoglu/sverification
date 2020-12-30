function uout = mixed_integer_CS(x0, N, A, B, D, z, Ec, Fc, c, H)

%dimensions
n = size(A,2); m = size(D,2); k = size(B,2);

%Define variables
for i = 1 : N+1
    x{i} = sdpvar(n,1);
    if i < N+1
        u{i} = sdpvar(k,1);
    end
end

%constraints
F = [];
F = [F, x{1} == x0];
for i = 1:N
   F = [F, x{i+1} == A * x{i} + B * u{i} + z];
   F = [F, u{i} <= 3, u{i} >= -3, -4*ones(2,1) <= x{i}, x{i} <= 4*ones(2,1)];
end

%objective
obj = 0;
P = 10*eye(2);
R = 1;
Pn = 10*eye(2);

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
    uu{i} = double(u{i});
    xout = [xout xx{i}];
    uout = [uout uu{i}];
end
xx{N+1} = double(x{N+1});
xout = [xout xx{N+1}];

uout = uout(1);

end

