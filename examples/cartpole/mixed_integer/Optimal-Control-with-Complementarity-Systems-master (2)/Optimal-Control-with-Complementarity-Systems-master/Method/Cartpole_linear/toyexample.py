from system import LinearSystem
from polt_helper import *
from model import *

plt.close('all')
A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m = toy_example()

toy = LinearSystem(A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m)

iter = 6

omg_0_all = np.zeros(((n+m+k)*N, m*N))
delta_0_all = np.zeros(((n+m+k)*N, m*N))

M1 = np.repeat(1000, m*N)
M2 = np.repeat(1000, m*N)


all_states, all_inputs, all_cost_solver, all_cost_real, all_cost_miqp = toy.simulate_system(omg_0_all, delta_0_all, iter, "ADMM_sep",  M1, M2)
# all_states, all_inputs, all_cost_solver, all_cost_real, all_cost_miqp= toy.simulate_system(omg_0_all, delta_0_all, iter, "MIQP", M1, M2)


print(all_cost_miqp)
print(all_cost_solver)
print(all_cost_real)

plot_cost(iter, all_cost_miqp, all_cost_solver, all_cost_real, "ADMM", N)
plot_all_states(all_states, n, N, "ADMM")

# # plot_all_inputs(all_states, k, "MIQP")
# plot_inputs(iter, all_inputs, N)