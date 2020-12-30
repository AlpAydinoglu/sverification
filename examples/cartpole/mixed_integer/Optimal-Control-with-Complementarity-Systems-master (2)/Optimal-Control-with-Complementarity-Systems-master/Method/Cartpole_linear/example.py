
from system import LinearSystem
from polt_helper import *
from model import *
from scipy import linalg

data = 5000

inputs = np.zeros(data)
states = np.zeros((data,4))

np.random.seed(1)

for z in range(0, data):

    A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m = Discrete_CartPole()
    
    cartpole = LinearSystem(A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m)
    
    iter = 1
    #delta_0 = np.random.randn((n+m+k)*N, 1)
    # omg_0 = np.random.randn((n+m+k)*N, 1)
    delta_0 = np.zeros(((n+m+k)*N, 1))
    omg_0 = np.zeros(((n+m+k)*N, 1))
    
    omg_0_all = np.zeros(((n+m+k)*N, m*N))
    delta_0_all = np.zeros(((n+m+k)*N, m*N))
    
    M1 = np.repeat(1000, m*N)
    M2 = np.repeat(1000, m*N)
    
    
    # all_states, all_inputs, all_cost_solver, all_cost_real, all_cost_miqp = cartpole.simulate_system(omg_0_all, delta_0_all, iter, "ADMM_sep",  M1, M2)
    # plot_all_states(all_states, n, N, "ADMM")
    
    
    # all_states, all_inputs, all_cost_solver, all_cost_real= cartpole.simulate_system(omg_0, delta_0, iter, "ADMM",  M1, M2)
    
    all_states, all_inputs, all_cost_solver, all_cost_real, all_cost_miqp = cartpole.simulate_system(omg_0, delta_0, iter, "MIQP", M1, M2)
    
    #plot_cost(iter, all_cost_miqp, all_cost_solver, all_cost_real, "MIQP", N)
    #plot_state(iter+1, all_states, N)
    #plot_inputs_plots(iter, all_inputs, N)
    
    states[z,:] = all_states[:,0]
    inputs[z] = all_inputs[0]
    
    print(z)


# n, k, m = 4, 1, 2
# mc, mp = 1, 1
# L, d = 1, 1
# k_1,  k_2= 1, 1
# g = 9.81
# Ts = 0.1    # 1ms

# A = np.array([[0, 0, 1, 0],
#               [0, 0, 0, 1],
#               [0, g * mp / mc, 0, 0],
#               [0, g * (mc + mp) / (L * mc), 0, 0]])
# B = np.array([[0, 0, 1 / mc, 1 / (L * mc)]])
# B = B.reshape((n, 1))
# D = np.array([[0, 0],
#               [0, 0],
#               [0, 0],
#               [-1 / (L * mp), 1 / (L * mp)]])
#     # D = np.zeros((n,m))
# A_d = np.identity(A.shape[0]) + Ts*A

# B_d = Ts*B


# calc = A_d @ all_states[:,0] + B_d @ all_inputs[0]