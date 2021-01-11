###Creates 5000 samples for the optimal control problem
###Saves the samples as states_ARE.npy and inputs_ARE.npy

from system import LinearSystem
from polt_helper import *
from model import *
from scipy import linalg

#number of samples
data = 5000

#save the samples in inputs and states
inputs = np.zeros(data)
states = np.zeros((data,4))

#seeding for repeatability
np.random.seed(1)

#sampling
for z in range(0, data):
    
    A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m = Discrete_CartPole()
    
    cartpole = LinearSystem(A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m)
    
    iter = 1
    delta_0 = np.zeros(((n+m+k)*N, 1))
    omg_0 = np.zeros(((n+m+k)*N, 1))
    
    omg_0_all = np.zeros(((n+m+k)*N, m*N))
    delta_0_all = np.zeros(((n+m+k)*N, m*N))
    
    M1 = np.repeat(1000, m*N)
    M2 = np.repeat(1000, m*N)
    
    all_states, all_inputs, all_cost_solver, all_cost_real, all_cost_miqp = cartpole.simulate_system(omg_0, delta_0, iter, "MIQP", M1, M2)
    
    states[z,:] = all_states[:,0]
    inputs[z] = all_inputs[0]
    
    print(z)
    
    #save the samples
    np.save('states_ARE', states)
    np.save('inputs_ARE', inputs)