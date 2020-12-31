import numpy as np
from scipy import linalg

def toy_example():
    n, k, m = 2, 1, 1
    A_d = np.array(([1, 1], [0, 1.01]))
    B_d = np.array(([0], [1]))
    D_d = np.array(([0], [1]))

    w = np.zeros((n, 1))
    E = np.array([[1,1]])
    F = np.array([[1]])
    c = np.array([[0]])
    N = 2

    x_0 = np.array(([0.01], [0.01]))

    Q = np.identity(n)
    R = np.identity(k)
    # P = linalg.solve_discrete_are(A, B, Q, R)
    P = 10000*np.identity(n)

    return A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m


def Discrete_CartPole():
    """
    n, k, m: num of states, inputs, compact forces
    mc, mp: mass of the cart and the pole
    L, d: length of the pole, s half of the distance between the walls
    k_1,  k_2: stiffness parameters of the soft walls
    g: gravitational acceleration
    Ts:  explicit Euler method time step
    A: (n, n)
    B: (n, k)
    D: (n, m)
    w: (n, 1)
    E: (m, n)
    F: (m, m)
    c: (m, 1)
    x: (n, 1)
    u: (k, 1)
    """
    n, k, m = 4, 1, 2
    mc, mp = 1, 1
    L, d = 1, 1
    k_1,  k_2= 1, 1
    g = 9.81
    Ts = 0.1    # 1ms

    A = np.array([[0, 0, 1, 0],
                  [0, 0, 0, 1],
                  [0, g * mp / mc, 0, 0],
                  [0, g * (mc + mp) / (L * mc), 0, 0]])
    B = np.array([[0, 0, 1 / mc, 1 / (L * mc)]])
    B = B.reshape((n, 1))
    D = np.array([[0, 0],
                  [0, 0],
                  [0, 0],
                  [-1 / (L * mp), 1 / (L * mp)]])
    # D = np.zeros((n,m))
    
    #A_d = np.identity(A.shape[0]) + Ts*A
    
    A_d = np.array([[1, 0, 0.1, 0],
                  [0, 1, 0, 0.1],
                  [0.1846, -4.4253, 1.4139, -1.4272],
                  [0.1846, -3.4443, 0.4139, -0.4272]])
    B_d = Ts*B

    D_d = Ts*D
    w = np.zeros((n, 1))
    E = np.array([[-1, L, 0, 0],
                  [1, -L, 0, 0]])
    F = np.array([[1 / k_1, 0],
                  [0, 1 / k_2]])
    c = np.array([[d, d]])
    c = c.reshape((m, 1))


    x_0 =  3 * np.random.random((n, 1)) - 1.5
    #x_0 = np.array(([0.1, 0, 0, 0.0]))
    x_0 = x_0.reshape((n, 1))
    # x_0[0] = np.random.random()/10
    
    N = 10
    Q = 10 * np.identity(n)
    R = 1*np.identity(k)
    P = linalg.solve_discrete_are(A_d, B_d, Q, R)
    #P = np.array([[  292.3727184 , -1623.31261745,   362.79039641,  -516.39416876],

    #P = 10*np.identity(n)
    return A_d, B_d, D_d, w, E, F, c, x_0, N, Q, P, R, n, k, m