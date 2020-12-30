import numpy as np
import control
from scipy.integrate import solve_ivp






class LinearSystem(object):
    def __init__(self, A, B):
        self.A = A
        self.B = B

    def cartpole_input(self, t, x):
        g = 9.8100
        Etilde = .5 * x[3] ** 2 - g * np.cos(x[1]) - g
        q1_ddot_des = x[3] * np.cos(x[1]) * Etilde - 5 * x[2] - x[0]
        u = (2 - np.cos(x[1]) ** 2) * q1_ddot_des - g * np.sin(x[1]) * np.cos(x[1]) - np.sin(x[1]) * x[3] ** 2
        return u

    def cart_pole(self, t, x):
        theta = x[1] % 2*np.pi
        g = 9.8100
        mc = 1
        mp = 1
        L = 1
        M = np.array([[mc + mp, mp*L*np.cos(theta)], [mp*L*np.cos(theta), mp*L**2]])
        C = np.array([[-mp*L*np.sin(theta)*x[3]**2], [mp*g*L*np.sin(theta)]])
        BB = np.array([[1], [0]])
        Etilde = .5 * x[3] ** 2 - g * np.cos(theta) - g
        angular_distance = x[3] ** 2 + (theta - np.pi) ** 2
        if np.abs(Etilde) < 1 and angular_distance < 1:
            # if we are close enough to desired state, let LQR take over
            K = control.lqr(self.A, self.B, np.diag([1, 1, 1, 10]), np.array([1]))
            tmp = np.array([0,np.pi,0,0])
            tmp.reshape((4,1))
            u = - K*(x-tmp)
        else:
            # else, use energy shaping
            u = self.cartpole_input(t, x)

        tmp = (np.linalg.inv(M).dot(BB*u - C))
        tmp.tolist()

        dx = [x[2],x[3], tmp[0][0], tmp[1][0]]
        return dx

    def simulate_cartpole(self):
        tf = 10
        # x_0 = np.array([[0], np.pi / 4, 0, 0])
        # x_0.reshape([4, 1])
        x_0 = [0, np.pi / 4, 0, 0]
        sol = solve_ivp(self.cart_pole, [0, tf], x_0, method='RK45', max_step=0.1, rtol=1e-4, atol=1e-4)

        return sol

