import numpy as np
from lemkelcp import lemkelcp
import cvxpy as cp
import scipy.linalg as linalg
from Solver import Solver
import time
from polt_helper import *


class LinearSystem(object):
    """
    x_0: the initial state for the system
    Q: the state cost for i = 0, ... N-1
    R: the input cost for all the u
    P: the terminal cost for x_N
    N: the MPC horizon

    min sum_{t=0} ^N    x_t ^T Q_t x_t + u_t ^T R_t u_t
    s.t. x_{t+1} = A x_t + B u_t + D lambda_t + w, t=0, ..., N-1
         E x_t + F lambda_t + c >= 0
         lambda_t >= 0
         lambda_t ^T  (E x_t + F lambda_t) =0
         x_{0}=x(0)

    """

    def __init__(self, A, B, D, w, E, F, c, x_0, N, Q, P, R, n, k, m,):

        self.A, self.B, self.D, self.w = A, B, D, w
        self.E, self.F, self.c = E, F, c
        self.x_0 = x_0
        self.N = N
        self.n, self.k, self.m = n, k, m
        self.Q, self.P, self.R = Q, P, R

    def simulate_system(self, omg_0, delta_0, iter, solver, M1=[], M2=[]):
        """

        :param x_0: the initial state for MPC problem
        :param delta_0: initial variable value for ADMM algorithm
        :param omg_0: initial variable value for ADMM algorithm
        :param iter: evaluation time for system
        :param solver: "ADMM" or "MIQP"
        :return: history of states and inputs
        """
        x_real = self.x_0
        omg = omg_0
        delta = delta_0
        n, k, m = self.n, self.k, self.m
        N = self.N

        all_states = np.zeros((n, iter+1))
        all_states[:, 0] = np.squeeze(np.asarray(x_real))

        all_inputs = np.zeros((k+m, iter))
        all_cost_miqp = np.zeros((iter, 1))
        all_cost_real = np.zeros((iter, 1))
        all_cost_solver = np.zeros((iter, 1))

        for i in range(iter):
            sol = Solver(x_real, n, k, m, N, i)

            if solver == "ADMM":
                H_wo_x0, M_1, M_2, M_3, m_1, m_2, m_3 = self.helper(x_real, "ADMM")
                x_seq, u_seq, lmb_seq, omg, delta = sol.ADMM_solver(M_1, M_2, M_3, m_1, m_2, m_3, H_wo_x0, omg, delta)
            elif solver == "ADMM_sep":
                H_wo_x0, G_1, G_2, G_3, H_3, g_1, g_2, g_3 = self.helper(x_real, "MIQP", M1, M2)
                x_seq, u_seq, lmb_seq = sol.MIQP_solver(H_wo_x0, G_1, G_2, G_3, H_3, g_1, g_2, g_3)
                all_cost_miqp[i] = self.solver_cost(x_seq, u_seq, lmb_seq, H_wo_x0)

                H_wo_x0, M_1, M_2, M_3, m_1, m_2, m_3 = self.helper(x_real, "ADMM")
                print("x_0", x_real)
                MAX_ITER = 50

                penalty_para = range(1, 2)
                length = len(penalty_para)

                cost_admm_test = np.zeros((length))
                cost_real_test = np.zeros((length))

                for para in range(length):
                    rho_0 = np.ones((m*N, 1)) * penalty_para[para]
                    # print("para", para)

                    x_seq, u_seq, lmb_seq = sol.ADMM_solver_sep(M_1, M_2, M_3, m_1, m_2, m_3, H_wo_x0, omg_0, delta_0, rho_0, MAX_ITER)
                    cost_admm_test[para] = self.solver_cost(x_seq, u_seq, lmb_seq, H_wo_x0)
                    cost_real_test[para] = self.real_cost(x_real, u_seq, lmb_seq, H_wo_x0)
                    # print("true cost: ", cost_real_test[para])
                    # print("ADMM cost: ", cost_admm_test[para])

                # constraint_vs_steps(MAX_ITER, cond3, N, penalty_para)
                # cost_vs_steps(MAX_ITER, cost_admm_test / N, N, penalty_para)

                cost_miqp_test = np.ones((length,1)) * all_cost_miqp[i]

                # plot_cost_vs_penalty(cost_miqp_test, cost_admm_test, cost_real_test, penalty_para, N)


            elif solver == "MIQP":
                H_wo_x0, G_1, G_2, G_3, H_3, g_1, g_2, g_3 = self.helper(x_real, "MIQP", M1, M2)
                x_seq, u_seq, lmb_seq = sol.MIQP_solver(H_wo_x0, G_1, G_2, G_3, H_3, g_1, g_2, g_3)

            all_cost_solver[i] = self.solver_cost(x_seq, u_seq, lmb_seq, H_wo_x0)
            all_cost_real[i] = self.real_cost(x_real, u_seq, lmb_seq, H_wo_x0)

            x_real, lam_real = self.next_state(x_real, u_seq[0:k])

            all_states[:, i+1] = np.squeeze(np.asarray(x_real))
            all_inputs[0:k, i] = u_seq[0:k]
            all_inputs[k:k+m, i] = np.squeeze(np.asarray(lam_real))

        return all_states, all_inputs, all_cost_solver, all_cost_real, all_cost_miqp

    def next_state(self, x, u):

        """
        :param x: current state
        :param u: optimal input for current state
        :return: next state after applying the controller
        """

        M = self.F
        q = self.E.dot(x)+self.c
        sol = lemkelcp(M, q)

        lmb = sol[0]
        dim = lmb.shape[0]
        lmb = lmb.reshape((dim, 1))

        # print('COND1: ', lam >= 0)
        # print('COND2: ', q+ M.dot(lam) >= 0)
        # print('COND3: ', np.multiply((q+ M.dot(lam)),lam)== 0 )

        # Here u is an element, using * to indicate the multiplication
        x_next = self.A @ x + self.B * u + self.D @ lmb + self.w

        return x_next, lmb

    def real_cost(self, x_0, u_seq, lmb_seq, H):
        n, N = self.n, self.N
        x_seq = np.zeros((n * (N+1), 1))
        x_seq[0 : n] = x_0
        for i in range(self.N):
            x_seq[(i+1)*n: (i+2)*n], lmb = self.next_state(x_seq[i*n: (i+1)*n], u_seq[i])

        z = np.concatenate([x_seq[n : n * (N+1)], u_seq, lmb_seq])
        cost = (z.T @ H @ z) / N

        return cost

    def solver_cost(self, x_seq, u_seq, lmb_seq, H):

        z = np.concatenate([x_seq, u_seq, lmb_seq])
        cost = (z.T @ H @ z) / self.N
        return cost


    def helper(self, x_0, solver, M1 = [], M2 = []):
        """
        :purpose: matrix concatenation
        :param x_0: the initial state for MPC problem
        :param Q: the state cost for i = 0, ... N-1
        :param R: the input cost for all the u
        :param P: the terminal cost for x_N
        :param N: the MPC horizon
        :param n: dim of states
        :param k: dim of inputs
        :param m: dim of compact forces
        :param solver: "ADMM" or "MIQP"
        :return: the compact matrix for optimization problem
        """
        n, k, m = self.n, self.k, self.m
        Q, P, R = self.Q, self.P, self.R
        N = self.N

        Q_all = np.kron(np.identity(N-1), Q)
        R_all = np.kron(np.identity(N), R)
        lmd_all = np.zeros((m*N, m*N))

        # cost without x_0 term
        H_wo_x0 = linalg.block_diag(Q_all, P, R_all, lmd_all)

        tmp1 = np.kron(np.identity(N - 1), -self.A)
        tmp1 = np.vstack([np.zeros([n, n * (N - 1)]), tmp1])
        tmp1 = np.hstack((tmp1, np.zeros([n * N, n])))
        M_1_part1 = tmp1 + np.identity(N * n)
        M_1_part2 = np.kron(np.identity(N), -self.B)
        M_1_part3 = np.kron(np.identity(N), -self.D)
        M_1 = np.hstack((M_1_part1, M_1_part2, M_1_part3))

        m_1 = np.kron(np.ones((N,1)), self.w)
        m_1[0:n] += self.A @ x_0

        tmp2 = np.kron(np.identity(N - 1), -self.E)

        tmp2 = np.vstack([np.zeros([m, n * (N - 1)]), tmp2])
        tmp2 = np.hstack((tmp2, np.zeros([m * (N), n])))
        M_2_eq1_part1 = tmp2
        M_2_eq1_part2 = np.kron(np.identity(N), np.zeros((m, k)))
        M_2_eq1_part3 = np.kron(np.identity(N), -self.F)
        M_2_eq1 = np.hstack((M_2_eq1_part1, M_2_eq1_part2, M_2_eq1_part3))

        M_2_eq2_part1 = np.zeros((m*N, (n+k)*N))
        M_2_eq2_part2 = -np.identity(m*N)
        M_2_eq2 = np.hstack((M_2_eq2_part1, M_2_eq2_part2))
        M_2 = np.vstack((M_2_eq1, M_2_eq2))

        m_2_eq1 = np.kron(np.ones((N,1)), self.c)
        m_2_eq1[0:m] += self.E @ x_0
        m_2_eq2 = np.zeros((m*N,1))
        m_2 = np.vstack((m_2_eq1, m_2_eq2))

        if solver == "ADMM":

            M_3_part1 = np.zeros(((n+k)*N,(n+k+m)*N))
            M_3_part2 = -M_2_eq1
            M_3 = np.vstack((M_3_part1, M_3_part2))

            m_3_eq1 = np.zeros(((n+k)*N, 1))
            m_3_eq2 = m_2_eq1
            m_3 = np.vstack((m_3_eq1, m_3_eq2))

            return H_wo_x0, M_1, M_2, M_3, m_1, m_2, m_3

        elif solver == "MIQP":
            G_1 = M_1
            g_1 = m_1
            g_1 = np.squeeze(np.asarray(g_1))
            G_2 = M_2
            g_2 = m_2
            g_2 = np.squeeze(np.asarray(g_2))
            G_3 = -G_2

            H_3_part1 = -np.diag(M1)
            H_3_part2 = np.diag(M2)
            H_3 = np.vstack((H_3_part1, H_3_part2))

            g_3_part1 = -m_2_eq1
            g_3_part2 = np.reshape(M2, (m*N, 1))
            g_3 = np.vstack((g_3_part1, g_3_part2))
            g_3 = np.squeeze(np.asarray(g_3))

            return H_wo_x0, G_1, G_2, G_3, H_3, g_1, g_2, g_3