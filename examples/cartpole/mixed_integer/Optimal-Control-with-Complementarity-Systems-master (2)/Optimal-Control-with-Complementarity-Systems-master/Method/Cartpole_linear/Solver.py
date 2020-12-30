import numpy as np
import cvxpy as cp
import scipy.linalg as linalg
from polt_helper import *

class Solver(object):

    def __init__(self, x_0, n, k, m, N, sys_iter):
        self.x_0 = x_0
        self.n = n
        self.k = k
        self.m = m
        self.N = N
        self.sys_iter = sys_iter

    def ADMM_solver_sep(self, M_1, M_2, M_3, m_1, m_2, m_3, H, omg_0, delta_0, rho_0, MAX_ITER=50):

        """ QCQP nonconvex problem
        min (z.T * H * z)
        s.t M_1 * z == m_1
            M_2 * z <= m_1
            delta.T * M_3 * delta + delta.T * m_3 = 0
            z = delta

        for each k:

            z: ((n+k+m)N, 1)
            delta: ((n+k+m)N, mN)
            omg: ((n+k+m)N, mN)
            rho: (mN, 1)

        for all ADMM k steps:

            z_all: ((n+k+m)N, k_iter)
                z_all[:,0] is meaningless, z starts from z^1
            delta_all: ((n+k+m)N, mN, k_iter)
            omg_all: ((n+k+m)N, mN, k_iter)
            rho_all: ((mN, k_iter)

         """

        n, k, m = self.n, self.k, self.m
        N = self.N

        m_1 = np.squeeze(np.asarray(m_1))
        m_2 = np.squeeze(np.asarray(m_2))

        mu, tau = 10, 1
        z_all, delta_all, omg_all, rho_all, s_all, r_all = self.admm_init(delta_0, omg_0, rho_0, MAX_ITER)
        M_3_all, m_3_all, lmb_P_all, Q_all = self.sep_constraint(M_3, m_3)

        for k_iter in range(MAX_ITER-1):

            """ z_update
              Variable:     z
              Minimize:     z^T H z + rho ||delta - z||_2^2
              Subject to:   M_1 z == m_1
                            M_2 z <= m_2
            """

            z = cp.Variable(((self.n + self.k + self.m) * self.N))
            z_obj = cp.quad_form(z, H)

            for i in range(m*N):
                z_obj = z_obj + rho_all[i, k_iter] * cp.sum_squares(z - delta_all[:, i, k_iter])

            z_cons = [M_1 @ z == m_1, M_2 @ z <= m_2]
            z_update = cp.Problem(cp.Minimize(z_obj), z_cons)
            try:
                z_update.solve(solver='MOSEK', verbose=False)
            except Exception as e:
                print("Error message: \n", e)
                plot_residual_vs_iter(s_all, r_all, rho_all, N, k_iter, self.sys_iter)
                break

            z_all[:, k_iter + 1] = z.value

            # print("status:", z_update.status)
            # print("H, delta and omg", H, delta, omg)
            # print("The optimal value is", z_update.value)
            # print("A solution z is")
            # print("z value in ADMM", z.value)

            for i in range(m*N):
                """ delta_update
                              Variable:     delta
                              Minimize:     ||delta - z - omg||_2^2
                              Subject to:   delta^T M_3 delta + delta^T m_3 = 0
                """

                delta_all[:, i, k_iter + 1] = self.onecons_qcqp(z_all[:, k_iter + 1] + omg_all[:, i, k_iter], M_3_all[:, :, i], m_3_all[:, i], 0, lmb_P_all[:, i], Q_all[:, :, i])

                """ omega_update
                    omg_i^{k+1} = omg_i^k + rho_i^k (z^{k+1} - delta_i^{k+1})
                """
                omg_all[:, i, k_iter + 1] = omg_all[:, i, k_iter] + rho_all[i, k_iter] * (z_all[:, k_iter + 1] - delta_all[:, i, k_iter + 1])

                """ penalty parameter (rho) adjusting
            
                    rho_i^{k+1} =
                        tau rho_i^k         if ||delta_i^{k+1} - z^{k+1}|| > mu ||delta_i^{k+1} - delta_i^k||
                        rho_i^k / tau       if ||delta_i^{k+1} - z^{k+1}|| < mu ||delta_i^{k+1} - delta_i^k||
                        rho^k               otherwise
                """

                """                
                dual residual:
                    s_i^{k+1} = rho_i^{k} ||(delta_i^{k+1} - delta_i^{k})||
                primal residual: 
                    r_i^{k+1} = ||(z^{k+1} - delta_i^{k+1})||
                """
                # abs(rho_all[i, k_iter])*
                s_all[i, k_iter+1] = abs(rho_all[i, k_iter]) * linalg.norm((delta_all[:, i, k_iter+1] - delta_all[:, i, k_iter]))
                r_all[i, k_iter+1] = linalg.norm((z_all[:, k_iter+1] - delta_all[:, i, k_iter+1]))

                if r_all[i, k_iter+1] > mu * s_all[i, k_iter+1]:
                    rho_all[i, k_iter+1] = tau * rho_all[i, k_iter]
                elif r_all[i, k_iter+1] < mu * s_all[i, k_iter+1]:
                    rho_all[i, k_iter+1] = rho_all[i, k_iter] / tau
                else:
                    rho_all[i, k_iter+1] = rho_all[i, k_iter]

            """ Stopping Criterion
                s^{k+1} < 10e-2
                r^{k+1} < 10e-2
            """
            if np.min(s_all[:, k_iter+1]) > 1e3 or np.min(r_all[:, k_iter+1]) > 1e3:
                print("ADMM doesn't converge")
                break

            if np.max(s_all[:, k_iter+1]) < 1e-3 and np.max(r_all[:, k_iter+1]) < 1e-3 or k_iter == MAX_ITER-2:
                if k_iter == MAX_ITER-2:
                    print("stopping criterion: reach max iteration")
                else:
                    print("stopping criterion: err under tolerance")

                plot_residual_vs_iter(s_all, r_all, rho_all, N, k_iter+1, self.sys_iter)
                break

            print("z^{} = ".format(k_iter+1), z.value)
        plot_all_var(z_all, n, k, m, N, "z", self.sys_iter, k_iter+1)
        for i in range(m*N):
            plot_all_var(delta_all[:, i, :], n, k, m, N, r"$\delta_{}$".format(i+1), self.sys_iter, k_iter+1)
            plot_all_var(omg_all[:, i, :], n, k, m, N, r"$\omega_{}$".format(i+1), self.sys_iter, k_iter+1)

        x_seq, u_seq, lmb_seq = self.var_decouple(z_all[:,k_iter+1])

        return x_seq, u_seq, lmb_seq

    def ADMM_solver(self, M_1, M_2, M_3, m_1, m_2, m_3, H, omg_0, delta_0, rho=0.1, MAX_ITER=10):

        """ QCQP nonconvex problem
        min (z.T * H * z)
        s.t M_1 * z == m_1
            M_2 * z <= m_1
            delta.T * M_3 * delta + delta.T * m_3 = 0
            z = delta
         """

        omg = omg_0
        delta = delta_0

        Psymm = (M_3 + M_3.T) / 2.0
        lmb_P, Q = linalg.eigh(np.asarray(Psymm))
        lmb_P = lmb_P.reshape(lmb_P.shape[0], 1)

        for i in range(MAX_ITER):
            """ z_update
              Variable:     z
              Minimize:     z^T H z + rho ||z - delta + omg||_2^2
              Subject to:   M_1 z == m_1
                            M_2 z <= m_2
            """
            z = cp.Variable(((self.n + self.k + self.m) * self.N, 1))

            z_obj = (cp.quad_form(z, H) + rho * cp.sum_squares(z - delta + omg))
            z_cons = [M_1 @ z == m_1, M_2 @ z <= m_2]
            z_update = cp.Problem(cp.Minimize(z_obj), z_cons)
            z_update.solve(solver='MOSEK', verbose=False)
            z = z.value

            # print("status:", z_update.status)
            # print("H, delta and omg", H, delta, omg)
            # print("The optimal value is", z_update.value)
            # print("A solution z is")
            # print("z value in ADMM", z.value)

            """ delta_update
              Variable:     delta
              Minimize:     ||delta - z - omg||_2^2
              Subject to:   delta^T M_3 delta + delta^T m_3 = 0
            """
            delta = self.onecons_qcqp(z + omg, M_3, m_3, 0, lmb_P, Q)

            """ omega_update
                omg = omg + z - delta
            """

            omg = omg + (z - delta)

        x_seq, u_seq, lmb_seq = self.var_decouple(z)

        return x_seq, u_seq, lmb_seq, omg, delta

    def onecons_qcqp(self, z, P, q, r, lmb, Q, tol=1e-9):
        """ Solving nonconvex QCQP with one constraint
          Variable:     x
          Minimize:     ||x-z||_2^2
          Subject to:    x^T P x + q^T x + r <= 0
        """

        if z.T @ P @ z + q.T @ z + r <= 0:
            return z

        else:

            z_hat = Q.T @ z
            q_hat = Q.T @ q
            q_hat = np.squeeze(np.asarray(q_hat))

            xhat = lambda nu: -np.divide(nu * q_hat - 2 * z_hat, 2 * (1 + nu * lmb))
            phi = lambda xhat: lmb.T.dot(np.power(xhat, 2)) + q_hat.T.dot(xhat) + r

            s = -np.inf
            e = np.inf
            for l in lmb:
                if l > 0: s = max(s, -1. / l)
                if l < 0: e = min(e, -1. / l)
            if s == -np.inf:
                s = -1.
                while phi(xhat(s)) <= 0: s *= 2.
            if e == np.inf:
                e = 1.
                while phi(xhat(e)) >= 0: e *= 2.
            while e - s > tol:
                m = (s + e) / 2.
                p = phi(xhat(m))
                if p > 0:
                    s = m
                elif p < 0:
                    e = m
                else:
                    s = e = m
                    break
            nu = (s + e) / 2.
            return Q.dot(xhat(nu))

    def MIQP_solver(self, H_wo_x0, G_1, G_2, G_3, H_3, g_1, g_2, g_3):
        """
        :param x_0: the initial state for MPC problem
        :return: optimal input value at i = 0

          Mixed-Integer Quadratic Program
             z = [x, lmd, u]
             s = {0, 1} in R^m
             min (z.T * H * z)
             s.t G_1 z == g_1
                 G_2 z <= g_2
                 G_3 z + H_3 s <= g3
                 x_0 = x(0)
        """
        n, k, m = self.n, self.k, self.m
        N = self.N

        z = cp.Variable((n + k + m) * N)
        s = cp.Variable(m * N, boolean=True)

        constraints = [G_1 @ z == g_1, G_2 @ z <= g_2, G_3 @ z + H_3 @ s <= g_3  ]
        #constraints = [G_1 @ z == g_1, G_2 @ z <= g_2, G_3 @ z + H_3 @ s <= g_3, z[n*(N-1)] == 0, z[n*(N-1)+1] == 0, z[n*(N-1)+2] == 0, z[n*(N-1)+3] == 0  ]
        # z[n*(N-1)] == 0, z[n*(N-1)+1] == 0, z[n*(N-1)+2] == 0, z[n*(N-1)+3] == 0]        # set x_N = 0
        obj = cp.Minimize(cp.quad_form(z, H_wo_x0))
        prob = cp.Problem(obj, constraints)
        prob.solve(solver='MOSEK')
        print("status:", prob.status)
        # print("optimal value", prob.value)
        # print("optimal var", z.value, s.value)

        x_seq, u_seq, lmb_seq = self.var_decouple(z.value)

        return x_seq, u_seq, lmb_seq

    def var_decouple(self, z):

        n, k, m = self.n, self.k, self.m
        N = self.N

        z = z.reshape((n + k + m) * N, 1)

        x = np.concatenate((z[0: n * N]))
        x = x.reshape((n*N, 1))

        u = z[n * N: (n + k) * N]
        lmb = z[(n + k) * N: (n + k + m) * N]
        return x, u, lmb

    def sep_constraint(self, M_3, m_3):
        n, k, m = self.n, self.k, self.m
        N = self.N

        M_3_all = np.repeat(M_3[:, :, np.newaxis], m*N, axis=2)
        m_3_all = np.repeat(m_3[:, np.newaxis], m*N, axis=1)
        Psymm_all = np.zeros((M_3.shape[0], M_3.shape[1], m*N))
        lmb_P_all = np.zeros((M_3.shape[0], m*N))
        Q_all = np.zeros((M_3.shape[0], M_3.shape[1], m*N))
        for i in range(m*N):
            M_3_all[:, :, i][:(n + k) * N + i] = 0
            M_3_all[:, :, i][(n + k) * N + (i + 1):] = 0
            m_3_all[:, i][:(n + k) * N + i] = 0
            m_3_all[:, i][(n + k) * N + (i + 1):] = 0
            Psymm_all[:, :, i] = (M_3_all[:, :, i] + M_3_all[:, :, i].T) / 2.0
            lmb_P, Q_all[:, :, i] = linalg.eigh(np.asarray(Psymm_all[:, :, i]))
            # lmb_P_all[:, i] = lmb_P.reshape(lmb_P.shape[0], 1)
            lmb_P_all[:, i] = lmb_P

        return M_3_all, m_3_all, lmb_P_all, Q_all

    def admm_init(self, delta_0, omg_0, rho_0, MAX_ITER):
        '''

        :param omg_0:
        :param delta_0:
        :param rho_0:
        :param MAX_ITER:
        :return:


        z starts from z^1
        z^k = z_all[:, k]
        delta_i^k = delta_all[:, i, k]
        omg_i^k = omg_all[:, i, k]
        rho_i^k = rho_all[i, k]
        '''

        n, k, m = self.n, self.k, self.m
        N = self.N

        z_all = np.zeros(((n+k+m)*N, MAX_ITER))
        delta_all = np.zeros(((n+k+m)*N, m*N, MAX_ITER))
        omg_all = np.zeros(((n+k+m)*N, m*N, MAX_ITER))
        rho_all = np.zeros((m*N, MAX_ITER))
        s_all = np.zeros((m*N, MAX_ITER))
        r_all = np.zeros((m*N, MAX_ITER))

        delta_all[:,:,0] = delta_0
        omg_all[:,:,0] = omg_0
        rho_all[:,0] = np.squeeze(np.asarray(rho_0))

        return z_all, delta_all, omg_all, rho_all, s_all, r_all
