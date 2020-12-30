import numpy as np
import scipy.linalg as linalg
import cvxpy as cp

# tic = time.perf_counter()
# toc = time.perf_counter()
# print(f"each MPC step takes {toc - tic:0.4f} seconds")

# E = np.array([[1, 1]])
# F = np.array([[1]])
# c = np.array([[0]])
#
# lam = z_all[(n + k) * N:, k_iter + 1].reshape((N, m))
# lam = lam.T
# x_fortest = z_all[:n * N, k_iter + 1].reshape((N, n))
# x_fortest = x_fortest.T
#
# cond1 = np.zeros((m, 1))
# cond2 = np.zeros((m, 1))
# for j in range(N):
#     cond1 += (E @ x_fortest[:, j] + F @ lam[:, j] + c).reshape((m, 1))
#     cond2 += lam[:, j].reshape((m, 1))
#
# print("cond1, ", cond1)
# print("cond2, ", cond2)
# print("cond3: product, ", cond3)


A_d = np.array(([1, 1], [0, 1.01]))
C = linalg.eigvals(A_d)
print(C)