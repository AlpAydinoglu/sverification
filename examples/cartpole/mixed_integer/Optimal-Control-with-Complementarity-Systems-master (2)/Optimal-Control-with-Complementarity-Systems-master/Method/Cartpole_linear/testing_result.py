# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 15:00:24 2020

@author: alp1a
"""

import numpy as np
import scipy.linalg as linalg

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
A_d = np.identity(A.shape[0]) + Ts*A

B_d = Ts*B
