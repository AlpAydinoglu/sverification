# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 22:56:06 2020

@author: alp
"""

# Import packages.
import cvxpy as cp
import numpy as np
import scipy.io

data = scipy.io.loadmat('sys_params.mat')
A = data['A']
D = data['D']
E = data['Ec']
F = data['Fc']
c = data['c']
z = data['z']

n = np.size(A,0)
m = np.size(F,0)
np.random.seed(1)
eps = 0.001

#Define Lyapunov function variables and the constraints set

P = cp.Variable((n,n), symmetric=True)
Q = cp.Variable((n,m), symmetric=False)
R = cp.Variable((m,m), symmetric=True)
c1 = cp.Variable((n,1), symmetric=False)
c2 = cp.Variable((m,1), symmetric=False)
c3 = cp.Variable((1,1), symmetric=False)

constraints = []

#V and DV

V = cp.bmat([ [P, Q, c1/2], [Q.T, R, c2/2 ], [c1.T / 2, c2.T / 2, c3]     ])
DV11 = (A.T @ P @ A) - P 
DV12 = (A.T @ P @ D) - Q
DV13 = A.T @ Q
DV14 = (A.T @ P @ z) - (c1/2)
DV21 = D.T @ P @ A - Q.T
DV22 = D.T @ P @ D - R
DV23 = D.T @ Q
DV24 = (D.T @ P @ z) + (D.T @ c1 /2) - c2/2
DV31 = Q.T @ A
DV32 = Q.T @ D
DV33 = R
DV34= Q.T @ z + c2/2
DV41 = (z.T @ P @ A - c1.T / 2)
DV42 = z.T @ P @ D + (c1.T @ D / 2) - (c2.T / 2)
DV43 = (z.T @ Q + c2.T / 2)
DV44 = z.T @ P @ z + c1.T @ z;
DV = cp.bmat([ [DV11, DV12, DV13, DV14], [DV21, DV22, DV23, DV24], [DV31, DV32, DV33, DV34], [DV41, DV42, DV43, DV44]])

#S-procedure terms

#(\lam)^2 >= 0 and (Ex + F \lam + c)^2 >= 0

Vquad = cp.bmat([ [E, F, c], [np.zeros([m,n]), np.eye(m), np. zeros([m,1])] , [np.zeros([1,n+m]) , np.eye(1)]       ])
DVquad = cp.bmat([ [E, F, np.zeros([m,m]),  c], [np.zeros([m,n]), np.eye(m), np. zeros([m,m+1])] , [np.zeros([1,n+2*m]) , np.eye(1)]]) 
W = cp.Variable((2*m+1,2*m+1), symmetric=False)
U = cp.Variable((2*m+1,2*m+1), symmetric=False)
W2 = cp.Variable((2*m+1,2*m+1), symmetric=False)
constraints += [W >= 0 , U >= 0, W2 >= 0]

#(Ex_plus + F \lamd + c)^2 >= 0, (\lamd)^2 >= 0
DVquad2 = cp.bmat([ [E @ A, E @ D, F,  E@z + c], [np.zeros([m,n+m]), np.eye(m), np.zeros([m,1]) ]    ]) 
U_V2 = cp.Variable((2*m,2*m), symmetric=False)
constraints += [U_V2 >= 0]

#x^2 >= bon
bon = 0.1;
bnd = np.zeros([n+m+1,n+m+1]) 
bnd[0:n, 0:n] = np.eye(n)
bnd[n+m,n+m] = -bon
bnde = np.zeros([n+2*m+1,n+2*m+1]) 
bnde[0:n, 0:n] = np.eye(n)
bnde[n+2*m,n+2*m] = -bon
t1 = cp.Variable((1,1))
t2 = cp.Variable((1,1))
t3 = cp.Variable((1,1))

#\lam^T (Ex + F\lam + c) = 0
tau1 = cp.Variable((m,1))
J1 = cp.diag(tau1)
tau2 = cp.Variable((m,1))
J2 = cp.diag(tau2)
tau3 = cp.Variable((m,1))
J3 = cp.diag(tau3)
S1 = cp.bmat([ [np.zeros([n,n+m+1])], [J1 @ E, J1 @ F, J1 @ c] , [np.zeros([1,n+m+1])]      ])
JJ1 = (S1.T + S1)/2
S2 = cp.bmat([ [np.zeros([n,n+m+1])], [J2 @ E, J2 @ F, J2 @ c] , [np.zeros([1,n+m+1])]      ])
JJ2 = (S2.T + S2)/2
S3 = cp.bmat([ [np.zeros([n,n+2*m+1])], [J3 @ E, J3 @ F, np.zeros([m,m]), J3 @ c] , [np.zeros([m+1,n+2*m+1])]      ])
JJ3 = (S3.T + S3)/2

#\lamd^T (Ex_plus + F\lamd + c) = 0
tau4 = cp.Variable((m,1))
J4 = cp.diag(tau4)
S4 = cp.bmat([ [np.zeros([n,n+2*m+1])], [J4 @ E @ A, J4 @ E @ D, J4 @ F, J4 @ E @ z + c] , [np.zeros([m+1,n+2*m+1])]      ])
JJ4 = (S4.T + S4)/2


#construct the LMI's
matl = np.zeros([n+m+1,n+m+1]) 
matl[0:n, 0:n] = eps*np.eye(n)
matk = np.zeros([n+2*m+1,n+2*m+1]) 
matk[0:n, 0:n] = eps*np.eye(n)
alpha = cp.Variable((1,1))
matz = np.zeros([n+m+1,n+m+1]) 
matz[0:n, 0:n] = np.eye(n)
constraints += [alpha >= eps]


#V > eps x^2
constraints += [ V - Vquad.T @ W @ Vquad - t1 * bnd - JJ1 -matl >> 0 ]

#DV < -eps x^2
constraints += [ -DV - DVquad.T @ U @ DVquad - DVquad2.T @ U_V2 @ DVquad2 - t2 * bnde - JJ3 - JJ4 - matk >> 0 ]

#V < alpha x^2
constraints += [-V - Vquad.T @ W2 @ Vquad - t3 * bnd - JJ2 + alpha * matz >> 0 ]

prob = cp.Problem(cp.Minimize(cp.square(c3)), constraints)

prob.solve(verbose = True)