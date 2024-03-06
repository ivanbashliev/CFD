# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:28:25 2024

@author: Ivan Bashliev
"""

import numpy as np


#######################################################################################
#                                   FUNCTIONS                                         #
#######################################################################################

def incidence_tE21(N):
    
    tE21 = np.zeros((N**2 + 4*N, 2*N*(N+1) + 4*N), dtype=int)
    alpha = 0
    beta = 0
    gamma = 0
    ind = N 
    
    
    for i in range(N**2+4*N):
        if N <= i < N**2 + 3*N:
            tE21[i,i-N+beta] = -1 
            tE21[i,i-N+1+beta] = 1
            alpha += 1
            if alpha == N + 2:
                beta += 1
                alpha = 0
        if N == 1:
            if i != N and i != 2*N + 1:
                tE21[i,gamma + N*(N+1) + 2*N] = -1
                tE21[i,gamma + N*(N+1) + 3*N] = 1
                gamma += 1
        elif N == 2:
            if i != N and i != 2*N + 1 and i != 2*N + 2 and i != 3*N + 3:
                tE21[i,gamma + N*(N+1) + 2*N] = -1
                tE21[i,gamma + N*(N+1) + 3*N] = 1
                gamma += 1
        elif N == 3:
            if i != N and i != 2*N + 1 and i!= 2*N+2 and i!= 3*N + 3 and i!=3*N+4 and i!=4*N+5:
                tE21[i,gamma + N*(N+1) + 2*N] = -1
                tE21[i,gamma + N*(N+1) + 3*N] = 1
                gamma += 1
        elif N == 4:
            if i != N and i != 2*N + 1 and i!= 2*N+2 and i!= 3*N + 3 and i!=3*N+4 and i!=4*N+5 and i!=4*N+6 and i!=5*N+7:
                tE21[i,gamma + N*(N+1) + 2*N] = -1
                tE21[i,gamma + N*(N+1) + 3*N] = 1
                gamma += 1
    return tE21


#######################################################################################
#                                     INPUT                                           #
#######################################################################################

N = 4               # Number of points per dimension
tol = 1e-6            # Convergence criterion
one = int(1)
mone = int(-1)
L = float(1.0)        # Domain length
Re = float(1000)      # Reynolds number

#  Boundary conditions for the lid driven acvity test case
U_wall_top = -1
U_wall_bot = 0
U_wall_left = 0
U_wall_right = 0
V_wall_top = 0
V_wall_bot = 0
V_wall_left = 0
V_wall_right = 0 

############################## Main body ##############################################


u = np.zeros([2*N*(N+1),1], dtype = float)
p = np.zeros([N*N+4*N,1], dtype = float)
tx = np.zeros([N+1,1], dtype = float)     # grid points on primal grid
x = np.zeros([N+2,1], dtype = float)      # grid points on dual grid
th = np.zeros([N], dtype = float)       # mesh width primal grid
h = np.zeros([N+1], dtype = float)      # mesh width dual grid 


# Set up the sparse incidence matrix tE21. Use the orientations described in the assignment.
tE21 = incidence_tE21(N)