import math
import numpy as np


def permutations_k(k):
    #p_1
    row_1 = [0, -1, 0, 0]
    row_2 = [1, 0, 0, 0]
    row_3 = [0, 0, 0, 1]
    row_4 = [0, 0, -1, 0]
    p_1 = np.array([row_1, row_2, row_3, row_4])

    #p_2
    row_1 = [0, 0, -1, 0]
    row_2 = [0, 0, 0, -1]
    row_3 = [1, 0, 0, 0]
    row_4 = [0, 1, 0, 0]
    p_2 = np.array([row_1, row_2, row_3, row_4])

    #p_3
    row_1 = [0, 0, 0, -1]
    row_2 = [0, 0, 1, 0]
    row_3 = [0, -1, 0, 0]
    row_4 = [1, 0, 0, 0]
    p_3 = np.array([row_1, row_2, row_3, row_4])

    p = np.array([p_1, p_2, p_3])

    return p[k - 1]


def zeta_k(k, m_inertia, p, q):
    return 1/(4*m_inertia[k])*(np.transpose(p) @ permutations_k(k) @ q)

def propagator_k(k, m_inertia, p, q, dt):
    #external functions
    zeta = zeta_k(k, m_inertia, p, q)
    p_k = permutations_k(k)
    #propagator analytic expansion
    q_new = math.cos(zeta*dt)*q + math.sin(zeta*dt)*(p_k @ q)
    p_new = math.cos(zeta*dt)*p + math.sin(zeta*dt)*(p_k @ p)

    return [q_new, p_new]

def propagator_p(p, q_torque, dt):
    return p + dt/2*q_torque
