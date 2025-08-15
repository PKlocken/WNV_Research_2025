import numpy as np
from numpy.linalg import eig, inv

from model_DE_Parameters_Vectors import model_DE_Parameters_Vectors
from model_DE_Parameters_Hosts import model_DE_Parameters_Hosts

def model_DE_R0_NextGen_Function(N, urban, transmit, extraInput=None):
    # Vector Parameters
    pV = model_DE_Parameters_Vectors()
    rs, ri, phi, qs, qi, m_E, m_L, mu_L, mu_V, b, c_L, k_L, p_mh = pV[:13]

    # Derived Vector Parameters
    d_l = ((rs * m_L * qs / mu_V) - mu_L - m_L)
    tau = m_L + mu_L + d_l
    M_star = m_L * c_L / mu_V  # DFE susceptible adult mosquitoes

    # Vector next-gen matrix blocks
    W_M = np.array([
        [m_E, 0, 0, 0],
        [-m_E * qi * phi, tau, 0, 0],
        [0, 0, k_L + mu_V, 0],
        [0, -m_L, -k_L, mu_V]
    ])

    F_MM = np.array([
        [0, 0, 0, ri],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]
    ])

    # Host Parameters
    if urban <= 2:
        hostParams = model_DE_Parameters_Hosts(urban, transmit)
    else:
        hostParams = model_DE_Parameters_Hosts(urban, transmit, extraInput)

    s = []
    omega_list = []
    f_MH = []
    f_HM = []

    NH = 0
    NHa = 0

    for j in range(N):
        params = hostParams[j, :]

        p_hm = params[0]
        omega_j = params[1]
        g = params[2]
        gamma = params[3]
        Lambda = params[4]
        mu_h = params[5]

        if urban == 2:
            c_h = params[7]
        else:
            c_h = params[6]

        alpha = params[8]

        d_h = Lambda - mu_h
        sj = gamma + g + mu_h + d_h

        s.append(sj)
        omega_list.append(omega_j)
        f_MH.append(b * p_mh * alpha * c_h)
        f_HM.append(b * p_hm * alpha * M_star)

        NH += alpha * c_h
        NHa += c_h

    # Dead-end host (row 5)
    params = hostParams[4, :]
    c_h_dead = params[7] if urban == 2 else params[6]
    alpha_dead = params[8]

    NH += alpha_dead * c_h_dead
    NHa += c_h_dead

    # Ratios
    MH_ratio = M_star / NH
    MH_ratio_a = M_star / NHa

    # Host matrices
    W_H = np.diag(s)
    F_HH = np.diag(omega_list)
    F_HM_mat = np.zeros((4, N))
    F_MH_mat = np.zeros((N, 4))

    F_HM_mat[2, :] = np.array(f_HM) / NH
    F_MH_mat[:, 3] = np.array(f_MH) / NH

    # Assemble full matrices
    W_H_inv = inv(W_H)
    W_M_inv = inv(W_M)

    top_left = F_HH @ W_H_inv
    top_right = F_MH_mat @ W_M_inv
    bottom_left = F_HM_mat @ W_H_inv
    bottom_right = F_MM @ W_M_inv

    FWinv = np.block([
        [top_left, top_right],
        [bottom_left, bottom_right]
    ])

    # Spectral radius of next generation matrix
    eigVals, eigVecs = eig(FWinv)
    eigVals_diag = np.real(eigVals)
    R0 = np.max(eigVals_diag)

    idx = np.argmax(eigVals_diag)
    v_R0 = eigVecs[:, idx]
    v_R0 = v_R0 / np.max(np.abs(v_R0))  # Normalize

    return R0, MH_ratio, MH_ratio_a, v_R0
