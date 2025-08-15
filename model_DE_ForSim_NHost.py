import numpy as np
from model_DE_Parameters_Vectors import model_DE_Parameters_Vectors
from model_DE_Parameters_Control import model_DE_Parameters_Control
from model_DE_Parameters_Hosts import model_DE_Parameters_Hosts

def model_DE_ForSim_NHost(t, x, N, urban, transmit, *args):
    # Check for urban > 2 requiring extra input
    if 2 < urban < 5:
        if not args:
            raise ValueError("Additional input required for urban > 2 cases.")
        tag_input = args[0]
        factor = 100
    else:
        tag_input = None

    # Vector parameters
    pV = model_DE_Parameters_Vectors()
    rs, ri, phi, qs, qi, m_E, m_L, mu_L, mu_V, b, c_L, kl, p_mh = pV[:13]

    d_l = ((rs * m_L * qs / mu_V) - mu_L - m_L)  # density-dependent larval death

    # State variables
    Es, Ei, Ls, Li, Vs, Ve, Vi = x[:7]

    # Control parameters
    pU = model_DE_Parameters_Control(1, 10)
    km1 = pU[0]
    km2 = pU[3]

    # Control state variables (set to 0)
    Ul = 0
    Ua = 0

    # Host parameters
    hostParams = model_DE_Parameters_Hosts(urban, transmit, *args)

    NH = np.zeros(N)
    Y = np.zeros(N)
    alpha = np.zeros(N)

    dH = []

    dVs = m_L * Ls
    dVe = 0

    # Loop through host types
    for j in range(N):
        pH = hostParams[j, :]

        alpha[j] = pH[8]  # biting preference

        Hs = x[7 + 3 * j]
        Hi = x[8 + 3 * j]
        Hr = x[9 + 3 * j]

        NH[j] = Hs + Hi + Hr
        Y[j] = alpha[j] * NH[j]

    NY = np.sum(Y)

    # Dead-end host contribution
    params = hostParams[4, :]
    if urban == 2:
        c_h = params[7]
    else:
        c_h = params[6]
    alphaD = params[8]
    NY += alphaD * c_h

    # Host ODEs and vector interactions
    for j in range(N):
        pH = hostParams[j, :]

        p_hm, omega, g, gamma, Lambda, mu_h = pH[:6]
        if urban == 2:
            c_h = pH[7]
        else:
            c_h = pH[6]
        p_hh = pH[9]

        d_h = Lambda - mu_h

        Hs = x[7 + 3 * j]
        Hi = x[8 + 3 * j]
        Hr = x[9 + 3 * j]
        NHj = NH[j]

        dHs = Lambda * NHj - b * p_mh * Vi * alpha[j] * Hs / NY \
              - omega * p_hh * Hi * Hs / NHj - d_h * NHj * Hs / c_h - mu_h * Hs

        dHi = b * p_mh * Vi * alpha[j] * Hs / NY + omega * p_hh * Hi * Hs / NHj \
              - (gamma + g) * Hi - d_h * NHj * Hi / c_h - mu_h * Hi

        dHr = g * Hi - d_h * NHj * Hr / c_h - mu_h * Hr

        dH.extend([dHs, dHi, dHr])

        dVs -= b * p_hm * Vs * alpha[j] * Hi / NY
        dVe += b * p_hm * Vs * alpha[j] * Hi / NY

    # Final vector equations
    dVs -= mu_V * Vs + km2 * Vs * Ua
    dVe -= kl * Ve + mu_V * Ve + km2 * Ve * Ua

    # Other vector dynamics
    dEs = rs * (Vs + Ve) - m_E * Es
    dEi = ri * Vi - m_E * Ei

    dLs = m_E * qs * Es + m_E * qi * (1 - phi) * Ei - mu_L * Ls - m_L * Ls \
          - d_l * Ls * (Ls + Li) / c_L - km1 * Ls * Ul

    dLi = m_E * qi * phi * Ei - mu_L * Li - m_L * Li \
          - d_l * Li * (Li + Ls) / c_L - km1 * Li * Ul

    dVi = m_L * Li + kl * Ve - mu_V * Vi - km2 * Vi * Ua

    # Combine all derivatives
    dxdt = np.array([dEs, dEi, dLs, dLi, dVs, dVe, dVi] + dH)

    return dxdt
