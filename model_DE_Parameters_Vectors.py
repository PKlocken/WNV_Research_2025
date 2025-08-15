import numpy as np

def model_DE_Parameters_Vectors():
    """
    Returns vector parameters as a NumPy array.
    """
    pV = np.zeros(13)

    # Egg laying rates
    pV[0] = 150 / (2 * 8)   # rs
    pV[1] = 100 / (2 * 8)   # ri

    # Egg infection and hatching
    pV[2] = 0.003           # phi
    pV[3] = 0.56            # qs
    pV[4] = 0.43            # qi

    # Development and death rates
    pV[5] = 1 / 2           # m_E
    pV[6] = 1 / 7           # m_L
    pV[7] = 0.16            # mu_L
    pV[8] = 1 / 10.4        # mu_V
    pV[9] = 1 / 5           # b

    # Carrying capacity and disease progression
    pV[10] = 100            # c_L (per ha)
    pV[11] = 1 / 10         # kl
    pV[12] = 1              # p_mh

    return pV
