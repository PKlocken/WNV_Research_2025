import numpy as np
from math import log

def model_DE_Parameters_Hosts(urban, transmit, *args):
    """
    Returns all specified host parameters as a NumPy array.
    """
    if urban > 2 and urban < 5:
        if len(args) == 0:
            raise ValueError("Additional input required for urban > 2 cases.")
        tag_input = args[0]
        factor = 100
    elif urban == 5:
        if len(args) == 0:
            raise ValueError("Additional input required for urban = 5.")
        factor = args[0]
    else:
        tag_input = 0
        factor = 1

    # Host type 1: vulnerable amplifiers with horizontal transmission
    pH1 = np.zeros(10)
    pH1[0] = 0.68  # p_hm1 – host-to-mosquito transmission probability
    pH1[2] = (1/3)*((1/4.7)+(np.log(0.6)/365)) + np.log(0.6)/365  # g1 – WNV recovery rate
    pH1[3] = (1/4.7)+(np.log(0.6)/365)  # gamma1 – WNV induced death rate
    pH1[4] = (0.15 - np.log(0.6))/365  # Lambda1 – per capita birth rate
    pH1[5] = -np.log(0.6)/365  # mu_h1 – natural death rate
    pH1[1] = 3.5 * (pH1[5] + pH1[3] + pH1[2])  # omega1 – maximal horizontal transmission rate
    pH1[6] = 1     # c_h1 – carrying capacity (less urban)
    pH1[7] = 0.25  # c_h1 – carrying capacity (more urban)
    pH1[8] = 1     # alpha1 – biting preference
    pH1[9] = 1     # p_hh1 – probability of horizontal transmission

    # Host type 2: vulnerable amplifiers without horizontal transmission
    pH2 = np.zeros(10)
    pH2[0] = 0.53  # p_hm2
    pH2[1] = 0     # omega2
    pH2[2] = (1/4.7) + 2 * (np.log(0.6)/365)  # g2 – WNV recovery rate
    pH2[3] = (1/4.7) + (np.log(0.6)/365)      # gamma2 – induced death rate
    pH2[4] = (0.3 - np.log(0.6))/365          # Lambda2 – birth rate
    pH2[5] = -np.log(0.6)/365                 # mu_h2 – natural death rate
    pH2[6] = 2     # c_h2 less urban
    pH2[7] = 5     # c_h2 more urban
    pH2[8] = 0.1   # alpha2 – biting preference
    pH2[9] = 1     # p_hh2

    # Host type 3: invulnerable amplifiers
    pH3 = np.zeros(10)
    pH3[0] = 0.36  # p_hm3
    pH3[1] = 0     # omega3
    pH3[2] = 0.3333  # g3 – recovery rate
    pH3[3] = 0       # gamma3 – no induced death
    pH3[4] = (0.1570 - np.log(0.6))/365  # Lambda3
    pH3[5] = -np.log(0.6)/365            # mu_h3
    pH3[6] = 2     # c_h3 less urban
    pH3[7] = 1     # c_h3 more urban
    pH3[8] = 1     # alpha3
    pH3[9] = 1     # p_hh3

    # Host type 4: diluters
    pH4 = np.zeros(10)
    pH4[0] = 0.1   # p_hm4
    pH4[1] = 0     # omega4
    pH4[2] = 1     # g4 – total recovery
    pH4[3] = 0     # gamma4 – no death
    pH4[4] = (0.3 - np.log(0.6))/365  # Lambda4
    pH4[5] = -np.log(0.6)/365         # mu_h4
    pH4[6] = 2     # c_h4 less urban
    pH4[7] = 1     # c_h4 more urban
    pH4[8] = 0.1   # alpha4
    pH4[9] = 1     # p_hh4

    # Host type 5: dead-end hosts
    pH5 = np.zeros(10)
    pH5[0] = 0     # p_hm5
    pH5[1] = 0     # omega5
    pH5[2] = 0     # g5
    pH5[3] = 0     # gamma5
    pH5[4] = 0     # Lambda5
    pH5[5] = 0     # mu_h5
    pH5[6] = 70    # c_h5 less urban
    pH5[7] = 70    # c_h5 more urban
    pH5[8] = 0.05  # alpha5
    pH5[9] = 0     # p_hh5

    # Disable horizontal transmission for host 1 if transmit == 1
    if transmit == 1:
        pH1[1] = 0  # omega1 = 0

    # Disable host-to-mosquito transmission for host 3 in case 34
    if transmit == 34:
        pH3[0] = 0  # p_hm3 = 0

    # Disable host-to-mosquito transmission for host 4 in case 34 or 4
    if transmit == 34 or transmit == 4:
        pH4[0] = 0  # p_hm4 = 0

    # Modify alpha or carrying capacity under special urban cases
    if urban == 3:
        if tag_input == 1:
            pH3[6] = pH3[6] / factor  # c_h3 (less urban)
        elif tag_input == 2:
            pH3[8] = pH3[8] / factor  # alpha3
    elif urban == 4:
        if tag_input == 1:
            pH4[6] = pH4[6] * factor  # c_h4
        elif tag_input == 2:
            pH4[8] = pH4[8] * factor  # alpha4
    elif urban == 5:
        pH4[8] = 1                   # alpha4
        pH1[8] = 1 / factor          # alpha1
        pH2[8] = (1 / factor) * pH2[8]  # alpha2
        pH3[8] = (1 / factor) * pH3[8]  # alpha3

    hostParams = np.array([pH1, pH2, pH3, pH4, pH5])
    return hostParams
