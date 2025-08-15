import numpy as np
from math import exp, log

from model_DE_Parameters_Vectors import model_DE_Parameters_Vectors

def model_DE_Parameters_Control(larvicide_type, Tf):
    """
    Returns control parameters as a NumPy array.
    """
    pU = np.zeros(12)
    pV = model_DE_Parameters_Vectors()

    # Larvicide: 1 = Methoprene
    if larvicide_type == 1:
        min_ef = 0.03
        min_ef_day = 150
        half_ef_day = 100

        max_ef = 1 / (1 + (min_ef / (1 - min_ef)) ** (half_ef_day / (min_ef_day - half_ef_day)))
        pU[0] = max_ef * (pV[7] + pV[6]) / (1 - max_ef)  # km1
        pU[1] = -log((1 - max_ef) / max_ef) / half_ef_day  # decay rate

    # Larvicide: 2 = Vectobac
    elif larvicide_type == 2:
        min_ef = 0.22
        min_ef_day = 42
        mid_ef_day = 34
        mid_ef = 0.57

        part1 = (min_ef / (1 - min_ef)) ** (mid_ef_day / (min_ef_day - mid_ef_day))
        part2 = ((1 - mid_ef) / mid_ef) ** (min_ef_day / (min_ef_day - mid_ef_day))

        max_ef = 1 / (1 + part1 * part2)
        pU[0] = max_ef * (pV[7] + pV[6]) / (1 - max_ef)  # km1
        pU[1] = -log(mid_ef * (1 - max_ef) / ((1 - mid_ef) * max_ef)) / mid_ef_day  # decay rate

    # Adulticide decay and kill rates
    pU[2] = 24  # decay rate of adulticide (per hour?)
    pU[3] = -log(0.1) * pU[2] * 2  # max kill rate of adulticide (km2)

    # Percent of adulticide remaining after 1 hour (not used directly here)
    _ = exp(pU[2] * 0.5 * exp(-pU[3] / 24) / pU[3] - pU[2] * 0.5 / pU[3])

    # Cost weights and time controls
    pU[4] = 5000     # weight of infected components
    pU[5] = 1        # cost of larvicide
    pU[6] = 10       # cost of adulticide
    pU[7] = 0.05     # cost of time
    pU[8] = 5000     # cost of eggs at final time
    pU[9] = -100000  # cost of hosts at final time
    pU[10] = Tf      # max time between controls
    pU[11] = 1       # min time between controls

    return pU
