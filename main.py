import os
from datetime import datetime
from itertools import product

from model_DE_Sim_Wrapper import run_all_simulations as run_normal
from model_DE_ShortTime_Tolerance_Wrapper import run_all_simulations as run_tolerance


def tolerance_experiment_runs(N, num_years, urban_list, transmit_list, origin_list, rtol_list, atol_list):
    timestamp = datetime.now().strftime("%Y_%m_%d_%H%M%S")
    base_dir = os.path.join('Plots', f'Run_{timestamp}_rtol_atol_experiments')
    os.makedirs(base_dir, exist_ok=True)

    # Loop over all combinations of rtol and atol
    for rtol, atol in product(rtol_list, atol_list):
        # Compose descriptive subfolder name (do NOT create this folder here!)
        subdir_name = f"rtol_{rtol}_atol_{atol}"

        print(f"Running simulation with rtol={rtol}, atol={atol}...\n")

        run_tolerance(
            N=N,
            num_years=num_years,
            urban_list=urban_list,
            transmit_list=transmit_list,
            origin_list=origin_list,
            rtol=rtol,
            atol=atol,
            base_dir=base_dir,
            extra_subdir=subdir_name  # Just a string, folder created inside run_all_simulations
        )


if __name__ == "__main__":

    """
    Scenario variable explanation:
    N = 4  # Number of host types
    num_years = 5  # Simulation length in years
    urban_list = [1, 2, 3, 4, 5]  # Urban settings
    transmit_list = [0, 1, 34, 4]  # Transmission restriction codes
    origin_list = ['e', 'm', '1', '2', '3', '4']  # Infection origin identifiers
    rtol_list = [0.1, 0.05, 0.01, 0.005, 0.001] # List of relative tolerances to try
    atol_list = [1e-6, 1e-7] # List of absolute tolerances to try
    
    """
    

    N = 4
    num_years = 5
    urban_list = [1]
    origin_list = ['m']
    transmit_list = [0]
    
    """ run_all_simulations(N, num_years, urban, transmit, origin) """
    
    """ rtol = 0.001
    atol = 1e-6
    run_all_simulations(N, num_years, urban_list, transmit_list, origin_list, rtol, atol) """    


    #rtol_list = [0.1, 0.05, 0.01, 0.005, 0.001] # List of relative tolerances to try
    #atol_list = [1e-4, 5*1e-5, 1e-5, 1e-6, 1e-7] # List of absolute tolerances to try
    rtol_list = [0.001] # List of relative tolerances to try
    atol_list = [1e-8, 1e-9, 1e-10] # List of absolute tolerances to try
    tolerance_experiment_runs(N, num_years, urban_list, transmit_list, origin_list, rtol_list, atol_list)