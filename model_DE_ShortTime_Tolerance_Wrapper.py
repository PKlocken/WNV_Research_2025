import os
from datetime import datetime
from itertools import product

from model_DE_ShortTimeOnly_Simulation_NHost_Function import model_DE_ShortTimeOnly_Simulation_NHost_Function

def run_all_simulations(N, num_years, urban_list, transmit_list, origin_list, rtol, atol, base_dir=None, extra_subdir=None):
    if base_dir is None:
        timestamp = datetime.now().strftime("%Y_%m_%d_%H%M%S")
        base_dir = os.path.join('Plots', f'Run_{timestamp}')
        os.makedirs(base_dir, exist_ok=True)

    if extra_subdir is not None:
        base_dir = os.path.join(base_dir, extra_subdir)
        os.makedirs(base_dir, exist_ok=True)
    
    print(f"Saving results to: {base_dir}\n")

    # Loop through all combinations of parameters
    print("Starting simulation runs...\n")
    for urban, transmit, origin in product(urban_list, transmit_list, origin_list):
        # Create a subdirectory for each parameter combination
        combo_subdir = os.path.join(base_dir, f'urban_{urban}__transmission_{transmit}__origin_{origin}')
        os.makedirs(combo_subdir, exist_ok=True)

        print(f'Running simulation with Urban={urban}, Transmission={transmit}, Origin={origin}\n')

        # <-- DEBUG PRINT: confirm the exact folder being passed to the model
        print("WRAPPER DEBUG: passing savePath =", os.path.abspath(combo_subdir), f"  rtol={rtol}  atol={atol}")

        try:
            # Call your simulation function with current parameters
            model_DE_ShortTimeOnly_Simulation_NHost_Function(
                N=N,
                urban=urban,
                transmit=transmit,
                infection_origin=origin,
                num_years=num_years,
                savePath=combo_subdir,
                rtol=rtol,
                atol=atol
            )
        except Exception as e:
            print(f"Error during simulation for Urban={urban}, Transmission={transmit}, Origin={origin}: {e}\n")

    print(f'\nAll model runs completed.\nResults saved in directory: {base_dir}')