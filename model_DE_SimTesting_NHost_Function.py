import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.integrate import solve_ivp

from model_DE_Parameters_Vectors import model_DE_Parameters_Vectors
from model_DE_Parameters_Hosts import model_DE_Parameters_Hosts
from model_DE_R0_NextGen_Function import model_DE_R0_NextGen_Function
from model_DE_ForSim_NHost import model_DE_ForSim_NHost

def model_DE_SimTesting_NHost_Function(N, urban, transmit, infection_origin, num_years, savePath=None, extraInput=None, rtol=1e-3, atol=1e-6):
    
    # If savePath is not provided or is empty, create it with timestamp
    if not savePath:
        timestamp = datetime.now().strftime("%Y_%m_%d_%H%M%S")
        savePath = os.path.join("Plots", f"Run_{timestamp}")
        os.makedirs(savePath, exist_ok=True)
    elif not os.path.exists(savePath):
        os.makedirs(savePath, exist_ok=True)
    
     # DEBUG: confirm where the model will save and which tolerances it received
    print(f"MODEL DEBUG: savePath={os.path.abspath(savePath)!r}  rtol={rtol}  atol={atol}")

    if 2 < urban < 5:
        if extraInput is None:
            raise ValueError("extraInput is required when urban = 3 or 4.")
        tag_input = extraInput
        if tag_input == 1:
            tag = "cd"  # capacity dilution
        elif tag_input == 2:
            tag = "cb"  # biting preference
        else:
            tag = "cn"  # control/no-change
    else:
        tag_input = None
        tag = "cn"

    T_short_initial = 0
    T_short_final   = 250
    T_long_initial  = 250
    T_long_final    = num_years * 365
    tspan = [T_short_initial, T_long_final]

    # Load vector parameters
    pV = model_DE_Parameters_Vectors()
    rs, ri, phi, qs, qi, m_E, m_L, mu_L, mu_V, b, c_L, kl, p_mh = pV[:13]
    d_l = ((rs * m_L * qs / mu_V) - mu_L - m_L)

    # Compute disease-free equilibrium
    Es_DFE = c_L * m_L * rs / (mu_V * m_E)
    Ls_DFE = c_L
    Vs_DFE = c_L * m_L / mu_V

    # Initial vector condition based on infection_origin
    if infection_origin == 'e':
        x0 = [0.999*Es_DFE, 0.001*Es_DFE, Ls_DFE, 0, Vs_DFE, 0, 0]
    elif infection_origin == 'l':
        x0 = [Es_DFE, 0, 0.999*Ls_DFE, 0.001*Ls_DFE, Vs_DFE, 0, 0]
    elif infection_origin == 'M':
        x0 = [Es_DFE, 0, Ls_DFE, 0, 0.999*Vs_DFE, 0.001*Vs_DFE, 0]
    elif infection_origin == 'm':
        x0 = [Es_DFE, 0, Ls_DFE, 0, 0.999*Vs_DFE, 0, 0.001*Vs_DFE]
    else:
        x0 = [Es_DFE, 0, Ls_DFE, 0, Vs_DFE, 0, 0]

    # Load host parameters
    if urban <= 2:
        hostParams = model_DE_Parameters_Hosts(urban, transmit)
    else:
        hostParams = model_DE_Parameters_Hosts(urban, transmit, extraInput)

    # Initialize host compartments
    for j in range(N):
        pH = hostParams[j, :]
        if urban == 2:
            c_h = pH[7]
        else:
            c_h = pH[6]
        Hs_DFE = c_h
        if infection_origin.isnumeric() and int(infection_origin) == j+1:
            x0.extend([0.999*Hs_DFE, 0.001*Hs_DFE, 0])
        else:
            x0.extend([Hs_DFE, 0, 0])

    # Calculate R0 & Solve ODE
    if urban <= 2:
        R0, MH_ratio, MH_ratio_a, v_R0 = model_DE_R0_NextGen_Function(N, urban, transmit, None)
        sol = solve_ivp(lambda t, y: model_DE_ForSim_NHost(t, y, N, urban, transmit),
                        tspan, x0,
                        method='RK45',      # ode45 equivalent
                        rtol=rtol, atol=atol,
                        dense_output=True,
                        max_step=1.0)
    else:
        R0, MH_ratio, MH_ratio_a, v_R0 = model_DE_R0_NextGen_Function(N, urban, transmit, extraInput)
        sol = solve_ivp(lambda t, y: model_DE_ForSim_NHost(t, y, N, urban, transmit, extraInput),
                        tspan, x0,
                        method='RK45',      # ode45 equivalent
                        rtol=rtol, atol=atol,
                        dense_output=True,
                        max_step=1.0)

    t = sol.t
    x = sol.y.T

    # Time slicing
    short_idx = t <= T_short_final
    long_idx = t >= T_long_initial

    t_short = t[short_idx]
    x_short = x[short_idx, :]
    t_long  = t[long_idx]
    x_long  = x[long_idx, :]

    # Extract vector data, short time
    Vs_short, Ve_short, Vi_short = x_short[:, 4], x_short[:, 5], x_short[:, 6]
    Mosquito_Prevalence = Vi_short / (Vi_short + Ve_short + Vs_short)
    Max_Prevalence = Mosquito_Prevalence.max()

    # Host data
    hostData_short = x_short[:, 7:]
    hostData_long  = x_long[:, 7:]
    hostData_long_final = x_long[-1, 7:]

    # ------------------- PLOTTING -------------------

    # Fonts
    figureTitleFontSize = 25
    plotTitleFontSize = 25
    plotIndexFontSize = 20

    # Determine base tag string for filenames
    if transmit == 4:
        base_tag = 'no4__'
    elif transmit == 34:
        base_tag = 'no34__'
    elif transmit == 1:
        base_tag = 'no1__'
    else:
        base_tag = 'all__'

    # Infected compartments (short time)
    fig, axs = plt.subplots(2, 1, figsize=(16, 12), constrained_layout=True)

    for j in range(N):
        Hi = hostData_short[:, 3 * j + 1]
        axs[0].plot(t_short, Hi, linewidth=4)
    axs[0].set_title('Infected Hosts', fontsize=plotTitleFontSize)
    axs[0].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[0].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[0].legend([f'Type {j+1}' for j in range(N)], fontsize=plotIndexFontSize)
    axs[0].grid(True)

    axs[1].plot(t_short, Ve_short, 'm', linewidth=4, label='Exposed Mosquitoes')
    axs[1].plot(t_short, Vi_short, 'r', linewidth=4, label='Infected Mosquitoes')
    axs[1].set_title('Infected Vectors', fontsize=plotTitleFontSize)
    axs[1].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[1].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[1].legend(fontsize=plotIndexFontSize)
    axs[1].grid(True)

    file_base = f'infected_compartments_T0={T_short_initial}_Tf={T_short_final}_N={N}_urban={urban}_{base_tag}R0={R0:.2f}_{infection_origin}_{tag}'
    fig.savefig(os.path.join(savePath, file_base + '.eps'), format='eps')
    fig.savefig(os.path.join(savePath, file_base + '.png'), format='png')
    plt.close(fig)

    # Infected compartments (long time)
    fig, axs = plt.subplots(2, 1, figsize=(16, 12), constrained_layout=True)

    for j in range(N):
        Hi = hostData_long[:, 3 * j + 1]
        axs[0].plot(t_long, Hi, linewidth=4)
    axs[0].set_xlim([T_long_initial, t_long[-1]])
    axs[0].set_ylim(bottom=0)
    axs[0].set_title('Infected Hosts', fontsize=plotTitleFontSize)
    axs[0].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[0].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[0].legend([f'Type {j+1}' for j in range(N)], fontsize=plotIndexFontSize)
    axs[0].grid(True)

    axs[1].plot(t_long, Ve_long := x_long[:, 5], 'm', linewidth=4, label='Exposed Mosquitoes')
    axs[1].plot(t_long, Vi_long := x_long[:, 6], 'r', linewidth=4, label='Infected Mosquitoes')
    axs[1].set_xlim([T_long_initial, t_long[-1]])
    axs[1].set_ylim(bottom=0)
    axs[1].set_title('Infected Vectors', fontsize=plotTitleFontSize)
    axs[1].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[1].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[1].legend(fontsize=plotIndexFontSize)
    axs[1].grid(True)

    file_base = f'infected_compartments_T0={T_long_initial}_Tf={T_long_final}_N={N}_urban={urban}_{base_tag}R0={R0:.2f}_{infection_origin}_{tag}'
    fig.savefig(os.path.join(savePath, file_base + '.eps'), format='eps')
    fig.savefig(os.path.join(savePath, file_base + '.png'), format='png')
    plt.close(fig)

    # Susceptible and Recovered hosts (short time)
    fig, axs = plt.subplots(2, 1, figsize=(16, 12), constrained_layout=True)

    for j in range(N):
        Hs = hostData_short[:, 3 * j]
        axs[0].plot(t_short, Hs, linewidth=4)
    axs[0].set_title('Susceptible Hosts', fontsize=plotTitleFontSize)
    axs[0].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[0].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[0].legend([f'Type {j+1}' for j in range(N)], fontsize=plotIndexFontSize)
    axs[0].grid(True)

    for j in range(N):
        Hr = hostData_short[:, 3 * j + 2]
        axs[1].plot(t_short, Hr, linewidth=4)
    axs[1].set_title('Recovered Hosts', fontsize=plotTitleFontSize)
    axs[1].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[1].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[1].legend([f'Type {j+1}' for j in range(N)], fontsize=plotIndexFontSize)
    axs[1].grid(True)

    file_base = f'recovered_compartments_T0={T_short_initial}_Tf={T_short_final}_N={N}_urban={urban}_{base_tag}R0={R0:.2f}_{infection_origin}_{tag}'
    fig.savefig(os.path.join(savePath, file_base + '.eps'), format='eps')
    fig.savefig(os.path.join(savePath, file_base + '.png'), format='png')
    plt.close(fig)

    # Susceptible and Recovered hosts (long time)
    fig, axs = plt.subplots(2, 1, figsize=(16, 12), constrained_layout=True)

    for j in range(N):
        Hs = hostData_long[:, 3 * j]
        axs[0].plot(t_long, Hs, linewidth=4)
    axs[0].set_xlim([T_long_initial, t_long[-1]])
    axs[0].set_ylim(bottom=0)
    axs[0].set_title('Susceptible Hosts', fontsize=plotTitleFontSize)
    axs[0].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[0].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[0].legend([f'Type {j+1}' for j in range(N)], fontsize=plotIndexFontSize)
    axs[0].grid(True)

    for j in range(N):
        Hr = hostData_long[:, 3 * j + 2]
        axs[1].plot(t_long, Hr, linewidth=4)
    axs[1].set_xlim([T_long_initial, t_long[-1]])
    axs[1].set_ylim(bottom=0)
    axs[1].set_title('Recovered Hosts', fontsize=plotTitleFontSize)
    axs[1].set_xlabel('Time (days)', fontsize=plotIndexFontSize)
    axs[1].set_ylabel('Population (ind/ha)', fontsize=plotIndexFontSize)
    axs[1].legend([f'Type {j+1}' for j in range(N)], fontsize=plotIndexFontSize)
    axs[1].grid(True)

    file_base = f'noninfected_compartments_T0={T_long_initial}_Tf={T_long_final}_N={N}_urban={urban}_{base_tag}R0={R0:.2f}_{infection_origin}_{tag}'
    fig.savefig(os.path.join(savePath, file_base + '.eps'), format='eps')
    fig.savefig(os.path.join(savePath, file_base + '.png'), format='png')
    plt.close(fig)

    # ------------------- TXT OUTPUT -------------------
    if isinstance(extraInput, (int, float)):
        extraStr = str(extraInput)
    else:
        extraStr = str(extraInput)

    filename = f"SimResults_Urban{urban}_Trans{transmit:.2f}_{tag}_Extra{extraStr}.txt"
    txtFileName = os.path.join(savePath, filename)
    with open(txtFileName, 'w') as f:
        f.write("Simulation Summary\n")
        f.write("==================\n\n")
        f.write(f"R0: {R0:.6f}\n")
        f.write("Corresponding Eigenvector: " + " ".join(f"{val:.6f}" for val in v_R0) + "\n")
        f.write(f"Max Mosquito Prevalence: {Max_Prevalence:.6f}\n")
        f.write(f"Mosquito-to-Host Ratio (MH_ratio): {MH_ratio:.6f}\n")
        f.write(f"Adjusted MH_ratio (MH_ratio_a): {MH_ratio_a:.6f}\n\n")
        f.write("Final Host Compartment Values:\n")
        for j, val in enumerate(hostData_long_final, start=1):
            f.write(f"  x({j}): {val:.6f}\n")

    return hostData_long_final, Max_Prevalence, R0, MH_ratio, MH_ratio_a, v_R0
