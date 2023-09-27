"""

Analyse data from Monte Carlo simulations. Takes in energy 
and magnetisation output files.

"""

import sys
import os
from matplotlib import pyplot as plt
import numpy as np
import scipy

# function definitions

def get_heat_caps_diff(temps, energies):

    """

    Find heat capacities by differentiation of internal energy values
    at different temperatures:
    Cv = dE/dT
    
    Parameters:
    temps: list of temperature values
    energies: list of energy values for each respective temperature

    Returns:
    heat_caps: list of heat capacities for each temperature

    """
    
    # initiate heat capacity array
    heat_caps = []
    # number of energy values, i.e. number of temps
    n_energies = len(energies)

    # iterate over energy values
    for i in range(n_energies):
        # difference in energy
        d_energy = 0
        # difference in temperature
        d_temp = 0
        
        # check if there is more than one value
        if n_energies > 1:
            # first value
            if i == 0:
                # use only temp one above
                d_energy = energies[0] - energies[1]
                d_temp = temps[0] - temps[1]
            # last value
            elif i == n_energies-1:
                # use only temp one below
                d_energy = energies[-2] - energies[-1]
                d_temp = temps[-2] - temps[-1]
            # somewhere in the middle
            else:
                d_energy = energies[i-1] - energies[i+1]
                d_temp = temps[i-1] - temps[i+1]
        # if only one or no energy value
        # set temp to 1 to avoid division by zero error
        # but let heat cap be zero
        else:
            d_temp = 1

        # heat capacity
        heat_cap = d_energy/d_temp
        heat_caps.append(heat_cap)

    return heat_caps


def get_heat_caps_var(temps, energies, energies_sq):

    """

    Find heat capacities by differentiation of internal energy values
    at different temperatures:
    Cv/k = 1/(kT)^2*(mean(E^2) - mean(E)^2)

    Parameters:
    temps: list of temperature values (kT)
    energies: list of energy values for each respective temperature
    energies_sq: list of energy values squared for each respective temperature

    Returns:
    heat_caps: list of heat capacities for each temperature

    """

    # initiate heat capacity array
    heat_caps = []
    # number of energy values, i.e. number of temps
    n_energies = len(energies)

    # iterate over energy values
    for i in range(n_energies):
        # heat capacity
        heat_cap = 1/(temps[i]**2)*(energies_sq[i] - energies[i]**2)
        heat_caps.append(heat_cap)

    return heat_caps


def get_trans_temp(heat_caps):

    """

    Find transition temp from peak in heat capacity vs temp.
    Returns index of highest peak in heat capacity.

    """

    # find list of peaks
    peaks = scipy.signal.find_peaks(heat_caps)

    # peaks found
    if len(peaks[0]) > 0:
        # find highest peak
        highest_peak = peaks[0][0]
        for p in peaks[0]:
            if heat_caps[p] > heat_caps[highest_peak]:
                highest_peak = p
    # no peaks found
    else:
        highest_peak = False

    return highest_peak


def main():

    # read in input dir
    in_dir = sys.argv[1]
    # read in u0 subdirectory prefix
    u_subdir_prefix = sys.argv[2]
    # read in output directory
    out_dir = sys.argv[3]
    # read in total lattice size (e.g. 8000 for 20x20x20)
    lattice_total = float(sys.argv[4])
    # read in data that won't be considered
    if len(sys.argv) == 6:
        excluded_dirs = sys.argv[5].split(",")
    else:
        excluded_dirs = []

    # dictionary for all u0 values at all temps
    # containing energy values and heat capacities
    # hierarchy:
    #   data = {
    #       "u0": {
    #           "temps": [ T1, T2, T3, ... ],
    #           "energies": [ E(T1), E(T2), ... ],
    #           "heat_capacities": [ C(T1), C(T2), ... ],
    #           "magnetisation": [ M(T1), M(T2), ... ]
    #       }
    #   }
    data = {}

    # set up pyplot
    plt.figure(figsize=(6,4), dpi=220)
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Computer Modern",
        "font.size": 14
    })

    # walk through u0 subdirectories
    for u_subdir in os.listdir(in_dir):
        # if directory is really a directory, name contains set prefix, and output exists
        if os.path.isdir(os.path.join(in_dir, u_subdir)) \
                and (u_subdir.find(u_subdir_prefix) == 0) \
                and os.path.exists(os.path.join(in_dir, u_subdir, "energy.txt")) \
                and (u_subdir.replace(u_subdir_prefix, "") not in excluded_dirs):
            # determine value of u0
            u_value_str = u_subdir.replace(u_subdir_prefix, "")
            u_value = float(u_value_str)
            
            data[u_value_str] = { "temps": [], "energies": [], "energies_sq":[], "heat_caps": [], 
                    "magnetisation": [], "last_magn": 0, "last_magn_stdev": 0  }

            # get energies from u0 directory at respective temps
            energies_file_path = os.path.join(in_dir, u_subdir, "energy.txt")
            energies_file = open(energies_file_path, "r")
            for line in energies_file.readlines():
                chunks = line.split()
                temp = float(chunks[0])
                energy = float(chunks[1])
                data[u_value_str]["temps"].append(temp)
                data[u_value_str]["energies"].append(energy)

            # get energies squared from u0 directory at respective temps
            energies_sq_file_path = os.path.join(in_dir, u_subdir, "energy2.txt")
            energies_sq_file = open(energies_sq_file_path, "r")
            for line in energies_sq_file.readlines():
                chunks = line.split()
                energy_sq = float(chunks[1])
                data[u_value_str]["energies_sq"].append(energy_sq)

            # get magnetisation from u0 directory at respective temps
            magn_file_path = os.path.join(in_dir, u_subdir, "total_magnetisation.txt")
            #magn_file_path = os.path.join(in_dir, u_subdir, "spin_total.txt")
            magn_file = open(magn_file_path, "r")
            for line in magn_file.readlines():
                chunks = line.split()
                if len(chunks) == 2:
                    # is 1D
                    dim = 1
                    #magn = abs(float(chunks[1]))/lattice_total
                    magn = abs(float(chunks[1]))
                elif len(chunks) == 4:
                    # is 3D, so take vector magnitude
                    dim = 3
                    magn = np.sqrt(float(chunks[1])**2 + float(chunks[2])**2 \
                            + float(chunks[3])**2)/lattice_total
                data[u_value_str]["magnetisation"].append(magn)
            print(u_value_str)
            print(len(data[u_value_str]["temps"]))
            # get average spin magnitude and deviation from spins_after.txt
            ''' spin_file_path = os.path.join(in_dir, u_subdir, "spins_after.txt")
            spin_file = open(spin_file_path, "r")
            spin_lines = spin_file.readlines()
            
            # walk through lines and extract spins
            spins = []
            i = 0
            while i < len(spin_lines):
                if dim == 1:
                    # read only one line per spin
                    spin = abs(float(spin_lines[i]))
                    i += 1
                elif dim == 3:
                    # read three lines per spin
                    #spin = np.sqrt(float(spin_lines[i])**2 + float(spin_lines[i+1])**2 
                    #      + float(spin_lines[i+2])**2)
                    #i += 3
                    spin = 0
                    i += 1
                else:
                    i += 1
                spins.append(spin)
                
            avg_spin = np.mean(spins)
            dev_spins = np.std(spins)
            '''
            avg_spin = 1
            dev_spins = 1
            data[u_value_str]["last_magn"] = avg_spin
            data[u_value_str]["last_magn_stdev"] = dev_spins

    # go through all u0 values
    for u, u_data in data.items():
        heat_caps = get_heat_caps_diff(u_data["temps"], u_data["energies"])
        heat_caps_var = get_heat_caps_var(u_data["temps"], u_data["energies"], u_data["energies_sq"])
        # store heat capacities and convert to per lattice site
        data[u]["heat_caps"] = [ i/lattice_total for i in heat_caps ]
        data[u]["heat_caps_var"] = [ i/lattice_total for i in heat_caps_var ]

        # plot heat capacity vs temp
        plt.clf()
        plt.xlabel("Temperature $kT$")
        plt.ylabel("Heat capacity $C_V/k$")
        
        # TEST
        # plt.xlim((1,1.5))
        
        plt.plot(data[u]["temps"], data[u]["heat_caps"], label="differentiation")
        plt.plot(data[u]["temps"], data[u]["heat_caps_var"], label="variance", color="#1a6ece")
        plt.legend()
        plt.savefig(os.path.join(out_dir, ("heat_cap_vs_temp_u0_"+u+".png")))
        
        # write heat caps to output file
        with open(os.path.join(out_dir, "heat_cap_vs_temp_u0_"+u+".csv"), "w") as cv_out:
            # CSV header
            cv_out.write("temp,heat_cap\n")
            # CSV data
            for i in range(len(data[u]["temps"])):
                temp = data[u]["temps"][i]
                heat_cap = data[u]["heat_caps"][i]
                cv_out.write("{},{}\n".format(temp, heat_cap))

        # write variance heat caps to file
        with open(os.path.join(out_dir, "heat_cap_vs_temp_u0_"+u+"_var.csv"), "w") as cv_out:
            # CSV header
            cv_out.write("temp,heat_cap\n")
            # CSV data
            for i in range(len(data[u]["temps"])):
                temp = data[u]["temps"][i]
                heat_cap = data[u]["heat_caps_var"][i]
                cv_out.write("{},{}\n".format(temp, heat_cap))

        # plot magnetisation vs temp
        plt.clf()
        plt.xlabel("Temperature kT")
        print(u)
        plt.ylabel("Magnetisation per site")
        plt.plot(data[u]["temps"], data[u]["magnetisation"], color="#1a6ece")
        plt.savefig(os.path.join(out_dir, ("magn_vs_temp_u0_"+u+".png")))

       	# write magnetisations to output file
        with open(os.path.join(out_dir, "magn_vs_temp_u0_"+u+".csv"), "w") as magn_out:
            # CSV header
            magn_out.write("temp,magn\n")
            # CSV data
            for i in range(len(data[u]["temps"])):
                temp = data[u]["temps"][i]
                magn = data[u]["magnetisation"][i]
                magn_out.write("{},{}\n".format(temp, magn)) 

        print("u0 = " + u)
        print(data[u]["last_magn"])
        print(data[u]["last_magn_stdev"])
        print()
    
    # subfigures
    fig, (ax_cv, ax_m) = plt.subplots(2, sharex=True, figsize=(6,8), dpi=220)

    # plot heat capacity vs temp for several u in one plot
    #plt.clf()
    plt.xlabel("Temperature $kT$")
    #plt.ylabel("Heat capacity $C_V/k$")
    ax_cv.set(ylabel="Heat capacity $C_V/k$")
    u_selected = [ "1.0", "2.0", "6.0", "30", "100" ]
    colours = [ "#1a6ece", "#9ecbff", "#e85959", "#bd2929", "#6d0b0b" ]
    for i in range(len(u_selected)):
        u = u_selected[i]
        if u in data:
            ax_cv.plot(data[u]["temps"], data[u]["heat_caps_var"], label="$u_0 = {}$".format(u).replace(".0",""), 
                    linewidth=1, color=colours[i])

    #plt.legend()
    #plt.tight_layout()
    #plt.savefig(os.path.join(out_dir, "heat_caps_vs_temp_several_u0_var.png"))

    # plot heat capacity vs temp for several u in one plot
    #plt.clf()
    #plt.xlabel("Temperature $kT$")
    #plt.
    ax_m.set(ylabel="Net magnetisation $\\vert\\langle m \\rangle\\vert$")
    u_selected = [ "1.0", "2.0", "6.0", "30", "100" ]
    colours = [ "#1a6ece", "#9ecbff", "#e85959", "#bd2929", "#6d0b0b" ]
    for i in range(len(u_selected)):
        u = u_selected[i]
        if u in data:
            ax_m.plot(data[u]["temps"], data[u]["magnetisation"], label="$u_0 = {}$".format(u).replace(".0",""), 
                    linewidth=1, color=colours[i])
    
    ax_cv.tick_params(axis="both", direction="in")
    ax_m.tick_params(axis="both", direction="in")
    fig.align_ylabels([ ax_cv, ax_m ])
    ax_cv.legend(loc="upper left")
    ax_cv.legend()
    fig.tight_layout()
    #plt.savefig(os.path.join(out_dir, "magn_vs_temp_several_u0.png")) 
    plt.savefig(os.path.join(out_dir, "heat_caps_magn_vs_temp_several_u0_var.png"), bbox_inches='tight')

    # plot transition temp vs u0
    # x values on a log 10 scale
    u_values = []
    # y values
    trans_temps = []
    trans_temps_var = []
    
    # collect Tc vs u0 data
    for u, u_data in data.items():
        u_num = float(u)
        trans_temp_index = get_trans_temp(u_data["heat_caps"])
        trans_temp_index_var = get_trans_temp(u_data["heat_caps_var"])
        
        # transition found
        if trans_temp_index:
            trans_temp = u_data["temps"][trans_temp_index]
            trans_temp_var = u_data["temps"][trans_temp_index_var]
            u_values.append(u_num)
            trans_temps.append(trans_temp)
            trans_temps_var.append(trans_temp_var)
    
    # initialise u0 vs Tc plot
    plt.clf()
    plt.figure(figsize=(6,4), dpi=220)
    plt.xlabel("Potential parameter $u_0$")
    plt.xscale("log")
    plt.ylabel("Transition temperature $kT_c$")
    
    # copy u_values for sorting
    u_values_var = u_values.copy()

    # sort x,y values (heat caps by diff)
    vals = zip(u_values, trans_temps)
    res = sorted(vals, key = lambda x: x[0])

    # unzip after sorting
    u_values = [ i for i, j in res ]
    trans_temps = [ j for i, j in res ]
    
    # sort x,y values (heat caps by variance)
    vals = zip(u_values_var, trans_temps_var)
    res = sorted(vals, key = lambda x: x[0])

    # unzip after sorting
    u_values_var = [ i for i, j in res ]
    trans_temps_var = [ j for i, j in res ]
    
    # write sorted data to CSV output file (diff)
    with open(os.path.join(out_dir, "trans_temp_vs_u0.csv"), "w") as tc_out:
        # CSV header
        tc_out.write("u0,tc\n")
        # CSV data
        for i in range(len(u_values)):
            tc_out.write("{},{}\n".format(u_values[i], trans_temps[i]))

    # plot data and save
    #plt.plot(u_values, trans_temps, linewidth=1, marker="o", markersize="2", label="differentiation")
    plt.plot(u_values_var, trans_temps_var, linewidth=1, marker="o", markersize="3", color="#1a6ece", label="variance")
    plt.tick_params(axis="both", direction="in")
    #plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "trans_temp_vs_u0.png"), bbox_inches='tight')

    # write sorted data to CSV output file (var)
    with open(os.path.join(out_dir, "trans_temp_vs_u0_var.csv"), "w") as tc_out:
        # CSV header
        tc_out.write("u0,tc\n")
        # CSV data
        for i in range(len(u_values_var)):
            tc_out.write("{},{}\n".format(u_values_var[i], trans_temps_var[i]))

# heat capacity at temp

# find peak in heat capacity at critical temp Tc

# plot Tc vs u0


if __name__ == "__main__":
    main()
