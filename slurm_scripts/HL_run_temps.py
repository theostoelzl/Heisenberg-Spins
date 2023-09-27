"""

Takes in u0 as argv and runs MC programme for Heisenberg-Landau
for given temperature range. (Temperature is actually kT!)

"""

import sys
import subprocess

# potential characteristic
u0 = sys.argv[1]

if float(u0) > 1.5:
    # high temperature to start from
    highT = 1.35
    # go to low temperature
    lowT = 1.25
else:
    highT = 1.57
    lowT = 1.47

# how many temps should be sampled
samplesT = 40

# find temperature steps
dT = (highT-lowT)/samplesT

# start loop over temps, starting from high T
for i in range(samplesT):
    # get temp
    T = highT - dT*i

    # start MC
    if i == 0:
        out = subprocess.run([ "./HL_MC_metro.o", "none", str(T), "1084", str(u0) ])
    else:
        out = subprocess.run([ "./HL_MC_metro.o", "spins_after.txt", str(T), "1084", str(u0) ])

    print(out, flush=True)

print("Done.", flush=True)
