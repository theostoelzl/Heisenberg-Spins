#!/bin/bash
#
# Heat-bath Monte Carlo of Landau-Ising spinds
#
# Each task will process one of the 6 one-line text input files,
# passing it into our (not very magic) rev program, which simply
# writes the line backwards. The result is stored in an output
# file named in the same way as the inputs.
#
#SBATCH --job-name=hl_metro_sampling_ultralight
#SBATCH --partition=long
#SBATCH --time=7-00:00:00
#SBATCH --array=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 
# --array=0,1,2,3,4,5,6,7,8,9,10,11
#
########################################################################

echo "*** Starting job "$SLURM_JOB_ID $SLURM_ARRAY_TASK_ID $SLURM_JOB_NAME " ***"
echo "Start date: " $(date)
echo "Start time: " $(date +%s)
echo "Nodes:      " $SLURMD_NODENAME
echo "CPUs:       " $SLURM_NTASKS

# Decide which input file this Slurm task will process.
# We use the special $SLURM_ARRAY_TASK_ID variable, which tells the
# task which one it is among the whole job array.
#u0_vals=(1.0 1.5 2.0 2.4 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0)
u0_vals=(20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000)
u0=${u0_vals[${SLURM_ARRAY_TASK_ID}]}

# Create new directory for u0 value
mkdir hl_u0_$u0

# Copy necessary files to u0 dir
cp HL_run_temps.py hl_u0_$u0
cp HL_MC_metro.o hl_u0_$u0
cp system.txt hl_u0_$u0

# Change to new directory
cd hl_u0_$u0

# Call Python script to run all temps
# because Bash can't do floats and I don't
# want to waste my life energy getting the
# bc cmd to work
python HL_run_temps.py $u0 > output_$u0.out

# Change back to parent dir
cd ..

echo "End  date: " $(date)
echo "End  time: " $(date +%s)
echo "*** Finished job "$SLURM_JOB_ID $SLURM_ARRAY_TASK_ID $SLURM_JOB_NAME " ***"
