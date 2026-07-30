#!/bin/bash
#SBATCH --job-name=mag_col_lammps
#SBATCH --array=0-39
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=4-00:00:00
#SBATCH --error=%x-%A_%a.err # This sets the output files with the name of the job + the ID (if arrays are employed, used "_%a" for the number of the element)
#SBATCH --output=%x-%A_%a.out

# Load necessary modules (Fortran compiler, CMake for the Makefile, OpenMPI for parallelizing...)
module load GCC/13.2.0
module load CMake/3.27.6-GCCcore-13.2.0
module load OpenMPI/4.1.6-GCC-13.2.0
module load FFTW/3.3.10-GCC-13.2.0

# Math managing
MY_LIST=($(seq 0.25 0.25 5.0))
MY_LIST2=(6 40)
#MY_LIST2=(0.5 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5)

LEN1=${#MY_LIST[@]}
task_id=${SLURM_ARRAY_TASK_ID}

index1=$((task_id % $LEN1))
index2=$((task_id / $LEN1))

MY_F=$(echo "${MY_LIST2[$index2]}" | tr ',' '.') #The variables in the list are somehow output with commas, I have to swap them to dots
#MY_MU=$(echo "${MY_LIST2[$index2]}" | tr ',' '.')

#MY_F=2
MY_MU=1.8
MY_VAR=$(echo "${MY_LIST[$index1]}" | tr ',' '.')

# Calculates the value of the corresponding coefficient and changes the output file names accordingly
GAMMA=$(python3 -c "print(f'{2*($MY_MU**2):.2f}'.replace('.', 'p'))")
LE=$(python3 -c "print(f'{$MY_F / (2*($MY_MU)**2):.2f}'.replace('.', 'p'))") 
PE=$(python3 -c "print(f'{$MY_F:.2f}'.replace('.', 'p'))") 
TAUR=$(python3 -c "print(f'{$MY_VAR:.2f}'.replace('.', 'p'))")

VAR_P_LABEL=gamma_${GAMMA}_Pe_${PE}_tauR_${TAUR} 

# Setup working directories
TMPDIR=/tmp/${SLURM_JOB_ID}
mkdir -p ${TMPDIR}
cp -pv ${SLURM_SUBMIT_DIR}/initial.lmp ${TMPDIR}
cd ${TMPDIR}
# MNT_DIR="/mnt/storage/lzhu/LAMMPS_magcol" # This is just for saving output, you could just make an output folder in your directory.

OUTDIR=${SLURM_SUBMIT_DIR}/output/
mkdir -p $OUTDIR

# Run LAMMPS calling the binary file compiled with the modified version 
srun -n4 --mpi=pmi2 /home/lzhu/lammps-22Jul2025/build/lmp_mpi -v myforce $MY_F -v mymoment $MY_MU -v myvar $MY_VAR -v my_str $VAR_P_LABEL -in $TMPDIR/initial.lmp

# This is for copying from the TMPDIR to the MNT or storage dir, if you made an output folder as mentioned above, skip this step and do the commented part before this.
# scp -i "${HOME}/.ssh/computacio" *.lammps* superfe:"${MNT_DIR}/"

# Copy results to output
cp -v $TMPDIR/*.lammpstrj $OUTDIR/
cp -v $TMPDIR/*.data $OUTDIR/

module purge

