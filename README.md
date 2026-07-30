NOTE: The LAMMPS script in this repository were designed to run on a **custom, modified build of LAMMPS** developed by the soft matter group at Universitat de Barcelona.

Specifically, this build includes customized fixes to **pair_lj_cut_dipole_cut**

If you are a collaborator or reviewer wishing to reproduce these exact trajectories, please contact carles.calero@ub.edu for access to the modified LAMMPS

The Python clustering analysis script (Count_Chains5.0.py) is completely standalone and process standard LAMMPS dump files (`.lammpstrj`). They can be tested using the sample file (dump_gamma_6p48_Pe_2p00_tauR_0p33.lammpstrj) provided.

### Count_Chains5.0.py
The code is designed to analyze the aggregation of active paramagnetic colloids. 
The code requires an input file generated via the dump function in LAMMPS as input and needs to be in the same directory as the code itself.
The input file must be named as dump_gamma_{mygammastr}_Pe_{myPestr}_tauR_{mytauRstr}.lammpstrj
Where gamma is the magnetic coupling parameter calculated as 2*m**2/T, where m is the magnetic moment of the particle and T the temperature
Pe is the active Peclet number equal to the self propelling force defined in "fix myactive" within the LAMMPS script
tauR is the rotational persistence time equalt to the gamma_r defined in "fix mybrownian" within the LAMMPS script
Note that the above are only valid when LJ reduced units are applied. For a detailed derivation of these parameters, please refer to the annex document

Within the code, the user is required to define the list of values of gamma, Pe and tauR spanned by the data files of interest
Note that, by default, the code iterates through all possible combinations of the three parameters. The user may adapt the iteration procedure to its needs

Once the list of parameters are defined, the user may choose the desired graphs to be output by the code
In plot_flags, the user determines the graph of interest that will be output for each specific case
In aggregate_flags, the user determines the graph of interest that will be output for the whole parameter space considered

Finally, the user may execute the code

### script_run.sh
The code initialises the parameter sweep in Slurm.
The parameters that can be swept are the magnetic dipole (MY_MU), the self-propelling force (MY_F) and a custom variable (MY_VAR). This last one is set to be the rotational drag coefficient gamma_r by default. This can be changed in the .lmp
The code sweeps through all possible combinations of elements in MYLIST and MYLIST2, which can be assigned to any of the previously mentioned parameters

### initial.lmp
The code initialises the 2D system in LAMMPS
The system is set to be in a 2D squared box of size 200x200 with periodic boundary conditions applied to both spatial coordinates containing 3000 particles
The code produces a .lammpstrj file containing the evolution of the system
