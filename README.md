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
