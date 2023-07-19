# Finite deformation nudged elastic band method


The method is reported in the following paper:
A. Ghasemi, P. Xiao, W. Gao, Nudged elastic band method for solid-solid transition under finite deformation , The Journal of Chemical Physics, 151, 054110, 2019. https://doi.org/10.1063/1.5113716

FD-NEB is implemented based on the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/index.html). The FD-NEB computation code is developed based on the G-SSNEB code. It is implemented based on an open source project [Transition State Library for ASE (TSASE)](https://theory.cm.utexas.edu/tsase/).

Once you have a local copy of ASE and this code, include the top directory to your $PYTHONPATH.
```
export PYTHONPATH=<path to code>/neb:$PYTHONPATH
export PYTHONPATH=<path to code>/ase:$PYTHONPATH
```

Within the 'example' folder, two cases are presented. One employs LAMMPS as the calculator while the other utilizes VASP. Both cases examine the phase transition of silicon from the diamond to the &beta-tin phase.

When using LAMMPS as calculator, set up the enviroment variable to point to the LAMMPS executable:
```
export ASE_LAMMPSRUN_COMMAND="<path to lammps executable>/lmp_serial"
```
or 
```
export ASE_LAMMPSRUN_COMMAND="mpirun <path to lammps executable>/lmp_mpi"
```

When using VASP as calculator, write a script called run_vasp.py containing something like this:
```
import os
exitcode = os.system('mpirun <path to vasp executable>/vasp_std')
```

then setup environment variables to run_vasp.py and vasp potential directory:

```
export VASP_SCRIPT=<path to run_vasp.py>/run_vasp.py
export VASP_PP_PATH=<path to vasp potential directory>/vasp_pot
```

To run the examples:

```
cd ./example/lammps
python neb_lammps.py
```
or 
```
cd ./example/vasp
python neb_vasp.py
```

For NEB calculations involving solid states, it's necessary to load the initial and final states to a target stress prior to executing the NEB. This can be done using the stressbox.py script.

The methodology behind this approach is outlined in the following publication:

Ghasemi, A., Gao, W., "A method to apply Piola-Kirchhoff stress in molecular statics simulations", DOI: 10.1016/j.commatsci.2021.110496.

We've included an illustrative script, opt_example.py, to demonstrate the application of either Piola-Kirchhoff or Cauchy stress.