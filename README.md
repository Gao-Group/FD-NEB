# Finite deformation nudged elastic band method


Solid-state nudged elastic band (SSNEB) methods can be used for finding solid-solid transition paths when solids are subjected to external stress fields. However, previous SSNEB methods may lead to inaccurate barriers and deviated reaction paths for transitions under finite/large deformation due to an inaccurate evaluation of the external work contributions in enthalpies. We report Finite deformation NEB method in 

A. Ghasemi, P. Xiao, W. Gao, Nudged elastic band method for solid-solid transition under finite deformation , The Journal of Chemical Physics, 151, 054110, 2019. [https://doi.org/10.1063/1.5113716](https://doi.org/10.1063/1.5113716)

FD-NEB is implemented based on open source project [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/index.html) and [Transition State Library for ASE (TSASE)](https://theory.cm.utexas.edu/tsase/).


## Usage



Once you have a local copy of ASE and this code, include the top directory to your $PYTHONPATH.
```
export PYTHONPATH=<path to FD-NEB code>/neb:$PYTHONPATH
export PYTHONPATH=<path to ASE code>/ase:$PYTHONPATH
```

Within the 'example' folder, two cases are presented. One employs LAMMPS as the calculator while the other utilizes VASP. Both cases examine the phase transition of silicon from the diamond to the &beta-tin phase.

When using LAMMPS as calculator, set up the enviroment variable to point to the LAMMPS executable:
```
export ASE_LAMMPSRUN_COMMAND="<path to lammps executable>/lmp_serial"
```
or for mpi version LAMMPS
```
export ASE_LAMMPSRUN_COMMAND="mpirun <path to lammps executable>/lmp_mpi"
```

When using VASP as calculator, write a script called run_vasp.py containing something like this:
```
import os
exitcode = os.system('mpirun <path to vasp executable>/vasp_std')
```

then setup environment variables to run_vasp.py and VASP potential directory:

```
export VASP_SCRIPT=<path to run_vasp.py>/run_vasp.py
export VASP_PP_PATH=<path to VASP potential directory>/VASP_potential
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

For solid state NEB calculations, it's necessary to pre-load the initial and final states to a target stress prior to executing the NEB. This can be done using the stressbox.py script.

The methodology behind stressbox approach is outlined in the following publication:

Ghasemi, A., Gao, W., "A method to apply Piola-Kirchhoff stress in molecular statics simulations", [DOI: 10.1016/j.commatsci.2021.110496](https://doi.org/10.1016/j.commatsci.2021.110496).

We've included an illustrative script, opt_example.py, to demonstrate the application of either Piola-Kirchhoff or Cauchy stress.


## Stress used in FD-NEB

For large or finite deformations, obtaining a more accurate enthalpy barrier requires conducting Nudged Elastic Band (NEB) calculations with Piola-Kirchhoff (PK) stresses, instead of Cauchy stress.  

In FD-NEB code, user can choose 'cauchy', 'PK1' or 'PK2'. For small deformaiton, they should yield very similar resutls. For finite deformaiton, 'PK1' or 'PK2' should be used. There is anohter option in stres type, which is called 'hydro', meaning hydrostatic pressure. 

Since PK stresses may not be familiar to nonexperts in mechanics, we brifely explain the difference in stress definition:  

Cauchy stress is most commonly used in atomistic simulations, because it measures the force per unit area in the deformed configuration and can be directly computed using Viral stress formula. On the other hand, PK stresses (including the first and second kind) have been widely used in solid mechanics for finite deformation problems. Briefly, the first PK stress tensor is also called engineering stress or nominal stress, because it measures the force per unit area in reference configuration. The second PK stress tensor (S) is defined entirely in the reference configuration: using a fictitious force pulled from the deformed configuration, which is then divided by the corresponding area in the reference configuration. One of the advantages of PK stresses is that they have well defined work conjugates, allowing accurate evaluation of the work done by a constant external stress. In fact, in experiments, it is the applied force (or hydrostatic presssure) that is easily controlled not the Cauchy stress due to the difficulty of tracking the deformed area. Therefore, it is difficult to have Cauchy stress (except hydrostatic pressure) stays constant during the transition process.


