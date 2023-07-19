#!/usr/bin/env python
'''

Cell optimization with stress applied on z direction of Silicon, use SW potential
'''

from ase.optimize.fire import FIRE
from ase.optimize import BFGS
from ase.optimize import MDMin
from ase import *
from ase.io import read,write
import numpy as np
from ase.calculators.lammpsrun import LAMMPS
from ase import units
from ase.stressbox import stressbox
import math

parameters = { 'units' : 'metal', 
               'atom_style':'atomic',
               'boundary': 'p p p', 
               'pair_style': 'sw',
               'pair_coeff': ['* * Si.sw Si'] }

lmp_calc=LAMMPS(parameters=parameters)

refatom = read('si_alpha_1000atoms_sw_relax.lmpdata',format='lammps-data',style="atomic")

patom = read('si_alpha_1000atoms_random.lmpdata',format='lammps-data',style="atomic")

patom.set_calculator(lmp_calc)

pstress = patom.get_cell()*0.0

fixstrain = np.ones((3,3))

ref_atom = refatom.copy()
ref_coord = refatom.get_cell()

pstress[2][2] = 1.34523705e+01 

# PK1 can be replaced by PK2 or Cauchy
pbox =  stressbox(patom,express=pstress, fixstrain = fixstrain, ref_atom=ref_atom, stress_type='PK1')

dyn=MDMin(pbox)
#dyn=FIRE(pbox)
#dyn=BFGS(pbox)

f_max = 1E-6

dyn.run(fmax=f_max,steps=10000)

print('pk1 stress=',pbox.get_stress_pk1()/units.GPa)
print('pk2 stress=',pbox.get_stress_pk2()/units.GPa)
print('cauchy stress=',pbox.get_stress(voigt=False)/units.GPa)
print('deformation gradient=',pbox.get_defgrad())

write("output.lmpdata",patom, format='lammps-data',atom_style='atomic')
