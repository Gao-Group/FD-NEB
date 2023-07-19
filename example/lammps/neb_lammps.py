#!/usr/bin/env python

'''
neb optimization example: silicon phase transition under stress, use LAMMPS as calculator
'''

from neb import fdneb, fire_ssneb, qm_ssneb
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read
import os
import sys
import numpy

# read initial and final images
p1 = read('0.CON',format='vasp')
p2 = read('6.CON',format='vasp')

# define reference image, usually use intial image
# reference image is only needed when using PK1 or PK2 stress, no need to define it when using Cauchy stress
ref_img = read('ref.CON',format='vasp')

root_dir=os.getcwd()

# define lammps input script keywords
parameters = { 'units' : 'metal', 'atom_style':'charge', 'pair_style': 'sw',
               'boundary': 'p p p', 'mass': ['1 28.0855'],
              'pair_coeff': ['* *  '+root_dir+'/Si.sw  Si']}

calc=LAMMPS(parameters=parameters)

# use this calculator to save all lammps files for debuging 
'''
calc = LAMMPS(parameters=parameters, tmp_dir="lmp_dir")
'''

p1.set_calculator(calc)
p2.set_calculator(calc)

# set fixstrain to zero to fix the simulation cell 
# e.g. fixstrain[2] *= 0.0 fix the cell deformation along xz,yz,zz directions
fixstrain = numpy.ones((3,3))

# define the target stress
# e.g., pstress[0][0] = -2.0 is to set sigma_xx = -2 GPa, negative means compression 
pstress = p1.get_cell()*0.0
pstress[0][0] = -2.0
pstress[1][1] = -2.0
pstress[2][2] = 4.0

# number of images along MEP
nim = 7



# p1.......... one endpoint of the path
# p2.......... the other endpoint of the path
# numImages... the total number of images in the path, including the endpoints
# method...... "ci" for the climbing image method, anything else for normal NEB method 
# ss.......... boolean, solid-state neb or regular neb 
# express..... external press, 3*3 lower triangular matrix in the unit of GPa
# fixstrain... 3*3 matrix as express. 0 fixes strain at the corresponding direction
# stress_type. 'cauchy' or 'PK1'(first Piola–Kirchhoff) or 'PK2' (second Piola–Kirchhoff) stress
# ref_img..... refence image for computing PK stress, not needed for cauchy stress      
band = fdneb(p1, p2, numImages=nim, method='ci', ss=True, fixstrain=fixstrain, express=pstress, stress_type='PK1', ref_img=ref_img)

         



# this commented code can be used to read a pre-existing band in order to continue the NEB calcualtion 
'''
for i in range(1,nim-1):
    filename = str(i)+'.CON'
    b = read(filename,format='vasp')
    band.path[i].set_positions(b.get_positions())
    band.path[i].set_cell(b.get_cell())
'''    

# define optimizer and parameters, use FIRE or quick_min for minimizaiton
opt = fire_ssneb(band, maxmove =0.1, dtmax = 0.1, dt=0.1)
'''
opt = qm_ssneb(band, maxmove = 0.100, dt = 0.08)
'''

# run the minimization, define convergence criteria and maximum interation step
opt.minimize(forceConverged=0.01, maxIterations = 250)


# create a movie to show MEP
os.system('touch movie.con')
for p in range(len(band.path)):
    os.system('cat '+str(p)+'.CON >>movie.con')
