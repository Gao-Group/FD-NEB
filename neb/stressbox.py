#!/usr/bin/env python
'''
The method implemented in this script is reported in the following reference:
Arman Ghasemi, Wei Gao. "A method to apply Piola-Kirchhoff stress in molecular statics simulations", submitted.

The script is modified based on the TSASE package developed in Henkelman Group at UT Austin. 

To use this script, please copy it to the directoy of ASE (tested on versiton 3.19.0): ase-3.19.0/ase/
'''

from ase import *
from ase.io import read,write
from ase import units
import numpy as np
import numpy
import copy 
from numpy.linalg import inv

class stressbox(Atoms):
    def __init__(self, atomsx, express=np.zeros((3,3)), fixstrain=np.ones((3,3)), ref_atom=[], stress_type= 'cauchy'):
        """box relaxation
        atomsx: an Atoms object
        fixstrain: 3*3 matrix as express. 0 fixes strain at the corresponding direction
        """
        self.atomsx = atomsx 
        self.stress_type = stress_type
        
        self.express= express *units.GPa # convert to ASE units: eV/A^3

        self.fixstrain = fixstrain
        
        self.ref_atom = ref_atom

        cell       = atomsx.get_cell()
        vol        = atomsx.get_volume()
        self.natom = atomsx.get_global_number_of_atoms()

        Atoms.__init__(self,atomsx)

    def get_positions(self):
        r    = self.atomsx.get_positions()*0.0
        Rc   = np.vstack((r, self.atomsx.get_cell()*0.0))
        return Rc

    def set_positions(self,dr):
        rcell  = self.atomsx.get_cell()        
        rcell += dr[-3:]*((self.natom)**(1.0/3.0))
        self.atomsx.set_cell(rcell)
        ratom  = self.atomsx.get_positions() + dr[:-3]
        self.atomsx.set_positions(ratom)

    def __len__(self):
        return self.natom+3

    def get_forces(self,apply_constraint=True):
        f    = self.atomsx.get_forces(apply_constraint)
        stt  = self.atomsx.get_stress()
        vol  = self.atomsx.get_volume()
        st   = np.zeros((3,3))
        
        if self.stress_type == 'PK1' or self.stress_type == 'PK2':                
            loc_cord = self.atomsx.get_cell()
            ref_cord = self.ref_atom.get_cell()
   
            deform_grad = np.dot(np.transpose(loc_cord),inv(np.transpose(ref_cord)))
            J_def_grad = np.linalg.det(deform_grad)

        #following the order of get_stress in vasp.py
        # (the order of stress in ase are the same for all calculators)
        st[0][0] = -stt[0]   
        st[1][1] = -stt[1] 
        st[2][2] = -stt[2] 
        st[2][1] = -stt[3] 
        st[2][0] = -stt[4] 
        st[1][0] = -stt[5] 

        if self.stress_type == 'PK1':
            express_cauchy = 1/(J_def_grad)*numpy.dot(self.express,numpy.transpose(deform_grad))

        if self.stress_type == 'PK2':
            express_cauchy = 1/(J_def_grad)*numpy.dot(numpy.dot(deform_grad,self.express),numpy.transpose(deform_grad))

        if self.stress_type == 'cauchy':
            express_cauchy = np.copy(self.express) # do a hard copy, so self.express does not change along with express_cauchy

        # upper triangular of stress matrix is not used, set them to zero
        express_cauchy[1][2] = 0
        express_cauchy[0][2] = 0
        express_cauchy[0][1] = 0

        st  -= express_cauchy
        # apply constrain
        st *= self.fixstrain

        Fc = np.vstack((f, st*((vol/(self.natom))**(2.0/3.0))))
        
        return Fc
    
    def get_potential_energy(self,force_consistent=True):
        return self.atomsx.get_potential_energy(force_consistent)

    def get_stress(self,voigt=True):
        return self.atomsx.get_stress(voigt)

    def get_defgrad(self):
        loc_cord = self.atomsx.get_cell()
        ref_cord = self.ref_atom.get_cell()
        deform_grad = np.dot(np.transpose(loc_cord),inv(np.transpose(ref_cord)))
        return deform_grad
         
    def get_stress_pk1(self):
        loc_cord = self.atomsx.get_cell()
        ref_cord = self.ref_atom.get_cell()
        deform_grad = np.dot(np.transpose(loc_cord),inv(np.transpose(ref_cord)))
        J_def_grad = np.linalg.det(deform_grad)
        return  J_def_grad*np.dot(self.atomsx.get_stress(voigt=False),np.transpose(inv(deform_grad)))

    def get_stress_pk2(self):
        loc_cord = self.atomsx.get_cell()
        ref_cord = self.ref_atom.get_cell()
        deform_grad = np.dot(np.transpose(loc_cord),inv(np.transpose(ref_cord)))
        J_def_grad = np.linalg.det(deform_grad)
        return  J_def_grad*np.dot(np.dot(inv(deform_grad),self.atomsx.get_stress(voigt=False)),np.transpose(inv(deform_grad)))
        
    def copy(self):
        """Return a copy."""
        import copy
        atomsy = self.atomsx.copy()
        atoms = self.__class__(atomsy, self.express)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()
        return atoms

