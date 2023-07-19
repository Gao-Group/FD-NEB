'''
Nudged elastic band quick-min optimizer module
'''

from .util import vproj,vdot,vmag,vunit
from .minimizer_ssneb import minimizer_ssneb
import numpy as np

class qm_ssneb(minimizer_ssneb):
    '''
    Class containing a quick-min optimizer for nudged elastic bands
    '''

    def __init__(self, path, maxmove=0.2, dt=0.05):
        '''
        Constructor
            path    - neb object to optimize
            pot     - potential energy surface to optimize on
            maxmove - maximum distance that the optimizer can move in one step
            dt      - differential timestep
        '''
        minimizer_ssneb.__init__(self, path)
        self.maxmove = maxmove
        self.dt=dt

        #self.v contains the velocity for both atoms and box for one image
        i = self.band.numImages-2
        j = self.band.natom+3
        self.v = np.zeros((i,j,3))


    def step(self):
        '''
        Take a step
        '''
        self.band.forces()
        for i in range(1, len(self.band.path) - 1):
            totalf = self.band.path[i].totalf
            Power = vdot(totalf,self.v[i-1])
            if Power > 0.0 :
                self.v[i-1]  = vproj(self.v[i-1], totalf)
            else:
                self.v[i-1] *= 0.0
                #self.dt *= 0.999
            # Euler step
            self.v[i-1] += self.dt * totalf

            # check for max step
            if vmag(self.v[i-1]) > self.maxmove/self.dt :
                self.v[i-1] = self.maxmove/self.dt * vunit(self.v[i-1])
            dR = self.dt * self.v[i-1]
            # move R
            rt  = self.band.path[i].get_positions()
            rt += dR[:-3]
            self.band.path[i].set_positions(rt)
            # move box and update cartesian coordinates
            ct  = self.band.path[i].get_cell()
            ct += np.dot(ct, dR[-3:]) / self.band.jacobian
            self.band.path[i].set_cell(ct, scale_atoms=True)


