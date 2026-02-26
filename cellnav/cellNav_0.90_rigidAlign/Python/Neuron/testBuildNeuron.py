# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 16:26:22 2021

@author: jlmorgan
"""

'''
import sys
sys.path.append("c:\\nrn\\lib\\python")
'''


import neuron
from neuron import h
from neuron.units import ms, mV


class Cell:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self.all = self.soma.wholetree()
        self._setup_biophysics()
    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)



class BallAndStick(Cell):
    name = 'BallAndStick'
    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.dend.connect(self.soma)
        self.soma.L = self.soma.diam = 12.6157
        self.dend.L = 200
        self.dend.diam = 1
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        self.soma.insert('hh')                                          
        for seg in self.soma:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        self.dend.insert('pas')                 
        for seg in self.dend:
            seg.pas.g = 0.001  # Passive conductance in S/cm2
            seg.pas.e = -65    # Leak reversal potential mV

my_cell = BallAndStick(0)



for sec in h.allsec():
    print('%s: %s' % (sec, ', '.join(sec.psection()['density_mechs'].keys())))

cell = BallAndStick('hi')
otherCell = BallAndStick('bye')
h.topology()

del otherCell

import matplotlib.pyplot as plt

h.PlotShape(False).plot(plt)





soma = h.Section(name= 'soma')
h.topology()
soma.psection()
soma.psection()['morphology']['L']
soma.L = 20;
soma.diam = 20;
soma.insert('hh') # Figure out best channels for amacrine cell
iclamp = h.IClamp(soma(0.5))
print([item for item in dir(iclamp) if not item.startswith('__')])
iclamp.delay = 2
iclamp.dur = 0.1
iclamp.amp = 0.9
soma.psection()
v = h.Vector().record(soma(0.5)._ref_v)
t = h.Vector().record(h._ref_t)

h.load_file('stdrun.hoc')
h.finitialize(-65 * mV)

h.continuerun(40 * ms)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(t,v)
plt.xlabel('t (ms)')
plt.ylabel('v (mV)')
plt.show()












'''
def load_json(morphfile):
    with open(morphfile, 'r') as f:
        secdata = json.load(morphfile)
        
    seclist = []
    for sd in secdata:
        # make section
        sec = h.Section(name = sd['name'])
        seclist.append(sec)
        
        #make 3d morphology
        for x,y,z,d in zip(sd['x'], sd['y'], sd['z'], sd('diam')):
            h.pt3dadd(x, y, z, d, sec = sec)
            
    # connect children to parent compartments
    for sec,sd in zip(seclist,secdata):
            if sd['parent_loc'] >= 0:
                parent_sec = sec_list[sd['parent']]
                sec.connect(parent_sec(sd['parent_loc']), sd['section_ori'])
                
    return seclsection

'''










