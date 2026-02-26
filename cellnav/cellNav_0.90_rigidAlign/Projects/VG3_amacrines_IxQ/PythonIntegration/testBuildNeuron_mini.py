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
#from neuron import gui
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from matplotlib.pyplot import cm

h.load_file('stdrun.hoc')




'''
class Cell:
    def __init__(self,name='neuron',soma=None,apic=None,dend=None,axon=None):
        self.soma = soma if soma is not None else []
        self.apic = apic if apic is not None else []
        self.dend = dend if dend is not None else []
        self.axon = axon if axon is not None else []

        self.all = self.soma + self.apic + self.dend + self.axon

    def delete(self):
        self.soma = None
        self.apic = None
        self.dend = None
        self.axon = None
        self.all = None

    def __str__(self):
        return self.name

'''


class Cell:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self.all = self.seed.wholetree()
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        h.define_shape()
        
        
    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)
    
    def _set_position(self, x, y, z):
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(i,
                               x - self.x + sec.x3d(i),
                               y - self.y + sec.y3d(i),
                               z - self.z + sec.z3d(i),
                              sec.diam3d(i))
        self.x, self.y, self.z = x, y, z
        
   
class skel(Cell):
    name = 'skel'
    
    def _setup_morphology(self):
        self.seed = h.Section(name='seed', cell=self)


    def _setup_biophysics(self):
        pass
        for sec in self.all:
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        self.seed.insert('hh')                                          
        for seg in self.seed:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        self.seed.insert('pas')                 
        for seg in self.seed:
            seg.pas.g = 0.001  # Passive conductance in S/cm2
            seg.pas.e = -65    # Leak reversal potential mV


cell = skel(0)


pos = [0 0 0]


## Create first section by hand
sec_list = []
sec = cell.seed
buffersize =  h.pt3dclear(sec=sec)
print(sec)
sec_list.append(sec)   
h.pt3dconst(1, sec=sec)
style = h.pt3dstyle(0, sec=sec)
h.pt3dadd(pos[0,1], pos[0,0], pos[0,2] , rad[0][0]*2, sec=sec) 




## add non-seed sections
for secC in range(1,sNum):
    sec = h.Section(cell=cell) #just creating a random neuron section
    #buffersize =  h.pt3dclear(sec=sec)
    print(sec)
    sec_list.append(sec)  
    style = h.pt3dstyle(0, sec=sec)
    x, y, z = pos[secC,:]
    d = rad[secC][0] * 2
    h.pt3dadd( x, y, z, d, sec=sec) 
    
for secC in range(0,sNum):
    sec.insert('pas')
    sec.Ra = 100000000#ex['Ra'][0] # 100    # Axial resistance in Ohm * cm
    sec.cm = 100000000#ex['cm'][0] #1 #Capacitance in microFarad/cm2
    sec.g_pas = ex['g_pas'][0] #0.001  # Passive conductance in S/cm2          
    sec.e_pas = ex['e_pas'][0] #-65

    #'''
    sec.insert('hh')                                          
    for seg in sec:
        seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
        seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
        seg.hh.gl = 1#0.0003    # Leak conductance in S/cm2
        seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        sec.insert('pas')                 
    for seg in sec:
        seg.pas.g = 1#0.001  # Passive conductance in S/cm2
        seg.pas.e = -65    # Leak reversal potential mV
    #'''
 
    #cell.all = cell.soma + cell.dend + cell.axon

## create connections
for secC in range(1,sNum):
    sec = sec_list[secC]    
    sec.connect(sec_list[parent[secC][0]]) 
    h.pt3dconst(1, sec=sec)

h.define_shape()


## make volt List
volt_list = []
for secC in range(0,sNum):   
    sec = sec_list[secC]
    dend_v = h.Vector().record(sec(0)._ref_v)
    volt_list.append(dend_v)

## make syn volt list
synVolt_list = []
for secS in range(0,len(synEdges)):
    sec = sec_list[secS]
    syn_v =  h.Vector().record(sec(0)._ref_v)
    synVolt_list.append(syn_v)


v_max = []
x_sec = []
y_sec = []
z_sec = []
diam_d = []
for secC in range(0,sNum):
    sec = sec_list[secC]
    #print(h.x3d(1,sec=sec))
    diam_d.append(2*sec.diam**2)
    x_sec.append(h.x3d(0,sec=sec))
    y_sec.append(h.y3d(0,sec=sec))
    z_sec.append(h.z3d(0,sec=sec))
    v_max.append(1)


if 0:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_axis_off()
    
    pl = ax.scatter(x_sec, y_sec, z_sec, zdir='z', s=0.5, c=v_max, cmap = cm.jet,vmin = -80, vmax = 20, depthshade=0) #cm.jet, depthshade=1
    plt.show()


## run model for each synapse
maxVolt = np.zeros((len(synEdges),len(sec_list)))
synSecPos = np.zeros((len(synEdges),3))

for syn_num in range(0,len(synEdges)):
    
    print('running {} of {} synapses'.format(syn_num+1,len(synEdges)))

    pickEdge = int(synEdges[syn_num])
    sec_syn = sec_list[pickEdge]
    t = h.Vector().record(h._ref_t)

    
    if 0:
        stim = h.NetStim()
        syn_ = h.ExpSyn(sec_syn(.5))
        stim.number = 1
        stim.start = 0
        stim.stop = 1
        ncstim = h.NetCon(stim, syn_)
        ncstim.delay = 1 * ms
        ncstim.weight[0] = 0.04
        syn_.tau = 2 * ms
        
    else:
        syn = h.AlphaSynapse(0,sec = sec_syn) #0.1 is the synapse position along cylindrical section 1610. The cylinder length is scaled to the range 0 - 1.
        #syn = h.ExpSyn(0,sec = sec_syn)
        syn.e = ex['e'][0] # reversal potential in mV
        syn.gmax = -.1#ex['gmax'][0] # max conductance in S/cm2
        syn.onset = 0 # Time delay in ms 
        syn.tau = .1 #ex['tau'][0]
        syn.i = 1   
        
        #Recording the simulation/We'll start out recording the membrane potential at the center of the soma and the time in two NEURON Vectors:
        syn_i_vec = h.Vector()
        syn_i_vec.record(syn._ref_i)
    
    h.finitialize(-65 * mV)
    h.continuerun(20 * ms)
    del(syn)
    
    if 1:
        plt.plot(t, volt_list[pickEdge-1], label='soma(0.5)')
        plt.plot(t, volt_list[0], label='soma(0.5)')
        plt.plot(t,syn_i_vec)
        plt.plot(t, synVolt_list[syn_num])
        plt.show()
    
    v_max = []
    for a in range(0,len(sec_list)):
        sec = sec_list[a]
        v_max.append(min(volt_list[a]))   
        
    if 0:
        plt.plot(v_max)
        plt.show()

    maxVolt[syn_num, :] = v_max
    
    
    x_sec.append(h.x3d(0,sec=sec))
    y_sec.append(h.y3d(0,sec=sec))
    z_sec.append(h.z3d(0,sec=sec))
    synSecPos[syn_num,0] = sec_syn.y3d(0)
    synSecPos[syn_num,1] = sec_syn.x3d(0)
    synSecPos[syn_num,2] = sec_syn.z3d(0)
    
## Save data for matlab       
secPos = np.zeros((len(y_sec),3))
secPos[:,0] = y_sec  
secPos[:,1] = x_sec
secPos[:,2] = z_sec


nrnFile = '{}nrn_cid{}.mat'.format(WPNswc,cid)
from scipy.io import savemat
mdic = {'maxVolt': maxVolt, 'secPos':secPos,'synSecPos':synSecPos}
savemat(nrnFile,mdic)




