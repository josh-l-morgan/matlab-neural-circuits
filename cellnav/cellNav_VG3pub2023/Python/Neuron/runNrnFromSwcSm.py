# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 13:51:36 2021

@author: jlmorgan
"""

from __future__ import division
import sys
sys.path.append("c:\\nrn8\\lib\\python")
#import neuron
from neuron import h
#from neuron import gui
import time

h.load_file('stdlib.hoc')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import string
#import matlab.engine
#eng = matlab.engine.start_matlab()
#import mat73
import scipy.io
h.load_file('import3d.hoc')
h.load_file('stdrun.hoc')
import numbers
from neuron.units import ms, mV
from mpl_toolkits.mplot3d import Axes3D

#Reconstructing morphology
cid = 3323;
WPNswc = 'E:\\IxQ_KarlsRetinaVG3_2019\\CellNavLibrary_IxQ\Volumes\\AprilMerge\\Analysis\\swc\\'
sm2nrnFile = '{}sm2nrn_cid{}.mat'.format(WPNswc,cid)
swcFile = '{}cid{}.swc'.format(WPNswc,cid)
readVar = scipy.io.loadmat(sm2nrnFile)
sm2nrn = readVar['sm2nrn']
synEdges = sm2nrn['synCloseNode'][0][0];


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

# load the data. Use Import3d_SWC_read for swc, Import3d_Neurolucida3 for
# Neurolucida V3, Import3d_MorphML for MorphML (level 1 of NeuroML), or
# Import3d_Eutectic_read for Eutectic.
def load_swc(filename, fileformat=None, cell=None, use_axon=True, xshift=0, yshift=0, zshift=0):
    """
    Load an SWC from filename and instantiate inside cell. Code kindly provided
    by @ramcdougal.

    Args:
        filename = .swc file containing morphology
        cell = Cell() object. (Default: None, creates new object)
        filename = the filename of the SWC file
        use_axon = include the axon? Default: True (yes)
        xshift, yshift, zshift = use to position the cell

    Returns:
        Cell() object with populated soma, axon, dend, & apic fields

    Minimal example:
        # pull the morphology for the demo from NeuroMorpho.Org
        from PyNeuronToolbox import neuromorphoorg, load
        with open('c91662.swc', 'w') as f:
            f.write(neuromorphoorg.morphology('c91662'))
        cell = load(filename)

    """
    #    https://www.neuron.yale.edu/neuron/static/py_doc/modelspec/programmatic/topology/geometry.html

    if cell is None:
        cell = Cell(name="".join(filename.split('.')[:-1]))

    if fileformat is None:
        fileformat = filename.split('.')[-1]

    name_form = {1: 'soma[%d]', 2: 'axon[%d]', 3: 'dend[%d]', 4: 'apic[%d]'}

    # load the data. Use Import3d_SWC_read for swc, Import3d_Neurolucida3 for
    # Neurolucida V3, Import3d_MorphML for MorphML (level 1 of NeuroML), or
    # Import3d_Eutectic_read for Eutectic.
    if fileformat == 'swc':
        morph = h.Import3d_SWC_read()
    elif fileformat == 'asc':
        morph = h.Import3d_Neurolucida3()
    else:
        raise Exception('file format `%s` not recognized'%(fileformat))
    morph.input(filename)

    # easiest to instantiate by passing the loaded morphology to the Import3d_GUI
    # tool; with a second argument of 0, it won't display the GUI, but it will allow
    # use of the GUI's features
    i3d = h.Import3d_GUI(morph, 0)

    # get a list of the swc section objects
    swc_secs = i3d.swc.sections
    swc_secs = [swc_secs.object(i) for i in range(int(swc_secs.count()))] #1681 sections are put in swc_secs. Each section is called Import3d_Section[i]
    # initialize the lists of sections
    sec_list = {1: cell.soma, 2: cell.axon, 3: cell.dend, 4: cell.apic}
    # name and create the sections
    real_secs = {}
    for swc_sec in swc_secs:
        cell_part = int(swc_sec.type) #axon,soma or dendrite...
        # skip everything else if it's an axon and we're not supposed to
        # use it... or if is_subsidiary
        if (not(use_axon) and cell_part == 2) or swc_sec.is_subsidiary:
            continue
        # figure out the name of the new section
        if cell_part not in name_form:
            raise Exception('unsupported point type')
        name = name_form[cell_part] % len(sec_list[cell_part])
        # create the section
        sec = h.Section(cell=cell) #just creating a random neuron section
        # connect to parent, if any
        if swc_sec.parentsec is not None: #swc_sec.parentsec finds parent of swc_sec. Then connects to parent
            sec.connect(real_secs[swc_sec.parentsec.hname()](swc_sec.parentx))
        #real_secs has list of section number and corresponding Neuron section ID. Like ..'Import3d_Section[215]': __nrnsec_0x55e69b406450,.. 
        # define shape
        if swc_sec.first == 1: #swc_sec.first == 1 only for first section. Maybe its the parent section so you define shape first? Or is it a non-spherical seciton? 
            h.pt3dstyle(1, swc_sec.raw.getval(0, 0), swc_sec.raw.getval(1, 0),
                        swc_sec.raw.getval(2, 0), sec=sec)
        j = swc_sec.first
        xx, yy, zz = [swc_sec.raw.getrow(i).c(j) for i in range(3)] # [Vector[22683], Vector[22685], Vector[22687]]
        dd = swc_sec.d.c(j)
        if swc_sec.iscontour_:
            # never happens in SWC files, but can happen in other formats supported
            # by NEURON's Import3D GUI
            raise Exception('Unsupported section style: contour')
        if dd.size() == 1: #xx,yy,zz,dd is a vector for location
            # single point soma; treat as sphere
            x, y, z, d = [dim.x[0] for dim in [xx, yy, zz, dd]]
            for xprime in [x - d / 2., x, x + d / 2.]:
                h.pt3dadd(xprime + xshift, y + yshift, z + zshift, d, sec=sec) 
        else:
            for x, y, z, d in zip(xx, yy, zz, dd):
                h.pt3dadd(x + xshift, y + yshift, z + zshift, d, sec=sec)

        # store the section in the appropriate list in the cell and lookup table               
        sec_list[cell_part].append(sec)    
        real_secs[swc_sec.hname()] = sec  
        #this_seg = sec
        sec.insert('pas')
        sec.Ra = 100    # Axial resistance in Ohm * cm
        sec.cm = 1 #Capacitance in microFarad/cm2
        sec.g_pas = 0.001  # Passive conductance in S/cm2          
        sec.e_pas = -65
        #print(sec.diam,d) 
        #sec.diam = d #look at diam being funny
    cell.all = cell.soma + cell.apic + cell.dend + cell.axon
    return cell,real_secs


cell1,secs_list = load_swc(swcFile)
#cell1,secs_list = load_swc('cid4.swc')

# a helper library, included with NEURON


maxVolt = np.zeros((len(synEdges),len(secs_list)))

for syn_num in range(0,len(synEdges)):
    
    print('running {} of {} synapses'.format(syn_num,len(synEdges)))

    #sec_syn = secs_list['Import3d_Section[%s]' %str(2000)]
    #secStr = 'Import3d_Section[{}]'.format(1000);
    secStr = 'Import3d_Section{}'.format(synEdges[syn_num])
    sec_syn = secs_list[secStr]
    #sec_syn = secs_list['Import3d_Section[%s]' %str(synEdges[syn_num])]
    syn = h.AlphaSynapse(0.5,sec = sec_syn) #0.1 is the synapse position along cylindrical section 1610. The cylinder length is scaled to the range 0 - 1.
    syn.e = 0 # reversal potential in mV
    syn.gmax = 0.8 # max conductance in S/cm2
    syn.onset = 0 # Time delay in ms 
    syn.tau = 0.1
    
    """
    sec_syn1 = secs_list['Import3d_Section[4460]']
    syn1 = h.AlphaSynapse(0.1,sec = sec_syn1) #0.1 is the synapse position along cylindrical section 1610. The cylinder length is scaled to the range 0 - 1. 
    syn1.e = 0 # reversal potential in mV
    syn1.gmax = 0.8 # max conductance in S/cm2
    syn1.onset = 5 # Time delay in ms 
    syn1.tau = 0.1 # Decay time constant in ms
    print(h.psection(sec_syn))
    print(syn.get_segment())
    print("asyn.tau = {}".format(syn.tau))
    """
    
    
    #Recording the simulation/We'll start out recording the membrane potential at the center of the soma and the time in two NEURON Vectors:
    t = h.Vector().record(h._ref_t)
    
    syn_i_vec = h.Vector()
    syn_i_vec.record(syn._ref_i)
    """
    
    soma = secs_list['Import3d_Section[1]']
    dend_sec_1040 = secs_list['Import3d_Section[1040]']
    dend_1610_v = h.Vector().record(sec_syn(0.5)._ref_v) #0.5 is posn of recording electrode along cylindrical section
    dend_1040_v = h.Vector().record(dend_sec_1040(0.5)._ref_v)
    soma_v = h.Vector().record(soma(0.5)._ref_v)
    """
    
    volt_list = []
    
    for w in range(0,len(secs_list)):
        dend_sec = secs_list['Import3d_Section[%s]' %str(w)]
        dend_v = h.Vector().record(dend_sec(0.5)._ref_v)
        volt_list.append(dend_v)
    
    h.finitialize(-65 * mV)
    h.continuerun(2 * ms)
    
    
    v_max = [] 
    for a in range(0,len(secs_list)):
        sec = secs_list['Import3d_Section[%s]' %str(a)]
        #print(h.x3d(1,sec=sec))
       # diam_d.append(2*sec.diam**2)
        #x_sec.append(h.x3d(1,sec=sec))
       # y_sec.append(h.y3d(1,sec=sec))
       # z_sec.append(h.z3d(1,sec=sec))
        v_max.append(max(volt_list[a]))   
        
    maxVolt[syn_num,:] = v_max;
            
## Save

#nrnFile = '{}cid{}.npy'.format(WPNswc,cid)
#np.save(nrnFile,maxVolt)

nrnFile = '{}nrn_cid{}.mat'.format(WPNswc,cid)
from scipy.io import savemat
mdic = {'maxVolt': maxVolt, 'label':'experiment'}
savemat(nrnFile,mdic)
