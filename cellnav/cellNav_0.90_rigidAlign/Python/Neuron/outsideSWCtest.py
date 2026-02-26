# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:36:57 2021

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



## Get experiment information from matlab
WPNswc = 'E:\\IxQ_KarlsRetinaVG3_2019\\CellNavLibrary_IxQ\Volumes\\AprilMerge\\Analysis\\swc\\'
exFile = '{}experiment.mat'.format(WPNswc)
exVar = scipy.io.loadmat(exFile)
ex = {'cids':0,'tau':0,'e':0,'gmax':0}
ex['cids'] = exVar['ex']['cids'][0][0][0]
ex['tau'] = exVar['ex']['tau'][0][0][0]
ex['e'] = exVar['ex']['e'][0][0][0]
ex['gmax'] = exVar['ex']['gmax'][0][0][0]
ex['Ra'] = exVar['ex']['Ra'][0][0][0]
ex['cm'] = exVar['ex']['cm'][0][0][0]
ex['g_pas'] = exVar['ex']['g_pas'][0][0][0]
ex['e_pas'] = exVar['ex']['e_pas'][0][0][0]



## Get cell  from matlab
cid = ex['cids'][0];
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
    c = 0;
    for swc_sec in swc_secs:
        cell_part = int(swc_sec.type) #axon,soma or dendrite...
        print('start section')
        print(swc_sec.id)
        c = c+1;
        print(c)
        #time.sleep(1)
        if swc_sec.id > 470:
            print('hi')
            
        
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
        print(swc_sec.hname())
        print(len(sec_list[3]))
        real_secs[swc_sec.hname()] = sec  
        #this_seg = sec
        sec.insert('pas')
        sec.Ra = ex['Ra'][0] # 100    # Axial resistance in Ohm * cm
        sec.cm = ex['cm'][0] #1 #Capacitance in microFarad/cm2
        sec.g_pas = ex['g_pas'][0] #0.001  # Passive conductance in S/cm2          
        sec.e_pas = ex['e_pas'][0] #-65
        #print(sec.diam,d) 
        #sec.diam = d #look at diam being funny
    cell.all = cell.soma + cell.apic + cell.dend + cell.axon
    return cell,real_secs


swcFile = 'C:\\Users\\jlmorgan\\Documents\\MATLAB\\Neuron\\neuron_code\\5_josh\\cid4.swc'

cell1,secs_list = load_swc(swcFile,'swc')
#cell1,secs_list = load_swc('cid4.swc')
secs_list.values()
# a helper library, included with NEURON

volt_list = []
key_list = []
oldID = []
#for w in range(0,len(secs_list)):
for key in secs_list:   
    #id = int(key.split('[]')[0])
    dend_sec = secs_list[key]
    #dend_sec = secs_list['Import3d_Section[%s]' %str(w)]
    dend_v = h.Vector().record(dend_sec(0)._ref_v)
    volt_list.append(dend_v)
    key_list.append(key)
    oldID.append(dend_sec.id)

maxVolt = np.zeros((len(synEdges),len(secs_list)))


print('found {} segments. Maximimum edge number is {}'.format(len(secs_list),max(synEdges)))


for syn_num in range(0,len(synEdges)):

    
    print('running {} of {} synapses'.format(syn_num,len(synEdges)))

    
    secStr = 'Import3d_Section{}'.format(synEdges[syn_num])
    secStr = 'Import3d_Section{}'.format([1]) #!!!!!!!!!!!!!!!
    pickEdge = int(synEdges[syn_num])
    if pickEdge>len(key_list):
        print('synEdge index {} excedes the length ({}) of sec_list'.format(pickEdge,len(key_list)))
        pickEdge = len(key_list)
        
    secStr = key_list[pickEdge-1]

    sec_syn = secs_list[secStr]
    syn = h.AlphaSynapse(0.5,sec = sec_syn) #0.1 is the synapse position along cylindrical section 1610. The cylinder length is scaled to the range 0 - 1.
    syn.e = ex['e'][0] # reversal potential in mV
    syn.gmax = ex['gmax'][0] # max conductance in S/cm2
    syn.onset = 0 # Time delay in ms 
    syn.tau = ex['tau'][0]
    
    
    #Recording the simulation/We'll start out recording the membrane potential at the center of the soma and the time in two NEURON Vectors:
    t = h.Vector().record(h._ref_t)
    
    syn_i_vec = h.Vector()
    syn_i_vec.record(syn._ref_i)

    
    h.finitialize(-65 * mV)
    h.continuerun(2 * ms)
    
    
   
        
    v_max = []

    for a in range(0,len(secs_list)):
        sec = secs_list[key_list[a]]
        v_max.append(max(volt_list[a]))   
    
    ##record section properties
    if syn_num == 0:
        x_sec = []
        y_sec = []
        z_sec = []
        diam_d = []
        for a in range(0,len(secs_list)):
            sec = secs_list[key_list[a]]   
            diam_d.append(2*sec.diam**2)
            x_sec.append(h.x3d(1, sec=sec))
            y_sec.append(h.y3d(1, sec=sec))
            z_sec.append(h.z3d(1, sec=sec))
            
            
            

    maxVolt[syn_num, :] = v_max
    
       
pos = np.zeros((len(y_sec),3))
pos[:,0] = y_sec  
pos[:,1] = x_sec
pos[:,2] = z_sec
## Save

#nrnFile = '{}cid{}.npy'.format(WPNswc,cid)
#np.save(nrnFile,maxVolt)

nrnFile = '{}nrn_cid{}.mat'.format(WPNswc,cid)
from scipy.io import savemat
mdic = {'maxVolt': maxVolt, 'pos':pos}
savemat(nrnFile,mdic)











