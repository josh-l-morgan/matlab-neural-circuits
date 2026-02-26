# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 16:26:22 2021

@author: jlmorgan
"""

'''
import sys
sys.path.append("c:\\nrn\\lib\\python")
'''

if 1:
    import neuron
    from neuron import h
    from neuron.units import ms, mV
    #from neuron import gui
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.io
    from matplotlib.pyplot import cm    
    h.load_file('stdrun.hoc')
    import random



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
        
    def __repr__(self):
        return '{}[{}]'.format(self.name, self._gid)
   
class skel(Cell):
    name = 'skel'
    
   

useDefaultShape = 0;
cell = skel(0)


WPNswc = 'Y:\\Active\morganLab\\karlsRetina\\CellNavLibrary_IxQ\Volumes\\AprilMerge\\Analysis\\swc\\'
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
synEdges = sm2nrn['synCloseNode'][0][0]-1;
nep = sm2nrn['nep'][0][0]
swc = nep['swcS'][0][0]
parent = swc['pred'][0][0]
#pos = swc['pos'][0][0]
pos = nep['pos'][0][0]
edges = nep['edges'][0][0]-1 #parent is [:,1]?
rad = nep['nodeRad'][0][0]
edgeRad = nep['edgeRad'][0][0]
sNum = len(rad)
eNum = edges.shape[0]
# L = np.sqrt(np.subtract(pos[:,0].reshape(sNum,1),pos[parent,0])**2 + \
#          np.subtract(pos[:,1].reshape(sNum,1),pos[parent,1])**2 + \
#              np.subtract(pos[:,2].reshape(sNum,1),pos[parent,2])**2 )
    
    
L = np.sqrt(np.subtract(pos[edges[:,0],0],pos[edges[:,1],0])**2 + \
        np.subtract(pos[edges[:,0],1],pos[edges[:,1],1])**2 + \
            np.subtract(pos[edges[:,0],2],pos[edges[:,1],2])**2 )
        

if 0:
    sNum = sNum
    #pos = pos[0:sNum,:]
    #parent = parent[0:sNum]
    #rad = rad[0:sNum] # * 0 + 20
    synEdges = [sNum-100]
    #L = L[0:sNum] #* 0 + 20
    #h.topology()


'''
##Setup temp properties for sections
sNum = 10;
pos = np.zeros((sNum,3))
posSize = pos.shape 
for r in range(posSize[0]):
    pos[r,:] = np.array([1,2,3]) + r * 10
print(pos) 
rad = np.ones((sNum))
parent = np.zeros((sNum,1),int) 
parent = np.array(range(0,sNum))
'''


sec_list = []

'''
## Create first section by hand
sec = cell.seed
buffersize =  h.pt3dclear(sec=sec)
print(sec)
sec_list.append(sec)   
h.pt3dconst(1, sec=sec)
style = h.pt3dstyle(0, sec=sec)
h.pt3dadd(pos[0,1], pos[0,0], pos[0,2] , rad[0][0]*2, sec=sec) 
sec.diam = rad[0][0]*2
'''

## create sections
for secC in range(0,sNum):
    sec = h.Section(cell=cell) #just creating a random neuron section
    #buffersize =  h.pt3dclear(sec=sec)
    sec_list.append(sec)  
    sec.nseg = 1
    style = h.pt3dstyle(0, sec=sec)
    #x, y, z = pos[secC,:]
    d = rad[secC][0] * 2
    #h.pt3dadd( x, y, z, d, sec=sec) 
    sec.diam = d;
    sec.L = 1#L[secC];

## set up biophysics    
for secC in range(0,sNum):
    sec.insert('pas')
    sec.Ra = ex['Ra'][0] # 100    # Axial resistance in Ohm * cm
    sec.cm = ex['cm'][0] #1 #Capacitance in microFarad/cm2
    sec.g_pas = ex['g_pas'][0] #0.001  # Passive conductance in S/cm2          
    sec.e_pas = ex['e_pas'][0] #-65

    '''
    sec.insert('hh')                                          
    for seg in sec:
        seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
        seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
        seg.hh.gl = 0.0003    # Leak conductance in S/cm2
        seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite
        sec.insert('pas')                 
    for seg in sec:
        seg.pas.g = 1#0.001  # Passive conductance in S/cm2
        seg.pas.e = -65    # Leak reversal potential mV
    #'''
 
    #cell.all = cell.soma + cell.dend + cell.axon

## create connections
for secC in range(0,eNum):
    sec = sec_list[edges[secC,0]]   
    #sec.connect(sec_list[random.randint(0,secC-1)],0,1) 
    parent = sec_list[edges[secC,1]]
    
    sec.connect(parent,0,1) 
    #sec.connect(sec_list[parent[secC][0]]) 
    #h.pt3dconst(1, sec=sec)
    #print(sec.parentseg())
    
if useDefaultShape:
    h.define_shape()


## make volt List
volt_list = []
for secC in range(0,sNum):   
    sec = sec_list[secC]
    dend_v = h.Vector().record(sec(1)._ref_v)
    volt_list.append(dend_v)

## make syn volt list
synVolt_list = []
for secS in range(0,len(synEdges)):
    sec = sec_list[secS]
    syn_v =  h.Vector().record(sec(0)._ref_v)
    synVolt_list.append(syn_v)


'''
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


if 1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_axis_off()
    
    pl = ax.scatter(x_sec, y_sec, z_sec, zdir='z', s=0.5, c=v_max, cmap = cm.jet,vmin = -80, vmax = 20, depthshade=0) #cm.jet, depthshade=1
    plt.show()
'''

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
        syn = h.ExpSyn(sec_syn(1))
        stim.number = 1
        stim.start = 0
        ncstim = h.NetCon(stim, syn)
        ncstim.delay = 10 * ms
        ncstim.weight[0] = 1
        syn.tau = 1 * ms
        
    else:
        syn = h.AlphaSynapse(.5,sec = sec_syn) #0.1 is the synapse position along cylindrical section 1610. The cylinder length is scaled to the range 0 - 1.
        #syn = h.ExpSyn(0,sec = sec_syn)
        syn.e = ex['e'][0] # reversal potential in mV
        syn.gmax = 1#ex['gmax'][0] # max conductance in S/cm2
        syn.onset = 1 # Time delay in ms 
        syn.tau = 1 * ms#ex['tau'][0]
        syn.i = 1  
        
    #Recording the simulation/We'll start out recording the membrane potential at the center of the soma and the time in two NEURON Vectors:
    syn_i_vec = h.Vector()
    syn_i_vec.record(syn._ref_i)
    
    h.finitialize(-65 * mV)
    h.continuerun(100 * ms)
    del(syn)
    #del(stim)
    
    if 1:
        #plt.plot(t, volt_list[0], label='soma(0.5)')
        #plt.plot(t, volt_list[-1], label='soma(0.5)')
        plt.plot(t, volt_list[pickEdge], label='soma(0.5)')

        #plt.plot(t,syn_i_vec) 
        plt.plot(t, synVolt_list[syn_num])
        plt.show()
    
    v_max = []
    for a in range(0,len(sec_list)):
        sec = sec_list[a]
        v_max.append(max(volt_list[a]))   
        
    if 0:
        plt.plot(v_max)
        plt.show()

    maxVolt[syn_num, :] = v_max
    
    ''' 
    x_sec.append(h.x3d(0,sec=sec))
    y_sec.append(h.y3d(0,sec=sec))
    z_sec.append(h.z3d(0,sec=sec))
    synSecPos[syn_num,0] = sec_syn.y3d(0)
    synSecPos[syn_num,1] = sec_syn.x3d(0)
    synSecPos[syn_num,2] = sec_syn.z3d(0)
    '''
    
## Save data for matlab    
secPos = []   
#secPos = np.zeros((len(y_sec),3))
#secPos[:,0] = y_sec  
#secPos[:,1] = x_sec
#secPos[:,2] = z_sec


nrnFile = '{}nrn_cid{}.mat'.format(WPNswc,cid)
from scipy.io import savemat
mdic = {'maxVolt': maxVolt, 'secPos':secPos,'synSecPos':synSecPos}
savemat(nrnFile,mdic)




