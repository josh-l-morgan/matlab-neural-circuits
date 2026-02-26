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
        
        # everything below here in this method is NEW
        self._spike_detector = h.NetCon(self.seed(0.5)._ref_v, None, sec=self.seed)
        self.spike_times = h.Vector()
        self._spike_detector.record(self.spike_times)
        
        self._ncs = []
        
        self.soma_v = h.Vector().record(self.seed(0.5)._ref_v)

        
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
synEdges = sm2nrn['synCloseNode'][0][0];
nep = sm2nrn['nep'][0][0]
swc = nep['swcS'][0][0]
parent = swc['pred'][0][0]
pos = swc['pos'][0][0]
rad = swc['rad'][0][0]
sNum = len(rad)




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
for secC in range(sNum):
    sec = h.Section(cell=cell) #just creating a random neuron section
    #buffersize =  h.pt3dclear(sec=sec)
    print(sec)
    sec_list.append(sec)  
    style = h.pt3dstyle(0, sec=sec)
    x, y, z = pos[secC,:]
    d = rad[secC][0] * 2
    h.pt3dadd( x, y, z, d, sec=sec) 

 
    #cell.all = cell.soma + cell.dend + cell.axon



## create connections
for secC in range(1,sNum):
    sec = sec_list[secC]    
    sec.connect(sec_list[parent[secC][0]]) 
    h.pt3dconst(1, sec=sec)

h.define_shape()


v_max = []
x_sec = []
y_sec = []
z_sec = []
diam_d = []
for secC in range(1,sNum):
    sec = sec_list[secC]
    #print(h.x3d(1,sec=sec))
    diam_d.append(2*sec.diam**2)
    x_sec.append(h.x3d(0,sec=sec))
    y_sec.append(h.y3d(0,sec=sec))
    z_sec.append(h.z3d(0,sec=sec))
    v_max.append(1)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_axis_off()

pl = ax.scatter(x_sec, y_sec, z_sec, zdir='z', s=0.5, c=v_max, cmap = cm.jet,vmin = -80, vmax = 20, depthshade=0) #cm.jet, depthshade=1
plt.show()


'''

#name = name_form[cell_part] % len(sec_list[cell_part])
# create the section

# connect to parent, if any


#real_secs has list of section number and corresponding Neuron section ID. Like ..'Import3d_Section[215]': __nrnsec_0x55e69b406450,.. 
# define shape
#if swc_sec.first == 1: #swc_sec.first == 1 only for first section. Maybe its the parent section so you define shape first? Or is it a non-spherical seciton? 
#    h.pt3dstyle(1, swc_sec.raw.getval(0, 0), swc_sec.raw.getval(1, 0),
#               swc_sec.raw.getval(2, 0), sec=sec)
    
j = swc_sec.first
# [Vector[22683], Vector[22685], Vector[22687]]








# store the section in the appropriate list in the cell and lookup table               

#print(swc_sec.hname())
#print(len(sec_list[3]))
real_secs[0] = sec  
#this_seg = sec
sec.insert('pas')
sec.Ra = ex['Ra'][0] # 100    # Axial resistance in Ohm * cm
sec.cm = ex['cm'][0] #1 #Capacitance in microFarad/cm2
sec.g_pas = ex['g_pas'][0] #0.001  # Passive conductance in S/cm2          
sec.e_pas = ex['e_pas'][0] #-65
#print(sec.diam,d) 
#sec.diam = d #look at diam being funny
return cell,real_secs


'''










'''
for secC in range(sNum):
    cell_part = 'dendrite' #axon,soma or dendrite...
    #print('start section')
    #print(swc_sec.id)
    c = c+1
    #print(c)
    #time.sleep(1)
    
    
    if cell_part not in name_form:
        raise Exception('unsupported point type')
    name = name_form[cell_part] % len(sec_list[cell_part])
    # create the section
    sec = h.Section(cell=cell) #just creating a random neuron section
    # connect to parent, if any
    if swc_sec.parentsec is not None: #swc_sec.parentsec finds parent of swc_sec. Then connects to parent
        sec.connect(real_secs[swc_sec.parentsec.hname()](swc_sec.parentx))
   
    sec.connect(startSec) 
   
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
    #print(swc_sec.hname())
    #print(len(sec_list[3]))
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
    
'''





## build sec
'''
sec = new Import3d_Section(first, i-first)
sec.append(0, first, i-first, x, y, z, d)

  pt2sec(pix2ix(first), sec.parentsec)
  psec = sec.parentsec

 sec.append(0, pix2ix(first), 1, x, y, z, d)
 sec.append(1, first, i-first, x, y, z, d)

 sec.type = type.x[first]
 sections.append(sec)
  sec.pid = sec.parentsec.id
 
 
 for i=0, id.size-1 {
f (point2sec.x[i] > isec ) { // in next section
ptadd(pix2ix(i), point2sec.x[i])

sec = point2sec.x[i]
tadd(i, isec)
 '''
 
 
 
 
 
 
 
 
 
 
'''
## Ball and stick tutorial 2

#mycell = BallAndStick(0, 0, 0, 0, 0)

my_cells = create_n_BallAndStick(7, 50)

#h.PlotShape(False).plot(plt)

ps = h.PlotShape(True)
ps.show(0)

my_cells = create_n_BallAndStick(5, 50)

h.topology()



##    Make synapses

stim = h.NetStim() # Make a new stimulator

# Attach it to a synapse in the middle of the dendrite
# of the first cell in the network. (Named 'syn_' to avoid
# being overwritten with the 'syn' var assigned later.)
syn_ = h.ExpSyn(my_cells[0].dend(0.5))

stim.number = 1
stim.start = 9
ncstim = h.NetCon(stim, syn_)
ncstim.delay = 1 * ms
ncstim.weight[0] = 0.04 # NetCon weight is a vector.
syn_.tau = 2 * ms


recording_cell = my_cells[1]
soma_v = h.Vector().record(recording_cell.soma(0.5)._ref_v)
dend_v = h.Vector().record(recording_cell.dend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)
syn_i = h.Vector().record(syn_._ref_i)
spike_times = [h.Vector() for nc in netcons]
for nc, spike_times_vec in zip(netcons, spike_times):

h.finitialize(-65 * mV)
h.continuerun(25 * ms)


import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(2, 1, 1)
soma_plot = ax1.plot(t, soma_v, color='black', label='soma(0.5)')
dend_plot = ax1.plot(t, dend_v, color='red', label='dend(0.5)')
rev_plot = ax1.plot([t[0], t[-1]], [syn_.e, syn_.e], label='syn reversal',
        color='blue', linestyle=':')
ax1.legend()
ax1.set_ylabel('mV')
ax1.set_xticks([]) # Use ax2's tick labels

ax2 = fig.add_subplot(2, 1, 2)
syn_plot = ax2.plot(t, syn_i, color='blue', label='synaptic current')
ax2.legend()
ax2.set_ylabel(h.units('ExpSyn.i'))
ax2.set_xlabel('time (ms)')
plt.show()


syns = []
netcons = []
for source, target in zip(my_cells, my_cells[1:] + [my_cells[0]]):
    syn = h.ExpSyn(target.dend(0.5))
    nc = h.NetCon(source.soma(0.5)._ref_v, syn, sec=source.soma)
    nc.weight[0] = 0.05
    nc.delay = 5
    netcons.append(nc)
    syns.append(syn)

h.finitialize(-65 * mV)
h.continuerun(100 * ms)
plt.plot(t, soma_v, label='soma(0.5)')
plt.plot(t, dend_v, label='dend(0.5)')
plt.legend()
plt.show()

for i, spike_times_vec in enumerate(spike_times):
    print('cell {}: {}'.format(i, list(spike_times_vec)))


spike_times = [h.Vector() for nc in netcons]
for nc, spike_times_vec in zip(netcons, spike_times):
    nc.record(spike_times_vec)

h.finitialize(-65 * mV)
h.continuerun(100 * ms)

for i, spike_times_vec in enumerate(spike_times):
    print('cell {}: {}'.format(i, list(spike_times_vec)))

'''



'''
## Ball and stick 1
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










