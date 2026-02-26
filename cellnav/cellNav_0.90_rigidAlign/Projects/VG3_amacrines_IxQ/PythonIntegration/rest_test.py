from __future__ import division
from neuron import h
h.load_file('stdrun.hoc')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import string

import numbers
from neuron.units import ms, mV
from mpl_toolkits.mplot3d import Axes3D

# a helper library, included with NEURON
h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

#Reconstructing morphology

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
        sec.Ra = 150    # Axial resistance in Ohm * cm
        sec.cm = 1 #Capacitance in microFarad/cm2
        sec.g_pas = 0.0001  # Passive conductance in S/cm2          
        sec.e_pas = -49.2
        #print(sec.diam,d) 
        #sec.diam = d #look at diam being funny
    cell.all = cell.soma + cell.apic + cell.dend + cell.axon
    return cell,real_secs
cell1,secs_list = load_swc('cid4.swc')



#stim = h.IClamp(0.1,sec = secc)
#stim.delay = 5
#stim.dur = 1
#stim.amp = 0.1


sec_syn = secs_list['Import3d_Section[4460]']
syn = h.AlphaSynapse(0.1,sec = sec_syn) #0.1 is the synapse position along cylindrical section 1610. The cylinder length is scaled to the range 0 - 1. 
syn.e = 0 # reversal potential in mV
syn.gmax = 0.0 # max conductance in S/cm2
syn.onset = 1 # Time delay in ms 
syn.tau = 2

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


soma = secs_list['Import3d_Section[1]']
dend_sec_1040 = secs_list['Import3d_Section[1040]']
dend_1610_v = h.Vector().record(sec_syn(0.5)._ref_v) #0.5 is posn of recording electrode along cylindrical section
dend_1040_v = h.Vector().record(dend_sec_1040(0.5)._ref_v)
soma_v = h.Vector().record(soma(0.5)._ref_v)


volt_list = []

for w in range(0,len(secs_list)):
    dend_sec = secs_list['Import3d_Section[%s]' %str(w)]
    dend_v = h.Vector().record(dend_sec(0.5)._ref_v)
    volt_list.append(dend_v)

h.finitialize(-65 * mV)
h.continuerun(105 * ms)
print(max(volt_list[1]))






#print(max(dend_1610_v))
"""
plt.plot(t, soma_v,'-', color = 'blue',label='soma_v')
plt.plot(t, dend_1610_v,'-' , color = 'black', label='dend_with_synapse')
plt.plot(t, dend_1040_v,'-' , color = 'red', label='far_away_dend')
plt.xlabel('time (ms)')
plt.ylabel('Voltage (mV)')
plt.legend()
plt.show()

"""
"""
fig = plt.figure(figsize=(8,4))
ax1 = fig.add_subplot(2,1,1)
soma_plot = ax1.plot(t, soma_v,'-', color = 'blue')
dend1_plot = ax1.plot(t, dend_1610_v,'-' , color = 'black')
dend2_plot = ax1.plot(t, dend_1040_v,'-' , color = 'red')

rev_plot = ax1.plot([t[0], t[-1]], [syn.e, syn.e],
        color='green', linestyle=':')
ax1.legend(soma_plot + dend1_plot + dend2_plot + rev_plot,
        ['soma_v', 'dend_with_synapse(sec1610)', 'far_away_dend(sec1040)', 'syn_reversal_potential'])
ax1.set_ylabel('Voltage (mV)')
ax1.set_xticks([]) # Use ax2's tick labels

ax2 = fig.add_subplot(2,1,2)
syn_plot = ax2.plot(t, syn_i_vec, color='green')
ax2.legend(syn_plot, ['synaptic current'])
ax2.set_ylabel('Current (%s)' %str(h.units('AlphaSynapse.i')))
ax2.set_xlabel('time (ms)')
plt.show()




#h.PlotShape(False).plot(plt)


#secs_list['Import3d_Section[1]']
sl = h.SectionList([secs_list['Import3d_Section[%s]'%str(j)] for j in range(1,100)]) #len(secs_list)
ps = h.PlotShape(sl,False)

a = 1
for sec in sl:
    sec.v = max(volt_list[a])
    a = a+1

ps.scale(-65, 0)
ps.variable('v')
ax = ps.plot(plt, cmap=cm.jet)
ps.show(0)
#ps.len_scale()
#ps.show(0)
#ps.point_mark_remove()
#plt.show()


"""



swc = np.loadtxt('cid4.swc')
#print(swc)
v_max = []
x_sec = []
y_sec = []
z_sec = []
diam_d = []
a = 0
print((volt_list[2][-1]))
for a in range(0,len(secs_list)):
    sec = secs_list['Import3d_Section[%s]' %str(a)]
    #print(h.x3d(1,sec=sec))
    diam_d.append(2*sec.diam**2)
    x_sec.append(h.x3d(1,sec=sec))
    y_sec.append(h.y3d(1,sec=sec))
    z_sec.append(h.z3d(1,sec=sec))
    v_max.append((volt_list[a][-1]))
    a = a+1
x_pos = []
y_pos = []
z_pos = []
diam = []
parent = []
index = []
color = []
#print(len(swc[0]))
for a in range(0,len(swc)):
    index.append(swc[a][0])
    x_pos.append(swc[a][2])
    y_pos.append(swc[a][3])
    z_pos.append(swc[a][4])
    diam.append(100*swc[a][5]**2)
    parent.append(swc[a][6])
    #color.append([255 0 0 0.2])
#print(index[4460])


#plt.plot([x1,x2],[y1,y2],linewidth=mylw,color=plt.cm.rainbow(intensity))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_axis_off()
ax.scatter(x_pos[0], y_pos[0], z_pos[0], zdir='z', s=2000, c='black') #x_pos,  cmap=cm.jet
ax.scatter(x_pos, y_pos, z_pos, zdir='z', s=diam, c=(0.1, 0.1, 0.1, 0.4)) 
ax.scatter(x_sec[4460], y_sec[4460]-1, z_sec[4460], zdir='z', s=3, c='red' )
#ax.scatter(x_sec[4460], y_sec[4460]-1, z_sec[4460], zdir='z', s=3, c='black')

#ax.scatter(-15, -15, 15, zdir='z', s=3, c='black')
#ax.scatter(-13, -15, 15, zdir='z', s=3, c='black')
#ax.scatter(-15, -13, 15, zdir='z', s=3, c='black')
#ax.scatter(-15, -15, 17, zdir='z', s=3, c='black')


pl = ax.scatter(x_sec, y_sec, z_sec, zdir='z', s=0.5, c=v_max, cmap = cm.jet,vmin = -80, vmax = 20, depthshade=0) #cm.jet, depthshade=1

#cm.jet, depthshade=1 s = 0.2
fig.colorbar(pl, shrink=0.2, aspect=5,label = 'Voltage (in mV)')

#print(dir(fig))
#fig.set_clim(-80,20)
plt.show()













#plt.show()




#Writing a tree plot

















