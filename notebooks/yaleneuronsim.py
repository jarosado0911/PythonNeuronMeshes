# import some libraries, not all  of these I use, but there are here for reference
import neuron
import os
print(neuron.__version__)
from neuron import h
from neuron import rxd, gui, gui2 
from neuron.units import ms, mV, um
import matplotlib.pyplot as plt
import numpy as np

# if you get a warning then rerun this cell
gui2.set_backend('jupyter')

h.load_file("import3d.hoc")
h.load_file("stdrun.hoc")

class yaleneuronsim:
    def __init__(self,filename):
        self.load_morphology(filename)
        
    def load_morphology(self,filename):
        cell = h.Import3d_SWC_read()
        cell.input(filename)
        i3d = h.Import3d_GUI(cell, 0)
        i3d.instantiate(None)
        
    def set_discretization(self,dx):
        # add Hodgkin Huxley dynamics to each section
        for s in h.allsec():
            s.insert('hh')

            # use delta to specify the number of segments on the section
            n=int(s.L/dx)
            if n==0:
                n=1
            s.nseg=n

    def get_swc_locs(self,print_data=False):
        # initialize an empty dictionary
        nlocs={}

        # this forloop collects all the information
        for s in h.allsec():
            # here I set each section to have axial resistance 1000
            s.Ra=5000
            n3d=int(h.n3d(sec=s))
            if n3d==0:
                continue

            if print_data:
                # print some information
                print('\nSection:        ',s)
                print('Section Length: ',s.L)
                print('Number of discretization points (not the same as points in .swc): ', s.nseg)
                print('Mechanism: ', s.psection()['density_mechs'].keys())
                print('Axial Resistance: ',s.psection()['Ra'])

            # add items to the dictionary and print the coordinates
            for i in range(n3d):
                info={}
                info['sec']=s; info['L']=s.L;
                # these coordinates are the points in the .swc file NOT the discretization points!
                xx=h.x3d(i,sec=s); yy=h.y3d(i,sec=s); zz=h.z3d(i,sec=s); pos=tuple([xx,yy,zz]); ratio=s.arc3d(i)/s.L;
                info['ratio']=ratio
                if print_data:
                    print('%:',"{:.4f}".format(ratio),' (x,y,z): ',str(xx),',',str(yy),',',str(zz))
                nlocs[pos]=info
                
        return nlocs
    
    def setup_run(self,v_initial,end_time,dur,amp,delay):
        # set up the current clamp
        iclamp = h.IClamp(h.soma[0](0.5))
        iclamp.delay = delay * ms
        iclamp.amp = amp
        iclamp.dur = dur * ms

        # make a vector the time variable
        t=h.Vector().record(h._ref_t)

        # initialize empty data dictionary
        vdata={};

        # iterate through all sections
        for s in h.allsec():
            n3d=int(h.n3d(sec=s))
            if n3d==0:
                continue
            for i in range(n3d):
                xx=h.x3d(i,sec=s); yy=h.y3d(i,sec=s); zz=h.z3d(i,sec=s); 
                pos=tuple([xx,yy,zz]); ratio=s.arc3d(i)/s.L;
                # make the dictionary with recording vector using the ratio
                vdata[pos]=h.Vector().record(s(ratio)._ref_v)

        # initialize the simulation
        h.finitialize(v_initial * mV)

        # run the simulation
        h.continuerun(end_time * ms)
        
        return vdata,t