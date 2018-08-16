# Import modules
from __future__ import division
import numpy as np 
import matplotlib as mpl
from brian2 import *
from brian2.core.variables import *
import yaml
import matplotlib.pyplot as plt 
import time
import networkx as nx 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


class netSimObj:    
    eqs_neuron = '''
        dv/dt = (1/C)*(k*(v-v_r)*(v-v_t)/mV - u + I_app + I_syn) : volt
        du/dt = int(1-is_fs)*a*(b*(v-v_r) - u) + int(is_fs)*a*(big_U-u): amp
        big_U = int(v>=v_b)*b*((v-v_b)**3)/(mV*mV) : amp
        
        I_app = i0*(t>=ti)*(t<=tf) : amp
        I_syn = I_elec + I_AMPA + I_GABA : amp 

        I_elec : amp
        I_AMPA : amp 
        I_GABA : amp 
        C : farad
        v_r : volt
        v_t : volt
        v_p : volt
        k : amp/volt
        a : 1/second
        b : amp/volt
        c : volt
        d : amp
        is_fs : 1
        v_b : volt
        e_vpu : volt/amp
        e_cu : volt/amp
        e_du : amp
        i0 : amp 
        ti : second
        tf : second
        
    '''
    eqs_AMPA = '''
        g_AMPA : siemens 
        I_AMPA_post = g_AMPA*(E_AMPA-v)*s_AMPA : amp (summed)
        ds_AMPA/dt = -s_AMPA/tau_AMPA : 1 (clock-driven) 
    '''
    eqs_pre_glut = '''
        s_AMPA += 1
    '''
    eqs_GABA = '''
        g_GABA : siemens 
        I_GABA_post = g_GABA*(E_GABA-v)*s_GABA : amp (summed)
        ds_GABA/dt = -s_GABA/tau_GABA : 1  (clock-driven) 
    '''
    eqs_pre_gaba = '''
        s_GABA += 1
    '''
    eqs_elec = '''
        g_elec:siemens
        I_elec_post = g_elec*(v_pre-v_post) : amp (summed)
    '''
    def __init__(self, prm_file, net_file, net_id): 
        
        netgraph = yaml.load(open(net_file))[net_id]
        param = yaml.load(open(prm_file))
        self.N = shape(netgraph['nodes'])[0]
        self.cell_names = np.asarray(netgraph['nodes'])[:,0]
        self.cell_types = np.asarray(netgraph['nodes'])[:,1]
        self.connections = np.asarray(netgraph['edges'])
        
        cell_dict = []
        for it in self.cell_types: 
            cell_dict.append(param[it])
        unit_dict = param['units']
        init_vals = {} 
        for k in cell_dict[0].iterkeys():
            init_vals[k] = [_cd[k] for _cd in cell_dict]  * eval(unit_dict[k])
        self.init_states = init_vals 
        self.init_states['i0'] = np.zeros(self.N) * pA
        self.init_states['ti'] = np.zeros(self.N) * ms
        self.init_states['tf'] = np.zeros(self.N) * ms
        self.data = {'t': [], 'v' : []}
    def drawNetwork(self): 
        loc_conn = self.connections != '0'
        edge_labels = dict(zip(zip(*np.where(loc_conn)), self.connections[loc_conn]))
        node_labels = dict(zip(xrange(self.N), self.cell_names))
        graph = nx.from_numpy_matrix(loc_conn, create_using=nx.DiGraph())
        pos = nx.kamada_kawai_layout(graph,dist=(0.05,0.05))
        nx.draw(graph,pos,node_size=2000,linewidths=0,width=1.5)
        nx.draw_networkx_labels(graph,pos,labels=node_labels)
        nx.draw_networkx_edge_labels(graph,pos,edge_labels=edge_labels,label_pos=0.4,font_size=10)
        show()
    def setGlobalVars(self, globalVars, constVars):
        self.glob = globalVars 
        self.const = constVars
    def setSynapseVars(self, synVars, changeDefaultSynCond): 
        self.syn = {'default': synVars, 'custom': changeDefaultSynCond} 
    def setDCStimVars(self, dcstimVars): 
        for cell in dcstimVars.keys(): 
            loc_ = self.indexOfCell(cell) 
            for att in dcstimVars[cell].keys():
                self.init_states[att][loc_] = dcstimVars[cell][att]
    def returnData(self, cell_name2find):
        return self.data['v'][self.indexOfCell(cell_name2find)]
    def indexOfCell(self, cell_name2find):
        return np.where(np.array(self.cell_names) == cell_name2find)[0][0]
    def runSim(self, reportElapsedTime):
        G = NeuronGroup(N=self.N, model=self.eqs_neuron, threshold='v >= v_p + e_vpu*u',\
                reset='v = c + e_cu*u; u = (u+d)*(u+d<e_du) + e_du*(u+d>=e_du)', method='euler')
        S_elec = Synapses(G, G, model=self.eqs_elec)
        S_AMPA = Synapses(G, G, model=self.eqs_AMPA, on_pre=self.eqs_pre_glut, method='euler') 
        S_GABA = Synapses(G, G, model=self.eqs_GABA, on_pre=self.eqs_pre_gaba, method='euler') 
        
        for key in self.const.keys(): 
            val = float(self.const[key])
            dim = self.const[key].dim
            G.variables.add_constant(key, val, dim)
            S_elec.variables.add_constant(key, val, dim)
            S_AMPA.variables.add_constant(key, val, dim)
            S_GABA.variables.add_constant(key, val, dim)
            
        for i in xrange(self.N): 
            for j in xrange(self.N): 
                conn_ij_raw = self.connections[i,j] 
                if conn_ij_raw != '0': 
                    for conn_ij in conn_ij_raw.split('-'): 
                        if conn_ij.upper() == 'ELEC': # SYMMETRICAL ELECTRICAL SYNAPSE 
                            S_elec.connect(i=i,j=j)
                            S_elec.connect(i=j,j=i)
                        elif conn_ij.upper() == 'AMPA':
                            S_AMPA.connect(i=i,j=j)
                        elif conn_ij.upper() == 'GABA':
                            S_GABA.connect(i=i,j=j)
                        else:
                            print "pair (%d, %d) - value = %s " %(i,j,conn_ij)
                            raise ValueError('Cannot accept any other values representing '\
                                             'synapses except "0", "ELEC", "AMPA" or "GABA"')
        
        if size(S_elec.g_elec) > 0:
            S_elec.g_elec[:] = self.syn['default']['g_elec']
        if size(S_AMPA.g_AMPA) > 0:
            S_AMPA.g_AMPA[:] = self.syn['default']['g_AMPA'] 
        if size(S_GABA.g_GABA) > 0:
            S_GABA.g_GABA[:] = self.syn['default']['g_GABA']
        
        cust_syn = self.syn['custom']
        if size(cust_syn) != 0: 
            for c_s in cust_syn.keys():
                if size(c_s) != 3 : 
                    raise ValueError('Each "customed-synapse dict" key needs be a 3-element tuple '\
                                     '("source_name", "synapse_type", "target_name")')
                src_loc = self.indexOfCell(c_s[0]) 
                tgt_loc = self.indexOfCell(c_s[2])
                val2replace = cust_syn[c_s]
                if c_s[1].upper() == 'ELEC': # SYMMETRICAL ELECTRICAL SYNAPSE 
                    S_elec.g_elec[S_elec.indices[src_loc,tgt_loc]] = val2replace 
                    S_elec.g_elec[S_elec.indices[tgt_loc,src_loc]] = val2replace
                elif c_s[1].upper() == 'AMPA':
                    S_AMPA.g_AMPA[S_AMPA.indices[src_loc,tgt_loc]] = val2replace
                elif c_s[1].upper() == 'GABA':
                    S_GABA.g_GABA[S_GABA.indices[src_loc,tgt_loc]] = val2replace
                else: 
                    raise ValueError('Cannot accept any other values like "%s" representing '\
                                     'synapses except "0", "ELEC", "AMPA" or "GABA"' %(c_s[1]))
        G.set_states(self.init_states)
        
        defaultclock.dt = self.glob['dt']
        monitors = StateMonitor(G,'v',record=True)
        net = Network() 
        net.add(G, S_elec, S_AMPA, S_GABA, monitors)    
        startT = time.time()
        net.run(self.glob['tstop'])
        endT = time.time()
        if reportElapsedTime: 
            print('Total time elapsed = %.3f s' %(endT-startT))        
        self.data['t'] = np.array(monitors.t/ms)
        self.data['v'] = np.array(monitors.v/mV)  
    def plotSimulatedData(self, org, lim_x):
        t = self.data['t']
        cnt = 1 
        if size(org) != 0:
            for i in xrange(shape(org)[0]):
                for j in xrange(shape(org)[1]): 
                    if org[i][j] != '': 
                        subplot(shape(org)[0], shape(org)[1], cnt)
                        loc_ = self.indexOfCell(org[i][j])
                        plot(t, np.transpose(self.data['v'][loc_,:]))
                        title(self.cell_names[loc_][0])
                        xlim(lim_x)
                    cnt += 1
        else: 
            for i in xrange(self.N): 
                subplot(np.ceil(self.N/2),2,i+1)
                plot(t,self.data['v'][i,:])
                title(self.cell_names[i])
                xlim(lim_x)
        tight_layout()
        show()