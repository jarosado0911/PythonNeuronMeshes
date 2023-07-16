import os
import numpy as np
import networkx as nx
import copy
import scipy
from scipy.sparse import csr_matrix
from scipy import interpolate
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib
import json
import random
print('JSON:           ',json.__version__)
print('networkx:       ', nx.__version__)
print('scipy:          ', scipy.__version__)
def split_edges(G):
    """
    This function takes in a graph object (G) and refines the geometry once by splitting all edges
    The refine geometry is return as a new graph object (H)
    """
    g_size=len(G.nodes)
    G2=nx.DiGraph()          # initialize empty digraph

    # new indices don't start at index = 0
    new_ind=g_size+1
    # iterate through the edges
    for e in G.edges:
        # get average of two nodes
        pos0 = (G.nodes[e[0]]['pos']); pos1 = (G.nodes[e[1]]['pos']);
        r0   = (G.nodes[e[0]]['r']);   r1 = (G.nodes[e[1]]['r']);
        t0   = (G.nodes[e[0]]['t']);   t1 = (G.nodes[e[1]]['t']);
        pos_ave=tuple([(pos0[0]+pos1[0])*0.5,(pos0[1]+pos1[1])*0.5,(pos0[2]+pos1[2])*0.5 ]);
        #pos_ave = (pos0+pos1)/2.0; 
        r_ave = (r0+r1)/2.0;
        t_ave=t1;
        
        # start building the new graph by adding the nodes
        G2.add_node(e[0],t=t0,pos=pos0,r=r0)
        G2.add_node(new_ind,t=t_ave,pos=pos_ave,r=r_ave)
        G2.add_node(e[1],t=t1, pos=pos1,r=r1)
        
        # then add the edges
        G2.add_edge(e[0],new_ind)
        G2.add_edge(new_ind,e[1])
        new_ind+=1

    # this is the topological sort call --> put into separate function?
    H = nx.DiGraph()
    H.add_nodes_from(sorted(G2.nodes(data=True)))
    H.add_edges_from(G2.edges(data=True))
    x=nx.topological_sort(H)
    mapping = dict(zip(list(x),list(H.nodes)))
    H=nx.relabel_nodes(H,mapping)
    return H

def split_refine(G,N,savetofile):
    """
    This function performs a sequence of split refinements on the geometry and save into .ugx and .swc
    """
    Gnxt=copy.deepcopy(G)
    for i in range(N):
        Gnxt=split_edges(Gnxt)
        swcfilename=savetofile+'split_refinement'+str(i)+'.swc'
        save_1d_neuron_swc(Gnxt,swcfilename)
        Gnxt=read_1d_neuron_swc(swcfilename)
        ugxfilename=swcfilename.replace('swc','ugx')
        write_1d_ugx(Gnxt,ugxfilename)


def spline_refine(G,dx,savetofile):
    """
    This function performs a sequence of spline refinements and saves to .swc and .ugx formats
    """
    for i in range(len(dx)):
        delta=dx[i];
        Gspline=spline_neuron(G,delta)
        swcfilename=savetofile+'spline_refinement'+str(i)+'.swc'
        save_1d_neuron_swc(Gspline,swcfilename)
        Gspline=read_1d_neuron_swc(swcfilename) # for some reason writing neuron directly to .ugx after spline is messed up
        ugxfilename=swcfilename.replace('swc','ugx')
        write_1d_ugx(Gspline,ugxfilename)

def spline_neuron(Gin,delta_x):
    """
    this function performs a spline refinement of the neuron graph geometry
    """
    G=copy.deepcopy(Gin)
    trunks,T=get_trunks(G)
    n_nodes=len(G.nodes())
    G2=G
    
    # we need to spline using the points in a trunk
    for j in trunks.keys():
        x=[]; y=[]; z=[]; r=[]; tp=[];
        lst=trunks[j]
        for n in lst:
            pos=G.nodes[n]['pos']; rad=G.nodes[n]['r']; t=G.nodes[n]['t'];
            x.append(pos[0]); y.append(pos[1]); z.append(pos[2]);
            r.append(rad); tp.append(t);

        new_r= [x/1.0 for x in r]; ave_t=max(tp)

        if len(lst)<=5:
            spl_deg=len(lst)-1
        else:
            spl_deg=3   # using a cubic spline for most of cell
            
        # get length of trunk 
        xc=np.array(x); yc=np.array(y); zc=np.array(z);
        d=np.diff(xc)**2+np.diff(yc)**2+np.diff(zc)**2
        d=list(d**(0.5)); d.insert(0,0);
        
        # get length, and number of points
        trunk_length=sum(d); num_true_pts=math.ceil(trunk_length/delta_x)
        if num_true_pts <=3:
            num_true_pts=4

        u_fine = np.linspace(0,1.0,num_true_pts)

        # do interpolation of points
        tck, u = interpolate.splprep([x,y,z],k=spl_deg,s=0)
        x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)

        # do interpolation of radii   
        if min(new_r)==max(new_r):
            r_fine=[new_r[0] for xx in x_fine]
        else:
            try:
                tck,u = interpolate.splprep([d,new_r],k=1,s=0)
            except:
                tck,u = interpolate.splprep([[max(d),min(d)],[max(new_r),min(new_r)]],k=1,s=0)
            d_fine,r_fine=interpolate.splev(u_fine,tck)
    
        # add none connected nodes to end of graph
        start_node=lst[0]
        for j in range(1,len(x_fine)-1):
            n_nodes+=1
            G2.add_node(n_nodes,t=ave_t,pos=(x_fine[j],y_fine[j],z_fine[j]),r=r_fine[j])
            G2.add_edge(start_node,n_nodes) # add edge from start_node to new_node
            start_node=n_nodes
        # add last edge and set last point type
        G2.add_edge(start_node,lst[-1])
        nx.set_node_attributes(G2,{lst[-1]:ave_t},'t')

        # remove old edges
        start_node=lst[0]
        for j in range(1,len(lst)):
            next_node=lst[j]; G2.remove_edge(start_node,next_node); start_node=next_node

    remove = [node for node,degree in dict(G2.degree()).items() if degree == 0]
    G2.remove_nodes_from(remove)
    
    # do a topological sort of the graph
    x=nx.topological_sort(G2)
    mapping = dict(zip(list(x),list(range(1,len(G2.nodes)+1))))
    G2=nx.relabel_nodes(G2,mapping)        
    return G2

def read_swc(filename):
    """
    This function takes in a filename (.swc file) and reads the lines of the file
    into a variable called lines. The lines are split accordingly.
    """
    with open(filename,'r') as f:
        lines=[]
        boolean_found_line=False
        for line in f:
            # the first line is the soma point
            if '1 1 ' in line: 
                # set this flag and the start reading lines in.
                boolean_found_line=True
            if boolean_found_line:
                ll=line.strip('\n').lstrip(' ') #remove newline and left space
                ll=ll.split(' ')                #split the string on spaces
                for i in range(len(ll)):        #iterate through the line to change strings to values
                    if i==2 or i==3 or i==4 or i==5:
                        ll[i]=float(ll[i])  # the radii and x,y,z need to be floats
                    else:
                        ll[i]=int(ll[i])
                lines.append(ll)
    return lines

def get_graph(lines):
    """
    This function turns the lines that were read in, into a digraph object
    """
    G=nx.DiGraph()
    for item in lines:
        G.add_node(item[0],t=item[1],pos=(item[2],item[3],item[4]),r=item[5])

    # Now let us add the edges
    for item in lines:
        if item[6] != -1: # add edges, but don't try to add edge that has pid = -1
            G.add_edge(item[6],item[0])
    return G

def read_1d_neuron_swc(filename):
    """
    This function reads a .swc file and returns the digraph of the neuron
    with radius, position, and type data on each node
    """
    lines=read_swc(filename)
    return get_graph(lines)

def save_1d_neuron_swc(G,filename):
    """
    This function saves a neuron in .swc format
    """
    with open(filename,'w') as f:
        for i in range(1,len(G)+1):
            ll=list(G.predecessors(i))
            if not ll:
                pid=-1
            else:
                pid=ll[0]
            x=G.nodes[i]['pos'][0]; y=G.nodes[i]['pos'][1]; z=G.nodes[i]['pos'][2];
            r=G.nodes[i]['r']
            t=G.nodes[i]['t']
            line=str(i)+' '+str(t)+' '+str(x)+' '+str(y)+' '+str(z)+' '+str(r)+' '+str(pid)+'\n'
            f.write(line)
    print('Saved to ',filename)

def save_1d_neuron_swc_soma_neurite(G,filename):
    """
    This function saves a neuron in .swc format but with only two subsets
    the soma and neurite subset. All non-soma subsets are put in neurite subset
    """
    with open(filename,'w') as f:
        for i in range(1,len(G)+1):
            ll=list(G.predecessors(i))
            if not ll:
                pid=-1
            else:
                pid=ll[0]
            x=G.nodes[i]['pos'][0]; y=G.nodes[i]['pos'][1]; z=G.nodes[i]['pos'][2];
            r=G.nodes[i]['r']
            t=G.nodes[i]['t']
            if t != 1: t=7; # if not soma save as neurite
            line=str(i)+' '+str(t)+' '+str(x)+' '+str(y)+' '+str(z)+' '+str(r)+' '+str(pid)+'\n'
            f.write(line)
    print('Saved to ',filename)
    
def get_trunks(G):
    """
    This function finds all the trunks (sections) of a neuron
    It returns a dictionary of trunks, which contain the nodes from the graph
    corresponding to that trunk. It also returns a graph structure that is same
    as the original Graph G, but the subsets are all assigned to different values
    """
    T=copy.deepcopy(G)
    # nodes which are branch points
    branch_nodes=[]; end_nodes=[];
    for i in range(1,len(T)+1):
        if T.nodes[int(i)]['t']==1 and T.degree(int(i))<=2: #what if soma has degree 2 or 1 only?
            branch_nodes.append(int(i))
        if T.degree(int(i))>2:
            branch_nodes.append(int(i))
        elif T.degree(int(i))<2:
            end_nodes.append(int(i))
    
    trunks={}; cnt=1;
    # iterate through branch nodes
    #print(branch_nodes)
    for i in branch_nodes:
        # find successors of branch node
        succ=list(T.successors(i))
        for k in range(len(succ)):
            ind=succ[k]; P=[i];
            while True:
                P.append(ind)
                ind=list(T.successors(ind))
                if len(ind)!= 1:
                    break
                ind=ind[0]
            trunks[cnt]=P; cnt+=1;
    tkeys=list(trunks.keys())
    for i in range(len(tkeys)):
        for j in range(1,len(trunks[tkeys[i]])):  #why did I skip the first point in the trunk? --> this causes soma to dissappear?
            nx.set_node_attributes(T,{trunks[tkeys[i]][j]: i+2},'t')  #doing i+2 avoid soma type why?
    return trunks,T

def get_cell_structure(G):
    """
    This function returns the .ugx structure of the Neuron Graph
    """
    # Let us get ugx structure, first get the subset types
    subsets=set()
    for n in list(G.nodes()):
        subsets.add(G.nodes[n]['t'])

    # now collect nodes in particular subset
    node_setdict={}; edge_setdict={};
    attachment_dict={}
    for i in subsets:
        node_setdict[i]=[]; edge_setdict[i]=[]
    for n in list(G.nodes()):
        t=G.nodes[n]['t']
        node_setdict[t].append(n) 
        rad=G.nodes[n]['r']; pos=G.nodes[n]['pos']
        Dict={'r': rad,'pos': pos}
        attachment_dict[n]=Dict   

    for i in subsets:
        node_setdict[i].sort()
    # now make a list of edges, edge set dictionary
    edge_list=[]
       
    for n in list(G.nodes()):
        pred_list=list(G.predecessors(n))
        if len(pred_list)!=0:
            from_to=tuple([pred_list[0], n]) 
            edge_list.append(from_to)
    edge_list.sort()   
     
    edge_number=0 #ugx edge numbers start at 0 
    for e in edge_list:
        from_to=e
        t=G.nodes[from_to[1]]['t']  
        Dict={'e':edge_number,'fromto':from_to}; edge_number+=1
        edge_setdict[t].append(Dict)
    
    # now assemble dictionary containing the subset structures
    cell={}; cell['subsets']=subsets; cell['nodesets']=node_setdict
    cell['edgesets']=edge_setdict; cell['edges']=edge_list; cell['attachments']=attachment_dict
    return cell

def get_cell_structure_by_trunks(G):
    #get the trunks of the geometry
    trunks,T=get_trunks(G)
    
    # Let us get ugx structure, first get the subset types
    subsets=set()
    for n in list(G.nodes()):
        subsets.add(G.nodes[n]['t'])
    
    # now collect nodes in particular subset
    node_setdict={}; edge_setdict={};
    attachment_dict={}
    for i in subsets:
        node_setdict[i]=[]; edge_setdict[i]=[]
    
    trunks_keys=list(trunks.keys())
    vertex_number=1;
    ugx_to_graph={};
    graph_to_ugx={};
    ugx1dcoords ={};
    
    for i in range(len(trunks_keys)): # should I skip the last??
        trunk_nodes=trunks[trunks_keys[i]];
        for j in range(len(trunk_nodes)):
            n=trunk_nodes[j]
            t=G.nodes[n]['t']
            if n in graph_to_ugx.keys():
                continue
            node_setdict[t].append(vertex_number) 
            ugx_to_graph[vertex_number]=n; graph_to_ugx[n]=vertex_number; vertex_number+=1;
            ugx1dcoords[vertex_number]=G.nodes[n]['pos']
            rad=G.nodes[n]['r']; pos=G.nodes[n]['pos']
            Dict={'r': rad,'pos': pos}
            attachment_dict[graph_to_ugx[n]]=Dict
            
    for i in subsets:
        node_setdict[i].sort()
        
    # now make a list of edges, edge set dictionary
    edge_list=[]
    for i in range(len(trunks_keys)): 
        trunk_nodes=trunks[trunks_keys[i]];
        for j in range(1,len(trunk_nodes)): #skip the first node
            n=trunk_nodes[j]
            pred_list=list(G.predecessors(n))
            if len(pred_list)!=0:
                from_to=tuple([graph_to_ugx[pred_list[0]], graph_to_ugx[n]]) 
                edge_list.append(from_to)
    edge_list.sort()
    
    edge_number=0 #ugx edge numbers start at 0 
    for e in edge_list:
        from_to=e
        t=G.nodes[ugx_to_graph[from_to[1]]]['t']  
        Dict={'e':edge_number,'fromto':from_to}; edge_number+=1
        edge_setdict[t].append(Dict)
        
    # now assemble dictionary containing the subset structures
    cell={}; cell['subsets']=subsets; cell['nodesets']=node_setdict
    cell['edgesets']=edge_setdict; cell['edges']=edge_list; cell['attachments']=attachment_dict
    cell['graphugx']=graph_to_ugx;
    cell['ugxgraph']=ugx_to_graph;
    cell['ugx1dcoords']=ugx1dcoords;
    return cell

def convert_to_soma_neurite(G):
    GSN=copy.deepcopy(G);
    for n in GSN.nodes():
        if GSN.nodes[n]['t']!=1: GSN.nodes[n]['t']=7  # assign all non soma subsets to neurite subset
    return GSN
    
def write_1d_ugx_soma_neurite(G,filename):
    """
    This function write the neuron graph to .ugx form but only with two subets
    the soma and neurite subset
    """
    GSN=convert_to_soma_neurite(G)
    write_1d_ugx(GSN,filename)
    print('Saved to ',filename)
    
def write_1d_ugx(G,filename):
    """
    This function writes the neuron graph to .ugx preserving all subsets by usings the
    get_cell_structure function call
    """
    cell = get_cell_structure_by_trunks(G)
    #cell=get_cell_structure(G)
    with open(filename,'w') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">');
        s=''
        for n in cell['attachments'].keys():
            pos=cell['attachments'][n]['pos'];
            s+=str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+' '
        s=s[:-1]; f.write(s);
        f.write('</vertices>\n');
        f.write('<edges>');
        s=''
        for e in cell['edges']:
            s+=str(e[0]-1)+' '+str(e[1]-1)+' '
        s=s[:-1]; f.write(s)
        f.write('</edges>\n');
        f.write('<vertex_attachment name="diameter" type="double" passOn="0" global="1">');
        s=''
        for n in cell['attachments'].keys():
            s+=str(cell['attachments'][n]['r'])+' ';
        s=s[:-1]; f.write(s);
        f.write('</vertex_attachment>\n');
        f.write('<subset_handler name="defSH">');
        for ss in cell['subsets']:
            name='subset'; clr=tuple([random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)])
            if ss==1: name='soma'
            if ss==2: name='axon'
            if ss==3: name='dend'
            if ss==4: name='apic'
            if ss==5: name='fork'
            if ss==6: name='end'
            if ss==7: name='neurite'
            f.write('<subset name="'+name+'" color="'+str(clr[0])+' '+str(clr[1])+' '+str(clr[2])+'" state="0">\n')
            f.write('<vertices>')
            s=''
            for n in cell['nodesets'][ss]: s+=str(n-1)+' '
            s=s[:-1]; f.write(s)  # remove the last space
            f.write('</vertices>\n')
            f.write('<edges>')
            s=''
            for e in cell['edgesets'][ss]: s+=str(e['e'])+' '
            s=s[:-1]; f.write(s)
            f.write('</edges>\n'); 
            f.write('</subset>\n')
        #end for loop
        f.write('</subset_handler>'); 
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n'); 
        f.write('</projection_handler>\n');
        f.write('</grid>');
    print('Saved to ',filename)
    
def remove_soma_line(G):
    # Let us make a routine for removing the soma line segment and only keep a single point
    # find all points whose type and succesor types are soma i.e. t=1
    soma_list=[]
    for i in range(1,len(G.nodes)+1):
        if G.nodes[i]['t']==1:
            pred=list(G.predecessors(i))
            if len(pred)!=0:
                pre=pred[0]
                if G.nodes[pre]['t']==1: soma_list.append(i)
    G.remove_nodes_from(soma_list)
    x=nx.topological_sort(G)
    mapping = dict(zip(list(x),list(range(1,len(G.nodes)+1))))
    G2=nx.relabel_nodes(G,mapping) 
    return G2