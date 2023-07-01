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
print('JSON:           ',json.__version__,'\n')
print('networkx:  ', nx.__version__)
print('scipy:     ', scipy.__version__)

def write_vrn(FOLDER):
    print('Making VRN ...')
    flist=os.listdir(FOLDER)
    ugx1d=[x for x in flist if '1d' in x]
    ugxtris=[x for x in flist if '3d' in x]

    jsonfilename=FOLDER+'/MetaInfo.json'
    jsonmet={}
    jsonmet['geom1d']=[]
    jsonmet['ARCHIVE']="Rosado"
    jsonmet['SPECIES']="Animal"
    jsonmet['STRAIN']="not reported"

    for i in range(len(ugx1d)):
        g1d={}
        g1d["name"]=ugx1d[i]
        g1d["description"]="1d mesh coarse mesh"
        g1d["refinement"]=str(i)
        g1d["inflations"]=[]
        infl={}
        infl["name"]=ugxtris[i]
        infl["description"]="2d surface mesh"
        infl["inflation"]=str(1.0)
        g1d["inflations"].append(infl)
        jsonmet['geom1d'].append(g1d)

    # Serializing json
    json_object = json.dumps(jsonmet,indent=1)
    # Writing to sample.json
    with open(jsonfilename, "w") as outfile:
        outfile.write(json_object)

    c='zip -j -r '+FOLDER+'.vrn '+FOLDER
    os.system(c)

def spline_neuron(Gin,delta_x):
    G=copy.deepcopy(Gin)
    trunks,T=get_trunks(G)
    n_nodes=len(G.nodes())
    G2=G
        
    for j in trunks.keys():
        x=[]; y=[]; z=[]; r=[]; tp=[];
        lst=trunks[j]
        for n in lst:
            pos=G.nodes[n]['pos']; rad=G.nodes[n]['r']; t=G.nodes[n]['t'];
            x.append(pos[0]); y.append(pos[1]); z.append(pos[2]);
            r.append(rad); tp.append(t);

        new_r= [x/1.0 for x in r]; ave_t= max(tp)

        if len(lst)<=5:
            spl_deg=len(lst)-1
        else:
            spl_deg=3
            
        # get length of trunk 
        xc=np.array(x); yc=np.array(y); zc=np.array(z);
        d=np.diff(xc)**2+np.diff(yc)**2+np.diff(zc)**2
        d=list(d**(0.5)); d.insert(0,0);
        
        # get length, and number of points
        trunk_length=sum(d); num_true_pts=math.ceil(trunk_length/delta_x)
        if num_true_pts <=3:
            num_true_pts=4

        u_fine = np.linspace(0,1,num_true_pts)

        # do interpolation of points
        tck, u = interpolate.splprep([x,y,z],k=spl_deg,s=0)
        x_fine, y_fine, z_fine = interpolate.splev(u_fine, tck)

        # do interpolation of radii    
        tck,u = interpolate.splprep([d,new_r],k=1,s=0)
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
            next_node=lst[j]
            G2.remove_edge(start_node,next_node)
            start_node=next_node

    remove = [node for node,degree in dict(G2.degree()).items() if degree == 0]
    G2.remove_nodes_from(remove)
    
    x=nx.topological_sort(G2)
    mapping = dict(zip(list(x),list(range(1,len(G2.nodes)+1))))
    G2=nx.relabel_nodes(G2,mapping)        
    return G2
    
    
def get_trunks(G):
    T=copy.deepcopy(G)
    # nodes which are branch points
    branch_nodes=[]; end_nodes=[];
    for i in range(1,len(T)+1):
        if T.degree(i)>2:
            branch_nodes.append(i)
        elif T.degree(i)<2:
            end_nodes.append(i)
    
    trunks={}; cnt=1;
    # iterate through branch nodes
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
    for tr in trunks.keys():
        for j in range(1,len(trunks[tr])):
            nx.set_node_attributes(T,{trunks[tr][j]:((tr % 6)+1)},'t')
    return trunks,T

def refine_and_save(G,path,n=5):
    GG=[G]
    if not os.path.exists(path):
        os.mkdir(path)
    save_to_swc(G,path+'/refinement0.swc')
    H=G
    print('refining round... ',end='')
    for i in range(1,n+1):
        filename=path+'/refinement'+str(i)+'.swc'
        print(str(i),'',end='')
        H=refine_once(H)
        GG.append(H)
        save_to_swc(H,filename)
    print('Done!')
    
    return GG

def remove_soma_line(G):
    # Let us make a routine for removing the soma line segment and only keep a single point
    # find all points whose type and succesor types are soma i.e. t=1
    soma_list=[]
    for i in range(1,len(G.nodes)+1):
        if G.nodes[i]['t']==1:
            #print("Node info: ",i,' : ',G.nodes[i])
            pred=list(G.predecessors(i))
            #print("Node succ", pred)
            if len(pred)!=0:
                pre=pred[0]
                if G.nodes[pre]['t']==1:
                    soma_list.append(i)
    G.remove_nodes_from(soma_list)
    x=nx.topological_sort(G)
    mapping = dict(zip(list(x),list(range(1,len(G.nodes)+1))))
    G2=nx.relabel_nodes(G,mapping) 
    return G2

def save_to_swc(G,filename):
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

def print_id_pid(G):
    for i in range(1,len(G)+1):
        ll=list(G.predecessors(i))
        if not ll:
            pid=-1
        else:
            pid=ll[0]
        print('(Id,Pid): (',i,', ',pid,')')

def get_graph(lines):
    """
    This function turns the lines that were read in, into a graph object
    """
    G=nx.DiGraph()
    for item in lines:
        G.add_node(item[0],t=item[1],pos=(item[2],item[3],item[4]),r=item[5])

    # Now let us add the edges
    for item in lines:
        if item[6] != -1: # add edges, but don't try to add edge that has pid = -1
            G.add_edge(item[6],item[0])
    return G

def refine_once(G):
    """
    This function takes in a graph object (G) and refines the geometry once.
    The refine geometry is return as a new graph object (H)
    """
    g_size=len(G.nodes)
    G2=nx.DiGraph()          # initialize empty digraph

    # new indices don't start at index = 0
    new_ind=g_size+1
    # iterate through the edges
    for e in G.edges:
        # get average of two nodes
        pos0 = np.array(G.nodes[e[0]]['pos']); pos1 = np.array(G.nodes[e[1]]['pos']);
        r0   = np.array(G.nodes[e[0]]['r']);   r1 = np.array(G.nodes[e[1]]['r']);
        t0   = np.array(G.nodes[e[0]]['t']);   t1 = np.array(G.nodes[e[1]]['t']);
        pos_ave = (pos0+pos1)/2.0; r_ave = (r0+r1)/2.0;
        t_ave=t1;
        
        # start building the new graph by adding the nodes
        G2.add_node(e[0],t=t0,pos=pos0,r=r0)
        G2.add_node(new_ind,t=t_ave,pos=tuple(pos_ave),r=r_ave)
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

def plot_hines_sparsity(G):
    """
    This function make a sparsity matrix that describes the connectivity of the neuron graph.
    This function takes in a graph object (G) and uses the edges information to
    generate the sparsity matrix.
    """
    rows=[]
    cols=[]
    for e in G.edges:
        rows.append(e[0]); rows.append(e[0]); rows.append(e[1]); rows.append(e[1]);
        cols.append(e[1]); cols.append(e[0]); cols.append(e[0]); cols.append(e[1]);

    rows=np.array(rows)
    cols=np.array(cols)
    data=np.ones(len(G.edges)*4)
    A=csr_matrix((data, (rows, cols)), shape=(len(G.nodes)+1, len(G.nodes)+1)).toarray()
    plt.spy(A,markersize=1.5)
    plt.show()
    
def plot_neuron_r(G):
    """
    This function takes a graph object (G) and generates the plot
    of the graph using the position attribute of the graph nodes and the radius attribute.
    The radius attribute is used to color and set the width of the segments.
    No return value, only a plot is generated. 
    """
    rr=nx.get_node_attributes(G,'r').values()
    maxima=max(rr)
    minima=min(rr)

    norm = matplotlib.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap='viridis')

    e=list(G.edges())

    maxR=G.nodes

    fig = plt.figure(figsize=(6,6))#,dpi=400)
    ax = fig.add_subplot(111, projection="3d")
    for (e1,e2) in e:
        p1=G.nodes[e1]['pos']; p2=G.nodes[e2]['pos'];
        r=G.nodes[e1]['r'];
        t=G.nodes[e1]['t'];
        x=[p1[0],p2[0]];y=[p1[1],p2[1]];z=[p1[2],p2[2]] 
        clr='b'
        if t==1:
            clr='r'
        ax.plot(x,y,z,c=mapper.to_rgba(r),linewidth=r)
    ax.set_axis_off()
    ax.set_zticks([])
    ax.patch.set_edgecolor('black')  
    ax.patch.set_linewidth(1)  
    ax.view_init(90,-90,0)

def plot_neuron(G):
    """
    This function takes a graph object (G) and generates the plot
    of the graph using the position attribute of the graph nodes.
    No return value, only a plot is generated. 
    """
    pos={i:ps for (i,ps) in nx.get_node_attributes(G,'pos').items()}
    rad=[r for r in nx.get_node_attributes(G,'r').values()]
    
    # Extract node and edge positions from the layout
    node_xyz = np.array([pos[v] for v in sorted(G)])
    edge_xyz = np.array([(pos[u], pos[v]) for u, v in G.edges()])
    
    # Create the 3D figure
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(111, projection="3d")

    # Plot the nodes - alpha is scaled by "depth" automatically
    ax.scatter(*node_xyz.T, s=1, c=rad)
 
    # Plot the edges
    for vizedge in edge_xyz:
        ax.plot(*vizedge.T, color="tab:gray")

    def _format_axes(ax):
        """Visualization options for the 3D axes."""
        # Turn gridlines off
        ax.grid(False)
        # Suppress tick labels
        for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
            dim.set_ticks([])
        # Set axes labels
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_axis_off()
        ax.patch.set_edgecolor('black')  
        ax.patch.set_linewidth(1)  
        ax.view_init(90,-90,0)


    _format_axes(ax)
    fig.tight_layout()
    plt.show()