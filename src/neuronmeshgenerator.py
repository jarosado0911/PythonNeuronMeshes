import numpy as np
import neuronswc as ns
import networkx as nx
import math
import random
import os
import json

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

def rotationMatrix(axis, theta):
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def get_UVT(G,tk):
    points=[] # gather points
    for i in tk:
        pos=list(G.nodes[i]['pos'])
        points.append(pos)

    n = len(points)  # Number of points
    T = np.apply_along_axis(np.gradient, axis=0, arr=points) # Calculate all tangents 

    # Normalize all tangents
    f = lambda m : m / np.linalg.norm(m)
    T = np.apply_along_axis(f, axis=1, arr=T)

    # Initialize the first parallel-transported normal vector V
    V = np.zeros(np.shape(points))
    V[0] = (T[0][1], -T[0][0], 0)
    V[0] = V[0] / np.linalg.norm(V[0])

    # Compute the values for V for each tangential vector from T
    for i in range(n - 1):
        b = np.cross(T[i], T[i + 1])
        if np.linalg.norm(b) < 0.00001:
            V[i + 1] = V[i]
        else:
            b = b / np.linalg.norm(b)
            phi = np.arccos(np.dot(T[i], T[i + 1]))
            R = rotationMatrix(b,phi)
            V[i + 1] = np.dot(R, V[i])

    # Calculate the second parallel-transported normal vector U
    U = np.array([np.cross(t, v) for (t, v) in zip(T, V)])
    
    return points,U,V,T

def get_pft_frames(G,npts):
    trunks,_=ns.get_trunks(G)
    t=np.linspace(0,360,npts+1); t=t[:-1]
    t=t*np.pi/180; contour={}
    for ky,tk in zip(trunks.keys(),trunks.values()):
        points,U,V,T=get_UVT(G,tk)
        circle_pts={}
        for i in range(len(points)):
            rr=G.nodes[tk[i]]['r']
            if G.degree(tk[i])==1:
                rr=0.00001
                
            posxyz=[]
            for rad in t:
                xx=points[i][0]+rr*(U[i][0]*np.cos(rad) + V[i][0]*np.sin(rad))
                yy=points[i][1]+rr*(U[i][1]*np.cos(rad) + V[i][1]*np.sin(rad))
                zz=points[i][2]+rr*(U[i][2]*np.cos(rad) + V[i][2]*np.sin(rad))
                posxyz.append(tuple([xx,yy,zz]))
            circle_pts[tk[i]]=posxyz
        contour[ky]=circle_pts 
    return contour

#################################################################
################### General Mesh Functions ######################
#################################################################

def make_meshes(neuron_filename,start_dx,num_ref,n_cir_points,sph_contours,sph_points):
    MESH_FOLDER=neuron_filename.split('/')
    MESH_FOLDER=MESH_FOLDER[-1]
    MESH_FOLDER=MESH_FOLDER.replace('swc','mesh')
    cell_name=MESH_FOLDER.replace('.mesh','')
    # make output folder
    if not os.path.exists(MESH_FOLDER):
        print('folder not found')
        os.makedirs(MESH_FOLDER)


    G=ns.read_1d_neuron_swc(neuron_filename)     # read the .swc file
    G=ns.remove_soma_line(G)              # remove soma line

    dx=[start_dx*pow(2.0,p) for p in range(0,-num_ref,-1)]
    for i in range(len(dx)):
        Grefine=ns.spline_neuron(G,dx[i])
        swcfilename=MESH_FOLDER+'/'+cell_name+'refinement'+str(i)+'.swc'
        ns.save_1d_neuron_swc_soma_neurite(Grefine,swcfilename)
        Grefine=ns.read_1d_neuron_swc(swcfilename)
        
        #reorder nodes into dfs ordering
        #old_labels=list(set(Grefine.nodes()))
        #new_labels=list(nx.dfs_preorder_nodes(Grefine,source=1))
        #mapping=dict(zip(old_labels,new_labels))
        #Grefine=nx.relabel_nodes(Grefine,mapping)
        
        os.remove(swcfilename)
        ugx1dfilename=swcfilename.replace('.swc','.ugx')
        ugx1dfilename=ugx1dfilename.replace('refinement','_1d_refinement')
        ns.write_1d_ugx(Grefine,ugx1dfilename)
        ugxmeshfilename=ugx1dfilename.replace('_1d_refinement','_3d_mesh_refinment')
        write_mesh_with_soma_sphere_ugx(Grefine,n_cir_points,sph_contours,sph_points,ugxmeshfilename)
    return MESH_FOLDER


def get_mesh_structure(G,ncirclepoints=6):
    """
    This function returns the .ugx surface mesh structure of the Neuron Graph
    a) returns the list of vertices of the surface mesh in their corresponding subsets
    Note when using ProMesh4 set the remove duplicates threshold to 0 if you set to anything else then it will remove vertices
    but what is the limit on how close the decimal values for coordinates should be??
    """
    # Let us get ugx structure, first get the subset types
    subsets=set()
    for n in list(G.nodes()):
        subsets.add(G.nodes[n]['t'])
        
    # now collect nodes in particular subset
    vertex_setdict={}; edge_setdict={}; face_setdict={};
    attachment_dict={}
    for i in subsets:
        vertex_setdict[i]=[]; edge_setdict[i]=[]; face_setdict[i]=[];
        
    trunks,_=ns.get_trunks(G)  # get the trunks
    contour=get_pft_frames(G,ncirclepoints)  # get the cont points
    
    coords=[]; vertex_number=1;
    for trunk in contour:
        keylist=list(contour[trunk].keys())
        for i in range(0,len(keylist)): #should I skip first contour or last, for now all? 
            for xyz in contour[trunk][keylist[i]]:
                coords.append(tuple([xyz[0],xyz[1],xyz[2]]))
                if i<len(keylist)-1:
                    t=G.nodes[keylist[i+1]]['t']
                else:
                    t=G.nodes[keylist[i]]['t']
                if t==1:
                    try:
                        t=G.nodes[keylist[i+1]]['t']
                    except IndexError:
                        t=G.nodes[keylist[i]]['t']
                vertex_setdict[t].append(vertex_number); vertex_number+=1
      
    edge_number=0 # ugx starts at 0
    edge_list=[];
    template=np.array(range(ncirclepoints))
    temp1=template; temp2=np.roll(template,1); temp3=temp2+ncirclepoints;
    for trunk in contour:
        keylist=list(contour[trunk].keys())
        for i in range(0,len(keylist)):  #should I skip first contour?
            if i<len(keylist)-1:
                t=G.nodes[keylist[i+1]]['t']
            else:
                t=G.nodes[keylist[i]]['t']
            if t==1:
                try:
                    t=G.nodes[keylist[i+1]]['t']
                except IndexError:
                    t=G.nodes[keylist[i]]['t']
            for e1,e2,e3 in zip(temp1,temp2,temp3):
                edge_list.append(tuple([e1,e2]))
                edge_setdict[t].append(edge_number); edge_number+=1
                if i != len(keylist)-1:
                    edge_list.append(tuple([e1,e1+ncirclepoints]))
                    edge_setdict[t].append(edge_number); edge_number+=1
                    edge_list.append(tuple([e1,e3]))
                    edge_setdict[t].append(edge_number); edge_number+=1
            temp1+=ncirclepoints; temp2+=ncirclepoints; temp3=temp2+ncirclepoints;
    
    face_number=0; face_list=[]
    template=np.array(range(ncirclepoints))
    temp1=template; temp2=temp1+ncirclepoints; temp3=np.roll(temp2,1); temp4=np.roll(temp1,1);
    
    face_mapping=[];
    
    for trunk in contour:
        keylist=list(contour[trunk].keys())
        for i in range(0,len(keylist)):
            for v1,v2,v3,v4 in zip(temp1,temp2,temp3,temp4):
                if i<len(keylist)-1:
                    nx=list(G.predecessors(keylist[i])); 
                    if len(nx) != 0:
                        nx=nx[0]
                    else:
                        nx=keylist[i+1]
                    nxt=G.nodes[keylist[i]]['pos']; cur=G.nodes[nx]['pos'];
                    f1=coords[v1]; f2=coords[v2]; f3=coords[v3]; LL=compute_lambda(cur,nxt,f1,f2,f3);
                    face_mapping.append(tuple([cur,nxt,LL]))
                if i==len(keylist)-1:
                    nx=list(G.predecessors(keylist[i])); nx=nx[0]; nxt=G.nodes[nx]['pos'];
                    cur=G.nodes[keylist[i]]['pos']
                    LL=compute_lambda(cur,nxt,f1,f2,f3);
                    face_mapping.append(tuple([cur,nxt,LL]))
                    
            temp1+=ncirclepoints; temp2=temp1+ncirclepoints; temp3=np.roll(temp2,1); temp4=np.roll(temp1,1);
            
    
    template=np.array(range(ncirclepoints))
    temp1=template; temp2=temp1+ncirclepoints; temp3=np.roll(temp2,1); temp4=np.roll(temp1,1);
    for trunk in contour:
        keylist=list(contour[trunk].keys())
        for i in range(0,len(keylist)):  #should I skip first contour?
            if i<len(keylist)-1:
                t=G.nodes[keylist[i+1]]['t']
            else:
                t=G.nodes[keylist[i]]['t']
            if t==1:
                try:
                    t=G.nodes[keylist[i+1]]['t']
                except IndexError:
                    t=G.nodes[keylist[i]]['t']
            for v1,v2,v3,v4 in zip(temp1,temp2,temp3,temp4):
                if i != len(keylist)-1:
                    nx=list(G.predecessors(keylist[i])); 
                    if len(nx) != 0:
                        nx=nx[0]
                    else:
                        nx=keylist[i+1]
                        
                    nxt=G.nodes[keylist[i]]['pos']; cur=G.nodes[nx]['pos'];
                    
                    #first face
                    #f1=coords[v1]; f2=coords[v2]; f3=coords[v3];
                    #LL=compute_lambda(cur,nxt,f1,f2,f3);
                    face_list.append(tuple([v3,v2,v1])); #face_mapping.append(tuple([cur,nxt,LL]))
                    face_setdict[t].append(face_number); face_number+=1
                    
                    #second face
                    #f1=coords[v1]; f2=coords[v3]; f3=coords[v4];
                    #LL=compute_lambda(cur,nxt,f1,f2,f3); 
                    face_list.append(tuple([v4,v3,v1])); #face_mapping.append(tuple([cur,nxt,LL]))
                    face_setdict[t].append(face_number); face_number+=1
            temp1+=ncirclepoints; temp2=temp1+ncirclepoints; temp3=np.roll(temp2,1); temp4=np.roll(temp1,1);
    
    surface={}
    surface['subsets']=subsets; surface['trunks']=trunks;
    surface['coords']=coords;   surface['vertexsets']=vertex_setdict
    surface['edges']=edge_list; surface['edgesets']=edge_setdict;
    surface['faces']=face_list; surface['facesets']=face_setdict;
    surface['fmapping']=face_mapping;
    return surface

def compute_lambda(cur,nxt,f1,f2,f3):
    cr=np.array(cur); nx=np.array(nxt);
    fc1=np.array(f1); fc2=np.array(f2); fc3=np.array(f3);
    
    centroid_point= (fc1+fc2+fc3)/3;
    dotprod1=np.dot((nx-cr),(centroid_point-cr));
    dotprod2=np.dot((nx-cr),(nx-cr));
    
    return (dotprod1/dotprod2);

def write_mesh_ugx_soma_neurite(G,ncirclepoints,filename):
    GSN=ns.convert_to_soma_neurite(G)
    write_mesh_ugx(GSN,ncirclepoints,filename)

def write_mesh_ugx(G,ncirclepoints,filename):
    """
    This function writes the neuron graph to .ugx preserving all subsets by usings the
    get_cell_structure function call
    """
    surface=get_mesh_structure(G,ncirclepoints)
    with open(filename,'w') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">');
        s=''
        for xyz in surface['coords']:
            for crd in xyz:
                s+=str(crd)+' '
        s=s[:-1]; f.write(s);
        f.write('</vertices>\n');
        f.write('<edges>');
        s=''
        for e in surface['edges']:
            s+=str(e[0])+' '+str(e[1])+' '
        s=s[:-1]; f.write(s)
        f.write('</edges>\n');
        f.write('<triangles>')
        s=''
        for fc in surface['faces']:
            s+=str(fc[0])+' '+str(fc[1])+' '+str(fc[2])+' '
        s=s[:-1]; f.write(s)
        f.write('</triangles>\n')
        f.write('<subset_handler name="defSH">');
        for ss in surface['subsets']:
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
            for n in surface['vertexsets'][ss]: s+=str(n-1)+' '  # ugx starts at 0
            s=s[:-1]; f.write(s)  # remove the last space
            f.write('</vertices>\n')
            f.write('<edges>')
            s=''
            for e in surface['edgesets'][ss]: s+=str(e)+' '
            s=s[:-1]; f.write(s)
            f.write('</edges>\n');
            f.write('<faces>')
            s=''
            for fc in surface['facesets'][ss]: s+=str(fc)+' '
            s=s[:-1]; f.write(s)
            f.write('</faces>\n')
            f.write('</subset>\n')
        f.write('</subset_handler>'); 
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n'); 
        f.write('</projection_handler>\n');
        f.write('</grid>');
    print('Saved to ',filename)
    
    
def write_mesh_with_soma_sphere_ugx(G,ncirclepoints,n_spcontours,n_spcircle_points,filename):
    """
    This function writes the neuron graph to .ugx preserving all subsets by usings the
    get_cell_structure function call
    """
    surface=get_mesh_structure(G,ncirclepoints); m=len(surface['coords']); l=len(surface['edges']); p=len(surface['faces']);
    # note this function assumes that there is NO soma line!
    radius=G.nodes[1]['r']*1.1
    soma_center=list(G.nodes[1]['pos'])
    spheresurface=sphere_mesh_structure(n_spcontours,n_spcircle_points,radius,soma_center)
    with open(filename,'w') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">');
        s=''
        for xyz in surface['coords']:
            for crd in xyz:
                s+=str(crd)+' '
        for xyz in spheresurface['coords']:
            for crd in xyz:
                s+=str(crd)+' '
        s=s[:-1]; f.write(s);
        f.write('</vertices>\n');
        f.write('<edges>');
        s=''
        for e in surface['edges']:
            s+=str(e[0])+' '+str(e[1])+' '
        for e in spheresurface['edges']:
            s+=str(e[0]+m)+' '+str(e[1]+m)+' '
        s=s[:-1]; f.write(s)
        f.write('</edges>\n');
        f.write('<triangles>')
        s=''
        for fc in surface['faces']:
            s+=str(fc[0])+' '+str(fc[1])+' '+str(fc[2])+' '
        for fc in spheresurface['faces']:
            s+=str(fc[0]+m)+' '+str(fc[1]+m)+' '+str(fc[2]+m)+' '
        s=s[:-1]; f.write(s)
        f.write('</triangles>\n')
        f.write('<vertex_attachment name="npMapping" type="Mapping" passOn="1" global="1">')
        s=''
        for mp in surface['fmapping']:
            pos1=mp[0]; pos2=mp[1]; LL=mp[2];
            s+=str(pos1[0])+' '+str(pos1[1])+' '+str(pos1[2])+' '+str(pos2[0])+' '+str(pos2[1])+' '+str(pos2[2])+' '+str(LL)+' ';
        for mp in spheresurface['fmapping']:
            pos1=mp[0]; pos2=mp[1]; LL=mp[2];
            s+=str(pos1[0])+' '+str(pos1[1])+' '+str(pos1[2])+' '+str(pos2[0])+' '+str(pos2[1])+' '+str(pos2[2])+' '+str(LL)+' ';
        s=s[:-1]; f.write(s)
        f.write('</vertex_attachment>\n')
        f.write('<subset_handler name="defSH">');
        for ss in surface['subsets']:
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
            s=''; 
            if ss==1:
                for n in spheresurface['vertexsets']: s+=str(n+m-1)+' '
            else:
                for n in surface['vertexsets'][ss]: s+=str(n-1)+' '  # ugx starts at 0
            s=s[:-1]; f.write(s)  # remove the last space
            f.write('</vertices>\n')
            f.write('<edges>')
            s=''
            if ss==1:
                for e in spheresurface['edgesets']: s+=str(e+l)+' '
            else:
                for e in surface['edgesets'][ss]: s+=str(e)+' '
            s=s[:-1]; f.write(s)
            f.write('</edges>\n');
            f.write('<faces>')
            s=''
            if ss==1:
                for fc in spheresurface['facesets']: s+=str(fc+p)+' '
            else:
                for fc in surface['facesets'][ss]: s+=str(fc)+' '
            s=s[:-1]; f.write(s)
            f.write('</faces>\n')
            f.write('</subset>\n')
        f.write('</subset_handler>'); 
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n'); 
        f.write('</projection_handler>\n');
        f.write('</grid>');
    print('Saved to ',filename)

#################################################################
###################### Soma Mesh Functions ######################
#################################################################
def sphere_contour_circles(PHI,R,theta_n,c):
    thetas=np.linspace(0,2*np.pi,theta_n+1)
    
    cx=c[0]; cy=c[1]; cz=c[2]
    
    points=[];
    for THETA in thetas:
        x = R * np.sin(PHI) * np.cos(THETA)+cx
        y = R * np.sin(PHI) * np.sin(THETA)+cy
        z = R * np.cos(PHI)+cz
        points.append(tuple([x,y,z]))
    return points

def get_sphere_contours(n_contours,n_circle_points,radius,soma_center):
    cont={}
    cont[0]=[tuple([soma_center[0],soma_center[1],soma_center[2]+radius])];
    for i in range(1,n_contours+1):
        cont[i]=sphere_contour_circles(i*np.pi/(n_contours+1),radius,n_circle_points,soma_center)
        cont[i]=cont[i][:-1]
    cont[n_contours+1]=[tuple([soma_center[0],soma_center[1],soma_center[2]-radius])];
    return cont
    
def sphere_mesh_structure(n_contours,n_circle_points,radius,soma_center):
    points=[];
    num_vertices=n_contours*n_circle_points+2;
    
    cont=get_sphere_contours(n_contours,n_circle_points,radius,soma_center)
    chk=0
        
    # add coordinates
    coords=[]; vertex_number=1;
    vertexset=[];
    for c in cont.keys():
        for xyz in cont[c]:
            coords.append(xyz)
            vertexset.append(vertex_number); vertex_number+=1
    
    # add edges
    template=np.array(range(n_circle_points))+1 #start template at 1 because the first index (0) corresponds to top point of sphere
    temp1=template; temp2=np.roll(template,1); temp3=temp2+n_circle_points;
    edge_list=[]; edge_setdict=[]; edgenumber=0;
    contour_keys=cont.keys()
    for i in range(0,len(contour_keys)):
        if i==0:
            for v2 in temp2:
                edge_list.append(tuple([0,v2])); edge_setdict.append(edgenumber); edgenumber+=1;
        elif i==len(contour_keys)-1:
            for v1 in temp1-n_circle_points:
                edge_list.append(tuple([v1,num_vertices-1])); edge_setdict.append(edgenumber); edgenumber+=1;
        else:
            for v1,v2,v3 in zip(temp1,temp2,temp3):
                edge_list.append(tuple([v1,v2])); edge_setdict.append(edgenumber); edgenumber+=1;
                if i != len(contour_keys)-2:
                    edge_list.append(tuple([v1,v1+n_circle_points])); edge_setdict.append(edgenumber); edgenumber+=1;
                    edge_list.append(tuple([v1,v3])); edge_setdict.append(edgenumber); edgenumber+=1;
            temp1+=n_circle_points; temp2=np.roll(temp1,1); temp3=temp2+n_circle_points;
    
    # add faces
    face_number=0; face_list=[]; face_setdict=[]; facenumber=0;
    template=np.array(range(n_circle_points))+1
    temp1=template; temp2=temp1+n_circle_points; temp3=np.roll(temp2,1); temp4=np.roll(temp1,1);
    face_mapping=[]; 
    LL=0.0; somapos=tuple([soma_center[0],soma_center[1],soma_center[2]]);
    for i in range(0,len(contour_keys)):
        if i==0:
            for v1,v4 in zip(temp1,temp4):
                face_list.append(tuple([v4,v1,0])); face_setdict.append(facenumber); facenumber+=1;
                face_mapping.append(tuple([somapos,somapos,LL]))
        elif i==len(contour_keys)-1:
            for v1,v4 in zip(temp1-n_circle_points,temp4-n_circle_points):
                face_list.append(tuple([num_vertices-1,v1,v4])); face_setdict.append(facenumber); facenumber+=1;
                face_mapping.append(tuple([somapos,somapos,LL]))
        else:
            for v1,v2,v3,v4 in zip(temp1,temp2,temp3,temp4):
                if i != len(contour_keys)-2:
                    face_list.append(tuple([v3,v2,v1]));face_setdict.append(facenumber); facenumber+=1;
                    face_list.append(tuple([v3,v1,v4]));face_setdict.append(facenumber); facenumber+=1;
                    face_mapping.append(tuple([somapos,somapos,LL]))
                    #face_mapping.append(tuple([somapos,somapos,LL]))
                    
            temp1+=n_circle_points; temp2=temp1+n_circle_points; temp3=np.roll(temp2,1); temp4=np.roll(temp1,1);
        
    surface={}
    surface['coords']=coords; surface['vertexsets']=vertexset;
    surface['edges']=edge_list; surface['edgesets']=edge_setdict;
    surface['faces']=face_list; surface['facesets']=face_setdict;
    surface['contours']=cont
    surface['fmapping']=face_mapping;
    return surface

def write_soma_ugx(filename,n_contours,n_circle_points,radius,soma_center):
    surface=sphere_mesh_structure(n_contours,n_circle_points,radius,soma_center)
    with open(filename,'w') as f:        
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">');
        s=''
        for xyz in surface['coords']:
            s+=str(xyz[0])+' '+str(xyz[1])+' '+str(xyz[2])+' '
        s=s[:-1];f.write(s);
        f.write('</vertices>\n');
        f.write('<edges>')
        s=''
        for e in surface['edges']:
            s+=str(e[0])+' '+str(e[1])+' '
        s=s[:-1];f.write(s);
        f.write('</edges>\n')
        f.write('<triangles>')
        s=''
        for fc in surface['faces']:
            s+=str(fc[0])+' '+str(fc[1])+' '+str(fc[2])+' '
        s=s[:-1]
        f.write(s)
        f.write('</triangles>\n')
        f.write('<subset_handler name="defSH">\n')
        f.write('<subset name="Soma" color="1 0 0" state="0">\n')
        f.write('<vertices>')
        s=''
        for v in surface['vertexsets']: s+=str(v-1)+' '
        s=s[:-1];f.write(s);
        f.write('</vertices>\n')
        f.write('<edges>')
        s=''
        for e in surface['edgesets']: s+=str(e)+' '
        s=s[:-1];f.write(s);
        f.write('</edges>\n')
        f.write('<faces>')
        s=''
        for fc in surface['facesets']: s+=str(fc)+' '
        s=s[:-1];f.write(s);
        f.write('</faces>\n')        
        f.write('</subset>\n')
        f.write('</subset_handler>\n')
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n');
        f.write('</projection_handler>\n');
        f.write('</grid>');
    