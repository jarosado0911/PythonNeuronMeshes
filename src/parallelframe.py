import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import math
import renderneuron as rn

print("numpy:     ",np.__version__)
print("maplotlib: ", matplotlib.__version__)

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
        
def get_circle(p):
    npts=8
    t=np.linspace(0,360,npts+1); t=t[:-1]
    t=t*np.pi/180

    x=[]; y=[]; z=[]; r=10;
    for ph in t:
        x.append(r *np.cos(ph) + p[0])
        y.append(r *np.sin(ph) + p[1])
        z.append(p[2])
    return x,y,z
    
def get_pft_frames(G,npts):
    trunks,_=rn.get_trunks(G)
    t=np.linspace(0,360,npts+1); t=t[:-1]
    t=t*np.pi/180
    contour={}
    #centerpt={}
    for ky,tk in zip(trunks.keys(),trunks.values()):
        points,U,V,T=get_UVT(G,tk)
        circle_pts={}
        #cnterpt={}
        for i in range(len(points)):
            rr=G.nodes[tk[i]]['r']
            posxyz=[]
            
            for rad in t:
                xx=points[i][0]+rr*(U[i][0]*np.cos(rad) + V[i][0]*np.sin(rad))
                yy=points[i][1]+rr*(U[i][1]*np.cos(rad) + V[i][1]*np.sin(rad))
                zz=points[i][2]+rr*(U[i][2]*np.cos(rad) + V[i][2]*np.sin(rad))
                posxyz.append(tuple([xx,yy,zz]))
            circle_pts[tk[i]]=posxyz
            #cnterpt[tk[i]]=tuple([points[i][0],points[i][1],points[i][2]])
        contour[ky]=circle_pts 
        #centerpt[ky]=cnterpt
    return contour#,centerpt

def write_1d_ugx(G,filename):
    with open(filename,'w') as f:
        soma_points=[]; neurite_points=[];
        soma_edges=[]; neurite_edges=[];        
        point_cnt=0; edge_cnt=0;
        
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">');
        for i in range(1,len(G)+1):
            if G.nodes[i]['t']==1:
                soma_points.append(point_cnt)
            else:
                neurite_points.append(point_cnt)
                
            x=G.nodes[i]['pos'][0]; y=G.nodes[i]['pos'][1]; z=G.nodes[i]['pos'][2];
            line=str(x)+' '+str(y)+' '+str(z)+' '
            f.write(line); point_cnt+=1
        f.write('</vertices>\n');
        f.write('<edges>');
        for i in range(1,len(G)+1):
            ll=list(G.predecessors(i))
            if not ll:
                continue #pid=-1
            else:
                pid=ll[0]
            line=str(i-1)+' '+str(pid-1)+' '
            
            if G.nodes[i]['t']==1:
                soma_edges.append(edge_cnt)
            else:
                neurite_edges.append(edge_cnt)
            edge_cnt+=1
            
            f.write(line)
        f.write('</edges>\n');
        f.write('<vertex_attachment name="diameter" type="double" passOn="0" global="1">');
        for i in range(1,len(G)+1):
            r=G.nodes[i]['r']
            line=str(r)+' '
            f.write(line)
        f.write('</vertex_attachment>\n');
        f.write('<subset_handler name="defSH">');
        f.write('<subset name="Neurites" color="0.588235 0.588235 1 1" state="0">\n')
        f.write('<vertices>')
        for i in neurite_points: 
            f.write(str(i)+' ')
        f.write('</vertices>\n')
        f.write('<edges>')
        for i in neurite_edges: 
            f.write(str(i)+' ')
        f.write('</edges>\n')
        f.write('</subset>\n')
        f.write('<subset name="Soma" color="1 0 0" state="0">\n')
        f.write('<vertices>')
        for i in soma_points: 
            f.write(str(i)+' ')
        f.write('</vertices>\n')
        f.write('<edges>')
        for i in soma_edges: 
            f.write(str(i)+' ')
        f.write('</edges>\n')
        f.write('</subset>\n')
        f.write('</subset_handler>');
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n');
        f.write('</projection_handler>\n');
        f.write('</grid>');
                
def write_ugx_subsets(cont,G,npts,filename):
    subsets=[]
    
    clr={1: [1,0,0], 2:[0,1,0], 3: [0,0,1], 4:[0.5,0.5,1],5:[1,1,0.5],6:[0.5,1,1] ,7:[1,0.5,0.5]}
    
    for i in range(1,len(G.nodes())+1):
           if G.nodes[i]['t'] not in subsets:
                   subsets.append(G.nodes[i]['t'])
                    
    subset_points={}; subset_edges={}; subset_faces={};
    for i in subsets:
        subset_points[i]=[]; subset_edges[i]=[]; subset_faces[i]=[];
    
    num_of_vertices=0; num_of_edges=0; num_of_faces=0;
    
    with open(filename,'w') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">\n');
        s=''
        
        # add the vertices to the .ugx
        for ky0 in cont.keys(): # iterate through trunks
            klst=list(cont[ky0].keys()) # list of trunk vertices
            for j in range(len(klst)): # iterate through trunk vertices
                for pos in cont[ky0][klst[j]]:
                    if j!=(len(klst)-1):
                        t=G.nodes[klst[j+1]]['t']
                    else:
                        t=G.nodes[klst[j]]['t']
                    subset_points[t].append(num_of_vertices)
                    s+=str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+' '; num_of_vertices+=1;   
        f.write(s); 
        f.write('</vertices>\n');
        
        cur=0; nxt=cur+1;
        ncirpts=npts;
        
        f.write('<edges>')
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                init=cur
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    s=str(cur)+' '+str(nxt)+' '; 
                    f.write(s);
                    cur=nxt; nxt+=1;
                    if j!=(len(klst)-1):
                        t=G.nodes[klst[j+1]]['t']
                    else:
                        t=G.nodes[klst[j]]['t']
                    subset_edges[t].append(num_of_edges)
                    num_of_edges+=1;
                f.write(str(nxt-1)+' '+str(init)+' '); 
                if j!=(len(klst)-1):
                    t=G.nodes[klst[j+1]]['t']
                else:
                    t=G.nodes[klst[j]]['t']
                subset_edges[t].append(num_of_edges)
                num_of_edges+=1;
                cur=nxt; nxt+=1;

        cur=0; nxt=cur+1;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    if j != (len(klst)-1):
                        t=G.nodes[klst[j+1]]['t']                        
                        subset_edges[t].append(num_of_edges)
                        subset_edges[t].append(num_of_edges+1)
                        s=str(cur)+' '+str(cur+ncirpts)+' '; num_of_edges+=1;
                        s+=str(cur)+' '+str(cur+ncirpts-1)+' '; num_of_edges+=1;
                        f.write(s); 
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    t=G.nodes[klst[j+1]]['t']
                    subset_edges[t].append(num_of_edges)
                    subset_edges[t].append(num_of_edges+1)
                    subset_edges[t].append(num_of_edges+2)
                    f.write(str(nxt-1)+' '+str(nxt-1+ncirpts)+' '); num_of_edges+=1;
                    f.write(str(nxt-1)+' '+str(nxt-2+ncirpts)+' '); num_of_edges+=1;
                    f.write(str(nxt-1)+' '+str(nxt)+' '); num_of_edges+=1;
                cur=nxt; nxt+=1;
        f.write('</edges>\n')    
        f.write('<triangles>')
        cur=0; nxt=cur+1;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    if j != (len(klst)-1):
                        t=G.nodes[klst[j+1]]['t']
                        subset_faces[t].append(num_of_faces)
                        subset_faces[t].append(num_of_faces+1)
                        s= str(cur+ncirpts)+' '+str(cur+1)+' '+str(cur)+' '; num_of_faces+=1;
                        s+=str(cur)+' '+str(cur+ncirpts-1)+' '+str(cur+ncirpts)+' '; num_of_faces+=1;
                        
                        f.write(s)
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    t=G.nodes[klst[j+1]]['t']
                    subset_faces[t].append(num_of_faces)
                    subset_faces[t].append(num_of_faces+1)
                    f.write( str(nxt-1+ncirpts)+' '+str(nxt)+' '+str(nxt-1)+' '); num_of_faces+=1;
                    f.write(str(nxt-1)+' '+str(nxt+ncirpts-2)+' '+str(nxt-1+ncirpts)+' '); num_of_faces+=1;
                cur=nxt; nxt+=1;
        f.write('</triangles>\n')
        f.write('<vertex_attachment name="npMapping" type="Mapping" passOn="1" global="1">')
        for ky0 in cont.keys():
            for ky1 in cont[ky0].keys():
                for pos in cont[ky0][ky1]:
                    xc=G.nodes('pos')[ky1][0]; yc=G.nodes('pos')[ky1][1]; zc=G.nodes('pos')[ky1][2];
                    f.write(str(xc)+' '+str(yc)+' '+str(zc)+' ')
        f.write('</vertex_attachment>\n')
        f.write('<subset_handler name="defSH">\n')
        for i in subsets:
            cl='"'+str(clr[i][0])+' '+str(clr[i][1])+' '+str(clr[i][2])+'"'
            f.write('<subset name="Neurites'+str(i)+'" color='+cl+' state="0">\n')
            f.write('<vertices>')
            s=''
            for p in subset_points[i]:
                s+=str(p)+' '
            f.write(s)
            f.write('</vertices>\n')
            f.write('<edges>')
            s=''
            for p in subset_edges[i]:
                s+=str(p)+' '
            f.write(s)
            f.write('</edges>\n')
            f.write('<faces>')
            s=''
            for p in subset_faces[i]:
                s+=str(p)+' '
            f.write(s)
            f.write('</faces>\n')
            f.write('</subset>\n')
        f.write('</subset_handler>\n')
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n');
        f.write('</projection_handler>\n');
        f.write('</grid>');
        
def write_ugx(cont,G,npts,filename):
    cur=0; nxt=cur+1;
    ncirpts=npts;
    
    soma_faces=[]; neurite_faces=[];
    soma_edges=[]; neurite_edges=[];
    soma_points=[]; neurite_points=[];
    
    with open(filename,'w') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">\n');
        s=''
        num_of_vertices=0;
        num_of_edges=0;
        num_of_faces=0;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
            #for ky1 in cont[ky0].keys():
                #for pos in cont[ky0][ky1]:
                 for pos in cont[ky0][klst[j]]:
                    #if G.nodes[ky1]['t']==1:
                    if j != (len(klst)-1):
                        if G.nodes[klst[j+1]]['t']==1:
                            soma_points.append(num_of_vertices);
                        else:
                            neurite_points.append(num_of_vertices);
                    else:
                        if G.nodes[klst[j]]['t']==1:
                            soma_points.append(num_of_vertices);
                        else:
                            neurite_points.append(num_of_vertices);
                    s+=str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+' '; num_of_vertices+=1;
                    
        f.write(s); 
        f.write('</vertices>\n');
        #####################################################################
        f.write('<edges>')
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
            #for ky1 in cont[ky0].keys():
                init=cur
                #npts=len(cont[ky0][ky1])
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    s=str(cur)+' '+str(nxt)+' '; 
                    f.write(s);
                    cur=nxt; nxt+=1;
                    #if G.nodes[ky1]['t']==1:
                    if j != (len(klst)-1):
                        if (G.nodes[klst[j+1]]['t']==1):
                            soma_edges.append(num_of_edges);
                        else:
                            neurite_edges.append(num_of_edges);
                    else:
                        if (G.nodes[klst[j]]['t']==1):
                            soma_edges.append(num_of_edges);
                        else:
                            neurite_edges.append(num_of_edges);
                    num_of_edges+=1;
                    
                f.write(str(nxt-1)+' '+str(init)+' '); 
                if j != (len(klst)-1):
                    if (G.nodes[klst[j+1]]['t']==1):
                        soma_edges.append(num_of_edges);
                    else:
                        neurite_edges.append(num_of_edges);
                else:
                    if (G.nodes[klst[j]]['t']==1):
                        soma_edges.append(num_of_edges);
                    else:
                        neurite_edges.append(num_of_edges);
                num_of_edges+=1;
                cur=nxt; nxt+=1;
           
        cur=0; nxt=cur+1;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    if j != (len(klst)-1):
                        if (G.nodes[klst[j+1]]['t']==1): #and (G.nodes[klst[j]]['t']==1):
                            soma_edges.append(num_of_edges);
                            soma_edges.append(num_of_edges+1);
                        else:
                            neurite_edges.append(num_of_edges);
                            neurite_edges.append(num_of_edges+1);
                        s=str(cur)+' '+str(cur+ncirpts)+' '; num_of_edges+=1;
                        s+=str(cur)+' '+str(cur+ncirpts-1)+' '; num_of_edges+=1;
                        
                        f.write(s); 
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    if G.nodes[klst[j+1]]['t']==1: #and (G.nodes[klst[j]]['t']==1):
                        soma_edges.append(num_of_edges);
                        soma_edges.append(num_of_edges+1);
                        soma_edges.append(num_of_edges+2);
                    else:
                        neurite_edges.append(num_of_edges);
                        neurite_edges.append(num_of_edges+1);
                        neurite_edges.append(num_of_edges+2);
                    f.write(str(nxt-1)+' '+str(nxt-1+ncirpts)+' '); num_of_edges+=1;
                    f.write(str(nxt-1)+' '+str(nxt-2+ncirpts)+' '); num_of_edges+=1;
                    f.write(str(nxt-1)+' '+str(nxt)+' '); num_of_edges+=1;
                cur=nxt; nxt+=1;
        f.write('</edges>\n')
        ##########################################################################################
        f.write('<triangles>')
        cur=0; nxt=cur+1;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    if j != (len(klst)-1):
                        if G.nodes[klst[j+1]]['t']==1:
                            soma_faces.append(num_of_faces);
                            soma_faces.append(num_of_faces+1);
                        else:
                            neurite_faces.append(num_of_faces);
                            neurite_faces.append(num_of_faces+1);
                        s= str(cur+ncirpts)+' '+str(cur+1)+' '+str(cur)+' '; num_of_faces+=1;
                        s+=str(cur)+' '+str(cur+ncirpts-1)+' '+str(cur+ncirpts)+' '; num_of_faces+=1;
                        
                        f.write(s)
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    if G.nodes[klst[j+1]]['t']==1:
                        soma_faces.append(num_of_faces);
                        soma_faces.append(num_of_faces+1);
                    else:
                        neurite_faces.append(num_of_faces);
                        neurite_faces.append(num_of_faces+1);
                    f.write( str(nxt-1+ncirpts)+' '+str(nxt)+' '+str(nxt-1)+' '); num_of_faces+=1;
                    f.write(str(nxt-1)+' '+str(nxt+ncirpts-2)+' '+str(nxt-1+ncirpts)+' '); num_of_faces+=1;
                cur=nxt; nxt+=1;
        f.write('</triangles>\n')
        f.write('<vertex_attachment name="npMapping" type="Mapping" passOn="1" global="1">')
        for ky0 in cont.keys():
            for ky1 in cont[ky0].keys():
                for pos in cont[ky0][ky1]:
                    xc=G.nodes('pos')[ky1][0]; yc=G.nodes('pos')[ky1][1]; zc=G.nodes('pos')[ky1][2];
                    f.write(str(xc)+' '+str(yc)+' '+str(zc)+' ')
        f.write('</vertex_attachment>\n')
        f.write('<subset_handler name="defSH">\n')
        f.write('<subset name="Neurites" color="0.588235 0.588235 1 1" state="0">\n')
        f.write('<vertices>')
        for i in neurite_points: #range(num_of_vertices):
            f.write(str(i)+' ')
        f.write('</vertices>\n')
        f.write('<edges>')
        for i in neurite_edges: #range(1,num_of_edges):
            f.write(str(i)+' ')
        f.write('</edges>\n')
        f.write('<faces>')
        for i in neurite_faces: #range(1,num_of_faces):
            f.write(str(i)+' ')
        f.write('</faces>\n')
        f.write('</subset>\n')
        f.write('<subset name="Soma" color="1 0 0" state="0">\n')
        f.write('<vertices>')
        for i in soma_points: #range(num_of_vertices):
            f.write(str(i)+' ')
        f.write('</vertices>\n')
        f.write('<edges>')
        for i in soma_edges: #range(num_of_vertices):
            f.write(str(i)+' ')
        f.write('</edges>\n')
        f.write('<faces>')
        for i in soma_faces: #range(num_of_vertices):
            f.write(str(i)+' ')
        f.write('</faces>\n')
        f.write('</subset>\n')
        f.write('</subset_handler>\n')
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n');
        f.write('</projection_handler>\n');
        f.write('</grid>');
    
def plot_pft_frames(G):
    trunks,_=rn.get_trunks(G)
    ax = plt.figure(figsize=(18,18)).add_subplot(projection='3d')
    
    xn=[]; yn=[]; zn=[];
    npts=10;
    t=np.linspace(0,360,npts+1); t=t[:-1]
    t=t*np.pi/180
    contour={}
    for ky,tk in zip(trunks.keys(),trunks.values()):
        points,U,V,T=get_UVT(G,tk)
        circle_pts={}
        for i in range(len(points)):
            rr=G.nodes[tk[i]]['r']
            xp=[]; yp=[]; zp=[];
            posxyz=[]
            
            for rad in t:
                xx=points[i][0]+rr*(U[i][0]*np.cos(rad) + V[i][0]*np.sin(rad))
                yy=points[i][1]+rr*(U[i][1]*np.cos(rad) + V[i][1]*np.sin(rad))
                zz=points[i][2]+rr*(U[i][2]*np.cos(rad) + V[i][2]*np.sin(rad))
                xp.append(xx); yp.append(yy); zp.append(zz);
                posxyz.append(tuple([xx,yy,zz]))
            
            circle_pts[tk[i]]=posxyz
            
            xn+=xp; yn+=yp; zn+=zp;
            xp.append(xp[0]); yp.append(yp[0]); zp.append(zp[0]);
            verts = [list(zip(xp,yp,zp))]
            
            ax.add_collection3d(Line3DCollection(verts,colors='b', linewidths=0.35,alpha=0.15))
        contour[ky]=circle_pts
        
    ax.scatter(xn,yn,zn,c='b',s=0.25,alpha=0.15)
    plt.show()
    return contour
    

def plot_pft_vectors(G):
    trunks,_=rn.get_trunks(G)
    ax = plt.figure(figsize=(18,18)).add_subplot(projection='3d')
    
    for tk in trunks.values():
        points,U,V,T=get_UVT(G,tk)
        
        x=[p[0] for p in points]; y=[p[1] for p in points]; z=[p[2] for p in points]
        dux=[dd[0] for dd in U];  duy=[dd[1] for dd in U];  duz=[dd[2] for dd in U];
        dvx=[dd[0] for dd in V];  dvy=[dd[1] for dd in V];  dvz=[dd[2] for dd in V];
        dtx=[dd[0] for dd in T];  dty=[dd[1] for dd in T];  dtz=[dd[2] for dd in T];
        
        ax.quiver(x, y, z, dux, duy, duz, length=1, normalize=True,color='g')
        ax.quiver(x, y, z, dvx, dvy, dvz, length=1, normalize=True,color='b')
        ax.quiver(x, y, z, dtx, dty, dtz, length=1, normalize=True,color='r')
    #ax.view_init(90,-90,0)
    plt.show()