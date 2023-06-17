import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
import math
import renderneuron as rn

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
    for ky,tk in zip(trunks.keys(),trunks.values()):
        points,U,V,T=get_UVT(G,tk)
        circle_pts={}
        for i in range(len(points)):
            rr=G.nodes[tk[i]]['r']
            posxyz=[]
            
            for rad in t:
                xx=points[i][0]+rr*(U[i][0]*np.cos(rad) + V[i][0]*np.sin(rad))
                yy=points[i][1]+rr*(U[i][1]*np.cos(rad) + V[i][1]*np.sin(rad))
                zz=points[i][2]+rr*(U[i][2]*np.cos(rad) + V[i][2]*np.sin(rad))
                posxyz.append(tuple([xx,yy,zz]))
            circle_pts[tk[i]]=posxyz
        contour[ky]=circle_pts 
    return contour

def write_ugx(cont,npts,filename):
    cur=0; nxt=cur+1;
    ncirpts=npts;

    with open(filename,'w') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">\n');
        s=''
        for ky0 in cont.keys():
            for ky1 in cont[ky0].keys():
                for pos in cont[ky0][ky1]:
                    s+=str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+' '
        f.write(s)
        f.write('</vertices>\n');
        f.write('<edges>')
        for ky0 in cont.keys():
            for ky1 in cont[ky0].keys():
                init=cur
                npts=len(cont[ky0][ky1])
                for i in range(npts-1):
                    s=str(cur)+' '+str(nxt)+' '
                    f.write(s)
                    cur=nxt; nxt+=1;
                f.write(str(nxt-1)+' '+str(init)+' ')
                cur=nxt; nxt+=1;

        cur=0; nxt=cur+1;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    s=str(cur)+' '+str(cur+ncirpts)+' '
                    if j != (len(klst)-1):
                        f.write(s)
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    f.write(str(nxt-1)+' '+str(nxt-1+ncirpts)+' ')
                cur=nxt; nxt+=1;
        f.write('</edges>\n')
        f.write('<triangles>')
        cur=0; nxt=cur+1;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    s=str(cur)+' '+str(cur+1)+' '+str(cur+ncirpts)+' '
                    s+=str(cur)+' '+str(cur+ncirpts-1)+' '+str(cur+ncirpts)+' '
                    if j != (len(klst)-1):
                        f.write(s)
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    f.write(str(nxt-1)+' '+str(nxt)+' '+str(nxt-1+ncirpts)+' ')
                    f.write(str(nxt-1)+' '+str(nxt+ncirpts-2)+' '+str(nxt-1+ncirpts)+' ')
                cur=nxt; nxt+=1;
        f.write('</triangles>\n')
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