import numpy as np

def contour_circles(PHI,R,theta_n,c):
    thetas=np.linspace(0,2*np.pi,theta_n+1)
    
    cx=c[0]; cy=c[1]; cz=c[2]
    
    points=[];
    for THETA in thetas:
        x = R * np.sin(PHI) * np.cos(THETA)+cx
        y = R * np.sin(PHI) * np.sin(THETA)+cy
        z = R * np.cos(PHI)+cz
        points.append(tuple([x,y,z]))
    return points

def write_ugx(cont,G,npts,filename,scs,nsphere_contours,nsphere_contour_pts):
    cur=0; nxt=cur+1;
    ncirpts=npts;
    
    soma_faces=[]; neurite_faces=[];
    soma_edges=[]; neurite_edges=[];
    soma_points=[]; neurite_points=[];
    
    somafilename='soma.ugx'
    for nd in G.nodes:
        if G.nodes[nd]['t']==1:
            soma_rad=G.nodes[nd]['r'];
            soma_center=G.nodes[nd]['pos'];
    points, edgelist, facelist=write_ugx_sphere(somafilename,nsphere_contours,nsphere_contour_pts,soma_rad*scs,soma_center)
    
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
                 for pos in cont[ky0][klst[j]]:
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
        
        print(num_of_vertices)
        start_index=num_of_vertices;
        for pos in points:
            s+=str(pos[0])+' '+str(pos[1])+' '+str(pos[2])+' '; num_of_vertices+=1
            soma_points.append(num_of_vertices-1);
        
        s=s[:-1]            
        f.write(s); 
        f.write('</vertices>\n');
        #####################################################################
        f.write('<edges>')
        s=''
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                init=cur
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    s+=str(cur)+' '+str(nxt)+' '; 
                    #f.write(s);
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
                    
                #f.write(str(nxt-1)+' '+str(init)+' '); 
                s+=str(nxt-1)+' '+str(init)+' '
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
                        if (G.nodes[klst[j+1]]['t']==1): 
                            soma_edges.append(num_of_edges);
                            soma_edges.append(num_of_edges+1);
                        else:
                            neurite_edges.append(num_of_edges);
                            neurite_edges.append(num_of_edges+1);
                        s+=str(cur)+' '+str(cur+ncirpts)+' '; num_of_edges+=1;
                        s+=str(cur)+' '+str(cur+ncirpts-1)+' '; num_of_edges+=1;
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    if G.nodes[klst[j+1]]['t']==1: 
                        soma_edges.append(num_of_edges);
                        soma_edges.append(num_of_edges+1);
                        soma_edges.append(num_of_edges+2);
                    else:
                        neurite_edges.append(num_of_edges);
                        neurite_edges.append(num_of_edges+1);
                        neurite_edges.append(num_of_edges+2);
                    s+=str(nxt-1)+' '+str(nxt-1+ncirpts)+' '; num_of_edges+=1;
                    s+=str(nxt-1)+' '+str(nxt-2+ncirpts)+' '; num_of_edges+=1;
                    s+=str(nxt-1)+' '+str(nxt)+' '; num_of_edges+=1;
                cur=nxt; nxt+=1;
        
        for ee in edgelist:
            e1=ee[0]+start_index; e2=ee[1]+start_index;
            s+=str(e1)+' '+str(e2)+' ';  soma_edges.append(num_of_edges); num_of_edges+=1;          
        
        s=s[:-1]
        f.write(s)
        f.write('</edges>\n')
        ##########################################################################################
        f.write('<triangles>')
        s=''
        mapping_points=[];
        cur=0; nxt=cur+1;
        for ky0 in cont.keys():
            klst=list(cont[ky0].keys())
            for j in range(len(klst)):
                npts=len(cont[ky0][klst[j]])
                for i in range(npts-1):
                    if j != (len(klst)-1):
                        ## collect face points here!!!!!!
                        if G.nodes[klst[j+1]]['t']==1:
                            soma_faces.append(num_of_faces);
                            soma_faces.append(num_of_faces+1);
                        else:
                            neurite_faces.append(num_of_faces);
                            neurite_faces.append(num_of_faces+1);
                        s+=str(cur+ncirpts)+' '+str(cur+1)+' '+str(cur)+' '; num_of_faces+=1;
                        s+=str(cur)+' '+str(cur+ncirpts-1)+' '+str(cur+ncirpts)+' '; num_of_faces+=1;
                        if G.nodes[klst[j]]['t']!=1:
                            mapping_points.append(G.nodes[klst[j]]['pos'])
                            mapping_points.append(G.nodes[klst[j+1]]['pos'])
                    cur=nxt; nxt+=1;
                if j != (len(klst)-1):
                    if G.nodes[klst[j+1]]['t']==1:
                        soma_faces.append(num_of_faces);
                        soma_faces.append(num_of_faces+1);
                    else:
                        neurite_faces.append(num_of_faces);
                        neurite_faces.append(num_of_faces+1);
                    s+= str(nxt-1+ncirpts)+' '+str(nxt)+' '+str(nxt-1)+' '; num_of_faces+=1;
                    s+= str(nxt-1)+' '+str(nxt+ncirpts-2)+' '+str(nxt-1+ncirpts)+' '; num_of_faces+=1;
                    if G.nodes[klst[j]]['t']!=1:
                        mapping_points.append(G.nodes[klst[j]]['pos'])
                        mapping_points.append(G.nodes[klst[j+1]]['pos'])
                cur=nxt; nxt+=1;
        
        for ff in facelist:
            f1=ff[0]+start_index; f2=ff[1]+start_index; f3=ff[2]+start_index;
            s+=str(f1)+' '+str(f2)+' '+str(f3)+' ';  soma_faces.append(num_of_faces); num_of_faces+=1;
         
        s=s[:-1]
        f.write(s)
        f.write('</triangles>\n')
        #f.write('<vertex_attachment name="npNormals" type="vector3" passOn="1" global="1">')
        #f.write('</vertex_attachment>\n')
        #f.write('<vertex_attachment name="npSurfParams" type="NeuriteProjectorSurfaceParams" passOn="1" global="1">')
        #f.write('</vertex_attachment>\n')
        f.write('<vertex_attachment name="npMapping" type="Mapping" passOn="1" global="1">')
        s=''
        cnt=0
        for mppos in mapping_points:
            s+=str(mppos[0])+' '+str(mppos[1])+' '+str(mppos[2])+' '
            if cnt != 0 and cnt % 2 ==1:
                s+=str(0.004)+' '
            cnt+=1
        
        # Soma face mapping        
        for ff in facelist:
            s+=str(soma_center[0])+' '+str(soma_center[1])+' '+str(soma_center[2])+' '+str(0.0)+' '
            s+=str(0.004)+' '
            s+=str(soma_center[0])+' '+str(soma_center[1])+' '+str(soma_center[2])+' '+str(0.0)+' '
            s+=str(0.004)+' '
        s=s[:-1]
        f.write(s)
        f.write('</vertex_attachment>\n')
        f.write('<subset_handler name="defSH">\n')
        f.write('<subset name="Neurites" color="0.588235 0.588235 1 1" state="0">\n')
        f.write('<vertices>')
        s=''
        for i in neurite_points: #range(num_of_vertices):
            s+=str(i)+' '
        s=s[:-1]
        f.write(s)
        f.write('</vertices>\n')
        f.write('<edges>')
        s=''
        for i in neurite_edges: #range(1,num_of_edges):
            s+=str(i)+' '
        s=s[:-1]
        f.write(s)
        f.write('</edges>\n')
        f.write('<faces>')
        s=''
        for i in neurite_faces: #range(1,num_of_faces):
            s+=str(i)+' '
        s=s[:-1]
        f.write(s)
        f.write('</faces>\n')
        f.write('</subset>\n')
        f.write('<subset name="Soma" color="1 0 0" state="0">\n')
        f.write('<vertices>')
        s=''
        for i in soma_points: #range(num_of_vertices):
            s+=str(i)+' '
        s=s[:-1]
        f.write(s)
        f.write('</vertices>\n')
        f.write('<edges>')
        s=''
        for i in soma_edges: #range(num_of_vertices):
            s+=str(i)+' '
        s=s[:-1]
        f.write(s)
        f.write('</edges>\n')
        f.write('<faces>')
        s=''
        for i in soma_faces: #range(num_of_vertices):
            s+=str(i)+' '
        s=s[:-1]
        f.write(s)
        f.write('</faces>\n')
        f.write('</subset>\n')
        f.write('</subset_handler>\n')
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n');
        f.write('</projection_handler>\n');
        f.write('</grid>');


def write_ugx_sphere(filename,n,cirpts,rad,soma_center):
    points=[];
    num_vertices=n*cirpts+2;

    cont={}
    for i in range(1,n+1):
        cont[i]=contour_circles(i*np.pi/(n+1),rad,cirpts,soma_center)
        cont[i]=cont[i][:-1]
    cont[n+1]=[tuple([soma_center[0],soma_center[1],soma_center[2]+rad])];
    cont[n+2]=[tuple([soma_center[0],soma_center[1],soma_center[2]-rad])];

    npoints=[]; nedges=[]; nfaces=[];
    edgelist=[]; facelist=[];
    num_points=0; num_edges=0; num_faces=0;

    with open(filename,'w') as f:        
        f.write('<?xml version="1.0" encoding="utf-8"?>\n');
        f.write('<grid name="defGrid">\n');
        f.write('<vertices coords="3">');
        line=''
        for ky in cont.keys():
            for pts in cont[ky]:
                line+=str(pts[0])+' '+str(pts[1])+' '+str(pts[2])+' '
                npoints.append(num_points); num_points+=1;
                points.append(pts)
                
        line=line[:-1]
        f.write(line)
        f.write('</vertices>\n');
        f.write('<edges>')
        cur=0; nxt=1;
        line='';

        for i in range(n):
            init=cur
            for j in range(cirpts-1):
                line+=str(cur)+' '+str(nxt)+' ';nedges.append(num_edges); num_edges+=1
                edgelist.append(tuple([cur,nxt]))
                cur=nxt; nxt+=1
            line+=str(cur)+' '+str(init)+ ' ';nedges.append(num_edges); num_edges+=1
            edgelist.append(tuple([cur,init]))
            cur+=1; nxt=cur+1   
        cur=0;

        for i in range(n-1):
            for j in range(cirpts):
                line+=str(cur)+' '+str(cur+cirpts)+' ';nedges.append(num_edges); num_edges+=1
                edgelist.append(tuple([cur,cur+cirpts]))
                cur+=1

        cur=0; cnt=1;
        for i in range(n-1):
            init=cur
            for j in range(cirpts):
                if cur % cirpts == 0:
                    ss=str(cur)+' '+str((cnt+1)*cirpts-1)+' ';nedges.append(num_edges); num_edges+=1
                    edgelist.append(tuple([cur,(cnt+1)*cirpts-1]))
                    cnt+=1
                else:
                    ss=str(cur)+' '+str(cur+cirpts-1)+' ';nedges.append(num_edges); num_edges+=1
                    edgelist.append(tuple([cur,cur+cirpts-1]))
                cur+=1; line+=ss;

        # bottom point edges
        for i in range(cirpts):
            line+=str(cur+i)+' '+str(num_vertices-1)+' ';nedges.append(num_edges); num_edges+=1
            edgelist.append(tuple([cur+i,num_vertices-1]))
        # top point edges
        for i in range(cirpts):
            line+=str(i)+' '+str(num_vertices-2)+' ';nedges.append(num_edges); num_edges+=1
            edgelist.append(tuple([i,num_vertices-2]))
        line=line[:-1]
        f.write(line)
        f.write('</edges>\n')
        f.write('<triangles>')
        cur=0; nxt=1;
        cnt=0;
        line='';
        for i in range(n-1):
            init=cur
            for j in range(cirpts):
                if j<(cirpts):
                    if nxt % cirpts==0:
                        line+=str(cur+cirpts)+' '+str(cnt*cirpts)+' '+str(cur)+' ';nfaces.append(num_faces); num_faces+=1;
                        facelist.append(tuple([cur+cirpts,cnt*cirpts,cur]))
                        cnt+=1
                    else:
                        line+=str(cur+cirpts)+' '+str(nxt)+' '+str(cur)+' ';nfaces.append(num_faces); num_faces+=1;
                        facelist.append(tuple([cur+cirpts,nxt,cur]))
                    cur=nxt; nxt+=1

        cur=0; nxt=1; cnt=1;
        for i in range(n-1):
            init=cur
            for j in range(cirpts):
                if cur % cirpts != 0:
                    line+=str(cur)+' '+str(cur+cirpts-1)+' '+str(cur+cirpts)+' ';nfaces.append(num_faces); num_faces+=1;
                    facelist.append(tuple([cur,cur+cirpts-1,cur+cirpts]))
                else:
                    line+=str(cur)+' '+str((cnt+1)*cirpts-1)+' '+str(cur+cirpts)+' ';nfaces.append(num_faces); num_faces+=1;
                    facelist.append(tuple([cur,(cnt+1)*cirpts-1,cur+cirpts]))
                    cnt+=1
                cur=nxt; nxt+=1

        # bottom point faces
        for i in range(cirpts):
            init=cur
            if i<(cirpts-1):
                line+=str(num_vertices-1)+' '+str(cur+i+1)+' '+str(cur+i)+' ';nfaces.append(num_faces); num_faces+=1;
                facelist.append(tuple([num_vertices-1,cur+i+1,cur+i]))
            else:
                line+=str(num_vertices-1)+' '+str(init)+' '+str(cur+i)+' ';nfaces.append(num_faces); num_faces+=1;
                facelist.append(tuple([num_vertices-1,init,cur+i]))

        # top point faces
        for i in range(cirpts):
            if i<(cirpts-1):
                line+=str(i)+' '+str(i+1)+' '+str(num_vertices-2)+' ';nfaces.append(num_faces); num_faces+=1;
                facelist.append(tuple([i,i+1,num_vertices-2]))
            else:
                line+=str(i)+' '+str(0)+' '+str(num_vertices-2)+' ';nfaces.append(num_faces); num_faces+=1;
                facelist.append(tuple([i,0,num_vertices-2]))

        line=line[:-1]
        f.write(line)
        f.write('</triangles>\n')
        f.write('<subset_handler name="defSH">\n')
        f.write('<subset name="Soma" color="1 0 0" state="0">\n')
        f.write('<vertices>')
        ss=''
        for i in npoints:
            ss+=str(i)+' '
        ss=ss[:-1]
        f.write(ss)
        f.write('</vertices>\n')
        f.write('<edges>')
        ss=''
        for i in nedges:
            ss+=str(i)+' '
        ss=ss[:-1]
        f.write(ss)
        f.write('</edges>\n')
        f.write('<faces>')
        ss=''
        for i in nfaces:
            ss+=str(i)+' '
        ss=ss[:-1]
        f.write(ss)
        f.write('</faces>\n')
        f.write('</subset>\n')
        f.write('</subset_handler>\n')
        f.write('<projection_handler name="defPH" subset_handler="0">\n');
        f.write('<default type="default">0 0</default>\n');
        f.write('</projection_handler>\n');
        f.write('</grid>');
        
        return points, edgelist, facelist