import argparse
import os
import sys
import copy
print("Python version: ",sys.version)
print("Version info:   ",sys.version_info)
print('Argparse:       ',argparse.__version__,'\n')

sys.path.insert(0, os.getcwd()+'/src/old/')

import renderneuron as rn
import parallelframe as pf
import somasphere as sph
print('\n')

parser = argparse.ArgumentParser(description='''This program will generate .swc refinements, .ugx 1d refinements, and .ugx surface meshes and zip all files into a .vrn file,
                                             usage: python3 generate_mesh.py -n 4 -c 6 -s 1.10 -p 10 -q 16 -i cells/<cellname> -o <outfoldername> --spline''')
parser.add_argument('-n', '--numrefine',help="Number of Refinements", required=True)
parser.add_argument('-c', '--numcontpts',help="Number of Contour points",required=True)
parser.add_argument('-s', '--scalesoma', help="Scale the soma", required=True)
parser.add_argument('-p', '--spherecontours',help="Number of sphere contours",required=True)
parser.add_argument('-q', '--spherepoints', help="Number of points per sphere contour",required=True)
parser.add_argument('-i', '--input',help="The input .swc file", required=True)
parser.add_argument('-o', '--output',help="The output folder name", required=True)
parser.add_argument('--spline', action='store_true',help="Use splines")
args = parser.parse_args()

print('Number of refinements:               ',args.numrefine)
print('Number of contour points:            ',args.numcontpts)
print('Soma Scale:                          ',args.scalesoma)
print('Number of Sphere contours:           ',args.spherecontours)
print('Number of Sphere points per contour: ',args.spherepoints)
print('Input file:                          ',args.input)
print('Output file:                         ',args.output)
print('Use splines:                         ',args.spline)

n=int(args.numrefine)
scale_soma=float(args.scalesoma);
print(scale_soma)
nsphere_contours=int(args.spherecontours);
nsphere_contour_pts=int(args.spherepoints);

OUTPUT_DIR=args.output
if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)

CELL_FILE=args.input
lines=rn.read_swc(CELL_FILE)
G=rn.get_graph(lines)
npts=int(args.numcontpts);

if (args.spline):
    pw=list(range(7,7-int(n),-1))
    dx=[pow(2,p) for p in pw]
    print('DX = ', dx)
    print('refining round...',end=' ')
    for i in range(len(dx)):
        GS=rn.spline_neuron(G,dx[i])
        print(str(dx[i]),' ',end='')
        GS=rn.remove_soma_line(GS)
        rn.save_to_swc(GS,OUTPUT_DIR+'/refinement'+str(i)+'.swc')
        pf.write_1d_ugx(GS,OUTPUT_DIR+'/cell_1d_ref_'+str(i)+'.ugx')
        cont=pf.get_pft_frames(GS,npts)
        outfilename=OUTPUT_DIR+'/cell_soma_sphere_3d_tris_x1_ref_'+str(i)+'.ugx'
        # write the .ugx to file
        sph.write_ugx(cont,GS,npts,outfilename,scale_soma,nsphere_contours,nsphere_contour_pts)        
    print(' Done!')
    rn.write_vrn(OUTPUT_DIR)
    
    G0=rn.remove_soma_line(G)
    trunks,T=rn.get_trunks(G0)
    rn.save_to_swc(T,OUTPUT_DIR+'/trunks.swc')
    
else:
    GG=rn.refine_and_save(G,OUTPUT_DIR+'/',int(n))
    for i in range(len(GG)):
        H=copy.deepcopy(GG[i])
        H=rn.remove_soma_line(H)
        cont=pf.get_pft_frames(H,npts)
        outfilename=OUTPUT_DIR+'/cell_3d_tris_x1_ref_'+str(i)+'.ugx'
        
        sph.write_ugx(cont,H,npts,outfilename,scale_soma,nsphere_contours,nsphere_contour_pts)
        outfilename=OUTPUT_DIR+'/mesh_subsets'+str(i)+'.ugx'
        pf.write_ugx_subsets(cont,H,npts,outfilename)
        pf.write_1d_ugx(H,OUTPUT_DIR+'/cell_1d_ref_'+str(i)+'.ugx')
    trunks,T=rn.get_trunks(G)
    rn.save_to_swc(T,OUTPUT_DIR+'/trunks.swc')
    rn.write_vrn(OUTPUT_DIR)