import argparse
import os
import shutil
import sys
import copy
print("Python version: ",sys.version)
print("Version info:   ",sys.version_info)
print('Argparse:       ',argparse.__version__)

sys.path.insert(0, os.getcwd()+'/src/')

import neuronswc as ns
import neuronmeshgenerator as nmg

parser = argparse.ArgumentParser(description='''This program will generate .swc refinements, .ugx 1d refinements, and .ugx surface meshes and zip all files into a .vrn file,
                                             usage: python3 NEW_generate_meshes.py -n 6 -c 6 -d 32.0 -p 10 -q 16 -i cells/<cellname>''')
parser.add_argument('-n', '--numrefine',help="Number of Refinements", required=False,default=6)
parser.add_argument('-d', '--startdx', help="Largest DeltaX",required=False,default=32)
parser.add_argument('-c', '--numcontpts',help="Number of Contour points",required=False,default=12)
parser.add_argument('-p', '--spherecontours',help="Number of sphere contours",required=False,default=10)
parser.add_argument('-q', '--spherepoints', help="Number of points per sphere contour",required=False,default=16)
parser.add_argument('-i', '--input',help="The input .swc file", required=True)
args = parser.parse_args()

print('Number of refinements:               ',args.numrefine)
print('Start DX:                            ',args.startdx)
print('Number of contour points:            ',args.numcontpts)
print('Number of Sphere contours:           ',args.spherecontours)
print('Number of Sphere points per contour: ',args.spherepoints)
print('Input file:                          ',args.input)

n=int(args.numrefine)
nsphere_contours=int(args.spherecontours);
nsphere_contour_pts=int(args.spherepoints);

CELL_FILE=args.input
npts=int(args.numcontpts);
start_dx=float(args.startdx);

MESHFOLDER=nmg.make_meshes(CELL_FILE,start_dx,n,npts,nsphere_contours,nsphere_contour_pts)
nmg.write_vrn(MESHFOLDER)
print('Done making VRN!...',MESHFOLDER,' is the folder')
print('Removing MESHFOLDER...')
shutil.rmtree(MESHFOLDER)
