import renderneuron as rn
import parallelframe as pf
import argparse
import os

parser = argparse.ArgumentParser(description='This program will generate .swc refinements')
parser.add_argument('-n', '--numrefine',help="Number of Refinements", required=True)
parser.add_argument('-c', '--numcontpts',help="Number of Contour points",required=True)
parser.add_argument('-i', '--input',help="The input .swc file", required=True)
parser.add_argument('-o', '--output',help="The output folder name", required=True)
parser.add_argument('--spline', action='store_true',help="Use splines")
args = parser.parse_args()

print('Number of refinements:    ',args.numrefine)
print('Number of contour points: ', args.numcontpts)
print('Input file:               ',args.input)
print('Output file:              ',args.output)

n=args.numrefine

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
        rn.save_to_swc(GS,OUTPUT_DIR+'/refinement'+str(i)+'.swc')
        cont=pf.get_pft_frames(GS,npts)
        outfilename=OUTPUT_DIR+'/mesh'+str(i)+'.ugx'
        pf.write_ugx(cont,npts,outfilename)
    print(' Done!')
else:
    GG=rn.refine_and_save(G,OUTPUT_DIR+'/',int(n))
    for i in range(len(GG)):
        H=GG[i]
        cont=pf.get_pft_frames(H,npts)
        outfilename=OUTPUT_DIR+'/mesh'+str(i)+'.ugx'
        pf.write_ugx(cont,npts,outfilename)
    trunks,T=rn.get_trunks(G)
    rn.save_to_swc(T,OUTPUT_DIR+'/trunks.swc')