# PythonNeuronMeshes

Contents
========

 * [What?](#what)
 * [Installation](#installation)
 * [Usage](#usage)
 
### What?
This project takes a `.swc` file downloaded from [`NeuroMorpho.org`](https://neuromorpho.org/) and generates:

+ 1D geometry refinements in .swc format
+ 3D surface mesh geometries in .ugx format

The code will either use a spline interpolation for regularizing the geometry and producing refinements or only does pure refinements of the original `.swc` file

### Installation
The project uses Python in particular I provide some version information that is currently used.
I recommend using Jupyter-Notebooks and downloading [`Promesh4`](https://promesh3d.com/) for view the geometries.
Other than that there is nothing else to install.
```
Python version:  3.8.10 (default, Mar 13 2023, 10:26:41)
[GCC 9.4.0]
Version info:    sys.version_info(major=3, minor=8, micro=10, releaselevel='final', serial=0)
Argparse:        1.1
networkx:   3.1
scipy:      1.10.1
numpy:      1.24.3
maplotlib:  3.7.1
```

### Usage
There is a file called `generate_meshes.py` if you execute in the commandline `python3 generate_mesh.py` you will receive the following output
```
usage: generate_meshes.py [-h] -n NUMREFINE -c NUMCONTPTS -i INPUT -o OUTPUT [--spline]
generate_meshes.py: error: the following arguments are required: -n/--numrefine, -c/--numcontpts, -i/--input, -o/--output
```
If you `python3 generate_meshes.py -h` you will get some more help information:
```
usage: generate_meshes.py [-h] -n NUMREFINE -c NUMCONTPTS -i INPUT -o OUTPUT [--spline]

This program will generate .swc refinements, usage: python3 generate_mesh.py -n 4 -c 6 -i cells/<cellname> -o
<outfoldername> --spline

optional arguments:
  -h, --help            show this help message and exit
  -n NUMREFINE, --numrefine NUMREFINE
                        Number of Refinements
  -c NUMCONTPTS, --numcontpts NUMCONTPTS
                        Number of Contour points
  -i INPUT, --input INPUT
                        The input .swc file
  -o OUTPUT, --output OUTPUT
                        The output folder name
  --spline              Use splines
```