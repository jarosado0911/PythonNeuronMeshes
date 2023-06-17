# PythonNeuronMeshes

Contents
========

 * [What?](#what)
 * [Installation](#installation)
 * [Usage](#usage)
 
### What?
This project takes a [`.swc`] file downloaded from [`NeuroMorpho.org`](https://neuromorpho.org/) and generates:

+ 1D geometry refinements in .swc format
+ 3D surface mesh geometries in .ugx format

The code will either use a spline interpolation for regularizing the geometry and producing refinements or only does pure refinements of the original [`.swc`] file

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