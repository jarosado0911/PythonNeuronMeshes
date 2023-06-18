# Source Files

Contents
========

 * [Files](#files)
 * [Render Neuron](#renderneuron)
 * [Parallel Frame](#parallelframe)
 
 ### Files
 This is the source directory which contains the following two python files:
+ `renderneuron.py` contains the functions for making and saving simple refinements, and spline refinements. The file also contains some plotting features
+ `parallelframe.py` contains the functions which uses [`Parallel Transport Frames`](https://legacy.cs.indiana.edu/ftp/techreports/TR425.pdf) to generate contours along a path, and then connects the contours to make a mesh. There is a function for writing the `.ugx` file which contains the mesh with surface triangulations.

### Render Neuron

### Parallel Frame Transport