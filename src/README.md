# Source Files

Contents
========

 * [Files](#files)
 * [Render Neuron](#renderneuron)
 * [Parallel Frame](#parallelframe)
 
 ### Files (NEED to UPDATE description)
 This is the source directory which contains the following two python files:
+ `neuronswc.py` contains the functions for making and saving simple refinements, and spline refinements.
+ `neuronmeshgenerator.py` contains the functions which uses [`Parallel Transport Frames`](https://legacy.cs.indiana.edu/ftp/techreports/TR425.pdf) to generate contours along a path, and then connects the contours to make a mesh. There is a function for writing the `.ugx` file which contains the mesh with surface triangulations.