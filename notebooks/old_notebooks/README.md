# Notebooks [OLD]

Contents
========

* [About](#about)
* [Files](#files)
* [Installation](#installation)

## About

This folder contains exploratory notebooks using the mesh generation modules and interfacing with Yale Neuron.
This link: [`Yale Neuron`](https://www.neuron.yale.edu/neuron/) is the main website for Yale's Neuron simulator. 
Neuron is a simulation environment for modeling individual and networks of neurons. It was primarily developed by Michael Hines, John W. Moore, and Ted Carnevale at Yale and Duke.

Neuron models individual neurons via the use of sections that are automatically subdivided into individual compartments, instead of requiring the user to manually create compartments. The primary scripting language is hoc but a Python interface is also available. 

There is also an example video I made called `movie.mp4` showing electrical signal propagation from a current clamp simulation inside the `Running_Yale_Neuron_Simulation` notebook.

## Files

- The `Exploring_Functionality` goes through the process of reading an `.swc` file, making direct refinements, making spline refinements, and generating `.ugx` surface geometries

- The `Running_Yale_Neuron_Simulation` goes through loading a `.swc` file, assigning Hodgkin Huxley mechanisms to the membrane, collecting cell information into a dictionary, setting up a current clamp simulation, running the simulation, and then making a video using a sequence of images. The notebook will make a directory called `images` where the images for the video are stored, the `images` directory is not part of this repo. Below is an example video
<p align="center">
  <img src="./../img/yaleneuron.gif" alt="Size Limit CLI" width="500">
</p>

## Installation
For this `Running_Yale_Neuron_Simulation` you will need Neuron 8.2.2, this is the version I used. For the the `Exploring_Functionality` notebook you will need the versions that are listed in the main `readme` of the project.

Here is the version of `ffmpeg` I used for making the movie:
```
ffmpeg version 4.2.7-0ubuntu0.1 Copyright (c) 2000-2022 the FFmpeg developers
  built with gcc 9 (Ubuntu 9.4.0-1ubuntu1~20.04.1)
```
Here is another video:

https://github.com/jarosado0911/PythonNeuronMeshes/assets/57325500/8202f368-b1c0-4aa8-a41d-b6f44b1be511
