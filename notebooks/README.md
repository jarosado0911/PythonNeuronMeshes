# Notebooks

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

## Installation
For this `Running_Yale_Neuron_Simulation` you will need Neuron 8.2.2, this is the version I used. For the the `Exploring_Functionality` notebook you will need the versions that are listed in the main `readme` of the project.

Here is the version of `ffmpeg` I used for making the movie:
```
ffmpeg version 4.2.7-0ubuntu0.1 Copyright (c) 2000-2022 the FFmpeg developers
  built with gcc 9 (Ubuntu 9.4.0-1ubuntu1~20.04.1)
```