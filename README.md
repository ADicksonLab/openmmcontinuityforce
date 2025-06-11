OpenMM Continuity Force
=====================

This is a plugin for [OpenMM](https://openmm.org) that implements a custom many-particle
force called a "ContinuityForce".  This force aims to ensure that all particles in a set
remain continuously connected in space.  More precisely, if we imagine the particles as
nodes of a network and draw an edge between particles that are closer than a cutoff distance ($d$),
then the nodes would all be in the same component of the network.  In other words, there exists
a path of some length between any two nodes $i$ and $j$.

ContinuityForce achieves this by classifying particles into components on the fly and adding
an attractive force between pairs of particles belonging to different components.  The particles
are chosen each dynamics step as the closest possible pairs between the two components.  The
attractive force is of the form $E(r)=k(r-d)^2$, where r is the separation between
the particles and $d$ is a user-defined cutoff distance.

This was originally designed to aid FlexibleTopology simulations ([paper](https://pubs.acs.org/doi/full/10.1021/acs.jctc.3c00409)).  This plugin was made using the ["OpenMMExamplePlugin" template](https://github.com/openmm/openmmexampleplugin).

Preparing the Environment
===================

Use a CONDA environment for easy installation. Before installing this plugin,
you should create this environment and install openmm:
```
conda install -c conda-forge openmm=8.1 cudatoolkit=11.8
```
Note: explicitly installing the cudatoolkit is now necessary for newer openmm versions.  At this point you
can also check if openmm can communicate with your GPU:
```
python -m openmm.testInstallation
```
where you should see that CUDA successfully computes forces.

You will also need swig:
```
conda install swig
```
And if you need to install cmake this can also be done using CONDA:
```
conda install cmake
```

Building/Installing The Plugin
===================

To build and install `openmmcontinuityforce`, follow these steps:

1. Clone this repository to your own machine.

2. Create a directory in which to build the plugin, e.g. `openmmcontinuityforce/build`.

3. Change to the build directory, and run cmake: `cmake ..` and open ccmake: `ccmake ..`

4. Set OPENMM_DIR to point to the directory where OpenMM is installed.  This is needed to locate
the OpenMM header files and libraries.  If you are using CONDA, this directory will be something
like: `/your/path/to/anaconda/envs/env_name`. 

5. Set CMAKE_INSTALL_PREFIX to the same directory.

6. If you plan to build the CUDA platform, make sure that CUDA_TOOLKIT_ROOT_DIR is set correctly
and that BUILD_CUDA_LIB is selected.

7. Press "Configure" again if necessary, then press "Generate".

8. Type `make install`, followed by `make PythonInstall`.  Note: there will likely be many warnings for PythonInstall related to `OpenMMSwigHeaders.i`.  These can be safely ignored and as long as the build reaches 100%.


Test Cases
==========

To run all the test cases build the "test" target by typing `make test`.

Accessing ContinuityForce in Python
==========

There are different ways to make the plugin accessible in a python script.
I recommend including the following line in your python script:
```
omm.Platform.loadPluginsFromDirectory('/path/to/your/anaconda/pkgs/openmm-7.7.0-py39h9717219_0/lib/plugins')
```
where again you will likely need to change the version and build of OpenMM to match that
on your machine.

You can then use the plugin in your Python scripts:
```
from openmmcontinuityforce import ContinuityForce
force = ContinuityForce()
```

