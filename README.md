seekr2_openmm_plugin
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/seekr2_openmm_plugin/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/seekr2_openmm_plugin/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/seekr2_openmm_plugin/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/seekr2_openmm_plugin/branch/master)


An OpenMM plugin for SEEKR2.

## TABLE OF CONTENTS:
1. INTRODUCTION
2. INSTALLATION OF OPENMM
3. INSTALLATION OF SEEKR2_OPENMM_PLUGIN
4. MMVT LANGEVIN INTEGRATOR
5. ELBER LANGEVIN INTEGRATOR:
6. MMVT/ELBER SURFACE DEFINITIONS
7. SAMPLES SCRIPTS
8. THE CROSSINGS FILE
9. CROSSING STATE ANALYSIS


## INTRODUCTION:
The SEEKR2 OpenMM plugin is a plugin for the molecular dynamics (MD) simulation 
software package OpenMM. Technically, the plugin defines new Integrator 
objects in OpenMM, which allow the user to define a series of surfaces using 
Force objects, and monitors crossing events on these surfaces. Normally, the 
integrators act just like an the LangevinIntegrator or other common integrators
that ships with OpenMM. But,for instance, in MMVT when a surface crossing event
is detected, the positions, velocities, and forces of the previous timestep are
restored, the velocities are reversed, the crossing event is logged, and 
(optionally) the locations and velocities of the atoms during the crossing are 
recorded. For Elber milestoning, crossing events of the various milestones
are recorded by the integrators.

## INSTALLATION OF OPENMM

The SEEKR2 OpenMM plugin requires OpenMM to run.

The instructions for installing OpenMM from source can be found on OpenMM's
main manual: http://docs.openmm.org/6.1.0/userguide/library.html

One may also follow the instructions below for installing OpenMM from source.

Download and install Conda for Python3 (We have tested python3.5, and 
python3.6). Conda can be downloaded from: 
https://conda.io/en/latest/miniconda.html

Install numpy, scipy, netcdf4, and mpi4py:

```
$ conda install numpy scipy netcdf4 mpi4py
```

Make sure git is installed. Make sure 'curses' is installed.

```
$ sudo apt-get install cmake-curses-gui
```

Make sure 'doxygen' is installed.

```
$ conda install -c conda-forge doxygen
```

You will also need 'swig'

```
$ conda install swig
```

Clone OpenMM and cd into OpenMM directory, then perform necessary build steps.

```
$ git clone https://github.com/openmm/openmm.git
$ cd openmm
$ mkdir build
$ cd build
$ ccmake ..
```

The ccmake gui should come up. Press 'c' and then 't'

You should modify the following variables:

CMAKE_INSTALL_PREFIX: change to a local directory that exists (example: 
/home/username/bin/openmm). If such a directory doesn't exist, then make one.

Press 'c'. When the configuration is successful, type 'g' to generate. Then 
ccmake should close on its own.

Type:

```
$ make
$ make install
$ make PythonInstall
$ make test
```

If the PythonInstall step fails, then install cython

```
$ pip install --upgrade cython
```

You can test your OpenMM installation by opening python and importing OpenMM:

```
>>> from simtk import openmm
```

If the import completes without error, the installation should have been 
successful.



## SEEKR2 PLUGIN INSTALLATION:

Clone the SEEKR2 Plugin and cd into the directory, then perform necessary 
build steps.

```
$ git clone https://github.com/seekrcentral/seekr2_openmm_plugin.git
$ cd seekr2_openmm_plugin/seekr2plugin
$ mkdir build
$ cd build
$ ccmake ..
```

The ccmake gui should come up. Press 'c' and then 't'

You should modify the following variables:

CMAKE_INSTALL_PREFIX: change to the same directory as the variable used in the 
OpenMM installation above.

OPENMM_DIR: Use the same directory that you used for CMAKE_INSTALL_PREFIX for 
OpenMM above.

Press 'c'. When the configuration is successful, type 'g' to generate. Then 
ccmake should close on its own.

Type:

```
$ make
$ make install
$ make PythonInstall
$ make test
```

You can test your installation by opening Python and typing:

```
>>> from seekr2plugin import MmvtLangevinIntegrator
```

If the import completes without error, the installation should have been 
successful.



## MMVT LANGEVIN INTEGRATOR:

Using SEEKR2 first requires one to define an integrator object, and we will 
start by describine the MmvtLangevinIntegrator object. An initialization of 
the object is called using the following syntax:

MmvtLangevinIntegrator(temperature, frictionCoefficient, timeStep, 
                       crossingsFileName)

Where temperature, frictionCoefficient, and timeStep are the same parameters 
that would be passed to OpenMM's LangevinIntegrator object.

The argument for crossingsFileName is a string that defines the location to 
write the crossing events. A description of this file output can be found in 
the "CROSSINGS FILE" section lower in this README.

The integrator object has a number of methods, the most important of which are:

 - addMilestoneGroup(milestoneGroup): allows one to define which force groups 
   SEEKR2 should monitor for crossing events. This method must be called for 
   each MMVT surface that is desired.
 - setSaveStateFileName(fileName): A directory-filename prefix that will be 
   used if the user chooses to write position/velocity states upon a milestone 
   crossing event. If this method is never called or a blank string is passed, 
   then the state will not be saved upon crossing. A method for extracting 
   atomic positions from these saved states is provided in the section below 
   labeled "CROSSING STATE ANALYSIS".
 - setSaveStatisticsFileName(fileName): The argument is a string defining the 
   location to write MMVT statistics directly.
 - setBounceCounter(counter): the argument is an integer that will define the
   starting number of bounces. This is used when restarting MMVT simulations.

## ELBER LANGEVIN INTEGRATOR:

An integrator designed to perform Elber milestoning is also provided as an
alternative to MMVT milestoning. Like the MMVT integrator, the Elber
integrator is initialized according to the following syntax:

ElberLangevinIntegrator(temperature, frictionCoefficient, timeStep, 
                       crossingsFileName)

Where temperature, frictionCoefficient, and timeStep are the same parameters 
that would be passed to OpenMM's LangevinIntegrator object.

The argument for crossingsFileName is a string that defines the location to 
write the crossing events. A description of this file output can be found in 
the "CROSSINGS FILE" section lower in this README.

The integrator object has a number of methods, the most important of which are:

 - addSrcMilestoneGroup(milestoneGroup): allows one to define which force groups 
   SEEKR2 should monitor for source milestone crossing events. This method must
   be called for each Elber source milestone surface that is desired. The
   source milestone is a surface that, if crossed in the reversal stage, will
   cause the system to be rejected from the first hitting point distribution 
   (FHPD).
 - addDestMilestoneGroup(milestoneGroup): allows one to define which force 
   groups SEEKR2 should monitor for source milestone crossing events. This 
   method must be called for each Elber destination milestone surface that is 
   desired. The destination milestone is a surface that, if crossed in the 
   forward stage, will be logged as a transition.
 - setSaveStateFileName(fileName): A directory-filename prefix that will be 
   used if the user chooses to write position/velocity states upon a milestone 
   crossing event. If this method is never called or a blank string is passed, 
   then the state will not be saved upon crossing. A method for extracting 
   atomic positions from these saved states is provided in the section below 
   labeled "CROSSING STATE ANALYSIS".
 - setCrossingCounter(counter): the argument is an integer that will define the
   starting number of crossings. This is used to reset each subsequent Elber
   reversal or forward trajectory.

## MMVT AND ELBER SURFACE DEFINITIONS:

The MMVT and Elber surfaces in SEEKR2 are defined using OpenMM Custom Force 
Objects. Before proceeding, the reader should become familiar with Forces, 
Custom Forces, Force groups, and Atom groups by carefully reading the 
following OpenMM docs:

http://docs.openmm.org/7.0.0/userguide/theory.html#standard-forces
http://docs.openmm.org/7.0.0/userguide/theory.html#custom-forces
http://docs.openmm.org/7.0.0/api-c++/generated/OpenMM.CustomNonbondedForce.html

In SEEKR2, a force group must be used to define the MMVT/Elber surface. Any 
OpenMM Force object can be used for this. The potential energy expression of the 
Force must be provided. Depending on the potential energy expression, the 
proper number of atom groups and Force parameters must be provided. Finally, 
SEEKR2 must be told to monitor this force group for crossing by using the 
MmvtLangevinIntegrator.addMilestoneGroup() function (or equivalent functions
in the ElberLangevinIntegrator).

SEEKR2 monitors the specified force groups. If the potential energy expression 
is negative or zero, then the system advances in time according to a typical 
Langevin integrator. However, once the potential energy of the Force is found 
to be positive, then the MMVT reflection mechanism is triggered, the crossing 
information is logged, and, optionally, the state is saved.

Therefore, the key step to defining these surfaces is to formulate an 
appropriate mathematical expression whose value is exactly zero all along the 
MMVT surface. Therefore, the MMVT surface is the zero isosurface of the 
provided expression.

A maximum of 32 force groups are allowed in OpenMM. Since group 0 should be 
reserved for non-SEEKR2 forces, that leaves 31 possible distinct SEEKR surfaces 
that can be defined per simulation.

## SAMPLE SCRIPTS:

A number of sample scripts are located in the examples/ directory to provide 
some examples of how to create MMVT/Elber surfaces of different types.

mmvt_demonstration1_argon_sphere.py: A system of 10 argon atoms interacting 
with the Lennard-Jones potential are bound within a spherical "Cartesian" 
surface ("Cartesian" surfaces are defined by an absolute location in space, 
and do not move with any particles in the system).

mmvt_demonstration2_argon_box.py: A system of 10 argon atoms are bound within 
a box defined by six planar "Cartesian" surfaces.

mmvt_demonstration3_tryp_benz_spherical.py: A system of trypsin and 
benzamidine, where the benzamidine is bound between two spheres at 11 and 13 
Angstroms, where the center of the spheres is defined as the center of mass of 
a number of key residues of the trypsin active site.

mmvt_demonstration4_trp_cage.py: The tryptophan cage peptide, where two MMVT 
surfaces are defined as the sum of the squared distances of three pairs of 
side chains that are interacting in the folded native structure. This is 
intended as a demonstration of potential protein-folding applications.

mmvt_demonstration5_pyrophosphate.py: Inorganic pyrophosphate in the active 
site of pyrophosphatase is bound between two "angular" MMVT surfaces.

mmvt_demonstration6_tryp_benz_savestate.py: This script demonstrates that 
crossing states can be saved for later analysis of atomic positions. Again, 
the trypsin/benzamidine system is used.

mmvt_demonstration7_tryp_benz_box.py: This script defines six spherical MMVT 
surfaces in the trypsin/benzamidine system that approximate a box-like shape 
near the location of the benzamidine. The centers of the spheres move along 
with the centers of mass of some hand-picked atoms on the trypsin.

elber_demonstration_1_tryp_benz_spherical.py: A system of trypsin and 
benzamidine, where the benzamidine is starting at a distance approximately
12 Angstroms, close to a source milestone at exactly 12 Angstroms. There are
also two spherical destination milestones at 11 and 13 Angstroms. For all
milestones, the center of the spheres is defined as the center of mass of 
a number of key residues of the trypsin active site.

make_pdb_from_mm_state.py: This script provides a sample for how position 
information can be extracted from the OpenMM state objects saved upon 
milestone crossing. This script is further discussed in the 
"CROSSING STATE ANALYSIS" section below.


## THE CROSSINGS FILE:

When the file defined by the crossingsFileName argument is opened, it contains 
two columns. The first column defines the force group representing the 
milestone surface that was crossed, and the second column is the time of 
crossing (in ps) since the start of the simulation. Notice that each of the 
crossing events are sequential. If two or more milestones are crossed in the 
same timestep, both/all crossings are logged as happening at the same time.

This file contains all the information necessary to recreate the crossing 
probabilities and times for later MMVT kinetics/thermodynamic analysis.


## CROSSING STATE ANALYSIS:

When a file prefix is provided to the setSaveStateFileName() method, then a 
new file will be created upon each crossing that will be named according to 
the provided prefix, but with the crossing surface force group and timestep 
appended to the name. These state files are written in XML, and contain 
plaintext information about atomic positions and velocities, as well as some 
other system information.

The make_pdb_from_mm_state.py script is provided as an example script for
extracting the positions from the state files. Note that the parmed package 
is needed to run this script.


### Copyright

Copyright (c) 2021, Lane Votapka


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
