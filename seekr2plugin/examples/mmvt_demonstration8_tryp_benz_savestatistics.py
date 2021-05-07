from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from seekr2plugin import MmvtLangevinIntegrator, vectori, vectord
import seekr2plugin
from time import time
import numpy as np
import os

lig_indices = [3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229]
#3221 3222 3223 3224 3225 3226 3227 3228 3229
rec_indices = [2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
# 2478 2489 2499 2535 2718 2745 2769 2787 2794 2867 2926

box_vector = Quantity(
    [[61.23940982865766, 0.0, 0.0], 
     [-20.413134962473007, 57.73706986993814, 0.0], 
     [-20.413134962473007, -28.868531440986094, 50.00177126469543]], 
                      unit=angstrom)

prmtop = AmberPrmtopFile('tryp_ben.prmtop')
inpcrd = AmberInpcrdFile('tryp_ben.inpcrd')

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

mypdb = PDBFile('tryp_ben_new.pdb')


myforce1 = CustomCentroidBondForce(2, "-k*(distance(g1,g2)^2-radius^2)")
mygroup1a = myforce1.addGroup(rec_indices)
mygroup2a = myforce1.addGroup(lig_indices)
myforce1.setForceGroup(1)
myforce1.addPerBondParameter('k')
myforce1.addPerBondParameter('radius')
myforce1.addBond([mygroup1a, mygroup2a], [1.0e-9*kilojoules_per_mole, 11.0*angstroms])
forcenum1 = system.addForce(myforce1)

myforce2 = CustomCentroidBondForce(2, "k*(distance(g1,g2)^2-radius^2)")
mygroup1b = myforce2.addGroup(rec_indices)
mygroup2b = myforce2.addGroup(lig_indices)
myforce2.setForceGroup(2)
myforce2.addPerBondParameter('k')
myforce2.addPerBondParameter('radius')
myforce2.addBond([mygroup1b, mygroup2b], [1.0e-9*kilojoules_per_mole, 13.0*angstroms])
forcenum2 = system.addForce(myforce2)

integrator = MmvtLangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, "tryp_test_filename.txt")
integrator.addMilestoneGroup(1)
integrator.addMilestoneGroup(2)
integrator.setSaveStatisticsFileName("statistics.txt")


platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(mypdb.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)

simulation.context.setPeriodicBoxVectors(*box_vector)
    
<<<<<<< HEAD:plugin/examples/demonstration8_tryp_benz_savestatistics.py
#simulation.minimizeEnergy()
=======
simulation.minimizeEnergy()
>>>>>>> 42eaa5868c5b72c82b4482d5b0a4c1b903a86bfe:examples/demonstration3_tryp_benz_spherical.py
#simulation.reporters.append(PDBReporter('/home/lvotapka/doc/mmvt_test/tryp_ben_output.pdb', 20))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.step(100000)

#state = simulation.context.getState()
#print(state.getPeriodicBoxVectors())
