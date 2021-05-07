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
rec1_indices = [2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
# 2478 2489 2499 2535 2718 2745 2769 2787 2794 2867 2926
rec2_indices = [1115, 1129, 1140, 1154, 1168]
# 1115 1129 1140 1154 1168
rec3_indices = [566, 576, 597, 619, 630]
# 566 576 597 619 630

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
mygroup1a = myforce1.addGroup(rec1_indices)
mygroup2a = myforce1.addGroup(lig_indices)
myforce1.setForceGroup(1)
myforce1.addPerBondParameter('k')
myforce1.addPerBondParameter('radius')
myforce1.addBond([mygroup1a, mygroup2a], [1.0e-9*kilojoules_per_mole, 11.0*angstroms])
forcenum1 = system.addForce(myforce1)

myforce2 = CustomCentroidBondForce(2, "k*(distance(g1,g2)^2-radius^2)")
mygroup1b = myforce2.addGroup(rec1_indices)
mygroup2b = myforce2.addGroup(lig_indices)
myforce2.setForceGroup(2)
myforce2.addPerBondParameter('k')
myforce2.addPerBondParameter('radius')
myforce2.addBond([mygroup1b, mygroup2b], [1.0e-9*kilojoules_per_mole, 13.0*angstroms])
forcenum2 = system.addForce(myforce2)

myforce3 = CustomCentroidBondForce(2, "-k*(distance(g1,g2)^2-radius^2)")
mygroup1c = myforce3.addGroup(rec2_indices)
mygroup2c = myforce3.addGroup(lig_indices)
myforce3.setForceGroup(3)
myforce3.addPerBondParameter('k')
myforce3.addPerBondParameter('radius')
myforce3.addBond([mygroup1c, mygroup2c], [1.0e-9*kilojoules_per_mole, 11.0*angstroms])
forcenum3 = system.addForce(myforce3)

myforce4 = CustomCentroidBondForce(2, "k*(distance(g1,g2)^2-radius^2)")
mygroup1d = myforce4.addGroup(rec2_indices)
mygroup2d = myforce4.addGroup(lig_indices)
myforce4.setForceGroup(4)
myforce4.addPerBondParameter('k')
myforce4.addPerBondParameter('radius')
myforce4.addBond([mygroup1d, mygroup2d], [1.0e-9*kilojoules_per_mole, 13.0*angstroms])
forcenum4 = system.addForce(myforce4)

myforce5 = CustomCentroidBondForce(2, "-k*(distance(g1,g2)^2-radius^2)")
mygroup1e = myforce5.addGroup(rec3_indices)
mygroup2e = myforce5.addGroup(lig_indices)
myforce5.setForceGroup(5)
myforce5.addPerBondParameter('k')
myforce5.addPerBondParameter('radius')
myforce5.addBond([mygroup1e, mygroup2e], [1.0e-9*kilojoules_per_mole, 17.0*angstroms])
forcenum5 = system.addForce(myforce5)

myforce6 = CustomCentroidBondForce(2, "k*(distance(g1,g2)^2-radius^2)")
mygroup1f = myforce6.addGroup(rec3_indices)
mygroup2f = myforce6.addGroup(lig_indices)
myforce6.setForceGroup(6)
myforce6.addPerBondParameter('k')
myforce6.addPerBondParameter('radius')
myforce6.addBond([mygroup1f, mygroup2f], [1.0e-9*kilojoules_per_mole, 20.0*angstroms])
forcenum6 = system.addForce(myforce6)

integrator = MmvtLangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, "tryp_test_filename.txt")
integrator.addMilestoneGroup(1)
integrator.addMilestoneGroup(2)
integrator.addMilestoneGroup(3)
integrator.addMilestoneGroup(4)
integrator.addMilestoneGroup(5)
integrator.addMilestoneGroup(6)
integrator.setSaveStateFileName("trp_benz")


platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(mypdb.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)

simulation.context.setPeriodicBoxVectors(*box_vector)
    
#simulation.minimizeEnergy()
#simulation.reporters.append(PDBReporter('/home/lvotapka/doc/mmvt_test/tryp_ben_output.pdb', 20))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.step(4000)

#state = simulation.context.getState()
#print(state.getPeriodicBoxVectors())
