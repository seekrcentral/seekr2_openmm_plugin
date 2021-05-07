from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from seekr2plugin import MmvtLangevinIntegrator, vectori, vectord
import seekr2plugin
from time import time
import numpy as np
import os


pair1 = [82, 228]
pair2 = [144, 221]
pair3 = [36, 89]

prmtop = AmberPrmtopFile('trp_cage_linear.parm7')
inpcrd1 = AmberInpcrdFile('trp_cage_linear.rst7')
#inpcrd2 = AmberInpcrdFile('trp_cage.inpcrd')

#system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
#system = prmtop.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
system = prmtop.createSystem(implicitSolvent=GBn2)

#mypdb = PDBFile('/home/lvotapka/doc/mmvt_test/trp_cage_rmsd.pdb')


myforce1 = CustomCentroidBondForce(6, "-k*(distance(g1,g2)^2 + distance(g3,g4)^2 + distance(g5,g6)^2 - radius^2)")
mygroup1a = myforce1.addGroup([pair1[0]])
mygroup2a = myforce1.addGroup([pair1[1]])
mygroup3a = myforce1.addGroup([pair2[0]])
mygroup4a = myforce1.addGroup([pair2[1]])
mygroup5a = myforce1.addGroup([pair3[0]])
mygroup6a = myforce1.addGroup([pair3[1]])
myforce1.setForceGroup(1)
myforce1.addPerBondParameter('k')
myforce1.addPerBondParameter('radius')
myforce1.addBond([mygroup1a, mygroup2a, mygroup3a, mygroup4a, mygroup5a, mygroup6a], [1.0e-9*kilojoules_per_mole, 6.0*angstroms])
forcenum1 = system.addForce(myforce1)

myforce2 = CustomCentroidBondForce(6, "-k*(distance(g1,g2)^2 + distance(g3,g4)^2 + distance(g5,g6)^2 - radius^2)")
mygroup1b = myforce2.addGroup([pair1[0]])
mygroup2b = myforce2.addGroup([pair1[1]])
mygroup3b = myforce2.addGroup([pair2[0]])
mygroup4b = myforce2.addGroup([pair2[1]])
mygroup5b = myforce2.addGroup([pair3[0]])
mygroup6b = myforce2.addGroup([pair3[1]])
myforce2.setForceGroup(2)
myforce2.addPerBondParameter('k')
myforce2.addPerBondParameter('radius')
myforce2.addBond([mygroup1b, mygroup2b, mygroup3b, mygroup4b, mygroup5b, mygroup6b], [1.0e-9*kilojoules_per_mole, 7.0*angstroms])
forcenum2 = system.addForce(myforce2)

integrator = MmvtLangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, "trp_cage_filename.txt")
integrator.addMilestoneGroup(1)
integrator.addMilestoneGroup(2)


platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
#simulation.context.setPositions(inpcrd2.positions)
simulation.context.setPositions(inpcrd1.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)

#if inpcrd.boxVectors is not None:
#    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('trp_cage_output.pdb', 20))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.step(20000)

