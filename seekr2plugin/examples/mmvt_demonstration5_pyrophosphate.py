from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from seekr2plugin import MmvtLangevinIntegrator, vectori, vectord
import seekr2plugin
from time import time
import numpy as np
import os


index_list = [4539, 4534, 2337]

prmtop = AmberPrmtopFile('pyro.parm7')
inpcrd = AmberInpcrdFile('pyro.rst7')

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
#system = prmtop.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
#system = prmtop.createSystem(implicitSolvent=GBn2)

mypdb = PDBFile('pyro_single_frame.pdb')


myforce1 = CustomCentroidBondForce(3, "-k*(angle(g1, g2, g3) - ref_angle)")
mygroup1a = myforce1.addGroup([index_list[0]])
mygroup2a = myforce1.addGroup([index_list[1]])
mygroup3a = myforce1.addGroup([index_list[2]])
myforce1.setForceGroup(1)
myforce1.addPerBondParameter('k')
myforce1.addPerBondParameter('ref_angle')
myforce1.addBond([mygroup1a, mygroup2a, mygroup3a], [1.0e-9*kilojoules_per_mole, radians*np.pi*8.0/9.0])
forcenum1 = system.addForce(myforce1)

'''
myforce2 = CustomCentroidBondForce(3, "k*(angle(g1, g2, g3) - ref_angle)")
mygroup1b = myforce2.addGroup([index_list[0]])
mygroup2b = myforce2.addGroup([index_list[1]])
mygroup3b = myforce2.addGroup([index_list[2]])
myforce2.setForceGroup(2)
myforce2.addPerBondParameter('k')
myforce2.addPerBondParameter('ref_angle')
myforce2.addBond([mygroup1b, mygroup2b, mygroup3b], [1.0e-9*kilojoules_per_mole, radians*np.pi*17.0/18.0])
forcenum2 = system.addForce(myforce2)
'''
integrator = MmvtLangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, "pyro_output.txt")
integrator.addMilestoneGroup(1)
#integrator.addMilestoneGroup(2)


platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(mypdb.positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
#simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('pyro_output.pdb', 10))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.step(5000)

