from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from seekr2plugin import MmvtLangevinIntegrator, vectori, vectord
import seekr2plugin
from time import time
import numpy as np
import os

def writePdbFrame(filename, frameNum, state):
    posInNm = state.getPositions()
    myfile = open(filename, 'a')
    # Use PDB MODEL cards to number trajectory frames
    myfile.write("MODEL     %d\n" % frameNum) # start of frame
    #for (int a = 0; a < (int)posInNm.size(); ++a)
    for a in range(0, len(posInNm)):
        myfile.write("ATOM  %5d  AR   AR     1    " % (a+1)) # atom number
        myfile.write("%8.3f%8.3f%8.3f  1.00  0.00\n" % (posInNm[a][0].value_in_unit(angstroms), posInNm[a][1].value_in_unit(angstroms), posInNm[a][2].value_in_unit(angstroms)))
        
    myfile.write("ENDMDL\n") # end of frame
    myfile.close()

filename = "argon_sphere_out.pdb"
box_edge = 25 * angstroms
nparticles = 10
mass = 39.9 * amu # mass
sigma = 3.4 * angstrom # Lennard-Jones sigma
epsilon = 0.238 * kilocalories_per_mole # Lennard-Jones well-depth
charge = 0.0 * elementary_charge # argon model has no charge

cutoff = 2.5*sigma # Compute cutoff
#pdb = PDBFile('argons.pdb')
#forcefield = ForceField('amber14-all.xml', 'argon.xml')
#system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
system = System()
system.setDefaultPeriodicBoxVectors(np.array([box_edge, 0, 0]), np.array([0, box_edge, 0]), np.array([0, 0, box_edge]))
force = CustomNonbondedForce("""4*epsilon*(1/(((r/sigma)^6)^2) - 1/((r/sigma)^6));
                                   sigma=0.5*(sigma1+sigma2);
                                   epsilon=sqrt(epsilon1*epsilon2);
                                   """)
force.addPerParticleParameter("sigma")
force.addPerParticleParameter("epsilon")
force.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
force.setCutoffDistance(cutoff) 
for particle_index in range(nparticles):
    system.addParticle(mass)
    force.addParticle([sigma, epsilon])

system.addForce(force)


positions = Quantity(numpy.random.uniform(high=box_edge/angstroms, size=[nparticles,3]), angstrom)

positions[0][0] = 12 * angstrom
positions[0][1] = 12 * angstrom
positions[0][2] = 12 * angstrom

myforce = CustomCentroidBondForce(1, "k*((x1-centerx)^2 + (y1-centery)^2 + (z1-centerz)^2 - radius^2)")
mygroup1 = myforce.addGroup([0]) # atom index 0
myforce.setForceGroup(1)
myforce.addPerBondParameter('k')
myforce.addPerBondParameter('centerx')
myforce.addPerBondParameter('centery')
myforce.addPerBondParameter('centerz')
myforce.addPerBondParameter('radius')
myforce.addBond([mygroup1], [1.0e-9*kilojoules_per_mole, 12.0*angstroms, 12.0*angstroms, 12.0*angstroms, 10.0*angstroms])
forcenum = system.addForce(myforce)

integrator = MmvtLangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, "test_filename.txt")
integrator.addMilestoneGroup(1)

platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'} # CUDA platform
context = Context(system, integrator)
context.setPositions(positions)
context.setVelocitiesToTemperature(300*kelvin)
LocalEnergyMinimizer.minimize(context)

if os.path.exists(filename):
    os.system('rm %s' % filename) # delete the existing transition data file if it exists

for iteration in range(2000):
    integrator.step(10)
    state = context.getState(getPositions=True)
    positions = state.getPositions(asNumpy=True)
    writePdbFrame(filename, iteration, state)

'''
simulation = Simulation(pdb.topology, system, integrator, platform, properties) 

#simulation = Simulation(pdb.topology, system, integrator) # Reference platform
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('demonstration1_argon_sphere.pdb', 10))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
starttime = time()
simulation.step(1000)
print "time:", time() - starttime

'''
