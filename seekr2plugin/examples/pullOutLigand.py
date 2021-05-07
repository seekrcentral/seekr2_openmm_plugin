from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from mmvtplugin import MmvtLangevinIntegrator, vectori, vectord
import mmvtplugin
from time import time
import numpy as np
import os
import parmed

lig_indices = [3221, 3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229]
#3221 3222 3223 3224 3225 3226 3227 3228 3229
rec_indices = [2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926]
# 2478 2489 2499 2535 2718 2745 2769 2787 2794 2867 2926

prmtop_filename = 'tryp_ben.prmtop'
inpcrd_filename = 'tryp_ben.inpcrd'

prmtop = AmberPrmtopFile(prmtop_filename)
inpcrd = AmberInpcrdFile(inpcrd_filename)



mypdb = PDBFile('tryp_ben_new.pdb')

max_cycles = 4000

# radii are the locations of the milestones
radii = [11.0, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0]
num_steps=100 # steps per cycle

def test_for_upward_crossings(crossings_file_name, save_state_filename):
    if not os.path.exists(crossings_file_name):
        return None
    my_file = open(crossings_file_name,'r')
    for line in my_file:
        line = line.strip().split()
        milestone = line[0]
        if milestone == '2': # upward transition!
            line_num = line[1]
            state_filename = save_state_filename+'_%s_%s' % (milestone, line_num)
            return state_filename
    return None
            

def run_steps(inner_radius, outer_radius, crossings_file_name, save_state_filename, start_state):
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    myforce1 = CustomCentroidBondForce(2, "step(-k*(distance(g1, g2)^2 - radius^2))")
    mygroup1a = myforce1.addGroup(rec_indices)
    mygroup2a = myforce1.addGroup(lig_indices)
    myforce1.setForceGroup(1)
    myforce1.addPerBondParameter('k')
    myforce1.addPerBondParameter('radius')
    myforce1.addBond([mygroup1a, mygroup2a], [1.0e-9*kilojoules_per_mole, inner_radius*angstroms])
    forcenum1 = system.addForce(myforce1)

    myforce2 = CustomCentroidBondForce(2, "step(k*(distance(g1,g2)^2-radius^2))")
    mygroup1b = myforce2.addGroup(rec_indices)
    mygroup2b = myforce2.addGroup(lig_indices)
    myforce2.setForceGroup(2)
    myforce2.addPerBondParameter('k')
    myforce2.addPerBondParameter('radius')
    myforce2.addBond([mygroup1b, mygroup2b], [1.0e-9*kilojoules_per_mole, outer_radius*angstroms])
    forcenum2 = system.addForce(myforce2)

    integrator = MmvtLangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, crossings_file_name)
    integrator.addMilestoneGroup(1)
    integrator.addMilestoneGroup(2)
    integrator.setSaveStateFileName(save_state_filename)

    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    
    if start_state == None:
        simulation.context.setPositions(mypdb.positions)
        simulation.context.setVelocitiesToTemperature(300*kelvin)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.minimizeEnergy()
    else:
        simulation.loadState(start_state)
        state = simulation.context.getState(getPositions = True, enforcePeriodicBox = False)
        positions = state.getPositions()
        amber_parm = parmed.amber.AmberParm(prmtop_filename, inpcrd_filename)
        amber_parm.positions = positions
        pdb_save_filename = save_state_filename+'_%f.pdb' % inner_radius
        amber_parm.save(pdb_save_filename, overwrite=True)
    
        
    
    simulation.reporters.append(PDBReporter('tryp_ben_output.pdb', 20))
    #simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
    #        potentialEnergy=True, temperature=True, volume=True))
    
    cycle_counter = 0
    while cycle_counter < max_cycles:
        simulation.step(num_steps)
        crossed_upward = test_for_upward_crossings(crossings_file_name, save_state_filename)
        if crossed_upward:
            
            return crossed_upward
        cycle_counter += 1
        
    return None

start_state = None
for i in range(len(radii)-1):
print('loading start_state:', start_state)
    inner_radius = radii[i]
    outer_radius = radii[i+1]
    crossings_file_name = 'crossing_data_%d.dat' % i
    if os.path.exists(crossings_file_name):
        os.system('rm %s' % crossings_file_name)
    save_state_filename = 'state%d' % i
    print("running simulation between", inner_radius, 'and', outer_radius)
    start_state = run_steps(inner_radius, outer_radius, crossings_file_name, save_state_filename, start_state)
    
    
