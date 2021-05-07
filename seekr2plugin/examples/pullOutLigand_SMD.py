from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import time
import numpy as np
import os
import parmed

# Feel free to modify this section below

rec_indices = [213, 215, 218, 232, 320, 322, 580, 582, 1107, 1109, 1444, 1446, 
    1449, 1483, 1486, 1489, 1492, 1519, 2266, 2268, 2271, 2274]
#213 215 218 232 320 322 580 582 1107 1109 1444 1446 1449 1483 1486 1489 1492 1519 2266 2268 2271 2274
#lig_indices = [4753, 4754, 4755, 4756, 4757, 4758, 4759, 4760, 4761, 4762, 
#    4763, 4764, 4765, 4766, 4767, 4768, 4769, 4784, 4785, 4786, 4787,]
#4753 4754 4755 4756 4757 4758 4759 4760 4761 4762 4763 4764 4765 4766 4767 4768 4769 4784 4785 4786 4787
lig_indices = list(range(4753, 4804))
prmtop_filename = "system_TP4EW.parm7"
inpcrd_filename = "system_TP4EW.inpcrd"
pdb_filename = "system.pdb"
temperature = 300*kelvin
spring_constant = 9000.0*kilojoules_per_mole

#optionally print a trajectory
trajectory_filename = "jak_stat_out.pdb"
trajectory_frequency = 100000

# total number of timesteps to take in the CMD simulation
total_num_steps = 10000000

# How many "windows" from the bound state to the unbound state
num_windows = 100
show_state_output = False

# Radii of anchor locations. Warning: don't include zero, nor the lowest anchor
#  because the structure defined by pdb_filename should already be the lowest 
#  bound state.
target_radii = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]

#### Don't modify below this point

state_filename = "state.xml"
basename = "radius"
steps_per_window = total_num_steps // num_windows

prmtop = AmberPrmtopFile(prmtop_filename)
inpcrd = AmberInpcrdFile(inpcrd_filename)
mypdb = PDBFile(pdb_filename)

parmed_struct = parmed.load_file(prmtop_filename, xyz=pdb_filename)

def get_lig_rec_distance(parmed_struct, positions, lig_atom_list, rec_atom_list):
    parmed_struct.coordinates = positions
    center_of_mass_1 = np.array([[0., 0., 0.]])
    total_mass = 0.0
    for atom_index in lig_atom_list:
        atom_pos = parmed_struct.coordinates[atom_index,:]
        atom_mass = parmed_struct.atoms[atom_index].mass
        center_of_mass_1 += atom_mass * atom_pos
        total_mass += atom_mass
    center_of_mass_1 = center_of_mass_1 / total_mass
    
    center_of_mass_2 = np.array([[0., 0., 0.]])
    total_mass = 0.0
    for atom_index in rec_atom_list:
        atom_pos = parmed_struct.coordinates[atom_index,:]
        atom_mass = parmed_struct.atoms[atom_index].mass
        center_of_mass_2 += atom_mass * atom_pos
        total_mass += atom_mass
    center_of_mass_2 = center_of_mass_2 / total_mass
    distance = np.linalg.norm(center_of_mass_2 - center_of_mass_1)
    return distance

def run_window(target_radius, save_state_filename, index):
    system = prmtop.createSystem(
        nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    myforce1 = CustomCentroidBondForce(
        2, "0.5*k*(distance(g1, g2) - radius)^2")
    mygroup1a = myforce1.addGroup(rec_indices)
    mygroup2a = myforce1.addGroup(lig_indices)
    myforce1.setForceGroup(1)
    myforce1.addPerBondParameter('k')
    myforce1.addPerBondParameter('radius')
    myforce1.addBond([mygroup1a, mygroup2a], [spring_constant, 
                                              target_radius*angstroms])
    forcenum1 = system.addForce(myforce1)

    integrator = LangevinIntegrator(temperature, 1/picosecond, 
        0.002*picoseconds)
    platform = Platform.getPlatformByName('CUDA')
    properties = {'CudaDeviceIndex': '1', 'CudaPrecision': 'mixed'}
    simulation = Simulation(prmtop.topology, system, integrator, platform, 
                            properties)
                            
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
    if trajectory_filename and trajectory_frequency:
        simulation.reporters.append(pdb_reporter)
    
    if index == 0:
        simulation.context.setPositions(mypdb.positions)
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.minimizeEnergy()
        
    else:
        simulation.loadState(save_state_filename)
        
    if show_state_output:
        simulation.reporters.append(StateDataReporter(stdout, steps_per_window, step=True,
            potentialEnergy=True, temperature=True, volume=True))
    
    
    simulation.step(steps_per_window)
    state = simulation.context.getState(getPositions = True, 
                                        getVelocities = True,
                                        enforcePeriodicBox = True)
    positions = state.getPositions()
    distance = get_lig_rec_distance(parmed_struct, positions, lig_indices, rec_indices)
    simulation.saveState(state_filename)
    return distance, positions
    
if trajectory_filename and trajectory_frequency:
     pdb_reporter = PDBReporter(trajectory_filename, trajectory_frequency)

start_radius = get_lig_rec_distance(parmed_struct, mypdb.positions, lig_indices, rec_indices)
print("start_radius:", start_radius)
start_window = start_radius
last_window = target_radii[-1] + 1.0
increment = (last_window - start_window)/num_windows
print("simulating steered MD in windows from", start_window, "to", last_window, 
    "in increments of", increment)
windows = np.arange(start_window, last_window, increment)
old_distance = start_radius
distance = start_radius
goal_radius_index = 0

for i, window_radius in enumerate(windows):
    print("running window:", window_radius)
    goal_radius = target_radii[goal_radius_index]
    assert goal_radius > start_radius, "Error: your system is starting with "\
        "a ligand-receptor distance of %f, but a target radius is "\
        "listed at a distance of %f. This may never be reached. Please choose "\
        "a list of target radii that are larger than the starting ligand-"\
        "receptor distance."
    old_distance = distance
    distance, positions = run_window(window_radius, state_filename, i)
    print("current ligand-receptor distance:", distance)
    if (distance-goal_radius)*(old_distance-goal_radius) < 0:
        # then a crossing event has taken place
        print("radius crossed:", goal_radius, "saving state")
        amber_parm = parmed.amber.AmberParm(prmtop_filename, inpcrd_filename)
        amber_parm.positions = positions
        pdb_save_filename = basename+'_%.2f.pdb' % distance
        amber_parm.save(pdb_save_filename, overwrite=True)
        goal_radius_index += 1
    
