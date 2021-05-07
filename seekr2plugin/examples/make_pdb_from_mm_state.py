from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import parmed
import sys

state_files = sys.argv[1:]

print("number of files:", len(state_files))

prmtop_filename = 'tryp_ben.prmtop'
inpcrd_filename = 'tryp_ben.inpcrd'

prmtop = AmberPrmtopFile(prmtop_filename)
inpcrd = AmberInpcrdFile(inpcrd_filename)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

counter = 0
for state_file in state_files:
    simulation.loadState(state_file)
    state = simulation.context.getState(getPositions = True, enforcePeriodicBox = False)
    positions = state.getPositions()
    amber_parm = parmed.amber.AmberParm(prmtop_filename, inpcrd_filename)
    amber_parm.positions = positions
    pdb_save_filename = state_file+'.pdb'
    amber_parm.save(pdb_save_filename, overwrite=True)
    counter += 1
    if counter%10 == 0:
        print("frame:", counter)
    
'''
# VMD commands:

ls trp_benz_0_*.pdb
mol load pdb trp_benz_0_1.pdb
foreach filename [lrange [ls trp_benz_0_*.pdb] 2 end] {puts $filename; animate read pdb $filename}

# align protein
set lig [atomselect top "index 3221 3222 3223 3224 3225 3226 3227 3228 3229"]
draw color orange
for {set i 0} {$i < [molinfo top get numframes]} {incr i} {puts $i; $lig frame $i; draw point [measure center $lig weight [$lig get mass]]}


mol load pdb trp_benz_1_0.pdb
'''
