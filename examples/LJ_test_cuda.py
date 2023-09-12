import simtk.openmm as omm
import simtk.openmm.app as omma
import openmm.unit as unit
import numpy as np

from contforceplugin import ContForce

openmm_path = sys.argv[1]
plugins_path = osp.join(openmm_path,'lib','plugins')

if not os.path.exists(openmm_path):
    print(f"Error! openmm_path: {openmm_path} does not exist")
elif not os.path.exists(plugins_path):
    print(f"Error! openmm_path: {openmm_path} does not contain lib/plugins directory")
    print("e.g.: /path/to/your/conda/pkgs/openmm-7.7.0-py39h9717219_0/lib/plugins")

omm.Platform.loadPluginsFromDirectory(plugins_path)

npart = 20
mass = 10
sigma = 0.1
d = 0.3

system = omm.System()
for i in range(npart):
    system.addParticle(mass)

force = ContForce()

force.addBond(list(range(npart)),npart,d,2000) # ind1, ind2, length, k
system.addForce(force)

nbforce = omm.NonbondedForce()
for i in range(npart):
    nbforce.addParticle(0.0, sigma, 1.0) # charge, sigma (nm), epsilon
system.addForce(nbforce)

integrator = omm.LangevinIntegrator(300,100.0,0.002)  # using a high friction coefficient
platform = omm.Platform.getPlatformByName('CUDA')
top = omma.Topology()
ch = top.addChain()
for i in range(npart):
    R = top.addResidue(f'R{i}',ch)
    top.addAtom(f'Ar{i}',omma.Element.getBySymbol('Ar'),R)

simulation = omma.Simulation(top,system,integrator,platform=platform)
pos = 2*np.random.random(size=(npart,3))
simulation.context.setPositions(pos)

simulation.minimizeEnergy()

simulation.reporters.append(omma.PDBReporter('output.pdb', 100))
simulation.step(20000)
state = simulation.context.getState(getPositions=True,getVelocities=True)

pos = np.array(state.getPositions().value_in_unit(unit.nanometer))

dmat = np.array([[np.sqrt(np.sum(np.square(pos[i]-pos[j]))) for i in range(npart)] for j in range(npart)])
connections = np.array(dmat < d+0.2, dtype=int)

print("Final connection matrix:")
print(connections)
