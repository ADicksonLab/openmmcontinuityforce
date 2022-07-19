#include "ContForceProxy.h"
#include "ContForce.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>

using namespace ContForcePlugin;
using namespace OpenMM;
using namespace std;

ContForceProxy::ContForceProxy() : SerializationProxy("ContForce") {
}

void ContForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const ContForce& force = *reinterpret_cast<const ContForce*>(object);
    SerializationNode& bonds = node.createChildNode("Bonds");
    for (int i = 0; i < force.getNumBonds(); i++) {
	  vector<int> idxs;
	  int npart;
	  double distance, k;
	  force.getBondParameters(i, idxs, npart, distance, k);
	  SerializationNode& bond = bonds.createChildNode("Bond");
	  bond.setIntProperty("npart", npart).setDoubleProperty("d", distance).setDoubleProperty("k", k);
	  for (int idx = 0; idx < npart; idx++) {
		bond.createChildNode("Index").setIntProperty("idx",idxs[idx]);
	  }
    }
}

void* ContForceProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    ContForce* force = new ContForce();
    try {
        const SerializationNode& bonds = node.getChildNode("Bonds");
        for (int i = 0; i < (int) bonds.getChildren().size(); i++) {
            const SerializationNode& bond = bonds.getChildren()[i];
			vector<int> idxs ((int) bond.getChildren().size(), -1);
			for (int idx = 0; idx < (int) bond.getChildren().size(); idx++) {
			  const SerializationNode& at_idx = bond.getChildren()[idx];
			  idxs[idx] = at_idx.getIntProperty("idx");
			}
            force->addBond(idxs, bond.getIntProperty("npart"), bond.getDoubleProperty("d"), bond.getDoubleProperty("k"));
        }
    }
    catch (...) {
        delete force;
        throw;
    }
    return force;
}
