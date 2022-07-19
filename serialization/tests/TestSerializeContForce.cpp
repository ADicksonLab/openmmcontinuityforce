#include "ContForce.h"
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace ContForcePlugin;
using namespace OpenMM;
using namespace std;

extern "C" void registerExampleSerializationProxies();

void testSerialization() {
    // Create a Force.

    ContForce force;
	vector<int> idxs1 = {0,1,2,3};
    force.addBond(idxs1, 4, 1.0, 2.0);
	vector<int> idxs2 = {10,11,12};
    force.addBond(idxs2, 3, 2.0, 2.1);
	vector<int> idxs3 = {19,4,16,12,7};
    force.addBond(idxs3, 5, 3.0, 2.2);
	vector<int> idxs4 = {29,300,301,5,3,0};
    force.addBond(idxs4, 6, 4.0, 2.3);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<ContForce>(&force, "Force", buffer);
    ContForce* copy = XmlSerializer::deserialize<ContForce>(buffer);

    // Compare the two forces to see if they are identical.

    ContForce& force2 = *copy;
    ASSERT_EQUAL(force.getNumBonds(), force2.getNumBonds());
    for (int i = 0; i < force.getNumBonds(); i++) {
	  vector<int> a1, b1;
	  int a2, b2;
	  double da, db, ka, kb;
	  force.getBondParameters(i, a1, a2, da, ka);
	  force2.getBondParameters(i, b1, b2, db, kb);

	  ASSERT_EQUAL(a2, b2);
	  for (int idx = 0; idx < a2; idx++) {
		ASSERT_EQUAL(a1[idx], b1[idx]);
	  }
	  ASSERT_EQUAL(da, db);
	  ASSERT_EQUAL(ka, kb);
    }
}

int main() {
    try {
        registerExampleSerializationProxies();
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
