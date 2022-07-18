/* -------------------------------------------------------------------------- *
 *                                OpenMMExample                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ExampleForce.h"
#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

extern "C" void registerExampleSerializationProxies();

void testSerialization() {
    // Create a Force.

    ExampleForce force;
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
    XmlSerializer::serialize<ExampleForce>(&force, "Force", buffer);
    ExampleForce* copy = XmlSerializer::deserialize<ExampleForce>(buffer);

    // Compare the two forces to see if they are identical.

    ExampleForce& force2 = *copy;
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
