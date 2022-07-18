/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
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

/**
 * This tests the Reference implementation of ExampleForce.
 */

#include "ExampleForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include <cmath>
#include <iostream>
#include <vector>

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

extern "C" OPENMM_EXPORT void registerExampleReferenceKernelFactories();

void testForce() {
    // Create a system of 10 atoms connected with a continuity force
    
    const int numParticles = 10;
    System system;
    vector<Vec3> positions(numParticles);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
		if (i != 9) {
		  positions[i] = Vec3(i, 0.7, 0.5);
		} else {
		  positions[i] = Vec3(i, 0.7, 2.1);
		}
    }
    ExampleForce* force = new ExampleForce();
    system.addForce(force);
	vector<int> idxs = {0,1,2,3,4,5,6,7,8,9};
	force->addBond(idxs, 10, 1.0, 17);
    
    // Compute the forces and energy.

    VerletIntegrator integ(1.0);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    
    // See if the energy is correct.
    
	double length = 1.0;
	double k = 17;
	Vec3 delta = positions[9]-positions[8];
	double dr = sqrt(delta.dot(delta))-length;
    double expectedEnergy = k*dr*dr;

    ASSERT_EQUAL_TOL(expectedEnergy, state.getPotentialEnergy(), 1e-5);

    // Validate the forces by moving each particle along each axis, and see if the energy changes by the correct amount.
    
    double offset = 1e-3;
    for (int i = 0; i < numParticles; i++)
        for (int j = 0; j < 3; j++) {
            vector<Vec3> offsetPos = positions;
            offsetPos[i][j] = positions[i][j]-offset;
            context.setPositions(offsetPos);
            double e1 = context.getState(State::Energy).getPotentialEnergy();
            offsetPos[i][j] = positions[i][j]+offset;
            context.setPositions(offsetPos);
            double e2 = context.getState(State::Energy).getPotentialEnergy();
            ASSERT_EQUAL_TOL(state.getForces()[i][j], (e1-e2)/(2*offset), 1e-2);
        }
}

void testChangingParameters() {
    const double k = 1.5;
    const double length = 0.5;
    Platform& platform = Platform::getPlatformByName("Reference");
    
    // Create a system with one bond.
    
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    ExampleForce* force = new ExampleForce();
	vector<int> idxs = {0,1};
    force->addBond(idxs, 2, length, k);
    system.addForce(force);
    vector<Vec3> positions(2);
    positions[0] = Vec3(1, 0, 0);
    positions[1] = Vec3(2, 0, 0);
    
    // Check the energy.
    
    VerletIntegrator integ(1.0);
    Context context(system, integ, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(k*pow(1.0-length, 2), state.getPotentialEnergy(), 1e-5);
    
    // Modify the parameters.
    
    const double k2 = 2.2;
    const double length2 = 0.9;
    force->setBondParameters(0, idxs, 2, length2, k2);
    force->updateParametersInContext(context);
    state = context.getState(State::Energy);
    ASSERT_EQUAL_TOL(k2*pow(1.0-length2, 2), state.getPotentialEnergy(), 1e-5);
}

int main() {
    try {
        registerExampleReferenceKernelFactories();
        testForce();
        testChangingParameters();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
