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

#include "ReferenceExampleKernels.h"
#include "ExampleForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

void ReferenceCalcExampleForceKernel::initialize(const System& system, const ExampleForce& force) {
    // Initialize bond parameters.
    
    int numBonds = force.getNumBonds();
    idxs.resize(numBonds);
    npart.resize(numBonds);
    length.resize(numBonds);
    k.resize(numBonds);
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(i, idxs[i], npart[i], length[i], k[i]);
}

double ReferenceCalcExampleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    int numBonds = npart.size();
    double energy = 0;

    // Compute the interactions.
    
    for (int i = 0; i < numBonds; i++) {
	  for (int at_idx = 0; at_idx < npart[i]; at_idx++) {
		// find the closest distance from at_idx to the others
		c_dist = -1;
		c_idx = -1;
		for (int at2 = 0; at2 < npart[i]; at2++) {
		  if (at_idx != at2) {
			// ignore periodic boundaries for now
			RealVec delta = pos[idxs[i][at_idx]]-pos[idxs[i][at2]];
			RealOpenMM r2 = delta.dot(delta);
			RealOpenMM r = sqrt(r2);
			if (r < c_dist || c_dist == -1) {
			  c_dist = r;
			  c_idx = at2;
			}
		  }
		}
		// if closest distance > length, apply an attractive force between at_idx and at2
		if (c_dist > length[i]) {
		  RealOpenMM dr = (c_dist-length[i]);
		  RealOpenMM dr2 = dr*dr;
		  if (includeEnergy) {
			energy += k[i]*dr2;
		  }
		  if (includeForces) {
			RealOpenMM dEdR = 2*k[i]*dr;
			dEdR = (c_dist > 0) ? (dEdR/c_dist) : 0;
			RealVec delta = pos[idxs[i][at_idx]]-pos[idxs[i][c_idx]];
					
			force[idxs[i][at_idx]] -= delta*dEdR;
			force[idxs[i][c_idx]] += delta*dEdR;
		  }
		}
	  }
    }
    return energy;
}

void ReferenceCalcExampleForceKernel::copyParametersToContext(ContextImpl& context, const ExampleForce& force) {
    if (force.getNumBonds() != npart.size())
        throw OpenMMException("updateParametersInContext: The number of Example bonds has changed");
    for (int i = 0; i < force.getNumBonds(); i++) {
        vector<int> test_idxs;
        int test_npart;
        force.getBondParameters(i, test_idxs, test_npart, length[i], k[i]);
        if (test_npart != npart[i] || test_idxs != idxs[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
}
