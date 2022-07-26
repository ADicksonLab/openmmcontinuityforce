/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2021 Stanford University and the Authors.      *
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

#include "CommonContForceKernels.h"
#include "CommonContForceKernelSources.h"
#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/internal/ContextImpl.h"

using namespace ContForcePlugin;
using namespace OpenMM;
using namespace std;

// macro for checking the result of synchronization operation on CUDA
// copied from `openmm/platforms/cuda/src/CudaParallelKernels.cpp`
#define CHECK_RESULT(result, prefix) \
if (result != CUDA_SUCCESS) { \
	std::stringstream m; \
	m<<prefix<<": "<<cu.getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
	throw OpenMMException(m.str());\
}

void CommonCalcContForceKernel::initialize(const System& system, const ContForce& force) {

  usePeriodic = force.usesPeriodicBoundaryConditions();
  int numContexts = cc.getNumContexts();  
  int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
  int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
  numBonds = endIndex-startIndex;
  if (numBonds == 0)
	return;
  vector<vector<int> > idxs(numBonds, vector<int>);
  vector<int> npart(numBonds);
  vector<mm_float2> paramVector(numBonds);
  
  for (int i = 0; i < numBonds; i++) {
	double length, k;
	force.getBondParameters(startIndex+i, idxs[i], npart[i], length, k);
	paramVector[i] = mm_float2((float) length, (float) k);
  }
  params.upload(paramVector);

  // Inititalize CUDA objects.

  cc.setAsCurrent();
  map<string, string> defines;
  CCmodule program = cc.createModule(CommonContForceKernelSources::ContForce, defines);
  copyInputsKernel = cc.getKernel(program, "copyInputs");
  addForcesKernel = cc.getKernel(program, "addForces");

  /*
    map<string, string> replacements;
    replacements["PARAMS"] = cc.getBondedUtilities().addArgument(params, "float2");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonExampleKernelSources::exampleForce, replacements), force.getForceGroup());
    cc.addForce(new CommonExampleForceInfo(force));
  */
}

double CommonCalcContForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

  int numParticles = cc.getNumAtoms();
  vector<real4> pos;
  context.getPositions(pos);
  
  int numBonds = npart.size();
  double energy = 0;

  for (int i = 0; i < numBonds; i++) {
	// initialize distance nxn matrix with zeros
	vector<vector<double>> dmat( npart[i], vector<double> ( npart[i], 0));

	// calculate distances
	for (int at1 = 0; at1 < npart[i]-1; at1++) {
	  for (int at2 = at1+1; at2 < npart[i]; at2++) {
		// ignore periodic boundaries for now
		real4 delta = pos[idxs[i][at1]]-pos[idxs[i][at2]];
		RealOpenMM r2 = delta.dot(delta);
		dmat[at1][at2] = dmat[at2][at1] = sqrt(r2);
	  }
	}

	// define list of component indexes
	vector<int> comp_idxs( npart[i], -1);

	int curr_comp = 0;
	while (*std::min_element(comp_idxs.begin(), comp_idxs.end()) == -1) {
	  // assign the first -1 index element to curr_comp
	  for (int at = 0; at < npart[i]; at++) {
		if (comp_idxs[at] == -1) {
		  comp_idxs[at] = curr_comp;
		  break;
		}
	  }

	  bool updated = true;
	  while (updated) {
		updated = false;
		  
		// see if any of the -1 elements can be assigned to curr_comp
		for (int at = 0; at < npart[i]; at++) {
		  if (comp_idxs[at] == -1) {
			for (int at2 = 0; at2 < npart[i]; at2++) {
			  if (comp_idxs[at2] == curr_comp) {
				if (dmat[at][at2] < length[i]) {
				  comp_idxs[at] = curr_comp;
				  updated = true;
				}
			  }
			}
		  }
		}
	  }		
	  curr_comp++;
	} // all comp_idxs assigned! (curr_comp holds the number of components)
		  
	  // if more than one component exists:  for each component, find the closest
	  //   outside node and introduce an attractive force between them
	if (curr_comp > 1) {
	  // force_mat will mark which atom pairs will be restrained, to
	  // avoid double counting
	  vector<vector<int>> force_mat( npart[i], vector<int> ( npart[i], 0));

	  for (int comp_idx = 0; comp_idx < curr_comp; comp_idx++) {
		int c_atom_in = -1;
		int c_atom_out = -1;
		double c_dist = -1.0;
		  
		// find closest atom pair
		for (int in_idx = 0; in_idx < npart[i]; in_idx++) {
		  if (comp_idxs[in_idx] == comp_idx) {
			for (int out_idx = 0; out_idx < npart[i]; out_idx++) {
			  if (comp_idxs[out_idx] != comp_idx) {
				if (dmat[in_idx][out_idx] < c_dist | c_dist == -1.0) {
				  c_dist = dmat[in_idx][out_idx];
				  c_atom_in = in_idx;
				  c_atom_out = out_idx;
				}
			  }
			}
		  }
		}
		// mark this pair to be restrained in force_mat
		force_mat[c_atom_in][c_atom_out] = force_mat[c_atom_out][c_atom_in] = 1;
	  }

	  int paddedNumAtoms = cc.getPaddedNumAtoms();
	  if (includeForces) {
		vector<real4> force( paddedNumAtoms, real4 (0));
	  }
	  // add restraint force to designated atom pairs
	  for (int at1 = 0; at1 < npart[i]-1; at1++) {
		for (int at2 = at1+1; at2 < npart[i]; at2++) {
		  if (force_mat[at1][at2] == 1) {
			RealOpenMM dr = (dmat[at1][at2]-length[i]);
			RealOpenMM dr2 = dr*dr;
			if (includeEnergy) {
			  energy += k[i]*dr2;
			}
			if (includeForces) {
			  RealOpenMM dEdR = 2*k[i]*dr;
			  dEdR = (dmat[at1][at2] > 0) ? (dEdR/dmat[at1][at2]) : 0;
			  RealVec delta = pos[idxs[i][at1]]-pos[idxs[i][at2]];
			  
			  force[idxs[i][at1]] -= delta*dEdR;
			  force[idxs[i][at2]] += delta*dEdR;
			}
		  }
		}
	  }
	  if (includeForces) {
		void* forceArgs[] = {&force, &cc.getForce().getDevicePointer(),
							 &cc.getAtomIndexArray().getDevicePointer(), &numParticles, &paddedNumAtoms};
		cc.executeKernel(addForcesKernel, forceArgs, numParticles);

	  }
	}
  }
  return energy;
}

void CommonCalcContForceKernel::copyParametersToContext(ContextImpl& context, const ContForce& force) {
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
	  throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
	  return;
    
    // Record the per-bond parameters.
    
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int npart;
		vector<int> idxs;
        double length, k;
        force.getBondParameters(startIndex+i, idxs, npart, length, k);
        paramVector[i] = mm_float2((float) length, (float) k);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules();
}

