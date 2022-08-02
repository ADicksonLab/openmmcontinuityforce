/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2018 Stanford University and the Authors.           *
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

#include "CudaContForceKernels.h"
#include "CudaContForceKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/reference/RealVec.h"
#include <map>
#include <cuda_runtime_api.h>
#include <algorithm>
#include <cstring>
using namespace ContForcePlugin;
using namespace OpenMM;
using namespace std;





class CudaCalcContForceKernel::CopyForcesTask : public ThreadPool::Task {
public:
  CopyForcesTask(CudaContext& cu, vector<Vec3>& forces) : cu(cu), forces(forces) {
  }
  void execute(ThreadPool& threads, int threadIndex) {
	// Copy the forces applied by PLUMED to a buffer for uploading.  This is done in parallel for speed.

	int numParticles = cu.getNumAtoms();
	int numThreads = threads.getNumThreads();
	int start = threadIndex*numParticles/numThreads;
	int end = (threadIndex+1)*numParticles/numThreads;
	if (cu.getUseDoublePrecision()) {
	  double* buffer = (double*) cu.getPinnedBuffer();
	  for (int i = start; i < end; ++i) {
		const Vec3& p = forces[i];
		buffer[3*i] = p[0];
		buffer[3*i+1] = p[1];
		buffer[3*i+2] = p[2];
	  }
	}
	else {
	  float* buffer = (float*) cu.getPinnedBuffer();
	  for (int i = start; i < end; ++i) {
		const Vec3& p = forces[i];
		buffer[3*i] = (float) p[0];
		buffer[3*i+1] = (float) p[1];
		buffer[3*i+2] = (float) p[2];
	  }
	}
  }
  CudaContext& cu;
  vector<Vec3>& forces;
};




CudaCalcContForceKernel::~CudaCalcContForceKernel() {
}

/**
 * @brief
 *
 * @param system
 * @param force
 * @param nnModule
 */
void CudaCalcContForceKernel::initialize(const System& system, const ContForce& force) {


	int numBonds = force.getNumBonds();
	idxs.resize(numBonds);
	npart.resize(numBonds);
	length.resize(numBonds);	
	k.resize(numBonds);
	for (int i = 0; i < numBonds; i++) {
	  idxs[i].resize(1);
	  force.getBondParameters(i, idxs[i], npart[i], length[i], k[i]);
	}
	// Inititalize CUDA objects.
	cu.setAsCurrent();
	cuStreamCreate(&stream, CU_STREAM_NON_BLOCKING);
	cuEventCreate(&syncEvent, CU_EVENT_DISABLE_TIMING);
	int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
	contForces = new CudaArray(cu, 3*system.getNumParticles(), elementSize, "contForces");
	map<string, string> defines;
	defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
	defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
	CUmodule module = cu.createModule(CudaContForceKernelSources::ContForce, defines);
	addForcesKernel = cu.getKernel(module, "addForces");

}

double CudaCalcContForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

	vector<Vec3> pos, forces;
	context.getPositions(pos);
	int numParticles = cu.getNumAtoms();
	//context.getForces(forces);
	forces.resize(numParticles);
	memset(&forces[0], 0, numParticles*sizeof(Vec3));

	int numBonds = npart.size();
	double energy = 0;

	for (int i = 0; i < numBonds; i++) {
	  // initialize distance nxn matrix with zeros
	  vector<vector<double>> dmat( npart[i], vector<double> ( npart[i], 0));

	  // calculate distances
	  for (int at1 = 0; at1 < npart[i]-1; at1++) {
		for (int at2 = at1+1; at2 < npart[i]; at2++) {
		  // ignore periodic boundaries for now
		  RealVec delta = pos[idxs[i][at1]]-pos[idxs[i][at2]];
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

		// add restraint force to designated atom pairs
		for (int at1 = 0; at1 < npart[i]-1; at1++) {
		  for (int at2 = at1+1; at2 < npart[i]; at2++) {
			if (force_mat[at1][at2] == 1) {
			  RealOpenMM dr = (dmat[at1][at2]-length[i]);
			  RealOpenMM dr2 = dr*dr;
			  if (includeEnergy) {
				energy += k[i]*dr2;
			  }
				RealOpenMM dEdR = 2*k[i]*dr;
				dEdR = (dmat[at1][at2] > 0) ? (dEdR/dmat[at1][at2]) : 0;
				RealVec delta = pos[idxs[i][at1]]-pos[idxs[i][at2]];

				forces[idxs[i][at1]] -= delta*dEdR;
				forces[idxs[i][at2]] += delta*dEdR;

			}
		  }
		}
	  }
	}

	if (includeForces) {
	  CopyForcesTask task(cu, forces);
	  cu.getPlatformData().threads.execute(task);
	  cu.getPlatformData().threads.waitForThreads();
	  cu.setAsCurrent();
	  cuMemcpyHtoDAsync(contForces->getDevicePointer(), cu.getPinnedBuffer(), contForces->getSize()*contForces->getElementSize(), stream);
	  cuEventRecord(syncEvent, stream);
	  // Wait until executeOnWorkerThread() is finished.

	  cu.getWorkThread().flush();
	  cuStreamWaitEvent(cu.getCurrentStream(), syncEvent, 0);
	  void* args[] = {&contForces->getDevicePointer(), &cu.getForce().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer()};
	  cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
	}
	return energy;
}

void CudaCalcContForceKernel::copyParametersToContext(ContextImpl& context, const ContForce& force) {
	if (force.getNumBonds() != npart.size())
		throw OpenMMException("updateParametersInContext: The number of bonds has changed");
	for (int i = 0; i < force.getNumBonds(); i++) {
	  vector<int> test_idxs(1,1);
	  int test_npart;
	  force.getBondParameters(i, test_idxs, test_npart, length[i], k[i]);
	  if (test_npart != npart[i] || test_idxs != idxs[i])
		throw OpenMMException("updateParametersInContext: A particle index has changed");
	}
}
