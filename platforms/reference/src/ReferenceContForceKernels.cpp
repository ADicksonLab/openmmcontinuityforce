#include "ReferenceContForceKernels.h"
#include "ContForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include<algorithm>

using namespace ContForcePlugin;
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

void ReferenceCalcContForceKernel::initialize(const System& system, const ContForce& force) {
    // Initialize bond parameters.
    
    int numBonds = force.getNumBonds();
    idxs.resize(numBonds);
    npart.resize(numBonds);
    length.resize(numBonds);
    k.resize(numBonds);
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(i, idxs[i], npart[i], length[i], k[i]);
}

double ReferenceCalcContForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
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
	  }
    }
    return energy;
}

void ReferenceCalcContForceKernel::copyParametersToContext(ContextImpl& context, const ContForce& force) {
    if (force.getNumBonds() != npart.size())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    for (int i = 0; i < force.getNumBonds(); i++) {
        vector<int> test_idxs;
        int test_npart;
        force.getBondParameters(i, test_idxs, test_npart, length[i], k[i]);
        if (test_npart != npart[i] || test_idxs != idxs[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
}
