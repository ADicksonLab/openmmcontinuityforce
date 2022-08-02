#ifndef OPENMM_CONTFORCE_H_
#define OPENMM_CONTFORCE_H_

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

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "internal/windowsExportExample.h"

namespace ContForcePlugin {

class OPENMM_EXPORT_EXAMPLE ContForce : public OpenMM::Force {
public:
    /**
     * Create a ContForce object.
     */
    ContForce();
    /**
     * Get the number of bond stretch terms in the potential function
     */
    int getNumBonds() const {
        return bonds.size();
    }
    /**
     * Add a bond term to the force.
     *
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the force constant for the bond, measured in kJ/mol/nm^4
     * @return the index of the bond that was added
     */
    int addBond(std::vector<int> idxs , int npart, double length, double k);
    /**
     * Get the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to get parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond, measured in kJ/mol/nm^4
     */
    void getBondParameters(int index, std::vector<int>& loc_idxs, int& loc_npart, double& loc_length, double& loc_k) const;
    /**
     * Set the force field parameters for a bond term.
     * 
     * @param index     the index of the bond for which to set parameters
     * @param particle1 the index of the first particle connected by the bond
     * @param particle2 the index of the second particle connected by the bond
     * @param length    the equilibrium length of the bond, measured in nm
     * @param k         the harmonic force constant for the bond, measured in kJ/mol/nm^4
     */
    void setBondParameters(int index, std::vector<int> loc_idxs, int loc_npart, double loc_length, double loc_k);
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
     */
    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Returns true if the force uses periodic boundary conditions and false otherwise. Your force should implement this
     * method appropriately to ensure that `System.usesPeriodicBoundaryConditions()` works for all systems containing
     * your force.
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    class BondInfo;
    std::vector<BondInfo> bonds;
};

/**
 * This is an internal class used to record information about a bond.
 * @private
 */
class ContForce::BondInfo {
public:
    std::vector<int> idxs;
    int npart;
    double length, k;
    BondInfo() {
	  idxs.clear();
	  npart = 0;
	  length = k = 0.0;
    }
  BondInfo(std::vector<int> bi_idxs, int bi_npart, double bi_length, double bi_k) :
	idxs(bi_idxs), npart(bi_npart), length(bi_length), k(bi_k) {
  }
};

} // namespace ContForcePlugin

#endif /*OPENMM_CONTFORCE_H_*/
