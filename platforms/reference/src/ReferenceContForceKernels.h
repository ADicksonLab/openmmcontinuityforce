#ifndef REFERENCE_CONTFORCE_KERNELS_H_
#define REFERENCE_CONTFORCE_KERNELS_H_

#include "ContForceKernels.h"
#include "openmm/Platform.h"
#include <vector>

namespace ContForcePlugin {

class ReferenceCalcContForceKernel : public CalcContForceKernel {
public:
    ReferenceCalcContForceKernel(std::string name, const OpenMM::Platform& platform) : CalcContForceKernel(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the ContForce this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const ContForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ContForce to copy the parameters from
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const ContForce& force);
private:
    int numBonds;
    std::vector<std::vector<int>> idxs;
    std::vector<int> npart;
    std::vector<double> length, k;
};

} // namespace ContForcePlugin

#endif /*REFERENCE_CONTFORCE_KERNELS_H_*/
