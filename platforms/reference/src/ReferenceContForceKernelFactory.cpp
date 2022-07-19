#include "ReferenceContForceKernelFactory.h"
#include "ReferenceContForceKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace ContForcePlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
        Platform& platform = Platform::getPlatform(i);
        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
            ReferenceContForceKernelFactory* factory = new ReferenceContForceKernelFactory();
            platform.registerKernelFactory(CalcContForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerExampleReferenceKernelFactories() {
    registerKernelFactories();
}

KernelImpl* ReferenceContForceKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& data = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    if (name == CalcContForceKernel::Name())
        return new ReferenceCalcContForceKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
