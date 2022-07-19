#ifdef WIN32
#include <windows.h>
#include <sstream>
#else
#include <dlfcn.h>
#include <dirent.h>
#include <cstdlib>
#endif

#include "ContForce.h"
#include "ContForceProxy.h"
#include "openmm/serialization/SerializationProxy.h"

#if defined(WIN32)
    #include <windows.h>
    extern "C" OPENMM_EXPORT_EXAMPLE void registerExampleSerializationProxies();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            registerExampleSerializationProxies();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) registerExampleSerializationProxies();
#endif

using namespace ContForcePlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT_EXAMPLE void registerExampleSerializationProxies() {
  SerializationProxy::registerProxy(typeid(ContForce), new ContForceProxy());
}
