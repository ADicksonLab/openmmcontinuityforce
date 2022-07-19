#ifndef OPENMM_CONT_FORCE_PROXY_H_
#define OPENMM_CONT_FORCE_PROXY_H_

#include "internal/windowsExportExample.h"
#include "openmm/serialization/SerializationProxy.h"

namespace OpenMM {

/**
 * This is a proxy for serializing ContForce objects.
 */

  class OPENMM_EXPORT_EXAMPLE ContForceProxy : public SerializationProxy {
public:
    ContForceProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

} // namespace OpenMM

#endif /*OPENMM_CONT_FORCE_PROXY_H_*/
