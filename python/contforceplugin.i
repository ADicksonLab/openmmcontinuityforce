%module contforceplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "ContForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Add units to function outputs.

%pythonappend ContForcePlugin::ContForce::getBondParameters(int index, std::vector<int>& idxs, int& npart,
                                                             double& length, double& k) const %{
    val[2] = unit.Quantity(val[2], unit.nanometer)
    val[3] = unit.Quantity(val[3], unit.kilojoule_per_mole/unit.nanometer**4)
%}


 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


namespace ContForcePlugin {

class ContForce : public OpenMM::Force {
public:
    ContForce();

    int getNumBonds() const;

    int addBond(std::vector<int> idxs, int npart, double length, double k);

    void setBondParameters(int index, std::vector<int> idxs, int npart, double length, double k);

    void updateParametersInContext(OpenMM::Context& context);

    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    %apply std::vector<int>& OUTPUT {std::vector<int>& idxs};
    %apply int& OUTPUT {int& npart};
    %apply double& OUTPUT {double& length};
    %apply double& OUTPUT {double& k};
    void getBondParameters(int index, std::vector<int>& idxs, int& npart, double& length, double& k) const;
    %clear std::vector<int>& idxs;
    %clear int& npart;
    %clear double& length;
    %clear double& k;

    /*
     * Add methods for casting a Force to an ContForce.
    */
    %extend {
        static ContForcePlugin::ContForce& cast(OpenMM::Force& force) {
            return dynamic_cast<ContForcePlugin::ContForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<ContForcePlugin::ContForce*>(&force) != NULL);
        }
    }
};

}
