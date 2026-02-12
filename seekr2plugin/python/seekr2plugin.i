/*
   Copyright 2019 by Lane Votapka
   All rights reserved
*/

%module seekr2plugin

%include "factory.i"
%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

namespace OpenMM {
    class Integrator;
}


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

%include "std_string.i"

%{
#include "MmvtLangevinMiddleIntegrator.h"
#include "ElberLangevinMiddleIntegrator.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"

using namespace OpenMM;
using namespace Seekr2Plugin;
%}

%exception {
    try {
        $action
    }
    catch (std::exception &e) {
        PyObject* mm = PyImport_AddModule("openmm");
        PyObject* openmm_exception = PyObject_GetAttrString(mm, "OpenMMException");
        PyErr_SetString(openmm_exception, const_cast<char*>(e.what()));
        return NULL;
    }
}

%pythoncode %{
import openmm
import openmm.app as openmm_app
import openmm.unit as unit

def testSeekr2Plugin():
    system = openmm.System()
    topology = openmm_app.Topology()
    mass = 1.0 * unit.dalton
    system.addParticle(mass)
    atom_name = "X0"
    if not atom_name in openmm_app.element.Element._elements_by_symbol:
        element = openmm_app.element.Element(0, atom_name, atom_name, mass)
    else:
        element = openmm_app.element.Element._elements_by_symbol[atom_name]
    chain = topology.addChain()
    residue = topology.addResidue("UNK", chain)
    topology.addAtom(atom_name, element, residue)
    integrator = MmvtLangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds, "/tmp/dummy.txt")
    simulation = openmm_app.Simulation(topology, system, integrator)
    # Successfully created a Ref simulation object and Context object.
    return
%}

/*
 * Add units to function outputs.
*/

%apply std::vector<int> &INPUT { std::vector<int> & particles };
%apply std::vector<double> &INPUT { std::vector<double> & weights };

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getTemperature() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setTemperature(double size) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getFriction() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setFriction(double coeff) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getStepSize() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setStepSize(double size) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getConstraintTolerance() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setConstraintTolerance(double tol) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::step(int steps) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getOutputFileName() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setOutputFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getSaveStateFileName() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setSaveStateFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getSaveStatisticsFileName() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setSaveStatisticsFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getMilestoneGroup(int index) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getNumMilestoneGroups() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::addMilestoneGroup(int milestoneGroup) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getRandomNumberSeed() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setRandomNumberSeed(int seed) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::getBounceCounter() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinMiddleIntegrator::setBounceCounter(int counter) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getTemperature() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setTemperature(double size) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getFriction() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setFriction(double coeff) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getStepSize() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setStepSize(double size) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getConstraintTolerance() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setConstraintTolerance(double tol) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::step(int steps) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getOutputFileName() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setOutputFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getSaveStateFileName() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setSaveStateFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getSrcMilestoneGroup(int index) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getNumSrcMilestoneGroups() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::addSrcMilestoneGroup(int milestoneGroup) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getDestMilestoneGroup(int index) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getNumDestMilestoneGroups() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::addDestMilestoneGroup(int milestoneGroup) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getRandomNumberSeed() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setRandomNumberSeed(int seed) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getCrossingCounter() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setCrossingCounter(int counter) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::getEndOnSrcMilestone() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinMiddleIntegrator::setEndOnSrcMilestone(bool endOnSrc) const %{

%}

namespace OpenMM {
    class Integrator;
}

namespace Seekr2Plugin {

class MmvtLangevinMiddleIntegrator : public OpenMM::Integrator {
public:
    MmvtLangevinMiddleIntegrator(double temperature, double frictionCoeff, 
        double stepSize, std::string fileName);
    
    double getTemperature() const;
    
    void setTemperature(double temperature);
    
    double getFriction() const;
    
    void setFriction(double coeff);
    
    double getStepSize() const;
    
    void setStepSize(double size);
    
    double getConstraintTolerance() const;
    
    void setConstraintTolerance(double tol);
    
    virtual void step(int steps);
    
    std::string getOutputFileName() const;
    
    void setOutputFileName(std::string fileName);
    
    std::string getSaveStateFileName() const;
    
    void setSaveStateFileName(std::string fileName);
    
    std::string getSaveStatisticsFileName() const;
    
    void setSaveStatisticsFileName(std::string fileName);
    
    int getMilestoneGroup(int index) const;
    
    int getNumMilestoneGroups() const;
    
    int addMilestoneGroup(int milestoneGroup);
    
    int getRandomNumberSeed() const;
    
    void setRandomNumberSeed(int seed);
    
    int getBounceCounter() const;
    
    void setBounceCounter(int counter);
};

class ElberLangevinMiddleIntegrator : public OpenMM::Integrator {
public:
    ElberLangevinMiddleIntegrator(double temperature, double frictionCoeff, 
        double stepSize, std::string fileName);
    
    double getTemperature() const;
    
    void setTemperature(double temperature);
    
    double getFriction() const;
    
    void setFriction(double coeff);
    
    double getStepSize() const;
    
    void setStepSize(double size);
    
    double getConstraintTolerance() const;
    
    void setConstraintTolerance(double tol);
    
    virtual void step(int steps);
    
    std::string getOutputFileName() const;
    
    void setOutputFileName(std::string fileName);
    
    std::string getSaveStateFileName() const;
    
    void setSaveStateFileName(std::string fileName);
    
    int getSrcMilestoneGroup(int index) const;
    
    int getNumSrcMilestoneGroups() const;
    
    int addSrcMilestoneGroup(int milestoneGroup);
    
    int getDestMilestoneGroup(int index) const;
    
    int getNumDestMilestoneGroups() const;
    
    int addDestMilestoneGroup(int milestoneGroup);
    
    int getRandomNumberSeed() const;
    
    void setRandomNumberSeed(int seed);
    
    int getCrossingCounter() const;
    
    void setCrossingCounter(int counter);
    
    bool getEndOnSrcMilestone() const;
    
    void setEndOnSrcMilestone(bool endOnSrc);
};

}
