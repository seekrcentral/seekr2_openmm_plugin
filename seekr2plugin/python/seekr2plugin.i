/*
   Copyright 2019 by Lane Votapka
   All rights reserved
*/

%module seekr2plugin

%include "factory.i"
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

%include "std_string.i"

%{
#include "MmvtLangevinIntegrator.h"
#include "ElberLangevinIntegrator.h"
#include "MmvtLangevinMiddleIntegrator.h"
#include "ElberLangevinMiddleIntegrator.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
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
try:
    import openmm as mm
    import openmm.unit as unit
except:
    import simtk.openmm as mm
    import simtk.unit as unit
    from testInstallation import testInstallation
%}

/*
 * Add units to function outputs.
*/

%apply std::vector<int> &INPUT { std::vector<int> & particles };
%apply std::vector<double> &INPUT { std::vector<double> & weights };

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getTemperature() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setTemperature(double temperature) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getFriction() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setFriction(double coeff) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getStepSize() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setStepSize(double size) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getConstraintTolerance() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setConstraintTolerance(double tol) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::step(int steps) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getOutputFileName() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setOutputFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getSaveStateFileName() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setSaveStateFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getSaveStatisticsFileName() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setSaveStatisticsFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getMilestoneGroup(int index) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getNumMilestoneGroups() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::addMilestoneGroup(int milestoneGroup) const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getRandomNumberSeed() const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setRandomNumberSeed(int seed) const %{
    
%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::getBounceCounter() const %{

%}

%pythonappend Seekr2Plugin::MmvtLangevinIntegrator::setBounceCounter(int counter) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getTemperature() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setTemperature(double temperature) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getFriction() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setFriction(double coeff) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getStepSize() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setStepSize(double size) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getConstraintTolerance() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setConstraintTolerance(double tol) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::step(int steps) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getOutputFileName() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setOutputFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getSaveStateFileName() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setSaveStateFileName(std::string fileName) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getSrcMilestoneGroup(int index) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getNumSrcMilestoneGroups() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::addSrcMilestoneGroup(int milestoneGroup) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getDestMilestoneGroup(int index) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getNumDestMilestoneGroups() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::addDestMilestoneGroup(int milestoneGroup) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getRandomNumberSeed() const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setRandomNumberSeed(int seed) const %{
    
%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getCrossingCounter() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setCrossingCounter(int counter) const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::getEndOnSrcMilestone() const %{

%}

%pythonappend Seekr2Plugin::ElberLangevinIntegrator::setEndOnSrcMilestone(bool endOnSrc) const %{

%}

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

namespace Seekr2Plugin {

class MmvtLangevinIntegrator : public OpenMM::Integrator {
public:
    MmvtLangevinIntegrator(double temperature, double frictionCoeff, 
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

class ElberLangevinIntegrator : public OpenMM::Integrator {
public:
    ElberLangevinIntegrator(double temperature, double frictionCoeff, 
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
