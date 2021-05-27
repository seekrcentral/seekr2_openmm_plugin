/*
   Copyright 2019 by Lane Votapka
   All rights reserved
*/

%module seekr2plugin

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
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
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
try:
    import openmm as mm
    import openmm.unit as unit
except ImportError:
    import simtk.openmm as mm
    import simtk.unit as unit
%}

/*
 * Add units to function outputs.
*/

%apply std::vector<int> &INPUT { std::vector<int> & particles };
%apply std::vector<double> &INPUT { std::vector<double> & weights };

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

namespace Seekr2Plugin {

class MmvtLangevinIntegrator : public OpenMM::Integrator {
public:
    MmvtLangevinIntegrator(double temperature, double frictionCoeff, 
        double stepSize, std::string fileName);
    
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
