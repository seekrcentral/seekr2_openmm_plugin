#ifndef OPENMM_MMVTLANGEVINMIDDLEINTEGRATOR_H_
#define OPENMM_MMVTLANGEVINMIDDLEINTEGRATOR_H_

/*
   Copyright 2019 by Lane Votapka
   All rights reserved
 * -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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

#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include "openmm/TabulatedFunction.h"
#include "openmm/Force.h"
#include "openmm/System.h"
#include "internal/windowsExportSeekr2.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include <utility>
#include <map>
#include <set>
#include <string>

namespace Seekr2Plugin {

/**
 * This is an Integrator which simulates a System using Langevin dynamics.
 * A modification has been made to enable MMVT calculations. This includes
 * returning to a previous step and reversing velocities.
 */

class OPENMM_EXPORT MmvtLangevinMiddleIntegrator : public OpenMM::Integrator {
public:
    /**
     * Create a MmvtLangevinMiddleIntegrator.
     * 
     * @param temperature    the temperature of the heat bath (in Kelvin)
     * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
     * @param stepSize       the step size with which to integrate the system (in picoseconds)
     * @param fileName       the name of the file to write milestone transitions to
     */
    MmvtLangevinMiddleIntegrator(double temperature, double frictionCoeff, double stepSize, std::string fileName);
    /**
     * Get the temperature of the heat bath (in Kelvin).
     *
     * @return the temperature of the heat bath, measured in Kelvin
     */
    double getTemperature() const {
        return temperature;
    }
    /**
     * Set the temperature of the heat bath (in Kelvin).
     *
     * @param temp    the temperature of the heat bath, measured in Kelvin
     */
    void setTemperature(double temp) {
        temperature = temp;
    }
    /**
     * Get the friction coefficient which determines how strongly the system is coupled to
     * the heat bath (in inverse ps).
     *
     * @return the friction coefficient, measured in 1/ps
     */
    double getFriction() const {
        return friction;
    }
    /**
     * Set the friction coefficient which determines how strongly the system is coupled to
     * the heat bath (in inverse ps).
     *
     * @param coeff    the friction coefficient, measured in 1/ps
     */
    void setFriction(double coeff) {
        friction = coeff;
    }
    /**
     * Get the random number seed.  See setRandomNumberSeed() for details.
     */
    int getRandomNumberSeed() const {
        return randomNumberSeed;
    }
    /**
     * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
     * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
     * are run with different random number seeds, the sequence of random forces will be different.  On
     * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
     * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
     * results on successive runs, even if those runs were initialized identically.
     *
     * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
     * is created from this Force. This is done to ensure that each Context receives unique random seeds
     * without you needing to set them explicitly.
     */
    void setRandomNumberSeed(int seed) {
        randomNumberSeed = seed;
    }
    /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    void step(int steps);
    
    /**
     * Get the file name that the integrator writes milestone transitions to
     */
    const std::string& getOutputFileName() const;
    
    /**
     * Set the file name that the integrator would write milestone transitions to
     *
     * @param fileName    the name of the output file
     */
    void setOutputFileName(const std::string& fileName);
    
    /**
     * Get the file name that the integrator writes positions/velocities to
     */
    const std::string& getSaveStateFileName() const;
    
    /**
     * Set the file name that the integrator would write positions/velocities to
     *
     * @param fileName    the string of the state file name upon crossing
     */
    void setSaveStateFileName(const std::string& fileName);
    
     /**
     * Get the base name of the file that the simulation state is being 
     * written to upon a bounce/crossing event
     */
    const std::string& getSaveStatisticsFileName() const;
    
    /**
     * Set the file name that the integrator would write precomputed
     * transition statistics to
     * @param fileName    the name of the output file to which transition 
     * statistics are written to
     */
    void setSaveStatisticsFileName(const std::string& filename);
    
     /**
     * Get the force group that describes a particular milestone
     * 
     * @param index    the index of the milestone whose force group to return
     */
    const int getMilestoneGroup(int index) const;
    
     /**
     * Get the number of milestone boundaries (force groups)
     */
    int getNumMilestoneGroups() const {
        return milestoneGroups.size();
    }
    
    /**
     * Add a force group for a milestone boundary
     *
     * @param milestoneGroup    the force group for a milestone
     */
    int addMilestoneGroup(const int milestoneGroup);
    
    const int getBounceCounter() const;
    
    void setBounceCounter(int counter);
    
protected:
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(OpenMM::ContextImpl& context);
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    void cleanup();
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();
private:
    double temperature, friction;
    int randomNumberSeed;
    OpenMM::Kernel kernel;
    std::string outputFileName;
    std::string saveStateFileName;
    std::string saveStatisticsFileName;
    std::vector<int> milestoneGroups;
    int bounceCounter;
    
};

} // namespace Seekr2Plugin

#endif /*OPENMM_MMVTLANGEVINMIDDLEINTEGRATOR_H_*/
