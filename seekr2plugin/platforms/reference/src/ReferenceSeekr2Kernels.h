#ifndef REFERENCE_SEEKR2_KERNELS_H_
#define REFERENCE_SEEKR2_KERNELS_H_

/*
   Copyright 2019 by Lane Votapka
   All rights reserved
   
   -------------------------------------------------------------------------- *
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

#include "openmm/reference/ReferencePlatform.h"
#include "openmm/reference/ReferenceForce.h"
#include "openmm/reference/ReferenceStochasticDynamics.h"
#include "openmm/reference/RealVec.h"
#include "Seekr2Kernels.h"
#include "openmm/Platform.h"
#include <vector>
#include <map>

namespace Seekr2Plugin {

/**
 * This kernel is invoked by MmvtLangevinIntegrator to update the positions and velocities of the system
 based on the forces and will monitor for crossing events to reverse direction
 */

class ReferenceIntegrateMmvtLangevinStepKernel : public IntegrateMmvtLangevinStepKernel {
public:
    ReferenceIntegrateMmvtLangevinStepKernel(std::string name, const OpenMM::Platform& platform, OpenMM::ReferencePlatform::PlatformData& data) : IntegrateMmvtLangevinStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateMmvtLangevinStepKernel();
    /**
     * Initialize the kernel, setting up the particle masses.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the MmvtLangevinIntegrator this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const MmvtLangevinIntegrator& integrator);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the MmvtLangevinIntegrator this kernel is being used for
     */
    void execute(OpenMM::ContextImpl& context, const MmvtLangevinIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the MmvtLangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(OpenMM::ContextImpl& context, const MmvtLangevinIntegrator& integrator);
    
    
private:
    OpenMM::ReferencePlatform::PlatformData& data;
    OpenMM::ReferenceStochasticDynamics* dynamics;
    std::vector<double> masses;
    double prevTemp, prevFriction, prevStepSize;
    std::vector<OpenMM::Vec3> oldPosData;
    std::vector<OpenMM::Vec3> oldVelData;
    std::vector<OpenMM::Vec3> oldForceData;
    
    std::vector<int> N_alpha_beta;
    std::vector<std::vector <int> > Nij_alpha;
    std::vector<double> Ri_alpha;
    double T_alpha;
    std::string outputFileName;
    std::vector<int> milestoneGroups;
    bool saveStateBool = false;
    std::string saveStateFileName;
    bool saveStatisticsBool = false;
    std::string saveStatisticsFileName;
    std::vector<int> bitvector;
    std::vector<std::string> globalParameterNames;
    int numMilestoneGroups, bounceCounter, previousMilestoneCrossed;
    double firstCrossingTime;
    double incubationTime;
    
};

/**
 * This kernel is invoked by ElberLangevinIntegrator to update the positions and velocities of the system
 based on the forces and will monitor for crossing events to reverse direction
 */

class ReferenceIntegrateElberLangevinStepKernel : public IntegrateElberLangevinStepKernel {
public:
    ReferenceIntegrateElberLangevinStepKernel(std::string name, const OpenMM::Platform& platform, OpenMM::ReferencePlatform::PlatformData& data) : IntegrateElberLangevinStepKernel(name, platform),
        data(data), dynamics(0) {
    }
    ~ReferenceIntegrateElberLangevinStepKernel();
    /**
     * Initialize the kernel, setting up the particle masses.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the ElberLangevinIntegrator this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const ElberLangevinIntegrator& integrator);
    /**
     * Execute the kernel.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the ElberLangevinIntegrator this kernel is being used for
     */
    void execute(OpenMM::ContextImpl& context, const ElberLangevinIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the MmvtLangevinIntegrator this kernel is being used for
     */
    double computeKineticEnergy(OpenMM::ContextImpl& context, const ElberLangevinIntegrator& integrator);
    
private:
    OpenMM::ReferencePlatform::PlatformData& data;
    OpenMM::ReferenceStochasticDynamics* dynamics;
    std::vector<double> masses;
    double prevTemp, prevFriction, prevStepSize;
    
    std::string outputFileName;
    std::vector<int> srcbitvector;
    std::vector<int> destbitvector;
    std::vector<int> srcMilestoneGroups;
    std::vector<int> destMilestoneGroups;
    std::vector<double> srcMilestoneValues; // the force group potential energy "value" to detect crossings
    std::vector<double> destMilestoneValues;
    bool endOnSrcMilestone = true; // whether to end on one of the source milestone
    bool crossedSrcMilestone = false; // need to see if source milestone was crossed - only way to have valid statistics
    bool endSimulation = false; // If an ending milestone was crossed, then don't log any more crossings
    bool saveStateBool = false;
    std::string saveStateFileName;
    std::vector<std::string> globalParameterNames;
    int numSrcMilestoneGroups, numDestMilestoneGroups;
    int crossingCounter;
};

} // namespace Seekr2Plugin

#endif /*REFERENCE_SEEKR2_KERNELS_H_*/
