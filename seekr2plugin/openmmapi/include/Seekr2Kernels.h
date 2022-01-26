#ifndef SEEKR2_KERNELS_H_
#define SEEKR2_KERNELS_H_

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


#include "MmvtLangevinIntegrator.h"
#include "ElberLangevinIntegrator.h"
#include "MmvtLangevinMiddleIntegrator.h"
#include "ElberLangevinMiddleIntegrator.h"
#include "openmm/KernelImpl.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <string>

namespace Seekr2Plugin {

/**
 * This kernel is invoked by MmvtLangevinIntegrator to take one time step.
 */
class IntegrateMmvtLangevinStepKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "IntegrateMmvtLangevinStep";
    }
    IntegrateMmvtLangevinStepKernel(std::string name, const OpenMM::Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the MmvtLangevinIntegrator this kernel will be used for
     * @param force      the Force to get particle parameters from
     */
    virtual void initialize(const OpenMM::System& system, const MmvtLangevinIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the MmvtLangevinIntegrator this kernel is being used for
     */
    virtual void execute(OpenMM::ContextImpl& context, const MmvtLangevinIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     */
    virtual double computeKineticEnergy(OpenMM::ContextImpl& context, const MmvtLangevinIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by ElberLangevinIntegrator to take one time step.
 */
class IntegrateElberLangevinStepKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "IntegrateElberLangevinStep";
    }
    IntegrateElberLangevinStepKernel(std::string name, const OpenMM::Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the ElberLangevinIntegrator this kernel will be used for
     * @param force      the Force to get particle parameters from
     */
    virtual void initialize(const OpenMM::System& system, const ElberLangevinIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the ElberLangevinIntegrator this kernel is being used for
     */
    virtual void execute(OpenMM::ContextImpl& context, const ElberLangevinIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     */
    virtual double computeKineticEnergy(OpenMM::ContextImpl& context, const ElberLangevinIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by MmvtLangevinMiddleIntegrator to take one time step.
 */
class IntegrateMmvtLangevinMiddleStepKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "IntegrateMmvtLangevinMiddleStep";
    }
    IntegrateMmvtLangevinMiddleStepKernel(std::string name, const OpenMM::Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the MmvtLangevinMiddleIntegrator this kernel will be used for
     * @param force      the Force to get particle parameters from
     */
    virtual void initialize(const OpenMM::System& system, const MmvtLangevinMiddleIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the MmvtLangevinMiddleIntegrator this kernel is being used for
     */
    virtual void execute(OpenMM::ContextImpl& context, const MmvtLangevinMiddleIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     */
    virtual double computeKineticEnergy(OpenMM::ContextImpl& context, const MmvtLangevinMiddleIntegrator& integrator) = 0;
};

/**
 * This kernel is invoked by ElberLangevinMiddleIntegrator to take one time step.
 */
class IntegrateElberLangevinMiddleStepKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "IntegrateElberLangevinMiddleStep";
    }
    IntegrateElberLangevinMiddleStepKernel(std::string name, const OpenMM::Platform& platform) : KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the ElberLangevinMiddleIntegrator this kernel will be used for
     * @param force      the Force to get particle parameters from
     */
    virtual void initialize(const OpenMM::System& system, const ElberLangevinMiddleIntegrator& integrator) = 0;
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the ElberLangevinMiddleIntegrator this kernel is being used for
     */
    virtual void execute(OpenMM::ContextImpl& context, const ElberLangevinMiddleIntegrator& integrator) = 0;
    /**
     * Compute the kinetic energy.
     */
    virtual double computeKineticEnergy(OpenMM::ContextImpl& context, const ElberLangevinMiddleIntegrator& integrator) = 0;
};

} // namespace Seekr2Plugin

#endif /*SEEKR2_KERNELS_H_*/
