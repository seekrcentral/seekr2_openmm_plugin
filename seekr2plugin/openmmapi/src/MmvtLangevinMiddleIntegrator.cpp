/*
 * Copyright 2019 by Lane Votapka
 * All rights reserved
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

#include "MmvtLangevinMiddleIntegrator.h"
#include "Seekr2Kernels.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <string>
#include <iostream>

using namespace Seekr2Plugin;
using namespace OpenMM;
using std::string;
using std::vector;
using std::stringstream;
using std::map;
using std::set;
using std::cout;

MmvtLangevinMiddleIntegrator::MmvtLangevinMiddleIntegrator(double temperature, 
                          double frictionCoeff, double stepSize, 
                          string fileName) : temperature(temperature), 
                          friction(frictionCoeff) {
    setStepSize(stepSize);
    setRandomNumberSeed(0);
    setOutputFileName(fileName);
    setSaveStateFileName("");
    setConstraintTolerance(1e-5);
    setSaveStatisticsFileName("");
    setBounceCounter(0);
}

void MmvtLangevinMiddleIntegrator::initialize(ContextImpl& contextRef) {
    context = &contextRef;
    const System& system = contextRef.getSystem();
    kernel = context->getPlatform().createKernel(IntegrateMmvtLangevinMiddleStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateMmvtLangevinMiddleStepKernel>().initialize(contextRef.getSystem(), *this);
}

void MmvtLangevinMiddleIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> MmvtLangevinMiddleIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateMmvtLangevinMiddleStepKernel::Name());
    return names;
}

double MmvtLangevinMiddleIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateMmvtLangevinMiddleStepKernel>().computeKineticEnergy(*context, *this);
}

void MmvtLangevinMiddleIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");  
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false);
        kernel.getAs<IntegrateMmvtLangevinMiddleStepKernel>().execute(*context, *this);
    }
}

const string& MmvtLangevinMiddleIntegrator::getOutputFileName() const {
    return outputFileName;
}

void MmvtLangevinMiddleIntegrator::setOutputFileName(const string& fileName) {
    outputFileName = fileName;
}

const string& MmvtLangevinMiddleIntegrator::getSaveStateFileName() const {
    return saveStateFileName;
}

void MmvtLangevinMiddleIntegrator::setSaveStateFileName(const string& fileName) {
    saveStateFileName = fileName;
}

const string& MmvtLangevinMiddleIntegrator::getSaveStatisticsFileName() const {
    return saveStatisticsFileName;
}

void MmvtLangevinMiddleIntegrator::setSaveStatisticsFileName(const string& fileName) {
    saveStatisticsFileName = fileName;
}

const int MmvtLangevinMiddleIntegrator::getMilestoneGroup(int index) const {
    ASSERT_VALID_INDEX(index, milestoneGroups)
    return milestoneGroups[index];
}

int MmvtLangevinMiddleIntegrator::addMilestoneGroup(int milestoneGroup) {
    milestoneGroups.push_back(milestoneGroup);
    return milestoneGroups.size()-1;
}

const int MmvtLangevinMiddleIntegrator::getBounceCounter() const {
    return bounceCounter;
}

void MmvtLangevinMiddleIntegrator::setBounceCounter(int counter) {
    bounceCounter = counter;
}
