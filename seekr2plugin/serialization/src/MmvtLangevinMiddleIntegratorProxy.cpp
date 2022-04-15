/* Copyright 2019 by Lane Votapka
 * All rights reserved
 * -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "MmvtLangevinMiddleIntegratorProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "MmvtLangevinMiddleIntegrator.h"
#include <sstream>

using namespace OpenMM;
using namespace Seekr2Plugin;
using namespace std;

MmvtLangevinMiddleIntegratorProxy::MmvtLangevinMiddleIntegratorProxy() : SerializationProxy("MmvtLangevinMiddleIntegrator") {
}

void MmvtLangevinMiddleIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const MmvtLangevinMiddleIntegrator& integrator = *reinterpret_cast<const MmvtLangevinMiddleIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    node.setDoubleProperty("temperature", integrator.getTemperature());
    node.setDoubleProperty("friction", integrator.getFriction());
    node.setIntProperty("randomSeed", integrator.getRandomNumberSeed());
    node.setIntProperty("bounceCounter", integrator.getBounceCounter());
    node.setStringProperty("outputFileName", integrator.getOutputFileName());
    node.setStringProperty("saveStateFileName", integrator.getSaveStateFileName());
    node.setStringProperty("saveStatisticsFileName", integrator.getSaveStatisticsFileName());
    SerializationNode& perMilestoneGroups = node.createChildNode("milestoneGroups");
    for (int i = 0; i < integrator.getNumMilestoneGroups(); i++) {
        perMilestoneGroups.createChildNode("milestoneGroup").setIntProperty("forceGroupNumber", integrator.getMilestoneGroup(i));
    }
}

void* MmvtLangevinMiddleIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    MmvtLangevinMiddleIntegrator *integrator = new MmvtLangevinMiddleIntegrator(node.getDoubleProperty("temperature"),
            node.getDoubleProperty("friction"), node.getDoubleProperty("stepSize"), node.getStringProperty("outputFileName"));
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    integrator->setRandomNumberSeed(node.getIntProperty("randomSeed"));
    integrator->setBounceCounter(node.getIntProperty("bounceCounter"));
    integrator->setSaveStateFileName(node.getStringProperty("saveStateFileName"));
    integrator->setSaveStatisticsFileName(node.getStringProperty("saveStatisticsFileName"));
    const SerializationNode& perMilestoneGroups = node.getChildNode("milestoneGroups");
    for (auto& group : perMilestoneGroups.getChildren())
        integrator->addMilestoneGroup(group.getIntProperty("forceGroupNumber"));
    
    return integrator;
}
