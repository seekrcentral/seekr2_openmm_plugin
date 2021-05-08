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

#include "ElberLangevinIntegratorProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "ElberLangevinIntegrator.h"
#include <sstream>

using namespace OpenMM;
using namespace Seekr2Plugin;
using namespace std;

ElberLangevinIntegratorProxy::ElberLangevinIntegratorProxy() : SerializationProxy("ElberLangevinIntegrator") {
}

void ElberLangevinIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const ElberLangevinIntegrator& integrator = *reinterpret_cast<const ElberLangevinIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    node.setDoubleProperty("temperature", integrator.getTemperature());
    node.setDoubleProperty("friction", integrator.getFriction());
    node.setIntProperty("randomSeed", integrator.getRandomNumberSeed());
    node.setStringProperty("outputFileName", integrator.getOutputFileName());
    node.setStringProperty("saveStateFileName", integrator.getSaveStateFileName());
    SerializationNode& perSrcMilestoneGroups = node.createChildNode("srcMilestoneGroups");
    for (int i = 0; i < integrator.getNumSrcMilestoneGroups(); i++) {
        perSrcMilestoneGroups.createChildNode("srcMilestoneGroup").setIntProperty("forceGroupNumber", integrator.getSrcMilestoneGroup(i));
    }
    SerializationNode& perDestMilestoneGroups = node.createChildNode("destMilestoneGroups");
    for (int i = 0; i < integrator.getNumDestMilestoneGroups(); i++) {
        perDestMilestoneGroups.createChildNode("destMilestoneGroup").setIntProperty("forceGroupNumber", integrator.getDestMilestoneGroup(i));
    }
}

void* ElberLangevinIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    ElberLangevinIntegrator *integrator = new ElberLangevinIntegrator(node.getDoubleProperty("temperature"),
            node.getDoubleProperty("friction"), node.getDoubleProperty("stepSize"), node.getStringProperty("outputFileName"));
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    integrator->setRandomNumberSeed(node.getIntProperty("randomSeed"));
    integrator->setSaveStateFileName(node.getStringProperty("saveStateFileName"));
    const SerializationNode& perSrcMilestoneGroups = node.getChildNode("srcMilestoneGroups");
    for (auto& group : perSrcMilestoneGroups.getChildren())
        integrator->addSrcMilestoneGroup(group.getIntProperty("forceGroupNumber"));
    const SerializationNode& perDestMilestoneGroups = node.getChildNode("destMilestoneGroups");
    for (auto& group : perDestMilestoneGroups.getChildren())
        integrator->addDestMilestoneGroup(group.getIntProperty("forceGroupNumber"));
    return integrator;
}
