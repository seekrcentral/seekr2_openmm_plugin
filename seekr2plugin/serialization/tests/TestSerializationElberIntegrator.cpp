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

#include "openmm/Platform.h"
#include "openmm/internal/AssertionUtilities.h"
#include "ElberLangevinIntegrator.h"
#include "ElberLangevinMiddleIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace Seekr2Plugin;
using namespace std;

extern "C" void registerElberSerializationProxies();

void testSerialization() {
    // Create an Integrator.

    ElberLangevinIntegrator integ1(301.1, 0.95, 0.001, "/tmp/dummy.txt");
    integ1.setRandomNumberSeed(18);
    integ1.addSrcMilestoneGroup(4);
    integ1.addDestMilestoneGroup(5);
    integ1.setSaveStateFileName("/tmp/dummyState.txt");

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<ElberLangevinIntegrator>(&integ1, "Integrator", buffer);
    ElberLangevinIntegrator* copy = XmlSerializer::deserialize<ElberLangevinIntegrator>(buffer);

    // Compare the two integrators to see if they are identical.
    
    ElberLangevinIntegrator& integ2 = *copy;
    ASSERT_EQUAL(integ1.getTemperature(), integ2.getTemperature());
    ASSERT_EQUAL(integ1.getFriction(), integ2.getFriction());
    ASSERT_EQUAL(integ1.getRandomNumberSeed(), integ2.getRandomNumberSeed());
    ASSERT_EQUAL(integ1.getSrcMilestoneGroup(0), integ2.getSrcMilestoneGroup(0));
    ASSERT_EQUAL(integ1.getDestMilestoneGroup(0), integ2.getDestMilestoneGroup(0));
    ASSERT_EQUAL(integ1.getSaveStateFileName(), integ2.getSaveStateFileName());
}

void testSerializationMiddle() {
    // Create an Integrator.

    ElberLangevinMiddleIntegrator integ1(301.1, 0.95, 0.001, "/tmp/dummy.txt");
    integ1.setRandomNumberSeed(18);
    integ1.addSrcMilestoneGroup(4);
    integ1.addDestMilestoneGroup(5);
    integ1.setSaveStateFileName("/tmp/dummyState.txt");

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<ElberLangevinMiddleIntegrator>(&integ1, "Integrator", buffer);
    ElberLangevinMiddleIntegrator* copy = XmlSerializer::deserialize<ElberLangevinMiddleIntegrator>(buffer);

    // Compare the two integrators to see if they are identical.
    
    ElberLangevinMiddleIntegrator& integ2 = *copy;
    ASSERT_EQUAL(integ1.getTemperature(), integ2.getTemperature());
    ASSERT_EQUAL(integ1.getFriction(), integ2.getFriction());
    ASSERT_EQUAL(integ1.getRandomNumberSeed(), integ2.getRandomNumberSeed());
    ASSERT_EQUAL(integ1.getSrcMilestoneGroup(0), integ2.getSrcMilestoneGroup(0));
    ASSERT_EQUAL(integ1.getDestMilestoneGroup(0), integ2.getDestMilestoneGroup(0));
    ASSERT_EQUAL(integ1.getSaveStateFileName(), integ2.getSaveStateFileName());
}

int main() {
    try {
        registerElberSerializationProxies();
        testSerialization();
        testSerializationMiddle();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

