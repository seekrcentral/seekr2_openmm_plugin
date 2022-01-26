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
#include "MmvtLangevinIntegrator.h"
#include "MmvtLangevinMiddleIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>

using namespace OpenMM;
using namespace Seekr2Plugin;
using namespace std;

extern "C" void registerMmvtSerializationProxies();

void testSerialization() {
    // Create an Integrator.

    MmvtLangevinIntegrator integ1(301.1, 0.95, 0.001, "/tmp/dummy.txt");
    integ1.setRandomNumberSeed(18);
    integ1.addMilestoneGroup(4);
    integ1.setSaveStateFileName("/tmp/dummyState.txt");
    integ1.setSaveStatisticsFileName("/tmp/dummyStatistics.txt");

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<MmvtLangevinIntegrator>(&integ1, "Integrator", buffer);
    MmvtLangevinIntegrator* copy = XmlSerializer::deserialize<MmvtLangevinIntegrator>(buffer);

    // Compare the two integrators to see if they are identical.
    
    MmvtLangevinIntegrator& integ2 = *copy;
    ASSERT_EQUAL(integ1.getTemperature(), integ2.getTemperature());
    ASSERT_EQUAL(integ1.getFriction(), integ2.getFriction());
    ASSERT_EQUAL(integ1.getRandomNumberSeed(), integ2.getRandomNumberSeed());
    ASSERT_EQUAL(integ1.getMilestoneGroup(0), integ2.getMilestoneGroup(0));
    ASSERT_EQUAL(integ1.getSaveStateFileName(), integ2.getSaveStateFileName());
    ASSERT_EQUAL(integ1.getSaveStatisticsFileName(), integ2.getSaveStatisticsFileName());
}

void testSerializationMiddle() {
    // Create an Integrator.

    MmvtLangevinMiddleIntegrator integ1(301.1, 0.95, 0.001, "/tmp/dummyMiddle.txt");
    integ1.setRandomNumberSeed(18);
    integ1.addMilestoneGroup(4);
    integ1.setSaveStateFileName("/tmp/dummyStateMiddle.txt");
    integ1.setSaveStatisticsFileName("/tmp/dummyStatisticsMiddle.txt");

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<MmvtLangevinMiddleIntegrator>(&integ1, "Integrator", buffer);
    MmvtLangevinMiddleIntegrator* copy = XmlSerializer::deserialize<MmvtLangevinMiddleIntegrator>(buffer);

    // Compare the two integrators to see if they are identical.
    
    MmvtLangevinMiddleIntegrator& integ2 = *copy;
    ASSERT_EQUAL(integ1.getTemperature(), integ2.getTemperature());
    ASSERT_EQUAL(integ1.getFriction(), integ2.getFriction());
    ASSERT_EQUAL(integ1.getRandomNumberSeed(), integ2.getRandomNumberSeed());
    ASSERT_EQUAL(integ1.getMilestoneGroup(0), integ2.getMilestoneGroup(0));
    ASSERT_EQUAL(integ1.getSaveStateFileName(), integ2.getSaveStateFileName());
    ASSERT_EQUAL(integ1.getSaveStatisticsFileName(), integ2.getSaveStatisticsFileName());
}

int main() {
    try {
        registerMmvtSerializationProxies();
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

