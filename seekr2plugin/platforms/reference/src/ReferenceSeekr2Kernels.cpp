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


#include "ReferenceSeekr2Kernels.h"
#include "MmvtLangevinIntegrator.h"
#include "ElberLangevinIntegrator.h"
#include "MmvtLangevinMiddleIntegrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/reference/ReferenceForce.h"
#include "openmm/Context.h"
#include "openmm/reference/SimTKOpenMMUtilities.h"
#include "openmm/reference/ReferenceConstraints.h"
#include "openmm/reference/ReferenceVirtualSites.h"
#include "openmm/reference/ReferenceTabulatedFunction.h"
#include "openmm/State.h"
#include "openmm/serialization/XmlSerializer.h"
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace Seekr2Plugin;
using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->velocities);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(ReferenceConstraints*) data->constraints;
}
/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();
    
    // Compute the shifted velocities.
    
    vector<Vec3> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        if (masses[i] > 0)
            shiftedVel[i] = velData[i]+forceData[i]*(timeShift/masses[i]);
        else
            shiftedVel[i] = velData[i];
    }
    
    // Apply constraints to them.
    
    vector<double> inverseMasses(numParticles);
    for (int i = 0; i < numParticles; i++)
        inverseMasses[i] = (masses[i] == 0 ? 0 : 1/masses[i]);
    extractConstraints(context).applyToVelocities(posData, shiftedVel, inverseMasses, 1e-4);
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}

ReferenceIntegrateMmvtLangevinStepKernel::~ReferenceIntegrateMmvtLangevinStepKernel() {
    
}

void ReferenceIntegrateMmvtLangevinStepKernel::initialize(const System& system, const MmvtLangevinIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    bool output_file_already_exists;
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
    
    N_alpha_beta = vector<int> (integrator.getNumMilestoneGroups());
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        N_alpha_beta[i] = 0;
    }
    Nij_alpha = vector<vector <int> > (integrator.getNumMilestoneGroups());
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        Nij_alpha[i] = vector<int>(integrator.getNumMilestoneGroups());
        for (int j=0; j<integrator.getNumMilestoneGroups(); j++) {
            Nij_alpha[i][j] = 0;
        }
    }
    Ri_alpha = vector<double> (integrator.getNumMilestoneGroups());
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        Ri_alpha[i] = 0.0;
    }
    T_alpha = 0.0;
    outputFileName = integrator.getOutputFileName();
    
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        // TODO: marked for removal
        /*
        bool foundForceGroup=false;
        for (int j=0; j<system.getNumForces(); j++) {
            if (system.getForce(j).getForceGroup() == integrator.getMilestoneGroup(i))
                foundForceGroup=true;
        }
        if (foundForceGroup == false)
            throw OpenMMException("System contains no force groups used to detect MMVT boundary crossings. Check for mismatches between force group assignments and the groups added to the MMVT integrator.");
        */
        milestoneGroups.push_back(integrator.getMilestoneGroup(i));
    }
    saveStateFileName = integrator.getSaveStateFileName();
    if (saveStateFileName.empty()) {
        saveStateBool = false;
    } else {
        saveStateBool = true;
    }
    saveStatisticsFileName = integrator.getSaveStatisticsFileName();
    if (saveStatisticsFileName.empty()) {
        saveStatisticsBool = false;
    } else {
        saveStatisticsBool = true;
    }
    bounceCounter = integrator.getBounceCounter();
    incubationTime = 0.0;
    firstCrossingTime = 0.0;
    previousMilestoneCrossed = -1;
    assert(data.stepCount == 0);
    assert(data.time == 0.0);
    ofstream datafile; // open datafile for writing
    datafile.open(outputFileName);
    if (datafile) {
        output_file_already_exists = true;
    } else {
        output_file_already_exists = false;
    }
    datafile.close();
    if (output_file_already_exists == false) {
        ofstream datafile; // open datafile for writing
        datafile.open(outputFileName, std::ios_base::app); // write new file
        datafile << "#\"Bounced boundary ID\",\"bounce index\",\"total time (ps)\"\n";
        datafile.close(); // close data file
    }
}

void ReferenceIntegrateMmvtLangevinStepKernel::execute(ContextImpl& context, const MmvtLangevinIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    bool includeForces = false;
    bool includeEnergy = true;
    double value = 0.0; // The value to monitor for crossing events
    
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3> oldPosData(posData); // this will create an exact copy of old positions
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3> oldVelData(velData); // this will create an exact copy of old velocities
    vector<Vec3>& forceData = extractForces(context);
    vector<Vec3> oldForceData(forceData); // this will create an exact copy of old forces
    
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    
    
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        if (dynamics) {
            delete dynamics;
        }
        dynamics = new ReferenceStochasticDynamics(
                context.getSystem().getNumParticles(), 
                stepSize, 
                friction, 
                temperature);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    
    if (data.stepCount <= 0) {
        value = context.calcForcesAndEnergy(includeForces, includeEnergy, 2);
        if (value > 0.0) {
            throw OpenMMException("MMVT simulation bouncing on first step: the system is trapped behind a boundary. Check and revise MMVT boundary definitions and atomic positions.");
        }
    }
    
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    // EXTRACT POSITIONS HERE AND TEST FOR CRITERIA
    // NOTE: context positions and velocities are changed by reference
    // test if criteria satisfied
    // then reverse velocities by reference
    // restore old positions
    int bitcode;
    bool bounced = false;
    int num_bounced_surfaces = 0;
    // Handle bit string here
    value = context.calcForcesAndEnergy(includeForces, includeEnergy, 2);
    if (value > 0.0) { // take a step back and reverse velocities
        bounced = true;
        //if (data.stepCount <= 0)
    	//    throw OpenMMException("MMVT simulation bouncing on first step: the system is trapped behind a boundary. Check and revise MMVT boundary definitions and atomic positions.");
        
        ofstream datafile; // open datafile for writing
        datafile.open(outputFileName, std::ios_base::app); // append to file
        datafile.setf(std::ios::fixed,std::ios::floatfield);
        datafile.precision(3);
        bitcode = static_cast<int>(value);
        // check for corner bounce so as not to save state
        for (int i=0; i<milestoneGroups.size(); i++) {
            if ((bitcode % 2) != 0) {
                // count the number of bounced surfaces
                num_bounced_surfaces++;
            }
            bitcode = bitcode >> 1;
        }
        bitcode = static_cast<int>(value);
        for (int i=0; i<milestoneGroups.size(); i++) {
            // Write to output file
            if ((bitcode % 2) == 0) {
                // if value is not showing this as crossed, continue
                bitcode = bitcode >> 1;
                continue;
            }
            bitcode = bitcode >> 1;
            datafile << milestoneGroups[i] << "," << bounceCounter << ","<< context.getTime() << "\n";
            if (saveStateBool == true && num_bounced_surfaces == 1) {
                State myState = context.getOwner().getState(State::Positions | State::Velocities);
                stringstream buffer;
                stringstream number_str;
                number_str << "_" << bounceCounter << "_" << milestoneGroups[i] ;
                string trueFileName = saveStateFileName + number_str.str();
                XmlSerializer::serialize<State>(&myState, "State", buffer);
                ofstream statefile; // open datafile for writing
                statefile.open(trueFileName, std::ios_base::out); // append to file
                statefile << buffer.rdbuf();
                statefile.close(); // close data file
            }
            
            N_alpha_beta[i] += 1;
            if (previousMilestoneCrossed != i) {
                if (previousMilestoneCrossed != -1) { // if this isn't the first time a bounce has occurred
                    Nij_alpha[previousMilestoneCrossed][i] += 1; 
                    Ri_alpha[previousMilestoneCrossed] += incubationTime;
                } else {
                    firstCrossingTime = data.time;
                }
                incubationTime = 0.0;
            }
            T_alpha = data.time - firstCrossingTime;
            
            previousMilestoneCrossed = i;
            bounceCounter++;
        }
        datafile.close(); // close data file
        if (saveStatisticsBool == true) {
            ofstream stats;
            stats.open(saveStatisticsFileName, ios_base::trunc);
            stats.setf(std::ios::fixed,std::ios::floatfield);
            stats.precision(3);
            for (int i = 0; i < N_alpha_beta.size(); i++) {
                stats << "N_alpha_" << milestoneGroups[i] << ": " << N_alpha_beta[i] << "\n";
            }
            int nij_index = 0;
            for (int i = 0; i < integrator.getNumMilestoneGroups(); i++) {
                for (int j = 0; j < integrator.getNumMilestoneGroups(); j++) {
                    stats << "N_" << milestoneGroups[i] << "_" << milestoneGroups[j] << "_alpha: " << Nij_alpha[i][j] << "\n"; 
                    nij_index += 1;
                }
            }
            for (int i = 0; i < integrator.getNumMilestoneGroups(); i++) {
                stats << "R_" << milestoneGroups[i] << "_alpha: " << Ri_alpha[i] << "\n";       
            }
            stats << "T_alpha: " << T_alpha << "\n";  
            stats.close();
        }
    }
    if (bounced == true) {
        posData = oldPosData;
        for (int j=0; j<velData.size(); j++) {
            velData[j] = oldVelData[j] * -1.0; // take a step back and reverse the velocities of the particles
        }
        forceData = oldForceData;
    }
    
    data.time += stepSize;
    data.stepCount++;
    incubationTime += stepSize;
}

double ReferenceIntegrateMmvtLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const MmvtLangevinIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ReferenceIntegrateElberLangevinStepKernel::~ReferenceIntegrateElberLangevinStepKernel() {
    
}

void ReferenceIntegrateElberLangevinStepKernel::initialize(const System& system, const ElberLangevinIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    bool output_file_already_exists;
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
    
    outputFileName = integrator.getOutputFileName();
    endOnSrcMilestone = integrator.getEndOnSrcMilestone();
    for (int i=0; i<integrator.getNumSrcMilestoneGroups(); i++) {
        bool foundForceGroup=false;
        for (int j=0; j<system.getNumForces(); j++) {
            if (system.getForce(j).getForceGroup() == integrator.getSrcMilestoneGroup(i))
                foundForceGroup=true;
        }
        if (foundForceGroup == false)
            throw OpenMMException("System contains no force groups used to detect Elber boundary crossings. Check for mismatches between force group assignments and the groups added to the MMVT integrator.");
        srcMilestoneGroups.push_back(integrator.getSrcMilestoneGroup(i));
        srcMilestoneValues.push_back(-INFINITY);
    }
    for (int i=0; i<integrator.getNumDestMilestoneGroups(); i++) {
        bool foundForceGroup=false;
        for (int j=0; j<system.getNumForces(); j++) {
            if (system.getForce(j).getForceGroup() == integrator.getDestMilestoneGroup(i))
                foundForceGroup=true;
        }
        if (foundForceGroup == false)
            throw OpenMMException("System contains no force groups used to detect Elber boundary crossings. Check for mismatches between force group assignments and the groups added to the MMVT integrator.");
        destMilestoneGroups.push_back(integrator.getDestMilestoneGroup(i));
        destMilestoneValues.push_back(-INFINITY);
    }
    
    saveStateFileName = integrator.getSaveStateFileName();
    if (saveStateFileName.empty()) {
        saveStateBool = false;
    } else {
        saveStateBool = true;
    }
    
    srcbitvector.clear();
    for (int i=0; i<srcMilestoneGroups.size(); i++) {
        srcbitvector.push_back(pow(2,srcMilestoneGroups[i]));
    }
    destbitvector.clear();
    for (int i=0; i<destMilestoneGroups.size(); i++) {
        destbitvector.push_back(pow(2,destMilestoneGroups[i]));
    }
    crossingCounter = integrator.getCrossingCounter();
    assert(data.stepCount == 0);
    assert(data.time == 0.0);
    ifstream datafile; // open datafile for writing
    datafile.open(outputFileName);
    if (datafile) {
        output_file_already_exists = true;
    } else {
        output_file_already_exists = false;
    }
    datafile.close();
    if (output_file_already_exists == false) {
        ofstream datafile; // open datafile for writing
        datafile.open(outputFileName, std::ios_base::app); // write new file
        datafile << "#\"Crossed boundary ID\",\"crossing counter\",\"total time (ps)\"\n";
        datafile << "# An asterisk(*) indicates that source milestone was never crossed - asterisked statistics are invalid and should be excluded.\n";
        datafile.close(); // close data file
    }
}

void ReferenceIntegrateElberLangevinStepKernel::execute(ContextImpl& context, const ElberLangevinIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);    
    
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    
    
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        if (dynamics) {
            delete dynamics;
        }
        dynamics = new ReferenceStochasticDynamics(
                context.getSystem().getNumParticles(), 
                stepSize, 
                friction, 
                temperature);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    // EXTRACT POSITIONS HERE AND TEST FOR CRITERIA
    // NOTE: context positions and velocities are changed by reference
    // test if criteria satisfied
    // then reverse velocities by reference
    // restore old positions
    
    float value = 0.0;
    float oldvalue = 0.0;
    bool includeForces = false;
    bool includeEnergy = true;
    int endMilestoneGroup = 0;
    int num_bounced_surfaces = 0;
    if (endSimulation == false) {
        // first check source milestone crossings
        for (int i=0; i<integrator.getNumSrcMilestoneGroups(); i++) {
            value = context.calcForcesAndEnergy(includeForces, includeEnergy, srcbitvector[i]);
            if (srcMilestoneValues[i] == -INFINITY) {
                // First timestep
                srcMilestoneValues[i] = value;
            }
            oldvalue = srcMilestoneValues[i];
            if ((value - oldvalue) != 0.0) {
                // The source milestone has been crossed
                if (endOnSrcMilestone == true) {
                    endSimulation = true;
                    num_bounced_surfaces++;
                    ofstream datafile;
                    datafile.open(outputFileName, std::ios_base::app);
                    datafile << integrator.getSrcMilestoneGroup(i) << "," << crossingCounter << "," << context.getTime() << "\n";
                    endMilestoneGroup = integrator.getSrcMilestoneGroup(i);
                    datafile.close();
                } else {
                    crossedSrcMilestone = true;
                    context.setTime(0.0); // reset the timer
                    srcMilestoneValues[i] = value;
                }
            } 
        }
        
        // then check destination milestone crossings
        for (int i=0; i<integrator.getNumDestMilestoneGroups(); i++) {
            value = context.calcForcesAndEnergy(includeForces, includeEnergy, destbitvector[i]);
            if (destMilestoneValues[i] == -INFINITY) {
                // First timestep
                destMilestoneValues[i] = value;
            }
            oldvalue = destMilestoneValues[i];
            if ((value - oldvalue) != 0.0) {
                // The destination milestone has been crossed
                endSimulation = true;
                num_bounced_surfaces++;
                ofstream datafile;
                datafile.open(outputFileName, std::ios_base::app);
                if ((crossedSrcMilestone == true) || (endOnSrcMilestone == true)) {
                    datafile << integrator.getDestMilestoneGroup(i) << "," << crossingCounter << "," << context.getTime() << "\n";
                } else {
                    datafile << integrator.getDestMilestoneGroup(i) << "*," << crossingCounter << "," << context.getTime() << "\n";
                }
                endMilestoneGroup = integrator.getDestMilestoneGroup(i);
                datafile.close();
            } 
        }
        
        if (endSimulation == true) {
            // Then a crossing event has just occurred.
            if (saveStateBool == true && num_bounced_surfaces == 1) {
                State myState = context.getOwner().getState(State::Positions | State::Velocities);
                stringstream buffer;
                stringstream number_str;
                //number_str << "_" << crossingCounter << "_" << crossingCounter;
                number_str << "_" << crossingCounter << "_" << endMilestoneGroup;
                string trueFileName = saveStateFileName + number_str.str();
                XmlSerializer::serialize<State>(&myState, "State", buffer);
                ofstream statefile; // open datafile for writing
                statefile.open(trueFileName, std::ios_base::trunc); // append to file
                statefile << buffer.rdbuf();
                statefile.close(); // close data file
            }
            crossingCounter ++;
        }
    }
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateElberLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const ElberLangevinIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ReferenceIntegrateMmvtLangevinMiddleStepKernel::~ReferenceIntegrateMmvtLangevinMiddleStepKernel() {
    
}

void ReferenceIntegrateMmvtLangevinMiddleStepKernel::initialize(const System& system, const MmvtLangevinMiddleIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    bool output_file_already_exists;
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
    
    N_alpha_beta = vector<int> (integrator.getNumMilestoneGroups());
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        N_alpha_beta[i] = 0;
    }
    Nij_alpha = vector<vector <int> > (integrator.getNumMilestoneGroups());
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        Nij_alpha[i] = vector<int>(integrator.getNumMilestoneGroups());
        for (int j=0; j<integrator.getNumMilestoneGroups(); j++) {
            Nij_alpha[i][j] = 0;
        }
    }
    Ri_alpha = vector<double> (integrator.getNumMilestoneGroups());
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        Ri_alpha[i] = 0.0;
    }
    T_alpha = 0.0;
    outputFileName = integrator.getOutputFileName();
    
    for (int i=0; i<integrator.getNumMilestoneGroups(); i++) {
        milestoneGroups.push_back(integrator.getMilestoneGroup(i));
    }
    saveStateFileName = integrator.getSaveStateFileName();
    if (saveStateFileName.empty()) {
        saveStateBool = false;
    } else {
        saveStateBool = true;
    }
    saveStatisticsFileName = integrator.getSaveStatisticsFileName();
    if (saveStatisticsFileName.empty()) {
        saveStatisticsBool = false;
    } else {
        saveStatisticsBool = true;
    }
    bounceCounter = integrator.getBounceCounter();
    incubationTime = 0.0;
    firstCrossingTime = 0.0;
    previousMilestoneCrossed = -1;
    assert(data.stepCount == 0);
    assert(data.time == 0.0);
    ofstream datafile; // open datafile for writing
    datafile.open(outputFileName);
    if (datafile) {
        output_file_already_exists = true;
    } else {
        output_file_already_exists = false;
    }
    datafile.close();
    if (output_file_already_exists == false) {
        ofstream datafile; // open datafile for writing
        datafile.open(outputFileName, std::ios_base::app); // write new file
        datafile << "#\"Bounced boundary ID\",\"bounce index\",\"total time (ps)\"\n";
        datafile.close(); // close data file
    }
}

void ReferenceIntegrateMmvtLangevinMiddleStepKernel::execute(ContextImpl& context, const MmvtLangevinMiddleIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    bool includeForces = false;
    bool includeEnergy = true;
    double value = 0.0; // The value to monitor for crossing events
    
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3> oldPosData(posData); // this will create an exact copy of old positions
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3> oldVelData(velData); // this will create an exact copy of old velocities
    vector<Vec3>& forceData = extractForces(context);
    vector<Vec3> oldForceData(forceData); // this will create an exact copy of old forces
    
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    
    
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        if (dynamics) {
            delete dynamics;
        }
        dynamics = new ReferenceLangevinMiddleDynamics(
                context.getSystem().getNumParticles(), 
                stepSize, 
                friction, 
                temperature);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    
    if (data.stepCount <= 0) {
        value = context.calcForcesAndEnergy(includeForces, includeEnergy, 2);
        if (value > 0.0) {
            throw OpenMMException("MMVT simulation bouncing on first step: the system is trapped behind a boundary. Check and revise MMVT boundary definitions and atomic positions.");
        }
    }
    
    dynamics->update(context, posData, velData, masses, integrator.getConstraintTolerance());
    // EXTRACT POSITIONS HERE AND TEST FOR CRITERIA
    // NOTE: context positions and velocities are changed by reference
    // test if criteria satisfied
    // then reverse velocities by reference
    // restore old positions
    int bitcode;
    bool bounced = false;
    int num_bounced_surfaces = 0;
    // Handle bit string here
    value = context.calcForcesAndEnergy(includeForces, includeEnergy, 2);
    if (value > 0.0) { // take a step back and reverse velocities
        ofstream datafile; // open datafile for writing
        datafile.open(outputFileName, std::ios_base::app); // append to file
        datafile.setf(std::ios::fixed,std::ios::floatfield);
        datafile.precision(3);
        bitcode = static_cast<int>(value);
        // check for corner bounce so as not to save state
        for (int i=0; i<milestoneGroups.size(); i++) {
            if ((bitcode % 2) != 0) {
                // count the number of bounced surfaces
                num_bounced_surfaces++;
            }
            bitcode = bitcode >> 1;
        }
        bitcode = static_cast<int>(value);
        for (int i=0; i<milestoneGroups.size(); i++) {
            bounced = true;
            // Write to output file
            if ((bitcode % 2) == 0) {
                // if value is not showing this as crossed, continue
                bitcode = bitcode >> 1;
                continue;
            }
            bitcode = bitcode >> 1;
            datafile << milestoneGroups[i] << "," << bounceCounter << ","<< context.getTime() << "\n";
            if (saveStateBool == true && num_bounced_surfaces == 1) {
                State myState = context.getOwner().getState(State::Positions | State::Velocities);
                stringstream buffer;
                stringstream number_str;
                number_str << "_" << bounceCounter << "_" << milestoneGroups[i] ;
                string trueFileName = saveStateFileName + number_str.str();
                XmlSerializer::serialize<State>(&myState, "State", buffer);
                ofstream statefile; // open datafile for writing
                statefile.open(trueFileName, std::ios_base::out); // append to file
                statefile << buffer.rdbuf();
                statefile.close(); // close data file
            }
            
            N_alpha_beta[i] += 1;
            if (previousMilestoneCrossed != i) {
                if (previousMilestoneCrossed != -1) { // if this isn't the first time a bounce has occurred
                    Nij_alpha[previousMilestoneCrossed][i] += 1; 
                    Ri_alpha[previousMilestoneCrossed] += incubationTime;
                } else {
                    firstCrossingTime = data.time;
                }
                incubationTime = 0.0;
            }
            T_alpha = data.time - firstCrossingTime;
            
            previousMilestoneCrossed = i;
            bounceCounter++;
        }
        datafile.close(); // close data file
        if (saveStatisticsBool == true) {
            ofstream stats;
            stats.open(saveStatisticsFileName, ios_base::trunc);
            stats.setf(std::ios::fixed,std::ios::floatfield);
            stats.precision(3);
            for (int i = 0; i < N_alpha_beta.size(); i++) {
                stats << "N_alpha_" << milestoneGroups[i] << ": " << N_alpha_beta[i] << "\n";
            }
            int nij_index = 0;
            for (int i = 0; i < integrator.getNumMilestoneGroups(); i++) {
                for (int j = 0; j < integrator.getNumMilestoneGroups(); j++) {
                    stats << "N_" << milestoneGroups[i] << "_" << milestoneGroups[j] << "_alpha: " << Nij_alpha[i][j] << "\n"; 
                    nij_index += 1;
                }
            }
            for (int i = 0; i < integrator.getNumMilestoneGroups(); i++) {
                stats << "R_" << milestoneGroups[i] << "_alpha: " << Ri_alpha[i] << "\n";       
            }
            stats << "T_alpha: " << T_alpha << "\n";  
            stats.close();
        }
    }
    if (bounced == true) {
        posData = oldPosData;
        for (int j=0; j<velData.size(); j++) {
            velData[j] = oldVelData[j] * -1.0; // take a step back and reverse the velocities of the particles
        }
        forceData = oldForceData;
    }

    data.time += stepSize;
    data.stepCount++;
    incubationTime += stepSize;
}

double ReferenceIntegrateMmvtLangevinMiddleStepKernel::computeKineticEnergy(ContextImpl& context, const MmvtLangevinMiddleIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ReferenceIntegrateElberLangevinMiddleStepKernel::~ReferenceIntegrateElberLangevinMiddleStepKernel() {
    
}

void ReferenceIntegrateElberLangevinMiddleStepKernel::initialize(const System& system, const ElberLangevinMiddleIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    bool output_file_already_exists;
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
    
    outputFileName = integrator.getOutputFileName();
    endOnSrcMilestone = integrator.getEndOnSrcMilestone();
    for (int i=0; i<integrator.getNumSrcMilestoneGroups(); i++) {
        bool foundForceGroup=false;
        for (int j=0; j<system.getNumForces(); j++) {
            if (system.getForce(j).getForceGroup() == integrator.getSrcMilestoneGroup(i))
                foundForceGroup=true;
        }
        if (foundForceGroup == false)
            throw OpenMMException("System contains no force groups used to detect Elber boundary crossings. Check for mismatches between force group assignments and the groups added to the MMVT integrator.");
        srcMilestoneGroups.push_back(integrator.getSrcMilestoneGroup(i));
        srcMilestoneValues.push_back(-INFINITY);
    }
    for (int i=0; i<integrator.getNumDestMilestoneGroups(); i++) {
        bool foundForceGroup=false;
        for (int j=0; j<system.getNumForces(); j++) {
            if (system.getForce(j).getForceGroup() == integrator.getDestMilestoneGroup(i))
                foundForceGroup=true;
        }
        if (foundForceGroup == false)
            throw OpenMMException("System contains no force groups used to detect Elber boundary crossings. Check for mismatches between force group assignments and the groups added to the MMVT integrator.");
        destMilestoneGroups.push_back(integrator.getDestMilestoneGroup(i));
        destMilestoneValues.push_back(-INFINITY);
    }
    
    saveStateFileName = integrator.getSaveStateFileName();
    if (saveStateFileName.empty()) {
        saveStateBool = false;
    } else {
        saveStateBool = true;
    }
    
    srcbitvector.clear();
    for (int i=0; i<srcMilestoneGroups.size(); i++) {
        srcbitvector.push_back(pow(2,srcMilestoneGroups[i]));
    }
    destbitvector.clear();
    for (int i=0; i<destMilestoneGroups.size(); i++) {
        destbitvector.push_back(pow(2,destMilestoneGroups[i]));
    }
    crossingCounter = integrator.getCrossingCounter();
    assert(data.stepCount == 0);
    assert(data.time == 0.0);
    ofstream datafile; // open datafile for writing
    datafile.open(outputFileName);
    if (datafile) {
        output_file_already_exists = true;
    } else {
        output_file_already_exists = false;
    }
    datafile.close();
    if (output_file_already_exists == false) {
        ofstream datafile; // open datafile for writing
        datafile.open(outputFileName, std::ios_base::app); // write new file
        datafile << "#\"Crossed boundary ID\",\"crossing counter\",\"total time (ps)\"\n";
        datafile << "# An asterisk(*) indicates that source milestone was never crossed - asterisked statistics are invalid and should be excluded.\n";
        datafile.close(); // close data file
    }
}

void ReferenceIntegrateElberLangevinMiddleStepKernel::execute(ContextImpl& context, const ElberLangevinMiddleIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);    
    
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    
    
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        if (dynamics) {
            delete dynamics;
        }
        dynamics = new ReferenceLangevinMiddleDynamics(
                context.getSystem().getNumParticles(), 
                stepSize, 
                friction, 
                temperature);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    
    dynamics->update(context, posData, velData, masses, integrator.getConstraintTolerance());
    // EXTRACT POSITIONS HERE AND TEST FOR CRITERIA
    // NOTE: context positions and velocities are changed by reference
    // test if criteria satisfied
    // then reverse velocities by reference
    // restore old positions
    
    int num_bounced_surfaces = 0;
    float value = 0.0;
    float oldvalue = 0.0;
    bool includeForces = false;
    bool includeEnergy = true;
    
    if (endSimulation == false) {
        // first check source milestone crossings
        for (int i=0; i<integrator.getNumSrcMilestoneGroups(); i++) {
            value = context.calcForcesAndEnergy(includeForces, includeEnergy, srcbitvector[i]);
            if (srcMilestoneValues[i] == -INFINITY) {
                // First timestep
                srcMilestoneValues[i] = value;
            }
            oldvalue = srcMilestoneValues[i];
            if ((value - oldvalue) != 0.0) {
                // The source milestone has been crossed
                if (endOnSrcMilestone == true) {
                    endSimulation = true;
                    num_bounced_surfaces++;
                    ofstream datafile;
                    datafile.open(outputFileName, std::ios_base::app);
                    datafile << integrator.getSrcMilestoneGroup(i) << "," << crossingCounter << "," << context.getTime() << "\n";
                    datafile.close();
                } else {
                    crossedSrcMilestone = true;
                    context.setTime(0.0); // reset the timer
                    srcMilestoneValues[i] = value;
                }
            } 
        }
        
        // then check destination milestone crossings
        for (int i=0; i<integrator.getNumDestMilestoneGroups(); i++) {
            value = context.calcForcesAndEnergy(includeForces, includeEnergy, destbitvector[i]);
            if (destMilestoneValues[i] == -INFINITY) {
                // First timestep
                destMilestoneValues[i] = value;
            }
            oldvalue = destMilestoneValues[i];
            if ((value - oldvalue) != 0.0) {
                // The destination milestone has been crossed
                endSimulation = true;
                num_bounced_surfaces++;
                ofstream datafile;
                datafile.open(outputFileName, std::ios_base::app);
                if ((crossedSrcMilestone == true) || (endOnSrcMilestone == true)) {
                    datafile << integrator.getDestMilestoneGroup(i) << "," << crossingCounter << "," << context.getTime() << "\n";
                } else {
                    datafile << integrator.getDestMilestoneGroup(i) << "*," << crossingCounter << "," << context.getTime() << "\n";
                }
                datafile.close();
            } 
        }
        
        if (endSimulation == true) {
            // Then a crossing event has just occurred.
            if (saveStateBool == true && num_bounced_surfaces == 1) {
                State myState = context.getOwner().getState(State::Positions | State::Velocities);
                stringstream buffer;
                stringstream number_str;
                number_str << "_" << crossingCounter << "_" << crossingCounter;
                string trueFileName = saveStateFileName + number_str.str();
                XmlSerializer::serialize<State>(&myState, "State", buffer);
                ofstream statefile; // open datafile for writing
                statefile.open(trueFileName, std::ios_base::trunc); // append to file
                statefile << buffer.rdbuf();
                statefile.close(); // close data file
            }
            crossingCounter ++;
        }
    }
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateElberLangevinMiddleStepKernel::computeKineticEnergy(ContextImpl& context, const ElberLangevinMiddleIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}
