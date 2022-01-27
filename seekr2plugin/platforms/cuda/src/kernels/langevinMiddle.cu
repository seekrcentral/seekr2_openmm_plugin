enum {VelScale, NoiseScale};

/**
 * Perform the first step of Langevin Middle integration: velocity step.
 */

extern "C" __global__ void integrateMmvtLangevinMiddlePart1(int numAtoms, 
            int paddedNumAtoms, mixed4* __restrict__ velm, 
            const long long* __restrict__ force, 
            const mixed2* __restrict__ dt) {
    mixed fscale = dt[0].y/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += fscale*velocity.w*force[index];
            velocity.y += fscale*velocity.w*force[index+paddedNumAtoms];
            velocity.z += fscale*velocity.w*force[index+paddedNumAtoms*2];
            velm[index] = velocity;
        }
    }
}

/**
 * Perform the second part of integration: position half step, then interact with heat bath,
 * then another position half step.
 */

extern "C" __global__ void integrateMmvtLangevinMiddlePart2(int numAtoms, 
        mixed4* __restrict__ velm, mixed4* __restrict__ posDelta,
        mixed4* __restrict__ oldDelta, const mixed* __restrict__ paramBuffer, 
        const mixed2* __restrict__ dt, mixed4* __restrict__ oldVelm,
        const float4* __restrict__ random, unsigned int randomIndex
        ) {
    mixed vscale = paramBuffer[VelScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed halfdt = 0.5f*dt[0].y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    randomIndex += index;
    while (index < numAtoms) {
        mixed4 velocity = velm[index];
        oldVelm[index] = velm[index];
        
        if (velocity.w != 0.0) {
            mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            mixed sqrtInvMass = SQRT(velocity.w);
            velocity.x = vscale*velocity.x + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
            //delta = delta + make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            delta.x += (mixed) halfdt*velocity.x;
            delta.y += (mixed) halfdt*velocity.y;
            delta.z += (mixed) halfdt*velocity.z;
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
        
        randomIndex += blockDim.x*gridDim.x;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Perform the third part of integration: apply constraint forces to velocities, then record
 * the constrained positions.
 */

extern "C" __global__ void integrateMmvtLangevinMiddlePart3(int numAtoms, 
         real4* __restrict__ posq, mixed4* __restrict__ velm,
         const mixed4* __restrict__ posDelta, mixed4* __restrict__ oldDelta, 
         const mixed2* __restrict__ dt, real4* __restrict__ oldPosq, 
         real4* __restrict__ posqCorrection) {
    mixed invDt = 1/dt[0].y;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        oldPosq[index] = posq[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            velocity.x += (delta.x-oldDelta[index].x)*invDt;
            velocity.y += (delta.y-oldDelta[index].y)*invDt;
            velocity.z += (delta.z-oldDelta[index].z)*invDt;
            velm[index] = velocity;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
    }
}


extern "C" __global__ void integrateElberLangevinMiddlePart1(int numAtoms, 
            int paddedNumAtoms, mixed4* __restrict__ velm, 
            const long long* __restrict__ force, 
            const mixed2* __restrict__ dt) {
    mixed fscale = dt[0].y/(mixed) 0x100000000;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            velocity.x += fscale*velocity.w*force[index];
            velocity.y += fscale*velocity.w*force[index+paddedNumAtoms];
            velocity.z += fscale*velocity.w*force[index+paddedNumAtoms*2];
            velm[index] = velocity;
        }
    }
}

/**
 * Perform the second part of integration: position half step, then interact with heat bath,
 * then another position half step.
 */

extern "C" __global__ void integrateElberLangevinMiddlePart2(int numAtoms, 
        mixed4* __restrict__ velm, mixed4* __restrict__ posDelta,
        mixed4* __restrict__ oldDelta, const mixed* __restrict__ paramBuffer, 
        const mixed2* __restrict__ dt,
        const float4* __restrict__ random, unsigned int randomIndex
        ) {
    mixed vscale = paramBuffer[VelScale];
    mixed noisescale = paramBuffer[NoiseScale];
    mixed halfdt = 0.5f*dt[0].y;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    randomIndex += index;
    while (index < numAtoms) {
        mixed4 velocity = velm[index];
        
        if (velocity.w != 0.0) {
            mixed4 delta = make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            mixed sqrtInvMass = SQRT(velocity.w);
            velocity.x = vscale*velocity.x + noisescale*sqrtInvMass*random[randomIndex].x;
            velocity.y = vscale*velocity.y + noisescale*sqrtInvMass*random[randomIndex].y;
            velocity.z = vscale*velocity.z + noisescale*sqrtInvMass*random[randomIndex].z;
            velm[index] = velocity;
            //delta = delta + make_mixed4(halfdt*velocity.x, halfdt*velocity.y, halfdt*velocity.z, 0);
            delta.x += (mixed) halfdt*velocity.x;
            delta.y += (mixed) halfdt*velocity.y;
            delta.z += (mixed) halfdt*velocity.z;
            posDelta[index] = delta;
            oldDelta[index] = delta;
        }
        
        randomIndex += blockDim.x*gridDim.x;
        index += blockDim.x*gridDim.x;
    }
}

/**
 * Perform the third part of integration: apply constraint forces to velocities, then record
 * the constrained positions.
 */

extern "C" __global__ void integrateElberLangevinMiddlePart3(int numAtoms, 
         real4* __restrict__ posq, mixed4* __restrict__ velm,
         const mixed4* __restrict__ posDelta, mixed4* __restrict__ oldDelta, 
         const mixed2* __restrict__ dt,
         real4* __restrict__ posqCorrection) {
    mixed invDt = 1/dt[0].y;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
            mixed4 delta = posDelta[index];
            velocity.x += (delta.x-oldDelta[index].x)*invDt;
            velocity.y += (delta.y-oldDelta[index].y)*invDt;
            velocity.z += (delta.z-oldDelta[index].z)*invDt;
            velm[index] = velocity;
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
        }
    }
}


/**
 * Take a step back in time and reverse velocities
 */
extern "C" __global__ void mmvtBounce(int numAtoms, int paddedNumAtoms, 
            real4* __restrict__ posq,
            mixed4* __restrict__ velm, 
            const real4* __restrict__ oldPosq, 
            const mixed4* __restrict__ oldVelm) { 
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    while (index < numAtoms) {
        posq[index] = oldPosq[index];
        velm[index] = make_mixed4(-oldVelm[index].x, -oldVelm[index].y,
                -oldVelm[index].z, oldVelm[index].w);
        index += blockDim.x*gridDim.x;
    }
}