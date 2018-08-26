//
//  File:       CLDistribute.cpp
//
//  Function:   Distribute uniform u32 in various ways
//
//  Copyright:  Andrew Willmott, 2018
//

#include "Distribute.h"

using namespace DistLib;

namespace
{
    inline float InvSqrtFast(float x)
    {
        float xhalf = 0.5f * x;
        int32_t i = (int32_t&) x;

        i = 0x5f375a86 - (i >> 1);
        x = (float&) i;
        x = x * (1.5f - xhalf * x * x);

        return x;
    }
}

// 1D

uint32_t DistLib::ToWeighted(uint32_t u, int numWeights, const float weights[])
{
    float weightSum = 0.0f;

    for (int i = 0; i < numWeights; i++)
        weightSum += weights[i];

    float s = ToFloat(u, weightSum);

    for (int i = 0; i < numWeights; i++)
    {
        if (s < weights[i])
        {
            float fs = float(i) / numWeights + s / (numWeights * weights[i]);
            return fs * UINT32_MAX;
        }
        
        s -= weights[i];
    }

    return UINT32_MAX;
}

uint32_t DistLib::ToWeighted(uint32_t s, int numWeights, const int weights[])
{
    CL_ASSERT(numWeights > 0);

    uint32_t weightSum = 0;

    for (int i = 0; i < numWeights; i++)
        weightSum += weights[i];

    CL_ASSERT(weightSum >= 1);

    uint32_t step = UINT32_MAX / weightSum;
    
    for (int i = 0; i < numWeights; i++)
    {
        if (s < weights[i] * step)
        {
            uint32_t outStep = UINT32_MAX / numWeights;

            // avoiding 64-bit:
            // s < weights[i] * (UINT32_MAX / weightSum)
            // s * weightSum < weights[i] * UINT32_MAX;
            // so track max weights[i] and prescale? or use weightSum 
            // divisor then becomes inaccurate :(
            // most we can multiply 's' by is weightSum / weights[i]

            
            return outStep * i + weightSum * (s / weights[i]) / numWeights;

            // return outStep * i + (uint64_t(weightSum) * s) / (weights[i] * numWeights);
        }
        
        s -= weights[i] * step;
    }

    return UINT32_MAX;
}

int32_t DistLib::ToInt32Weighted(uint32_t u, int32_t numWeights, const float weights[])
{
    float weightSum = 0.0f;

    for (int i = 0; i < numWeights; i++)
        weightSum += weights[i];

    float s = ToFloat(u, weightSum);
    weightSum = 0.0f;

    for (int i = 0; i < numWeights; i++)
    {
        weightSum += weights[i];

        if (s < weightSum)
            return i;
    }

    return numWeights - 1;
}

uint32_t DistLib::ToGaussLike(uint32_t u)
{
    uint32_t x0 = u;
    uint32_t x1 = Internal::Next(x0);
    uint32_t x2 = Internal::Next(x1);
    uint32_t x3 = Internal::Next(x2);

#if 0
    uint32_t c0 = (x0 & x1 & 1);
    uint32_t c1 = (x2 & x3 & 1);

    uint32_t carry = c0 + c1 + (c0 & c1 & 1);
#else
    uint32_t carry = 1;  // because this situation is symmetrical we can use a constant at the expense of never getting 0 or UINT32_MAX
#endif

    uint32_t r = (x0 >> 2) + (x1 >> 2) + (x2 >> 2) + (x3 >> 2) + carry;

    return r;
}

float DistLib::ToFloatGaussian(uint32_t u)
{
    uint32_t ux = u;
    uint32_t uy = Internal::Next(ux);

    float x = ToFloatSigned(ux);
    float y = ToFloatSigned(uy);
    float s = x * x + y * y;

    while (s >= 1.0f || s == 0.0f)
    {
        // For convenenience we run a local LCG rather than going back to the source distribution.
        ux = Internal::Next(ux);
        uy = Internal::Next(uy);

        x = ToFloatSigned(ux);
        y = ToFloatSigned(uy);
        s = x * x + y * y;
    }

    float f = sqrtf(-logf(s) / s);

    return f * (x + y);
}


// 2D

DistLib::Vec2f DistLib::ToVec2fDirection(uint32_t u)
{
    // Reject points outside unit sphere and inside epsilon ball. May seem clunky but expected iteration count is < 2
    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);

        Vec2f p = Vec2f(ToFloatSigned(u0), ToFloatSigned(u1));
        float p2 = sqrlen(p);

        if (p2 <= 1.0f && p2 >= 1e-4f)
            return p * InvSqrtFast(p2);

        u = Internal::Next(u1);
    }
}

DistLib::Vec2f DistLib::ToVec2fCircle(uint32_t u)
{
    // Reject points outside unit sphere. May seem clunky but expected iteration count is < 2
    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);

        Vec2f p = Vec2f(ToFloatSigned(u0), ToFloatSigned(u1));

        if (sqrlen(p) <= 1.0f)
            return p;

        u = Internal::Next(u1);
    }
}

DistLib::Vec2f DistLib::ToVec2fRing(uint32_t u, float r)
{
    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);

        Vec2f v = Vec2f(ToFloatSigned(u0), ToFloatSigned(u1));

        float Vec2f = sqrlen(v);

        if (Vec2f >= 1.0f || Vec2f < 1e-4f)
        {
            u = Internal::Next(u1);
            continue;
        }

        float r1 = sqrtf(Vec2f);

        return v * (r1 * r + 1.0f - r) / r1;
    }
}

DistLib::Vec2f DistLib::ToVec2fCircle(uint32_t u0, uint32_t u1)
{
    float theta = ToFloatSigned(u0, vl_pi);
    float cs = cosf(theta);
    // float ss = sinf(theta); 
    float ss = copysignf(sqrtf(1.0f - cs * cs), theta);
    float r2 = ToFloat(u1);

    return Vec2f(cs, ss) * sqrtf(r2);
}

DistLib::Vec2f DistLib::ToVec2fRing(uint32_t u0, uint32_t u1, float r)
{
    float theta = ToFloatSigned(u0, vl_pi);
    float cs = cosf(theta);
    float ss = copysignf(sqrtf(1.0f - cs * cs), theta);
    float r0 = (1.0f - r) * (1.0f - r);
    float r2 = ToFloat(u1, r0, 1.0f);

    return Vec2f(cs, ss) * sqrtf(r2);
}


// 3D

DistLib::Vec3f DistLib::ToVec3fDirection(uint32_t u)
{
    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);
        uint32_t u2 = Internal::Next(u1);

        Vec3f p = Vec3f(ToFloatSigned(u0), ToFloatSigned(u1), ToFloatSigned(u2));
        float p2 = sqrlen(p);

        if (p2 <= 1.0f && p2 >= 1e-4f)
            return p * InvSqrtFast(p2);

        u = Internal::Next(u2);
    }
}

DistLib::Vec3f DistLib::ToVec3fSphere(uint32_t u)
{
    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);
        uint32_t u2 = Internal::Next(u1);

        Vec3f p = Vec3f(ToFloatSigned(u0), ToFloatSigned(u1), ToFloatSigned(u2));

        if (sqrlen(p) <= 1.0f)
            return p;

        u = Internal::Next(u2);
    }
}

DistLib::Vec3f DistLib::ToVec3fEllipsoid(uint32_t u, Vec3f min, Vec3f max)
{
    // TODO: rework
    Vec3f c = 0.5f * (min + max);
    Vec3f r = c - min;
    Vec3f r2 = r * r;

    Vec3f invR2 =
    {
        r2[0] > 0.0f ? 1.0f / r2[0] : 0.0f,
        r2[1] > 0.0f ? 1.0f / r2[1] : 0.0f,
        r2[2] > 0.0f ? 1.0f / r2[2] : 0.0f
    };

    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);
        uint32_t u2 = Internal::Next(u1);

        Vec3f p(ToFloat(u0, min.x, max.x), ToFloat(u1, min.y, max.y), ToFloat(u2, min.z, max.z));
        Vec3f v = (p - c);

        // is this point inside the ellipsoid?
        // pt[0]^2 / vec[0]^2 + ... <= 1
        if (dot(v * v, invR2) <= 1.0f)
            return p;

        u = Internal::Next(u2);
    }
}

DistLib::Vec3f DistLib::ToVec3fTorus(uint32_t u, Vec3f min, Vec3f max, float r)
{
    // Unit torus is r = 0.5, width = (0.5, 1)
    // The volume is 2pi * 0.5 * (pi * (0.5 * 1)) = pi^2 / 2 =~ 4.9. So the average
    // iteration count is a little less than for the sphere test.

    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);
        uint32_t u2 = Internal::Next(u1);

        Vec3f p = Vec3f(ToFloatSigned(u0), ToFloatSigned(u1), ToFloatSigned(u2));

        // is this point inside a r = 0.5 torus?
        float len2xy = sqrlen(p.AsVec2());

        if (len2xy < 1e-8f)
        {
            u = Internal::Next(u2);
            continue;
        }

        float len2z  = sqr(p.z);
        float r1 = sqrtf(len2xy);

        if (len2z + len2xy > r1)
        {
            u = Internal::Next(u2);
            continue;
        }

        // adjust according to true radius: r1' = r1 * r + (1 - r)
        float r1Dash = r + (1.0f - r) / r1;

        p *= Vec3f(r1Dash, r1Dash, r);

        return min + (max - min) * (0.5f * (p + vl_one));
    }
}

DistLib::Vec3f DistLib::ToVec3fSphere(uint32_t u0, uint32_t u1, uint32_t u2)
{
    float theta = ToFloatSigned(u0, vl_pi);
    float cs = cosf(theta);
    // float ss = sinf(theta); 
    float ss = copysignf(sqrtf(1.0f - cs * cs), theta);

    float cz = ToFloatSigned(u1);
    float sz = sqrtf(1.0f - cz * cz);

    float r3 = ToFloat(u2);
    
    return Vec3f(cs * sz, ss * sz, cz) * cbrt(r3);
}


DistLib::Vec3f DistLib::ToVec3fEllipsoid(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max)
{
    return min + ToVec3fSphere(u0, u1, u2) * (max - min);
}

DistLib::Vec3f DistLib::ToVec3fTorus(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max, float r)
{
    // can be seen as a disc of a certain r'?

    float cz = ToFloatSigned(u2, r);
    float sz = sqrtf(1.0f - cz * cz);


    float theta = ToFloatSigned(u0, vl_pi);
    float cs = cosf(theta);
    float ss = copysignf(sqrtf(1.0f - cs * cs), theta);
    float r0 = (1.0f - r) * (1.0f - r);
    float r2 = ToFloat(u1, r0, 1.0f);

    return Vec3f(cs * sqrtf(r2), ss * sqrtf(r2), cz * r);
}

