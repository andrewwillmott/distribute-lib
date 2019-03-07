//
//  File:       Distribute.cpp
//
//  Function:   Distribute uniform u32 in various ways
//
//  Copyright:  Andrew Willmott, 2018
//

#include "Distribute.h"

using namespace DistLib;

namespace
{
    constexpr float vl_pi = 3.14159265358979323846f;

    inline void sincosf_local(float theta, float* sOut, float* cOut)
    {
        float c = cosf(theta);
        
        *sOut = copysignf(sqrtf(1.0f - c * c), theta);
        *cOut = c;
    }

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

uint32_t DistLib::ModGaussLike(uint32_t u)
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

uint32_t DistLib::ModWeighted(uint32_t u, int numWeights, const float weights[])
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
            return uint32_t(fs * UINT32_MAX);
        }
        
        s -= weights[i];
    }

    return UINT32_MAX;
}

uint32_t DistLib::ModWeighted(uint32_t u, int numWeights, const int weights[])
{
    DL_ASSERT(numWeights > 0);

    uint32_t weightSum = 0;

    for (int i = 0; i < numWeights; i++)
        weightSum += weights[i];

    DL_ASSERT(weightSum >= 1);

    uint32_t step = UINT32_MAX / weightSum;
    
    for (int i = 0; i < numWeights; i++)
    {
        if (u < weights[i] * step)
        {
            uint32_t outStep = UINT32_MAX / numWeights;

            return outStep * i + (weightSum * (u / weights[i])) / numWeights; // re-arranged to stick to 32-bit arithmetic at the expense of a few bits
        }
        
        u -= weights[i] * step;
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

float DistLib::ToFloatGaussian(uint32_t u)
{
    float x, y, s;

    while (true)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);

        x = ToFloatSigned(u0);
        y = ToFloatSigned(u1);
        s = x * x + y * y;
        
        if (s <= 1.0f)
            break;

        u = Internal::Next(u1);
    }

    float f = sqrtf(-2.0f * logf(s) / s);

    return f * (x + y) * sqrtf(0.5f);
}

Vec2f DistLib::ToFloatGaussian(uint32_t u0, uint32_t u1)
{
    // Full Box-Muller transform
    Vec2f v2 = ToDirection(u0);

    float r = ToFloat(u1);
    float f = sqrtf(-2.0f * logf(r));

    return f * v2;
}


// 2D

Vec2f DistLib::ToDirection2(uint32_t u)
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

Vec2f DistLib::ToCircle(uint32_t u)
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

Vec2f DistLib::ToRing(uint32_t u, float r)
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

        return v * ((r1 * r + 1.0f - r) / r1);
    }
}

Vec2f DistLib::ToDirection(uint32_t u0)
{
    float theta = ToFloatSigned(u0, vl_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    return Vec2f(cs, ss);
}

Vec2f DistLib::ToCircle(uint32_t u0, uint32_t u1)
{
    float theta = ToFloatSigned(u0, vl_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float r2 = ToFloat(u1);

    return Vec2f(cs, ss) * sqrtf(r2);
}

Vec2f DistLib::ToRing(uint32_t u0, uint32_t u1, float r)
{
    float theta = ToFloatSigned(u0, vl_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float r0 = (1.0f - r) * (1.0f - r);
    float r2 = ToFloat(u1, r0, 1.0f);

    return Vec2f(cs, ss) * sqrtf(r2);
}


// 3D

Vec3f DistLib::ToDirection3(uint32_t u)
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

Vec3f DistLib::ToSphere(uint32_t u)
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

Vec3f DistLib::ToEllipsoid(uint32_t u, Vec3f min, Vec3f max)
{
    return (min + max + ToSphere(u) * (max - min)) * 0.5f;
}

Vec3f DistLib::ToTorus(uint32_t u, float r)
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
        float len2xy = p.x * p.x + p.y * p.y;

        if (len2xy < 1e-8f)
        {
            u = Internal::Next(u2);
            continue;
        }

        float len2z  = p.z * p.z;
        float r1 = sqrtf(len2xy);

        if (len2z + len2xy > r1)
        {
            u = Internal::Next(u2);
            continue;
        }

        // adjust according to true radius: r1' = r1 * r + (1 - r)
        float r1Dash = r + (1.0f - r) / r1;

        p = p * Vec3f(r1Dash, r1Dash, r);

        return p;
    }
}

Vec3f DistLib::ToDirection(uint32_t u0, uint32_t u1)
{
    float theta = ToFloatSigned(u0, vl_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float cz = ToFloatSigned(u1);
    float sz = sqrtf(1.0f - cz * cz);

    return Vec3f(cs * sz, ss * sz, cz);
}

Vec3f DistLib::ToSphere(uint32_t u0, uint32_t u1, uint32_t u2)
{
    float theta = ToFloatSigned(u0, vl_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float cz = ToFloatSigned(u1);
    float sz = sqrtf(1.0f - cz * cz);

    float r3 = ToFloat(u2);
    
    return Vec3f(cs * sz, ss * sz, cz) * cbrtf(r3);
}

Vec3f DistLib::ToEllipsoid(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max)
{
    return (min + max + ToSphere(u0, u1, u2) * (max - min)) * 0.5f;
}

Vec3f DistLib::ToTorus(uint32_t u0, uint32_t u1, uint32_t u2, float r)
{
    r *= 0.5f;
    
    float theta = ToFloatSigned(u0, vl_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float cz = ToFloatSigned(u2);
    float sz = sqrtf(1.0f - cz * cz);

    float r0 = (1.0f - r - sz * r);
    float r1 = (1.0f - r + sz * r);
    float r2 = ToFloat(u1, r0 * r0, r1 * r1);

    float rr = sqrtf(r2);
    
    return Vec3f(cs * rr, ss * rr, cz * r);
}

Vec2f DistLib::ToTriangleHier(uint32_t u)
{
    // Current tri defined by c, (cx + w, cy), (c.x, c.y + w)
    float cx = 0.0f, cy = 0.0f;
    float w = 1;

    for (int i = 0; i < 16; i++)
    {
        if (!u)
            break;  // picking middle for the remaining points

        w *= 0.5f;

        switch (u & 0x3)
        {
        case 0:     // middle (inverted)
            cx += w;
            cy += w;
            w = -w;
            break;
        case 1:     // bottom-right
            cx += w;
            break;
        case 2:     // top-left
            cy += w;
            break;
        case 3:     // bottom-left
            break;
        }

        u >>= 2;
    }

    return Vec2f(cx + w / 3.0f, cy + w / 3.0f);
}

Vec2f DistLib::ToTriangleHierRev(uint32_t u)
{
    // Current tri defined by c, (cx + w, cy), (c.x, c.y + w)
    float cx = 0.0f, cy = 0.0f;
    float w = 1;

    for (int i = 0; i < 16; i++)
    {
        if (!u)
            break;  // picking middle for the remaining points

        w *= 0.5f;

        switch ((u >> 30) & 0x3)
        {
        case 0:     // middle (inverted)
            cx += w;
            cy += w;
            w = -w;
            break;
        case 1:     // bottom-right
            cx += w;
            break;
        case 2:     // top-left
            cy += w;
            break;
        case 3:     // bottom-left
            break;
        }

        u <<= 2;
    }

    return Vec2f(cx + w / 3.0f, cy + w / 3.0f);
}
