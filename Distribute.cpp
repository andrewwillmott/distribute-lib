//
// Distribute.cpp
//
// Distribute uniform u32 in various ways
//
// Andrew Willmott
//

#include "Distribute.hpp"

using namespace DL;

namespace
{
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

uint32_t DL::ModWeighted(uint32_t u, int numWeights, const float weights[])
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

uint32_t DL::ModWeighted(uint32_t u, int numWeights, const int weights[])
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

int32_t DL::ToInt32Weighted(uint32_t u, int32_t numWeights, const float weights[])
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

float DL::ToFloatGaussian(uint32_t u)
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

Vec2f DL::ToFloatGaussian(uint32_t u0, uint32_t u1)
{
    // Full Box-Muller transform
    Vec2f v2 = ToCircle(u0);

    float r = ToFloat(u1);
    float f = sqrtf(-2.0f * logf(r));

    return f * v2;
}


// 2D

Vec2f DL::ToDir2(uint32_t u)
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

Vec2f DL::ToDir2(uint32_t u, float s)
{
    float theta = ToFloatSigned(u, vlf_pi * s);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);
    return Vec2f(-ss, cs);
}

Vec2f DL::ToDisc(uint32_t u)
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

Vec2f DL::ToCircle(uint32_t u0)
{
    float theta = ToFloatSigned(u0, vlf_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    return Vec2f(cs, ss);
}

Vec2f DL::ToDisc(uint32_t u0, uint32_t u1)
{
    float theta = ToFloatSigned(u0, vlf_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float r2 = ToFloat(u1);

    return Vec2f(cs, ss) * sqrtf(r2);
}

Vec2f DL::ToRing(uint32_t u0, uint32_t u1, float s)
{
    float phi = ToFloatSigned(u0, vlf_pi);

    float r2 = (1.0f - s) * (1.0f - s);
    r2 += ToFloat(u1) * (1.0f - r2);

    float y, x;
    sincosf_local(phi, &y, &x);

    return sqrtf(r2) * Vec2f(x, y);
}


// 3D

Vec3f DL::ToDir3(uint32_t u)
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

Vec3f DL::ToBall(uint32_t u)
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

Vec3f DL::ToEllipsoid(uint32_t u, Vec3f min, Vec3f max)
{
    return (min + max + ToBall(u) * (max - min)) * 0.5f;
}

Vec3f DL::ToTorus(uint32_t u, float r)
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
        float len2xy = sqr(p.x) + sqr(p.y);

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

        return p;
    }
}

Vec3f DL::ToDir3(uint32_t u0, uint32_t u1)
{
    float theta = ToFloatSigned(u0, vlf_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float cz = ToFloatSigned(u1);
    float sz = sqrtf(1.0f - cz * cz);

    return Vec3f(cs * sz, ss * sz, cz);
}

Vec3f DL::ToDir3(uint32_t u0, uint32_t u1, float s)
{
    float theta = ToFloatSigned(u0, vlf_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float cz = ToFloat(u1, 1.0f - 2.0f * s, 1.0f);
    float sz = sqrtf(1.0f - cz * cz);

    return Vec3f(cs * sz, ss * sz, cz);
}

Vec3f DL::ToBall(uint32_t u0, uint32_t u1, uint32_t u2)
{
    float theta = ToFloatSigned(u0, vlf_pi);

    float ss, cs;
    sincosf_local(theta, &ss, &cs);

    float cz = ToFloatSigned(u1);
    float sz = sqrtf(1.0f - cz * cz);

    float r3 = ToFloat(u2);

    return Vec3f(cs * sz, ss * sz, cz) * cbrtf(r3);
}

Vec3f DL::ToEllipsoid(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max)
{
    return (min + max + ToBall(u0, u1, u2) * (max - min)) * 0.5f;
}

Vec3f DL::ToTorus(uint32_t u0, uint32_t u1, uint32_t u2, float r)
{
    r *= 0.5f;

    float theta = ToFloatSigned(u0, vlf_pi);

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

Mat3f DL::ToRot3(uint32_t u0, uint32_t u1, uint32_t u2)
{
    // https://doc.lagout.org/Others/Game%20Development/Programming/Graphics%20Gems%203.pdf
    float theta = ToFloat(u0, vlf_twoPi);  // Rotation about the pole (Z).
    float phi   = ToFloat(u1, vlf_twoPi);  // For direction of pole deflection.
    float z     = ToFloat(u2, 2.0f);       // For magnitude of pole deflection.

    float r  = sqrtf(z);
    float vx = sinf(phi) * r;
    float vy = cosf(phi) * r;
    float vz = sqrtf(2.0f - z);

    float st = sinf(theta);
    float ct = cosf(theta);
    float sx = vx * ct - vy * st;
    float sy = vx * st + vy * ct;

    Mat3f m
    (
        { vx * sx - ct, vx * sy - st, vx * vz  },
        { vy * sx + st, vy * sy - ct, vy * vz  },
        { vz * sx     , vz * sy     , 1.0f - z }
    );

    DL_ASSERT(IsOrthoNormal(m, 1e-5f));
    return m;
}

Vec4f DL::ToQuat(uint32_t u0, uint32_t u1, uint32_t u2, uint32_t u3)
{
    float th0 = ToFloatSigned(u0, vlf_pi);
    float th1 = ToFloatSigned(u2, vlf_pi);

    float d0 = ToFloat(u1);
    float d1 = ToFloat(u3);

    float ss, cs;
    sincosf_local(th0, &ss, &cs);

    Vec2f p0(cs, ss);
    p0 *= sqrtf(d0);

    sincosf_local(th1, &ss, &cs);
    Vec2f p1(cs, ss);

    float r1 = sqrtf(d1);
    p1 *= r1;

    float s0 = sqrtf(1.0f - d0);
    float s1 = 1.0f / (r1 + 1e-12f);
    float s  = s0 * s1;

    DL_ASSERT(!IsNAN(s));
    DL_ASSERT(!HasNAN(p0));
    DL_ASSERT(!HasNAN(p1));

    Vec4f q(p0.y, s * p1.x, s * p1.y, p0.x);
    return q;
}

#define BRANCHLESS

#ifdef BRANCHLESS

Vec2f DL::ToTriangleBakuOwen(uint32_t u)
{
    float cx = 0.0f, cy = 0.0f;
    float w = 0.5f;

    for (int i = 0; i < 16; i++)
    {
        bool flip = (u & 3) == 0;

        cy += ((u & 1) == 0) * w;
        cx += ((u & 2) == 0) * w;

        w *= flip ? -0.5f : 0.5f;

        u >>= 2;
    }

    return Vec2f(cx + w / 3.0f, cy + w / 3.0f);
}

Vec2f DL::ToTriangleBakuOwenRev(uint32_t u)
{
    float cx = 0.0f, cy = 0.0f;
    float w = 0.5f;

    for (int i = 0; i < 16; i++)
    {
        uint32_t d = (u >> 30);

        cy += ((d & 1) == 0) * w;
        cx += ((d & 2) == 0) * w;

        w *= (d == 0) ? -0.5f : 0.5f;

        u <<= 2;
    }

    return Vec2f(cx + w / 3.0f, cy + w / 3.0f);
}

#else

Vec2f DL::ToTriangleBakuOwen(uint32_t u)
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

Vec2f DL::ToTriangleBakuOwenRev(uint32_t u)
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

#endif


