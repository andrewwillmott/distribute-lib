//
//  File:       Distribute.h
//
//  Function:   Distribute uniform u32 in various ways
//
//  Copyright:  Andrew Willmott, 2018
//

#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

#include "VLVec234f.h"

namespace DistLib
{
    using ::Vec2f;
    using ::Vec3f;

    // --------------------------------------------------------------------------
    // Helpers for remapping a uint32_t value to the desired distribution type and range.
    //
    // These can be conveniently used with any RNG that implements operator uint32_t(), e.g.,
    // cRNG rng;
    // float r = ToFloat(rng, 1.0f, 20.0f);
    // --------------------------------------------------------------------------

    uint32_t ToUInt32         (uint32_t u);                         ///< Returns x (here for consistency)
    uint32_t ToUInt32         (uint32_t u, uint32_t limit);         ///< Returns [0, limit - 1]. Must have limit >= 0.
    uint32_t ToUInt32         (uint32_t u, uint32_t a, uint32_t b); ///< Returns [a, b - 1]. Must have a <= b.
    uint32_t ToUInt32Inclusive(uint32_t u, uint32_t limit);         ///< Returns [0, limit]. Must have limit >= 0.
    uint32_t ToUInt32Inclusive   (uint32_t u, uint32_t a, uint32_t b); ///< Returns [a, b]. Must have a <= b.

    int32_t  ToInt32          (uint32_t u);                         ///< Returns a positive int32_t
    int32_t  ToInt32Signed    (uint32_t u);                         ///< Returns an int32_t (can be -ve)
    int32_t  ToInt32          (uint32_t u, int32_t limit);          ///< Returns [0, limit - 1]. Must have limit >= 0 -- asserts on this.
    int32_t  ToInt32Inclusive (uint32_t u, int32_t limit);          ///< Returns [0, limit]. Must have limit >= 0.
    int32_t  ToInt32Signed    (uint32_t u, int32_t limit);          ///< Returns [-limit, limit]. Must have limit >= 0.
    int32_t  ToInt32          (uint32_t u, int32_t a, int32_t b);   ///< Returns [a, b - 1]. Must have a <= b.
    int32_t  ToInt32Inclusive (uint32_t u, int32_t a, int32_t b);   ///< Returns [a, b]. Must have a <= b.

    // float variants
    float ToFloat      (uint32_t u);                    ///< Returns a float in the range [ 0, 1]
    float ToFloatSigned(uint32_t u);                    ///< Returns a float in the range [-1, 1]
    float ToFloat      (uint32_t u, float a);           ///< Returns a float in the range [ 0, a]
    float ToFloatSigned(uint32_t u, float a);           ///< Returns a float in the range [-a, a]
    float ToFloat      (uint32_t u, float a, float b);  ///< Returns a float in the range [ a, b]

    // Warp
    uint32_t ToTriangle(uint32_t u);
    ///< Returns a triangle-shaped distribution with the max at zero. Most useful in combination with other distributors, e.g., ToFloat(ToUInt32Triangle(x))

    uint32_t ToGaussLike(uint32_t u);
    ///< Returns a gaussian-like distribution centred on zero.
    ///< This is much cheaper than a proper gaussian. The distribution tails are
    ///< limited to the range, so ToFloat(ToUInt32_gauss_like(...), a, b) will always
    ///< return a number between a and b, which can be useful.

    uint32_t ToWeighted(uint32_t u, int numWeights, const int   weights[]);    ///< Returns x (here for consistency)
    uint32_t ToWeighted(uint32_t u, int numWeights, const float weights[]);    ///< Returns x (here for consistency)

    uint32_t ToHalfUp  (uint32_t u);    ///< Expands the left side of a symmetric distribution, e.g., a triangle distribution is converted to a ramp up.
    uint32_t ToHalfDown(uint32_t u);    ///< Expands the right side of a symmetric distribution, e.g., a triangle distribution is converted to a ramp down.

    // Specialty
    float ToFloatTriangle(uint32_t u, float a, float b);  ///< Returns a triangle-shaped distribution in [a, b]. Shorthand for ToFloat(ToUInt32Triangle(u), a, b).
    float ToFloatGaussLike(uint32_t u, float a, float b); ///< Returns a gaussian-like distribution in [a, b].

    float ToFloatGaussian(uint32_t u);
    ///< Returns a proper gaussian distribution given two inputs. Somewhat expensive, requires a log() and sqrt().
    float ToFloatGaussian(uint32_t u, float mean, float std_dev);
    ///< Returns a gaussian distribution with the given parameters

    int32_t ToInt32Weighted(uint32_t u, int numWeights, const float weights[]); ///< Returns [0, numWeights - 1], distributed according to 'weights'. Cheaper version of ToInt32(ToWeighted()) 

    float ToFloatGaussian(uint32_t u0, uint32_t u1);
    ///< Returns a gaussian distribution with the given parameters

    // 2D
    Vec2f ToVec2f         (uint32_t u);                         ///< Returns random vector in the unit square
    Vec2f ToVec2f         (uint32_t u, Vec2f min, Vec2f max);   ///< Return random vector in the given AABB
    Vec2f ToVec2fDirection(uint32_t u);                         ///< Returns random direction vector
    Vec2f ToVec2fCircle   (uint32_t u);                         ///< Returns random vector in the unit circle
    Vec2f ToVec2fRing     (uint32_t u, float r);                ///< Returns random vector in the unit ring with given width. r=1 -> full circle.

    // Full-u versions. These are costlier but higher quality, and work when XXX
    Vec2f ToVec2f         (uint32_t u0, uint32_t u1);               ///< Returns random vector in the unit square
    Vec2f ToVec2f         (uint32_t u0, uint32_t u1, Vec2f min, Vec2f max);   ///< Return random vector in the given AABB
    Vec2f ToVec2fDirection(uint32_t u0, uint32_t u1);               ///< Returns random direction vector
    Vec2f ToVec2fCircle(uint32_t u0, uint32_t u1);                  //!< u0 affects angular distribute, u1 the radial. E.g., ToVec2fCircle(r1, ToTriangle(r2)) gives more samples in the centre.
    Vec2f ToVec2fRing     (uint32_t u0, uint32_t u1, float r);      ///< Returns random vector in the unit ring with given width. r=1 -> full circle.

    // 3D
    Vec3f ToVec3f         (uint32_t u);                         ///< Returns random vector in the unit cube
    Vec3f ToVec3f         (uint32_t u, Vec3f min, Vec3f max);   ///< Return random vector in the given AABB
    Vec3f ToVec3fDirection(uint32_t u);                         ///< Returns random direction vector
    Vec3f ToVec3fSphere   (uint32_t u);                         ///< Returns random vector in the unit sphere
    Vec3f ToVec3fEllipsoid(uint32_t u, Vec3f min, Vec3f max);   ///< Generalised version of Tosphere, handles any axis-aligned ellipsoid with the given bounds
    Vec3f ToVec3fTorus    (uint32_t u, Vec3f min, Vec3f max, float r);  ///< Returns random point from a torus with the given radius and bounds

    // full-u versions
    Vec3f ToVec3f         (uint32_t u0, uint32_t u1, uint32_t u2);  ///< Returns random vector in the unit cube
    Vec3f ToVec3f         (uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max);    ///< Return random vector in the given range
    Vec3f ToVec3fDirection(uint32_t u0, uint32_t u1, uint32_t u2);  ///< Returns random direction vector
    Vec3f ToVec3fSphere   (uint32_t u0, uint32_t u1, uint32_t u2);  ///< Returns random vector in the unit sphere
    Vec3f ToVec3fEllipsoid(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max);    ///< Returns random vector in the unit sphere
    Vec3f ToVec3fTorus    (uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max, float r);  ///< Returns random point from a torus with the given radius and bounds
}


// --------------------------------------------------------------------------
// Inlines
// --------------------------------------------------------------------------

// uint32_t variants

inline uint32_t DistLib::ToUInt32(uint32_t x)
{
    return x;
}

inline uint32_t DistLib::ToUInt32(uint32_t u, uint32_t limit)
{
    return (u * (uint64_t(limit))) >> 32;
}

inline uint32_t DistLib::ToUInt32(uint32_t u, uint32_t a, uint32_t b)
{
    CL_ASSERT(a <= b);
    return a + ToUInt32(u, b - a);
}

inline uint32_t DistLib::ToUInt32Inclusive(uint32_t u, uint32_t limit)
{
    return (u * (uint64_t(limit) + 1)) >> 32;
}

inline uint32_t DistLib::ToUInt32Inclusive(uint32_t u, uint32_t a, uint32_t b)
{
    CL_ASSERT(a <= b);
    return a + DistLib::ToUInt32Inclusive(u, b - a);
}

// int32_t variants

inline int32_t DistLib::ToInt32(uint32_t u)
{
    return int32_t(u & 0x7FFFFFFF);
}

inline int32_t DistLib::ToInt32Signed(uint32_t u)
{
    return int32_t(u ^ 0x80000000);
}

inline int32_t DistLib::ToInt32Inclusive(uint32_t u, int32_t limit)
{
    return int32_t((u * (uint64_t(limit) + 1)) >> 32);
}

inline int32_t DistLib::ToInt32(uint32_t u, int32_t limit)
{
    CL_ASSERT(limit >= 0);
    return int32_t((u * (uint64_t(limit))) >> 32);
}

inline int32_t DistLib::ToInt32Signed(uint32_t u, int32_t limit)
{
    CL_ASSERT(limit >= 0);
    return ((u * uint64_t(2 * limit + 1)) >> 32) - limit;
}

inline int32_t DistLib::ToInt32(uint32_t u, int32_t a, int32_t b)
{
    CL_ASSERT(a <= b);
    return a + ToUInt32(u, b - a);
}

inline int32_t DistLib::ToInt32Inclusive(uint32_t u, int32_t a, int32_t b)
{
    CL_ASSERT(a <= b);
    return a + ToUInt32Inclusive(u, b - a);
}

// float variants
namespace DistLib { namespace Internal
{
    constexpr float kFloatFromUInt32Scale = float(1.0 / 0xFFFFFFFF);
    constexpr float kFloatFromInt32Scale = float(2.0 / 0xFFFFFFFF);

    inline uint32_t Next(uint32_t u)
    {
        return u * 1664525 + 1013904223;
    }
}}

inline float DistLib::ToFloat(uint32_t u)
{
    return Internal::kFloatFromUInt32Scale * u;
}

inline float DistLib::ToFloatSigned(uint32_t u)
{
    return Internal::kFloatFromInt32Scale * u - float(1);
}

inline float DistLib::ToFloat(uint32_t u, float a)
{
    return a * (Internal::kFloatFromUInt32Scale * u);
}

inline float DistLib::ToFloatSigned(uint32_t u, float a)
{
    return a * (Internal::kFloatFromInt32Scale * u - float(1));
}

inline float DistLib::ToFloat(uint32_t u, float a, float b)
{
    return a + (b - a) * (Internal::kFloatFromUInt32Scale * u);
}

// Warp
inline uint32_t DistLib::ToTriangle(uint32_t u)
{
    uint32_t x0 = u;
    uint32_t x1 = Internal::Next(x0);

    // Variance is 1/6
    uint32_t carry = (x0 & x1 & 1);

    return ((x0 >> 1) + (x1 >> 1) + carry);
}

inline uint32_t DistLib::ToHalfDown(uint32_t u)
{
    return (u & 0x80000000) ? u << 1 : ~u << 1;
}
inline uint32_t DistLib::ToHalfUp(uint32_t u)
{
    return (u & 0x80000000) ? ~u << 1 : u << 1;
}

// Specialty
inline float DistLib::ToFloatTriangle(uint32_t u, float a, float b)
{
    return ToFloat(ToTriangle(u), a, b);
}

inline float DistLib::ToFloatGaussLike(uint32_t u, float a, float b)
{
    return ToFloat(ToGaussLike(u), a, b);
}

inline float DistLib::ToFloatGaussian(uint32_t u, float mean, float std_dev)
{
    return ToFloatGaussian(u) * std_dev + mean;
}

// 2D

inline Vec2f DistLib::ToVec2f(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);

    return ToVec2f(u0, u1);
}

inline Vec2f DistLib::ToVec2f(uint32_t u, Vec2f min, Vec2f max)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);

    return ToVec2f(u0, u1, min, max);
}

inline Vec2f DistLib::ToVec2f(uint32_t u0, uint32_t u1)
{
    return Vec2f(ToFloatSigned(u0), ToFloatSigned(u1));
}

inline Vec2f DistLib::ToVec2f(uint32_t u0, uint32_t u1, Vec2f min, Vec2f max)
{
    return Vec2f(ToFloat(u0, min.x, max.x), ToFloat(u1, min.y, max.y));
}

// Vec3f

inline Vec3f DistLib::ToVec3f(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return ToVec3f(u0, u1, u2);
}

inline Vec3f DistLib::ToVec3f(uint32_t u, Vec3f min, Vec3f max)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return ToVec3f(u0, u1, u2, min, max);
}

inline Vec3f DistLib::ToVec3f(uint32_t u0, uint32_t u1, uint32_t u2)
{
    return Vec3f(ToFloatSigned(u0), ToFloatSigned(u1), ToFloatSigned(u2));
}

inline Vec3f DistLib::ToVec3f(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max)
{
    return Vec3f(ToFloat(u0, min.x, max.x), ToFloat(u1, min.y, max.y), ToFloat(u2, min.z, max.z));
}

#endif
