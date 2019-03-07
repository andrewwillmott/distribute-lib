//
//  File:       Distribute.h
//
//  Function:   Distribute uniform u32 values in various ways
//
//  Copyright:  Andrew Willmott, 2018
//

#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

#include <math.h>
#include <stdint.h>

#ifndef DL_ASSERT
    #define DL_ASSERT(X)
#endif

#ifndef DL_VEC2F_CONVERT
    #define DL_VEC2F_CONVERT    // E.g, #define DL_VEC2F_CONVERT Vec2f(MyV2 v2) : x(v2.x), y(v2.y) {}; operator MyV2() const { return { x, y }; }
    #define DL_VEC3F_CONVERT
#endif

namespace DistLib
{
    struct Vec2f { float x; float y;                    Vec2f() {}; Vec2f(float xi, float yi)                     : x(xi), y(yi)               {}; DL_VEC2F_CONVERT };
    struct Vec3f { float x; float y; float z;           Vec3f() {}; Vec3f(float xi, float yi, float zi)           : x(xi), y(yi), z(zi)        {}; DL_VEC3F_CONVERT };


    // --------------------------------------------------------------------------
    // Helpers for remapping an evenly distributed uint32_t value to the desired
    // distribution shape and range.
    //
    // These functions can be used with any RNG that returns a uint32_t or
    // implements operator uint32_t(), e.g.,
    //
    //   int i = ToInt32(my_rand(), 10);
    //
    //   cRNG rng;
    //   float r = ToFloat(rng, 1.0f, 20.0f);
    //
    // --------------------------------------------------------------------------

    // 1D
    uint32_t ToUInt32         (uint32_t u);                         ///< Returns x (here for consistency)
    uint32_t ToUInt32         (uint32_t u, uint32_t limit);         ///< Returns [0, limit - 1]. Must have limit >= 0.
    uint32_t ToUInt32Inclusive(uint32_t u, uint32_t limit);         ///< Returns [0, limit]. Must have limit >= 0.
    uint32_t ToUInt32         (uint32_t u, uint32_t a, uint32_t b); ///< Returns [a, b - 1]. Must have a <= b.
    uint32_t ToUInt32Inclusive(uint32_t u, uint32_t a, uint32_t b); ///< Returns [a, b]. Must have a <= b.

    int32_t  ToInt32          (uint32_t u);                         ///< Returns a positive int32_t
    int32_t  ToInt32Signed    (uint32_t u);                         ///< Returns an int32_t (can be -ve)
    int32_t  ToInt32          (uint32_t u, int32_t limit);          ///< Returns [0, limit - 1]. Must have limit >= 0 -- asserts on this.
    int32_t  ToInt32Inclusive (uint32_t u, int32_t limit);          ///< Returns [0, limit]. Must have limit >= 0.
    int32_t  ToInt32Signed    (uint32_t u, int32_t limit);          ///< Returns [-limit, limit]. Must have limit >= 0.
    int32_t  ToInt32          (uint32_t u, int32_t a, int32_t b);   ///< Returns [a, b - 1]. Must have a <= b.
    int32_t  ToInt32Inclusive (uint32_t u, int32_t a, int32_t b);   ///< Returns [a, b]. Must have a <= b.

    float    ToFloat      (uint32_t u);                             ///< Returns a float in the range [ 0, 1]
    float    ToFloatSigned(uint32_t u);                             ///< Returns a float in the range [-1, 1]
    float    ToFloat      (uint32_t u, float a);                    ///< Returns a float in the range [ 0, a]
    float    ToFloatSigned(uint32_t u, float a);                    ///< Returns a float in the range [-a, a]
    float    ToFloat      (uint32_t u, float a, float b);           ///< Returns a float in the range [ a, b]

    // 2D
    Vec2f ToSquare    (uint32_t u);                                 ///< Returns random point in the unit square
    Vec2f ToRectangle (uint32_t u, Vec2f min, Vec2f max);           ///< Returns random point between min and max
    Vec2f ToTriangle  (uint32_t u);                                 ///< Returns random point in triangle with vertices (1, 0), (0, 0), (1, 0).
    Vec2f ToTriangle  (uint32_t u, Vec2f v0, Vec2f v1, Vec2f v2);   ///< Returns random point in triangle with vertices v0, v1, v2.
    Vec2f ToDirection2(uint32_t u);                                 ///< Returns random direction vector
    Vec2f ToCircle    (uint32_t u);                                 ///< Returns random point in the unit circle
    Vec2f ToRing      (uint32_t u, float r);                        ///< Returns random point in the unit ring with given width. r=1 -> full circle.
    
    // 3D
    Vec3f ToCube      (uint32_t u);                                 ///< Returns random point in the unit cube
    Vec3f ToCuboid    (uint32_t u, Vec3f min, Vec3f max);           ///< Returns random point between min and max
    Vec3f ToTriangle  (uint32_t u, Vec3f v0, Vec3f v1, Vec3f v2);   ///< Returns random point in triangle with vertices v0, v1, v2.
    Vec3f ToDirection3(uint32_t u);                                 ///< Returns random direction vector
    Vec3f ToSphere    (uint32_t u);                                 ///< Returns random point in the unit sphere
    Vec3f ToEllipsoid (uint32_t u, Vec3f min, Vec3f max);           ///< Returns random point in the ellipsoid defined by min and max.
    Vec3f ToTorus     (uint32_t u, float r);                        ///< Returns random point from a torus of the given radial width. r=0 -> circle, r=1 -> no hole in the middle.

    // 2D Full-u versions. These are costlier and less convenient (requiring multiple rng inputs), but higher quality, and ensure the character of the generator is reflected across all dimensions.
    Vec2f ToSquare    (uint32_t u0, uint32_t u1);                               ///< Returns random point in the unit square
    Vec2f ToRectangle (uint32_t u0, uint32_t u1, Vec2f min, Vec2f max);         ///< Returns random point between min and max
    Vec2f ToTriangle  (uint32_t u0, uint32_t u1);                               ///< Returns random point in triangle with vertices (1, 0), (0, 0), (1, 0).
    Vec2f ToTriangle  (uint32_t u0, uint32_t u1, Vec2f v0, Vec2f v1, Vec2f v2); ///< Returns random point in triangle with vertices v0, v1, v2.
    Vec2f ToDirection (uint32_t u0);                                            ///< Returns random direction vector.
    Vec2f ToCircle    (uint32_t u0, uint32_t u1);                               ///< Returns random point in the unit circle. u0 affects angular distribute, u1 the radial. E.g., ToCircle(r1, ModTriangle(r2)) gives more samples in the centre.
    Vec2f ToRing      (uint32_t u0, uint32_t u1, float r);                      ///< Returns random point in the unit ring with given width. r=1 -> full circle.

    // 3D Full-u versions. These are costlier and less convenient (requiring multiple rng inputs), but higher quality, and ensure the character of the generator is reflected across all dimensions.
    Vec3f ToCube      (uint32_t u0, uint32_t u1, uint32_t u2);                          ///< Returns point vector in the unit cube
    Vec3f ToCuboid    (uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max);    ///< Returns random point between min and max
    Vec3f ToTriangle  (uint32_t u0, uint32_t u1, Vec3f v0, Vec3f v1, Vec3f v2);         ///< Returns random point in triangle with vertices v0, v1, v2.
    Vec3f ToDirection (uint32_t u0, uint32_t u1);                                       ///< Returns random direction vector
    Vec3f ToSphere    (uint32_t u0, uint32_t u1, uint32_t u2);                          ///< Returns random point in the unit sphere
    Vec3f ToEllipsoid (uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max);    ///< Returns random point in the ellipsoid defined by min and max.
    Vec3f ToTorus     (uint32_t u0, uint32_t u1, uint32_t u2, float r);                 ///< Returns random point from a torus with the given radius and bounds

    // Modifiers. These can be used to pre-warp the input distribution before calling the above.
    uint32_t ModTriangle (uint32_t u);      ///< Returns a triangle-shaped distribution. Useful in combination with other distributors, e.g., ToFloat(ModTriangle(u)), or ToSquare(ModTriangle(u0), ModTriangle(u1)).
    uint32_t ModGaussLike(uint32_t u);      ///< Returns a gaussian-like distribution with variance of 1/6. This is much cheaper than a proper gaussian. The tails are limited to the range, so ToFloat(ToGaussLike(u), a, b) will always return a number between a and b, which can be useful.

    uint32_t ModWeighted (uint32_t u, int numWeights, const int   weights[]);   ///< Returns distribution weighted to 'numWeights' bins, according to 'weights'. A weight can be zero (no samples), but at least one weight must be non-zero.
    uint32_t ModWeighted (uint32_t u, int numWeights, const float weights[]);   ///< Version of ModWeighted that takes float weights for convenience.

    uint32_t ModInvert   (uint32_t u);      ///< Inverts distribution
    uint32_t ModHalfUp   (uint32_t u);      ///< Expands the left side of a symmetric distribution, e.g., a triangle distribution is converted to a ramp up.
    uint32_t ModHalfDown (uint32_t u);      ///< Expands the right side of a symmetric distribution, e.g., a triangle distribution is converted to a ramp down.
    uint32_t ModSymmUp   (uint32_t u);      ///< Converts a ramp back to a symmetric distribution, inverse of ModHalfUp.
    uint32_t ModSymmDown (uint32_t u);      ///< Converts a ramp back to a symmetric distribution, inverse of ModHalfDown.

    // Specialty or convenience operations.
    float   ToFloatTriangle (uint32_t u, float a, float b); ///< Most useful variant of ModTriangle -- returns a triangle-shaped distribution in [a, b]. Shorthand for ToFloat(ModTriangle(u), a, b).
    float   ToFloatGaussLike(uint32_t u, float a, float b); ///< Most useful variant of ModGaussLike -- returns a gaussian-like distribution in [a, b].

    int32_t ToInt32Weighted(uint32_t u, int numWeights, const float weights[]); ///< Returns [0, numWeights - 1], distributed according to 'weights'. Cheaper version of ToInt32(ToWeighted(u, ...), numWeights)

    float   ToFloatGaussian(uint32_t u);                    ///< Returns a proper gaussian ("normal") distribution.
    float   ToFloatGaussian(uint32_t u, float mean, float std_dev);   ///< Returns a gaussian distribution with the given mean and standard deviation
    Vec2f   ToFloatGaussian(uint32_t u0, uint32_t u1);      ///< Full version of ToFloatGaussian via the Box-Muller transform, returns two normally-distributed samples for two inputs.

    Vec2f   ToTriangleHier   (uint32_t u);                  ///< Distributes u hierarchically in the triangle. (See Basu & Owen '14.) Can be more evenly spaced than default version but is more expensive.
    Vec2f   ToTriangleHierRev(uint32_t u);                  ///< Distributes u taking high bits first rather than low bits.
}


// --------------------------------------------------------------------------
// Inlines
// --------------------------------------------------------------------------

namespace DistLib
{
    // Internal helpers
    namespace Internal
    {
        constexpr float kFloatFromUInt32Scale = float(1.0 / 0xFFFFFFFF);
        constexpr float kFloatFromInt32Scale  = float(2.0 / 0xFFFFFFFF);

        inline uint32_t Next(uint32_t u)
        {
            return u * 1664525 + 1013904223;
        }
    }

    // uint32_t variants

    inline uint32_t ToUInt32(uint32_t x)
    {
        return x;
    }

    inline uint32_t ToUInt32(uint32_t u, uint32_t limit)
    {
        return (u * (uint64_t(limit))) >> 32;
    }

    inline uint32_t ToUInt32(uint32_t u, uint32_t a, uint32_t b)
    {
        DL_ASSERT(a <= b);
        return a + ToUInt32(u, b - a);
    }

    inline uint32_t ToUInt32Inclusive(uint32_t u, uint32_t limit)
    {
        return (u * (uint64_t(limit) + 1)) >> 32;
    }

    inline uint32_t ToUInt32Inclusive(uint32_t u, uint32_t a, uint32_t b)
    {
        DL_ASSERT(a <= b);
        return a + ToUInt32Inclusive(u, b - a);
    }

    // int32_t variants

    inline int32_t ToInt32(uint32_t u)
    {
        return int32_t(u & 0x7FFFFFFF);
    }

    inline int32_t ToInt32Signed(uint32_t u)
    {
        return int32_t(u ^ 0x80000000);
    }

    inline int32_t ToInt32Inclusive(uint32_t u, int32_t limit)
    {
        return int32_t((u * (uint64_t(limit) + 1)) >> 32);
    }

    inline int32_t ToInt32(uint32_t u, int32_t limit)
    {
        DL_ASSERT(limit >= 0);
        return int32_t((u * (uint64_t(limit))) >> 32);
    }

    inline int32_t ToInt32Signed(uint32_t u, int32_t limit)
    {
        DL_ASSERT(limit >= 0);
        return ((u * uint64_t(2 * limit + 1)) >> 32) - limit;
    }

    inline int32_t ToInt32(uint32_t u, int32_t a, int32_t b)
    {
        DL_ASSERT(a <= b);
        return a + ToUInt32(u, b - a);
    }

    inline int32_t ToInt32Inclusive(uint32_t u, int32_t a, int32_t b)
    {
        DL_ASSERT(a <= b);
        return a + ToUInt32Inclusive(u, b - a);
    }

    // float variants
    inline float ToFloat(uint32_t u)
    {
        return Internal::kFloatFromUInt32Scale * u;
    }

    inline float ToFloatSigned(uint32_t u)
    {
        return Internal::kFloatFromInt32Scale * u - float(1);
    }

    inline float ToFloat(uint32_t u, float a)
    {
        return a * (Internal::kFloatFromUInt32Scale * u);
    }

    inline float ToFloatSigned(uint32_t u, float a)
    {
        return a * (Internal::kFloatFromInt32Scale * u - float(1));
    }

    inline float ToFloat(uint32_t u, float a, float b)
    {
        return a + (b - a) * (Internal::kFloatFromUInt32Scale * u);
    }

    // Specialty

    inline float ToFloatTriangle(uint32_t u, float a, float b)
    {
        return ToFloat(ModTriangle(u), a, b);
    }

    inline float ToFloatGaussLike(uint32_t u, float a, float b)
    {
        return ToFloat(ModGaussLike(u), a, b);
    }

    inline float ToFloatGaussian(uint32_t u, float mean, float std_dev)
    {
        return ToFloatGaussian(u) * std_dev + mean;
    }

    // 2D

    inline float sqrlen(Vec2f v)             { return v.x * v.x + v.y * v.y; }
    inline Vec2f operator+(Vec2f a, Vec2f b) { return { a.x + b.x, a.y + b.y }; }
    inline Vec2f operator-(Vec2f a, Vec2f b) { return { a.x - b.x, a.y - b.y }; }
    inline Vec2f operator*(Vec2f a, Vec2f b) { return { a.x * b.x, a.y * b.y }; }
    inline Vec2f operator*(float s, Vec2f a) { return { s   * a.x, s   * a.y }; }
    inline Vec2f operator*(Vec2f a, float s) { return { s   * a.x, s   * a.y }; }

    inline Vec2f ToSquare(uint32_t u)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);

        return ToSquare(u0, u1);
    }

    inline Vec2f ToRectangle(uint32_t u, Vec2f min, Vec2f max)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);

        return ToRectangle(u0, u1, min, max);
    }

    inline Vec2f ToTriangle(uint32_t u)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);

        float x = ToFloat(u0);
        float y = ToFloat(u1);
        
        if (x + y > 1.0f)
            return { 1.0f - x, 1.0f - y };
        else
            return { x, y };
    }

    inline Vec2f ToTriangle(uint32_t u, Vec2f v0, Vec2f v1, Vec2f v2)
    {
        Vec2f c = ToTriangle(u);
        return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
    }

    // 3D

    inline float sqrlen(Vec3f v)             { return v.x * v.x + v.y * v.y + v.z * v.z; }
    inline Vec3f operator+(Vec3f a, Vec3f b) { return { a.x + b.x, a.y + b.y, a.z + b.z}; }
    inline Vec3f operator-(Vec3f a, Vec3f b) { return { a.x - b.x, a.y - b.y, a.z - b.z}; }
    inline Vec3f operator*(Vec3f a, Vec3f b) { return { a.x * b.x, a.y * b.y, a.z * b.z}; }
    inline Vec3f operator*(float s, Vec3f a) { return { s   * a.x, s   * a.y, s   * a.z}; }
    inline Vec3f operator*(Vec3f a, float s) { return { s   * a.x, s   * a.y, s   * a.z}; }

    inline Vec3f ToCube(uint32_t u)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);
        uint32_t u2 = Internal::Next(u1);

        return ToCube(u0, u1, u2);
    }

    inline Vec3f ToCuboid(uint32_t u, Vec3f min, Vec3f max)
    {
        uint32_t u0 = u;
        uint32_t u1 = Internal::Next(u0);
        uint32_t u2 = Internal::Next(u1);

        return ToCuboid(u0, u1, u2, min, max);
    }

    inline Vec3f ToTriangle(uint32_t u, Vec3f v0, Vec3f v1, Vec3f v2)
    {
        Vec2f c = ToTriangle(u);
        return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
    }

    // 2D Full

    inline Vec2f ToSquare(uint32_t u0, uint32_t u1)
    {
        return Vec2f(ToFloatSigned(u0), ToFloatSigned(u1));
    }

    inline Vec2f ToRectangle(uint32_t u0, uint32_t u1, Vec2f min, Vec2f max)
    {
        return Vec2f(ToFloat(u0, min.x, max.x), ToFloat(u1, min.y, max.y));
    }

    inline Vec2f ToTriangle(uint32_t u0, uint32_t u1)
    {
        float x = ToFloat(u0);
        float t = sqrtf(ToFloat(u1));

        return { t * x, 1.0f - t};
    }

    inline Vec2f ToTriangle(uint32_t u0, uint32_t u1, Vec2f v0, Vec2f v1, Vec2f v2)
    {
        Vec2f c = ToTriangle(u0, u1);
        return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
    }

    // 3D Full

    inline Vec3f ToCube(uint32_t u0, uint32_t u1, uint32_t u2)
    {
        return Vec3f(ToFloatSigned(u0), ToFloatSigned(u1), ToFloatSigned(u2));
    }

    inline Vec3f ToCuboid(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max)
    {
        return Vec3f(ToFloat(u0, min.x, max.x), ToFloat(u1, min.y, max.y), ToFloat(u2, min.z, max.z));
    }

    inline Vec3f ToTriangle(uint32_t u0, uint32_t u1, Vec3f v0, Vec3f v1, Vec3f v2)
    {
        Vec2f c = ToTriangle(u0, u1);
        return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
    }

    // Modifiers
    inline uint32_t ModTriangle(uint32_t u)
    {
        uint32_t x0 = u;
        uint32_t x1 = Internal::Next(x0);

        // Variance is 1/6
        uint32_t carry = (x0 & x1 & 1);

        return ((x0 >> 1) + (x1 >> 1) + carry);
    }

    inline uint32_t ModInvert(uint32_t u)
    {
        return ~u;
    }

    inline uint32_t ModHalfDown(uint32_t u)
    {
        return (u & 0x80000000) ? (u << 1) : ~(u << 1);
    }

    inline uint32_t ModHalfUp(uint32_t u)
    {
        return (u & 0x80000000) ? ~(u << 1) : (u << 1);
    }

    inline uint32_t ModSymmUp(uint32_t u)
    {
        return (u & 1) ? ~(u >> 1) : (u >> 1);
    }

    inline uint32_t ModSymmDown(uint32_t u)
    {
        return (u & 1) ? (u << 1) : ~(u << 1);
    }
}

#endif

