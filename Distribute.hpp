//
// Distribute.hpp
//
// Distribute uniform u32 in various ways
//
// Andrew Willmott
//

#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

#include "VLMini.hpp"
#include <stdint.h>

#ifndef DL_ASSERT
    #define DL_ASSERT(X)
#endif

namespace DL
{
    // --------------------------------------------------------------------------
    // Helpers for remapping an evenly distributed uint32_t value to the desired
    // distribution shape and range.
    //
    // These functions can be used with any RNG that returns a uint32_t or
    // implements operator uint32_t(), e.g.,
    //
    //   int i = ToInt32(my_rand(), 10);
    //
    //   RNG rng;
    //   float r = ToFloat(rng, 1.0f, 20.0f);
    //
    // --------------------------------------------------------------------------

    // 1D
    uint32_t ToUInt32         (uint32_t u);                         // Returns x (here for consistency)
    uint32_t ToUInt32         (uint32_t u, uint32_t limit);         // Returns [0, limit - 1]. Must have limit >= 0.
    uint32_t ToUInt32Inclusive(uint32_t u, uint32_t limit);         // Returns [0, limit]. Must have limit >= 0.
    uint32_t ToUInt32         (uint32_t u, uint32_t a, uint32_t b); // Returns [a, b - 1]. Must have a <= b.
    uint32_t ToUInt32Inclusive(uint32_t u, uint32_t a, uint32_t b); // Returns [a, b]. Must have a <= b.

    int32_t  ToInt32          (uint32_t u);                         // Returns a positive int32_t
    int32_t  ToInt32Signed    (uint32_t u);                         // Returns an int32_t (can be -ve)
    int32_t  ToInt32          (uint32_t u, int32_t limit);          // Returns [0, limit - 1]. Must have limit >= 0 -- asserts on this.
    int32_t  ToInt32Inclusive (uint32_t u, int32_t limit);          // Returns [0, limit]. Must have limit >= 0.
    int32_t  ToInt32Signed    (uint32_t u, int32_t limit);          // Returns [-limit, limit]. Must have limit >= 0.
    int32_t  ToInt32          (uint32_t u, int32_t a, int32_t b);   // Returns [a, b - 1]. Must have a <= b.
    int32_t  ToInt32Inclusive (uint32_t u, int32_t a, int32_t b);   // Returns [a, b]. Must have a <= b.

    float    ToFloat      (uint32_t u);                             // Returns a float in the range [ 0, 1]
    float    ToFloatSigned(uint32_t u);                             // Returns a float in the range [-1, 1]
    float    ToFloat      (uint32_t u, float a);                    // Returns a float in the range [ 0, a]
    float    ToFloatSigned(uint32_t u, float a);                    // Returns a float in the range [-a, a]
    float    ToFloat      (uint32_t u, float a, float b);           // Returns a float in the range [ a, b]

    double   ToDouble      (uint32_t u);                            // Returns a double in the range [ 0, 1]
    double   ToDoubleSigned(uint32_t u);                            // Returns a double in the range [-1, 1]
    double   ToDouble      (uint32_t u, double a);                  // Returns a double in the range [ 0, a]
    double   ToDoubleSigned(uint32_t u, double a);                  // Returns a double in the range [-a, a]
    double   ToDouble      (uint32_t u, double a, double b);        // Returns a double in the range [ a, b]

    // 2D
    Vec2f ToVec2      (uint32_t u);                                 // Returns random [0, 1]^2 vector
    Vec2f ToVec2Signed(uint32_t u);                                 // Returns random [-1, 1]^2 vector
    Vec2f ToDir2      (uint32_t u);                                 // Returns random 2d direction
    Vec2f ToDir2      (uint32_t u, float s);                        // Returns random 2d direction around y, spread determined by 's'. (0 = y, 0.5 = positive, 1 = all directions)
    Vec2f ToSquare    (uint32_t u);                                 // Returns random point in the unit square: [0, 1]^2
    Vec2f ToRectangle (uint32_t u, Vec2f min, Vec2f max);           // Returns random point between min and max
    Vec2f ToTriangle  (uint32_t u);                                 // Returns random point in triangle with vertices (1, 0), (0, 0), (1, 0).
    Vec2f ToTriangle  (uint32_t u, Vec2f v0, Vec2f v1, Vec2f v2);   // Returns random point in triangle with vertices v0, v1, v2.
    Vec2f ToCircle    (uint32_t u);                                 // Returns random point on the unit circle
    Vec2f ToDisc      (uint32_t u);                                 // Returns random point on or inside the unit circle
    Vec2f ToRing      (uint32_t u, float r);                        // Returns random point in the unit ring with given width. r=1 -> full circle.
    Mat2f ToMat2      (uint32_t u);                                 // Returns random [0, 1]^2^2 matrix
    Mat2f ToMat2Signed(uint32_t u);                                 // Returns random [-1, 1]^2^2 matrix
    Mat2f ToRot2      (uint32_t u);                                 // Returns random rotation, can be used with either row or column vectors.

    // 3D
    Vec3f ToVec3      (uint32_t u);                                 // Returns random [0, 1]^3 vector
    Vec3f ToVec3Signed(uint32_t u);                                 // Returns random [-1, 1]^3 vector
    Vec3f ToDir3      (uint32_t u);                                 // Returns random 3d direction
    Vec3f ToDir3      (uint32_t u, float s);                        // Returns random direction around z, with spread determined by 's'. (0 = z, 0.5 = hemisphere, 1 = sphere)
    Vec3f ToCube      (uint32_t u);                                 // Returns random point in the unit cube ([-1, 1]^3)
    Vec3f ToBox       (uint32_t u, Vec3f min, Vec3f max);           // Returns random point between min and max
    Vec3f ToTriangle  (uint32_t u, Vec3f v0, Vec3f v1, Vec3f v2);   // Returns random point in triangle with vertices v0, v1, v2.
    Vec3f ToSphere    (uint32_t u);                                 // Returns random point on the unit sphere
    Vec3f ToBall      (uint32_t u);                                 // Returns random point on or inside the unit sphere
    Vec3f ToEllipsoid (uint32_t u, Vec3f min, Vec3f max);           // Returns random point in the ellipsoid defined by min and max.
    Vec3f ToTorus     (uint32_t u, float r);                        // Returns random point from a torus of the given radial width. r=0 -> circle, r=1 -> no hole in the middle.
    Mat3f ToMat3      (uint32_t u);                                 // Returns random [0, 1]^3^3 matrix
    Mat3f ToMat3Signed(uint32_t u);                                 // Returns random [-1, 1]^3^3 matrix
    Mat3f ToRot3      (uint32_t u);                                 // Returns random rotation, can be used with either row or column vectors.
    Vec4f ToQuat      (uint32_t u);                                 // Returns random quaternion rotation

    // 2D u x 2 versions. These are costlier and less convenient (requiring multiple rng inputs), but higher quality, and ensure the character of the generator is reflected across all dimensions.
    Vec2f ToVec2      (uint32_t u0, uint32_t u1);                                // Returns random [0, 1]^3 vector
    Vec2f ToVec2Signed(uint32_t u0, uint32_t u1);                                // Returns random [-1, 1]^3 vector
    Vec2f ToSquare    (uint32_t u0, uint32_t u1);                                // Returns random point in the unit square [0, 1]^2
    Vec2f ToRectangle (uint32_t u0, uint32_t u1, Vec2f min, Vec2f max);          // Returns random point between min and max
    Vec2f ToTriangle  (uint32_t u0, uint32_t u1);                                // Returns random point in triangle with vertices (1, 0), (0, 0), (1, 0).
    Vec2f ToTriangle  (uint32_t u0, uint32_t u1, Vec2f v0, Vec2f v1, Vec2f v2);  // Returns random point in triangle with vertices v0, v1, v2.
    Vec2f ToCircle    (uint32_t u0);                                             // Returns random direction vector (i.e., point on the unit sphere surface).
    Vec2f ToDisc      (uint32_t u0, uint32_t u1);                                // Returns random point in the unit circle. u0 affects angular distribute, u1 the radial. E.g., ToDisc(r1, ModTriangle(r2)) gives more samples in the centre.
    Vec2f ToRing      (uint32_t u0, uint32_t u1, float r);                       // Returns random point in the unit ring with given width. r=1 -> full circle.

    // 3D full-u versions. These are costlier and less convenient (requiring multiple rng inputs), but higher quality, and ensure the character of the generator is reflected across all dimensions.
    Vec3f ToVec3      (uint32_t u0, uint32_t u1, uint32_t u2);                        // Returns random [0, 1]^3 vector
    Vec3f ToVec3Signed(uint32_t u0, uint32_t u1, uint32_t u2);                        // Returns random [-1, 1]^3 vector
    Vec3f ToDir3      (uint32_t u0, uint32_t u1);                                     // Returns random 3d direction
    Vec3f ToDir3      (uint32_t u0, uint32_t u1, float s);                            // Returns random direction around z, with spread determined by 's'. (0 = z, 0.5 = hemisphere, 1 = sphere)
    Vec3f ToCube      (uint32_t u0, uint32_t u1, uint32_t u2);                        // Returns point vector in the unit cube ([-1, 1]^3)
    Vec3f ToBox       (uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max);  // Returns random point between min and max
    Vec3f ToTriangle  (uint32_t u0, uint32_t u1, Vec3f v0, Vec3f v1, Vec3f v2);       // Returns random point in triangle with vertices v0, v1, v2.
    Vec3f ToSphere    (uint32_t u0, uint32_t u1);                                     // Returns random point on the unit sphere
    Vec3f ToBall      (uint32_t u0, uint32_t u1, uint32_t u2);                        // Returns random point on or inside the unit sphere
    Vec3f ToEllipsoid (uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max);  // Returns random point in the ellipsoid defined by min and max.
    Vec3f ToTorus     (uint32_t u0, uint32_t u1, uint32_t u2, float r);               // Returns random point from a torus with the given radius and bounds
    Mat3f ToRot3      (uint32_t u0, uint32_t u1, uint32_t u2);                        // Returns random rotation, can be used with either row or column vectors.
    Vec4f ToQuat      (uint32_t u0, uint32_t u1, uint32_t u2, uint32_t u3);           // Returns random quaternion rotation. Use FastRenormalize to improve error if necessary.

    // Modifiers. These can be used to pre-warp the input distribution before calling the above.
    uint32_t ModTriangle (uint32_t u);      // Returns a triangle-shaped distribution. Useful in combination with other distributors, e.g., ToFloat(ModTriangle(u)), or ToSquare(ModTriangle(u0), ModTriangle(u1)).
    uint32_t ModGaussLike(uint32_t u);      // Returns a Gaussian-like distribution with range=[0,1], mean=0.5, var=1/12. This is much cheaper than a proper Gaussian (for which see ToFloatGaussian). The tails are limited to the range, so ToFloat(ToGaussLike(u), a, b) will always return a number between a and b, which can be useful.

    uint32_t ModWeighted (uint32_t u, int numWeights, const int   weights[]);   // Returns distribution weighted to 'numWeights' bins, according to 'weights'. A weight can be zero (no samples), but at least one weight must be non-zero.
    uint32_t ModWeighted (uint32_t u, int numWeights, const float weights[]);   // Version of ModWeighted that takes float weights for convenience.

    uint32_t ModInvert   (uint32_t u);      // Inverts distribution
    uint32_t ModHalfUp   (uint32_t u);      // Expands the left side of a symmetric distribution, e.g., a triangle distribution is converted to a ramp up.
    uint32_t ModHalfDown (uint32_t u);      // Expands the right side of a symmetric distribution, e.g., a triangle distribution is converted to a ramp down.
    uint32_t ModSymmUp   (uint32_t u);      // Converts a ramp back to a symmetric distribution, inverse of ModHalfUp.
    uint32_t ModSymmDown (uint32_t u);      // Converts a ramp back to a symmetric distribution, inverse of ModHalfDown.

    // Modifier to add jitter to input sequence
    struct ModJitter
    {
        uint32_t mS;  // scale
        uint32_t mJ;  // LCG state

        ModJitter(float lambda, uint32_t seed = 12345);  // lambda = max displacement from 0-1.
        uint32_t operator()(uint32_t u);
    };

    // Helpers for specific distribution types
    ModJitter GridJitter(int dims, int count, float lambda);  // For use with general grid-like inputs
    ModJitter RdJitter  (int dims, int count, float lambda);  // For use with R[N]U

    // Specialty or convenience operations.
    float   ToFloatTriangle (uint32_t u, float a, float b); // Most useful variant of ModTriangle -- returns a triangle-shaped distribution in [a, b]. Shorthand for ToFloat(ModTriangle(u), a, b).
    float   ToFloatGaussLike(uint32_t u, float a, float b); // Most useful variant of ModGaussLike -- returns a Gaussian-like distribution in [a, b].

    int32_t ToInt32Weighted(uint32_t u, int numWeights, const float weights[]); // Returns [0, numWeights - 1], distributed according to 'weights'. Cheaper version of ToInt32(ToWeighted(u, ...), numWeights)

    float   ToFloatGaussian(uint32_t u);                    // Returns the full Gaussian ("normal") distribution with infinite tails.
    float   ToFloatGaussian(uint32_t u, float mean, float stdDev);   // Returns a Gaussian distribution with the given mean and standard deviation
    Vec2f   ToFloatGaussian(uint32_t u0, uint32_t u1);      // Full version of ToFloatGaussian via the Box-Muller transform, returns two normally-distributed samples for two inputs.

    float   ToFloatCauchy(uint32_t u, float mean = 0, float gamma = 1.0f);  // Returns a Cauchy distribution parameterised by 'mean' and 'gamma'. Has longer tails than Gaussian.

    float   GaussianWalk(uint32_t u, float previous, float mean, float stdDev, float theta = 0.0f);  // These return the next value in a random walk (Wiener process). If theta is non-zero, it controls the strength of mean-reversion, as per the Ornstein–Uhlenbeck process
    Vec2f   GaussianWalk(uint32_t u0, uint32_t u1, Vec2f previous, Vec2f mean, float stdDev, float theta = 0.0f);

    float   NextEventTime(uint32_t u, float rate);  // Returns random time until next event, given an average event rate per unit time of 'rate'. This uses the exponential distribution, assuming independent events.

    // Alternative mapping functions, whose different characteristics might be more useful in some situations
    Vec2f   ToTriangleHeitz(uint32_t u0, uint32_t u1);  // Uses the low distortion square->tri map from https://eheitzresearch.wordpress.com/749-2/.
                                                        // The only reason this isn't the default is that the map's discontinuity can be visible with some structured inputs. If this isn't an issue then it should be preferred
    Vec2f   ToTriangleBakuOwen   (uint32_t u);          // Distributes u hierarchically in the triangle. (See Basu & Owen '14.) Can be more evenly spaced than default version but is more expensive. Relies on randomness in all bits of 'u', so doesn't play well with stratified sequences (Golden/Halton/Rd).
    Vec2f   ToTriangleBakuOwenRev(uint32_t u);          // Distributes u taking high bits first rather than low bits. Works better with some stratified inputs.
}


// --------------------------------------------------------------------------
// Inlines
// --------------------------------------------------------------------------


// Internal helpers

namespace DL { namespace Internal
{
    constexpr float  kFloatFromUInt32Scale  = float (1.0 / 0xFFFFFFFF);
    constexpr float  kFloatFromInt32Scale   = float (2.0 / 0xFFFFFFFF);
    constexpr double kDoubleFromUInt32Scale = double(1.0 / 0xFFFFFFFF);
    constexpr double kDoubleFromInt32Scale  = double(2.0 / 0xFFFFFFFF);

    inline uint32_t Next(uint32_t u)
    {
        u = u * 1664525 + 1013904223;  // use LCG/XOR hash combo to avoid correlation with structured inputs
        u ^= (u << 13);                // LCG only is okay (and faster) for just pseudo-random sequences
        u ^= (u >> 17);
        u ^= (u << 5);
        return u;
    }
}}

// uint32_t variants

inline uint32_t DL::ToUInt32(uint32_t x)
{
    return x;
}

inline uint32_t DL::ToUInt32(uint32_t u, uint32_t limit)
{
    return (u * (uint64_t(limit))) >> 32;
}

inline uint32_t DL::ToUInt32(uint32_t u, uint32_t a, uint32_t b)
{
    DL_ASSERT(a <= b);
    return a + ToUInt32(u, b - a);
}

inline uint32_t DL::ToUInt32Inclusive(uint32_t u, uint32_t limit)
{
    return (u * (uint64_t(limit) + 1)) >> 32;
}

inline uint32_t DL::ToUInt32Inclusive(uint32_t u, uint32_t a, uint32_t b)
{
    DL_ASSERT(a <= b);
    return a + DL::ToUInt32Inclusive(u, b - a);
}

// int32_t variants

inline int32_t DL::ToInt32(uint32_t u)
{
    return int32_t(u & 0x7FFFFFFF);
}

inline int32_t DL::ToInt32Signed(uint32_t u)
{
    return int32_t(u ^ 0x80000000);
}

inline int32_t DL::ToInt32Inclusive(uint32_t u, int32_t limit)
{
    return int32_t((u * (uint64_t(limit) + 1)) >> 32);
}

inline int32_t DL::ToInt32(uint32_t u, int32_t limit)
{
    DL_ASSERT(limit >= 0);
    return int32_t((u * (uint64_t(limit))) >> 32);
}

inline int32_t DL::ToInt32Signed(uint32_t u, int32_t limit)
{
    DL_ASSERT(limit >= 0);
    return ((u * uint64_t(2 * limit + 1)) >> 32) - limit;
}

inline int32_t DL::ToInt32(uint32_t u, int32_t a, int32_t b)
{
    DL_ASSERT(a <= b);
    return a + ToUInt32(u, b - a);
}

inline int32_t DL::ToInt32Inclusive(uint32_t u, int32_t a, int32_t b)
{
    DL_ASSERT(a <= b);
    return a + ToUInt32Inclusive(u, b - a);
}

// float variants
inline float DL::ToFloat(uint32_t u)
{
    return Internal::kFloatFromUInt32Scale * u;
}

inline float DL::ToFloatSigned(uint32_t u)
{
    return Internal::kFloatFromInt32Scale * u - float(1);
}

inline float DL::ToFloat(uint32_t u, float a)
{
    return a * (Internal::kFloatFromUInt32Scale * u);
}

inline float DL::ToFloatSigned(uint32_t u, float a)
{
    return a * (Internal::kFloatFromInt32Scale * u - float(1));
}

inline float DL::ToFloat(uint32_t u, float a, float b)
{
    return a + (b - a) * (Internal::kFloatFromUInt32Scale * u);
}

// double variants
inline double DL::ToDouble(uint32_t u)
{
    return Internal::kDoubleFromUInt32Scale * u;
}

inline double DL::ToDoubleSigned(uint32_t u)
{
    return Internal::kDoubleFromInt32Scale * u - 1.0;
}

inline double DL::ToDouble(uint32_t u, double a)
{
    return a * (Internal::kDoubleFromUInt32Scale * u);
}

inline double DL::ToDoubleSigned(uint32_t u, double a)
{
    return a * (Internal::kDoubleFromInt32Scale * u - 1.0);
}

inline double DL::ToDouble(uint32_t u, double a, double b)
{
    return a + (b - a) * (Internal::kDoubleFromUInt32Scale * u);
}


// Specialty

inline float DL::ToFloatTriangle(uint32_t u, float a, float b)
{
    return ToFloat(ModTriangle(u), a, b);
}

inline float DL::ToFloatGaussLike(uint32_t u, float a, float b)
{
    return ToFloat(ModGaussLike(u), a, b);
}

inline float DL::ToFloatGaussian(uint32_t u, float mean, float stdDev)
{
    return ToFloatGaussian(u) * stdDev + mean;
}

inline float DL::ToFloatCauchy(uint32_t u, float x0, float g)
{
    return x0 + g * tanf(vlf_halfPi * (ToFloatSigned(u)));
}

inline float DL::GaussianWalk(uint32_t u, float previous, float mean, float stdDev, float theta)
{
    return previous + ToFloatGaussian(u) * stdDev + (mean - previous) * theta;
}

inline Vec2f DL::GaussianWalk(uint32_t u0, uint32_t u1, Vec2f previous, Vec2f mean, float stdDev, float theta)
{
    return previous + ToFloatGaussian(u0, u1) * stdDev + (mean - previous) * theta;
}

inline float DL::NextEventTime(uint32_t u, float rate)
{
    float uf = ToFloat(u);
    return uf > 0.0f ? -logf(uf) / rate : 50.0f / rate;
}

// 2D

inline Vec2f DL::ToVec2(uint32_t u)
{
    return Vec2f(ToFloat(u), ToFloat(Internal::Next(u)));
}

inline Vec2f DL::ToVec2Signed(uint32_t u)
{
    return Vec2f(ToFloatSigned(u), ToFloatSigned(Internal::Next(u)));
}

inline Vec2f DL::ToSquare(uint32_t u)
{
    return Vec2f(ToFloat(u), ToFloat(Internal::Next(u)));
}

inline Vec2f DL::ToRectangle(uint32_t u, Vec2f min, Vec2f max)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);

    return ToRectangle(u0, u1, min, max);
}

inline Vec2f DL::ToTriangle(uint32_t u)
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

inline Vec2f DL::ToTriangle(uint32_t u, Vec2f v0, Vec2f v1, Vec2f v2)
{
    Vec2f c = ToTriangle(u);
    return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
}

inline Vec2f DL::ToRing(uint32_t u, float r)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);

    return ToRing(u0, u1, r);
}

inline Mat2f DL::ToMat2(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);

    return Mat2f(ToVec2(u0), ToVec2(u1));
}

inline Mat2f DL::ToMat2Signed(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);

    return Mat2f(ToVec2Signed(u0), ToVec2Signed(u1));
}

inline Mat2f DL::ToRot2(uint32_t u)
{
    Mat2f result;
    result.x = ToCircle(u);
    result.y = cross(result.x);
    return result;
}

// 3D

inline Vec3f DL::ToVec3(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return Vec3f(ToFloat(u0), ToFloat(u1), ToFloat(u2));
}

inline Vec3f DL::ToVec3Signed(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return Vec3f(ToFloatSigned(u0), ToFloatSigned(u1), ToFloatSigned(u2));
}

inline Vec3f DL::ToDir3(uint32_t u, float s)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);

    return ToDir3(u0, u1, s);
}

inline Vec3f DL::ToCube(uint32_t u)
{
    return ToVec3Signed(u);
}

inline Vec3f DL::ToBox(uint32_t u, Vec3f min, Vec3f max)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return ToBox(u0, u1, u2, min, max);
}

inline Vec3f DL::ToTriangle(uint32_t u, Vec3f v0, Vec3f v1, Vec3f v2)
{
    Vec2f c = ToTriangle(u);
    return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
}

inline Vec3f DL::ToSphere(uint32_t u)
{
    return ToDir3(u);
}

inline Mat3f DL::ToMat3(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return Mat3f(ToVec3(u0), ToVec3(u1), ToVec3(u2));
}

inline Mat3f DL::ToMat3Signed(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return Mat3f(ToVec3Signed(u0), ToVec3Signed(u1), ToVec3Signed(u2));
}

inline Mat3f DL::ToRot3(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);

    return ToRot3(u0, u1, u2);
}

inline Vec4f DL::ToQuat(uint32_t u)
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);
    uint32_t u3 = Internal::Next(u2);

    return ToQuat(u0, u1, u2, u3);
}

// 2D Full

inline Vec2f DL::ToVec2(uint32_t u0, uint32_t u1)
{
    return Vec2f(ToFloat(u0), ToFloat(u1));
}

inline Vec2f DL::ToVec2Signed(uint32_t u0, uint32_t u1)
{
    return Vec2f(ToFloatSigned(u0), ToFloatSigned(u1));
}

inline Vec2f DL::ToSquare(uint32_t u0, uint32_t u1)
{
    return ToVec2Signed(u0, u1);
}

inline Vec2f DL::ToRectangle(uint32_t u0, uint32_t u1, Vec2f min, Vec2f max)
{
    return Vec2f(ToFloat(u0, min.x, max.x), ToFloat(u1, min.y, max.y));
}

inline Vec2f DL::ToTriangle(uint32_t u0, uint32_t u1)
{
    float x = ToFloat(u0);
    float y = ToFloat(u1);

    if (x + y > 1.0f)
        return { 1.0f - x, 1.0f - y };
    else
        return { x, y };
}

inline Vec2f DL::ToTriangle(uint32_t u0, uint32_t u1, Vec2f v0, Vec2f v1, Vec2f v2)
{
    Vec2f c = ToTriangle(u0, u1);
    return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
}

// 3D Full

inline Vec3f DL::ToVec3(uint32_t u0, uint32_t u1, uint32_t u2)
{
    return Vec3f(ToFloat(u0), ToFloat(u1), ToFloat(u2));
}

inline Vec3f DL::ToVec3Signed(uint32_t u0, uint32_t u1, uint32_t u2)
{
    return Vec3f(ToFloatSigned(u0), ToFloatSigned(u1), ToFloatSigned(u2));
}

inline Vec3f DL::ToCube(uint32_t u0, uint32_t u1, uint32_t u2)
{
    return ToVec3Signed(u0, u1, u2);
}

inline Vec3f DL::ToBox(uint32_t u0, uint32_t u1, uint32_t u2, Vec3f min, Vec3f max)
{
    return Vec3f(ToFloat(u0, min.x, max.x), ToFloat(u1, min.y, max.y), ToFloat(u2, min.z, max.z));
}

inline Vec3f DL::ToTriangle(uint32_t u0, uint32_t u1, Vec3f v0, Vec3f v1, Vec3f v2)
{
    Vec2f c = ToTriangle(u0, u1);
    return (1.0f - c.x - c.y) * v0 + c.x * v1 + c.y * v2;
}

inline Vec3f DL::ToSphere(uint32_t u0, uint32_t u1)
{
    return ToDir3(u0, u1);
}

// Modifiers

#if defined(DIST_MT_NO_SQRT)
inline uint32_t DL::ModTriangle(uint32_t u)  // Classic v0 + v1 approach, but, doesn't play well with structured 'u', and often sqrt() approach is no slower.
{
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t carry = (u0 & u1 & 1);  // Variance is 1/6

    return ((u0 >> 1) + (u1 >> 1) + carry);
}
#elif defined(DIST_MT_USE_SQRTI)
inline uint32_t DL::ModTriangle(uint32_t u)
{
    u = sqrti(u);  // only worth it if fast (hardware) sqrti available e.g. arm64, and willing to drop some precision
    u = (u << 16) | u;
    return ModSymmUp(u);
}
#else
inline uint32_t DL::ModTriangle(uint32_t u)  // sqrt approach using double to avoid dropping bits
{
    double uf = double(u) / UINT32_MAX;
    u = uint32_t(UINT32_MAX * sqrt(uf));
    return ModSymmUp(u);
}
#endif

inline uint32_t DL::ModGaussLike(uint32_t u)
{
    // Always use sum approach as inverse CDF for Irwin Hall distribution gets more complicated. For structured 'u'
    // ToFloatGaussian is an alternative.
    uint32_t u0 = u;
    uint32_t u1 = Internal::Next(u0);
    uint32_t u2 = Internal::Next(u1);
    uint32_t u3 = Internal::Next(u2);
#if 0
    uint32_t c0 = (u0 & u1 & 1);
    uint32_t c1 = (u2 & u3 & 1);

    uint32_t carry = c0 + c1 + (c0 & c1 & 1);
#else
    uint32_t carry = 1;  // because this situation is symmetrical we can use a constant at the expense of never getting 0 or UINT32_MAX
#endif
    uint32_t r = (u0 >> 2) + (u1 >> 2) + (u2 >> 2) + (u3 >> 2) + carry;

    return r;
}

inline uint32_t DL::ModInvert(uint32_t u)
{
    return ~u;
}

inline uint32_t DL::ModHalfDown(uint32_t u)
{
    return (u & 0x80000000) ? (u << 1) : ~(u << 1);
}

inline uint32_t DL::ModHalfUp(uint32_t u)
{
    return (u & 0x80000000) ? ~(u << 1) : (u << 1);
}

inline uint32_t DL::ModSymmUp(uint32_t u)
{
    return (u & 1) ? ~(u >> 1) : (u >> 1);
}

inline uint32_t DL::ModSymmDown(uint32_t u)
{
    return (u & 1) ? (u << 1) : ~(u << 1);
}

inline DL::ModJitter DL::GridJitter(int dim, int n, float lambda)
{
    // Note: this dynamic N-agnostic approach mentioned in link above doesn't really work that well,
    // early sample offsets are far too large
    // uint64_t si = uint64_t(mS / sqrtf(mJ.mIndex / n + 0.3f));

    return ModJitter(lambda / pow((n < 1 ? 1 : n), 1.0 / dim));
}

inline DL::ModJitter DL::RdJitter(int dim, int n, float lambda)
{
    // See http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/
    constexpr double kRdF32 = 0.76;  // our setup assumes lambda=0.5 is ideal so we don't fold that in here
    lambda *= kRdF32 / pow(n - 0.7, 1.0 / dim);
    return ModJitter(lambda < 1.0f ? lambda : 1.0f);
}

inline DL::ModJitter::ModJitter(float lambda, uint32_t seed) :
    mS(uint32_t(UINT32_MAX * lambda)),
    mJ(seed)
{
    DL_ASSERT(lambda <= 1.0f);
}

inline uint32_t DL::ModJitter::operator()(uint32_t u)
{
    mJ = Internal::Next(mJ);
    return u + ((uint64_t(mS) * mJ) >> 32) - (mS >> 1);
}

// Alternatives

inline Vec2f DL::ToTriangleHeitz(uint32_t u0, uint32_t u1)
{
    Vec2f t = 0.5f * ToVec2(u0, u1);
    float offset = t.y - t.x;

    if (offset > 0)
        t.y += offset;
    else
        t.x -= offset;

    return t;
}

#endif
