//
// Generate.cpp
//
// Various sample/number generators
//
// Andrew Willmott
//

#include "Generate.hpp"

#include <stdio.h>

#define HL_ERROR(X)
using namespace DL;

namespace
{
    constexpr float kOneOverThree = float(1.0 / 3.0);
    constexpr float kOneOverFive  = float(1.0 / 5.0);

    const uint32_t kMaxBase3 = 43046721 - 1;  // 3^16, 16 = max 2-bit digits. -1 as need to avoid last trailing carry
    const uint32_t kMaxBase5 = 9765625 - 1;   // 5^10, 10 = max 3-bit digits. -1 as need to avoid last trailing carry

    const uint64_t kMaxBase3d = 1853020188851841 - 1;  // 3^32, 16 = max 2-bit digits. -1 as need to avoid last trailing carry
    const uint64_t kMaxBase5d = 476837158203125 - 1;   // 5^21, 10 = max 3-bit digits. -1 as need to avoid last trailing carry
}

void PCG::Advance(int deltaInt)
{
    uint64_t delta = deltaInt;
    uint64_t cur_mult = HL_PCG_MULTIPLIER_64;
    uint64_t cur_plus = mInc;

    uint64_t acc_mult = 1u;
    uint64_t acc_plus = 0u;

    while (delta > 0)
    {
        if (delta & 1)
        {
            acc_mult *= cur_mult;
            acc_plus = acc_plus * cur_mult + cur_plus;
        }

        cur_plus = (cur_mult + 1) * cur_plus;
        cur_mult *= cur_mult;

        delta >>= 1;
    }

    mState = acc_mult * mState + acc_plus;
}

float DL::HaltonFloat(int n, int b)
// Return term i of the base b Halton sequence.
// This is just a generalization of Heckbert's bit reversal distribution trick.
// E.g., when b=3, write n as a base 3 number, digit 0 -> which third of interval the
// sample is in, 1 -> which third of that, 2 -> which third of that, etc.
{
    float result = 0;
    float ip = 1.0f / b;
    float p = ip;

    while (n > 0)
    {
        result += (n % b) * p;
        n = n / b;
        p *= ip;
    }

    return result;
}

uint32_t DL::HaltonUInt32(int n, int b)
{
    uint32_t result = 0;
    uint32_t p  = UINT32_MAX / b;

    while (n > 0)
    {
        result += (n % b) * p;
        n /= b;
        p /= b;
    }

    return result;
}

double DL::HaltonDouble(int n, int b)
{
    double result = 0;
    double ip = 1.0 / b;
    double p = ip;

    while (n > 0)
    {
        result += (n % b) * p;
        n = n / b;
        p *= ip;
    }

    return result;
}

double DL::Halton2Double(int i)
{
    uint32_t u = Halton2UInt32(i);

    uint64_t resultU(0x3ff0000000000000 | (uint64_t(u) << 20));
    double result = ((double&) resultU);
    return result - 1.0;
}


// --- Halton2 --------------------------------------------------------

int Halton2::Next()
{
    if (mBase2 >= kMaxBase3)
    {
        HL_ERROR("Overflow");
        return mBase2;
    }

    /////////////////////////////////////
    // base 2

    uint32_t oldBase2 = mBase2;
    mBase2++;
    uint32_t diff = mBase2 ^ oldBase2;

    // bottom bit always changes, higher bits
    // change less frequently.
    float s = 0.5f;

    // diff will be of the form 0 * 1 +, i.e. one bits up until the last carry.
    // expected iterations = 1 + 0.5 + 0.25 + ... = 2
    do
    {
        if (oldBase2 & 1)
            mU[0] -= s;
        else
            mU[0] += s;

        s *= 0.5f;

        diff = diff >> 1;
        oldBase2 = oldBase2 >> 1;
    }
    while (diff);


    /////////////////////////////////////
    // base 3: use 2 bits for each base 3 digit.

    uint32_t mask = 0x3;  // also the max base 3 digit
    uint32_t add  = 0x1;  // amount to add to force carry once digit==3
    s = kOneOverThree;

    mBase3++;

    // expected iterations: 1.5
    while (1)
    {
        if ((mBase3 & mask) == mask)
        {
            mBase3 += add;          // force carry into next 2-bit digit
            mU[1] -= 2 * s;

            mask = mask << 2;
            add  = add  << 2;

            s *= kOneOverThree;
        }
        else
        {
            mU[1] += s;     // we know digit n has gone from a to a + 1
            break;
        }
    };

    return mBase2; // return the index of this sequence point
}

void Halton2::Set(int n)
{
    if (mBase2 >= kMaxBase3)
    {
        HL_ERROR("Overflow");
        return;
    }

    mU[0] = Halton2Float(n);
    mBase2 = n;

    mU[1] = 0;
    mBase3 = 0;

    float ip = kOneOverThree;
    float p = ip;

    for (int i = 0, k = n; k; i += 2, k /= 3)
    {
        int d = (k % 3);
        mBase3 |= d << i;
        mU[1] += d * p;
        p *= ip;
    }
}


// --- Halton3 --------------------------------------------------------

int Halton3::Next()
{
    if (mBase2 >= kMaxBase5)
    {
        HL_ERROR("Overflow");
        return mBase2;
    }

    // base 2: 1 bit per digit
    uint32_t oldBase2 = mBase2;
    mBase2++;
    uint32_t diff = mBase2 ^ oldBase2;

    // bottom bit always changes, higher bits
    // change less frequently.
    float s = 0.5f;

    // diff will be of the form 0 * 1 + , i.e. one bits up until the last carry.
    // expected iterations = 1 + 0.5 + 0.25 + ... = 2
    do
    {
        if (oldBase2 & 1)
            mU[0] -= s;
        else
            mU[0] += s;

        s *= 0.5f;

        diff = diff >> 1;
        oldBase2 = oldBase2 >> 1;
    }
    while (diff);


    // base 3: use 2 bits for each base 3 digit.
    uint32_t mask = 0x3;  // also the max base 3 digit
    uint32_t add  = 0x1;  // amount to add to force carry once digit==3
    s = kOneOverThree;

    mBase3++;

    // expected iterations: 1.5
    while (1)
    {
        if ((mBase3 & mask) == mask)
        {
            mBase3 += add;          // force carry into next 2-bit digit
            mU[1] -= 2 * s;

            mask = mask << 2;
            add  = add  << 2;

            s *= kOneOverThree;
        }
        else
        {
            mU[1] += s;     // we know digit n has gone from a to a + 1
            break;
        }
    };

    // base 5: use 3 bits for each base 5 digit.
    mask = 0x7;
    add  = 0x3;  // amount to add to force carry once digit==dmax
    uint32_t dmax = 0x5;  // max digit

    s = kOneOverFive;

    mBase5++;

    // expected iterations: 1.25
    while (1)
    {
        if ((mBase5 & mask) == dmax)
        {
            mBase5 += add;          // force carry into next 3-bit digit
            mU[2] -= 4 * s;

            mask = mask << 3;
            dmax = dmax << 3;
            add  = add  << 3;

            s *= kOneOverFive;
        }
        else
        {
            mU[2] += s;     // we know digit n has gone from a to a + 1
            break;
        }
    };

    return mBase2; // return the index of this sequence point
}

void Halton3::Set(int n)
{
    if (n >= (int) kMaxBase5)
    {
        HL_ERROR("Overflow");
        return;
    }

    mU[0] = Halton2Float(n);
    mBase2 = n;

    mU[1] = 0.0f;
    mBase3 = 0;

    float p = kOneOverThree;

    for (int i = 0, k = n; k; i += 2, k /= 3)
    {
        int d = (k % 3);
        mBase3 |= d << i;
        mU[1] += d * p;
        p *= kOneOverThree;
    }

    mU[2] = 0.0f;
    mBase5 = 0;

    p = kOneOverFive;

    for (int i = 0, k = n; k; i += 3, k /= 5)
    {
        int d = (k % 5);
        mBase5 |= d << i;
        mU[2] += d * p;
        p *= kOneOverFive;
    }
}

// Halton2D

uint64_t Halton2D::Next()
{
    if (mBase2 >= kMaxBase3d)
    {
        HL_ERROR("Overflow");
        return mBase2;
    }

    /////////////////////////////////////
    // base 2

    uint64_t oldBase2 = mBase2;
    mBase2++;
    uint64_t diff = mBase2 ^ oldBase2;

    // bottom bit always changes, higher bits
    // change less frequently.
    double s = 0.5;

    // diff will be of the form 0 * 1 +, i.e. one bits up until the last carry.
    // expected iterations = 1 + 0.5 + 0.25 + ... = 2
    do
    {
        if (oldBase2 & 1)
            mU[0] -= s;
        else
            mU[0] += s;

        s *= 0.5;

        diff = diff >> 1;
        oldBase2 = oldBase2 >> 1;
    }
    while (diff);


    /////////////////////////////////////
    // base 3: use 2 bits for each base 3 digit.

    uint32_t mask = 0x3;  // also the max base 3 digit
    uint32_t add  = 0x1;  // amount to add to force carry once digit==3
    s = kOneOverThree;

    mBase3++;

    // expected iterations: 1.5
    while (1)
    {
        if ((mBase3 & mask) == mask)
        {
            mBase3 += add;          // force carry into next 2-bit digit
            mU[1] -= 2 * s;

            mask = mask << 2;
            add  = add  << 2;

            s *= kOneOverThree;
        }
        else
        {
            mU[1] += s;     // we know digit n has gone from a to a + 1
            break;
        }
    };

    return mBase2; // return the index of this sequence point
}

void Halton2D::Set(int n)
{
    if (mBase2 >= kMaxBase3d)
    {
        HL_ERROR("Overflow");
        return;
    }

    mU[0] = Halton2Double(n);
    mBase2 = n;

    mU[1] = 0;
    mBase3 = 0;

    double ip = 1.0 / 3.0;
    double p = ip;

    for (int i = 0, k = n; k; i += 2, k /= 3)
    {
        uint64_t d = (k % 3);
        mBase3 |= d << i;
        mU[1] += d * p;
        p *= ip;
    }
}

// Halton3d

uint64_t Halton3D::Next()
{
    if (mBase2 >= kMaxBase5d)
    {
        HL_ERROR("Overflow");
        return mBase2;
    }

    // base 2: 1 bit per digit
    uint64_t oldBase2 = mBase2;
    mBase2++;
    uint64_t diff = mBase2 ^ oldBase2;

    // bottom bit always changes, higher bits
    // change less frequently.
    double s = 0.5;

    // diff will be of the form 0 * 1 +, i.e. one bits up until the last carry.
    // expected iterations = 1 + 0.5 + 0.25 + ... = 2
    do
    {
        if (oldBase2 & 1)
            mU[0] -= s;
        else
            mU[0] += s;

        s *= 0.5;

        diff = diff >> 1;
        oldBase2 = oldBase2 >> 1;
    }
    while (diff);

    // base 3: use 2 bits for each base 3 digit.
    uint64_t mask = 0x3;  // also the max base 3 digit
    uint64_t add  = 0x1;  // amount to add to force carry once digit==dmax
    s = 1.0 / 3.0;

    mBase3++;

    // expected iterations: 1.5
    while (1)
    {
        if ((mBase3 & mask) == mask)
        {
            mBase3 += add;          // force carry into next 2-bit digit
            mU[1] -= 2.0 * s;

            mask = mask << 2;
            add  = add  << 2;

            s *= 1.0 / 3.0;
        }
        else
        {
            mU[1] += s;     // we know digit n has gone from a to a + 1
            break;
        }
    };

    // base 5: use 3 bits for each base 5 digit.
    uint64_t dmax;
    mask = 0x7;
    dmax = 0x5;  // max digit
    add  = 0x3;  // amount to add to force carry once digit==dmax

    s = 1.0 / 5.0;

    mBase5++;

    // expected iterations: 1.25
    while (1)
    {
        if ((mBase5 & mask) == dmax)
        {
            mBase5 += add;          // force carry into next 3-bit digit
            mU[2] -= 4.0 * s;

            mask = mask << 3;
            dmax = dmax << 3;
            add  = add  << 3;

            s *= 1.0 / 5.0;
        }
        else
        {
            mU[2] += s;     // we know digit n has gone from a to a + 1
            break;
        }
    };

    return mBase2; // return the index of this sequence point
}

void Halton3D::Set(int n)
{
    if (mBase2 >= kMaxBase5d)
    {
        HL_ERROR("Overflow");
        return;
    }

    mU[0] = Halton2Double(n);
    mBase2 = n;

    mU[1] = 0.0;
    mBase3 = 0;

    double p = 1.0 / 3.0;

    for (int i = 0, k = n; k; i += 2, k /= 3)
    {
        uint64_t d = (k % 3);
        mBase3 |= d << i;
        mU[1] += d * p;
        p *= 1.0 / 3.0;
    }

    mU[2] = 0.0;
    mBase5 = 0;

    p = 1.0 / 5.0;

    for (int i = 0, k = n; k; i += 3, k /= 5)
    {
        uint64_t d = (k % 5);
        mBase5 |= d << i;
        mU[2] += d * p;
        p *= 1.0 / 5.0;
    }
}


// --- Halton2U --------------------------------------------------------

int Halton2U::Next()
{
    if (mBase2 >= kMaxBase3)
    {
        HL_ERROR("Overflow");
        return mBase2;
    }

    mBase2++;

    uint32_t mask = 0x1;
    uint32_t s    = UINT32_MAX / 2;

    while (1) // expected iterations: 2
    {
        if ((mBase2 & mask) == 0)
        {
            mU[0] -= s;

            mask = mask << 1;

            s /= 2;
        }
        else
        {
            mU[0] += s;     // we know digit n has gone from 0 to 1
            break;
        }
    };

    mBase3++;

             mask = 0x3;  // also the max base 3 digit
    uint32_t add  = 0x1;  // amount to add to force carry once digit==3
             s    = UINT32_MAX / 3;

    while (1) // expected iterations: 1.5
    {
        if ((mBase3 & mask) == mask)
        {
            mBase3 += add;          // force carry into next 2-bit digit
            mU[1] -= 2 * s;

            mask = mask << 2;
            add  = add  << 2;

            s /= 3;
        }
        else
        {
            mU[1] += s;     // we know digit n has gone from x to x + 1
            break;
        }
    };

    return mBase2;
}

void Halton2U::Set(int n)
{
    if (mBase2 >= kMaxBase3)
    {
        HL_ERROR("Overflow");
        return;
    }

    mU[0] = Halton2UInt32(n);
    mBase2 = n;

    mU[1] = 0;
    mBase3 = 0;

    uint32_t p = UINT32_MAX / 3;

    for (int i = 0, k = n; k; i += 2, k /= 3)
    {
        int d = (k % 3);
        mBase3 |= d << i;
        mU[1] += d * p;
        p /= 3;
    }
}


// --- Halton3U --------------------------------------------------------

int Halton3U::Next()
{
    if (mBase2 >= kMaxBase5)
    {
        HL_ERROR("Overflow");
        return mBase2;
    }

    mBase2++;

    uint32_t mask = 0x1;
    uint32_t s    = UINT32_MAX / 2;

    while (1) // expected iterations: 2
    {
        if ((mBase2 & mask) == 0)
        {
            mU[0] -= s;

            mask = mask << 1;

            s /= 2;
        }
        else
        {
            mU[0] += s;     // we know digit n has gone from 0 to 1
            break;
        }
    };

    mBase3++;

             mask = 0x3;  // also the max base 3 digit
    uint32_t add  = 0x1;  // amount to add to force carry once digit==3
             s    = UINT32_MAX / 3;

    while (1) // expected iterations: 1.5
    {
        if ((mBase3 & mask) == mask)
        {
            mBase3 += add;          // force carry into next 2-bit digit
            mU[1] -= 2 * s;

            mask = mask << 2;
            add  = add  << 2;

            s /= 3;
        }
        else
        {
            mU[1] += s;     // we know digit n has gone from x to x + 1
            break;
        }
    };

    mBase5++;

    mask = 0x7;
    add  = 0x3;  // amount to add to force carry once digit==5
    uint32_t dmax = 5;
    s = UINT32_MAX / 5;

    while (1)   // expected iterations: 1.25
    {
        if ((mBase5 & mask) == dmax)
        {
            mBase5 += add;          // force carry into next 3-bit digit
            mU[2] -= 4 * s;

            mask = mask << 3;
            dmax = dmax << 3;
            add  = add  << 3;

            s /= 5;
        }
        else
        {
            mU[2] += s;     // we know digit n has gone from a to a + 1
            break;
        }
    };

    return mBase2; // return the index of this sequence point
}

void Halton3U::Set(int n)
{
    if (mBase2 >= kMaxBase5)
    {
        HL_ERROR("Overflow");
        return;
    }

    mU[0] = Halton2UInt32(n);
    mBase2 = n;

    mU[1] = 0;
    mBase3 = 0;

    uint32_t p = UINT32_MAX / 3;

    for (int i = 0, k = n; k; i += 2, k /= 3)
    {
        int d = (k % 3);
        mBase3 |= d << i;
        mU[1] += d * p;
        p /= 3;
    }

    mU[2] = 0;
    mBase5 = 0;

    p = UINT32_MAX / 5;

    for (int i = 0, k = n; k; i += 3, k /= 5)
    {
        int d = (k % 5);
        mBase5 |= d << i;
        mU[2] += d * p;
        p /= 5;
    }
}


//
// Sobel2f via Leonhard Gruenschloss: https://github.com/lgruen/sobol/tree/main/single-precision
//

namespace
{
    const uint32_t kSobolMatrices[kSobolSize * kSobolMaxDim] =
    {
        0x80000000U,
        0x40000000U,
        0x20000000U,
        0x10000000U,
        0x8000000U,
        0x4000000U,
        0x2000000U,
        0x1000000U,
        0x800000U,
        0x400000U,
        0x200000U,
        0x100000U,
        0x80000U,
        0x40000U,
        0x20000U,
        0x10000U,
        0x8000U,
        0x4000U,
        0x2000U,
        0x1000U,
        0x800U,
        0x400U,
        0x200U,
        0x100U,
        0x80U,
        0x40U,
        0x20U,
        0x10U,
        0x8U,
        0x4U,
        0x2U,
        0x1U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x0U,
        0x80000000U,
        0xc0000000U,
        0xa0000000U,
        0xf0000000U,
        0x88000000U,
        0xcc000000U,
        0xaa000000U,
        0xff000000U,
        0x80800000U,
        0xc0c00000U,
        0xa0a00000U,
        0xf0f00000U,
        0x88880000U,
        0xcccc0000U,
        0xaaaa0000U,
        0xffff0000U,
        0x80008000U,
        0xc000c000U,
        0xa000a000U,
        0xf000f000U,
        0x88008800U,
        0xcc00cc00U,
        0xaa00aa00U,
        0xff00ff00U,
        0x80808080U,
        0xc0c0c0c0U,
        0xa0a0a0a0U,
        0xf0f0f0f0U,
        0x88888888U,
        0xccccccccU,
        0xaaaaaaaaU,
        0xffffffffU,
        0x80000000U,
        0xc0000000U,
        0xa0000000U,
        0xf0000000U,
        0x88000000U,
        0xcc000000U,
        0xaa000000U,
        0xff000000U,
        0x80800000U,
        0xc0c00000U,
        0xa0a00000U,
        0xf0f00000U,
        0x88880000U,
        0xcccc0000U,
        0xaaaa0000U,
        0xffff0000U,
        0x80008000U,
        0xc000c000U,
        0xa000a000U,
        0xf000f000U,
        0x80000000U,
        0xc0000000U,
        0x60000000U,
        0x90000000U,
        0xe8000000U,
        0x5c000000U,
        0x8e000000U,
        0xc5000000U,
        0x68800000U,
        0x9cc00000U,
        0xee600000U,
        0x55900000U,
        0x80680000U,
        0xc09c0000U,
        0x60ee0000U,
        0x90550000U,
        0xe8808000U,
        0x5cc0c000U,
        0x8e606000U,
        0xc5909000U,
        0x6868e800U,
        0x9c9c5c00U,
        0xeeee8e00U,
        0x5555c500U,
        0x8000e880U,
        0xc0005cc0U,
        0x60008e60U,
        0x9000c590U,
        0xe8006868U,
        0x5c009c9cU,
        0x8e00eeeeU,
        0xc5005555U,
        0x68808000U,
        0x9cc0c000U,
        0xee606000U,
        0x55909000U,
        0x8068e800U,
        0xc09c5c00U,
        0x60ee8e00U,
        0x9055c500U,
        0xe880e880U,
        0x5cc05cc0U,
        0x8e608e60U,
        0xc590c590U,
        0x68686868U,
        0x9c9c9c9cU,
        0xeeeeeeeeU,
        0x55555555U,
        0x80000000U,
        0xc0000000U,
        0x60000000U,
        0x90000000U,
        0x80000000U,
        0xc0000000U,
        0x20000000U,
        0x50000000U,
        0xf8000000U,
        0x74000000U,
        0xa2000000U,
        0x93000000U,
        0xd8800000U,
        0x25400000U,
        0x59e00000U,
        0xe6d00000U,
        0x78080000U,
        0xb40c0000U,
        0x82020000U,
        0xc3050000U,
        0x208f8000U,
        0x51474000U,
        0xfbea2000U,
        0x75d93000U,
        0xa0858800U,
        0x914e5400U,
        0xdbe79e00U,
        0x25db6d00U,
        0x58800080U,
        0xe54000c0U,
        0x79e00020U,
        0xb6d00050U,
        0x800800f8U,
        0xc00c0074U,
        0x200200a2U,
        0x50050093U,
        0xf80f80d8U,
        0x74074025U,
        0xa20a2059U,
        0x930930e6U,
        0xd88d8878U,
        0x254254b4U,
        0x59e59e82U,
        0xe6de6dc3U,
        0x780f80a0U,
        0xb4074091U,
        0x820a20dbU,
        0xc3093025U,
        0x208d8858U,
        0x514254e5U,
        0xfbe59e79U,
        0x75de6db6U,
        0xa08f8000U,
        0x91474000U,
        0xdbea2000U,
        0x25d93000U,
        0x80000000U,
        0x40000000U,
        0x20000000U,
        0xb0000000U,
        0xf8000000U,
        0xdc000000U,
        0x7a000000U,
        0x9d000000U,
        0x5a800000U,
        0x2fc00000U,
        0xa1600000U,
        0xf0b00000U,
        0xda880000U,
        0x6fc40000U,
        0x81620000U,
        0x40bb0000U,
        0x22878000U,
        0xb3c9c000U,
        0xfb65a000U,
        0xddb2d000U,
        0x78022800U,
        0x9c0b3c00U,
        0x5a0fb600U,
        0x2d0ddb00U,
        0xa2878080U,
        0xf3c9c040U,
        0xdb65a020U,
        0x6db2d0b0U,
        0x800228f8U,
        0x400b3cdcU,
        0x200fb67aU,
        0xb00ddb9dU,
        0xf80780daU,
        0xdc09c06fU,
        0x7a05a081U,
        0x9d02d040U,
        0x5a8a2822U,
        0x2fcf3cb3U,
        0xa16db6fbU,
        0xf0b6dbddU,
        0xda8000f8U,
        0x6fc000dcU,
        0x8160007aU,
        0x40b0009dU,
        0x2288005aU,
        0xb3c4002fU,
        0xfb6200a1U,
        0xddbb00f0U,
        0x780780daU,
        0x9c09c06fU,
        0x5a05a081U,
        0x2d02d040U,
        0x80000000U,
        0x40000000U,
        0x60000000U,
        0x30000000U,
        0xc8000000U,
        0x24000000U,
        0x56000000U,
        0xfb000000U,
        0xe0800000U,
        0x70400000U,
        0xa8600000U,
        0x14300000U,
        0x9ec80000U,
        0xdf240000U,
        0xb6d60000U,
        0x8bbb0000U,
        0x48008000U,
        0x64004000U,
        0x36006000U,
        0xcb003000U,
        0x2880c800U,
        0x54402400U,
        0xfe605600U,
        0xef30fb00U,
        0x7e48e080U,
        0xaf647040U,
        0x1eb6a860U,
        0x9f8b1430U,
        0xd6c81ec8U,
        0xbb249f24U,
        0x80d6d6d6U,
        0x40bbbbbbU,
        0x60800000U,
        0x30400000U,
        0xc8600000U,
        0x24300000U,
        0x56c80000U,
        0xfb240000U,
        0xe0d60000U,
        0x70bb0000U,
        0xa8808000U,
        0x14404000U,
        0x9e606000U,
        0xdf303000U,
        0xb648c800U,
        0x8b642400U,
        0x48b65600U,
        0x648bfb00U,
        0x36486080U,
        0xcb643040U,
        0x28b6c860U,
        0x548b2430U,
        0x80000000U,
        0xc0000000U,
        0xa0000000U,
        0xd0000000U,
        0x58000000U,
        0x94000000U,
        0x3e000000U,
        0xe3000000U,
        0xbe800000U,
        0x23c00000U,
        0x1e200000U,
        0xf3100000U,
        0x46780000U,
        0x67840000U,
        0x78460000U,
        0x84670000U,
        0xc6788000U,
        0xa784c000U,
        0xd846a000U,
        0x5467d000U,
        0x9e78d800U,
        0x33845400U,
        0xe6469e00U,
        0xb7673300U,
        0x20f86680U,
        0x104477c0U,
        0xf8668020U,
        0x4477c010U,
        0x668020f8U,
        0x77c01044U,
        0x8020f866U,
        0xc0104477U,
        0xa0f86680U,
        0xd04477c0U,
        0x58668020U,
        0x9477c010U,
        0x3e8020f8U,
        0xe3c01044U,
        0xbe20f866U,
        0x23104477U,
        0x1e786680U,
        0xf38477c0U,
        0x46468020U,
        0x6767c010U,
        0x78f820f8U,
        0x84441044U,
        0xc666f866U,
        0xa7774477U,
        0xd800e680U,
        0x5400b7c0U,
        0x9e002020U,
        0x33001010U,
        0x80000000U,
        0x40000000U,
        0xa0000000U,
        0x50000000U,
        0x88000000U,
        0x24000000U,
        0x12000000U,
        0x2d000000U,
        0x76800000U,
        0x9e400000U,
        0x08200000U,
        0x64100000U,
        0xb2280000U,
        0x7d140000U,
        0xfea20000U,
        0xba490000U,
        0x1a248000U,
        0x491b4000U,
        0xc4b5a000U,
        0xe3739000U,
        0xf6800800U,
        0xde400400U,
        0xa8200a00U,
        0x34100500U,
        0x3a280880U,
        0x59140240U,
        0xeca20120U,
        0x974902d0U,
        0x6ca48768U,
        0xd75b49e4U,
        0xcc95a082U,
        0x87639641U,
        0x44a80322U,
        0xa35403d1U,
        0x568205eaU,
        0x8e590ea4U,
        0x200c8922U,
        0x100f46d1U,
        0x2817ad6bU,
        0x743a9ce7U,
        0x9a248000U,
        0x091b4000U,
        0x64b5a000U,
        0xb3739000U,
        0x7e800800U,
        0xfa400400U,
        0xba200a00U,
        0x19100500U,
        0x4ca80880U,
        0xc7540240U,
        0xe4820120U,
        0xf35902d0U,
        0x80000000U,
        0x40000000U,
        0xa0000000U,
        0x50000000U,
        0x28000000U,
        0xd4000000U,
        0x6a000000U,
        0x71000000U,
        0x38800000U,
        0x58400000U,
        0xea200000U,
        0x31100000U,
        0x98a80000U,
        0x08540000U,
        0xc22a0000U,
        0xe5250000U,
        0xf2b28000U,
        0x79484000U,
        0xfaa42000U,
        0xbd731000U,
        0x18a80800U,
        0x48540400U,
        0x622a0a00U,
        0xb5250500U,
        0xdab28280U,
        0xad484d40U,
        0x90a426a0U,
        0xcc731710U,
        0x20280b88U,
        0x10140184U,
        0x880a04a2U,
        0x84350611U,
        0x421a8b0aU,
        0xa51c4dc5U,
        0x528e2a82U,
        0x29561942U,
        0xd29a84a3U,
        0x695c4610U,
        0x72ae2b08U,
        0x39461dc6U,
        0x5ab28280U,
        0xed484d40U,
        0x30a426a0U,
        0x9c731710U,
        0x08280b88U,
        0xc4140184U,
        0xe20a04a2U,
        0xf5350611U,
        0x7a9a8b0aU,
        0xfd5c4dc5U,
        0xb8ae2a82U,
        0x18461942U,
        0x80000000U,
        0x40000000U,
        0xe0000000U,
        0xb0000000U,
        0x98000000U,
        0x94000000U,
        0x8a000000U,
        0x5b000000U,
        0x33800000U,
        0xd9c00000U,
        0x72200000U,
        0x3f100000U,
        0xc1b80000U,
        0xa6ec0000U,
        0x53860000U,
        0x29f50000U,
        0x0a3a8000U,
        0x1b2ac000U,
        0xd392e000U,
        0x69ff7000U,
        0xea380800U,
        0xab2c0400U,
        0x4ba60e00U,
        0xfde50b00U,
        0x60028980U,
        0xf006c940U,
        0x7834e8a0U,
        0x241a75b0U,
        0x123a8b38U,
        0xcf2ac99cU,
        0xb992e922U,
        0x82ff78f1U,
        0x41b80d9bU,
        0xe6ec072eU,
        0xb3860398U,
        0x99f50c2fU,
        0x923a8a1bU,
        0x8f2ac56eU,
        0x5992e2bbU,
        0x32ff70deU,
        0xd9b80980U,
        0x72ec0940U,
        0x398608a0U,
        0xc2f505b0U,
        0xa1ba8338U,
        0x56eacd9cU,
        0x2bb2e722U,
        0x0def73f1U,
        0x1800041bU,
        0xd4000e6eU,
        0x6a000b38U,
        0xeb00099fU,
        0x80000000U,
        0x40000000U,
        0xa0000000U,
        0x10000000U,
        0x08000000U,
        0x6c000000U,
        0x9e000000U,
        0x23000000U,
        0x57800000U,
        0xadc00000U,
        0x7fa00000U,
        0x91d00000U,
        0x49880000U,
        0xced40000U,
        0x880a0000U,
        0x2c0f0000U,
        0x3e0d8000U,
        0x3317c000U,
        0x5fb06000U,
        0xc1f8b000U,
        0xe18d8800U,
        0xb2d7c400U,
        0x1e106a00U,
        0x6328b100U,
        0xf7858880U,
        0xbdc3c2c0U,
        0x77ba63e0U,
        0xfdf7b330U,
        0xd7800df8U,
        0xedc0081cU,
        0xdfa0041aU,
        0x81d00a2dU,
        0x41880160U,
        0xa2d400f1U,
        0x160a069aU,
        0x0f0f09edU,
        0x698d8200U,
        0x9ed7c500U,
        0x20106a81U,
        0x5028b7c2U,
        0xa8058160U,
        0x7c03c0f1U,
        0x961a669aU,
        0x4f27b9edU,
        0xc9880a00U,
        0x8ed40100U,
        0x280a0081U,
        0x3c0f06c2U,
        0x360d89e0U,
        0x5f17c231U,
        0xc1b0657aU,
        0xe2f8baddU,
        0x80000000U,
        0x40000000U,
        0x20000000U,
        0x30000000U,
        0x58000000U,
        0xac000000U,
        0x96000000U,
        0x2b000000U,
        0xd4800000U,
        0x09400000U,
        0xe2a00000U,
        0x52500000U,
        0x4e280000U,
        0xc71c0000U,
        0x629e0000U,
        0x12670000U,
        0x6e138000U,
        0xf731c000U,
        0x3a98a000U,
        0xbe449000U,
        0xf83b8800U,
        0xdc2dc400U,
        0xee06a200U,
        0xb7239300U,
        0x1aa80d80U,
        0x8e5c0ec0U,
        0xa03e0b60U,
        0x703701b0U,
        0x783b88c8U,
        0x9c2dca54U,
        0xce06a74aU,
        0x87239795U,
        0x42a801aaU,
        0x225c08e5U,
        0x363e0a03U,
        0x5b370703U,
        0xacbb8783U,
        0x956dc9c2U,
        0x2ca6ace0U,
        0xd5739872U,
        0x0c800c2aU,
        0xe5400625U,
        0x54a00163U,
        0x495006b3U,
        0xc2a80f4bU,
        0x625c0396U,
        0x163e0baaU,
        0x6b370fe7U,
        0xf4bb8d80U,
        0x396dcec0U,
        0xbaa6ab60U,
        0xfe7391b0U,
        0x80000000U,
        0xc0000000U,
        0xa0000000U,
        0x50000000U,
        0xf8000000U,
        0x8c000000U,
        0xe2000000U,
        0x33000000U,
        0x0f800000U,
        0x21400000U,
        0x95a00000U,
        0x5e700000U,
        0xd8080000U,
        0x1c240000U,
        0xba160000U,
        0xef370000U,
        0x15868000U,
        0x9e6fc000U,
        0x781b6000U,
        0x4c349000U,
        0x420e8800U,
        0x630bcc00U,
        0xf7ad6a00U,
        0xad739500U,
        0x77800780U,
        0x6d4004c0U,
        0xd7a00420U,
        0x3d700630U,
        0x2f880f78U,
        0xb1640ad4U,
        0xcdb6077aU,
        0x824706d7U,
        0xc20e8d78U,
        0xa30bc3d6U,
        0x57ad62fbU,
        0xfd739b14U,
        0x8f8004d8U,
        0xe1400424U,
        0x35a00620U,
        0x0e700f30U,
        0x20080af8U,
        0x90240716U,
        0x581606dbU,
        0xdc370d24U,
        0x1a0683a0U,
        0xbf2fc2f0U,
        0xedbb6b5aU,
        0x12449ce7U,
        0x9a068000U,
        0x7f2fc000U,
        0x4dbb6000U,
        0x42449000U,
        0x80000000U,
        0xc0000000U,
        0x60000000U,
        0x90000000U,
        0x38000000U,
        0xc4000000U,
        0x42000000U,
        0xa3000000U,
        0xf1800000U,
        0xaa400000U,
        0xfce00000U,
        0x85100000U,
        0xe0080000U,
        0x500c0000U,
        0x58060000U,
        0x54090000U,
        0x7a038000U,
        0x670c4000U,
        0xb3842000U,
        0x094a3000U,
        0x0d6f1800U,
        0x2f5aa400U,
        0x1ce7ce00U,
        0xd5145100U,
        0xb8000080U,
        0x040000c0U,
        0x22000060U,
        0x33000090U,
        0xc9800038U,
        0x6e4000c4U,
        0xbee00042U,
        0x261000a3U,
        0x118800f1U,
        0xfa4c00aaU,
        0xa4e600fcU,
        0xd1190085U,
        0x9a0b80e0U,
        0x37004050U,
        0xeb822058U,
        0x5d433054U,
        0x776c987aU,
        0x4856e467U,
        0xaf63eeb3U,
        0xdc5e6109U,
        0xb56f188dU,
        0x2b5aa4efU,
        0x3ee7ce7cU,
        0xe6145145U,
        0x71800000U,
        0x6a400000U,
        0x9ce00000U,
        0x15100000U,
        0x80000000U,
        0x40000000U,
        0x20000000U,
        0xf0000000U,
        0xa8000000U,
        0x54000000U,
        0x9a000000U,
        0x9d000000U,
        0x1e800000U,
        0x5cc00000U,
        0x7d200000U,
        0x8d100000U,
        0x24880000U,
        0x71c40000U,
        0xeba20000U,
        0x75df0000U,
        0x6ba28000U,
        0x35d14000U,
        0x4ba3a000U,
        0xc5d2d000U,
        0xe3a16800U,
        0x91db8c00U,
        0x79aef200U,
        0x0cdf4100U,
        0x672a8080U,
        0x50154040U,
        0x1a01a020U,
        0xdd0dd0f0U,
        0x3e83e8a8U,
        0xaccacc54U,
        0xd52d529aU,
        0xd91d919dU,
        0xbe83e89eU,
        0xeccacc1cU,
        0xf52d525dU,
        0x291d917dU,
        0x1683e80cU,
        0xb8cacc65U,
        0x6f2d5251U,
        0xb41d9118U,
        0x0803e85dU,
        0xe40acc7dU,
        0x120d528cU,
        0x390d9125U,
        0x2c8be8f1U,
        0x95cecca8U,
        0xf9af5255U,
        0x4cd29199U,
        0x4729681eU,
        0xa01f8c5cU,
        0xb20cf27dU,
        0x8900418dU,
        0x80000000U,
        0xc0000000U,
        0x20000000U,
        0xd0000000U,
        0xd8000000U,
        0xc4000000U,
        0x46000000U,
        0x85000000U,
        0xa5800000U,
        0x76c00000U,
        0xada00000U,
        0x6ab00000U,
        0x2da80000U,
        0xaabc0000U,
        0x0daa0000U,
        0x7ab10000U,
        0xd5a78000U,
        0xbebd4000U,
        0x93a3e000U,
        0x3bb51000U,
        0x3629b800U,
        0x4d727c00U,
        0x9b836200U,
        0x27c4d700U,
        0xb629b880U,
        0x8d727cc0U,
        0xbb836220U,
        0xf7c4d7d0U,
        0x6e29b858U,
        0x49727c04U,
        0xfd836266U,
        0x72c4d755U,
        0xcba9b8fdU,
        0x3fb27c72U,
        0x502362cbU,
        0x1874d73fU,
        0xe601b8d0U,
        0x950e7cd8U,
        0x5d8962c6U,
        0x62c5d745U,
        0x33a63805U,
        0x2bb33c66U,
        0xce2a8255U,
        0x5970c77eU,
        0x058f8033U,
        0x66c1402bU,
        0x55a9e0ceU,
        0x7eb41059U,
        0xb3a63805U,
        0xebb33c66U,
        0xee2a8255U,
        0x8970c77eU,
    };
}

uint32_t DL::SobolUInt32(uint32_t index, uint32_t dim, uint32_t scramble)
{
    DL_ASSERT(dim < kSobolMaxDim);

    uint32_t r = scramble;
    uint32_t i = dim * kSobolSize;

    for ( ; index; index >>= 1)
    {
        if (index & 1)
            r ^= kSobolMatrices[i];

        i++;
    }

    return r ;
}

int Sobol2U::Next()
{
    int prevIndex = mIndex++;

    uint32_t flipped = mIndex ^ prevIndex;
    uint32_t is = mIBase;

    do
    {
        mU[0] ^= kSobolMatrices[is];
        mU[1] ^= kSobolMatrices[is + kSobolSize];

        flipped = flipped >> 1;
        is++;
    }
    while (flipped);

    return mIndex;
}

void Sobol2U::Set(int index, int scramble)
{
    mU[0] = scramble;
    mU[1] = scramble;
    mIndex = index;

    uint32_t is = mIBase;

    for ( ; index; index >>= 1, is++ )
    {
        if (index & 1)
        {
            mU[0] ^= kSobolMatrices[is];
            mU[1] ^= kSobolMatrices[is + kSobolSize];
        }
    }
}

int Sobol3U::Next()
{
    int prevIndex = mIndex++;

    uint32_t flipped = mIndex ^ prevIndex;
    uint32_t is = mIBase;

    do
    {
        mU[0] ^= kSobolMatrices[is];
        mU[1] ^= kSobolMatrices[is + kSobolSize];
        mU[2] ^= kSobolMatrices[is + kSobolSize * 2];

        flipped = flipped >> 1;
        is++;
    }
    while (flipped);

    return mIndex;
}

void Sobol3U::Set(int index, int scramble)
{
    mU[0] = scramble;
    mU[1] = scramble;
    mU[2] = scramble;
    mIndex = index;

    uint32_t is = mIBase;

    for ( ; index; index >>= 1, is++ )
    {
        if (index & 1)
        {
            mU[0] ^= kSobolMatrices[is];
            mU[1] ^= kSobolMatrices[is + kSobolSize];
            mU[2] ^= kSobolMatrices[is + kSobolSize * 2];
        }
    }
}

// JSF (Bob Jenkins)

JSF32::JSF32(uint32_t seed) : a(0xf1ea5eed), b(seed), c(seed), d(seed)
{
    for (int i = 0; i < 20; i++)
        (void) Next();
}
