//
//  File:       Generate.cpp
//
//  Function:   Various sample/number generators
//
//  Copyright:  Andrew Willmott, 2018
//

#include "Generate.h"

using namespace DistLib;

namespace
{
    constexpr float kOneOverThree = float(1.0 / 3.0);
    constexpr float kOneOverFive  = float(1.0 / 5.0);
}

float DistLib::HaltonFloat(int n, int b)
/// return term i of the base b Halton sequence
/// In fact, this is just a generalization of Heckbert's bit reversal distribution trick.
/// E.g., when b=3, write n as a base 3 number, digit 0 -> which third of interval the
/// sample is in, 1 -> which third of that, 2 -> which third of that, etc.
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

uint32_t DistLib::HaltonUInt32(int n, int b)
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


// --- cHalton2 --------------------------------------------------------

int cHalton2::Next()
{
    /////////////////////////////////////
    // base 2
    
    uint32_t oldBase2 = mBase2;
    mBase2++;
    uint32_t diff = mBase2 ^ oldBase2;

    // bottom bit always changes, higher bits
    // change less frequently.
    float s = 0.5;

    // diff will be of the form 0*1+, i.e. one bits up until the last carry.
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

void cHalton2::Set(int n)
{
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


// --- cHalton3 --------------------------------------------------------

int cHalton3::Next()
{
    // base 2: 1 bit per digit
    uint32_t oldBase2 = mBase2;
    mBase2++;
    uint32_t diff = mBase2 ^ oldBase2;

    // bottom bit always changes, higher bits
    // change less frequently.
    float s = 0.5;

    // diff will be of the form 0*1+, i.e. one bits up until the last carry.
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

void cHalton3::Set(int n)
{
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

// --- cHalton2U --------------------------------------------------------

int cHalton2U::Next()
{
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

void cHalton2U::Set(int n)
{
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


// --- cHalton3U --------------------------------------------------------

int cHalton3U::Next()
{
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

void cHalton3U::Set(int n)
{
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

