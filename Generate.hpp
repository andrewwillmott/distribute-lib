//
// Generate.hpp
//
// Various sample/number generators
//
// Andrew Willmott
//

#ifndef GENERATE_H
#define GENERATE_H

#include <math.h>
#include <stdint.h>

#ifndef DL_ASSERT
    #define DL_ASSERT(X)
#endif

namespace DL
{
    //--------------------------------------------------------------------------
    // LCG
    //--------------------------------------------------------------------------

    struct LCG
    {
        uint32_t mState = 0x12345678;

        LCG() {};
        LCG(uint32_t seed) : mState(seed) {}

        uint32_t Next();        // Explicit next number in sequence
        operator uint32_t();    // Return next number in sequence in uint32_t context
    };


    //--------------------------------------------------------------------------
    // PCG
    //--------------------------------------------------------------------------

    constexpr uint64_t kPCGDefaultState  = 0x853c49e6748fea9bULL;
    constexpr uint64_t kPCGDefaultStream = 0xda3e39cb94b95bdbULL;

    struct PCG
    // See http://www.pcg-random.org
    {
        uint64_t mState = kPCGDefaultState;
        uint64_t mInc   = kPCGDefaultStream;

        PCG() {}
        PCG(uint64_t seed, uint64_t stream = kPCGDefaultStream >> 1);

        uint32_t Next();        // Explicit next number in sequence
        operator uint32_t();    // Return next number in sequence in uint32_t context

        void Advance(int delta);  // Advance state by 'delta' steps.
    };

    PCG PCGStream(uint32_t seed, uint32_t streamID, int index = 0);
    PCG PCGStream(uint64_t seed, uint32_t streamID, int index = 0);


    //--------------------------------------------------------------------------
    // XorShift
    //--------------------------------------------------------------------------

    struct XORShift
    // See https://en.wikipedia.org/wiki/Xorshift
    {
        uint32_t mState = 0x12345678;

        XORShift() {}
        XORShift(uint32_t seed);

        uint32_t Next();        // Explicit next number in sequence
        operator uint32_t();    // Return next number in sequence in uint32_t context
    };


    //--------------------------------------------------------------------------
    // HashGen
    //--------------------------------------------------------------------------

    struct HashGen
    // Works by hashing a counter, which means it's random access, which may
    // be necessary in some applications. Applies an LCG and then xorshift to
    // the counter to avoid patterns.
    {
        uint32_t mIndex = 0;

        uint32_t Next();        // Advance to next point in the sequence. Returns the index of this point.
        void     Set (int n);   // Jump directly to term 'n' of the sequence
        operator uint32_t();    // Return next number in sequence in uint32_t context
    };
    uint32_t HashInt(uint32_t i);   // Hash used by above, in case direct use is more convenient


    //--------------------------------------------------------------------------
    // Halton
    //--------------------------------------------------------------------------

    float    HaltonFloat (int i, int b);  // return term i of the base b Halton sequence
    uint32_t HaltonUInt32(int i, int b);  // return term i of the base b Halton sequence as fraction of UINT32_MAX
    double   HaltonDouble(int i, int b);  // return term i of the base b Halton sequence

    float    Halton2Float (int i);   // return term i of the base 2 Halton sequence
    uint32_t Halton2UInt32(int i);   // return term i of the base 2 Halton sequence as fraction of UINT32_MAX
    double   Halton2Double(int i);   // return term i of the base 2 Halton sequence

    struct Halton2
    // This calculates the 2D Halton sequence incrementally, faster per-call than HaltonFloat(i, 2/3).
    {
        float    mU[2] = { 0.0f, 0.0f };

        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;

        int  Next();       // Advance to next point in the sequence. Returns the index of this point.
        void Set (int n);  // Jump directly to term 'n' of the sequence
    };

    struct Halton3
    // This calculates the 3D Halton sequence incrementally, faster per-call than HaltonFloat(i, 2/3/5).
    {
        float    mU[3] = { 0.0f, 0.0f, 0.0f };

        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;    // holds 16 base 3 digits, so max 'i' is 43,046,721
        uint32_t mBase5 = 0;    // holds 10 base 5 digits, so max 'i' is 9,765,625

        int  Next();       // Advance to next point in the sequence. Returns the index of this point.
        void Set (int n);  // Jump directly to term 'n' of the sequence
    };

    struct Halton2D
    // Double-precision variant of Halton2
    {
        double   mU[2] = { 0.0, 0.0 };

        uint64_t mBase2 = 0;
        uint64_t mBase3 = 0;

        uint64_t Next();       // Advance to next point in the sequence. Returns the index of this point.
        void     Set (int n);  // Jump directly to term 'n' of the sequence
    };

    struct Halton3D
    // Double-precision variant of Halton3
    {
        double   mU[3] = { 0.0, 0.0, 0.0 };

        uint64_t mBase2 = 0;
        uint64_t mBase3 = 0;
        uint64_t mBase5 = 0;

        uint64_t Next();       // Advance to next point in the sequence. Returns the index of this point.
        void     Set (int n);  // Jump directly to term 'n' of the sequence
    };

    struct Halton2U
    // U32 variant of Halton2 for use with Distribute functions.
    {
        uint32_t mU[2] = { 0 };     // Current sample point

        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;

        int  Next();       // Advance to next point in the sequence. Returns the index of this point.
        void Set (int n);  // Jump directly to term 'n' of the sequence
    };

    struct Halton3U
    // U32 variant of Halton3 for use with Distribute functions.
    {
        uint32_t mU[3] = { 0 };    // Current sample point

        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;
        uint32_t mBase5 = 0;

        int  Next();       // Advance to next point in the sequence. Returns the index of this point.
        void Set (int n);  // Jump directly to term 'n' of the sequence
    };


    //--------------------------------------------------------------------------
    // 2D Sobol
    //--------------------------------------------------------------------------

    float    SobolFloat (uint32_t index, uint32_t dim = 0, uint32_t scramble = 0);  // return term 'index' for dimension 'dim' of the Sobol sequence as 0-1 float.
    uint32_t SobolUInt32(uint32_t index, uint32_t dim = 0, uint32_t scramble = 0);  // return term 'index' for dimension 'dim' of the Sobol sequence as uint32 suitable for use with distribute calls.

    struct Sobol2U
    // Calculates Sobol sequence incrementally, faster than standalone calls.
    {
        uint32_t mU[2] = { 0 };

        uint32_t mIBase;
        uint32_t mIndex;

        Sobol2U(int startDim = 0, int scramble = 0);

        int  Next();       // Advance to next point in the sequence. Returns the index of this point.
        void Set (int n, int scramble = 0);  // Jump directly to term 'n' of the sequence
    };

    struct Sobol3U
    // Calculates Sobol sequence incrementally, faster than standalone calls.
    {
        uint32_t mU[3] = { 0 };

        uint32_t mIBase;
        uint32_t mIndex;

        Sobol3U(int startDim = 0, int scramble = 0);

        int  Next();       // Advance to next point in the sequence. Returns the index of this point.
        void Set (int n, int scramble = 0);  // Jump directly to term 'n' of the sequence
    };


    //--------------------------------------------------------------------------
    // Fibonacci/Golden Spiral
    //--------------------------------------------------------------------------

    float    GoldenFloat (int i); // return term i of the golden angle sequence as float
    uint32_t GoldenUInt32(int i); // return term i of the golden angle sequence as fraction of UINT32_MAX

    struct Golden2U
    {
        uint32_t mU[2] = { 0 };     // Current sample point
        uint32_t mStep;

        Golden2U(int numSamples);

        void Next();        // Advance to next point in the sequence.
        void Set(int n);    // Jump directly to term 'n' of the sequence
    };


    //--------------------------------------------------------------------------
    // Generalised Kronecker/golden ratio combination, Rd.
    // See http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    //--------------------------------------------------------------------------

    struct R1U
    {
        uint32_t mU[1] = { UINT32_MAX / 2 }; // Current sample point

        uint32_t Next();      // Advance to next point in the sequence.
        void     Set(int n);  // Jump directly to term 'n' of the sequence

        operator uint32_t();  // Return next number in sequence in uint32_t context
    };

    struct R2U
    {
        uint32_t mU[2] = { UINT32_MAX / 2, UINT32_MAX / 2 }; // Current sample point

        void Next();        // Advance to next point in the sequence.
        void Set(int n);    // Jump directly to term 'n' of the sequence
    };

    struct R3U  // Note: this is not nearly as good as R2U, very clumpy with obvious striding.
    {
        uint32_t mU[3] = { UINT32_MAX / 2, UINT32_MAX / 2, UINT32_MAX / 2 }; // Current sample point

        void Next();        // Advance to next point in the sequence.
        void Set(int n);    // Jump directly to term 'n' of the sequence
    };


    //--------------------------------------------------------------------------
    // JSF (Bob Jenkins)
    //--------------------------------------------------------------------------

    struct JSF32
    {
        uint32_t a;  // State
        uint32_t b;
        uint32_t c;
        uint32_t d;

        JSF32(uint32_t seed = 0x12345678);
        uint32_t Next();    // Advance to next point in the sequence.

        operator uint32_t();
    };


    // --- Inlines -------------------------------------------------------------

    inline float Halton2Float(int i)
    {
        uint32_t u = Halton2UInt32(i);

        uint32_t resultU(0x3f800000 | (u >> 9));
        float result = ((float&) resultU);
        return result - 1.0f;
    }

    inline uint32_t Halton2UInt32(int a)
    {
        uint32_t b;

        b = ((a & 0x55555555) << 1)  | ((a & 0xAAAAAAAA) >> 1);
        a = ((b & 0x33333333) << 2)  | ((b & 0xCCCCCCCC) >> 2);
        b = ((a & 0x0F0F0F0F) << 4)  | ((a & 0xF0F0F0F0) >> 4);
        a = ((b & 0x00FF00FF) << 8)  | ((b & 0xFF00FF00) >> 8);
        b = ((a & 0x0000FFFF) << 16) | ((a & 0xFFFF0000) >> 16);

        return b;
    }

    inline uint32_t LCG::Next()
    {
        uint32_t oldState = mState;
        mState = mState * 1103515245 + 12345;
        return oldState;
    }

    inline LCG::operator uint32_t()
    {
        return Next();
    }

    inline PCG::PCG(uint64_t state, uint64_t stream) : mState(0), mInc((stream << 1u) | 1u)
    {
        Next();
        mState += state;
        Next();
    }

    #define HL_PCG_MULTIPLIER_64  6364136223846793005ULL

    inline uint32_t PCG::Next()
    {
        uint64_t oldState = mState;
        mState = oldState * HL_PCG_MULTIPLIER_64 + mInc;

        uint32_t xorShifted = uint32_t(((oldState >> 18u) ^ oldState) >> 27u);
        uint32_t rot = oldState >> 59u;

        return (xorShifted >> rot) | (xorShifted << ((-int32_t(rot)) & 31)); // int32_t cast added as latest VS treats -u as error by default \o/
    }

    inline PCG::operator uint32_t()
    {
        return Next();
    }

    inline PCG PCGStream(uint32_t seed, uint32_t streamID, int index)
    {
      return PCG((uint64_t(seed) << 32) | seed, (uint64_t(streamID) << 31) | index);
    }

    inline PCG PCGStream(uint64_t seed, uint32_t streamID, int index)
    {
      return PCG(seed, (uint64_t(streamID) << 31) | index);
    }

    inline XORShift::XORShift(uint32_t seed) : mState(seed)
    {}

    inline uint32_t XORShift::Next()
    {
        mState ^= (mState << 13);
        mState ^= (mState >> 17);
        mState ^= (mState << 5);
        return mState;
    }

    inline XORShift::operator uint32_t()
    {
        return Next();
    }

    inline uint32_t HashInt(uint32_t i)
    {
        uint32_t hash = i * 1103515245 + 12345;
        hash ^= (hash << 13);
        hash ^= (hash >> 17);
        hash ^= (hash << 5);
        return hash;
    }

    inline void HashGen::Set(int i)
    {
        mIndex = i;
    }

    inline uint32_t HashGen::Next()
    {
        return HashInt(mIndex++);
    }

    inline HashGen::operator uint32_t()
    {
        return Next();
    }

    const uint32_t kSobolSize = 52;
    const uint32_t kSobolMaxDim = 16;

    inline float SobolFloat(uint32_t index, uint32_t dim, uint32_t scramble)
    {
        return float(1.0 / 0xFFFFFFFF) * SobolUInt32(index, dim, scramble);
    }

    inline Sobol2U::Sobol2U(int startDim, int scramble) :
        mIBase(startDim * kSobolSize),
        mIndex(0)
    {
        DL_ASSERT(startDim + 1 < kSobolMaxDim);
        mU[0] = scramble;
        mU[1] = scramble;
    }

    inline Sobol3U::Sobol3U(int startDim, int scramble) :
        mIBase(startDim * kSobolSize),
        mIndex(0)
    {
        DL_ASSERT(startDim + 2 < kSobolMaxDim);
        mU[0] = scramble;
        mU[1] = scramble;
        mU[2] = scramble;
    }

    constexpr float    kGoldenF32 = 0.6180339887498949f;   // = 0.5f * (sqrtf(5.0f) - 1.0f). sqrt _still_ isn't constexpr
    constexpr uint32_t kGoldenU32 = uint32_t(double(UINT32_MAX) * kGoldenF32);

    inline float GoldenFloat(int i)
    {
        float f = i * kGoldenF32;
        return f - floorf(f);
    }

    inline uint32_t GoldenUInt32(int i)
    {
        return uint32_t(i * kGoldenU32);
    }


    inline Golden2U::Golden2U(int numSamples) :
        mStep(UINT32_MAX / numSamples)
    {
        mU[1] = mStep / 2;
    }

    inline void Golden2U::Next()
    {
        mU[0] += kGoldenU32;
        mU[1] += mStep;
    }

    inline void Golden2U::Set(int i)
    {
        mU[0] = uint32_t(i * kGoldenU32);
        mU[1] = i * mStep + mStep / 2;
    }


    constexpr double   kG1 = 1.6180339887498948482;
    constexpr double   kR1xF32 = 1.0 / kG1;
    constexpr uint32_t kR1xU32 = uint32_t(UINT32_MAX * kR1xF32);    // Same as kGoldenU32, the Rn sequences are a generalisation of this

    inline uint32_t R1U::Next()
    {
        mU[0] += kR1xU32;
        return mU[0];
    }

    inline void R1U::Set(int i)
    {
        mU[0] = uint32_t(UINT32_MAX / 2 + i * kR1xU32);
    }

    inline R1U::operator uint32_t()
    {
        uint32_t u = mU[0];
        mU[0] += kR1xU32;
        return u;
    }

    constexpr double   kG2 = 1.32471795724474602596;
    constexpr double   kR2xF32 = 1.0f / kG2;
    constexpr double   kR2yF32 = 1.0f / (kG2 * kG2);
    constexpr uint32_t kR2xU32 = uint32_t(UINT32_MAX * kR2xF32);
    constexpr uint32_t kR2yU32 = uint32_t(UINT32_MAX * kR2yF32);

    inline void R2U::Next()
    {
        mU[0] += kR2xU32;
        mU[1] += kR2yU32;
    }

    inline void R2U::Set(int i)
    {
        mU[0] = uint32_t(UINT32_MAX / 2 + i * kR2xU32);
        mU[1] = uint32_t(UINT32_MAX / 2 + i * kR2yU32);
    }

    constexpr double kG3 = 1.22074408460575947536f;
    constexpr double kR3xF32 = 1.0 / kG3;
    constexpr double kR3yF32 = 1.0 / (kG3 * kG3);
    constexpr double kR3zF32 = 1.0 / (kG3 * kG3 * kG3);
    constexpr uint32_t kR3xU32 = uint32_t(UINT32_MAX * kR3xF32);
    constexpr uint32_t kR3yU32 = uint32_t(UINT32_MAX * kR3yF32);
    constexpr uint32_t kR3zU32 = uint32_t(UINT32_MAX * kR3zF32);

    inline void R3U::Next()
    {
        mU[0] += kR3xU32;
        mU[1] += kR3yU32;
        mU[2] += kR3zU32;
    }

    inline void R3U::Set(int i)
    {
        mU[0] = uint32_t(UINT32_MAX / 2 + i * kR3xU32);
        mU[1] = uint32_t(UINT32_MAX / 2 + i * kR3yU32);
        mU[2] = uint32_t(UINT32_MAX / 2 + i * kR3zU32);
    }

    inline uint32_t JSF32::Next()
    {
        uint32_t e;

        #define rot32(x, k) (((x) << (k)) | ((x) >> (32 - (k))))
        e = a - rot32(b, 27);
        a = b ^ rot32(c, 17);
        b = c + d;
        c = d + e;
        d = e + a;
        #undef rot32

        return d;
    }

    inline JSF32::operator uint32_t()
    {
        return Next();
    }
}

#endif
