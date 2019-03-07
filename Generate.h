//
//  File:       Generate.h
//
//  Function:   Various sample/number generators
//
//  Copyright:  Andrew Willmott, 2018
//

#ifndef GENERATE_H
#define GENERATE_H

#include <math.h>
#include <stdint.h>

namespace DistLib
{
    //--------------------------------------------------------------------------
    // LCG
    //--------------------------------------------------------------------------

    struct cLCG
    {
        uint32_t mState = 0x12345678;

        cLCG() {};
        cLCG(uint32_t seed) : mState(seed) {}

        uint32_t Next();        ///< Explicit next number in sequence
        operator uint32_t();    ///< Return next number in sequence in uint32_t context
    };


    //--------------------------------------------------------------------------
    // PCG
    //--------------------------------------------------------------------------

    struct cPCG
    /// See http://www.pcg-random.org
    {
        uint64_t mState = 0x853c49e6748fea9bULL;
        uint64_t mInc   = 0xda3e39cb94b95bdbULL;

        cPCG() {}
        cPCG(uint64_t initstate, uint64_t initseq);

        uint32_t Next();        ///< Explicit next number in sequence
        operator uint32_t();    ///< Return next number in sequence in uint32_t context
    };


    //--------------------------------------------------------------------------
    // XorShift
    //--------------------------------------------------------------------------

    struct cXORShift
    /// See https://en.wikipedia.org/wiki/Xorshift
    {
        uint32_t mState = 0x12345678;

        cXORShift() {}
        cXORShift(uint32_t seed);

        uint32_t Next();        ///< Explicit next number in sequence
        operator uint32_t();    ///< Return next number in sequence in uint32_t context
    };


    //--------------------------------------------------------------------------
    // cHashGen
    //--------------------------------------------------------------------------

    struct cHashGen
    /// Works by hashing a counter, which means it's random access, which may
    /// be necessary in some applications. Applies an LCG and then xorshift to
    /// the counter to avoid patterns.
    {
        uint32_t mIndex = 0;

        uint32_t Next();        ///< Advance to next point in the sequence. Returns the index of this point.
        void     Set (int n);   ///< Jump directly to term 'n' of the sequence
        operator uint32_t();    ///< Return next number in sequence in uint32_t context
    };
    uint32_t HashInt(uint32_t i);   ///< Hash used by above, in case direct use is more convenient


    //--------------------------------------------------------------------------
    // Halton
    //--------------------------------------------------------------------------

    float    HaltonFloat (int i, int b); ///< return term i of the base b Halton sequence
    uint32_t HaltonUInt32(int i, int b); ///< return term i of the base b Halton sequence as fraction of UINT32_MAX

    float    Halton2Float (int i);   ///< return term i of the base 2 Halton sequence
    uint32_t Halton2UInt32(int i);   ///< return term i of the base 2 Halton sequence as fraction of UINT32_MAX

    struct cHalton2
    /// This calculates the 2D Halton sequence incrementally, faster per-call than HaltonFloat(i, 2/3).
    {
        float    mU[2] = { 0.0f, 0.0f };

        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;
        
        int  Next ();       ///< Advance to next point in the sequence. Returns the index of this point.
        void Set  (int n);  ///< Jump directly to term 'n' of the sequence
    };
    
    struct cHalton3
    /// This calculates the 3D Halton sequence incrementally, faster per-call than HaltonFloat(i, 2/3/5).
    {
        float    mU[3] = { 0.0f, 0.0f, 0.0f };

        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;
        uint32_t mBase5 = 0;
        
        int  Next ();       ///< Advance to next point in the sequence. Returns the index of this point.
        void Set  (int n);  ///< Jump directly to term 'n' of the sequence
    };

    struct cHalton2U
    /// This calculates the 2D Halton sequence incrementally, faster per call than HaltonUInt32.
    {
        uint32_t mU[2] = { 0 };     ///< Current sample point

        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;
        
        int  Next ();       ///< Advance to next point in the sequence. Returns the index of this point.
        void Set  (int n);  ///< Jump directly to term 'n' of the sequence
    };

    struct cHalton3U
    /// This calculates the 3D Halton sequence incrementally, faster per-call than HaltonUInt32(i, 2/3/5).
    {
        uint32_t mU[3] = { 0 };     ///< Current sample point
        
        uint32_t mBase2 = 0;
        uint32_t mBase3 = 0;
        uint32_t mBase5 = 0;
        
        int  Next ();       ///< Advance to next point in the sequence. Returns the index of this point.
        void Set  (int n);  ///< Jump directly to term 'n' of the sequence
    };

    
    //--------------------------------------------------------------------------
    // Fibonacci/Golden Spiral
    //--------------------------------------------------------------------------
    
    float    GoldenFloat (int i); ///< return term i of the golden angle sequence as float
    uint32_t GoldenUInt32(int i); ///< return term i of the golden angle sequence as fraction of UINT32_MAX

    struct cGolden2U
    {
        uint32_t mU[2] = { 0 };     ///< Current sample point
        uint32_t mStep;

        cGolden2U(int numSamples);
        
        void Next();        ///< Advance to next point in the sequence.
        void Set(int n);    ///< Jump directly to term 'n' of the sequence
    };


    //--------------------------------------------------------------------------
    // Generalised Kronecker/golden ratio combination, Rd.
    // See http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    //--------------------------------------------------------------------------

    struct cR1U
    {
        uint32_t mU = UINT32_MAX / 2; ///< Current sample point

        void Next();        ///< Advance to next point in the sequence.
        void Set(int n);    ///< Jump directly to term 'n' of the sequence

        operator uint32_t();    ///< Return next number in sequence in uint32_t context
    };

    struct cR2U
    {
        uint32_t mU[2] = { UINT32_MAX / 2, UINT32_MAX / 2 }; ///< Current sample point

        void Next();        ///< Advance to next point in the sequence.
        void Set(int n);    ///< Jump directly to term 'n' of the sequence
    };

    struct cR3U
    {
        uint32_t mU[3] = { UINT32_MAX / 2, UINT32_MAX / 2, UINT32_MAX / 2 }; ///< Current sample point

        void Next();        ///< Advance to next point in the sequence.
        void Set(int n);    ///< Jump directly to term 'n' of the sequence
    };


    //--------------------------------------------------------------------------
    // Jittered version of R2.
    // See http://extremelearning.com.au/a-simple-method-to-construct-isotropic-quasirandom-blue-noise-point-sequences/
    //--------------------------------------------------------------------------

    struct cR2JitterU
    {
        uint32_t mU[2]  = { 0 };     ///< Current sample point

        cR2U     mR;        ///< R2 sequence
        cHashGen mJ;        ///< Jitter
        float    mS;        ///< Scale term

        cR2JitterU(float lambda = 1.0f);    ///< Lambda controls jitter amount

        void Next();        ///< Advance to next point in sequence
        void Set(int n);    ///< Jump directly to term 'n' of the sequence
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

    inline uint32_t cLCG::Next()
    {
        uint32_t oldState = mState;
        mState = uint32_t(mState * uint64_t(1103515245) + 12345);
        return oldState;
    }

    inline cLCG::operator uint32_t()
    {
        return Next();
    }

    inline cPCG::cPCG(uint64_t initstate, uint64_t initseq) : mState(0), mInc((initseq << 1u) | 1u)
    {
        Next();
        mState += initstate;
        Next();
    }

    inline uint32_t cPCG::Next()
    {
        uint64_t oldState = mState;
        mState = oldState * 6364136223846793005ULL + mInc;

        uint32_t xorShifted = uint32_t(((oldState >> 18u) ^ oldState) >> 27u);
        uint32_t rot = oldState >> 59u;

        return (xorShifted >> rot) | (xorShifted << ((-int32_t(rot)) & 31)); // int32_t cast added as latest VS treats -u as error by default \o/
    }

    inline cPCG::operator uint32_t()
    {
        return Next();
    }

    inline cXORShift::cXORShift(uint32_t seed) : mState(seed)
    {}

    inline uint32_t cXORShift::Next()
    {
        mState ^= (mState << 13);
        mState ^= (mState >> 17);
        mState ^= (mState << 5);
        return mState;
    }

    inline cXORShift::operator uint32_t()
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

    inline void cHashGen::Set(int i)
    {
        mIndex = i;
    }

    inline uint32_t cHashGen::Next()
    {
        return HashInt(mIndex++);
    }

    inline cHashGen::operator uint32_t()
    {
        return Next();
    }

    const float    kGoldenF32 = 0.5f * (sqrtf(5.0f) - 1.0f);
    const uint32_t kGoldenU32 = uint32_t(UINT32_MAX * kGoldenF32);

    inline float GoldenFloat(int i)
    {
        float f = i * kGoldenF32;
        return f - floorf(f);
    }

    inline uint32_t GoldenUInt32(int i)
    {
        return uint32_t(i * kGoldenU32);
    }


    inline cGolden2U::cGolden2U(int numSamples) :
        mStep(UINT32_MAX / numSamples)
    {
        mU[1] = mStep / 2;
    }

    inline void cGolden2U::Next()
    {
        mU[0] += kGoldenU32;
        mU[1] += mStep;
    }

    inline void cGolden2U::Set(int i)
    {
        mU[0] = uint32_t(i * kGoldenU32);
        mU[1] = i * mStep + mStep / 2;
    }


    constexpr float    kG1 = 1.6180339887498948482f;
    constexpr float    kR1xF32 = 1.0f / kG1;
    constexpr uint32_t kR1xU32 = uint32_t(UINT32_MAX * kR1xF32);    // Same as kGoldenU32, the Rn sequences are a generalisation of this

    inline void cR1U::Next()
    {
        mU += kR1xU32;
    }

    inline void cR1U::Set(int i)
    {
        mU = uint32_t(UINT32_MAX / 2 + i * kR1xU32);
    }

    inline cR1U::operator uint32_t()
    {
        uint32_t u = mU;
        mU += kR1xU32;
        return u;
    }

    constexpr float    kG2 = 1.32471795724474602596f;
    constexpr float    kR2xF32 = 1.0f / kG2;
    constexpr float    kR2yF32 = 1.0f / (kG2 * kG2);
    constexpr uint32_t kR2xU32 = uint32_t(UINT32_MAX * kR2xF32);
    constexpr uint32_t kR2yU32 = uint32_t(UINT32_MAX * kR2yF32);

    inline void cR2U::Next()
    {
        mU[0] += kR2xU32;
        mU[1] += kR2yU32;
    }

    inline void cR2U::Set(int i)
    {
        mU[0] = uint32_t(UINT32_MAX / 2 + i * kR2xU32);
        mU[1] = uint32_t(UINT32_MAX / 2 + i * kR2yU32);
    }

    constexpr float    kG3 = 1.22074408460575947536f;
    constexpr float    kR3xF32 = 1.0f / kG3;
    constexpr float    kR3yF32 = 1.0f / (kG3 * kG3);
    constexpr float    kR3zF32 = 1.0f / (kG3 * kG3 * kG3);
    constexpr uint32_t kR3xU32 = uint32_t(UINT32_MAX * kR3xF32);
    constexpr uint32_t kR3yU32 = uint32_t(UINT32_MAX * kR3yF32);
    constexpr uint32_t kR3zU32 = uint32_t(UINT32_MAX * kR3zF32);

    inline void cR3U::Next()
    {
        mU[0] += kR3xU32;
        mU[1] += kR3yU32;
        mU[2] += kR3zU32;
    }

    inline void cR3U::Set(int i)
    {
        mU[0] = uint32_t(UINT32_MAX / 2 + i * kR3xU32);
        mU[1] = uint32_t(UINT32_MAX / 2 + i * kR3yU32);
        mU[2] = uint32_t(UINT32_MAX / 2 + i * kR3zU32);
    }

    constexpr float    kR2dF32 = 0.76f * 1.772453851f / 4.0f;  // 1.77245 is sqrt(pi), as C++11 didn't make sqrt constexpr \o/
    constexpr uint32_t kR2dU32 = uint32_t(UINT32_MAX * kR2dF32);

    inline cR2JitterU::cR2JitterU(float lambda) :
        mS(kR2dU32 * lambda)
    {
        uint64_t si = uint64_t(mS / sqrtf(0.3f));

        mU[0] = mR.mU[0] + ((si * mJ.Next()) >> 32);
        mU[1] = mR.mU[1] + ((si * mJ.Next()) >> 32);
    }

    inline void cR2JitterU::Next()
    {
        mR.Next();
        mJ.mIndex += 2;

        uint64_t si = uint64_t(mS / sqrtf(mJ.mIndex / 2 + 0.3f));

        mU[0] = mR.mU[0] + ((si * mJ.Next()) >> 32);
        mU[1] = mR.mU[1] + ((si * mJ.Next()) >> 32);
    }

    inline void cR2JitterU::Set(int n)
    {
        mR.Set(n);
        mJ.Set(2 * n);

        uint64_t si = uint64_t(mS / sqrtf(n + 0.3f));

        mU[0] = mR.mU[0] + ((si * mJ.Next()) >> 32);
        mU[1] = mR.mU[1] + ((si * mJ.Next()) >> 32);
    }
}

#endif
