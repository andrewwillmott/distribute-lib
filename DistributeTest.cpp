#include "Distribute.h"

#include <stdio.h>
#include <stdlib.h>

using namespace DistLib;

namespace
{
    inline uint32_t LCG(uint32_t* state)
    {
        uint32_t current(*state);
        *state = uint32_t(*state * uint64_t(1103515245) + 12345);
        return current;
    }
    
    const int kWeights[] = { 2, 1, 0, 1 };
}

int main(int argc, const char* argv[])
{
    int numSamples = 100;
    
    if (argc > 1)
        numSamples = atoi(argv[1]);
    
    uint32_t state = 0x12456789;

    for (int i = 0; i < numSamples; i++)
    {
        uint32_t u = LCG(&state);

        printf("-- u: 0x%08x ------------------------------------------------\n", u);
        
        printf("  ToInt32(u, 10)    = %2d\n", ToInt32(u, 10));
        printf("  ToFloat(u, 100)   = %5.2f\n", ToFloat(u, 100.0f));

        Vec2f v2 = ToCircle(u);
        printf("  ToCircle(u)       = %5.2f, %5.2f\n", v2.x, v2.y);

        Vec3f v3 = ToSphere(u);
        printf("  ToSphere(u)       = %5.2f, %5.2f, %5.2f\n", v3.x, v3.y, v3.z);

        uint32_t u0 = u;

        printf("  u = ModTriangle(u):\n");
        u = ModTriangle(u0);
        printf("    ToInt32(u, 10)  = %2d\n", ToInt32(u, 10));
        printf("    ToFloat(u, 100) = %5.2f\n", ToFloat(u, 100.0f));

        printf("  u = ModWeighted(u, [2, 1, 0, 1]):\n");
        u = ModWeighted(u0, sizeof(kWeights) / sizeof(kWeights[0]), kWeights);
        printf("    ToInt32(u, 10)  = %2d\n", ToInt32(u, 10));
        printf("    ToFloat(u, 100) = %5.2f\n", ToFloat(u, 100.0f));
    }
    
    return 0;
}
