#include "Distribute.h"

#include <stdio.h>

using namespace DistLib;

inline uint32_t LCG(uint32_t* state)
{
    uint32_t current(*state);
    *state = uint32_t(*state * uint64_t(1103515245) + 12345);
    return current;
}

int main(int argc, const char* argv[])
{
    uint32_t state = 0x12456789;

    for (int i = 0; i < 100; i++)
    {
        uint32_t u = LCG(&state);

        printf("u: 0x%08x\n", u); 
        
        printf("  ToInt32(u, 10) = %d\n", ToInt32(u, 10)); 
        printf("  ToFloat(u)   = %g\n", ToFloat(u)); 
        //printf("  ToIntWeighted(u, 10) = %d\n", ToInt32Weighted(u, 10));
        
        Vec2f v2 = ToVec2fCircle(u); 
        printf("  ToVec2fCircle(u)   = %g, %g\n", v2.x, v2.y); 
    }
    
    return 0;
}
