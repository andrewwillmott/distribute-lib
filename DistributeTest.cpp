#include "Generate.hpp"
#include "Distribute.hpp"

#include <stdio.h>
#include <stdlib.h>

using namespace DL;

#ifdef _MSC_VER
    #pragma warning (disable: 4996)     // let's not make fopen an error.
#endif

namespace
{
    const int kWeights[] = { 2, 1, 0, 1 };

    // Mini SVG api
    void svg_header(FILE* out, int size)
    {
        fprintf(out, "<svg xmlns=\"http://www.w3.org/2000/svg\" height=\"%d\" width=\"%d\">\n", size, size);
    }

    void svg_trailer(FILE* out)
    {
        fprintf(out, "</svg>\n");
    }

    void svg_circle(FILE* out, Vec2f p0, float r, const char* colour, int outline = 0)
    {
        fprintf(out, "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"black\" stroke-width=\"%d\" fill=\"%s\" />\n",
            p0.x, p0.y, r,
            outline, colour
        );
    }

    void svg_line(FILE* out, Vec2f p0, Vec2f p1)
    {
        fprintf(out, "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" style=\"stroke:rgb(0,0,0);stroke-width:2\" />\n",
            p0.x, p0.y, p1.x, p1.y
        );
    }
    void svg_text(FILE* out, Vec2f p0, const char* text, const char* colour)
    {
        fprintf(out, "<text x=\"%f\" y=\"%f\" fill=\"%s\">%s</text>\n", p0.x, p0.y, colour, text);
    }
}

int main(int argc, const char* argv[])
{
    int numSamples = 1000;
    
    if (argc > 1)
        numSamples = atoi(argv[1]);

    uint32_t seed = 0x12345678;

    if (argc > 2)
        seed = atoi(argv[2]);

    // Exercise various basics over a range of numSamples samples using a LCG generator.
    
    LCG generator(seed);

    for (int i = 0; i < numSamples; i++)
    {
        uint32_t u = generator;

        printf("-- u: 0x%08x ------------------------------------------------\n", u);
        
        printf("  ToInt32(u, 10)    = %2d\n", ToInt32(u, 10));
        printf("  ToFloat(u, 100)   = %5.2f\n", ToFloat(u, 100.0f));

        Vec2f v2 = ToDisc(u);
        printf("  ToDisc(u)       = %5.2f, %5.2f\n", v2.x, v2.y);

        Vec3f v3 = ToBall(u);
        printf("  ToBall(u)       = %5.2f, %5.2f, %5.2f\n", v3.x, v3.y, v3.z);

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
    printf("\n");
    
    // Test examples
    float score     = ToFloat(generator, 1.0f, 100.0f);
    int   modifier  = ToInt32Signed(generator, 5);
    int   dayOfYear = ToInt32Inclusive(generator, 1, 365);
    float weightKG  = ToFloat(ModGaussLike(generator), 50.0f, 130.0f);

    printf("score = %4.1f, modifier = %2d, dayOfYear = %3d, weightKG = %5.1f\n",
        score     ,
        modifier  ,
        dayOfYear ,
        weightKG  
    );
    
    Vec2f discLoc       = ToDisc(generator);
    Vec2f pixTentSample = ToSquare(ModTriangle(generator), ModTriangle(generator));
    Vec3f rayDir        = ToDir3(generator);

    printf("discLoc       = [%g, %g]\n", discLoc.x, discLoc.y);
    printf("pixTentSample = [%g, %g]\n", pixTentSample.x, pixTentSample.y);
    printf("rayDir        = [%g, %g, %g]\n", rayDir.x, rayDir.y, rayDir.z);
    
    // Create an SVG showing samples distributed over a circle with triangle-ramp-down
    // radial distribution.    
    FILE* svgFile = fopen("distribute.svg", "w");

    int wh = 256;
    svg_header(svgFile, wh);
    float size = wh;
    svg_text(svgFile, Vec2f(0.0f, 20.0f), "ToDisc(u, ModTriangle(u))", "black");

    svg_line(svgFile, Vec2f(0.0f, 0.0f), Vec2f(size, 0.0f));
    svg_line(svgFile, Vec2f(size, 0.0f), Vec2f(size, size));
    svg_line(svgFile, Vec2f(size, size), Vec2f(0.0f, size));
    svg_line(svgFile, Vec2f(0.0f, size), Vec2f(0.0f, 0.0f));

    for (int i = 0; i < numSamples; i++)
    {
        Vec2f v = ToDisc(generator, ModHalfDown(ModTriangle(generator))); // use theta/r version, so we can modulate just r
        v = 0.5f * (v + Vec2f(1.0f, 1.0f));
        v = v * size;
        svg_circle(svgFile, v, 2.0f, "green");
    }

    svg_trailer(svgFile);
    fclose(svgFile);

    return 0;
}
