//------------------------------------------------------------------------------
// SH Library Test
//------------------------------------------------------------------------------

#define VL_DEBUG

#include "SHLib.h"
#include "ZHLib.h"

using namespace SHL;
using namespace ZHL;

int main(int argc, const char* argv[])
{
    float sample;

    sample = SH(1, -1, vl_y);
    printf("SH(1, -1, e_y) = %f\n", sample);
    VL_ASSERT(sample == kZH_Y[1]);

    sample = SH(1, 0, vl_z);
    printf("SH(1,  0, e_z) = %f\n", sample);
    VL_ASSERT(sample == kZH_Y[1]);

    sample = SH(1, 1, vl_x);
    printf("SH(1, +1, e_x) = %f\n", sample);
    VL_ASSERT(sample == kZH_Y[1]);

    printf("\n");
    float zcoeffs[7];
    CalcCosPowerZH7(3, zcoeffs);  // saturate(cos(theta)^3)

    float coeffs[49];
    RotateZHToSH(vl_x, 7, zcoeffs, coeffs);

    sample = SampleSH(vl_x, 7, coeffs);
    printf("Sample(cos^3(x), +e_x) = %4.1f\n", sample);
    VL_ASSERT(fabsf(sample - 1.0f) < 1e-2f);

    sample = SampleSH(vl_minus_x, 7, coeffs);
    printf("Sample(cos^3(x), -e_x) = %4.1f\n", sample);
    VL_ASSERT(fabsf(sample - -1.0f) < 1e-2f);

    sample = SampleSH(vl_y, 7, coeffs);
    printf("Sample(cos^3(x), +e_y) = %4.1f\n", sample);
    VL_ASSERT(fabsf(sample - 0.0f) < 1e-2f);

	return 0;
}