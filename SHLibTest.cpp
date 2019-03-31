//------------------------------------------------------------------------------
// SH Library Test
//------------------------------------------------------------------------------

#include "SHLib.h"

using namespace SHL;

int main(int argc, const char* argv[])
{
    float zcoeffs[7];

    CalcCosPowerSatZH7(7, zcoeffs);

    float shcoeffs[49];
    RotateZHToSH(Vec3f(1.0f, 1.0f, 1.0f), 7, zcoeffs, shcoeffs);

    // Vec4f SampleColour(const Vec4f* coeffs, int numBands, const Vec3f& v);
    // SampleColour();

	return 0;
}