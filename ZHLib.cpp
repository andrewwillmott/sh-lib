//------------------------------------------------------------------------------
// Zonal Harmonics Utility Library
//------------------------------------------------------------------------------

#include "ZHLib.h"

using namespace ZHL;

namespace
{
    const float kFourPi = 4 * vl_pi;

    // The 'central' (z-symmetric, m=0) basis function set is a legendre poly in z:
    //   1, z, (3z^2 - 1), (5z^3 - 3z), 3 - 30z^2 + 35z^4 ...
    // Sometimes referred to as zonal harmonics.

    // The ZH normalization constants are just sqrtf( (2l + 1) / (4 pi) )
    // We also roll in any constant factor from the polynomial.
    // The constant factor for these is easy -- just evaluate with z = 1
    //                                                      N      P(z) / N
    const float kZH_Y_0 = sqrtf( 1 / (      kFourPi));  //         1
    const float kZH_Y_1 = sqrtf( 3 / (      kFourPi));  //         z
    const float kZH_Y_2 = sqrtf( 5 / (  4 * kFourPi));  // 1/2     (3 z^2 - 1)
    const float kZH_Y_3 = sqrtf( 7 / (  4 * kFourPi));  // 1/2     (5 z^3 - 3 z)
    const float kZH_Y_4 = sqrtf( 9 / ( 64 * kFourPi));  // 1/8     (35 z^4 - 30 z^2 + 3)
    const float kZH_Y_5 = sqrtf(11 / ( 64 * kFourPi));  // 1/8     (63 z^5 - 70 z^3 + 15 z)
    const float kZH_Y_6 = sqrtf(13 / (256 * kFourPi));  // 1/16    (231 z^6 - 315 z^4 + 105 z^2 - 5)

    const float kZH_Ym_0 = 1.0    / kZH_Y_0;
    const float kZH_Ym_1 = 1.0    / kZH_Y_1;
    const float kZH_Ym_2 = 0.5    / kZH_Y_2;
    const float kZH_Ym_3 = 0.5    / kZH_Y_3;
    const float kZH_Ym_4 = 0.125  / kZH_Y_4;
    const float kZH_Ym_5 = 0.125  / kZH_Y_5;
    const float kZH_Ym_6 = 0.0625 / kZH_Y_6;

    // diffuse reflectance convolution weights.
    // Multiply bands by these constants to convert incident radiance to exit
    // radiance for a diffuse surface.
    //const float kSH_A0 = 1.0f;
    const float kSH_A1 = (2.0f / 3.0f);
    const float kSH_A2 = (1.0f / 4.0f);
    //const float kSH_A3 = 0.0f;
    const float kSH_A4 = -(1.0f / 24.0f);

    // Integrating f(z) over the sphere is:
    //   integral[-1;1]  f(z) . 2pi dz
}

namespace
{
    inline float ZH_K(int l)
    {
        return sqrtf(float(2 * l + 1) / (4.0f * vl_pi));
    }

    float ZH_P(int l, float x)
    {
        float pm0 = 1.0f;

        if (l == 0)
            return pm0;

        float pm1 = x * pm0;

        if (l == 1)
            return pm1;

        float pm2 = 0.0f;

        for (int i = 2; i <= l; i++)
        {
            pm2 = (float(2 * i - 1) * x * pm1 - float(i - 1) * pm0) / float(i);
            pm0 = pm1;
            pm1 = pm2;
        }

        return pm2;
    }
}

float ZHL::ZH(int l, float z)   // dir must be normalised
{
    return ZH_K(l) * ZH_P(l, z);
}

const float ZHL::kZH_Y[7] =
{
    sqrtf( 1 / (      kFourPi)),  //         1
    sqrtf( 3 / (      kFourPi)),  //         z
    sqrtf( 5 / (  4 * kFourPi)),  // 1/2     (3 z^2 - 1)
    sqrtf( 7 / (  4 * kFourPi)),  // 1/2     (5 z^3 - 3 z)
    sqrtf( 9 / ( 64 * kFourPi)),  // 1/8     (35 z^4 - 30 z^2 + 3)
    sqrtf(11 / ( 64 * kFourPi)),  // 1/8     (63 z^5 - 70 z^3 + 15 z)
    sqrtf(13 / (256 * kFourPi)),  // 1/16    (231 z^6 - 315 z^4 + 105 z^2 - 5)
};

// Main ZH utils

void ZHL::CalcZHWeights(float z, int numBands, float w[7])
{
    float z2, z3, z4, z5, z6;

    z2 = z *  z;
    z3 = z * z2;
    if (numBands > 4)
    {
        z4 = z2 * z2;
        z5 = z2 * z3;
        z6 = z3 * z3;
    }

    switch (numBands)
    {
    case 7:
        w[6] = kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5);
    case 6:
        w[5] = kZH_Y_5 * (63 * z5 - 70 * z3 + 15 * z);
    case 5:
        w[4] = kZH_Y_4 * (35 * z4 - 30 * z2 + 3);
    case 4:
        w[3] = kZH_Y_3 * (5 * z3 - 3 * z);
    case 3:
        w[2] = kZH_Y_2 * (3 * z2 - 1);
    case 2:
        w[1] = kZH_Y_1 * z;
    case 1:
        w[0] = kZH_Y_0;
    }
}

float ZHL::SampleZH(float z, int numBands, const float zcoeffs[])
{
    float z2, z3, z4, z5, z6;

    z2 = z *  z;
    z3 = z * z2;
    if (numBands > 4)
    {
        z4 = z2 * z2;
        z5 = z2 * z3;
        z6 = z3 * z3;
    }

    float result = kZH_Y_0 * zcoeffs[0];

    switch (numBands)
    {
    case 7:
        result += zcoeffs[6] * kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5);
    case 6:
        result += zcoeffs[5] * kZH_Y_5 * (63  * z5 -  70 * z3 +  15 * z);
    case 5:
        result += zcoeffs[4] * kZH_Y_4 * (35  * z4 -  30 * z2 +  3);
    case 4:
        result += zcoeffs[3] * kZH_Y_3 * (5 * z3 - 3 * z);
    case 3:
        result += zcoeffs[2] * kZH_Y_2 * (3 * z2 - 1);
    case 2:
        result += zcoeffs[1] * kZH_Y_1 * z;
    }

    return result;
}

float ZHL::SampleZH_p1(int numBands, const float zcoeffs[7])
{
    float result = kZH_Y_0 * zcoeffs[0];

    switch (numBands)
    {
    case 7:
        result += kZH_Y_6 * 16 * zcoeffs[6];
    case 6:
        result += kZH_Y_5 * 8  * zcoeffs[5];
    case 5:
        result += kZH_Y_4 * 8  * zcoeffs[4];
    case 4:
        result += kZH_Y_3 * 2  * zcoeffs[3];
    case 3:
        result += kZH_Y_2 * 2  * zcoeffs[2];
    case 2:
        result += kZH_Y_1      * zcoeffs[1];
    }

    return result;
}

void ZHL::AddZHSample(float x, float z, int numBands, float zcoeffs[])
{
    float z2, z3, z4, z5, z6;

    z2 = z *  z;
    z3 = z * z2;
    if (numBands > 4)
    {
        z4 = z2 * z2;
        z5 = z2 * z3;
        z6 = z3 * z3;
    }

    switch (numBands)
    {
    case 7:
        zcoeffs[6] += x * kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5);
    case 6:
        zcoeffs[5] += x * kZH_Y_5 * (63  * z5 -  70 * z3 +  15 * z);
    case 5:
        zcoeffs[4] += x * kZH_Y_4 * (35  * z4 -  30 * z2 +  3);
    case 4:
        zcoeffs[3] += x * kZH_Y_3 * (5 * z3 - 3 * z);
    case 3:
        zcoeffs[2] += x * kZH_Y_2 * (3 * z2 - 1);
    case 2:
        zcoeffs[1] += x * kZH_Y_1 * z;
    case 1:
        zcoeffs[0] += x * kZH_Y_0;
    }
}

void ZHL::CalcZH7Weights(float z, float w[7])
{
    float z2 = z * z;
    float z3 = z2 * z;
    float z4 = z2 * z2;
    float z5 = z2 * z3;
    float z6 = z3 * z3;

    w[0] = kZH_Y_0;
    w[1] = kZH_Y_1 * z;
    w[2] = kZH_Y_2 * (3 * z2 - 1);
    w[3] = kZH_Y_3 * (5 * z3 - 3 * z);
    w[4] = kZH_Y_4 * (35 * z4 - 30 * z2 + 3);
    w[5] = kZH_Y_5 * (63 * z5 - 70 * z3 + 15 * z);
    w[6] = kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5);
}

float ZHL::SampleZH7(float z, const float zcoeffs[7])
{
    float z2 = z * z;
    float z3 = z2 * z;
    float z4 = z2 * z2;
    float z5 = z3 * z2;
    float z6 = z3 * z3;

    float result;

    result  = kZH_Y_0 * zcoeffs[0];
    result += kZH_Y_1 * zcoeffs[1] * z;
    result += kZH_Y_2 * zcoeffs[2] * (3 * z2 - 1);
    result += kZH_Y_3 * zcoeffs[3] * (5 * z3 - 3 * z);
    result += kZH_Y_4 * zcoeffs[4] * (35  * z4 -  30 * z2 +  3);
    result += kZH_Y_5 * zcoeffs[5] * (63  * z5 -  70 * z3 +  15 * z);
    result += kZH_Y_6 * zcoeffs[6] * (231 * z6 - 315 * z4 + 105 * z2 - 5);

    return result;
}

float ZHL::SampleZH7_p1(const float zcoeffs[7])
{
    float result;

    result  = kZH_Y_0      * zcoeffs[0];
    result += kZH_Y_1      * zcoeffs[1];
    result += kZH_Y_2 * 2  * zcoeffs[2];
    result += kZH_Y_3 * 2  * zcoeffs[3];
    result += kZH_Y_4 * 8  * zcoeffs[4];
    result += kZH_Y_5 * 8  * zcoeffs[5];
    result += kZH_Y_6 * 16 * zcoeffs[6];

    return result;
}

void ZHL::AddZH7Sample(float x, float z, float zcoeffs[7])
{
    float z2 = z * z;
    float z3 = z2 * z;
    float z4 = z2 * z2;
    float z5 = z3 * z2;
    float z6 = z3 * z3;

    zcoeffs[0] += kZH_Y_0 * x;
    zcoeffs[1] += kZH_Y_1 * z * x;
    zcoeffs[2] += kZH_Y_2 * (3 * z2 - 1) * x;
    zcoeffs[3] += kZH_Y_3 * (5 * z3 - 3 * z) * x;
    zcoeffs[4] += kZH_Y_4 * (35  * z4 -  30 * z2 +  3) * x;
    zcoeffs[5] += kZH_Y_5 * (63  * z5 -  70 * z3 +  15 * z) * x;
    zcoeffs[6] += kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5) * x;
}

void ZHL::CalcCosPowerZH7(int n, float zcoeffs[7])
{
    // We require n as an integer in this case so we can optimise -- either the
    // odd or even set of bands will be zero.
    if (n & 1)
    {
        zcoeffs[0] = 0.0f;
        zcoeffs[1] =   2.0f / (n + 2);
        zcoeffs[2] = 0.0f;
        zcoeffs[3] =  10.0f / (n + 4) -   6.0f / (n + 2);
        zcoeffs[4] = 0.0f;
        zcoeffs[5] = 126.0f / (n + 6) - 140.0f / (n + 4) +  30.0f / (n + 2);
        zcoeffs[6] = 0.0f;
    }
    else
    {
        zcoeffs[0] =   2.0f / (n + 1);
        zcoeffs[1] = 0.0f;
        zcoeffs[2] =   6.0f / (n + 3) -   2.0f / (n + 1);
        zcoeffs[3] = 0.0f;
        zcoeffs[4] =  70.0f / (n + 5) -  60.0f / (n + 3) +   6.0f / (n + 1);
        zcoeffs[5] = 0.0f;
        zcoeffs[6] = 462.0f / (n + 7) - 630.0f / (n + 5) + 210.0f / (n + 3) - 10.0f / (n + 1);
    }

    // apply norm constants
    zcoeffs[0] *= vl_twoPi * kZH_Y_0;
    zcoeffs[1] *= vl_twoPi * kZH_Y_1;
    zcoeffs[2] *= vl_twoPi * kZH_Y_2;
    zcoeffs[3] *= vl_twoPi * kZH_Y_3;
    zcoeffs[4] *= vl_twoPi * kZH_Y_4;
    zcoeffs[5] *= vl_twoPi * kZH_Y_5;
    zcoeffs[6] *= vl_twoPi * kZH_Y_6;
}

void ZHL::CalcCosPowerSatZH7(float n, float zcoeffs[7])
{
    zcoeffs[0] =   1.0f / (n + 1);
    zcoeffs[1] =   1.0f / (n + 2);
    zcoeffs[2] =   3.0f / (n + 3) -   1.0f / (n + 1);
    zcoeffs[3] =   5.0f / (n + 4) -   3.0f / (n + 2);
    zcoeffs[4] =  35.0f / (n + 5) -  30.0f / (n + 3) +   3.0f / (n + 1);
    zcoeffs[5] =  63.0f / (n + 6) -  70.0f / (n + 4) +  15.0f / (n + 2);
    zcoeffs[6] = 231.0f / (n + 7) - 315.0f / (n + 5) + 105.0f / (n + 3) - 5.0f / (n + 1);

    // apply norm constants
    zcoeffs[0] *= vl_twoPi * kZH_Y_0;
    zcoeffs[1] *= vl_twoPi * kZH_Y_1;
    zcoeffs[2] *= vl_twoPi * kZH_Y_2;
    zcoeffs[3] *= vl_twoPi * kZH_Y_3;
    zcoeffs[4] *= vl_twoPi * kZH_Y_4;
    zcoeffs[5] *= vl_twoPi * kZH_Y_5;
    zcoeffs[6] *= vl_twoPi * kZH_Y_6;
}

void ZHL::CalcHemisphereZH7(float zcoeffs[7])
{
    // projection of cos(theta) >= 0 ? 1 : 0
    // by definition only anti-symmetrical z terms are non-zero.
    // before weighting, they're basically +- 0.5
    zcoeffs[0] =  1.0f  * vl_twoPi * kZH_Y_0;
    zcoeffs[1] =  0.5f  * vl_twoPi * kZH_Y_1;
    zcoeffs[2] =  0.0f;
    zcoeffs[3] = -0.25f * vl_twoPi * kZH_Y_3;
    zcoeffs[4] =  0.0f;
    zcoeffs[5] =  0.5f  * vl_twoPi * kZH_Y_5;
    zcoeffs[6] =  0.0f;
}

void ZHL::CalcGatedSpotZH7(float t, float zcoeffs[7])
{
    // Just the integral of P_i(z) 2pi dz over [t, 1].
    float t2 = sqr(t);
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;
    float t6 = t3 * t3;
    float t7 = t4 * t3;

    zcoeffs[0] = 1 - t;
    zcoeffs[1] = 0.5f * (1.0f - t2); // [z^2 /2, t, 1] = 1/2 - t^2/2
    zcoeffs[2] = t - t3;
    zcoeffs[3] = 0.25f * (-1.0f + t2 * 6.0f - 5.0f * t4);
    zcoeffs[4] = -3.0f * t + 10.0f * t3 - 7.0f * t5;
    zcoeffs[5] =  0.5f * (1.0f - 15.0f * t2 + 35.0f * t4 - 21.0f * t6);
    zcoeffs[6] =  5.0f * t - 35.0f * t3 + 63.0f * t5 - 33.0f * t7;

    // apply norm constants
    zcoeffs[0] *= vl_twoPi * kZH_Y_0;
    zcoeffs[1] *= vl_twoPi * kZH_Y_1;
    zcoeffs[2] *= vl_twoPi * kZH_Y_2;
    zcoeffs[3] *= vl_twoPi * kZH_Y_3;
    zcoeffs[4] *= vl_twoPi * kZH_Y_4;
    zcoeffs[5] *= vl_twoPi * kZH_Y_5;
    zcoeffs[6] *= vl_twoPi * kZH_Y_6;
}

void ZHL::CalcGatedCosZH7(float t, float zcoeffs[7])
{
    // The integral of z P_i(z) 2pi dz over [t, 1].
    float t2 = sqr(t);
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;
    float t6 = t3 * t3;
    float t7 = t4 * t3;
    float t8 = t4 * t4;

    zcoeffs[0] = 0.5f * (1.0f - t2);
    zcoeffs[1] = (1.0f - t3) / 3.0f;
    zcoeffs[2] = (1.0f + 2.0f * t2 - 3.0f * t4) / 4.0f;
    zcoeffs[3] = t3 - t5;
    zcoeffs[4] = -(1.0f + 9.0f * t2 - 45.0f * t4 + 35.0f * t6) / 6.0f;
    zcoeffs[5] = -5.0f * t3 + 14.0f * t5 - 9.0f * t7;
    zcoeffs[6] = 0.125f * (1.0f + 20.0f * t2 - 210.0f * t4 + 420.0f * t6 - 231.0f * t8);

    // apply norm constants
    zcoeffs[0] *= vl_twoPi * kZH_Y_0;
    zcoeffs[1] *= vl_twoPi * kZH_Y_1;
    zcoeffs[2] *= vl_twoPi * kZH_Y_2;
    zcoeffs[3] *= vl_twoPi * kZH_Y_3;
    zcoeffs[4] *= vl_twoPi * kZH_Y_4;
    zcoeffs[5] *= vl_twoPi * kZH_Y_5;
    zcoeffs[6] *= vl_twoPi * kZH_Y_6;
}

void ZHL::CalcGatedCosBands5(float t, float bandScale[5])
{
    // produces band weights for lerp(cos, 1, t) convolved with target.
    // the normalisation coefficients drop out.
    if (t >= 1.0f)
    {
        for (int i = 0; i < 5; i++)
          bandScale[i] = 1.0f;
        return;
    }

    float t2 = sqr(t);
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;
    float t6 = t3 * t3;

    bandScale[0] = 0.5f * (1.0f - t2);
    bandScale[1] = (1.0f - t3) / 3.0f;
    bandScale[2] = (1.0f + 2.0f * t2 - 3.0f * t4) / 4.0f;
    bandScale[3] = t3 - t5;
    bandScale[4] = -(1.0f + 9.0f * t2 - 45.0f * t4 + 35.0f * t6) / 6.0f;

    float invB0 = 1.0f / bandScale[0];

    bandScale[0] = 1.0f;
    bandScale[1] *= invB0;
    bandScale[2] *= invB0 * 0.5f;
    bandScale[3] *= invB0 * 0.5f;
    bandScale[4] *= invB0 * 0.125f;
}

void ZHL::MirrorZH(int n, float* coeffs)
{
    for (int l = 1; l < n; l += 2)
        coeffs[l] = -coeffs[l];
}

void ZHL::MirrorZH(int n, Vec4f* coeffs)
{
    for (int l = 1; l < n; l += 2)
        coeffs[l] = -coeffs[l];
}

void ZHL::ApplyDiffuseBRDFZH7(float zcoeffs[7])
{
    zcoeffs[1] *= kSH_A1;
    zcoeffs[2] *= kSH_A2;
    zcoeffs[3] = 0;
    zcoeffs[4] *= kSH_A4;
    zcoeffs[5] = 0;
    zcoeffs[6] = 0;
}

void ZHL::ApplyDiffuseBRDFZH7(const float zcoeffsIn[], float zcoeffsOut[])
{
    zcoeffsOut[0] = zcoeffsIn[0];
    zcoeffsOut[1] = kSH_A1 * zcoeffsIn[1];
    zcoeffsOut[2] = kSH_A2 * zcoeffsIn[2];
    zcoeffsOut[3] = 0;
    zcoeffsOut[4] = kSH_A4 * zcoeffsIn[4];
    zcoeffsOut[5] = 0;
    zcoeffsOut[6] = 0;
}

void ZHL::CalcHGPhaseZH(float g, float weight, int n, float zcoeffs[])
{
    // The Heyney-Greenstein phase function turns out to be a
    // simple power series in g when decomposed into zonal harmonics.
    for (int i = 0; i < n; i++)
    {
        zcoeffs[i] = sqrtf(4.0f * vl_pi * float(i * 2 + 1)) * weight;
        weight *= g;
    }
}

void ZHL::CalcNormalizedHGPhaseZH(float g, float weight, int n, float zcoeffs[])
{
    // normalize so that max of function is always 1.
    weight *= sqr(1 - fabsf(g)) / (1 + fabsf(g));

    for (int i = 0; i < n; i++)
    {
        zcoeffs[i] = sqrtf(4.0f * vl_pi * float(i * 2 + 1)) * weight;
        weight *= g;
    }
}

void ZHL::ConvolveZH(int n, const float brdfCoeffs[], const float zhCoeffsIn[], float zhCoeffsOut[])
{
    for (int i = 0; i < n; i++)
    {
        int n = (2 * i + 1);
        float alpha = sqrtf(4.0f * vl_pi / n);

        zhCoeffsOut[i] = zhCoeffsIn[i] * brdfCoeffs[i] * alpha;
    }
}

void ZHL::ApplyZHMaxScale(int n, float coeffs[])
{
    VL_ASSERT(n > 0);

    (*coeffs++) *= kZH_Ym_0;

    if (n < 2)
        return;

    (*coeffs++) *= kZH_Ym_1;

    if (n < 3)
        return;

    (*coeffs++) *= kZH_Ym_2;

    if (n < 4)
        return;

    (*coeffs++) *= kZH_Ym_3;

    if (n < 5)
        return;

    (*coeffs++) *= kZH_Ym_4;

    if (n < 6)
        return;

    (*coeffs++) *= kZH_Ym_5;

    if (n < 7)
        return;

    (*coeffs++) *= kZH_Ym_6;

    VL_ASSERT(n <= 7);
}


void ZHL::ApplyZHWindowing(float gamma, int n, float* coeffs)
{
    for (int i = 0; i < n; i++)
        coeffs[i] *= WindowScale(i, gamma);
}

void ZHL::ApplyZHWindowing(float gamma, int n, Vec4f* coeffs)
{
    for (int i = 0; i < n; i++)
        coeffs[i] *= WindowScale(i, gamma);
}

namespace
{
    struct ZHBasisTriple
    {
        int8_t i;
        int8_t j;
        int8_t k;
        int8_t m;  // mode
        float  s;
    };

    const ZHBasisTriple kZHBasisTriples2[] =
    {
        { 0, 0, 0, 0, 0.28209479 },
        { 0, 1, 1, 1, 0.28209479 },
        { -1, -1, -1, -1, 0.0 }
    };

    const ZHBasisTriple kZHBasisTriples3[] =
    {
        { 0, 0, 0, 0, 0.28209479 },
        { 0, 1, 1, 1, 0.28209479 },
        { 0, 2, 2, 1, 0.28209479 },
        { 1, 1, 2, 2, 0.25231325 },
        { 2, 2, 2, 0, 0.18022375 },
        { -1, -1, -1, -1, 0.0 }
    };

    const ZHBasisTriple kZHBasisTriples4[] =
    {
        { 0, 0, 0, 0, 0.28209479 },
        { 0, 1, 1, 1, 0.28209479 },
        { 0, 2, 2, 1, 0.28209479 },
        { 0, 3, 3, 1, 0.28209479 },
        { 1, 1, 2, 2, 0.25231325 },
        { 1, 2, 3, 3, 0.24776670 },
        { 2, 2, 2, 0, 0.18022375 },
        { 2, 3, 3, 1, 0.16820883 },
        { -1, -1, -1, -1, 0.0 }
    };

    const ZHBasisTriple kZHBasisTriples5[] =
    {
        { 0, 0, 0, 0, 0.28209479 },
        { 0, 1, 1, 1, 0.28209479 },
        { 0, 2, 2, 1, 0.28209479 },
        { 0, 3, 3, 1, 0.28209479 },
        { 0, 4, 4, 1, 0.28209479 },
        { 1, 1, 2, 2, 0.25231325 },
        { 1, 2, 3, 3, 0.24776670 },
        { 1, 3, 4, 3, 0.24623252 },
        { 2, 2, 2, 0, 0.18022375 },
        { 2, 2, 4, 2, 0.24179554 },
        { 2, 3, 3, 1, 0.16820883 },
        { 2, 4, 4, 1, 0.16383977 },
        { 3, 3, 4, 2, 0.15386989 },
        { 4, 4, 4, 0, 0.13696111 },
        { -1, -1, -1, -1, 0.0 }
    };

    const ZHBasisTriple kZHBasisTriples6[] =
    {
        { 0, 0, 0, 0, 0.28209479 },
        { 0, 1, 1, 1, 0.28209479 },
        { 0, 2, 2, 1, 0.28209479 },
        { 0, 3, 3, 1, 0.28209479 },
        { 0, 4, 4, 1, 0.28209479 },
        { 0, 5, 5, 1, 0.28209479 },
        { 1, 1, 2, 2, 0.25231325 },
        { 1, 2, 3, 3, 0.24776670 },
        { 1, 3, 4, 3, 0.24623252 },
        { 1, 4, 5, 3, 0.24553200 },
        { 2, 2, 2, 0, 0.18022375 },
        { 2, 2, 4, 2, 0.24179554 },
        { 2, 3, 3, 1, 0.16820883 },
        { 2, 3, 5, 3, 0.23961470 },
        { 2, 4, 4, 1, 0.16383977 },
        { 2, 5, 5, 1, 0.16173926 },
        { 3, 3, 4, 2, 0.15386989 },
        { 3, 4, 5, 3, 0.14837393 },
        { 4, 4, 4, 0, 0.13696111 },
        { 4, 5, 5, 1, 0.13019760 },
        { -1, -1, -1, -1, 0.0 }
    };

    const ZHBasisTriple kZHBasisTriples7[] =
    {
        { 0, 0, 0, 0, 0.28209479 },
        { 0, 1, 1, 1, 0.28209479 },
        { 0, 2, 2, 1, 0.28209479 },
        { 0, 3, 3, 1, 0.28209479 },
        { 0, 4, 4, 1, 0.28209479 },
        { 0, 5, 5, 1, 0.28209479 },
        { 0, 6, 6, 1, 0.28209479 },
        { 1, 1, 2, 2, 0.25231325 },
        { 1, 2, 3, 3, 0.24776670 },
        { 1, 3, 4, 3, 0.24623252 },
        { 1, 4, 5, 3, 0.24553200 },
        { 1, 5, 6, 3, 0.24515397 },
        { 2, 2, 2, 0, 0.18022375 },
        { 2, 2, 4, 2, 0.24179554 },
        { 2, 3, 3, 1, 0.16820883 },
        { 2, 3, 5, 3, 0.23961470 },
        { 2, 4, 4, 1, 0.16383977 },
        { 2, 4, 6, 3, 0.23856513 },
        { 2, 5, 5, 1, 0.16173926 },
        { 2, 6, 6, 1, 0.16056298 },
        { 3, 3, 4, 2, 0.15386989 },
        { 3, 3, 6, 2, 0.23708793 },
        { 3, 4, 5, 3, 0.14837393 },
        { 3, 5, 6, 3, 0.14563067 },
        { 4, 4, 4, 0, 0.13696111 },
        { 4, 4, 6, 2, 0.14225276 },
        { 4, 5, 5, 1, 0.13019760 },
        { 4, 6, 6, 1, 0.12671638 },
        { 5, 5, 6, 2, 0.12272787 },
        { 6, 6, 6, 0, 0.11450687 },
        { -1, -1, -1, -1, 0.0 }
    };

    template<typename T_SRC, typename T_DST> inline void MultiplyZH(int n, const T_SRC* a, const T_DST* b, T_DST* c, const ZHBasisTriple* triples)
    /// Multiplies two sets of coeffs together using given triple table.
    {
        VL_ASSERT((void*) a != c && b != c);

        for (int i = 0; i < n; i++)
            c[i] = vl_0;

        while (triples->i >= 0)
        {
            VL_ASSERT(triples->i < n && triples->j < n && triples->k < n);

            switch (triples->m)
            {
            case 0: // i == j == k
                c[triples->k] += triples->s * a[triples->i] * b[triples->j];
                break;
            case 1: // j == k
                c[triples->i] += triples->s * a[triples->j] * b[triples->j];
                c[triples->j] += triples->s * (a[triples->k] * b[triples->i] + a[triples->i] * b[triples->k]);
                break;
            case 2: // i == j
                c[triples->k] += triples->s * a[triples->i] * b[triples->j];
                c[triples->i] += triples->s * (a[triples->j] * b[triples->k] + a[triples->k] * b[triples->j]);
                break;
            case 3: // i != j != k
                c[triples->i] += triples->s * (a[triples->j] * b[triples->k] + a[triples->k] * b[triples->j]);
                c[triples->j] += triples->s * (a[triples->k] * b[triples->i] + a[triples->i] * b[triples->k]);
                c[triples->k] += triples->s * (a[triples->i] * b[triples->j] + a[triples->j] * b[triples->i]);
                break;
            }

            triples++;
        }
    }
}

void ZHL::MultiplyZH(int numBands, const float* a, const float* b, float* c)
{
    switch (numBands)
    {
    case 1: c[0] = a[0] * b[0] * kZH_Y_0; break;
    case 2: MultiplyZH(2, a, b, c, kZHBasisTriples2); break;
    case 3: MultiplyZH(3, a, b, c, kZHBasisTriples3); break;
    case 4: MultiplyZH(4, a, b, c, kZHBasisTriples4); break;
    case 5: MultiplyZH(5, a, b, c, kZHBasisTriples5); break;
    case 6: MultiplyZH(6, a, b, c, kZHBasisTriples6); break;
    case 7: MultiplyZH(7, a, b, c, kZHBasisTriples7); break;
    default:
        VL_ERROR("Unhandled band count\n");
    }
}

void ZHL::MultiplyZH(int numBands, const Vec4f* a, const Vec4f* b, Vec4f* c)
{
    switch (numBands)
    {
    case 1: c[0] = a[0] * b[0] * kZH_Y_0; break;
    case 2: MultiplyZH(2, a, b, c, kZHBasisTriples2); break;
    case 3: MultiplyZH(3, a, b, c, kZHBasisTriples3); break;
    case 4: MultiplyZH(4, a, b, c, kZHBasisTriples4); break;
    case 5: MultiplyZH(5, a, b, c, kZHBasisTriples5); break;
    case 6: MultiplyZH(6, a, b, c, kZHBasisTriples6); break;
    case 7: MultiplyZH(7, a, b, c, kZHBasisTriples7); break;
    default:
        VL_ERROR("Unhandled band count\n");
    }
}

// Misc lighting

namespace
{
    void GetTerminatorZH7(float* zcoeffs)
    {
        // projection of z = cos(theta); f(z) = ( z > cutoff ) ? 1.0f : max<float>( 0.0f, z + cutoff ) / (2.0f * cutoff);
        zcoeffs[0] =  1.77245322546612f;
        zcoeffs[1] =  1.51447090327194f;
        zcoeffs[2] = -6.08337227802301e-05f;
        zcoeffs[3] = -0.540183395337511f;
        zcoeffs[4] =  8.55535266573624e-05f;
        zcoeffs[5] =  0.297859938699253f;
        zcoeffs[6] = -0.000439377751693618f;
    }
}

void ZHL::CalcAtmosphereZH
(
    int          numPhases,
    const Vec4f* colourPhases,
    int          bands,
    Vec4f*       zcoeffs
)
{
    VL_ASSERT(bands >=1 && bands <= 5);

    Vec4f phaseCoeffs[5] = {};
    float zh_coeffs[5];

    for (int i = 1; i < numPhases; i++)
    {
        CalcHGPhaseZH(colourPhases[i][3], 1.0f, bands, zh_coeffs);

        for (int i = 0; i < bands; i++)
            phaseCoeffs[i] += colourPhases[i] * zh_coeffs[i];
    }

    float termZHCoeffs[7];
    GetTerminatorZH7(termZHCoeffs);

    switch (bands)
    {
    case 3:
        ::MultiplyZH<float, Vec4f>(bands, termZHCoeffs, phaseCoeffs, zcoeffs, kZHBasisTriples3);
        break;
    case 4:
        ::MultiplyZH<float, Vec4f>(bands, termZHCoeffs, phaseCoeffs, zcoeffs, kZHBasisTriples4);
        break;
    case 5:
        ::MultiplyZH<float, Vec4f>(bands, termZHCoeffs, phaseCoeffs, zcoeffs, kZHBasisTriples5);
        break;
    default:
        VL_ERROR("unhandled bands\n");
    }
}

