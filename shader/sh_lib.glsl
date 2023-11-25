//
// sh.sh
//
// Spherical Harmonics utilities
//

#ifndef SH_SH
#define SH_SH

// Unfortunately we can't pass subsets of arrays to GLSL functions,
// so this has to be fixed. It should be set to match the size of the uniforms
// array you have available.
#ifndef SH_BANDS
    #define SH_BANDS 5
#endif
#ifndef SH_TYPE
    #define SH_TYPE vec4
#endif
#ifndef CONST
    #define CONST(X) const X
#endif

#define SH_COEFFS (SH_BANDS * SH_BANDS)
#define ZH_COEFFS SH_BANDS

////////////////////////////////////////////////////////////////////////////////
// SH sampling
////////////////////////////////////////////////////////////////////////////////

CONST(float kSH_A0 = 1.0);
CONST(float kSH_A1 = 2.0 / 3.0);
CONST(float kSH_A2 = 1.0 / 4.0);
CONST(float kSH_A3 = 0.0);
CONST(float kSH_A4 = -1.0 / 24.0);

CONST(float kPi     = 3.14159265358979323846f);
CONST(float kTwoPi  = 2.0 * kPi);
CONST(float kFourPi = 4.0 * kPi);


//
// SH eval with normalisation constants already applied to the coefficients on
// the CPU side.
//

#define SH_PN_BAND_0(coeffs) (coeffs[0].xyz)
#define SH_PN_BAND_1(coeffs) \
    ( coeffs[1].xyz * x0.y   \
    + coeffs[2].xyz * x0.z   \
    + coeffs[3].xyz * x0.x   )

#define SH_PN_BAND_2(coeffs) \
    ( coeffs[4].xyz  * x2.x               \
    + coeffs[5].xyz  * x2.y               \
    + coeffs[6].xyz  * (3.0 * x1.z - 1.0) \
    + coeffs[7].xyz  * x2.z               \
    + coeffs[8].xyz  * (x1.x - x1.y)      )

#define SH_PN_BAND_3(coeffs) \
    ( coeffs[ 9].xyz * (3.0 * x1.x - x1.y)   * x0.y \
    + coeffs[10].xyz * x2.x                  * x0.z \
    + coeffs[11].xyz * (5.0 * x1.z - 1.0)    * x0.y \
    + coeffs[12].xyz * (5.0 * x1.z - 3.0)    * x0.z \
    + coeffs[13].xyz * (5.0 * x1.z - 1.0)    * x0.x \
    + coeffs[14].xyz * (x1.x - x1.y)         * x0.z \
    + coeffs[15].xyz * (x1.x - 3.0 * x1.y)   * x0.x )

#define SH_PN_BAND_4(coeffs) \
    ( coeffs[16].xyz * (x1.x - x1.y) * x2.x                            \
    + coeffs[17].xyz * (3.0 * x1.x - x1.y) * x2.y                      \
    + coeffs[18].xyz * (7.0 * x1.z - 1.0) * x2.x                       \
    + coeffs[19].xyz * (7.0 * x1.z - 3.0) * x2.y                       \
    + coeffs[20].xyz * (3.0 - 30.0 * x1.z + 35.0 * x1.z * x1.z)        \
    + coeffs[21].xyz * (7.0 * x1.z - 3.0) * x2.z                       \
    + coeffs[22].xyz * (7.0 * x1.z - 1.0) * (x1.x - x1.y)              \
    + coeffs[23].xyz * (x1.x - 3.0 * x1.y) * x2.z                      \
    + coeffs[24].xyz * (x1.x * x1.x + x1.y * x1.y - 6.0 * x1.x * x1.y) )

vec3 SH_PN(int n, vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    vec3 x1 = x0 * x0;
    vec3 x2 = x0 * x0.yzx;

    vec3       sh  = SH_PN_BAND_0(coeffs);
    if (n > 1) sh += SH_PN_BAND_1(coeffs);
#if SH_BANDS >= 3
    if (n > 2) sh += SH_PN_BAND_2(coeffs);
#endif
#if SH_BANDS >= 4
    if (n > 3) sh += SH_PN_BAND_3(coeffs);
#endif
#if SH_BANDS >= 5
    if (n > 4) sh += SH_PN_BAND_4(coeffs);
#endif
    return sh;
}

vec3 SH_PN_Diffuse(int n, vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    vec3 x1 = x0 * x0;
    vec3 x2 = x0 * x0.yzx;

    vec3       sh  = kSH_A0 * SH_PN_BAND_0(coeffs);
    if (n > 1) sh += kSH_A1 * SH_PN_BAND_1(coeffs);
#if SH_BANDS >= 3
    if (n > 2) sh += kSH_A2 * SH_PN_BAND_2(coeffs);
#endif
#if SH_BANDS >= 5
    if (n > 4) sh += kSH_A4 * SH_PN_BAND_4(coeffs);
#endif
    return sh;
}

vec3 SH_PN(vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    return SH_PN(SH_BANDS, x0, coeffs);
}

vec3 SH_PN_Diffuse(vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    return SH_PN_Diffuse(SH_BANDS, x0, coeffs);
}


//
// Normalisation constants in shader
//

CONST(float kSH_Y_00 = sqrt(  1.0 / (       kFourPi)));   // 1
CONST(float kSH_Y_10 = sqrt(  3.0 / (       kFourPi)));   // y
CONST(float kSH_Y_11 = sqrt(  3.0 / (       kFourPi)));   // z
CONST(float kSH_Y_12 = sqrt(  3.0 / (       kFourPi)));   // x
CONST(float kSH_Y_20 = sqrt( 15.0 / (       kFourPi)));   // xy
CONST(float kSH_Y_21 = sqrt( 15.0 / (       kFourPi)));   // yz
CONST(float kSH_Y_22 = sqrt(  5.0 / ( 4.0 * kFourPi)));   // (3z^2 - 1)
CONST(float kSH_Y_23 = sqrt( 15.0 / (       kFourPi)));   // xz
CONST(float kSH_Y_24 = sqrt( 15.0 / ( 4.0 * kFourPi)));   // (x^2 - y^2)
CONST(float kSH_Y_30 = sqrt( 35.0 / ( 8.0 * kFourPi)));   // y (3x^2 - y^2)
CONST(float kSH_Y_31 = sqrt(105.0 / ( 1.0 * kFourPi)));   // z xy
CONST(float kSH_Y_32 = sqrt( 21.0 / ( 8.0 * kFourPi)));   // y (5z^2 - 1)
CONST(float kSH_Y_33 = sqrt(  7.0 / ( 4.0 * kFourPi)));   // z (5z^2 - 3)
CONST(float kSH_Y_34 = sqrt( 21.0 / ( 8.0 * kFourPi)));   // x (5z^2 - 1)
CONST(float kSH_Y_35 = sqrt(105.0 / ( 4.0 * kFourPi)));   // z (x^2 - y^2)
CONST(float kSH_Y_36 = sqrt( 35.0 / ( 8.0 * kFourPi)));   // x (x^2 - 3y^2)
CONST(float kSH_Y_40 = sqrt(315.0 / ( 4.0 * kFourPi)));   // (x^2 - y^2)  xy
CONST(float kSH_Y_41 = sqrt(315.0 / ( 8.0 * kFourPi)));   // (3x^2 - y^2) yz
CONST(float kSH_Y_42 = sqrt( 45.0 / ( 4.0 * kFourPi)));   // (7z^2 - 1)   xy
CONST(float kSH_Y_43 = sqrt( 45.0 / ( 8.0 * kFourPi)));   // (7z^2 - 3)   yz
CONST(float kSH_Y_44 = sqrt(  9.0 / (64.0 * kFourPi)));   // (3 - 30z^2 + 35z^4)
CONST(float kSH_Y_45 = sqrt( 45.0 / ( 8.0 * kFourPi)));   // (7z^2 - 3)   xz
CONST(float kSH_Y_46 = sqrt( 45.0 / (16.0 * kFourPi)));   // (7z^2 - 1)   (x^2 -y^2)
CONST(float kSH_Y_47 = sqrt(315.0 / ( 8.0 * kFourPi)));   // (x^2 - 3y^2) xz
CONST(float kSH_Y_48 = sqrt(315.0 / (64.0 * kFourPi)));   // (x^4 - 6 x^2 y^2 + y^4)

#define SH_BAND_0(coeffs) (coeffs[0].xyz * kSH_Y_00)
#define SH_BAND_1(coeffs) \
    ( coeffs[1].xyz * kSH_Y_10 * x0.y   \
    + coeffs[2].xyz * kSH_Y_11 * x0.z   \
    + coeffs[3].xyz * kSH_Y_12 * x0.x   )

#define SH_BAND_2(coeffs) \
    ( coeffs[4].xyz * kSH_Y_20 * x2.x               \
    + coeffs[5].xyz * kSH_Y_21 * x2.y               \
    + coeffs[6].xyz * kSH_Y_22 * (3.0 * x1.z - 1.0) \
    + coeffs[7].xyz * kSH_Y_23 * x2.z               \
    + coeffs[8].xyz * kSH_Y_24 * (x1.x - x1.y)      )

#define SH_BAND_3(coeffs) \
    ( coeffs[ 9].xyz * kSH_Y_30 * (3.0 * x1.x - x1.y) * x0.y \
    + coeffs[10].xyz * kSH_Y_31 * x2.x                * x0.z \
    + coeffs[11].xyz * kSH_Y_32 * (5.0 * x1.z - 1.0)  * x0.y \
    + coeffs[12].xyz * kSH_Y_33 * (5.0 * x1.z - 3.0)  * x0.z \
    + coeffs[13].xyz * kSH_Y_34 * (5.0 * x1.z - 1.0)  * x0.x \
    + coeffs[14].xyz * kSH_Y_35 * (x1.x - x1.y)       * x0.z \
    + coeffs[15].xyz * kSH_Y_36 * (x1.x - 3.0 * x1.y) * x0.x )

#define SH_BAND_4(coeffs) \
    ( coeffs[16].xyz * kSH_Y_40 * (      x1.x - x1.y) * x2.x                     \
    + coeffs[17].xyz * kSH_Y_41 * (3.0 * x1.x - x1.y) * x2.y                     \
    + coeffs[18].xyz * kSH_Y_42 * (7.0 * x1.z - 1.0)  * x2.x                     \
    + coeffs[19].xyz * kSH_Y_43 * (7.0 * x1.z - 3.0)  * x2.y                     \
    + coeffs[20].xyz * kSH_Y_44 * (3.0 - 30.0 * x1.z + 35.0 * x1.z * x1.z)       \
    + coeffs[21].xyz * kSH_Y_45 * (7.0 * x1.z - 3.0)  * x2.z                     \
    + coeffs[22].xyz * kSH_Y_46 * (7.0 * x1.z - 1.0)  * (x1.x - x1.y)            \
    + coeffs[23].xyz * kSH_Y_47 * (x1.x - 3.0 * x1.y) * x2.z                     \
    + coeffs[24].xyz * kSH_Y_48 * (x1.x * x1.x + x1.y * x1.y - 6.0 * x1.x * x1.y))

vec3 SH(int n, vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    vec3 x1 = x0 * x0;
    vec3 x2 = x0 * x0.yzx;

    vec3       sh  = SH_BAND_0(coeffs);
    if (n > 1) sh += SH_BAND_1(coeffs);
#if SH_BANDS >= 3
    if (n > 2) sh += SH_BAND_2(coeffs);
#endif
#if SH_BANDS >= 4
    if (n > 3) sh += SH_BAND_3(coeffs);
#endif
#if SH_BANDS >= 5
    if (n > 4) sh += SH_BAND_4(coeffs);
#endif
    return sh;
}

vec3 SH_Diffuse(int n, vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    vec3 x1 = x0 * x0;
    vec3 x2 = x0 * x0.yzx;

    vec3       sh  = kSH_A0 * SH_BAND_0(coeffs);
    if (n > 1) sh += kSH_A1 * SH_BAND_1(coeffs);
#if SH_BANDS >= 3
    if (n > 2) sh += kSH_A2 * SH_BAND_2(coeffs);
#endif
#if SH_BANDS >= 5
    if (n > 4) sh += kSH_A4 * SH_BAND_4(coeffs);
#endif
    return sh;
}

vec3 SH(vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    return SH(SH_BANDS, x0, coeffs);
}

vec3 SH_Diffuse(vec3 x0, SH_TYPE coeffs[SH_COEFFS])
{
    return SH_Diffuse(SH_BANDS, x0, coeffs);
}


////////////////////////////////////////////////////////////////////////////////
// ZH construction/sample utilities
////////////////////////////////////////////////////////////////////////////////

CONST(float kZH_Y_0 = sqrt( 1.0 / (        kFourPi)));  //         1
CONST(float kZH_Y_1 = sqrt( 3.0 / (        kFourPi)));  //         z
CONST(float kZH_Y_2 = sqrt( 5.0 / (  4.0 * kFourPi)));  // 1/2     (3 z^2 - 1)
CONST(float kZH_Y_3 = sqrt( 7.0 / (  4.0 * kFourPi)));  // 1/2     (5 z^3 - 3 z)
CONST(float kZH_Y_4 = sqrt( 9.0 / ( 64.0 * kFourPi)));  // 1/8     (35 z^4 - 30 z^2 + 3)

float SampleZH(float z, const float zcoeffs[SH_BANDS])
{
    float z2 = z *  z;
    float z3 = z * z2;
    float z4 = z2 * z2;

    float result;
    result  = zcoeffs[0] * kZH_Y_0;
    result += zcoeffs[1] * kZH_Y_1 * z;
#if ZH_COEFFS >= 3
    result += zcoeffs[2] * kZH_Y_2 * (3.0 * z2 - 1.0);
#endif
#if ZH_COEFFS >= 4
    result += zcoeffs[3] * kZH_Y_3 * (5.0 * z3 - 3.0 * z);
#endif
#if ZH_COEFFS >= 5
    result += zcoeffs[4] * kZH_Y_4 * (35.0  * z4 -  30.0 * z2 +  3.0);
#endif

    return result;
}


struct ZH
{
    float c[ZH_COEFFS];
};

// Corresponds to constant light from z = cos(theta) = 't' to z = 1, 0 otherwise. Basically a spot light.
ZH GatedSpotZH(float t)
{
    // Just the integral of P_i(z) 2pi dz over [t, 1].
    float t2 = t * t;
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;

    ZH zh;
    zh.c[0] = kTwoPi * kZH_Y_0 * (1.0 - t);
    zh.c[1] = kTwoPi * kZH_Y_1 * 0.5 * (1.0 - t2); // [z^2 /2, t, 1] = 1/2 - t^2/2
#if ZH_COEFFS >= 3
    zh.c[2] = kTwoPi * kZH_Y_2 * (t - t3);
#endif
#if ZH_COEFFS >= 4
    zh.c[3] = kTwoPi * kZH_Y_3 * 0.25 * (-1.0 + t2 * 6.0 - 5.0 * t4);
#endif
#if ZH_COEFFS >= 5
    zh.c[4] = kTwoPi * kZH_Y_4 * (-3.0 * t + 10.0 * t3 - 7.0 * t5);
#endif

    return zh;
}

// Corresponds to cos(theta) from z = cos(theta) = 't' to z = 1, 0 otherwise. Basically a cos-weighted spotlight.
ZH GatedCosZH(float t)
{
    float t2 = t * t;
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;
    float t6 = t3 * t3;

    ZH zh;
    zh.c[0] = kTwoPi * kZH_Y_0 * 0.5f * (1.0f - t2);
    zh.c[1] = kTwoPi * kZH_Y_1 * (1.0f - t3) / 3.0f;
#if ZH_COEFFS >= 3
    zh.c[2] = kTwoPi * kZH_Y_2 * (1.0f + 2.0f * t2 - 3.0f * t4) / 4.0f;
#endif
#if ZH_COEFFS >= 4
    zh.c[3] = kTwoPi * kZH_Y_3 * (t3 - t5);
#endif
#if ZH_COEFFS >= 5
    zh.c[4] = kTwoPi * kZH_Y_4 * -(1.0f + 9.0f * t2 - 45.0f * t4 + 35.0f * t6) / 6.0f;
#endif

    return zh;
}

ZH HGPhaseZH(float g, float weight)
{
    // The Heyney-Greenstein phase function turns out to be a
    // simple power series in g when decomposed into zonal harmonics.
    ZH zh;
    for (int i = 0; i < ZH_COEFFS; i++)
    {
        zh.c[i] = sqrt(kFourPi * float(i * 2 + 1)) * weight;
        weight *= g;
    }
    return zh;
}

ZH NormalizedHGPhaseZH(float g, float weight)
{
    // normalize so that max of function is always 1.
    weight *= (1.0 - abs(g)) * (1.0 - abs(g)) / (1.0 + abs(g));

    ZH zh;
    for (int i = 0; i < ZH_COEFFS; i++)
    {
        zh.c[i] = sqrt(kFourPi * float(i * 2 + 1)) * weight;
        weight *= g;
    }
    return zh;
}

float WindowScale(int n, float gamma)
{
    float nt = float(n) * (float(n) + 1.0);
    return 1.0 / (1.0 + gamma * nt * nt);
}

void ApplyWindowingZH(float gamma, inout float zh_coeffs[ZH_COEFFS])
{
    for (int i = 0; i < ZH_COEFFS; i++)
        zh_coeffs[i] *= WindowScale(i, gamma);
}




////////////////////////////////////////////////////////////////////////////////
// ZH/SH Rotation utilities
////////////////////////////////////////////////////////////////////////////////

CONST(float kSqrt02_01  = sqrt( 2.0 /  1.0));
CONST(float kSqrt01_02  = sqrt( 1.0 /  2.0));
CONST(float kSqrt03_02  = sqrt( 3.0 /  2.0));
CONST(float kSqrt01_03  = sqrt( 1.0 /  3.0));
CONST(float kSqrt02_03  = sqrt( 2.0 /  3.0));
CONST(float kSqrt04_03  = sqrt( 4.0 /  3.0));
CONST(float kSqrt01_04  = sqrt( 1.0 /  4.0));
CONST(float kSqrt03_04  = sqrt( 3.0 /  4.0));
CONST(float kSqrt05_04  = sqrt( 5.0 /  4.0));
CONST(float kSqrt01_05  = sqrt( 1.0 /  5.0));
CONST(float kSqrt02_05  = sqrt( 2.0 /  5.0));
CONST(float kSqrt03_05  = sqrt( 3.0 /  5.0));
CONST(float kSqrt04_05  = sqrt( 4.0 /  5.0));
CONST(float kSqrt06_05  = sqrt( 6.0 /  5.0));
CONST(float kSqrt08_05  = sqrt( 8.0 /  5.0));
CONST(float kSqrt09_05  = sqrt( 9.0 /  5.0));
CONST(float kSqrt01_06  = sqrt( 1.0 /  6.0));
CONST(float kSqrt05_06  = sqrt( 5.0 /  6.0));
CONST(float kSqrt07_06  = sqrt( 7.0 /  6.0));
CONST(float kSqrt02_07  = sqrt(02.0 /  7.0));
CONST(float kSqrt06_07  = sqrt( 6.0 /  7.0));
CONST(float kSqrt10_07  = sqrt(10.0 /  7.0));
CONST(float kSqrt12_07  = sqrt(12.0 /  7.0));
CONST(float kSqrt15_07  = sqrt(15.0 /  7.0));
CONST(float kSqrt16_07  = sqrt(16.0 /  7.0));
CONST(float kSqrt01_08  = sqrt( 1.0 /  8.0));
CONST(float kSqrt03_08  = sqrt( 3.0 /  8.0));
CONST(float kSqrt05_08  = sqrt( 5.0 /  8.0));
CONST(float kSqrt07_08  = sqrt( 7.0 /  8.0));
CONST(float kSqrt09_08  = sqrt( 9.0 /  8.0));
CONST(float kSqrt05_09  = sqrt( 5.0 /  9.0));
CONST(float kSqrt08_09  = sqrt( 8.0 /  9.0));
CONST(float kSqrt01_10  = sqrt( 1.0 / 10.0));
CONST(float kSqrt03_10  = sqrt( 3.0 / 10.0));
CONST(float kSqrt07_10  = sqrt( 7.0 / 10.0));
CONST(float kSqrt01_12  = sqrt( 1.0 / 12.0));
CONST(float kSqrt07_12  = sqrt( 7.0 / 12.0));
CONST(float kSqrt01_14  = sqrt( 1.0 / 14.0));
CONST(float kSqrt03_14  = sqrt( 3.0 / 14.0));
CONST(float kSqrt15_14  = sqrt(15.0 / 14.0));
CONST(float kSqrt04_15  = sqrt( 4.0 / 15.0));
CONST(float kSqrt07_15  = sqrt( 7.0 / 10.0));
CONST(float kSqrt14_15  = sqrt(14.0 / 15.0));
CONST(float kSqrt16_15  = sqrt(16.0 / 15.0));
CONST(float kSqrt01_16  = sqrt( 1.0 / 16.0));
CONST(float kSqrt03_16  = sqrt( 3.0 / 16.0));
CONST(float kSqrt07_16  = sqrt( 7.0 / 16.0));
CONST(float kSqrt15_16  = sqrt(15.0 / 16.0));
CONST(float kSqrt01_18  = sqrt( 1.0 / 18.0));
CONST(float kSqrt01_24  = sqrt( 1.0 / 24.0));
CONST(float kSqrt03_25  = sqrt( 3.0 / 25.0));
CONST(float kSqrt14_25  = sqrt(14.0 / 25.0));
CONST(float kSqrt15_25  = sqrt(15.0 / 25.0));
CONST(float kSqrt18_25  = sqrt(18.0 / 25.0));
CONST(float kSqrt03_28  = sqrt( 3.0 / 28.0));
CONST(float kSqrt05_28  = sqrt( 5.0 / 28.0));
CONST(float kSqrt01_30  = sqrt( 1.0 / 30.0));
CONST(float kSqrt01_32  = sqrt( 1.0 / 32.0));
CONST(float kSqrt03_32  = sqrt( 3.0 / 32.0));
CONST(float kSqrt15_32  = sqrt(15.0 / 32.0));
CONST(float kSqrt21_32  = sqrt(21.0 / 32.0));
CONST(float kSqrt01_50  = sqrt( 1.0 / 50.0));
CONST(float kSqrt03_50  = sqrt( 3.0 / 50.0));
CONST(float kSqrt21_50  = sqrt(21.0 / 50.0));
CONST(float kSqrt15_56  = sqrt(15.0 / 56.0));
CONST(float kSqrt01_60  = sqrt( 1.0 / 60.0));
CONST(float kSqrt01_112 = sqrt( 1.0 / 112.0));
CONST(float kSqrt03_112 = sqrt( 3.0 / 112.0));
CONST(float kSqrt15_112 = sqrt(15.0 / 112.0));

void RotateZHToSHAdd(vec3 dir, SH_TYPE c, int n, const float zcoeffs[SH_BANDS], inout SH_TYPE coeffs[SH_COEFFS])
{
    coeffs[0] += c * zcoeffs[0];

    if (n < 2)
        return;

    vec3 sh1 = dir.yzx;
    coeffs[1] += c * zcoeffs[1] * sh1[0];
    coeffs[2] += c * zcoeffs[1] * sh1[1];
    coeffs[3] += c * zcoeffs[1] * sh1[2];

#if SH_BANDS >= 3
    if (n < 3)
        return;

    float sh2[5];

    sh2[0] =  kSqrt03_04 * (sh1[2] * sh1[0] + sh1[0] * sh1[2]);
    sh2[1] =  kSqrt03_04 * (sh1[1] * sh1[0] + sh1[0] * sh1[1]);
    sh2[2] = -kSqrt01_04 * (sh1[2] * sh1[2] + sh1[0] * sh1[0]) + sh1[1] * sh1[1];
    sh2[3] =  kSqrt03_04 * (sh1[1] * sh1[2] + sh1[2] * sh1[1]);
    sh2[4] =  kSqrt03_04 * (sh1[2] * sh1[2] - sh1[0] * sh1[0]);

    for (int i = 0; i < 5; i++)
        coeffs[i + 4] += c * zcoeffs[2] * sh2[i];
#endif
#if SH_BANDS >= 4
    if (n < 4)
        return;

    float sh3[7];

    sh3[0] =  kSqrt05_06 * (sh1[2] * sh2[0] + sh1[0] * sh2[4]);
    sh3[1] =  kSqrt05_09 * (sh1[1] * sh2[0] + (sh1[2] * sh2[1] + sh1[0] * sh2[3]));
    sh3[2] =  kSqrt08_09 *  sh1[1] * sh2[1] + kSqrt02_03 *  sh1[0] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[0] - sh1[0] * sh2[4]);
    sh3[3] =                sh1[1] * sh2[2] - kSqrt01_03 * (sh1[2] * sh2[3] + sh1[0] * sh2[1]);
    sh3[4] =  kSqrt08_09 *  sh1[1] * sh2[3] + kSqrt02_03 *  sh1[2] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[4] + sh1[0] * sh2[0]);
    sh3[5] =  kSqrt05_09 * (sh1[1] * sh2[4] + (sh1[2] * sh2[3] - sh1[0] * sh2[1]));
    sh3[6] =  kSqrt05_06 * (sh1[2] * sh2[4] - sh1[0] * sh2[0]);

    for (int i = 0; i < 7; i++)
        coeffs[i + 9] += c * zcoeffs[3] * sh3[i];
#endif
#if SH_BANDS >= 5
    if (n < 5)
        return;

    float sh4[9];

    sh4[0] =  kSqrt07_08 * (sh1[2] * sh3[0] + sh1[0] * sh3[6]);
    sh4[1] =  kSqrt07_16 *  sh1[1] * sh3[0] + kSqrt21_32 * (sh1[2] * sh3[1] + sh1[0] * sh3[5]);
    sh4[2] =  kSqrt03_04 *  sh1[1] * sh3[1] + kSqrt15_32 * (sh1[2] * sh3[2] + sh1[0] * sh3[4]) - kSqrt01_32 * (sh1[2] * sh3[0] - sh1[0] * sh3[6]);
    sh4[3] =  kSqrt15_16 *  sh1[1] * sh3[2] + kSqrt05_08 *  sh1[0] * sh3[3] - kSqrt03_32 * (sh1[2] * sh3[1] - sh1[0] * sh3[5]);
    sh4[4] =                sh1[1] * sh3[3] - kSqrt03_08 * (sh1[2] * sh3[4] + sh1[0] * sh3[2]);
    sh4[5] =  kSqrt15_16 *  sh1[1] * sh3[4] + kSqrt05_08 *  sh1[2] * sh3[3] - kSqrt03_32 * (sh1[2] * sh3[5] + sh1[0] * sh3[1]);
    sh4[6] =  kSqrt03_04 *  sh1[1] * sh3[5] + kSqrt15_32 * (sh1[2] * sh3[4] - sh1[0] * sh3[2]) - kSqrt01_32 * (sh1[2] * sh3[6] + sh1[0] * sh3[0]);
    sh4[7] =  kSqrt07_16 *  sh1[1] * sh3[6] + kSqrt21_32 * (sh1[2] * sh3[5] - sh1[0] * sh3[1]);
    sh4[8] =  kSqrt07_08 * (sh1[2] * sh3[6] - sh1[0] * sh3[0]);

    for (int i = 0; i < 9; i++)
        coeffs[i + 16] += c * zcoeffs[4] * sh4[i];
#endif
}

void RotateZHToSH(vec3 dir, SH_TYPE c, int n, const float zcoeffs[SH_BANDS], out SH_TYPE coeffs[SH_COEFFS])
{
    coeffs[0] = c * zcoeffs[0];

    if (n < 2)
        return;

    vec3 sh1 = dir.yzx;
    coeffs[1] = c * zcoeffs[1] * sh1[0];
    coeffs[2] = c * zcoeffs[1] * sh1[1];
    coeffs[3] = c * zcoeffs[1] * sh1[2];

#if SH_BANDS >= 3
    if (n < 3)
        return;

    float sh2[5];

    sh2[0] =  kSqrt03_04 * (sh1[2] * sh1[0] + sh1[0] * sh1[2]);
    sh2[1] =  kSqrt03_04 * (sh1[1] * sh1[0] + sh1[0] * sh1[1]);
    sh2[2] = -kSqrt01_04 * (sh1[2] * sh1[2] + sh1[0] * sh1[0]) + sh1[1] * sh1[1];
    sh2[3] =  kSqrt03_04 * (sh1[1] * sh1[2] + sh1[2] * sh1[1]);
    sh2[4] =  kSqrt03_04 * (sh1[2] * sh1[2] - sh1[0] * sh1[0]);

    for (int i = 0; i < 5; i++)
        coeffs[i + 4] = c * zcoeffs[2] * sh2[i];
#endif
#if SH_BANDS >= 4
    if (n < 4)
        return;

    float sh3[7];

    sh3[0] =  kSqrt05_06 * (sh1[2] * sh2[0] + sh1[0] * sh2[4]);
    sh3[1] =  kSqrt05_09 * (sh1[1] * sh2[0] + (sh1[2] * sh2[1] + sh1[0] * sh2[3]));
    sh3[2] =  kSqrt08_09 *  sh1[1] * sh2[1] + kSqrt02_03 *  sh1[0] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[0] - sh1[0] * sh2[4]);
    sh3[3] =                sh1[1] * sh2[2] - kSqrt01_03 * (sh1[2] * sh2[3] + sh1[0] * sh2[1]);
    sh3[4] =  kSqrt08_09 *  sh1[1] * sh2[3] + kSqrt02_03 *  sh1[2] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[4] + sh1[0] * sh2[0]);
    sh3[5] =  kSqrt05_09 * (sh1[1] * sh2[4] + (sh1[2] * sh2[3] - sh1[0] * sh2[1]));
    sh3[6] =  kSqrt05_06 * (sh1[2] * sh2[4] - sh1[0] * sh2[0]);

    for (int i = 0; i < 7; i++)
        coeffs[i + 9] = c * zcoeffs[3] * sh3[i];
#endif
#if SH_BANDS >= 5
    if (n < 5)
        return;

    float sh4[9];

    sh4[0] =  kSqrt07_08 * (sh1[2] * sh3[0] + sh1[0] * sh3[6]);
    sh4[1] =  kSqrt07_16 *  sh1[1] * sh3[0] + kSqrt21_32 * (sh1[2] * sh3[1] + sh1[0] * sh3[5]);
    sh4[2] =  kSqrt03_04 *  sh1[1] * sh3[1] + kSqrt15_32 * (sh1[2] * sh3[2] + sh1[0] * sh3[4]) - kSqrt01_32 * (sh1[2] * sh3[0] - sh1[0] * sh3[6]);
    sh4[3] =  kSqrt15_16 *  sh1[1] * sh3[2] + kSqrt05_08 *  sh1[0] * sh3[3] - kSqrt03_32 * (sh1[2] * sh3[1] - sh1[0] * sh3[5]);
    sh4[4] =                sh1[1] * sh3[3] - kSqrt03_08 * (sh1[2] * sh3[4] + sh1[0] * sh3[2]);
    sh4[5] =  kSqrt15_16 *  sh1[1] * sh3[4] + kSqrt05_08 *  sh1[2] * sh3[3] - kSqrt03_32 * (sh1[2] * sh3[5] + sh1[0] * sh3[1]);
    sh4[6] =  kSqrt03_04 *  sh1[1] * sh3[5] + kSqrt15_32 * (sh1[2] * sh3[4] - sh1[0] * sh3[2]) - kSqrt01_32 * (sh1[2] * sh3[6] + sh1[0] * sh3[0]);
    sh4[7] =  kSqrt07_16 *  sh1[1] * sh3[6] + kSqrt21_32 * (sh1[2] * sh3[5] - sh1[0] * sh3[1]);
    sh4[8] =  kSqrt07_08 * (sh1[2] * sh3[6] - sh1[0] * sh3[0]);

    for (int i = 0; i < 9; i++)
        coeffs[i + 16] = c * zcoeffs[4] * sh4[i];
#endif
}

void RotateZHToSH(vec3 dir, int n, const in SH_TYPE zcoeffs[SH_BANDS], out SH_TYPE coeffs[SH_COEFFS])
{
    coeffs[0] = zcoeffs[0];

    if (n < 2)
        return;

    vec3 sh1 = dir.yzx;
    coeffs[1] = zcoeffs[1] * sh1[0];
    coeffs[2] = zcoeffs[1] * sh1[1];
    coeffs[3] = zcoeffs[1] * sh1[2];

#if SH_BANDS >= 3
    if (n < 3)
        return;

    float sh2[5];

    sh2[0] =  kSqrt03_04 * (sh1[2] * sh1[0] + sh1[0] * sh1[2]);
    sh2[1] =  kSqrt03_04 * (sh1[1] * sh1[0] + sh1[0] * sh1[1]);
    sh2[2] = -kSqrt01_04 * (sh1[2] * sh1[2] + sh1[0] * sh1[0]) + sh1[1] * sh1[1];
    sh2[3] =  kSqrt03_04 * (sh1[1] * sh1[2] + sh1[2] * sh1[1]);
    sh2[4] =  kSqrt03_04 * (sh1[2] * sh1[2] - sh1[0] * sh1[0]);

    for (int i = 0; i < 5; i++)
        coeffs[i + 4] = zcoeffs[2] * sh2[i];
#endif
#if SH_BANDS >= 4
    if (n < 4)
        return;

    float sh3[7];

    sh3[0] =  kSqrt05_06 * (sh1[2] * sh2[0] + sh1[0] * sh2[4]);
    sh3[1] =  kSqrt05_09 * (sh1[1] * sh2[0] + (sh1[2] * sh2[1] + sh1[0] * sh2[3]));
    sh3[2] =  kSqrt08_09 *  sh1[1] * sh2[1] + kSqrt02_03 *  sh1[0] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[0] - sh1[0] * sh2[4]);
    sh3[3] =                sh1[1] * sh2[2] - kSqrt01_03 * (sh1[2] * sh2[3] + sh1[0] * sh2[1]);
    sh3[4] =  kSqrt08_09 *  sh1[1] * sh2[3] + kSqrt02_03 *  sh1[2] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[4] + sh1[0] * sh2[0]);
    sh3[5] =  kSqrt05_09 * (sh1[1] * sh2[4] + (sh1[2] * sh2[3] - sh1[0] * sh2[1]));
    sh3[6] =  kSqrt05_06 * (sh1[2] * sh2[4] - sh1[0] * sh2[0]);

    for (int i = 0; i < 7; i++)
        coeffs[i + 9] = zcoeffs[3] * sh3[i];
#endif
#if SH_BANDS >= 5
    if (n < 5)
        return;

    float sh4[9];

    sh4[0] =  kSqrt07_08 * (sh1[2] * sh3[0] + sh1[0] * sh3[6]);
    sh4[1] =  kSqrt07_16 *  sh1[1] * sh3[0] + kSqrt21_32 * (sh1[2] * sh3[1] + sh1[0] * sh3[5]);
    sh4[2] =  kSqrt03_04 *  sh1[1] * sh3[1] + kSqrt15_32 * (sh1[2] * sh3[2] + sh1[0] * sh3[4]) - kSqrt01_32 * (sh1[2] * sh3[0] - sh1[0] * sh3[6]);
    sh4[3] =  kSqrt15_16 *  sh1[1] * sh3[2] + kSqrt05_08 *  sh1[0] * sh3[3] - kSqrt03_32 * (sh1[2] * sh3[1] - sh1[0] * sh3[5]);
    sh4[4] =                sh1[1] * sh3[3] - kSqrt03_08 * (sh1[2] * sh3[4] + sh1[0] * sh3[2]);
    sh4[5] =  kSqrt15_16 *  sh1[1] * sh3[4] + kSqrt05_08 *  sh1[2] * sh3[3] - kSqrt03_32 * (sh1[2] * sh3[5] + sh1[0] * sh3[1]);
    sh4[6] =  kSqrt03_04 *  sh1[1] * sh3[5] + kSqrt15_32 * (sh1[2] * sh3[4] - sh1[0] * sh3[2]) - kSqrt01_32 * (sh1[2] * sh3[6] + sh1[0] * sh3[0]);
    sh4[7] =  kSqrt07_16 *  sh1[1] * sh3[6] + kSqrt21_32 * (sh1[2] * sh3[5] - sh1[0] * sh3[1]);
    sh4[8] =  kSqrt07_08 * (sh1[2] * sh3[6] - sh1[0] * sh3[0]);

    for (int i = 0; i < 9; i++)
        coeffs[i + 16] = zcoeffs[4] * sh4[i];
#endif
}

void RotateSHAboutZ(float theta, int n, inout SH_TYPE coeffs[SH_COEFFS])
{
    if (n < 2)
        return;

    float rc[2 * SH_BANDS - 2];

    rc[0] = cos(theta);
    rc[1] = sin(theta);

    SH_TYPE tmp = coeffs[1] * rc[0] + coeffs[3] * rc[1];
    coeffs[3]   = coeffs[3] * rc[0] - coeffs[1] * rc[1];
    coeffs[1]   = tmp;

    for (int band = 2; band < n; band++)
    {
        // on to the next band
        int ri = 2 * (band - 1);
        int ci = band * band;

        // find rot coeffs for outermost band coeffs
        rc[ri + 0] = rc[0] * rc[ri - 2] - rc[1] * rc[ri - 1];
        rc[ri + 1] = rc[0] * rc[ri - 1] + rc[1] * rc[ri - 2];

        // rotate outer bands
        for (int i = 0; i < band; i++)
        {
            int j = 2 * band - i;
            int k = 2 * (band - 1) - 2 * i;

            SH_TYPE tmp    = coeffs[ci + i] * rc[k] + coeffs[ci + j] * rc[k + 1];
            coeffs[ci + j] = coeffs[ci + j] * rc[k] - coeffs[ci + i] * rc[k + 1];
            coeffs[ci + i] = tmp;
        }
    }
}

void RotateSHAboutZ(float theta, inout SH_TYPE coeffs[SH_COEFFS])
{
    RotateSHAboutZ(theta, SH_BANDS, coeffs);
}

// General SH rotation

SH_TYPE dp3(SH_TYPE c[SH_COEFFS], float s[3])
{
    return c[1] * s[0] + c[2] * s[1] + c[3] * s[2];
}
#if SH_BANDS >= 3
SH_TYPE dp5(SH_TYPE c[SH_COEFFS], float s[5])
{
    return c[4] * s[0] + c[5] * s[1] + c[6] * s[2] + c[7] * s[3] + c[8] * s[4];
}
#endif
#if SH_BANDS >= 4
SH_TYPE dp7(SH_TYPE c[SH_COEFFS], float s[7])
{
    return c[9] * s[0] + c[10] * s[1] + c[11] * s[2] + c[12] * s[3] + c[13] * s[4] + c[14] * s[5] + c[15] * s[6];
}
#endif
#if SH_BANDS >= 5
SH_TYPE dp9(SH_TYPE c[SH_COEFFS], float s[9])
{
    SH_TYPE result = c[16] * s[0];
    for (int i = 1; i < 9; i++)
        result += c[16 + i] * s[i];
    return result;
}
#endif

void RotateSH(mat3 orient, int n, SH_TYPE coeffsIn[SH_COEFFS], out SH_TYPE coeffs[SH_COEFFS])
{
    int dst = 0;

    coeffs[dst++] = coeffsIn[0];
    if (n < 2)
        return;

#ifdef BGFX_SHADER_LANGUAGE_GLSL
    float sh10[3] = { orient[1][1], orient[2][1], orient[0][1] };
    float sh11[3] = { orient[1][2], orient[2][2], orient[0][2] };
    float sh12[3] = { orient[1][0], orient[2][0], orient[0][0] };
#else
    float sh10[3] = float[3](orient[1][1], orient[2][1], orient[0][1]);
    float sh11[3] = float[3](orient[1][2], orient[2][2], orient[0][2]);
    float sh12[3] = float[3](orient[1][0], orient[2][0], orient[0][0]);
#endif

    coeffs[dst++] = dp3(coeffsIn, sh10);
    coeffs[dst++] = dp3(coeffsIn, sh11);
    coeffs[dst++] = dp3(coeffsIn, sh12);

#if SH_BANDS >= 3
    if (n < 3)
        return;

    float sh20[5];
    sh20[0] = kSqrt01_04 * ((sh12[2] * sh10[0] + sh12[0] * sh10[2]) + (sh10[2] * sh12[0] + sh10[0] * sh12[2]));
    sh20[1] =               (sh12[1] * sh10[0] + sh10[1] * sh12[0]);
    sh20[2] = kSqrt03_04 *  (sh12[1] * sh10[1] + sh10[1] * sh12[1]);
    sh20[3] =               (sh12[1] * sh10[2] + sh10[1] * sh12[2]);
    sh20[4] = kSqrt01_04 * ((sh12[2] * sh10[2] - sh12[0] * sh10[0]) + (sh10[2] * sh12[2] - sh10[0] * sh12[0]));

    coeffs[dst++] = dp5(coeffsIn, sh20);

    float sh21[5];
    sh21[0] = kSqrt01_04 * ((sh11[2] * sh10[0] + sh11[0] * sh10[2]) + (sh10[2] * sh11[0] + sh10[0] * sh11[2]));
    sh21[1] =                sh11[1] * sh10[0] + sh10[1] * sh11[0];
    sh21[2] = kSqrt03_04 *  (sh11[1] * sh10[1] + sh10[1] * sh11[1]);
    sh21[3] =                sh11[1] * sh10[2] + sh10[1] * sh11[2];
    sh21[4] = kSqrt01_04 * ((sh11[2] * sh10[2] - sh11[0] * sh10[0]) + (sh10[2] * sh11[2] - sh10[0] * sh11[0]));

    coeffs[dst++] = dp5(coeffsIn, sh21);

    float sh22[5];
    sh22[0] = kSqrt01_03 * (sh11[2] * sh11[0] + sh11[0] * sh11[2]) - kSqrt01_12 * ((sh12[2] * sh12[0] + sh12[0] * sh12[2]) + (sh10[2] * sh10[0] + sh10[0] * sh10[2]));
    sh22[1] = kSqrt04_03 *  sh11[1] * sh11[0]                      - kSqrt01_03 *  (sh12[1] * sh12[0] + sh10[1] * sh10[0]);
    sh22[2] =               sh11[1] * sh11[1]                      - kSqrt01_04 *  (sh12[1] * sh12[1] + sh10[1] * sh10[1]);
    sh22[3] = kSqrt04_03 *  sh11[1] * sh11[2]                      - kSqrt01_03 *  (sh12[1] * sh12[2] + sh10[1] * sh10[2]);
    sh22[4] = kSqrt01_03 * (sh11[2] * sh11[2] - sh11[0] * sh11[0]) - kSqrt01_12 * ((sh12[2] * sh12[2] - sh12[0] * sh12[0]) + (sh10[2] * sh10[2] - sh10[0] * sh10[0]));

    coeffs[dst++] = dp5(coeffsIn, sh22);

    float sh23[5];
    sh23[0] = kSqrt01_04 * ((sh11[2] * sh12[0] + sh11[0] * sh12[2]) + (sh12[2] * sh11[0] + sh12[0] * sh11[2]));
    sh23[1] =                sh11[1] * sh12[0] + sh12[1] * sh11[0];
    sh23[2] = kSqrt03_04 *  (sh11[1] * sh12[1] + sh12[1] * sh11[1]);
    sh23[3] =                sh11[1] * sh12[2] + sh12[1] * sh11[2];
    sh23[4] = kSqrt01_04 * ((sh11[2] * sh12[2] - sh11[0] * sh12[0]) + (sh12[2] * sh11[2] - sh12[0] * sh11[0]));

    coeffs[dst++] = dp5(coeffsIn, sh23);

    float sh24[5];
    sh24[0] = kSqrt01_04 * ((sh12[2] * sh12[0] + sh12[0] * sh12[2]) - (sh10[2] * sh10[0] + sh10[0] * sh10[2]));
    sh24[1] =               (sh12[1] * sh12[0] - sh10[1] * sh10[0]);
    sh24[2] = kSqrt03_04 *  (sh12[1] * sh12[1] - sh10[1] * sh10[1]);
    sh24[3] =               (sh12[1] * sh12[2] - sh10[1] * sh10[2]);
    sh24[4] = kSqrt01_04 * ((sh12[2] * sh12[2] - sh12[0] * sh12[0]) - (sh10[2] * sh10[2] - sh10[0] * sh10[0]));

    coeffs[dst++] = dp5(coeffsIn, sh24);
#endif

#if SH_BANDS >= 4
    if (n < 4)
        return;

    float sh30[7];
    sh30[0] = kSqrt01_04 * ((sh12[2] * sh20[0] + sh12[0] * sh20[4]) + (sh10[2] * sh24[0] + sh10[0] * sh24[4]));
    sh30[1] = kSqrt03_02 * (sh12[1] * sh20[0] + sh10[1] * sh24[0]);
    sh30[2] = kSqrt15_16 * (sh12[1] * sh20[1] + sh10[1] * sh24[1]);
    sh30[3] = kSqrt05_06 * (sh12[1] * sh20[2] + sh10[1] * sh24[2]);
    sh30[4] = kSqrt15_16 * (sh12[1] * sh20[3] + sh10[1] * sh24[3]);
    sh30[5] = kSqrt03_02 * (sh12[1] * sh20[4] + sh10[1] * sh24[4]);
    sh30[6] = kSqrt01_04 * ((sh12[2] * sh20[4] - sh12[0] * sh20[0]) + (sh10[2] * sh24[4] - sh10[0] * sh24[0]));

    coeffs[dst++] = dp7(coeffsIn, sh30);

    float sh31[7];
    sh31[0] = kSqrt01_06 * (sh11[2] * sh20[0] + sh11[0] * sh20[4]) + kSqrt01_06 * ((sh12[2] * sh21[0] + sh12[0] * sh21[4]) + (sh10[2] * sh23[0] + sh10[0] * sh23[4]));
    sh31[1] = sh11[1] * sh20[0] + (sh12[1] * sh21[0] + sh10[1] * sh23[0]);
    sh31[2] = kSqrt05_08 * sh11[1] * sh20[1] + kSqrt05_08 * (sh12[1] * sh21[1] + sh10[1] * sh23[1]);
    sh31[3] = kSqrt05_09 * sh11[1] * sh20[2] + kSqrt05_09 * (sh12[1] * sh21[2] + sh10[1] * sh23[2]);
    sh31[4] = kSqrt05_08 * sh11[1] * sh20[3] + kSqrt05_08 * (sh12[1] * sh21[3] + sh10[1] * sh23[3]);
    sh31[5] = sh11[1] * sh20[4] + (sh12[1] * sh21[4] + sh10[1] * sh23[4]);
    sh31[6] = kSqrt01_06 * (sh11[2] * sh20[4] - sh11[0] * sh20[0]) + kSqrt01_06 * ((sh12[2] * sh21[4] - sh12[0] * sh21[0]) + (sh10[2] * sh23[4] - sh10[0] * sh23[0]));

    coeffs[dst++] = dp7(coeffsIn, sh31);

    float sh32[7];
    sh32[0] = kSqrt04_15 * (sh11[2] * sh21[0] + sh11[0] * sh21[4]) + kSqrt01_05 * (sh10[2] * sh22[0] + sh10[0] * sh22[4]) - kSqrt01_60 * ((sh12[2] * sh20[0] + sh12[0] * sh20[4]) - (sh10[2] * sh24[0] + sh10[0] * sh24[4]));
    sh32[1] = kSqrt08_05 * sh11[1] * sh21[0] + kSqrt06_05 * sh10[1] * sh22[0] - kSqrt01_10 * (sh12[1] * sh20[0] - sh10[1] * sh24[0]);
    sh32[2] = sh11[1] * sh21[1] + kSqrt03_04 * sh10[1] * sh22[1] - kSqrt01_16 * (sh12[1] * sh20[1] - sh10[1] * sh24[1]);
    sh32[3] = kSqrt08_09 * sh11[1] * sh21[2] + kSqrt02_03 * sh10[1] * sh22[2] - kSqrt01_18 * (sh12[1] * sh20[2] - sh10[1] * sh24[2]);
    sh32[4] = sh11[1] * sh21[3] + kSqrt03_04 * sh10[1] * sh22[3] - kSqrt01_16 * (sh12[1] * sh20[3] - sh10[1] * sh24[3]);
    sh32[5] = kSqrt08_05 * sh11[1] * sh21[4] + kSqrt06_05 * sh10[1] * sh22[4] - kSqrt01_10 * (sh12[1] * sh20[4] - sh10[1] * sh24[4]);
    sh32[6] = kSqrt04_15 * (sh11[2] * sh21[4] - sh11[0] * sh21[0]) + kSqrt01_05 * (sh10[2] * sh22[4] - sh10[0] * sh22[0]) - kSqrt01_60 * ((sh12[2] * sh20[4] - sh12[0] * sh20[0]) - (sh10[2] * sh24[4] - sh10[0] * sh24[0]));

    coeffs[dst++] = dp7(coeffsIn, sh32);

    float sh33[7];
    sh33[0] = kSqrt03_10 * (sh11[2] * sh22[0] + sh11[0] * sh22[4]) - kSqrt01_10 * ((sh12[2] * sh23[0] + sh12[0] * sh23[4]) + (sh10[2] * sh21[0] + sh10[0] * sh21[4]));
    sh33[1] = kSqrt09_05 * sh11[1] * sh22[0] - kSqrt03_05 * (sh12[1] * sh23[0] + sh10[1] * sh21[0]);
    sh33[2] = kSqrt09_08 * sh11[1] * sh22[1] - kSqrt03_08 * (sh12[1] * sh23[1] + sh10[1] * sh21[1]);
    sh33[3] = sh11[1] * sh22[2] - kSqrt01_03 * (sh12[1] * sh23[2] + sh10[1] * sh21[2]);
    sh33[4] = kSqrt09_08 * sh11[1] * sh22[3] - kSqrt03_08 * (sh12[1] * sh23[3] + sh10[1] * sh21[3]);
    sh33[5] = kSqrt09_05 * sh11[1] * sh22[4] - kSqrt03_05 * (sh12[1] * sh23[4] + sh10[1] * sh21[4]);
    sh33[6] = kSqrt03_10 * (sh11[2] * sh22[4] - sh11[0] * sh22[0]) - kSqrt01_10 * ((sh12[2] * sh23[4] - sh12[0] * sh23[0]) + (sh10[2] * sh21[4] - sh10[0] * sh21[0]));

    coeffs[dst++] = dp7(coeffsIn, sh33);

    float sh34[7];
    sh34[0] = kSqrt04_15 * (sh11[2] * sh23[0] + sh11[0] * sh23[4]) + kSqrt01_05 * (sh12[2] * sh22[0] + sh12[0] * sh22[4]) - kSqrt01_60 * ((sh12[2] * sh24[0] + sh12[0] * sh24[4]) + (sh10[2] * sh20[0] + sh10[0] * sh20[4]));
    sh34[1] = kSqrt08_05 * sh11[1] * sh23[0] + kSqrt06_05 * sh12[1] * sh22[0] - kSqrt01_10 * (sh12[1] * sh24[0] + sh10[1] * sh20[0]);
    sh34[2] = sh11[1] * sh23[1] + kSqrt03_04 * sh12[1] * sh22[1] - kSqrt01_16 * (sh12[1] * sh24[1] + sh10[1] * sh20[1]);
    sh34[3] = kSqrt08_09 * sh11[1] * sh23[2] + kSqrt02_03 * sh12[1] * sh22[2] - kSqrt01_18 * (sh12[1] * sh24[2] + sh10[1] * sh20[2]);
    sh34[4] = sh11[1] * sh23[3] + kSqrt03_04 * sh12[1] * sh22[3] - kSqrt01_16 * (sh12[1] * sh24[3] + sh10[1] * sh20[3]);
    sh34[5] = kSqrt08_05 * sh11[1] * sh23[4] + kSqrt06_05 * sh12[1] * sh22[4] - kSqrt01_10 * (sh12[1] * sh24[4] + sh10[1] * sh20[4]);
    sh34[6] = kSqrt04_15 * (sh11[2] * sh23[4] - sh11[0] * sh23[0]) + kSqrt01_05 * (sh12[2] * sh22[4] - sh12[0] * sh22[0]) - kSqrt01_60 * ((sh12[2] * sh24[4] - sh12[0] * sh24[0]) + (sh10[2] * sh20[4] - sh10[0] * sh20[0]));

    coeffs[dst++] = dp7(coeffsIn, sh34);

    float sh35[7];
    sh35[0] = kSqrt01_06 * (sh11[2] * sh24[0] + sh11[0] * sh24[4]) + kSqrt01_06 * ((sh12[2] * sh23[0] + sh12[0] * sh23[4]) - (sh10[2] * sh21[0] + sh10[0] * sh21[4]));
    sh35[1] = sh11[1] * sh24[0] + (sh12[1] * sh23[0] - sh10[1] * sh21[0]);
    sh35[2] = kSqrt05_08 * sh11[1] * sh24[1] + kSqrt05_08 * (sh12[1] * sh23[1] - sh10[1] * sh21[1]);
    sh35[3] = kSqrt05_09 * sh11[1] * sh24[2] + kSqrt05_09 * (sh12[1] * sh23[2] - sh10[1] * sh21[2]);
    sh35[4] = kSqrt05_08 * sh11[1] * sh24[3] + kSqrt05_08 * (sh12[1] * sh23[3] - sh10[1] * sh21[3]);
    sh35[5] = sh11[1] * sh24[4] + (sh12[1] * sh23[4] - sh10[1] * sh21[4]);
    sh35[6] = kSqrt01_06 * (sh11[2] * sh24[4] - sh11[0] * sh24[0]) + kSqrt01_06 * ((sh12[2] * sh23[4] - sh12[0] * sh23[0]) - (sh10[2] * sh21[4] - sh10[0] * sh21[0]));

    coeffs[dst++] = dp7(coeffsIn, sh35);

    float sh36[7];
    sh36[0] = kSqrt01_04 * ((sh12[2] * sh24[0] + sh12[0] * sh24[4]) - (sh10[2] * sh20[0] + sh10[0] * sh20[4]));
    sh36[1] = kSqrt03_02 * (sh12[1] * sh24[0] - sh10[1] * sh20[0]);
    sh36[2] = kSqrt15_16 * (sh12[1] * sh24[1] - sh10[1] * sh20[1]);
    sh36[3] = kSqrt05_06 * (sh12[1] * sh24[2] - sh10[1] * sh20[2]);
    sh36[4] = kSqrt15_16 * (sh12[1] * sh24[3] - sh10[1] * sh20[3]);
    sh36[5] = kSqrt03_02 * (sh12[1] * sh24[4] - sh10[1] * sh20[4]);
    sh36[6] = kSqrt01_04 * ((sh12[2] * sh24[4] - sh12[0] * sh24[0]) - (sh10[2] * sh20[4] - sh10[0] * sh20[0]));

    coeffs[dst++] = dp7(coeffsIn, sh36);
#endif

#if SH_BANDS >= 5
    if (n < 5)
        return;

    float sh40[9];

    sh40[0] = kSqrt01_04 * ((sh12[2] * sh30[0] + sh12[0] * sh30[6]) + (sh10[2] * sh36[0] + sh10[0] * sh36[6]));
    sh40[1] = kSqrt02_01 *  (sh12[1] * sh30[0] + sh10[1] * sh36[0]);
    sh40[2] = kSqrt07_06 *  (sh12[1] * sh30[1] + sh10[1] * sh36[1]);
    sh40[3] = kSqrt14_15 *  (sh12[1] * sh30[2] + sh10[1] * sh36[2]);
    sh40[4] = kSqrt07_08 *  (sh12[1] * sh30[3] + sh10[1] * sh36[3]);
    sh40[5] = kSqrt14_15 *  (sh12[1] * sh30[4] + sh10[1] * sh36[4]);
    sh40[6] = kSqrt07_06 *  (sh12[1] * sh30[5] + sh10[1] * sh36[5]);
    sh40[7] = kSqrt02_01 *  (sh12[1] * sh30[6] + sh10[1] * sh36[6]);
    sh40[8] = kSqrt01_04 * ((sh12[2] * sh30[6] - sh12[0] * sh30[0]) + (sh10[2] * sh36[6] - sh10[0] * sh36[0]));

    coeffs[dst++] = dp9(coeffsIn, sh40);

    float sh41[9];
    sh41[0] = kSqrt01_08 * (sh11[2] * sh30[0] + sh11[0] * sh30[6]) + kSqrt03_16 * ((sh12[2] * sh31[0] + sh12[0] * sh31[6]) + (sh10[2] * sh35[0] + sh10[0] * sh35[6]));
    sh41[1] =               sh11[1] * sh30[0]                      + kSqrt03_02 *  (sh12[1] * sh31[0] + sh10[1] * sh35[0]);
    sh41[2] = kSqrt07_12 *  sh11[1] * sh30[1]                      + kSqrt07_08 *  (sh12[1] * sh31[1] + sh10[1] * sh35[1]);
    sh41[3] = kSqrt07_15 *  sh11[1] * sh30[2]                      + kSqrt07_10 *  (sh12[1] * sh31[2] + sh10[1] * sh35[2]);
    sh41[4] = kSqrt07_16 *  sh11[1] * sh30[3]                      + kSqrt21_32 *  (sh12[1] * sh31[3] + sh10[1] * sh35[3]);
    sh41[5] = kSqrt07_15 *  sh11[1] * sh30[4]                      + kSqrt07_10 *  (sh12[1] * sh31[4] + sh10[1] * sh35[4]);
    sh41[6] = kSqrt07_12 *  sh11[1] * sh30[5]                      + kSqrt07_08 *  (sh12[1] * sh31[5] + sh10[1] * sh35[5]);
    sh41[7] =               sh11[1] * sh30[6]                      + kSqrt03_02 *  (sh12[1] * sh31[6] + sh10[1] * sh35[6]);
    sh41[8] = kSqrt01_08 * (sh11[2] * sh30[6] - sh11[0] * sh30[0]) + kSqrt03_16 * ((sh12[2] * sh31[6] - sh12[0] * sh31[0]) + (sh10[2] * sh35[6] - sh10[0] * sh35[0]));

    coeffs[dst++] = dp9(coeffsIn, sh41);

    float sh42[9];
    sh42[0] = kSqrt03_14 * (sh11[2] * sh31[0] + sh11[0] * sh31[6]) + kSqrt15_112 * ((sh12[2] * sh32[0] + sh12[0] * sh32[6]) + (sh10[2] * sh34[0] + sh10[0] * sh34[6])) - kSqrt01_112 * ((sh12[2] * sh30[0] + sh12[0] * sh30[6]) - (sh10[2] * sh36[0] + sh10[0] * sh36[6]));
    sh42[1] = kSqrt12_07 *  sh11[1] * sh31[0]                      + kSqrt15_14  *  (sh12[1] * sh32[0] + sh10[1] * sh34[0]) - kSqrt01_14 * (sh12[1] * sh30[0] - sh10[1] * sh36[0]);
    sh42[2] =               sh11[1] * sh31[1]                      + kSqrt05_08  *  (sh12[1] * sh32[1] + sh10[1] * sh34[1]) - kSqrt01_24 * (sh12[1] * sh30[1] - sh10[1] * sh36[1]);
    sh42[3] = kSqrt04_05 *  sh11[1] * sh31[2]                      + kSqrt01_02  *  (sh12[1] * sh32[2] + sh10[1] * sh34[2]) - kSqrt01_30 * (sh12[1] * sh30[2] - sh10[1] * sh36[2]);
    sh42[4] = kSqrt03_04 *  sh11[1] * sh31[3]                      + kSqrt15_32  *  (sh12[1] * sh32[3] + sh10[1] * sh34[3]) - kSqrt01_32 * (sh12[1] * sh30[3] - sh10[1] * sh36[3]);
    sh42[5] = kSqrt04_05 *  sh11[1] * sh31[4]                      + kSqrt01_02  *  (sh12[1] * sh32[4] + sh10[1] * sh34[4]) - kSqrt01_30 * (sh12[1] * sh30[4] - sh10[1] * sh36[4]);
    sh42[6] =               sh11[1] * sh31[5]                      + kSqrt05_08  *  (sh12[1] * sh32[5] + sh10[1] * sh34[5]) - kSqrt01_24 * (sh12[1] * sh30[5] - sh10[1] * sh36[5]);
    sh42[7] = kSqrt12_07 *  sh11[1] * sh31[6]                      + kSqrt15_14  *  (sh12[1] * sh32[6] + sh10[1] * sh34[6]) - kSqrt01_14 * (sh12[1] * sh30[6] - sh10[1] * sh36[6]);
    sh42[8] = kSqrt03_14 * (sh11[2] * sh31[6] - sh11[0] * sh31[0]) + kSqrt15_112 * ((sh12[2] * sh32[6] - sh12[0] * sh32[0]) + (sh10[2] * sh34[6] - sh10[0] * sh34[0])) - kSqrt01_112 * ((sh12[2] * sh30[6] - sh12[0] * sh30[0]) - (sh10[2] * sh36[6] - sh10[0] * sh36[0]));

    coeffs[dst++] = dp9(coeffsIn, sh42);

    float sh43[9];
    sh43[0] = kSqrt15_56 * (sh11[2] * sh32[0] + sh11[0] * sh32[6]) + kSqrt05_28 * (sh10[2] * sh33[0] + sh10[0] * sh33[6]) - kSqrt03_112 * ((sh12[2] * sh31[0] + sh12[0] * sh31[6]) - (sh10[2] * sh35[0] + sh10[0] * sh35[6]));
    sh43[1] = kSqrt15_07 *  sh11[1] * sh32[0]                      + kSqrt10_07 *  sh10[1] * sh33[0] - kSqrt03_14 * (sh12[1] * sh31[0] - sh10[1] * sh35[0]);
    sh43[2] = kSqrt05_04 *  sh11[1] * sh32[1]                      + kSqrt05_06 *  sh10[1] * sh33[1] - kSqrt01_08 * (sh12[1] * sh31[1] - sh10[1] * sh35[1]);
    sh43[3] =               sh11[1] * sh32[2]                      + kSqrt02_03 *  sh10[1] * sh33[2] - kSqrt01_10 * (sh12[1] * sh31[2] - sh10[1] * sh35[2]);
    sh43[4] = kSqrt15_16 *  sh11[1] * sh32[3]                      + kSqrt05_08 *  sh10[1] * sh33[3] - kSqrt03_32 * (sh12[1] * sh31[3] - sh10[1] * sh35[3]);
    sh43[5] =               sh11[1] * sh32[4]                      + kSqrt02_03 *  sh10[1] * sh33[4] - kSqrt01_10 * (sh12[1] * sh31[4] - sh10[1] * sh35[4]);
    sh43[6] = kSqrt05_04 *  sh11[1] * sh32[5]                      + kSqrt05_06 *  sh10[1] * sh33[5] - kSqrt01_08 * (sh12[1] * sh31[5] - sh10[1] * sh35[5]);
    sh43[7] = kSqrt15_07 *  sh11[1] * sh32[6]                      + kSqrt10_07 *  sh10[1] * sh33[6] - kSqrt03_14 * (sh12[1] * sh31[6] - sh10[1] * sh35[6]);
    sh43[8] = kSqrt15_56 * (sh11[2] * sh32[6] - sh11[0] * sh32[0]) + kSqrt05_28 * (sh10[2] * sh33[6] - sh10[0] * sh33[0]) - kSqrt03_112 * ((sh12[2] * sh31[6] - sh12[0] * sh31[0]) - (sh10[2] * sh35[6] - sh10[0] * sh35[0]));

    coeffs[dst++] = dp9(coeffsIn, sh43);

    float sh44[9];
    sh44[0] = kSqrt02_07 * (sh11[2] * sh33[0] + sh11[0] * sh33[6]) - kSqrt03_28 * ((sh12[2] * sh34[0] + sh12[0] * sh34[6]) + (sh10[2] * sh32[0] + sh10[0] * sh32[6]));
    sh44[1] = kSqrt16_07 *  sh11[1] * sh33[0]                      - kSqrt06_07 *  (sh12[1] * sh34[0] + sh10[1] * sh32[0]);
    sh44[2] = kSqrt04_03 *  sh11[1] * sh33[1]                      - kSqrt01_02 *  (sh12[1] * sh34[1] + sh10[1] * sh32[1]);
    sh44[3] = kSqrt16_15 *  sh11[1] * sh33[2]                      - kSqrt02_05 *  (sh12[1] * sh34[2] + sh10[1] * sh32[2]);
    sh44[4] =               sh11[1] * sh33[3]                      - kSqrt03_08 *  (sh12[1] * sh34[3] + sh10[1] * sh32[3]);
    sh44[5] = kSqrt16_15 *  sh11[1] * sh33[4]                      - kSqrt02_05 *  (sh12[1] * sh34[4] + sh10[1] * sh32[4]);
    sh44[6] = kSqrt04_03 *  sh11[1] * sh33[5]                      - kSqrt01_02 *  (sh12[1] * sh34[5] + sh10[1] * sh32[5]);
    sh44[7] = kSqrt16_07 *  sh11[1] * sh33[6]                      - kSqrt06_07 *  (sh12[1] * sh34[6] + sh10[1] * sh32[6]);
    sh44[8] = kSqrt02_07 * (sh11[2] * sh33[6] - sh11[0] * sh33[0]) - kSqrt03_28 * ((sh12[2] * sh34[6] - sh12[0] * sh34[0]) + (sh10[2] * sh32[6] - sh10[0] * sh32[0]));

    coeffs[dst++] = dp9(coeffsIn, sh44);

    float sh45[9];
    sh45[0] = kSqrt15_56 * (sh11[2] * sh34[0] + sh11[0] * sh34[6]) + kSqrt05_28 * (sh12[2] * sh33[0] + sh12[0] * sh33[6]) - kSqrt03_112 * ((sh12[2] * sh35[0] + sh12[0] * sh35[6]) + (sh10[2] * sh31[0] + sh10[0] * sh31[6]));
    sh45[1] = kSqrt15_07 *  sh11[1] * sh34[0]                      + kSqrt10_07 *  sh12[1] * sh33[0] - kSqrt03_14 * (sh12[1] * sh35[0] + sh10[1] * sh31[0]);
    sh45[2] = kSqrt05_04 *  sh11[1] * sh34[1]                      + kSqrt05_06 *  sh12[1] * sh33[1] - kSqrt01_08 * (sh12[1] * sh35[1] + sh10[1] * sh31[1]);
    sh45[3] =               sh11[1] * sh34[2]                      + kSqrt02_03 *  sh12[1] * sh33[2] - kSqrt01_10 * (sh12[1] * sh35[2] + sh10[1] * sh31[2]);
    sh45[4] = kSqrt15_16 *  sh11[1] * sh34[3]                      + kSqrt05_08 *  sh12[1] * sh33[3] - kSqrt03_32 * (sh12[1] * sh35[3] + sh10[1] * sh31[3]);
    sh45[5] =               sh11[1] * sh34[4]                      + kSqrt02_03 *  sh12[1] * sh33[4] - kSqrt01_10 * (sh12[1] * sh35[4] + sh10[1] * sh31[4]);
    sh45[6] = kSqrt05_04 *  sh11[1] * sh34[5]                      + kSqrt05_06 *  sh12[1] * sh33[5] - kSqrt01_08 * (sh12[1] * sh35[5] + sh10[1] * sh31[5]);
    sh45[7] = kSqrt15_07 *  sh11[1] * sh34[6]                      + kSqrt10_07 *  sh12[1] * sh33[6] - kSqrt03_14 * (sh12[1] * sh35[6] + sh10[1] * sh31[6]);
    sh45[8] = kSqrt15_56 * (sh11[2] * sh34[6] - sh11[0] * sh34[0]) + kSqrt05_28 * (sh12[2] * sh33[6] - sh12[0] * sh33[0]) - kSqrt03_112 * ((sh12[2] * sh35[6] - sh12[0] * sh35[0]) + (sh10[2] * sh31[6] - sh10[0] * sh31[0]));

    coeffs[dst++] = dp9(coeffsIn, sh45);

    float sh46[9];
    sh46[0] = kSqrt03_14 * (sh11[2] * sh35[0] + sh11[0] * sh35[6]) + kSqrt15_112 * ((sh12[2] * sh34[0] + sh12[0] * sh34[6]) - (sh10[2] * sh32[0] + sh10[0] * sh32[6])) - kSqrt01_112 * ((sh12[2] * sh36[0] + sh12[0] * sh36[6]) + (sh10[2] * sh30[0] + sh10[0] * sh30[6]));
    sh46[1] = kSqrt12_07 *  sh11[1] * sh35[0]                      + kSqrt15_14 *  (sh12[1] * sh34[0] - sh10[1] * sh32[0]) - kSqrt01_14 * (sh12[1] * sh36[0] + sh10[1] * sh30[0]);
    sh46[2] =               sh11[1] * sh35[1]                      + kSqrt05_08 *  (sh12[1] * sh34[1] - sh10[1] * sh32[1]) - kSqrt01_24 * (sh12[1] * sh36[1] + sh10[1] * sh30[1]);
    sh46[3] = kSqrt04_05 *  sh11[1] * sh35[2]                      + kSqrt01_02 *  (sh12[1] * sh34[2] - sh10[1] * sh32[2]) - kSqrt01_30 * (sh12[1] * sh36[2] + sh10[1] * sh30[2]);
    sh46[4] = kSqrt03_04 *  sh11[1] * sh35[3]                      + kSqrt15_32 *  (sh12[1] * sh34[3] - sh10[1] * sh32[3]) - kSqrt01_32 * (sh12[1] * sh36[3] + sh10[1] * sh30[3]);
    sh46[5] = kSqrt04_05 *  sh11[1] * sh35[4]                      + kSqrt01_02 *  (sh12[1] * sh34[4] - sh10[1] * sh32[4]) - kSqrt01_30 * (sh12[1] * sh36[4] + sh10[1] * sh30[4]);
    sh46[6] =               sh11[1] * sh35[5]                      + kSqrt05_08 *  (sh12[1] * sh34[5] - sh10[1] * sh32[5]) - kSqrt01_24 * (sh12[1] * sh36[5] + sh10[1] * sh30[5]);
    sh46[7] = kSqrt12_07 *  sh11[1] * sh35[6]                      + kSqrt15_14 *  (sh12[1] * sh34[6] - sh10[1] * sh32[6]) - kSqrt01_14 * (sh12[1] * sh36[6] + sh10[1] * sh30[6]);
    sh46[8] = kSqrt03_14 * (sh11[2] * sh35[6] - sh11[0] * sh35[0]) + kSqrt15_112 * ((sh12[2] * sh34[6] - sh12[0] * sh34[0]) - (sh10[2] * sh32[6] - sh10[0] * sh32[0])) - kSqrt01_112 * ((sh12[2] * sh36[6] - sh12[0] * sh36[0]) + (sh10[2] * sh30[6] - sh10[0] * sh30[0]));

    coeffs[dst++] = dp9(coeffsIn, sh46);

    float sh47[9];
    sh47[0] = kSqrt01_08 * (sh11[2] * sh36[0] + sh11[0] * sh36[6]) + kSqrt03_16 * ((sh12[2] * sh35[0] + sh12[0] * sh35[6]) - (sh10[2] * sh31[0] + sh10[0] * sh31[6]));
    sh47[1] =               sh11[1] * sh36[0]                      + kSqrt03_02 *  (sh12[1] * sh35[0] - sh10[1] * sh31[0]);
    sh47[2] = kSqrt07_12 *  sh11[1] * sh36[1]                      + kSqrt07_08 *  (sh12[1] * sh35[1] - sh10[1] * sh31[1]);
    sh47[3] = kSqrt07_15 *  sh11[1] * sh36[2]                      + kSqrt07_10 *  (sh12[1] * sh35[2] - sh10[1] * sh31[2]);
    sh47[4] = kSqrt07_16 *  sh11[1] * sh36[3]                      + kSqrt21_32 *  (sh12[1] * sh35[3] - sh10[1] * sh31[3]);
    sh47[5] = kSqrt07_15 *  sh11[1] * sh36[4]                      + kSqrt07_10 *  (sh12[1] * sh35[4] - sh10[1] * sh31[4]);
    sh47[6] = kSqrt07_12 *  sh11[1] * sh36[5]                      + kSqrt07_08 *  (sh12[1] * sh35[5] - sh10[1] * sh31[5]);
    sh47[7] =               sh11[1] * sh36[6]                      + kSqrt03_02 *  (sh12[1] * sh35[6] - sh10[1] * sh31[6]);
    sh47[8] = kSqrt01_08 * (sh11[2] * sh36[6] - sh11[0] * sh36[0]) + kSqrt03_16 * ((sh12[2] * sh35[6] - sh12[0] * sh35[0]) - (sh10[2] * sh31[6] - sh10[0] * sh31[0]));

    coeffs[dst++] = dp9(coeffsIn, sh47);

    float sh48[9];
    sh48[0] = kSqrt01_04 * ((sh12[2] * sh36[0] + sh12[0] * sh36[6]) - (sh10[2] * sh30[0] + sh10[0] * sh30[6]));
    sh48[1] = kSqrt02_01 *  (sh12[1] * sh36[0] - sh10[1] * sh30[0]);
    sh48[2] = kSqrt07_06 *  (sh12[1] * sh36[1] - sh10[1] * sh30[1]);
    sh48[3] = kSqrt14_15 *  (sh12[1] * sh36[2] - sh10[1] * sh30[2]);
    sh48[4] = kSqrt07_08 *  (sh12[1] * sh36[3] - sh10[1] * sh30[3]);
    sh48[5] = kSqrt14_15 *  (sh12[1] * sh36[4] - sh10[1] * sh30[4]);
    sh48[6] = kSqrt07_06 *  (sh12[1] * sh36[5] - sh10[1] * sh30[5]);
    sh48[7] = kSqrt02_01 *  (sh12[1] * sh36[6] - sh10[1] * sh30[6]);
    sh48[8] = kSqrt01_04 * ((sh12[2] * sh36[6] - sh12[0] * sh36[0]) - (sh10[2] * sh30[6] - sh10[0] * sh30[0]));

    coeffs[dst++] = dp9(coeffsIn, sh48);
#endif
}

void RotateSH(mat3 orient, SH_TYPE srcCoeffs[SH_COEFFS], out SH_TYPE coeffs[SH_COEFFS])
{
    RotateSH(orient, SH_BANDS, srcCoeffs, coeffs);
}


////////////////////////////////////////////////////////////////////////////////
// Lighting utilities
////////////////////////////////////////////////////////////////////////////////

// simple lerp between diffuse and pass-through, which is the full lighting
// environment bandpass limited by the SH order.
void ApplyGlossyBRDF(float gloss, int n, inout SH_TYPE coeffs[SH_COEFFS])
{
    if (n < 2)
        return;

    if (gloss >= 1.0)
        return;

    float diff = 1.0 - gloss;

    float w = gloss + diff * kSH_A1;
    coeffs[1] *= w;
    coeffs[2] *= w;
    coeffs[3] *= w;

    if (n < 3)
        return;

    w = gloss + diff * kSH_A2;
    for (int i = 4; i < 9; i++)
        coeffs[i] *= w;

    if (n < 4)
        return;

    for (int i = 9; i < 16; i++)
        coeffs[i] *= gloss;

    if (n < 5)
        return;

    w = gloss + diff * kSH_A4;
    for (int i = 16; i < 25; i++)
        coeffs[i] *= w;

    if (n < 6)
        return;

    int nc = n * n;
    for (int i = 25; i < nc; i++)
        coeffs[i] *= gloss;
}

void ApplyGlossyBRDF(float gloss, inout SH_TYPE coeffs[SH_COEFFS])
{
    ApplyGlossyBRDF(gloss, SH_BANDS, coeffs);
}

// AddSphereLight: Uses ZH to represent a pseudo spherical area light
void AddSphereLight(vec3 Ld, SH_TYPE Li, float Lsize, inout SH_TYPE sh_coeffs[SH_COEFFS])
{
    float r2 = dot(Ld, Ld);

    float a = Lsize;
    float a2 = a * a;

    float t = sqrt(r2 / (a2 + r2));
    Ld *= inversesqrt(r2);

    RotateZHToSHAdd(Ld, Li, SH_BANDS, GatedSpotZH(t).c, sh_coeffs);
}

// Add in normalised Henyey-Greenstein phase function.
void AddHGPhase(float g, float w, SH_TYPE c, inout SH_TYPE sh_coeffs[SH_COEFFS])
{
    ZH zh = NormalizedHGPhaseZH(g, w);

    for (int i = 0; i < SH_BANDS; i++)
        sh_coeffs[i * (i + 1)] += zh.c[i] * c;
}

void AddHGPhase(float g, float w, SH_TYPE c, vec3 d, inout SH_TYPE sh_coeffs[SH_COEFFS])
{
    RotateZHToSHAdd(d, c, SH_BANDS, NormalizedHGPhaseZH(g, w).c, sh_coeffs);
}

#endif


#ifdef DEMO

//
// Quick shader toy-style demo
//

void mainImage(out vec4 fragColor, in vec2 pixel)
{
    vec2 p2 = pixel / iResolution.xy;


    p2 *= vec2(3.0, 2.0);
    vec2 tile = floor(p2);
    p2 -= tile;

    p2 = 2.0 * p2 - vec2(1.0);

    vec2 tile_size = iResolution.xy / vec2(3.0, 2.0);
    p2 *= tile_size.xy / min(tile_size.x, tile_size.y);  // squarify
    p2 /= 0.9;  // border

    float r = length(p2);

    if (r > 1.0)
    {
        fragColor = vec4(0.4, 0.4, 0.5, 1.0);
        return;
    }

    vec3 p = vec3(p2, sqrt(1.0 - r * r));

    vec4 sh_coeffs[SH_COEFFS];

    if (tile.x == 2.0)
    {
        vec3 Lsun = normalize(vec3(1.0, 1.0, 1.0));  // sun direction

        AddHGPhase(0.15, 0.9, vec4(0.0, 0.5, 1.0, 1.0)      , sh_coeffs);  // sky
        AddHGPhase(0.80, 1.2, vec4(1.0, 0.3, 0.3, 1.0), Lsun, sh_coeffs);  // sun halo
        AddHGPhase(0.92, 7.0, vec4(1.0, 1.0, 1.0, 1.0), Lsun, sh_coeffs);  // sun
    }
    else if (tile.x == 1.0)
    {
        const vec4 yellow_blue = vec4(0.5, 0.5, -0.5, 1.0);

        sh_coeffs[ 0] = vec4(0.8);
        sh_coeffs[ 2] = yellow_blue;
        sh_coeffs[ 4] = yellow_blue;
        sh_coeffs[10] = yellow_blue;
        sh_coeffs[21] = yellow_blue;
        sh_coeffs[17] = yellow_blue;
    }
    else
    {
        const vec4 red_green = vec4(-1.0, 1.0, 0.0, 1.0);

        sh_coeffs[ 1] = 1.0 * red_green;
        sh_coeffs[ 7] = 1.0 * red_green;
        sh_coeffs[10] = 1.0 * red_green;
        sh_coeffs[21] = 0.5 * red_green;
        sh_coeffs[24] = 1.0 * red_green;
    }
    if (tile.y == 0.0)
    {
        float t = iTime;

        if (tile.x == 0.0)
        {
            float c = cos(t);
            float s = sin(t);

            mat3 rot_y = mat3
            (
                 c , 0.0, -s,
                0.0, 1.0,  0.0,
                 s , 0.0,  c
            );

            vec4 sh_coeffs_src[SH_COEFFS] = sh_coeffs;
            RotateSH(rot_y, sh_coeffs_src, sh_coeffs);
        }
        else if (tile.x == 1.0)
        {
            ApplyGlossyBRDF(abs(2.0 * fract(t * 0.5) - 1.0), sh_coeffs);
        }
        else
            RotateSHAboutZ(t, sh_coeffs);
    }

    fragColor.rgb = SH(normalize(p), sh_coeffs);
    fragColor.a = 1.0;
}

#endif
