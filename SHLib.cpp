//------------------------------------------------------------------------------
// SH Utility Library
//------------------------------------------------------------------------------

#include "SHLib.h"

#include "ZHLib.h"  // only for WindowScale, kZH_Y, so could decouple.

using namespace SHL;
using namespace ZHL;

namespace
{
    const float kFourPi = 4 * vl_pi;

    // The normalization constants themselves are given by
    //   sqrtf(((2l + 1) (l - m)!)  / (4 pi (l + m)!))
    // However, we also roll in any overall multipliers from the
    // basis function itself. (E.g., the b.f. for Y_22 is (x^2 - y^2) / 2.)

    // normalization constants                            basis function:
    const float kSH_Y_00 = sqrtf( 1 / ( kFourPi));          // 1

    // 3, 3 -> 1, 1
    const float kSH_Y_10 = sqrtf( 3 / ( kFourPi));          // y
    const float kSH_Y_11 = sqrtf( 3 / ( kFourPi));          // z
    const float kSH_Y_12 = sqrtf( 3 / ( kFourPi));          // x

    // 15, 15, 5 -> 3, 3, 1
    const float kSH_Y_20 = sqrtf(15 / (     kFourPi));      // xy
    const float kSH_Y_21 = sqrtf(15 / (     kFourPi));      // yz
    const float kSH_Y_22 = sqrtf( 5 / ( 4 * kFourPi));      // (3z^2 - 1)
    const float kSH_Y_23 = sqrtf(15 / (     kFourPi));      // xz
    const float kSH_Y_24 = sqrtf(15 / ( 4 * kFourPi));      // (x^2 - y^2)

    // 35/2  105 21/2  7 -> 5/2, 15, 3/2, 1
    const float kSH_Y_30 = sqrtf( 35 / ( 8 * kFourPi));     // y (3x^2 - y^2)
    const float kSH_Y_31 = sqrtf(105 / ( 1 * kFourPi));     // z xy
    const float kSH_Y_32 = sqrtf( 21 / ( 8 * kFourPi));     // y (5z^2 - 1)
    const float kSH_Y_33 = sqrtf(  7 / ( 4 * kFourPi));     // z (5z^2 - 3)
    const float kSH_Y_34 = sqrtf( 21 / ( 8 * kFourPi));     // x (5z^2 - 1)
    const float kSH_Y_35 = sqrtf(105 / ( 4 * kFourPi));     // z (x^2 - y^2)
    const float kSH_Y_36 = sqrtf( 35 / ( 8 * kFourPi));     // x (x^2 - 3y^2)

    // 315, 315/2, 45, 45/2, 9 -> 35, 35/2, 5, 5/2, 1
    const float kSH_Y_40 = sqrtf(315 / ( 4 * kFourPi));     // (x^2 - y^2)  xy
    const float kSH_Y_41 = sqrtf(315 / ( 8 * kFourPi));     // (3x^2 - y^2) yz
    const float kSH_Y_42 = sqrtf( 45 / ( 4 * kFourPi));     // (7z^2 - 1)   xy
    const float kSH_Y_43 = sqrtf( 45 / ( 8 * kFourPi));     // (7z^2 - 3)   yz
    const float kSH_Y_44 = sqrtf(  9 / (64 * kFourPi));     // (3 - 30z^2 + 35z^4)
    const float kSH_Y_45 = sqrtf( 45 / ( 8 * kFourPi));     // (7z^2 - 3)   xz
    const float kSH_Y_46 = sqrtf( 45 / (16 * kFourPi));     // (7z^2 - 1)   (x^2 -y^2)
    const float kSH_Y_47 = sqrtf(315 / ( 8 * kFourPi));     // (x^2 - 3y^2) xz
    const float kSH_Y_48 = sqrtf(315 / (64 * kFourPi));     // (x^4 - 6 x^2 y^2 + y^4)

    // diffuse reflectance convolution weights.
    // Multiply bands by these constants to convert incident radiance to exit
    // radiance for a diffuse surface.
    //const float kSH_A0 = 1.0f;
    const float kSH_A1 = (2.0f / 3.0f);
    const float kSH_A2 = (1.0f / 4.0f);
    //const float kSH_A3 = 0.0f;
    const float kSH_A4 = -(1.0f / 24.0f);
}

namespace
{
    /*
        SH function:
            m > 0       sqrt(2) * K(l, m) * cos(|m| * phi) * P(l, |m|, cos(theta))
            m < 0       sqrt(2) * K(l, m) * sin(|m| * phi) * P(l, |m|, cos(theta))
            m = 0       K(l, 0) * P(l, 0, cos(theta))

        If r^ = (x, y, z), then

            z = cos(theta)
            x = cos(phi)
            y = sin(phi)

        sin(n phi) = sum(k=0..n) binom(k, n) cos(phi)^k sin(phi)^(n-k) sin(1/2(n-k)pi)
        cos(n phi) = sum(k=0..n) binom(k, n) cos(phi)^k sin(phi)^(n-k) cos(1/2(n-k)pi)
        binom(k, n) = n! / k!(n-k)!

        So for m = 2,
            sin(2 phi) = 1 y^2 0 + 2 xy 1 + 1 x^2 0 = 2xy
            cos(2 phi) = 1 y^2 1 + 2 xy 0 + 1 x^2 1 = y^2 + x^2

        and so on.

    */

    uint32_t kFactorialTable[] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600 };

    inline float fact(int n)
    {
        if (n < 13)
            return float(kFactorialTable[n]);

        float result = float(kFactorialTable[12]);
        while (n >= 13)
            result *= n--;

        return result;
    }

    inline float SH_K(int l, int pm)
    {
        if (pm == 0)    // factorials cancel out
            return sqrtf(float(2 * l + 1) / kFourPi);

        float top    = float(2 * l + 1) * fact(l - pm);
        float bottom = kFourPi          * fact(l + pm);

        return sqrtf(top / bottom);
    }

    float SH_P(int l, int m, float x)
    {
        float pmm = 1.0f;

        if (m > 0)
        {
            float tm1 = sqrtf((1.0f - x) * (1.0f + x));
            float fact = 1.0f;

            for (int i = 1; i <= m; i++)
            {
                pmm *= (-fact) * tm1;
                fact += 2.0f;
            }
        }

        if (l == m)
            return pmm;

        float pmmp1 = x * float(2 * m + 1) * pmm;

        if (l == m + 1)
            return pmmp1;

        float pll = 0.0f;

        for (int ll = m + 2; ll <= l; ll++)
        {
            pll = (float(2 * ll - 1) * x * pmmp1 - float(ll + m - 1) * pmm);
            pll /= float(ll - m);
            pmm = pmmp1;
            pmmp1 = pll;
        }

        return pll;
    }
}

float SHL::SH(int l, int m, float theta, float phi)
{
    const float kRoot2 = sqrtf(2.0f);

    int pm = abs(m);
    float term1 = SH_K(l, pm) * SH_P(l, pm, cosf(theta));

    if (m == 0)
        return term1;

    phi += vl_pi;

    if (m > 0)
        return term1 * kRoot2 * cosf(pm * phi);
    else
        return term1 * kRoot2 * sinf(pm * phi);
}

namespace
{
    // Chebyshev recurrence relations for sin(nx) / cos(nx) given sin(x) / cos(x)
    inline float ChebyshevN(int n, float c, float n1, float n2)
    {
        if (n == 0)
          return n1;

        for (int i = 1; i < n; i++)
        {
            float n0 = n1;
            n1 = n2;
            n2 = 2 * c * n1 - n0;
        }

        return n2;
    }

    inline float SinN(int n, float s, float c)
    {
        return ChebyshevN(n, c, 0.0f, s);
    }

    inline float CosN(int n, float c)
    {
        return ChebyshevN(n, c, 1.0f, c);
    }

    const float kRoot2 = sqrtf(2.0f);
}

float SHL::SH(int l, int m, Vec3f dir)   // dir must be normalised
{
    float cosTheta = dir[2];

    int pm = abs(m);
    float term1 = SH_K(l, pm) * SH_P(l, pm, cosTheta);

    if (m == 0)
        return term1;

    float rxy = sqrtf(1.0f - cosTheta * cosTheta);
    float cosPhi, sinPhi;

    rxy += 1e-16f;
    float invSinTheta = 1.0f / rxy; // epsilon to avoid 0 divide

    // minus signs recapitulate += pi from SH(polar).
    cosPhi = -dir[0] * invSinTheta;

    if (m > 0)
        return term1 * kRoot2 * CosN(pm, cosPhi);

    sinPhi = -dir[1] * invSinTheta;

    return term1 * kRoot2 * SinN(pm, sinPhi, cosPhi);
}

const float SHL::kSH_Y[] =
{
    sqrtf(  1 / (     kFourPi)),  // 1
    sqrtf(  3 / (     kFourPi)),  // y
    sqrtf(  3 / (     kFourPi)),  // z
    sqrtf(  3 / (     kFourPi)),  // x
    sqrtf( 15 / (     kFourPi)),  // xy
    sqrtf( 15 / (     kFourPi)),  // yz
    sqrtf(  5 / ( 4 * kFourPi)),  // (3z^2 - 1)
    sqrtf( 15 / (     kFourPi)),  // xz
    sqrtf( 15 / ( 4 * kFourPi)),  // (x^2 - y^2)
    sqrtf( 35 / ( 8 * kFourPi)),  // y (3x^2 - y^2)
    sqrtf(105 / ( 1 * kFourPi)),  // z xy
    sqrtf( 21 / ( 8 * kFourPi)),  // y (5z^2 - 1)
    sqrtf(  7 / ( 4 * kFourPi)),  // z (5z^2 - 3)
    sqrtf( 21 / ( 8 * kFourPi)),  // x (5z^2 - 1)
    sqrtf(105 / ( 4 * kFourPi)),  // z (x^2 - y^2)
    sqrtf( 35 / ( 8 * kFourPi)),  // x (x^2 - 3y^2)
    sqrtf(315 / ( 4 * kFourPi)),  // (x^2 - y^2)  xy
    sqrtf(315 / ( 8 * kFourPi)),  // (3x^2 - y^2) yz
    sqrtf( 45 / ( 4 * kFourPi)),  // (7z^2 - 1)   xy
    sqrtf( 45 / ( 8 * kFourPi)),  // (7z^2 - 3)   yz
    sqrtf(  9 / (64 * kFourPi)),  // (3 - 30z^2 + 35z^4)
    sqrtf( 45 / ( 8 * kFourPi)),  // (7z^2 - 3)   xz
    sqrtf( 45 / (16 * kFourPi)),  // (7z^2 - 1)   (x^2 -y^2)
    sqrtf(315 / ( 8 * kFourPi)),  // (x^2 - 3y^2) xz
    sqrtf(315 / (64 * kFourPi)),  // (x^4 - 6 x^2 y^2 + y^4)
};

const float* SHL::kSH_Y0 = kSH_Y + 0;
const float* SHL::kSH_Y1 = kSH_Y + 1;
const float* SHL::kSH_Y2 = kSH_Y + 4;
const float* SHL::kSH_Y3 = kSH_Y + 9;
const float* SHL::kSH_Y4 = kSH_Y + 16;


// Main SH utils

namespace
{
    template<typename T> inline T SampleSH(Vec3f v, int numBands, const T coeffs[])
    {
        VL_ASSERT(numBands >= 1);

        T s(coeffs[0] * kSH_Y_00);

        if (numBands < 2)
            return s;

        const float& x  = v[0];
        const float& y  = v[1];
        const float& z  = v[2];

        float  x2 = sqr(x);
        float  y2 = sqr(y);
        float  z2 = sqr(z);

        s += coeffs[1] * kSH_Y_10 * y;
        s += coeffs[2] * kSH_Y_11 * z;
        s += coeffs[3] * kSH_Y_12 * x;

        if (numBands < 3)
            return s;

        float  xy = x * y;
        float  yz = y * z;
        float  xz = x * z;

        s += coeffs[4] * kSH_Y_20 * xy;
        s += coeffs[5] * kSH_Y_21 * yz;
        s += coeffs[6] * kSH_Y_22 * (3.0f * z2 - 1.0f);
        s += coeffs[7] * kSH_Y_23 * xz;
        s += coeffs[8] * kSH_Y_24 * (x2 - y2);

        if (numBands < 4)
            return s;

        s += coeffs[ 9] * kSH_Y_30 * (3.0f * x2 - y2)   * y;
        s += coeffs[10] * kSH_Y_31 * xy                 * z;
        s += coeffs[11] * kSH_Y_32 * (5.0f * z2 - 1.0f) * y;
        s += coeffs[12] * kSH_Y_33 * (5.0f * z2 - 3.0f) * z;
        s += coeffs[13] * kSH_Y_34 * (5.0f * z2 - 1.0f) * x;
        s += coeffs[14] * kSH_Y_35 * (x2 - y2)          * z;
        s += coeffs[15] * kSH_Y_36 * (x2 - 3.0f * y2)   * x;

        if (numBands < 5)
            return s;

        s += coeffs[16] * kSH_Y_40 * (x2 - y2)          * xy;
        s += coeffs[17] * kSH_Y_41 * (3.0f * x2 - y2)   * yz;
        s += coeffs[18] * kSH_Y_42 * (7.0f * z2 - 1.0f) * xy;
        s += coeffs[19] * kSH_Y_43 * (7.0f * z2 - 3.0f) * yz;
        s += coeffs[20] * kSH_Y_44 * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
        s += coeffs[21] * kSH_Y_45 * (7.0f * z2 - 3.0f) * xz;
        s += coeffs[22] * kSH_Y_46 * (7.0f * z2 - 1.0f) * (x2 - y2);
        s += coeffs[23] * kSH_Y_47 * (x2 - 3.0f * y2)   * xz;
        s += coeffs[24] * kSH_Y_48 * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);

        if (numBands < 6)
            return s;

        coeffs += 25;

        for (int l = 5; l < numBands; l++)
            for (int m = -l; m <= l; m++)
                s += (*coeffs++) * SH(l, m, v);

        return s;
    }
}

float SHL::SampleSH(Vec3f v, int numBands, const float coeffs[])
{
    return ::SampleSH<float>(v, numBands, coeffs);
}

Vec4f SHL::SampleSH(Vec3f v, int numBands, const Vec4f coeffs[])
{
    return ::SampleSH<Vec4f>(v, numBands, coeffs);
}


void SHL::CalcSHWeights(Vec3f v, int numBands, float w[])
{
    VL_ASSERT(numBands >= 1);

    w[0] = kSH_Y_00;

    if (numBands < 2)
        return;

    const float& x  = v[0];
    const float& y  = v[1];
    const float& z  = v[2];

    float  x2 = sqr(x);
    float  y2 = sqr(y);
    float  z2 = sqr(z);

    w[1] = kSH_Y_10 * y;
    w[2] = kSH_Y_11 * z;
    w[3] = kSH_Y_12 * x;

    if (numBands < 3)
        return;

    float  xy = x * y;
    float  yz = y * z;
    float  xz = x * z;

    w[4] = kSH_Y_20 * xy;
    w[5] = kSH_Y_21 * yz;
    w[6] = kSH_Y_22 * (3.0f * z2 - 1.0f);
    w[7] = kSH_Y_23 * xz;
    w[8] = kSH_Y_24 * (x2 - y2);

    if (numBands < 4)
        return;

    w[ 9] = kSH_Y_30 * (3.0f * x2 - y2)   * y;
    w[10] = kSH_Y_31 * xy                 * z;
    w[11] = kSH_Y_32 * (5.0f * z2 - 1.0f) * y;
    w[12] = kSH_Y_33 * (5.0f * z2 - 3.0f) * z;
    w[13] = kSH_Y_34 * (5.0f * z2 - 1.0f) * x;
    w[14] = kSH_Y_35 * (x2 - y2)          * z;
    w[15] = kSH_Y_36 * (x2 - 3.0f * y2)   * x;

    if (numBands < 5)
        return;

    w[16] = kSH_Y_40 * (x2 - y2)          * xy;
    w[17] = kSH_Y_41 * (3.0f * x2 - y2)   * yz;
    w[18] = kSH_Y_42 * (7.0f * z2 - 1.0f) * xy;
    w[19] = kSH_Y_43 * (7.0f * z2 - 3.0f) * yz;
    w[20] = kSH_Y_44 * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
    w[21] = kSH_Y_45 * (7.0f * z2 - 3.0f) * xz;
    w[22] = kSH_Y_46 * (7.0f * z2 - 1.0f) * (x2 - y2);
    w[23] = kSH_Y_47 * (x2 - 3.0f * y2)   * xz;
    w[24] = kSH_Y_48 * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);

    if (numBands < 6)
        return;

    w += 25;

    for (int l = 5; l < numBands; l++)
        for (int m = -l; m <= l; m++)
            *w++ = SH(l, m, v);
}


void SHL::AddSHSample(float s, Vec3f v, int numBands, float coeffs[])
{
    VL_ASSERT(numBands >= 1);

    coeffs[0] += s * kSH_Y_00;

    if (numBands < 2)
        return;

    const float& x = v[0];
    const float& y = v[1];
    const float& z = v[2];

    coeffs[1] += s * kSH_Y_10 * y;
    coeffs[2] += s * kSH_Y_11 * z;
    coeffs[3] += s * kSH_Y_12 * x;

    if (numBands < 3)
        return;

    float x2 = sqr(x);
    float y2 = sqr(y);
    float z2 = sqr(z);
    float xy = x * y;
    float yz = y * z;
    float xz = x * z;

    coeffs[4] += s * kSH_Y_20 * xy;
    coeffs[5] += s * kSH_Y_21 * yz;
    coeffs[6] += s * kSH_Y_22 * (3.0f * z2 - 1.0f);
    coeffs[7] += s * kSH_Y_23 * xz;
    coeffs[8] += s * kSH_Y_24 * (x2 - y2);

    if (numBands < 4)
        return;

    coeffs[ 9] += s * kSH_Y_30 * (3.0f * x2 - y2) * y;
    coeffs[10] += s * kSH_Y_31 * xy            * z;
    coeffs[11] += s * kSH_Y_32 * (5.0f * z2 - 1.0f)  * y;
    coeffs[12] += s * kSH_Y_33 * (5.0f * z2 - 3.0f)  * z;
    coeffs[13] += s * kSH_Y_34 * (5.0f * z2 - 1.0f)  * x;
    coeffs[14] += s * kSH_Y_35 * (x2 - y2)     * z;
    coeffs[15] += s * kSH_Y_36 * (x2 - 3.0f * y2) * x;

    if (numBands < 5)
        return;

    coeffs[16] += s * kSH_Y_40 * (x2 - y2)     * xy;
    coeffs[17] += s * kSH_Y_41 * (3.0f * x2 - y2) * yz;
    coeffs[18] += s * kSH_Y_42 * (7.0f * z2 - 1.0f)  * xy;
    coeffs[19] += s * kSH_Y_43 * (7.0f * z2 - 3.0f)  * yz;
    coeffs[20] += s * kSH_Y_44 * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
    coeffs[21] += s * kSH_Y_45 * (7.0f * z2 - 3.0f)  * xz;
    coeffs[22] += s * kSH_Y_46 * (7.0f * z2 - 1.0f)  * (x2 - y2);
    coeffs[23] += s * kSH_Y_47 * (x2 - 3.0f * y2) * xz;
    coeffs[24] += s * kSH_Y_48 * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);

    if (numBands < 6)
        return;

    coeffs += 25;

    for (int l = 5; l < numBands; l++)
        for (int m = -l; m <= l; m++)
            *coeffs++ = s * SH(l, m, v);
}

void SHL::AddSHSample(Vec4f s, Vec3f v, int numBands, Vec4f coeffs[])
{
    VL_ASSERT(numBands >= 1);

    coeffs[0] += s * kSH_Y_00;

    if (numBands < 2)
        return;

    const float& x = v[0];
    const float& y = v[1];
    const float& z = v[2];

    coeffs[1] += s * kSH_Y_10 * y;
    coeffs[2] += s * kSH_Y_11 * z;
    coeffs[3] += s * kSH_Y_12 * x;

    if (numBands < 3)
        return;

    float x2 = sqr(x);
    float y2 = sqr(y);
    float z2 = sqr(z);
    float xy = x * y;
    float yz = y * z;
    float xz = x * z;

    coeffs[4] += s * kSH_Y_20 * xy;
    coeffs[5] += s * kSH_Y_21 * yz;
    coeffs[6] += s * kSH_Y_22 * (3.0f * z2 - 1.0f);
    coeffs[7] += s * kSH_Y_23 * xz;
    coeffs[8] += s * kSH_Y_24 * (x2 - y2);

    if (numBands < 4)
        return;

    coeffs[ 9] += s * kSH_Y_30 * (3.0f * x2 - y2)   * y;
    coeffs[10] += s * kSH_Y_31 * xy                 * z;
    coeffs[11] += s * kSH_Y_32 * (5.0f * z2 - 1.0f) * y;
    coeffs[12] += s * kSH_Y_33 * (5.0f * z2 - 3.0f) * z;
    coeffs[13] += s * kSH_Y_34 * (5.0f * z2 - 1.0f) * x;
    coeffs[14] += s * kSH_Y_35 * (x2 - y2)          * z;
    coeffs[15] += s * kSH_Y_36 * (x2 - 3.0f * y2)   * x;

    if (numBands < 5)
        return;

    coeffs[16] += s * kSH_Y_40 * (x2 - y2)          * xy;
    coeffs[17] += s * kSH_Y_41 * (3.0f * x2 - y2)   * yz;
    coeffs[18] += s * kSH_Y_42 * (7.0f * z2 - 1.0f) * xy;
    coeffs[19] += s * kSH_Y_43 * (7.0f * z2 - 3.0f) * yz;
    coeffs[20] += s * kSH_Y_44 * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
    coeffs[21] += s * kSH_Y_45 * (7.0f * z2 - 3.0f) * xz;
    coeffs[22] += s * kSH_Y_46 * (7.0f * z2 - 1.0f) * (x2 - y2);
    coeffs[23] += s * kSH_Y_47 * (x2 - 3.0f * y2)   * xz;
    coeffs[24] += s * kSH_Y_48 * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);

    if (numBands < 6)
        return;

    coeffs += 25;

    for (int l = 5; l < numBands; l++)
        for (int m = -l; m <= l; m++)
            *coeffs++ = s * SH(l, m, v);
}

namespace
{
    template<typename T> inline void AddSHSampleHemi(T sIn, Vec3f v, int numBands, T coeffs[])
    {
        VL_ASSERT(numBands >= 1);

        T s(2.0f * sIn);

        coeffs[0] += s * kSH_Y_00;

        if (numBands < 2)
            return;

        const float& x  = v[0];
        const float& y  = v[1];

        coeffs[1] += s * kSH_Y_10 * y;
        coeffs[3] += s * kSH_Y_12 * x;

        if (numBands < 3)
            return;

        const float& z  = v[2];
        float x2 = sqr(x);
        float y2 = sqr(y);
        float z2 = sqr(z);
        float xy = x * y;

        coeffs[4] += s * kSH_Y_20 * xy;
        coeffs[6] += s * kSH_Y_22 * (3.0f * z2 - 1.0f);
        coeffs[8] += s * kSH_Y_24 * (x2 - y2);

        if (numBands < 4)
            return;

        coeffs[ 9] += s * kSH_Y_30 * (3.0f * x2 -        y2) * y;
        coeffs[11] += s * kSH_Y_32 * (5.0f * z2 - 1.0f     ) * y;
        coeffs[13] += s * kSH_Y_34 * (5.0f * z2 - 1.0f     ) * x;
        coeffs[15] += s * kSH_Y_36 * (       x2 - 3.0f * y2) * x;

        if (numBands < 5)
            return;

        coeffs[16] += s * kSH_Y_40 * (x2 - y2)           * xy;
        coeffs[18] += s * kSH_Y_42 * (7.0f * z2 - 1.0f)  * xy;
        coeffs[20] += s * kSH_Y_44 * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
        coeffs[22] += s * kSH_Y_46 * (7.0f * z2 - 1.0f)  * (x2 - y2);
        coeffs[24] += s * kSH_Y_48 * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);

        if (numBands < 6)
            return;

        coeffs += 25;

        for (int l = 5; l < numBands; l++)
        {
            for (int m = -l; m <= l; m += 2)
                coeffs[m + l] += s * SH(l, m, v);
            coeffs += 2 * l + 1;
        }
    }
}

void SHL::AddSHSampleHemi(float s, Vec3f v, int numBands, float coeffs[])
{
    ::AddSHSampleHemi<float>(s, v, numBands, coeffs);
}

void SHL::AddSHSampleHemi(Vec4f s, Vec3f v, int numBands, Vec4f coeffs[])
{
    ::AddSHSampleHemi<Vec4f>(s, v, numBands, coeffs);
}


namespace
{
    // The VS compiler amongst others doesn't convert sqrt to constants at compile time,
    // so we pull out these constants to avoid lots of sqrt calls in Rotate*.
    // You'd think constexpr would solve this, but not in over a decade.
    const float kSqrt02_01  = sqrt( 2.0 /  1.0);
    const float kSqrt01_02  = sqrt( 1.0 /  2.0);
    const float kSqrt03_02  = sqrt( 3.0 /  2.0);
    const float kSqrt01_03  = sqrt( 1.0 /  3.0);
    const float kSqrt02_03  = sqrt( 2.0 /  3.0);
    const float kSqrt04_03  = sqrt( 4.0 /  3.0);
    const float kSqrt01_04  = sqrt( 1.0 /  4.0);
    const float kSqrt03_04  = sqrt( 3.0 /  4.0);
    const float kSqrt05_04  = sqrt( 5.0 /  4.0);
    const float kSqrt01_05  = sqrt( 1.0 /  5.0);
    const float kSqrt02_05  = sqrt( 2.0 /  5.0);
    const float kSqrt03_05  = sqrt( 3.0 /  5.0);
    const float kSqrt04_05  = sqrt( 4.0 /  5.0);
    const float kSqrt06_05  = sqrt( 6.0 /  5.0);
    const float kSqrt08_05  = sqrt( 8.0 /  5.0);
    const float kSqrt09_05  = sqrt( 9.0 /  5.0);
    const float kSqrt01_06  = sqrt( 1.0 /  6.0);
    const float kSqrt05_06  = sqrt( 5.0 /  6.0);
    const float kSqrt07_06  = sqrt( 7.0 /  6.0);
    const float kSqrt02_07  = sqrt(02.0 /  7.0);
    const float kSqrt06_07  = sqrt( 6.0 /  7.0);
    const float kSqrt10_07  = sqrt(10.0 /  7.0);
    const float kSqrt12_07  = sqrt(12.0 /  7.0);
    const float kSqrt15_07  = sqrt(15.0 /  7.0);
    const float kSqrt16_07  = sqrt(16.0 /  7.0);
    const float kSqrt01_08  = sqrt( 1.0 /  8.0);
    const float kSqrt03_08  = sqrt( 3.0 /  8.0);
    const float kSqrt05_08  = sqrt( 5.0 /  8.0);
    const float kSqrt07_08  = sqrt( 7.0 /  8.0);
    const float kSqrt09_08  = sqrt( 9.0 /  8.0);
    const float kSqrt05_09  = sqrt( 5.0 /  9.0);
    const float kSqrt08_09  = sqrt( 8.0 /  9.0);
    const float kSqrt01_10  = sqrt( 1.0 / 10.0);
    const float kSqrt03_10  = sqrt( 3.0 / 10.0);
    const float kSqrt07_10  = sqrt( 7.0 / 10.0);
    const float kSqrt09_10  = sqrt( 9.0 / 10.0);
    const float kSqrt01_12  = sqrt( 1.0 / 12.0);
    const float kSqrt07_12  = sqrt( 7.0 / 12.0);
    const float kSqrt11_12  = sqrt(11.0 / 12.0);
    const float kSqrt01_14  = sqrt( 1.0 / 14.0);
    const float kSqrt03_14  = sqrt( 3.0 / 14.0);
    const float kSqrt15_14  = sqrt(15.0 / 14.0);
    const float kSqrt04_15  = sqrt( 4.0 / 15.0);
    const float kSqrt07_15  = sqrt( 7.0 / 10.0);
    const float kSqrt14_15  = sqrt(14.0 / 15.0);
    const float kSqrt16_15  = sqrt(16.0 / 15.0);
    const float kSqrt01_16  = sqrt( 1.0 / 16.0);
    const float kSqrt03_16  = sqrt( 3.0 / 16.0);
    const float kSqrt07_16  = sqrt( 7.0 / 16.0);
    const float kSqrt15_16  = sqrt(15.0 / 16.0);
    const float kSqrt01_18  = sqrt( 1.0 / 18.0);
    const float kSqrt01_24  = sqrt( 1.0 / 24.0);
    const float kSqrt03_25  = sqrt( 3.0 / 25.0);
    const float kSqrt09_25  = sqrt( 9.0 / 25.0);
    const float kSqrt14_25  = sqrt(14.0 / 25.0);
    const float kSqrt16_25  = sqrt(16.0 / 25.0);
    const float kSqrt18_25  = sqrt(18.0 / 25.0);
    const float kSqrt21_25  = sqrt(21.0 / 25.0);
    const float kSqrt24_25  = sqrt(24.0 / 25.0);
    const float kSqrt03_28  = sqrt( 3.0 / 28.0);
    const float kSqrt05_28  = sqrt( 5.0 / 28.0);
    const float kSqrt01_30  = sqrt( 1.0 / 30.0);
    const float kSqrt01_32  = sqrt( 1.0 / 32.0);
    const float kSqrt03_32  = sqrt( 3.0 / 32.0);
    const float kSqrt15_32  = sqrt(15.0 / 32.0);
    const float kSqrt21_32  = sqrt(21.0 / 32.0);
    const float kSqrt11_36  = sqrt(11.0 / 36.0);
    const float kSqrt35_36  = sqrt(35.0 / 36.0);
    const float kSqrt01_50  = sqrt( 1.0 / 50.0);
    const float kSqrt03_50  = sqrt( 3.0 / 50.0);
    const float kSqrt21_50  = sqrt(21.0 / 50.0);
    const float kSqrt15_56  = sqrt(15.0 / 56.0);
    const float kSqrt01_60  = sqrt( 1.0 / 60.0);
    const float kSqrt01_112 = sqrt( 1.0 / 112.0);
    const float kSqrt03_112 = sqrt( 3.0 / 112.0);
    const float kSqrt15_112 = sqrt(15.0 / 112.0);

    template<typename T> void RotateZHToSH(const Vec3f& dir, int n, const T* zcoeffs, T* coeffs)
    {
        VL_ASSERT(zcoeffs != coeffs);
        VL_ASSERT(n >= 2 && n <= 8);

        coeffs[0] = zcoeffs[0];

        if (n < 2)
            return;

        coeffs++;

        float sh1[3] = { dir[1], dir[2], dir[0] };
        coeffs[0] = zcoeffs[1] * sh1[0];
        coeffs[1] = zcoeffs[1] * sh1[1];
        coeffs[2] = zcoeffs[1] * sh1[2];

        if (n < 3)
            return;

        coeffs += 3;
        float sh2[5];

        sh2[0] =  kSqrt03_04 * (sh1[2] * sh1[0] + sh1[0] * sh1[2]);
        sh2[1] =  kSqrt03_04 * (sh1[1] * sh1[0] + sh1[0] * sh1[1]);
        sh2[2] = -kSqrt01_04 * (sh1[2] * sh1[2] + sh1[0] * sh1[0]) + sh1[1] * sh1[1];
        sh2[3] =  kSqrt03_04 * (sh1[1] * sh1[2] + sh1[2] * sh1[1]);
        sh2[4] =  kSqrt03_04 * (sh1[2] * sh1[2] - sh1[0] * sh1[0]);

        for (int i = 0; i < 5; i++)
            coeffs[i] = zcoeffs[2] * sh2[i];

        if (n < 4)
            return;

        coeffs += 5;
        float sh3[7];

        sh3[0] =  kSqrt05_06 * (sh1[2] * sh2[0] + sh1[0] * sh2[4]);
        sh3[1] =  kSqrt05_09 * (sh1[1] * sh2[0] + (sh1[2] * sh2[1] + sh1[0] * sh2[3]));
        sh3[2] =  kSqrt08_09 *  sh1[1] * sh2[1] + kSqrt02_03 *  sh1[0] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[0] - sh1[0] * sh2[4]);
        sh3[3] =                sh1[1] * sh2[2] - kSqrt01_03 * (sh1[2] * sh2[3] + sh1[0] * sh2[1]);
        sh3[4] =  kSqrt08_09 *  sh1[1] * sh2[3] + kSqrt02_03 *  sh1[2] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[4] + sh1[0] * sh2[0]);
        sh3[5] =  kSqrt05_09 * (sh1[1] * sh2[4] + (sh1[2] * sh2[3] - sh1[0] * sh2[1]));
        sh3[6] =  kSqrt05_06 * (sh1[2] * sh2[4] - sh1[0] * sh2[0]);

        for (int i = 0; i < 7; i++)
            coeffs[i] = zcoeffs[3] * sh3[i];

        if (n < 5)
            return;

        coeffs += 7;
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
            coeffs[i] = zcoeffs[4] * sh4[i];

        if (n < 6)
            return;

        coeffs += 9;
        float sh5[11];

        sh5[ 0] = kSqrt09_10 * (sh1[2] * sh4[0] + sh1[0] * sh4[8]);
        sh5[ 1] = kSqrt09_25 * sh1[1] * sh4[0] + kSqrt18_25 * (sh1[2] * sh4[1] + sh1[0] * sh4[7]);
        sh5[ 2] = kSqrt16_25 * sh1[1] * sh4[1] + kSqrt14_25 * (sh1[2] * sh4[2] + sh1[0] * sh4[6]) - kSqrt01_50 * (sh1[2] * sh4[0] - sh1[0] * sh4[8]);
        sh5[ 3] = kSqrt21_25 * sh1[1] * sh4[2] + kSqrt21_50 * (sh1[2] * sh4[3] + sh1[0] * sh4[5]) - kSqrt03_50 * (sh1[2] * sh4[1] - sh1[0] * sh4[7]);
        sh5[ 4] = kSqrt24_25 * sh1[1] * sh4[3] + kSqrt03_05 *  sh1[0] * sh4[4] - kSqrt03_25 * (sh1[2] * sh4[2] - sh1[0] * sh4[6]);
        sh5[ 5] =              sh1[1] * sh4[4] - kSqrt02_05 * (sh1[2] * sh4[5] + sh1[0] * sh4[3]);
        sh5[ 6] = kSqrt24_25 * sh1[1] * sh4[5] + kSqrt03_05 *  sh1[2] * sh4[4] - kSqrt03_25 * (sh1[2] * sh4[6] + sh1[0] * sh4[2]);
        sh5[ 7] = kSqrt21_25 * sh1[1] * sh4[6] + kSqrt21_50 * (sh1[2] * sh4[5] - sh1[0] * sh4[3]) - kSqrt03_50 * (sh1[2] * sh4[7] + sh1[0] * sh4[1]);
        sh5[ 8] = kSqrt16_25 * sh1[1] * sh4[7] + kSqrt14_25 * (sh1[2] * sh4[6] - sh1[0] * sh4[2]) - kSqrt01_50 * (sh1[2] * sh4[8] + sh1[0] * sh4[0]);
        sh5[ 9] = kSqrt09_25 * sh1[1] * sh4[8] + kSqrt18_25 * (sh1[2] * sh4[7] - sh1[0] * sh4[1]);
        sh5[10] = kSqrt09_10 * (sh1[2] * sh4[8] - sh1[0] * sh4[0]);

        for (int i = 0; i < 11; i++)
            coeffs[i] = zcoeffs[5] * sh5[i];

        if (n < 7)
            return;

        coeffs += 11;
        float sh6[13];

        sh6[ 0] = kSqrt11_12 * (sh1[2] * sh5[0] + sh1[0] * sh5[10]);
        sh6[ 1] = kSqrt11_36 * sh1[1] * sh5[0] + sqrtf(55.0 / 72.0) * (sh1[2] * sh5[1] + sh1[0] * sh5[9]);
        sh6[ 2] = kSqrt05_09 * sh1[1] * sh5[1] + kSqrt05_08 * (sh1[2] * sh5[2] + sh1[0] * sh5[8]) - sqrtf(1.0 / 72.0) * (sh1[2] * sh5[0] - sh1[0] * sh5[10]);
        sh6[ 3] = kSqrt03_04 * sh1[1] * sh5[2] + kSqrt01_02 * (sh1[2] * sh5[3] + sh1[0] * sh5[7]) - kSqrt01_24 * (sh1[2] * sh5[1] - sh1[0] * sh5[9]);
        sh6[ 4] = kSqrt08_09 * sh1[1] * sh5[3] + sqrtf(7.0 / 18.0) * (sh1[2] * sh5[4] + sh1[0] * sh5[6]) - kSqrt01_12 * (sh1[2] * sh5[2] - sh1[0] * sh5[8]);
        sh6[ 5] = kSqrt35_36 * sh1[1] * sh5[4] + kSqrt07_12 * sh1[0] * sh5[5] - sqrtf(5.0 / 36.0) * (sh1[2] * sh5[3] - sh1[0] * sh5[7]);
        sh6[ 6] =              sh1[1] * sh5[5] - sqrtf(5.0 / 12.0) * (sh1[2] * sh5[6] + sh1[0] * sh5[4]);
        sh6[ 7] = kSqrt35_36 * sh1[1] * sh5[6] + kSqrt07_12 * sh1[2] * sh5[5] - sqrtf(5.0 / 36.0) * (sh1[2] * sh5[7] + sh1[0] * sh5[3]);
        sh6[ 8] = kSqrt08_09 * sh1[1] * sh5[7] + sqrtf(7.0 / 18.0) * (sh1[2] * sh5[6] - sh1[0] * sh5[4]) - kSqrt01_12 * (sh1[2] * sh5[8] + sh1[0] * sh5[2]);
        sh6[ 9] = kSqrt03_04 * sh1[1] * sh5[8] + kSqrt01_02 * (sh1[2] * sh5[7] - sh1[0] * sh5[3]) - kSqrt01_24 * (sh1[2] * sh5[9] + sh1[0] * sh5[1]);
        sh6[10] = kSqrt05_09 * sh1[1] * sh5[9] + kSqrt05_08 * (sh1[2] * sh5[8] - sh1[0] * sh5[2]) - sqrtf(1.0 / 72.0) * (sh1[2] * sh5[10] + sh1[0] * sh5[0]);
        sh6[11] = kSqrt11_36 * sh1[1] * sh5[10] + sqrtf(55.0 / 72.0) * (sh1[2] * sh5[9] - sh1[0] * sh5[1]);
        sh6[12] = kSqrt11_12 * (sh1[2] * sh5[10] - sh1[0] * sh5[0]);

        for (int i = 0; i < 13; i++)
            coeffs[i] = zcoeffs[6] * sh6[i];

        if (n < 8)
            return;

        coeffs += 13;
        zcoeffs += 7;

        for (int l = 8; l < n; l++)
        {
            T zl = *zcoeffs++ * sqrtf(kFourPi / (2 * l + 1));

            for (int m = -l; m <= +l; m++)
                *coeffs++ = zl * SH(l, m, dir);
        }
    }
}

void SHL::RotateZHToSH(const Vec3f& dir, int n, const float zcoeffs[], float coeffs[])
{
#ifdef REFERENCE
    for (int l = 0; l < n; l++)
    {
        float zl = *zcoeffs++ * sqrtf(kFourPi / (2 * l + 1));

        for (int m = -l; m <= +l; m++)
            *coeffs++ = zl * SH(l, m, dir);
    }
#else
    ::RotateZHToSH<float>(dir, n, zcoeffs, coeffs);
#endif
}

void SHL::RotateZHToSH(const Vec3f& dir, int n, const Vec4f* zcoeffs, Vec4f* coeffs)
{
    ::RotateZHToSH<Vec4f>(dir, n, zcoeffs, coeffs);
}

namespace
{
    template<typename T> void RotateZHToSHAdd(const Vec3f& dir, int n, const T* zcoeffs, T* coeffs)
    {
        VL_ASSERT(zcoeffs != coeffs);
        VL_ASSERT(n >= 2 && n <= 8);

        coeffs[0] += zcoeffs[0];

        if (n < 2)
            return;

        coeffs++;

        float sh1[3] = { dir[1], dir[2], dir[0] };
        coeffs[0] += zcoeffs[1] * sh1[0];
        coeffs[1] += zcoeffs[1] * sh1[1];
        coeffs[2] += zcoeffs[1] * sh1[2];

        if (n < 3)
            return;

        coeffs += 3;
        float sh2[5];

        sh2[0] =  kSqrt03_04 * (sh1[2] * sh1[0] + sh1[0] * sh1[2]);
        sh2[1] =  kSqrt03_04 * (sh1[1] * sh1[0] + sh1[0] * sh1[1]);
        sh2[2] = -kSqrt01_04 * (sh1[2] * sh1[2] + sh1[0] * sh1[0]) + sh1[1] * sh1[1];
        sh2[3] =  kSqrt03_04 * (sh1[1] * sh1[2] + sh1[2] * sh1[1]);
        sh2[4] =  kSqrt03_04 * (sh1[2] * sh1[2] - sh1[0] * sh1[0]);

        for (int i = 0; i < 5; i++)
            coeffs[i] += zcoeffs[2] * sh2[i];

        if (n < 4)
            return;

        coeffs += 5;
        float sh3[7];

        sh3[0] =  kSqrt05_06 * (sh1[2] * sh2[0] + sh1[0] * sh2[4]);
        sh3[1] =  kSqrt05_09 * (sh1[1] * sh2[0] + (sh1[2] * sh2[1] + sh1[0] * sh2[3]));
        sh3[2] =  kSqrt08_09 *  sh1[1] * sh2[1] + kSqrt02_03 *  sh1[0] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[0] - sh1[0] * sh2[4]);
        sh3[3] =                sh1[1] * sh2[2] - kSqrt01_03 * (sh1[2] * sh2[3] + sh1[0] * sh2[1]);
        sh3[4] =  kSqrt08_09 *  sh1[1] * sh2[3] + kSqrt02_03 *  sh1[2] * sh2[2] - kSqrt01_18 * (sh1[2] * sh2[4] + sh1[0] * sh2[0]);
        sh3[5] =  kSqrt05_09 * (sh1[1] * sh2[4] + (sh1[2] * sh2[3] - sh1[0] * sh2[1]));
        sh3[6] =  kSqrt05_06 * (sh1[2] * sh2[4] - sh1[0] * sh2[0]);

        for (int i = 0; i < 7; i++)
            coeffs[i] += zcoeffs[3] * sh3[i];

        if (n < 5)
            return;

        coeffs += 7;
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
            coeffs[i] += zcoeffs[4] * sh4[i];

        if (n < 6)
            return;

        coeffs += 9;
        float sh5[11];

        sh5[ 0] = kSqrt09_10 * (sh1[2] * sh4[0] + sh1[0] * sh4[8]);
        sh5[ 1] = kSqrt09_25 * sh1[1] * sh4[0] + kSqrt18_25 * (sh1[2] * sh4[1] + sh1[0] * sh4[7]);
        sh5[ 2] = kSqrt16_25 * sh1[1] * sh4[1] + kSqrt14_25 * (sh1[2] * sh4[2] + sh1[0] * sh4[6]) - kSqrt01_50 * (sh1[2] * sh4[0] - sh1[0] * sh4[8]);
        sh5[ 3] = kSqrt21_25 * sh1[1] * sh4[2] + kSqrt21_50 * (sh1[2] * sh4[3] + sh1[0] * sh4[5]) - kSqrt03_50 * (sh1[2] * sh4[1] - sh1[0] * sh4[7]);
        sh5[ 4] = kSqrt24_25 * sh1[1] * sh4[3] + kSqrt03_05 *  sh1[0] * sh4[4] - kSqrt03_25 * (sh1[2] * sh4[2] - sh1[0] * sh4[6]);
        sh5[ 5] =              sh1[1] * sh4[4] - kSqrt02_05 * (sh1[2] * sh4[5] + sh1[0] * sh4[3]);
        sh5[ 6] = kSqrt24_25 * sh1[1] * sh4[5] + kSqrt03_05 *  sh1[2] * sh4[4] - kSqrt03_25 * (sh1[2] * sh4[6] + sh1[0] * sh4[2]);
        sh5[ 7] = kSqrt21_25 * sh1[1] * sh4[6] + kSqrt21_50 * (sh1[2] * sh4[5] - sh1[0] * sh4[3]) - kSqrt03_50 * (sh1[2] * sh4[7] + sh1[0] * sh4[1]);
        sh5[ 8] = kSqrt16_25 * sh1[1] * sh4[7] + kSqrt14_25 * (sh1[2] * sh4[6] - sh1[0] * sh4[2]) - kSqrt01_50 * (sh1[2] * sh4[8] + sh1[0] * sh4[0]);
        sh5[ 9] = kSqrt09_25 * sh1[1] * sh4[8] + kSqrt18_25 * (sh1[2] * sh4[7] - sh1[0] * sh4[1]);
        sh5[10] = kSqrt09_10 * (sh1[2] * sh4[8] - sh1[0] * sh4[0]);

        for (int i = 0; i < 11; i++)
            coeffs[i] += zcoeffs[5] * sh5[i];

        if (n < 7)
            return;

        coeffs += 11;
        float sh6[13];

        sh6[ 0] = kSqrt11_12 * (sh1[2] * sh5[0] + sh1[0] * sh5[10]);
        sh6[ 1] = kSqrt11_36 * sh1[1] * sh5[0] + sqrtf(55.0 / 72.0) * (sh1[2] * sh5[1] + sh1[0] * sh5[9]);
        sh6[ 2] = kSqrt05_09 * sh1[1] * sh5[1] + kSqrt05_08 * (sh1[2] * sh5[2] + sh1[0] * sh5[8]) - sqrtf(1.0 / 72.0) * (sh1[2] * sh5[0] - sh1[0] * sh5[10]);
        sh6[ 3] = kSqrt03_04 * sh1[1] * sh5[2] + kSqrt01_02 * (sh1[2] * sh5[3] + sh1[0] * sh5[7]) - kSqrt01_24 * (sh1[2] * sh5[1] - sh1[0] * sh5[9]);
        sh6[ 4] = kSqrt08_09 * sh1[1] * sh5[3] + sqrtf(7.0 / 18.0) * (sh1[2] * sh5[4] + sh1[0] * sh5[6]) - kSqrt01_12 * (sh1[2] * sh5[2] - sh1[0] * sh5[8]);
        sh6[ 5] = kSqrt35_36 * sh1[1] * sh5[4] + kSqrt07_12 * sh1[0] * sh5[5] - sqrtf(5.0 / 36.0) * (sh1[2] * sh5[3] - sh1[0] * sh5[7]);
        sh6[ 6] =              sh1[1] * sh5[5] - sqrtf(5.0 / 12.0) * (sh1[2] * sh5[6] + sh1[0] * sh5[4]);
        sh6[ 7] = kSqrt35_36 * sh1[1] * sh5[6] + kSqrt07_12 * sh1[2] * sh5[5] - sqrtf(5.0 / 36.0) * (sh1[2] * sh5[7] + sh1[0] * sh5[3]);
        sh6[ 8] = kSqrt08_09 * sh1[1] * sh5[7] + sqrtf(7.0 / 18.0) * (sh1[2] * sh5[6] - sh1[0] * sh5[4]) - kSqrt01_12 * (sh1[2] * sh5[8] + sh1[0] * sh5[2]);
        sh6[ 9] = kSqrt03_04 * sh1[1] * sh5[8] + kSqrt01_02 * (sh1[2] * sh5[7] - sh1[0] * sh5[3]) - kSqrt01_24 * (sh1[2] * sh5[9] + sh1[0] * sh5[1]);
        sh6[10] = kSqrt05_09 * sh1[1] * sh5[9] + kSqrt05_08 * (sh1[2] * sh5[8] - sh1[0] * sh5[2]) - sqrtf(1.0 / 72.0) * (sh1[2] * sh5[10] + sh1[0] * sh5[0]);
        sh6[11] = kSqrt11_36 * sh1[1] * sh5[10] + sqrtf(55.0 / 72.0) * (sh1[2] * sh5[9] - sh1[0] * sh5[1]);
        sh6[12] = kSqrt11_12 * (sh1[2] * sh5[10] - sh1[0] * sh5[0]);

        for (int i = 0; i < 13; i++)
            coeffs[i] += zcoeffs[6] * sh6[i];

        if (n < 8)
            return;

        coeffs += 13;
        zcoeffs += 7;

        for (int l = 8; l < n; l++)
        {
            T zl = *zcoeffs++ * sqrtf(kFourPi / (2 * l + 1));

            for (int m = -l; m <= +l; m++)
                *coeffs++ += zl * SH(l, m, dir);
        }
    }
}

void SHL::RotateZHToSHAdd(const Vec3f& dir, int n, const float zcoeffs[], float coeffs[])
{
    ::RotateZHToSHAdd<float>(dir, n, zcoeffs, coeffs);
}

void SHL::RotateZHToSHAdd(const Vec3f& dir, int n, const Vec4f zcoeffs[], Vec4f coeffs[])
{
    ::RotateZHToSHAdd<Vec4f>(dir, n, zcoeffs, coeffs);
}

namespace
{
    template<typename T> void RotateSHAboutZ(float theta, int n, T* coeffs)
    {
        if (n < 2)
            return;

        // float rc[2 * n - 2];
        VL_ASSERT(n < 64);
        float rc[128];  // 64 bands should be enough for anyone

        rc[0] = cosf(theta);
        rc[1] = sinf(theta);

        T tmp     = coeffs[1] * rc[0] + coeffs[3] * rc[1];
        coeffs[3] = coeffs[3] * rc[0] - coeffs[1] * rc[1];
        coeffs[1] = tmp;

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

                T tmp          = coeffs[ci + i] * rc[k] + coeffs[ci + j] * rc[k + 1];
                coeffs[ci + j] = coeffs[ci + j] * rc[k] - coeffs[ci + i] * rc[k + 1];
                coeffs[ci + i] = tmp;
            }
        }
    }
}

void SHL::RotateSHAboutZ(float theta, int n, float* coeffs)
{
    ::RotateSHAboutZ<float>(theta, n, coeffs);
}

void SHL::RotateSHAboutZ(float theta, int n, Vec4f* coeffs)
{
    ::RotateSHAboutZ<Vec4f>(theta, n, coeffs);
}

void SHL::MirrorSHInZ(int n, float* coeffs)
{
    for (int l = 0; l < n; l++)
        for (int m = -l; m <= l; m++)
            if ((l + m) & 1)
            {
                int index = l * (l + 1) + m;
                coeffs[index] = -coeffs[index];
            }
}

void SHL::MirrorSHInZ(int n, Vec4f* coeffs)
{
    for (int l = 0; l < n; l++)
        for (int m = -l; m <= l; m++)
            if ((l + m) & 1)
            {
                int index = l * (l + 1) + m;
                coeffs[index] = -coeffs[index];
            }
}

namespace
{
    template<typename T> inline T dp(int n, const T* a, const float* b)
    {
        T result = (*a) * (*b);

        while (--n > 0)
        {
            a++;
            b++;
            result += (*a) * (*b);
        }

        return result;
    }

    template<typename T> void RotateSH(const Mat3f& rot, int n, const T* coeffsIn, T* coeffs)
    /// We could make faster "specialized" versions, as for the last band we don't
    /// need the full sh matrix, just a row.
    {
        VL_ASSERT(coeffsIn != coeffs);
        VL_ASSERT(n >= 1 && n <= 8);

        (*coeffs++) = coeffsIn[0];

        if (n < 2)
            return;

        coeffsIn += 1;

        float sh1[3][3] =
        {
            { rot.y.y, rot.z.y, rot.x.y },
            { rot.y.z, rot.z.z, rot.x.z },
            { rot.y.x, rot.z.x, rot.x.x }
        };

        (*coeffs++) = dp(3, coeffsIn, sh1[0]);
        (*coeffs++) = dp(3, coeffsIn, sh1[1]);
        (*coeffs++) = dp(3, coeffsIn, sh1[2]);

        if (n < 3)
            return;

        coeffsIn += 3;
        float sh2[5][5];

        sh2[0][0] = kSqrt01_04 * ((sh1[2][2] * sh1[0][0] + sh1[2][0] * sh1[0][2]) + (sh1[0][2] * sh1[2][0] + sh1[0][0] * sh1[2][2]));
        sh2[0][1] =               (sh1[2][1] * sh1[0][0] + sh1[0][1] * sh1[2][0]);
        sh2[0][2] = kSqrt03_04 *  (sh1[2][1] * sh1[0][1] + sh1[0][1] * sh1[2][1]);
        sh2[0][3] =               (sh1[2][1] * sh1[0][2] + sh1[0][1] * sh1[2][2]);
        sh2[0][4] = kSqrt01_04 * ((sh1[2][2] * sh1[0][2] - sh1[2][0] * sh1[0][0]) + (sh1[0][2] * sh1[2][2] - sh1[0][0] * sh1[2][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[0]);

        sh2[1][0] = kSqrt01_04 * ((sh1[1][2] * sh1[0][0] + sh1[1][0] * sh1[0][2]) + (sh1[0][2] * sh1[1][0] + sh1[0][0] * sh1[1][2]));
        sh2[1][1] =                sh1[1][1] * sh1[0][0] + sh1[0][1] * sh1[1][0];
        sh2[1][2] = kSqrt03_04 *  (sh1[1][1] * sh1[0][1] + sh1[0][1] * sh1[1][1]);
        sh2[1][3] =                sh1[1][1] * sh1[0][2] + sh1[0][1] * sh1[1][2];
        sh2[1][4] = kSqrt01_04 * ((sh1[1][2] * sh1[0][2] - sh1[1][0] * sh1[0][0]) + (sh1[0][2] * sh1[1][2] - sh1[0][0] * sh1[1][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[1]);

        sh2[2][0] = kSqrt01_03 * (sh1[1][2] * sh1[1][0] + sh1[1][0] * sh1[1][2]) - kSqrt01_12 * ((sh1[2][2] * sh1[2][0] + sh1[2][0] * sh1[2][2]) + (sh1[0][2] * sh1[0][0] + sh1[0][0] * sh1[0][2]));
        sh2[2][1] = kSqrt04_03 *  sh1[1][1] * sh1[1][0] - kSqrt01_03 * (sh1[2][1] * sh1[2][0] + sh1[0][1] * sh1[0][0]);
        sh2[2][2] =               sh1[1][1] * sh1[1][1] - kSqrt01_04 * (sh1[2][1] * sh1[2][1] + sh1[0][1] * sh1[0][1]);
        sh2[2][3] = kSqrt04_03 *  sh1[1][1] * sh1[1][2] - kSqrt01_03 * (sh1[2][1] * sh1[2][2] + sh1[0][1] * sh1[0][2]);
        sh2[2][4] = kSqrt01_03 * (sh1[1][2] * sh1[1][2] - sh1[1][0] * sh1[1][0]) - kSqrt01_12 * ((sh1[2][2] * sh1[2][2] - sh1[2][0] * sh1[2][0]) + (sh1[0][2] * sh1[0][2] - sh1[0][0] * sh1[0][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[2]);

        sh2[3][0] = kSqrt01_04 * ((sh1[1][2] * sh1[2][0] + sh1[1][0] * sh1[2][2]) + (sh1[2][2] * sh1[1][0] + sh1[2][0] * sh1[1][2]));
        sh2[3][1] =                sh1[1][1] * sh1[2][0] + sh1[2][1] * sh1[1][0];
        sh2[3][2] = kSqrt03_04 *  (sh1[1][1] * sh1[2][1] + sh1[2][1] * sh1[1][1]);
        sh2[3][3] =                sh1[1][1] * sh1[2][2] + sh1[2][1] * sh1[1][2];
        sh2[3][4] = kSqrt01_04 * ((sh1[1][2] * sh1[2][2] - sh1[1][0] * sh1[2][0]) + (sh1[2][2] * sh1[1][2] - sh1[2][0] * sh1[1][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[3]);

        sh2[4][0] = kSqrt01_04 * ((sh1[2][2] * sh1[2][0] + sh1[2][0] * sh1[2][2]) - (sh1[0][2] * sh1[0][0] + sh1[0][0] * sh1[0][2]));
        sh2[4][1] =               (sh1[2][1] * sh1[2][0] - sh1[0][1] * sh1[0][0]);
        sh2[4][2] = kSqrt03_04 *  (sh1[2][1] * sh1[2][1] - sh1[0][1] * sh1[0][1]);
        sh2[4][3] =               (sh1[2][1] * sh1[2][2] - sh1[0][1] * sh1[0][2]);
        sh2[4][4] = kSqrt01_04 * ((sh1[2][2] * sh1[2][2] - sh1[2][0] * sh1[2][0]) - (sh1[0][2] * sh1[0][2] - sh1[0][0] * sh1[0][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[4]);

        if (n < 4)
            return;

        coeffsIn += 5;
        float sh3[7][7];

        sh3[0][0] = kSqrt01_04 * ((sh1[2][2] * sh2[0][0] + sh1[2][0] * sh2[0][4]) + (sh1[0][2] * sh2[4][0] + sh1[0][0] * sh2[4][4]));
        sh3[0][1] = kSqrt03_02 *  (sh1[2][1] * sh2[0][0] + sh1[0][1] * sh2[4][0]);
        sh3[0][2] = kSqrt15_16 *  (sh1[2][1] * sh2[0][1] + sh1[0][1] * sh2[4][1]);
        sh3[0][3] = kSqrt05_06 *  (sh1[2][1] * sh2[0][2] + sh1[0][1] * sh2[4][2]);
        sh3[0][4] = kSqrt15_16 *  (sh1[2][1] * sh2[0][3] + sh1[0][1] * sh2[4][3]);
        sh3[0][5] = kSqrt03_02 *  (sh1[2][1] * sh2[0][4] + sh1[0][1] * sh2[4][4]);
        sh3[0][6] = kSqrt01_04 * ((sh1[2][2] * sh2[0][4] - sh1[2][0] * sh2[0][0]) + (sh1[0][2] * sh2[4][4] - sh1[0][0] * sh2[4][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[0]);

        sh3[1][0] = kSqrt01_06 * (sh1[1][2] * sh2[0][0] + sh1[1][0] * sh2[0][4]) + kSqrt01_06 * ((sh1[2][2] * sh2[1][0] + sh1[2][0] * sh2[1][4]) + (sh1[0][2] * sh2[3][0] + sh1[0][0] * sh2[3][4]));
        sh3[1][1] =               sh1[1][1] * sh2[0][0]                          +               (sh1[2][1] * sh2[1][0] + sh1[0][1] * sh2[3][0]);
        sh3[1][2] = kSqrt05_08 *  sh1[1][1] * sh2[0][1]                          + kSqrt05_08 *  (sh1[2][1] * sh2[1][1] + sh1[0][1] * sh2[3][1]);
        sh3[1][3] = kSqrt05_09 *  sh1[1][1] * sh2[0][2]                          + kSqrt05_09 *  (sh1[2][1] * sh2[1][2] + sh1[0][1] * sh2[3][2]);
        sh3[1][4] = kSqrt05_08 *  sh1[1][1] * sh2[0][3]                          + kSqrt05_08 *  (sh1[2][1] * sh2[1][3] + sh1[0][1] * sh2[3][3]);
        sh3[1][5] =               sh1[1][1] * sh2[0][4]                          +               (sh1[2][1] * sh2[1][4] + sh1[0][1] * sh2[3][4]);
        sh3[1][6] = kSqrt01_06 * (sh1[1][2] * sh2[0][4] - sh1[1][0] * sh2[0][0]) + kSqrt01_06 * ((sh1[2][2] * sh2[1][4] - sh1[2][0] * sh2[1][0]) + (sh1[0][2] * sh2[3][4] - sh1[0][0] * sh2[3][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[1]);

        sh3[2][0] = kSqrt04_15 * (sh1[1][2] * sh2[1][0] + sh1[1][0] * sh2[1][4]) + kSqrt01_05 * (sh1[0][2] * sh2[2][0] + sh1[0][0] * sh2[2][4]) - kSqrt01_60 * ((sh1[2][2] * sh2[0][0] + sh1[2][0] * sh2[0][4]) - (sh1[0][2] * sh2[4][0] + sh1[0][0] * sh2[4][4]));
        sh3[2][1] = kSqrt08_05 *  sh1[1][1] * sh2[1][0]                          + kSqrt06_05 *  sh1[0][1] * sh2[2][0] - kSqrt01_10 * (sh1[2][1] * sh2[0][0] - sh1[0][1] * sh2[4][0]);
        sh3[2][2] =               sh1[1][1] * sh2[1][1]                          + kSqrt03_04 *  sh1[0][1] * sh2[2][1] - kSqrt01_16 * (sh1[2][1] * sh2[0][1] - sh1[0][1] * sh2[4][1]);
        sh3[2][3] = kSqrt08_09 *  sh1[1][1] * sh2[1][2]                          + kSqrt02_03 *  sh1[0][1] * sh2[2][2] - kSqrt01_18 * (sh1[2][1] * sh2[0][2] - sh1[0][1] * sh2[4][2]);
        sh3[2][4] =               sh1[1][1] * sh2[1][3]                          + kSqrt03_04 *  sh1[0][1] * sh2[2][3] - kSqrt01_16 * (sh1[2][1] * sh2[0][3] - sh1[0][1] * sh2[4][3]);
        sh3[2][5] = kSqrt08_05 *  sh1[1][1] * sh2[1][4]                          + kSqrt06_05 *  sh1[0][1] * sh2[2][4] - kSqrt01_10 * (sh1[2][1] * sh2[0][4] - sh1[0][1] * sh2[4][4]);
        sh3[2][6] = kSqrt04_15 * (sh1[1][2] * sh2[1][4] - sh1[1][0] * sh2[1][0]) + kSqrt01_05 * (sh1[0][2] * sh2[2][4] - sh1[0][0] * sh2[2][0]) - kSqrt01_60 * ((sh1[2][2] * sh2[0][4] - sh1[2][0] * sh2[0][0]) - (sh1[0][2] * sh2[4][4] - sh1[0][0] * sh2[4][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[2]);

        sh3[3][0] = kSqrt03_10 * (sh1[1][2] * sh2[2][0] + sh1[1][0] * sh2[2][4]) - kSqrt01_10 * ((sh1[2][2] * sh2[3][0] + sh1[2][0] * sh2[3][4]) + (sh1[0][2] * sh2[1][0] + sh1[0][0] * sh2[1][4]));
        sh3[3][1] = kSqrt09_05 *  sh1[1][1] * sh2[2][0]                          - kSqrt03_05 *  (sh1[2][1] * sh2[3][0] + sh1[0][1] * sh2[1][0]);
        sh3[3][2] = kSqrt09_08 *  sh1[1][1] * sh2[2][1]                          - kSqrt03_08 *  (sh1[2][1] * sh2[3][1] + sh1[0][1] * sh2[1][1]);
        sh3[3][3] =               sh1[1][1] * sh2[2][2]                          - kSqrt01_03 *  (sh1[2][1] * sh2[3][2] + sh1[0][1] * sh2[1][2]);
        sh3[3][4] = kSqrt09_08 *  sh1[1][1] * sh2[2][3]                          - kSqrt03_08 *  (sh1[2][1] * sh2[3][3] + sh1[0][1] * sh2[1][3]);
        sh3[3][5] = kSqrt09_05 *  sh1[1][1] * sh2[2][4]                          - kSqrt03_05 *  (sh1[2][1] * sh2[3][4] + sh1[0][1] * sh2[1][4]);
        sh3[3][6] = kSqrt03_10 * (sh1[1][2] * sh2[2][4] - sh1[1][0] * sh2[2][0]) - kSqrt01_10 * ((sh1[2][2] * sh2[3][4] - sh1[2][0] * sh2[3][0]) + (sh1[0][2] * sh2[1][4] - sh1[0][0] * sh2[1][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[3]);

        sh3[4][0] = kSqrt04_15 * (sh1[1][2] * sh2[3][0] + sh1[1][0] * sh2[3][4]) + kSqrt01_05 * (sh1[2][2] * sh2[2][0] + sh1[2][0] * sh2[2][4]) - kSqrt01_60 * ((sh1[2][2] * sh2[4][0] + sh1[2][0] * sh2[4][4]) + (sh1[0][2] * sh2[0][0] + sh1[0][0] * sh2[0][4]));
        sh3[4][1] = kSqrt08_05 *  sh1[1][1] * sh2[3][0]                          + kSqrt06_05 *  sh1[2][1] * sh2[2][0] - kSqrt01_10 * (sh1[2][1] * sh2[4][0] + sh1[0][1] * sh2[0][0]);
        sh3[4][2] =               sh1[1][1] * sh2[3][1]                          + kSqrt03_04 *  sh1[2][1] * sh2[2][1] - kSqrt01_16 * (sh1[2][1] * sh2[4][1] + sh1[0][1] * sh2[0][1]);
        sh3[4][3] = kSqrt08_09 *  sh1[1][1] * sh2[3][2]                          + kSqrt02_03 *  sh1[2][1] * sh2[2][2] - kSqrt01_18 * (sh1[2][1] * sh2[4][2] + sh1[0][1] * sh2[0][2]);
        sh3[4][4] =               sh1[1][1] * sh2[3][3]                          + kSqrt03_04 *  sh1[2][1] * sh2[2][3] - kSqrt01_16 * (sh1[2][1] * sh2[4][3] + sh1[0][1] * sh2[0][3]);
        sh3[4][5] = kSqrt08_05 *  sh1[1][1] * sh2[3][4]                          + kSqrt06_05 *  sh1[2][1] * sh2[2][4] - kSqrt01_10 * (sh1[2][1] * sh2[4][4] + sh1[0][1] * sh2[0][4]);
        sh3[4][6] = kSqrt04_15 * (sh1[1][2] * sh2[3][4] - sh1[1][0] * sh2[3][0]) + kSqrt01_05 * (sh1[2][2] * sh2[2][4] - sh1[2][0] * sh2[2][0]) - kSqrt01_60 * ((sh1[2][2] * sh2[4][4] - sh1[2][0] * sh2[4][0]) + (sh1[0][2] * sh2[0][4] - sh1[0][0] * sh2[0][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[4]);

        sh3[5][0] = kSqrt01_06 * (sh1[1][2] * sh2[4][0] + sh1[1][0] * sh2[4][4]) + kSqrt01_06 * ((sh1[2][2] * sh2[3][0] + sh1[2][0] * sh2[3][4]) - (sh1[0][2] * sh2[1][0] + sh1[0][0] * sh2[1][4]));
        sh3[5][1] =               sh1[1][1] * sh2[4][0]                          +               (sh1[2][1] * sh2[3][0] - sh1[0][1] * sh2[1][0]);
        sh3[5][2] = kSqrt05_08 *  sh1[1][1] * sh2[4][1]                          + kSqrt05_08 *  (sh1[2][1] * sh2[3][1] - sh1[0][1] * sh2[1][1]);
        sh3[5][3] = kSqrt05_09 *  sh1[1][1] * sh2[4][2]                          + kSqrt05_09 *  (sh1[2][1] * sh2[3][2] - sh1[0][1] * sh2[1][2]);
        sh3[5][4] = kSqrt05_08 *  sh1[1][1] * sh2[4][3]                          + kSqrt05_08 *  (sh1[2][1] * sh2[3][3] - sh1[0][1] * sh2[1][3]);
        sh3[5][5] =               sh1[1][1] * sh2[4][4]                          +               (sh1[2][1] * sh2[3][4] - sh1[0][1] * sh2[1][4]);
        sh3[5][6] = kSqrt01_06 * (sh1[1][2] * sh2[4][4] - sh1[1][0] * sh2[4][0]) + kSqrt01_06 * ((sh1[2][2] * sh2[3][4] - sh1[2][0] * sh2[3][0]) - (sh1[0][2] * sh2[1][4] - sh1[0][0] * sh2[1][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[5]);

        sh3[6][0] = kSqrt01_04 * ((sh1[2][2] * sh2[4][0] + sh1[2][0] * sh2[4][4]) - (sh1[0][2] * sh2[0][0] + sh1[0][0] * sh2[0][4]));
        sh3[6][1] = kSqrt03_02 *  (sh1[2][1] * sh2[4][0] - sh1[0][1] * sh2[0][0]);
        sh3[6][2] = kSqrt15_16 *  (sh1[2][1] * sh2[4][1] - sh1[0][1] * sh2[0][1]);
        sh3[6][3] = kSqrt05_06 *  (sh1[2][1] * sh2[4][2] - sh1[0][1] * sh2[0][2]);
        sh3[6][4] = kSqrt15_16 *  (sh1[2][1] * sh2[4][3] - sh1[0][1] * sh2[0][3]);
        sh3[6][5] = kSqrt03_02 *  (sh1[2][1] * sh2[4][4] - sh1[0][1] * sh2[0][4]);
        sh3[6][6] = kSqrt01_04 * ((sh1[2][2] * sh2[4][4] - sh1[2][0] * sh2[4][0]) - (sh1[0][2] * sh2[0][4] - sh1[0][0] * sh2[0][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[6]);

        if (n < 5)
            return;

        coeffsIn += 7;
        float sh4[9][9];

        sh4[0][0] = kSqrt01_04 * ((sh1[2][2] * sh3[0][0] + sh1[2][0] * sh3[0][6]) + (sh1[0][2] * sh3[6][0] + sh1[0][0] * sh3[6][6]));
        sh4[0][1] = kSqrt02_01 *  (sh1[2][1] * sh3[0][0] + sh1[0][1] * sh3[6][0]);
        sh4[0][2] = kSqrt07_06 *  (sh1[2][1] * sh3[0][1] + sh1[0][1] * sh3[6][1]);
        sh4[0][3] = kSqrt14_15 *  (sh1[2][1] * sh3[0][2] + sh1[0][1] * sh3[6][2]);
        sh4[0][4] = kSqrt07_08 *  (sh1[2][1] * sh3[0][3] + sh1[0][1] * sh3[6][3]);
        sh4[0][5] = kSqrt14_15 *  (sh1[2][1] * sh3[0][4] + sh1[0][1] * sh3[6][4]);
        sh4[0][6] = kSqrt07_06 *  (sh1[2][1] * sh3[0][5] + sh1[0][1] * sh3[6][5]);
        sh4[0][7] = kSqrt02_01 *  (sh1[2][1] * sh3[0][6] + sh1[0][1] * sh3[6][6]);
        sh4[0][8] = kSqrt01_04 * ((sh1[2][2] * sh3[0][6] - sh1[2][0] * sh3[0][0]) + (sh1[0][2] * sh3[6][6] - sh1[0][0] * sh3[6][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[0]);

        sh4[1][0] = kSqrt01_08 * (sh1[1][2] * sh3[0][0] + sh1[1][0] * sh3[0][6]) + kSqrt03_16 * ((sh1[2][2] * sh3[1][0] + sh1[2][0] * sh3[1][6]) + (sh1[0][2] * sh3[5][0] + sh1[0][0] * sh3[5][6]));
        sh4[1][1] =               sh1[1][1] * sh3[0][0]                          + kSqrt03_02 *  (sh1[2][1] * sh3[1][0] + sh1[0][1] * sh3[5][0]);
        sh4[1][2] = kSqrt07_12 *  sh1[1][1] * sh3[0][1]                          + kSqrt07_08 *  (sh1[2][1] * sh3[1][1] + sh1[0][1] * sh3[5][1]);
        sh4[1][3] = kSqrt07_15 *  sh1[1][1] * sh3[0][2]                          + kSqrt07_10 *  (sh1[2][1] * sh3[1][2] + sh1[0][1] * sh3[5][2]);
        sh4[1][4] = kSqrt07_16 *  sh1[1][1] * sh3[0][3]                          + kSqrt21_32 *  (sh1[2][1] * sh3[1][3] + sh1[0][1] * sh3[5][3]);
        sh4[1][5] = kSqrt07_15 *  sh1[1][1] * sh3[0][4]                          + kSqrt07_10 *  (sh1[2][1] * sh3[1][4] + sh1[0][1] * sh3[5][4]);
        sh4[1][6] = kSqrt07_12 *  sh1[1][1] * sh3[0][5]                          + kSqrt07_08 *  (sh1[2][1] * sh3[1][5] + sh1[0][1] * sh3[5][5]);
        sh4[1][7] =               sh1[1][1] * sh3[0][6]                          + kSqrt03_02 *  (sh1[2][1] * sh3[1][6] + sh1[0][1] * sh3[5][6]);
        sh4[1][8] = kSqrt01_08 * (sh1[1][2] * sh3[0][6] - sh1[1][0] * sh3[0][0]) + kSqrt03_16 * ((sh1[2][2] * sh3[1][6] - sh1[2][0] * sh3[1][0]) + (sh1[0][2] * sh3[5][6] - sh1[0][0] * sh3[5][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[1]);

        sh4[2][0] = kSqrt03_14 * (sh1[1][2] * sh3[1][0] + sh1[1][0] * sh3[1][6]) + kSqrt15_112 * ((sh1[2][2] * sh3[2][0] + sh1[2][0] * sh3[2][6]) + (sh1[0][2] * sh3[4][0] + sh1[0][0] * sh3[4][6])) - kSqrt01_112 * ((sh1[2][2] * sh3[0][0] + sh1[2][0] * sh3[0][6]) - (sh1[0][2] * sh3[6][0] + sh1[0][0] * sh3[6][6]));
        sh4[2][1] = kSqrt12_07 *  sh1[1][1] * sh3[1][0]                          + kSqrt15_14  *  (sh1[2][1] * sh3[2][0] + sh1[0][1] * sh3[4][0]) - kSqrt01_14 * (sh1[2][1] * sh3[0][0] - sh1[0][1] * sh3[6][0]);
        sh4[2][2] =               sh1[1][1] * sh3[1][1]                          + kSqrt05_08  *  (sh1[2][1] * sh3[2][1] + sh1[0][1] * sh3[4][1]) - kSqrt01_24 * (sh1[2][1] * sh3[0][1] - sh1[0][1] * sh3[6][1]);
        sh4[2][3] = kSqrt04_05 *  sh1[1][1] * sh3[1][2]                          + kSqrt01_02  *  (sh1[2][1] * sh3[2][2] + sh1[0][1] * sh3[4][2]) - kSqrt01_30 * (sh1[2][1] * sh3[0][2] - sh1[0][1] * sh3[6][2]);
        sh4[2][4] = kSqrt03_04 *  sh1[1][1] * sh3[1][3]                          + kSqrt15_32  *  (sh1[2][1] * sh3[2][3] + sh1[0][1] * sh3[4][3]) - kSqrt01_32 * (sh1[2][1] * sh3[0][3] - sh1[0][1] * sh3[6][3]);
        sh4[2][5] = kSqrt04_05 *  sh1[1][1] * sh3[1][4]                          + kSqrt01_02  *  (sh1[2][1] * sh3[2][4] + sh1[0][1] * sh3[4][4]) - kSqrt01_30 * (sh1[2][1] * sh3[0][4] - sh1[0][1] * sh3[6][4]);
        sh4[2][6] =               sh1[1][1] * sh3[1][5]                          + kSqrt05_08  *  (sh1[2][1] * sh3[2][5] + sh1[0][1] * sh3[4][5]) - kSqrt01_24 * (sh1[2][1] * sh3[0][5] - sh1[0][1] * sh3[6][5]);
        sh4[2][7] = kSqrt12_07 *  sh1[1][1] * sh3[1][6]                          + kSqrt15_14  *  (sh1[2][1] * sh3[2][6] + sh1[0][1] * sh3[4][6]) - kSqrt01_14 * (sh1[2][1] * sh3[0][6] - sh1[0][1] * sh3[6][6]);
        sh4[2][8] = kSqrt03_14 * (sh1[1][2] * sh3[1][6] - sh1[1][0] * sh3[1][0]) + kSqrt15_112 * ((sh1[2][2] * sh3[2][6] - sh1[2][0] * sh3[2][0]) + (sh1[0][2] * sh3[4][6] - sh1[0][0] * sh3[4][0])) - kSqrt01_112 * ((sh1[2][2] * sh3[0][6] - sh1[2][0] * sh3[0][0]) - (sh1[0][2] * sh3[6][6] - sh1[0][0] * sh3[6][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[2]);

        sh4[3][0] = kSqrt15_56 * (sh1[1][2] * sh3[2][0] + sh1[1][0] * sh3[2][6]) + kSqrt05_28 * (sh1[0][2] * sh3[3][0] + sh1[0][0] * sh3[3][6]) - kSqrt03_112 * ((sh1[2][2] * sh3[1][0] + sh1[2][0] * sh3[1][6]) - (sh1[0][2] * sh3[5][0] + sh1[0][0] * sh3[5][6]));
        sh4[3][1] = kSqrt15_07 *  sh1[1][1] * sh3[2][0]                          + kSqrt10_07 *  sh1[0][1] * sh3[3][0] - kSqrt03_14 * (sh1[2][1] * sh3[1][0] - sh1[0][1] * sh3[5][0]);
        sh4[3][2] = kSqrt05_04 *  sh1[1][1] * sh3[2][1]                          + kSqrt05_06 *  sh1[0][1] * sh3[3][1] - kSqrt01_08 * (sh1[2][1] * sh3[1][1] - sh1[0][1] * sh3[5][1]);
        sh4[3][3] =               sh1[1][1] * sh3[2][2]                          + kSqrt02_03 *  sh1[0][1] * sh3[3][2] - kSqrt01_10 * (sh1[2][1] * sh3[1][2] - sh1[0][1] * sh3[5][2]);
        sh4[3][4] = kSqrt15_16 *  sh1[1][1] * sh3[2][3]                          + kSqrt05_08 *  sh1[0][1] * sh3[3][3] - kSqrt03_32 * (sh1[2][1] * sh3[1][3] - sh1[0][1] * sh3[5][3]);
        sh4[3][5] =               sh1[1][1] * sh3[2][4]                          + kSqrt02_03 *  sh1[0][1] * sh3[3][4] - kSqrt01_10 * (sh1[2][1] * sh3[1][4] - sh1[0][1] * sh3[5][4]);
        sh4[3][6] = kSqrt05_04 *  sh1[1][1] * sh3[2][5]                          + kSqrt05_06 *  sh1[0][1] * sh3[3][5] - kSqrt01_08 * (sh1[2][1] * sh3[1][5] - sh1[0][1] * sh3[5][5]);
        sh4[3][7] = kSqrt15_07 *  sh1[1][1] * sh3[2][6]                          + kSqrt10_07 *  sh1[0][1] * sh3[3][6] - kSqrt03_14 * (sh1[2][1] * sh3[1][6] - sh1[0][1] * sh3[5][6]);
        sh4[3][8] = kSqrt15_56 * (sh1[1][2] * sh3[2][6] - sh1[1][0] * sh3[2][0]) + kSqrt05_28 * (sh1[0][2] * sh3[3][6] - sh1[0][0] * sh3[3][0]) - kSqrt03_112 * ((sh1[2][2] * sh3[1][6] - sh1[2][0] * sh3[1][0]) - (sh1[0][2] * sh3[5][6] - sh1[0][0] * sh3[5][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[3]);

        sh4[4][0] = kSqrt02_07 * (sh1[1][2] * sh3[3][0] + sh1[1][0] * sh3[3][6]) - kSqrt03_28 * ((sh1[2][2] * sh3[4][0] + sh1[2][0] * sh3[4][6]) + (sh1[0][2] * sh3[2][0] + sh1[0][0] * sh3[2][6]));
        sh4[4][1] = kSqrt16_07 *  sh1[1][1] * sh3[3][0]                          - kSqrt06_07 *  (sh1[2][1] * sh3[4][0] + sh1[0][1] * sh3[2][0]);
        sh4[4][2] = kSqrt04_03 *  sh1[1][1] * sh3[3][1]                          - kSqrt01_02 *  (sh1[2][1] * sh3[4][1] + sh1[0][1] * sh3[2][1]);
        sh4[4][3] = kSqrt16_15 *  sh1[1][1] * sh3[3][2]                          - kSqrt02_05 *  (sh1[2][1] * sh3[4][2] + sh1[0][1] * sh3[2][2]);
        sh4[4][4] =               sh1[1][1] * sh3[3][3]                          - kSqrt03_08 *  (sh1[2][1] * sh3[4][3] + sh1[0][1] * sh3[2][3]);
        sh4[4][5] = kSqrt16_15 *  sh1[1][1] * sh3[3][4]                          - kSqrt02_05 *  (sh1[2][1] * sh3[4][4] + sh1[0][1] * sh3[2][4]);
        sh4[4][6] = kSqrt04_03 *  sh1[1][1] * sh3[3][5]                          - kSqrt01_02 *  (sh1[2][1] * sh3[4][5] + sh1[0][1] * sh3[2][5]);
        sh4[4][7] = kSqrt16_07 *  sh1[1][1] * sh3[3][6]                          - kSqrt06_07 *  (sh1[2][1] * sh3[4][6] + sh1[0][1] * sh3[2][6]);
        sh4[4][8] = kSqrt02_07 * (sh1[1][2] * sh3[3][6] - sh1[1][0] * sh3[3][0]) - kSqrt03_28 * ((sh1[2][2] * sh3[4][6] - sh1[2][0] * sh3[4][0]) + (sh1[0][2] * sh3[2][6] - sh1[0][0] * sh3[2][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[4]);

        sh4[5][0] = kSqrt15_56 * (sh1[1][2] * sh3[4][0] + sh1[1][0] * sh3[4][6]) + kSqrt05_28 * (sh1[2][2] * sh3[3][0] + sh1[2][0] * sh3[3][6]) - kSqrt03_112 * ((sh1[2][2] * sh3[5][0] + sh1[2][0] * sh3[5][6]) + (sh1[0][2] * sh3[1][0] + sh1[0][0] * sh3[1][6]));
        sh4[5][1] = kSqrt15_07 *  sh1[1][1] * sh3[4][0]                          + kSqrt10_07 *  sh1[2][1] * sh3[3][0] - kSqrt03_14 * (sh1[2][1] * sh3[5][0] + sh1[0][1] * sh3[1][0]);
        sh4[5][2] = kSqrt05_04 *  sh1[1][1] * sh3[4][1]                          + kSqrt05_06 *  sh1[2][1] * sh3[3][1] - kSqrt01_08 * (sh1[2][1] * sh3[5][1] + sh1[0][1] * sh3[1][1]);
        sh4[5][3] =               sh1[1][1] * sh3[4][2]                          + kSqrt02_03 *  sh1[2][1] * sh3[3][2] - kSqrt01_10 * (sh1[2][1] * sh3[5][2] + sh1[0][1] * sh3[1][2]);
        sh4[5][4] = kSqrt15_16 *  sh1[1][1] * sh3[4][3]                          + kSqrt05_08 *  sh1[2][1] * sh3[3][3] - kSqrt03_32 * (sh1[2][1] * sh3[5][3] + sh1[0][1] * sh3[1][3]);
        sh4[5][5] =               sh1[1][1] * sh3[4][4]                          + kSqrt02_03 *  sh1[2][1] * sh3[3][4] - kSqrt01_10 * (sh1[2][1] * sh3[5][4] + sh1[0][1] * sh3[1][4]);
        sh4[5][6] = kSqrt05_04 *  sh1[1][1] * sh3[4][5]                          + kSqrt05_06 *  sh1[2][1] * sh3[3][5] - kSqrt01_08 * (sh1[2][1] * sh3[5][5] + sh1[0][1] * sh3[1][5]);
        sh4[5][7] = kSqrt15_07 *  sh1[1][1] * sh3[4][6]                          + kSqrt10_07 *  sh1[2][1] * sh3[3][6] - kSqrt03_14 * (sh1[2][1] * sh3[5][6] + sh1[0][1] * sh3[1][6]);
        sh4[5][8] = kSqrt15_56 * (sh1[1][2] * sh3[4][6] - sh1[1][0] * sh3[4][0]) + kSqrt05_28 * (sh1[2][2] * sh3[3][6] - sh1[2][0] * sh3[3][0]) - kSqrt03_112 * ((sh1[2][2] * sh3[5][6] - sh1[2][0] * sh3[5][0]) + (sh1[0][2] * sh3[1][6] - sh1[0][0] * sh3[1][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[5]);

        sh4[6][0] = kSqrt03_14 * (sh1[1][2] * sh3[5][0] + sh1[1][0] * sh3[5][6]) + kSqrt15_112 * ((sh1[2][2] * sh3[4][0] + sh1[2][0] * sh3[4][6]) - (sh1[0][2] * sh3[2][0] + sh1[0][0] * sh3[2][6])) - kSqrt01_112 * ((sh1[2][2] * sh3[6][0] + sh1[2][0] * sh3[6][6]) + (sh1[0][2] * sh3[0][0] + sh1[0][0] * sh3[0][6]));
        sh4[6][1] = kSqrt12_07 *  sh1[1][1] * sh3[5][0]                          + kSqrt15_14 *  (sh1[2][1] * sh3[4][0] - sh1[0][1] * sh3[2][0]) - kSqrt01_14 * (sh1[2][1] * sh3[6][0] + sh1[0][1] * sh3[0][0]);
        sh4[6][2] =               sh1[1][1] * sh3[5][1]                          + kSqrt05_08 *  (sh1[2][1] * sh3[4][1] - sh1[0][1] * sh3[2][1]) - kSqrt01_24 * (sh1[2][1] * sh3[6][1] + sh1[0][1] * sh3[0][1]);
        sh4[6][3] = kSqrt04_05 *  sh1[1][1] * sh3[5][2]                          + kSqrt01_02 *  (sh1[2][1] * sh3[4][2] - sh1[0][1] * sh3[2][2]) - kSqrt01_30 * (sh1[2][1] * sh3[6][2] + sh1[0][1] * sh3[0][2]);
        sh4[6][4] = kSqrt03_04 *  sh1[1][1] * sh3[5][3]                          + kSqrt15_32 *  (sh1[2][1] * sh3[4][3] - sh1[0][1] * sh3[2][3]) - kSqrt01_32 * (sh1[2][1] * sh3[6][3] + sh1[0][1] * sh3[0][3]);
        sh4[6][5] = kSqrt04_05 *  sh1[1][1] * sh3[5][4]                          + kSqrt01_02 *  (sh1[2][1] * sh3[4][4] - sh1[0][1] * sh3[2][4]) - kSqrt01_30 * (sh1[2][1] * sh3[6][4] + sh1[0][1] * sh3[0][4]);
        sh4[6][6] =               sh1[1][1] * sh3[5][5]                          + kSqrt05_08 *  (sh1[2][1] * sh3[4][5] - sh1[0][1] * sh3[2][5]) - kSqrt01_24 * (sh1[2][1] * sh3[6][5] + sh1[0][1] * sh3[0][5]);
        sh4[6][7] = kSqrt12_07 *  sh1[1][1] * sh3[5][6]                          + kSqrt15_14 *  (sh1[2][1] * sh3[4][6] - sh1[0][1] * sh3[2][6]) - kSqrt01_14 * (sh1[2][1] * sh3[6][6] + sh1[0][1] * sh3[0][6]);
        sh4[6][8] = kSqrt03_14 * (sh1[1][2] * sh3[5][6] - sh1[1][0] * sh3[5][0]) + kSqrt15_112 * ((sh1[2][2] * sh3[4][6] - sh1[2][0] * sh3[4][0]) - (sh1[0][2] * sh3[2][6] - sh1[0][0] * sh3[2][0])) - kSqrt01_112 * ((sh1[2][2] * sh3[6][6] - sh1[2][0] * sh3[6][0]) + (sh1[0][2] * sh3[0][6] - sh1[0][0] * sh3[0][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[6]);

        sh4[7][0] = kSqrt01_08 * (sh1[1][2] * sh3[6][0] + sh1[1][0] * sh3[6][6]) + kSqrt03_16 * ((sh1[2][2] * sh3[5][0] + sh1[2][0] * sh3[5][6]) - (sh1[0][2] * sh3[1][0] + sh1[0][0] * sh3[1][6]));
        sh4[7][1] =               sh1[1][1] * sh3[6][0]                          + kSqrt03_02 *  (sh1[2][1] * sh3[5][0] - sh1[0][1] * sh3[1][0]);
        sh4[7][2] = kSqrt07_12 *  sh1[1][1] * sh3[6][1]                          + kSqrt07_08 *  (sh1[2][1] * sh3[5][1] - sh1[0][1] * sh3[1][1]);
        sh4[7][3] = kSqrt07_15 *  sh1[1][1] * sh3[6][2]                          + kSqrt07_10 *  (sh1[2][1] * sh3[5][2] - sh1[0][1] * sh3[1][2]);
        sh4[7][4] = kSqrt07_16 *  sh1[1][1] * sh3[6][3]                          + kSqrt21_32 *  (sh1[2][1] * sh3[5][3] - sh1[0][1] * sh3[1][3]);
        sh4[7][5] = kSqrt07_15 *  sh1[1][1] * sh3[6][4]                          + kSqrt07_10 *  (sh1[2][1] * sh3[5][4] - sh1[0][1] * sh3[1][4]);
        sh4[7][6] = kSqrt07_12 *  sh1[1][1] * sh3[6][5]                          + kSqrt07_08 *  (sh1[2][1] * sh3[5][5] - sh1[0][1] * sh3[1][5]);
        sh4[7][7] =               sh1[1][1] * sh3[6][6]                          + kSqrt03_02 *  (sh1[2][1] * sh3[5][6] - sh1[0][1] * sh3[1][6]);
        sh4[7][8] = kSqrt01_08 * (sh1[1][2] * sh3[6][6] - sh1[1][0] * sh3[6][0]) + kSqrt03_16 * ((sh1[2][2] * sh3[5][6] - sh1[2][0] * sh3[5][0]) - (sh1[0][2] * sh3[1][6] - sh1[0][0] * sh3[1][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[7]);

        sh4[8][0] = kSqrt01_04 * ((sh1[2][2] * sh3[6][0] + sh1[2][0] * sh3[6][6]) - (sh1[0][2] * sh3[0][0] + sh1[0][0] * sh3[0][6]));
        sh4[8][1] = kSqrt02_01 *  (sh1[2][1] * sh3[6][0] - sh1[0][1] * sh3[0][0]);
        sh4[8][2] = kSqrt07_06 *  (sh1[2][1] * sh3[6][1] - sh1[0][1] * sh3[0][1]);
        sh4[8][3] = kSqrt14_15 *  (sh1[2][1] * sh3[6][2] - sh1[0][1] * sh3[0][2]);
        sh4[8][4] = kSqrt07_08 *  (sh1[2][1] * sh3[6][3] - sh1[0][1] * sh3[0][3]);
        sh4[8][5] = kSqrt14_15 *  (sh1[2][1] * sh3[6][4] - sh1[0][1] * sh3[0][4]);
        sh4[8][6] = kSqrt07_06 *  (sh1[2][1] * sh3[6][5] - sh1[0][1] * sh3[0][5]);
        sh4[8][7] = kSqrt02_01 *  (sh1[2][1] * sh3[6][6] - sh1[0][1] * sh3[0][6]);
        sh4[8][8] = kSqrt01_04 * ((sh1[2][2] * sh3[6][6] - sh1[2][0] * sh3[6][0]) - (sh1[0][2] * sh3[0][6] - sh1[0][0] * sh3[0][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[8]);

        if (n < 6)
            return;

        coeffsIn += 9;
        float sh5[11][11];

        sh5[0][0] = kSqrt01_04 * ((sh1[2][2] * sh4[0][0] + sh1[2][0] * sh4[0][8]) + (sh1[0][2] * sh4[8][0] + sh1[0][0] * sh4[8][8]));
        sh5[0][1] = sqrtf(5.0 / 2.0) * (sh1[2][1] * sh4[0][0] + sh1[0][1] * sh4[8][0]);
        sh5[0][2] = sqrtf(45.0 / 32.0) * (sh1[2][1] * sh4[0][1] + sh1[0][1] * sh4[8][1]);
        sh5[0][3] = kSqrt15_14 * (sh1[2][1] * sh4[0][2] + sh1[0][1] * sh4[8][2]);
        sh5[0][4] = kSqrt15_16 * (sh1[2][1] * sh4[0][3] + sh1[0][1] * sh4[8][3]);
        sh5[0][5] = kSqrt09_10 * (sh1[2][1] * sh4[0][4] + sh1[0][1] * sh4[8][4]);
        sh5[0][6] = kSqrt15_16 * (sh1[2][1] * sh4[0][5] + sh1[0][1] * sh4[8][5]);
        sh5[0][7] = kSqrt15_14 * (sh1[2][1] * sh4[0][6] + sh1[0][1] * sh4[8][6]);
        sh5[0][8] = sqrtf(45.0 / 32.0) * (sh1[2][1] * sh4[0][7] + sh1[0][1] * sh4[8][7]);
        sh5[0][9] = sqrtf(5.0 / 2.0) * (sh1[2][1] * sh4[0][8] + sh1[0][1] * sh4[8][8]);
        sh5[0][10] = kSqrt01_04 * ((sh1[2][2] * sh4[0][8] - sh1[2][0] * sh4[0][0]) + (sh1[0][2] * sh4[8][8] - sh1[0][0] * sh4[8][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[0]);

        sh5[1][0] = kSqrt01_10 * (sh1[1][2] * sh4[0][0] + sh1[1][0] * sh4[0][8]) + kSqrt01_05 * ((sh1[2][2] * sh4[1][0] + sh1[2][0] * sh4[1][8]) + (sh1[0][2] * sh4[7][0] + sh1[0][0] * sh4[7][8]));
        sh5[1][1] = sh1[1][1] * sh4[0][0] + kSqrt02_01 * (sh1[2][1] * sh4[1][0] + sh1[0][1] * sh4[7][0]);
        sh5[1][2] = sqrtf(9.0 / 16.0) * sh1[1][1] * sh4[0][1] + kSqrt09_08 * (sh1[2][1] * sh4[1][1] + sh1[0][1] * sh4[7][1]);
        sh5[1][3] = sqrtf(3.0 / 7.0) * sh1[1][1] * sh4[0][2] + kSqrt06_07 * (sh1[2][1] * sh4[1][2] + sh1[0][1] * sh4[7][2]);
        sh5[1][4] = kSqrt03_08 * sh1[1][1] * sh4[0][3] + kSqrt03_04 * (sh1[2][1] * sh4[1][3] + sh1[0][1] * sh4[7][3]);
        sh5[1][5] = kSqrt09_25 * sh1[1][1] * sh4[0][4] + kSqrt18_25 * (sh1[2][1] * sh4[1][4] + sh1[0][1] * sh4[7][4]);
        sh5[1][6] = kSqrt03_08 * sh1[1][1] * sh4[0][5] + kSqrt03_04 * (sh1[2][1] * sh4[1][5] + sh1[0][1] * sh4[7][5]);
        sh5[1][7] = sqrtf(3.0 / 7.0) * sh1[1][1] * sh4[0][6] + kSqrt06_07 * (sh1[2][1] * sh4[1][6] + sh1[0][1] * sh4[7][6]);
        sh5[1][8] = sqrtf(9.0 / 16.0) * sh1[1][1] * sh4[0][7] + kSqrt09_08 * (sh1[2][1] * sh4[1][7] + sh1[0][1] * sh4[7][7]);
        sh5[1][9] = sh1[1][1] * sh4[0][8] + kSqrt02_01 * (sh1[2][1] * sh4[1][8] + sh1[0][1] * sh4[7][8]);
        sh5[1][10] = kSqrt01_10 * (sh1[1][2] * sh4[0][8] - sh1[1][0] * sh4[0][0]) + kSqrt01_05 * ((sh1[2][2] * sh4[1][8] - sh1[2][0] * sh4[1][0]) + (sh1[0][2] * sh4[7][8] - sh1[0][0] * sh4[7][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[1]);

        sh5[2][0] = sqrtf(8.0 / 45.0) * (sh1[1][2] * sh4[1][0] + sh1[1][0] * sh4[1][8]) + sqrtf(7.0 / 45.0) * ((sh1[2][2] * sh4[2][0] + sh1[2][0] * sh4[2][8]) + (sh1[0][2] * sh4[6][0] + sh1[0][0] * sh4[6][8])) - sqrtf(1.0 / 180.0) * ((sh1[2][2] * sh4[0][0] + sh1[2][0] * sh4[0][8]) - (sh1[0][2] * sh4[8][0] + sh1[0][0] * sh4[8][8]));
        sh5[2][1] = sqrtf(16.0 / 9.0) * sh1[1][1] * sh4[1][0] + sqrtf(14.0 / 9.0) * (sh1[2][1] * sh4[2][0] + sh1[0][1] * sh4[6][0]) - kSqrt01_18 * (sh1[2][1] * sh4[0][0] - sh1[0][1] * sh4[8][0]);
        sh5[2][2] = sh1[1][1] * sh4[1][1] + kSqrt07_08 * (sh1[2][1] * sh4[2][1] + sh1[0][1] * sh4[6][1]) - kSqrt01_32 * (sh1[2][1] * sh4[0][1] - sh1[0][1] * sh4[8][1]);
        sh5[2][3] = sqrtf(16.0 / 21.0) * sh1[1][1] * sh4[1][2] + kSqrt02_03 * (sh1[2][1] * sh4[2][2] + sh1[0][1] * sh4[6][2]) - sqrtf(1.0 / 42.0) * (sh1[2][1] * sh4[0][2] - sh1[0][1] * sh4[8][2]);
        sh5[2][4] = kSqrt02_03 * sh1[1][1] * sh4[1][3] + kSqrt07_12 * (sh1[2][1] * sh4[2][3] + sh1[0][1] * sh4[6][3]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh4[0][3] - sh1[0][1] * sh4[8][3]);
        sh5[2][5] = sqrtf(16.0 / 25.0) * sh1[1][1] * sh4[1][4] + kSqrt14_25 * (sh1[2][1] * sh4[2][4] + sh1[0][1] * sh4[6][4]) - kSqrt01_50 * (sh1[2][1] * sh4[0][4] - sh1[0][1] * sh4[8][4]);
        sh5[2][6] = kSqrt02_03 * sh1[1][1] * sh4[1][5] + kSqrt07_12 * (sh1[2][1] * sh4[2][5] + sh1[0][1] * sh4[6][5]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh4[0][5] - sh1[0][1] * sh4[8][5]);
        sh5[2][7] = sqrtf(16.0 / 21.0) * sh1[1][1] * sh4[1][6] + kSqrt02_03 * (sh1[2][1] * sh4[2][6] + sh1[0][1] * sh4[6][6]) - sqrtf(1.0 / 42.0) * (sh1[2][1] * sh4[0][6] - sh1[0][1] * sh4[8][6]);
        sh5[2][8] = sh1[1][1] * sh4[1][7] + kSqrt07_08 * (sh1[2][1] * sh4[2][7] + sh1[0][1] * sh4[6][7]) - kSqrt01_32 * (sh1[2][1] * sh4[0][7] - sh1[0][1] * sh4[8][7]);
        sh5[2][9] = sqrtf(16.0 / 9.0) * sh1[1][1] * sh4[1][8] + sqrtf(14.0 / 9.0) * (sh1[2][1] * sh4[2][8] + sh1[0][1] * sh4[6][8]) - kSqrt01_18 * (sh1[2][1] * sh4[0][8] - sh1[0][1] * sh4[8][8]);
        sh5[2][10] = sqrtf(8.0 / 45.0) * (sh1[1][2] * sh4[1][8] - sh1[1][0] * sh4[1][0]) + sqrtf(7.0 / 45.0) * ((sh1[2][2] * sh4[2][8] - sh1[2][0] * sh4[2][0]) + (sh1[0][2] * sh4[6][8] - sh1[0][0] * sh4[6][0])) - sqrtf(1.0 / 180.0) * ((sh1[2][2] * sh4[0][8] - sh1[2][0] * sh4[0][0]) - (sh1[0][2] * sh4[8][8] - sh1[0][0] * sh4[8][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[2]);

        sh5[3][0] = sqrtf(7.0 / 30.0) * (sh1[1][2] * sh4[2][0] + sh1[1][0] * sh4[2][8]) + sqrtf(7.0 / 60.0) * ((sh1[2][2] * sh4[3][0] + sh1[2][0] * sh4[3][8]) + (sh1[0][2] * sh4[5][0] + sh1[0][0] * sh4[5][8])) - kSqrt01_60 * ((sh1[2][2] * sh4[1][0] + sh1[2][0] * sh4[1][8]) - (sh1[0][2] * sh4[7][0] + sh1[0][0] * sh4[7][8]));
        sh5[3][1] = sqrtf(7.0 / 3.0) * sh1[1][1] * sh4[2][0] + kSqrt07_06 * (sh1[2][1] * sh4[3][0] + sh1[0][1] * sh4[5][0]) - kSqrt01_06 * (sh1[2][1] * sh4[1][0] - sh1[0][1] * sh4[7][0]);
        sh5[3][2] = sqrtf(21.0 / 16.0) * sh1[1][1] * sh4[2][1] + kSqrt21_32 * (sh1[2][1] * sh4[3][1] + sh1[0][1] * sh4[5][1]) - kSqrt03_32 * (sh1[2][1] * sh4[1][1] - sh1[0][1] * sh4[7][1]);
        sh5[3][3] = sh1[1][1] * sh4[2][2] + kSqrt01_02 * (sh1[2][1] * sh4[3][2] + sh1[0][1] * sh4[5][2]) - kSqrt01_14 * (sh1[2][1] * sh4[1][2] - sh1[0][1] * sh4[7][2]);
        sh5[3][4] = kSqrt07_08 * sh1[1][1] * sh4[2][3] + kSqrt07_16 * (sh1[2][1] * sh4[3][3] + sh1[0][1] * sh4[5][3]) - kSqrt01_16 * (sh1[2][1] * sh4[1][3] - sh1[0][1] * sh4[7][3]);
        sh5[3][5] = sqrtf(21.0 / 25.0) * sh1[1][1] * sh4[2][4] + kSqrt21_50 * (sh1[2][1] * sh4[3][4] + sh1[0][1] * sh4[5][4]) - kSqrt03_50 * (sh1[2][1] * sh4[1][4] - sh1[0][1] * sh4[7][4]);
        sh5[3][6] = kSqrt07_08 * sh1[1][1] * sh4[2][5] + kSqrt07_16 * (sh1[2][1] * sh4[3][5] + sh1[0][1] * sh4[5][5]) - kSqrt01_16 * (sh1[2][1] * sh4[1][5] - sh1[0][1] * sh4[7][5]);
        sh5[3][7] = sh1[1][1] * sh4[2][6] + kSqrt01_02 * (sh1[2][1] * sh4[3][6] + sh1[0][1] * sh4[5][6]) - kSqrt01_14 * (sh1[2][1] * sh4[1][6] - sh1[0][1] * sh4[7][6]);
        sh5[3][8] = sqrtf(21.0 / 16.0) * sh1[1][1] * sh4[2][7] + kSqrt21_32 * (sh1[2][1] * sh4[3][7] + sh1[0][1] * sh4[5][7]) - kSqrt03_32 * (sh1[2][1] * sh4[1][7] - sh1[0][1] * sh4[7][7]);
        sh5[3][9] = sqrtf(7.0 / 3.0) * sh1[1][1] * sh4[2][8] + kSqrt07_06 * (sh1[2][1] * sh4[3][8] + sh1[0][1] * sh4[5][8]) - kSqrt01_06 * (sh1[2][1] * sh4[1][8] - sh1[0][1] * sh4[7][8]);
        sh5[3][10] = sqrtf(7.0 / 30.0) * (sh1[1][2] * sh4[2][8] - sh1[1][0] * sh4[2][0]) + sqrtf(7.0 / 60.0) * ((sh1[2][2] * sh4[3][8] - sh1[2][0] * sh4[3][0]) + (sh1[0][2] * sh4[5][8] - sh1[0][0] * sh4[5][0])) - kSqrt01_60 * ((sh1[2][2] * sh4[1][8] - sh1[2][0] * sh4[1][0]) - (sh1[0][2] * sh4[7][8] - sh1[0][0] * sh4[7][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[3]);

        sh5[4][0] = kSqrt04_15 * (sh1[1][2] * sh4[3][0] + sh1[1][0] * sh4[3][8]) + kSqrt01_06 * (sh1[0][2] * sh4[4][0] + sh1[0][0] * sh4[4][8]) - kSqrt01_30 * ((sh1[2][2] * sh4[2][0] + sh1[2][0] * sh4[2][8]) - (sh1[0][2] * sh4[6][0] + sh1[0][0] * sh4[6][8]));
        sh5[4][1] = sqrtf(8.0 / 3.0) * sh1[1][1] * sh4[3][0] + sqrtf(5.0 / 3.0) * sh1[0][1] * sh4[4][0] - kSqrt01_03 * (sh1[2][1] * sh4[2][0] - sh1[0][1] * sh4[6][0]);
        sh5[4][2] = kSqrt03_02 * sh1[1][1] * sh4[3][1] + kSqrt15_16 * sh1[0][1] * sh4[4][1] - kSqrt03_16 * (sh1[2][1] * sh4[2][1] - sh1[0][1] * sh4[6][1]);
        sh5[4][3] = sqrtf(8.0 / 7.0) * sh1[1][1] * sh4[3][2] + sqrtf(5.0 / 7.0) * sh1[0][1] * sh4[4][2] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh4[2][2] - sh1[0][1] * sh4[6][2]);
        sh5[4][4] = sh1[1][1] * sh4[3][3] + kSqrt05_08 * sh1[0][1] * sh4[4][3] - kSqrt01_08 * (sh1[2][1] * sh4[2][3] - sh1[0][1] * sh4[6][3]);
        sh5[4][5] = sqrtf(24.0 / 25.0) * sh1[1][1] * sh4[3][4] + kSqrt03_05 * sh1[0][1] * sh4[4][4] - kSqrt03_25 * (sh1[2][1] * sh4[2][4] - sh1[0][1] * sh4[6][4]);
        sh5[4][6] = sh1[1][1] * sh4[3][5] + kSqrt05_08 * sh1[0][1] * sh4[4][5] - kSqrt01_08 * (sh1[2][1] * sh4[2][5] - sh1[0][1] * sh4[6][5]);
        sh5[4][7] = sqrtf(8.0 / 7.0) * sh1[1][1] * sh4[3][6] + sqrtf(5.0 / 7.0) * sh1[0][1] * sh4[4][6] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh4[2][6] - sh1[0][1] * sh4[6][6]);
        sh5[4][8] = kSqrt03_02 * sh1[1][1] * sh4[3][7] + kSqrt15_16 * sh1[0][1] * sh4[4][7] - kSqrt03_16 * (sh1[2][1] * sh4[2][7] - sh1[0][1] * sh4[6][7]);
        sh5[4][9] = sqrtf(8.0 / 3.0) * sh1[1][1] * sh4[3][8] + sqrtf(5.0 / 3.0) * sh1[0][1] * sh4[4][8] - kSqrt01_03 * (sh1[2][1] * sh4[2][8] - sh1[0][1] * sh4[6][8]);
        sh5[4][10] = kSqrt04_15 * (sh1[1][2] * sh4[3][8] - sh1[1][0] * sh4[3][0]) + kSqrt01_06 * (sh1[0][2] * sh4[4][8] - sh1[0][0] * sh4[4][0]) - kSqrt01_30 * ((sh1[2][2] * sh4[2][8] - sh1[2][0] * sh4[2][0]) - (sh1[0][2] * sh4[6][8] - sh1[0][0] * sh4[6][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[4]);

        sh5[5][0] = sqrtf(5.0 / 18.0) * (sh1[1][2] * sh4[4][0] + sh1[1][0] * sh4[4][8]) - sqrtf(1.0 / 9.0) * ((sh1[2][2] * sh4[5][0] + sh1[2][0] * sh4[5][8]) + (sh1[0][2] * sh4[3][0] + sh1[0][0] * sh4[3][8]));
        sh5[5][1] = sqrtf(25.0 / 9.0) * sh1[1][1] * sh4[4][0] - sqrtf(10.0 / 9.0) * (sh1[2][1] * sh4[5][0] + sh1[0][1] * sh4[3][0]);
        sh5[5][2] = sqrtf(25.0 / 16.0) * sh1[1][1] * sh4[4][1] - kSqrt05_08 * (sh1[2][1] * sh4[5][1] + sh1[0][1] * sh4[3][1]);
        sh5[5][3] = sqrtf(25.0 / 21.0) * sh1[1][1] * sh4[4][2] - sqrtf(10.0 / 21.0) * (sh1[2][1] * sh4[5][2] + sh1[0][1] * sh4[3][2]);
        sh5[5][4] = sqrtf(25.0 / 24.0) * sh1[1][1] * sh4[4][3] - sqrtf(5.0 / 12.0) * (sh1[2][1] * sh4[5][3] + sh1[0][1] * sh4[3][3]);
        sh5[5][5] = sh1[1][1] * sh4[4][4] - kSqrt02_05 * (sh1[2][1] * sh4[5][4] + sh1[0][1] * sh4[3][4]);
        sh5[5][6] = sqrtf(25.0 / 24.0) * sh1[1][1] * sh4[4][5] - sqrtf(5.0 / 12.0) * (sh1[2][1] * sh4[5][5] + sh1[0][1] * sh4[3][5]);
        sh5[5][7] = sqrtf(25.0 / 21.0) * sh1[1][1] * sh4[4][6] - sqrtf(10.0 / 21.0) * (sh1[2][1] * sh4[5][6] + sh1[0][1] * sh4[3][6]);
        sh5[5][8] = sqrtf(25.0 / 16.0) * sh1[1][1] * sh4[4][7] - kSqrt05_08 * (sh1[2][1] * sh4[5][7] + sh1[0][1] * sh4[3][7]);
        sh5[5][9] = sqrtf(25.0 / 9.0) * sh1[1][1] * sh4[4][8] - sqrtf(10.0 / 9.0) * (sh1[2][1] * sh4[5][8] + sh1[0][1] * sh4[3][8]);
        sh5[5][10] = sqrtf(5.0 / 18.0) * (sh1[1][2] * sh4[4][8] - sh1[1][0] * sh4[4][0]) - sqrtf(1.0 / 9.0) * ((sh1[2][2] * sh4[5][8] - sh1[2][0] * sh4[5][0]) + (sh1[0][2] * sh4[3][8] - sh1[0][0] * sh4[3][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[5]);

        sh5[6][0] = kSqrt04_15 * (sh1[1][2] * sh4[5][0] + sh1[1][0] * sh4[5][8]) + kSqrt01_06 * (sh1[2][2] * sh4[4][0] + sh1[2][0] * sh4[4][8]) - kSqrt01_30 * ((sh1[2][2] * sh4[6][0] + sh1[2][0] * sh4[6][8]) + (sh1[0][2] * sh4[2][0] + sh1[0][0] * sh4[2][8]));
        sh5[6][1] = sqrtf(8.0 / 3.0) * sh1[1][1] * sh4[5][0] + sqrtf(5.0 / 3.0) * sh1[2][1] * sh4[4][0] - kSqrt01_03 * (sh1[2][1] * sh4[6][0] + sh1[0][1] * sh4[2][0]);
        sh5[6][2] = kSqrt03_02 * sh1[1][1] * sh4[5][1] + kSqrt15_16 * sh1[2][1] * sh4[4][1] - kSqrt03_16 * (sh1[2][1] * sh4[6][1] + sh1[0][1] * sh4[2][1]);
        sh5[6][3] = sqrtf(8.0 / 7.0) * sh1[1][1] * sh4[5][2] + sqrtf(5.0 / 7.0) * sh1[2][1] * sh4[4][2] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh4[6][2] + sh1[0][1] * sh4[2][2]);
        sh5[6][4] = sh1[1][1] * sh4[5][3] + kSqrt05_08 * sh1[2][1] * sh4[4][3] - kSqrt01_08 * (sh1[2][1] * sh4[6][3] + sh1[0][1] * sh4[2][3]);
        sh5[6][5] = sqrtf(24.0 / 25.0) * sh1[1][1] * sh4[5][4] + kSqrt03_05 * sh1[2][1] * sh4[4][4] - kSqrt03_25 * (sh1[2][1] * sh4[6][4] + sh1[0][1] * sh4[2][4]);
        sh5[6][6] = sh1[1][1] * sh4[5][5] + kSqrt05_08 * sh1[2][1] * sh4[4][5] - kSqrt01_08 * (sh1[2][1] * sh4[6][5] + sh1[0][1] * sh4[2][5]);
        sh5[6][7] = sqrtf(8.0 / 7.0) * sh1[1][1] * sh4[5][6] + sqrtf(5.0 / 7.0) * sh1[2][1] * sh4[4][6] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh4[6][6] + sh1[0][1] * sh4[2][6]);
        sh5[6][8] = kSqrt03_02 * sh1[1][1] * sh4[5][7] + kSqrt15_16 * sh1[2][1] * sh4[4][7] - kSqrt03_16 * (sh1[2][1] * sh4[6][7] + sh1[0][1] * sh4[2][7]);
        sh5[6][9] = sqrtf(8.0 / 3.0) * sh1[1][1] * sh4[5][8] + sqrtf(5.0 / 3.0) * sh1[2][1] * sh4[4][8] - kSqrt01_03 * (sh1[2][1] * sh4[6][8] + sh1[0][1] * sh4[2][8]);
        sh5[6][10] = kSqrt04_15 * (sh1[1][2] * sh4[5][8] - sh1[1][0] * sh4[5][0]) + kSqrt01_06 * (sh1[2][2] * sh4[4][8] - sh1[2][0] * sh4[4][0]) - kSqrt01_30 * ((sh1[2][2] * sh4[6][8] - sh1[2][0] * sh4[6][0]) + (sh1[0][2] * sh4[2][8] - sh1[0][0] * sh4[2][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[6]);

        sh5[7][0] = sqrtf(7.0 / 30.0) * (sh1[1][2] * sh4[6][0] + sh1[1][0] * sh4[6][8]) + sqrtf(7.0 / 60.0) * ((sh1[2][2] * sh4[5][0] + sh1[2][0] * sh4[5][8]) - (sh1[0][2] * sh4[3][0] + sh1[0][0] * sh4[3][8])) - kSqrt01_60 * ((sh1[2][2] * sh4[7][0] + sh1[2][0] * sh4[7][8]) + (sh1[0][2] * sh4[1][0] + sh1[0][0] * sh4[1][8]));
        sh5[7][1] = sqrtf(7.0 / 3.0) * sh1[1][1] * sh4[6][0] + kSqrt07_06 * (sh1[2][1] * sh4[5][0] - sh1[0][1] * sh4[3][0]) - kSqrt01_06 * (sh1[2][1] * sh4[7][0] + sh1[0][1] * sh4[1][0]);
        sh5[7][2] = sqrtf(21.0 / 16.0) * sh1[1][1] * sh4[6][1] + kSqrt21_32 * (sh1[2][1] * sh4[5][1] - sh1[0][1] * sh4[3][1]) - kSqrt03_32 * (sh1[2][1] * sh4[7][1] + sh1[0][1] * sh4[1][1]);
        sh5[7][3] = sh1[1][1] * sh4[6][2] + kSqrt01_02 * (sh1[2][1] * sh4[5][2] - sh1[0][1] * sh4[3][2]) - kSqrt01_14 * (sh1[2][1] * sh4[7][2] + sh1[0][1] * sh4[1][2]);
        sh5[7][4] = kSqrt07_08 * sh1[1][1] * sh4[6][3] + kSqrt07_16 * (sh1[2][1] * sh4[5][3] - sh1[0][1] * sh4[3][3]) - kSqrt01_16 * (sh1[2][1] * sh4[7][3] + sh1[0][1] * sh4[1][3]);
        sh5[7][5] = sqrtf(21.0 / 25.0) * sh1[1][1] * sh4[6][4] + kSqrt21_50 * (sh1[2][1] * sh4[5][4] - sh1[0][1] * sh4[3][4]) - kSqrt03_50 * (sh1[2][1] * sh4[7][4] + sh1[0][1] * sh4[1][4]);
        sh5[7][6] = kSqrt07_08 * sh1[1][1] * sh4[6][5] + kSqrt07_16 * (sh1[2][1] * sh4[5][5] - sh1[0][1] * sh4[3][5]) - kSqrt01_16 * (sh1[2][1] * sh4[7][5] + sh1[0][1] * sh4[1][5]);
        sh5[7][7] = sh1[1][1] * sh4[6][6] + kSqrt01_02 * (sh1[2][1] * sh4[5][6] - sh1[0][1] * sh4[3][6]) - kSqrt01_14 * (sh1[2][1] * sh4[7][6] + sh1[0][1] * sh4[1][6]);
        sh5[7][8] = sqrtf(21.0 / 16.0) * sh1[1][1] * sh4[6][7] + kSqrt21_32 * (sh1[2][1] * sh4[5][7] - sh1[0][1] * sh4[3][7]) - kSqrt03_32 * (sh1[2][1] * sh4[7][7] + sh1[0][1] * sh4[1][7]);
        sh5[7][9] = sqrtf(7.0 / 3.0) * sh1[1][1] * sh4[6][8] + kSqrt07_06 * (sh1[2][1] * sh4[5][8] - sh1[0][1] * sh4[3][8]) - kSqrt01_06 * (sh1[2][1] * sh4[7][8] + sh1[0][1] * sh4[1][8]);
        sh5[7][10] = sqrtf(7.0 / 30.0) * (sh1[1][2] * sh4[6][8] - sh1[1][0] * sh4[6][0]) + sqrtf(7.0 / 60.0) * ((sh1[2][2] * sh4[5][8] - sh1[2][0] * sh4[5][0]) - (sh1[0][2] * sh4[3][8] - sh1[0][0] * sh4[3][0])) - kSqrt01_60 * ((sh1[2][2] * sh4[7][8] - sh1[2][0] * sh4[7][0]) + (sh1[0][2] * sh4[1][8] - sh1[0][0] * sh4[1][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[7]);

        sh5[8][0] = sqrtf(8.0 / 45.0) * (sh1[1][2] * sh4[7][0] + sh1[1][0] * sh4[7][8]) + sqrtf(7.0 / 45.0) * ((sh1[2][2] * sh4[6][0] + sh1[2][0] * sh4[6][8]) - (sh1[0][2] * sh4[2][0] + sh1[0][0] * sh4[2][8])) - sqrtf(1.0 / 180.0) * ((sh1[2][2] * sh4[8][0] + sh1[2][0] * sh4[8][8]) + (sh1[0][2] * sh4[0][0] + sh1[0][0] * sh4[0][8]));
        sh5[8][1] = sqrtf(16.0 / 9.0) * sh1[1][1] * sh4[7][0] + sqrtf(14.0 / 9.0) * (sh1[2][1] * sh4[6][0] - sh1[0][1] * sh4[2][0]) - kSqrt01_18 * (sh1[2][1] * sh4[8][0] + sh1[0][1] * sh4[0][0]);
        sh5[8][2] = sh1[1][1] * sh4[7][1] + kSqrt07_08 * (sh1[2][1] * sh4[6][1] - sh1[0][1] * sh4[2][1]) - kSqrt01_32 * (sh1[2][1] * sh4[8][1] + sh1[0][1] * sh4[0][1]);
        sh5[8][3] = sqrtf(16.0 / 21.0) * sh1[1][1] * sh4[7][2] + kSqrt02_03 * (sh1[2][1] * sh4[6][2] - sh1[0][1] * sh4[2][2]) - sqrtf(1.0 / 42.0) * (sh1[2][1] * sh4[8][2] + sh1[0][1] * sh4[0][2]);
        sh5[8][4] = kSqrt02_03 * sh1[1][1] * sh4[7][3] + kSqrt07_12 * (sh1[2][1] * sh4[6][3] - sh1[0][1] * sh4[2][3]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh4[8][3] + sh1[0][1] * sh4[0][3]);
        sh5[8][5] = sqrtf(16.0 / 25.0) * sh1[1][1] * sh4[7][4] + kSqrt14_25 * (sh1[2][1] * sh4[6][4] - sh1[0][1] * sh4[2][4]) - kSqrt01_50 * (sh1[2][1] * sh4[8][4] + sh1[0][1] * sh4[0][4]);
        sh5[8][6] = kSqrt02_03 * sh1[1][1] * sh4[7][5] + kSqrt07_12 * (sh1[2][1] * sh4[6][5] - sh1[0][1] * sh4[2][5]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh4[8][5] + sh1[0][1] * sh4[0][5]);
        sh5[8][7] = sqrtf(16.0 / 21.0) * sh1[1][1] * sh4[7][6] + kSqrt02_03 * (sh1[2][1] * sh4[6][6] - sh1[0][1] * sh4[2][6]) - sqrtf(1.0 / 42.0) * (sh1[2][1] * sh4[8][6] + sh1[0][1] * sh4[0][6]);
        sh5[8][8] = sh1[1][1] * sh4[7][7] + kSqrt07_08 * (sh1[2][1] * sh4[6][7] - sh1[0][1] * sh4[2][7]) - kSqrt01_32 * (sh1[2][1] * sh4[8][7] + sh1[0][1] * sh4[0][7]);
        sh5[8][9] = sqrtf(16.0 / 9.0) * sh1[1][1] * sh4[7][8] + sqrtf(14.0 / 9.0) * (sh1[2][1] * sh4[6][8] - sh1[0][1] * sh4[2][8]) - kSqrt01_18 * (sh1[2][1] * sh4[8][8] + sh1[0][1] * sh4[0][8]);
        sh5[8][10] = sqrtf(8.0 / 45.0) * (sh1[1][2] * sh4[7][8] - sh1[1][0] * sh4[7][0]) + sqrtf(7.0 / 45.0) * ((sh1[2][2] * sh4[6][8] - sh1[2][0] * sh4[6][0]) - (sh1[0][2] * sh4[2][8] - sh1[0][0] * sh4[2][0])) - sqrtf(1.0 / 180.0) * ((sh1[2][2] * sh4[8][8] - sh1[2][0] * sh4[8][0]) + (sh1[0][2] * sh4[0][8] - sh1[0][0] * sh4[0][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[8]);

        sh5[9][0] = kSqrt01_10 * (sh1[1][2] * sh4[8][0] + sh1[1][0] * sh4[8][8]) + kSqrt01_05 * ((sh1[2][2] * sh4[7][0] + sh1[2][0] * sh4[7][8]) - (sh1[0][2] * sh4[1][0] + sh1[0][0] * sh4[1][8]));
        sh5[9][1] = sh1[1][1] * sh4[8][0] + kSqrt02_01 * (sh1[2][1] * sh4[7][0] - sh1[0][1] * sh4[1][0]);
        sh5[9][2] = sqrtf(9.0 / 16.0) * sh1[1][1] * sh4[8][1] + kSqrt09_08 * (sh1[2][1] * sh4[7][1] - sh1[0][1] * sh4[1][1]);
        sh5[9][3] = sqrtf(3.0 / 7.0) * sh1[1][1] * sh4[8][2] + kSqrt06_07 * (sh1[2][1] * sh4[7][2] - sh1[0][1] * sh4[1][2]);
        sh5[9][4] = kSqrt03_08 * sh1[1][1] * sh4[8][3] + kSqrt03_04 * (sh1[2][1] * sh4[7][3] - sh1[0][1] * sh4[1][3]);
        sh5[9][5] = kSqrt09_25 * sh1[1][1] * sh4[8][4] + kSqrt18_25 * (sh1[2][1] * sh4[7][4] - sh1[0][1] * sh4[1][4]);
        sh5[9][6] = kSqrt03_08 * sh1[1][1] * sh4[8][5] + kSqrt03_04 * (sh1[2][1] * sh4[7][5] - sh1[0][1] * sh4[1][5]);
        sh5[9][7] = sqrtf(3.0 / 7.0) * sh1[1][1] * sh4[8][6] + kSqrt06_07 * (sh1[2][1] * sh4[7][6] - sh1[0][1] * sh4[1][6]);
        sh5[9][8] = sqrtf(9.0 / 16.0) * sh1[1][1] * sh4[8][7] + kSqrt09_08 * (sh1[2][1] * sh4[7][7] - sh1[0][1] * sh4[1][7]);
        sh5[9][9] = sh1[1][1] * sh4[8][8] + kSqrt02_01 * (sh1[2][1] * sh4[7][8] - sh1[0][1] * sh4[1][8]);
        sh5[9][10] = kSqrt01_10 * (sh1[1][2] * sh4[8][8] - sh1[1][0] * sh4[8][0]) + kSqrt01_05 * ((sh1[2][2] * sh4[7][8] - sh1[2][0] * sh4[7][0]) - (sh1[0][2] * sh4[1][8] - sh1[0][0] * sh4[1][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[9]);

        sh5[10][0] = kSqrt01_04 * ((sh1[2][2] * sh4[8][0] + sh1[2][0] * sh4[8][8]) - (sh1[0][2] * sh4[0][0] + sh1[0][0] * sh4[0][8]));
        sh5[10][1] = sqrtf(5.0 / 2.0) * (sh1[2][1] * sh4[8][0] - sh1[0][1] * sh4[0][0]);
        sh5[10][2] = sqrtf(45.0 / 32.0) * (sh1[2][1] * sh4[8][1] - sh1[0][1] * sh4[0][1]);
        sh5[10][3] = kSqrt15_14 * (sh1[2][1] * sh4[8][2] - sh1[0][1] * sh4[0][2]);
        sh5[10][4] = kSqrt15_16 * (sh1[2][1] * sh4[8][3] - sh1[0][1] * sh4[0][3]);
        sh5[10][5] = kSqrt09_10 * (sh1[2][1] * sh4[8][4] - sh1[0][1] * sh4[0][4]);
        sh5[10][6] = kSqrt15_16 * (sh1[2][1] * sh4[8][5] - sh1[0][1] * sh4[0][5]);
        sh5[10][7] = kSqrt15_14 * (sh1[2][1] * sh4[8][6] - sh1[0][1] * sh4[0][6]);
        sh5[10][8] = sqrtf(45.0 / 32.0) * (sh1[2][1] * sh4[8][7] - sh1[0][1] * sh4[0][7]);
        sh5[10][9] = sqrtf(5.0 / 2.0) * (sh1[2][1] * sh4[8][8] - sh1[0][1] * sh4[0][8]);
        sh5[10][10] = kSqrt01_04 * ((sh1[2][2] * sh4[8][8] - sh1[2][0] * sh4[8][0]) - (sh1[0][2] * sh4[0][8] - sh1[0][0] * sh4[0][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[10]);

        if (n < 7)
            return;

        coeffsIn += 11;
        float sh6[13][13];

        sh6[0][0] = kSqrt01_04 * ((sh1[2][2] * sh5[0][0] + sh1[2][0] * sh5[0][10]) + (sh1[0][2] * sh5[10][0] + sh1[0][0] * sh5[10][10]));
        sh6[0][1] = sqrtf(3.0 / 1.0) * (sh1[2][1] * sh5[0][0] + sh1[0][1] * sh5[10][0]);
        sh6[0][2] = sqrtf(33.0 / 20.0) * (sh1[2][1] * sh5[0][1] + sh1[0][1] * sh5[10][1]);
        sh6[0][3] = sqrtf(11.0 / 9.0) * (sh1[2][1] * sh5[0][2] + sh1[0][1] * sh5[10][2]);
        sh6[0][4] = sqrtf(33.0 / 32.0) * (sh1[2][1] * sh5[0][3] + sh1[0][1] * sh5[10][3]);
        sh6[0][5] = sqrtf(33.0 / 35.0) * (sh1[2][1] * sh5[0][4] + sh1[0][1] * sh5[10][4]);
        sh6[0][6] = sqrtf(11.0 / 12.0) * (sh1[2][1] * sh5[0][5] + sh1[0][1] * sh5[10][5]);
        sh6[0][7] = sqrtf(33.0 / 35.0) * (sh1[2][1] * sh5[0][6] + sh1[0][1] * sh5[10][6]);
        sh6[0][8] = sqrtf(33.0 / 32.0) * (sh1[2][1] * sh5[0][7] + sh1[0][1] * sh5[10][7]);
        sh6[0][9] = sqrtf(11.0 / 9.0) * (sh1[2][1] * sh5[0][8] + sh1[0][1] * sh5[10][8]);
        sh6[0][10] = sqrtf(33.0 / 20.0) * (sh1[2][1] * sh5[0][9] + sh1[0][1] * sh5[10][9]);
        sh6[0][11] = sqrtf(3.0 / 1.0) * (sh1[2][1] * sh5[0][10] + sh1[0][1] * sh5[10][10]);
        sh6[0][12] = kSqrt01_04 * ((sh1[2][2] * sh5[0][10] - sh1[2][0] * sh5[0][0]) + (sh1[0][2] * sh5[10][10] - sh1[0][0] * sh5[10][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[0]);

        sh6[1][0] = kSqrt01_12 * (sh1[1][2] * sh5[0][0] + sh1[1][0] * sh5[0][10]) + sqrtf(5.0 / 24.0) * ((sh1[2][2] * sh5[1][0] + sh1[2][0] * sh5[1][10]) + (sh1[0][2] * sh5[9][0] + sh1[0][0] * sh5[9][10]));
        sh6[1][1] = sh1[1][1] * sh5[0][0] + sqrtf(5.0 / 2.0) * (sh1[2][1] * sh5[1][0] + sh1[0][1] * sh5[9][0]);
        sh6[1][2] = sqrtf(11.0 / 20.0) * sh1[1][1] * sh5[0][1] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh5[1][1] + sh1[0][1] * sh5[9][1]);
        sh6[1][3] = sqrtf(11.0 / 27.0) * sh1[1][1] * sh5[0][2] + sqrtf(55.0 / 54.0) * (sh1[2][1] * sh5[1][2] + sh1[0][1] * sh5[9][2]);
        sh6[1][4] = sqrtf(11.0 / 32.0) * sh1[1][1] * sh5[0][3] + sqrtf(55.0 / 64.0) * (sh1[2][1] * sh5[1][3] + sh1[0][1] * sh5[9][3]);
        sh6[1][5] = sqrtf(11.0 / 35.0) * sh1[1][1] * sh5[0][4] + sqrtf(11.0 / 14.0) * (sh1[2][1] * sh5[1][4] + sh1[0][1] * sh5[9][4]);
        sh6[1][6] = sqrtf(11.0 / 36.0) * sh1[1][1] * sh5[0][5] + sqrtf(55.0 / 72.0) * (sh1[2][1] * sh5[1][5] + sh1[0][1] * sh5[9][5]);
        sh6[1][7] = sqrtf(11.0 / 35.0) * sh1[1][1] * sh5[0][6] + sqrtf(11.0 / 14.0) * (sh1[2][1] * sh5[1][6] + sh1[0][1] * sh5[9][6]);
        sh6[1][8] = sqrtf(11.0 / 32.0) * sh1[1][1] * sh5[0][7] + sqrtf(55.0 / 64.0) * (sh1[2][1] * sh5[1][7] + sh1[0][1] * sh5[9][7]);
        sh6[1][9] = sqrtf(11.0 / 27.0) * sh1[1][1] * sh5[0][8] + sqrtf(55.0 / 54.0) * (sh1[2][1] * sh5[1][8] + sh1[0][1] * sh5[9][8]);
        sh6[1][10] = sqrtf(11.0 / 20.0) * sh1[1][1] * sh5[0][9] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh5[1][9] + sh1[0][1] * sh5[9][9]);
        sh6[1][11] = sh1[1][1] * sh5[0][10] + sqrtf(5.0 / 2.0) * (sh1[2][1] * sh5[1][10] + sh1[0][1] * sh5[9][10]);
        sh6[1][12] = kSqrt01_12 * (sh1[1][2] * sh5[0][10] - sh1[1][0] * sh5[0][0]) + sqrtf(5.0 / 24.0) * ((sh1[2][2] * sh5[1][10] - sh1[2][0] * sh5[1][0]) + (sh1[0][2] * sh5[9][10] - sh1[0][0] * sh5[9][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[1]);

        sh6[2][0] = sqrtf(5.0 / 33.0) * (sh1[1][2] * sh5[1][0] + sh1[1][0] * sh5[1][10]) + sqrtf(15.0 / 88.0) * ((sh1[2][2] * sh5[2][0] + sh1[2][0] * sh5[2][10]) + (sh1[0][2] * sh5[8][0] + sh1[0][0] * sh5[8][10])) - sqrtf(1.0 / 264.0) * ((sh1[2][2] * sh5[0][0] + sh1[2][0] * sh5[0][10]) - (sh1[0][2] * sh5[10][0] + sh1[0][0] * sh5[10][10]));
        sh6[2][1] = sqrtf(20.0 / 11.0) * sh1[1][1] * sh5[1][0] + sqrtf(45.0 / 22.0) * (sh1[2][1] * sh5[2][0] + sh1[0][1] * sh5[8][0]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh5[0][0] - sh1[0][1] * sh5[10][0]);
        sh6[2][2] = sh1[1][1] * sh5[1][1] + kSqrt09_08 * (sh1[2][1] * sh5[2][1] + sh1[0][1] * sh5[8][1]) - sqrtf(1.0 / 40.0) * (sh1[2][1] * sh5[0][1] - sh1[0][1] * sh5[10][1]);
        sh6[2][3] = sqrtf(20.0 / 27.0) * sh1[1][1] * sh5[1][2] + kSqrt05_06 * (sh1[2][1] * sh5[2][2] + sh1[0][1] * sh5[8][2]) - sqrtf(1.0 / 54.0) * (sh1[2][1] * sh5[0][2] - sh1[0][1] * sh5[10][2]);
        sh6[2][4] = kSqrt05_08 * sh1[1][1] * sh5[1][3] + sqrtf(45.0 / 64.0) * (sh1[2][1] * sh5[2][3] + sh1[0][1] * sh5[8][3]) - sqrtf(1.0 / 64.0) * (sh1[2][1] * sh5[0][3] - sh1[0][1] * sh5[10][3]);
        sh6[2][5] = sqrtf(4.0 / 7.0) * sh1[1][1] * sh5[1][4] + sqrtf(9.0 / 14.0) * (sh1[2][1] * sh5[2][4] + sh1[0][1] * sh5[8][4]) - sqrtf(1.0 / 70.0) * (sh1[2][1] * sh5[0][4] - sh1[0][1] * sh5[10][4]);
        sh6[2][6] = kSqrt05_09 * sh1[1][1] * sh5[1][5] + kSqrt05_08 * (sh1[2][1] * sh5[2][5] + sh1[0][1] * sh5[8][5]) - sqrtf(1.0 / 72.0) * (sh1[2][1] * sh5[0][5] - sh1[0][1] * sh5[10][5]);
        sh6[2][7] = sqrtf(4.0 / 7.0) * sh1[1][1] * sh5[1][6] + sqrtf(9.0 / 14.0) * (sh1[2][1] * sh5[2][6] + sh1[0][1] * sh5[8][6]) - sqrtf(1.0 / 70.0) * (sh1[2][1] * sh5[0][6] - sh1[0][1] * sh5[10][6]);
        sh6[2][8] = kSqrt05_08 * sh1[1][1] * sh5[1][7] + sqrtf(45.0 / 64.0) * (sh1[2][1] * sh5[2][7] + sh1[0][1] * sh5[8][7]) - sqrtf(1.0 / 64.0) * (sh1[2][1] * sh5[0][7] - sh1[0][1] * sh5[10][7]);
        sh6[2][9] = sqrtf(20.0 / 27.0) * sh1[1][1] * sh5[1][8] + kSqrt05_06 * (sh1[2][1] * sh5[2][8] + sh1[0][1] * sh5[8][8]) - sqrtf(1.0 / 54.0) * (sh1[2][1] * sh5[0][8] - sh1[0][1] * sh5[10][8]);
        sh6[2][10] = sh1[1][1] * sh5[1][9] + kSqrt09_08 * (sh1[2][1] * sh5[2][9] + sh1[0][1] * sh5[8][9]) - sqrtf(1.0 / 40.0) * (sh1[2][1] * sh5[0][9] - sh1[0][1] * sh5[10][9]);
        sh6[2][11] = sqrtf(20.0 / 11.0) * sh1[1][1] * sh5[1][10] + sqrtf(45.0 / 22.0) * (sh1[2][1] * sh5[2][10] + sh1[0][1] * sh5[8][10]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh5[0][10] - sh1[0][1] * sh5[10][10]);
        sh6[2][12] = sqrtf(5.0 / 33.0) * (sh1[1][2] * sh5[1][10] - sh1[1][0] * sh5[1][0]) + sqrtf(15.0 / 88.0) * ((sh1[2][2] * sh5[2][10] - sh1[2][0] * sh5[2][0]) + (sh1[0][2] * sh5[8][10] - sh1[0][0] * sh5[8][0])) - sqrtf(1.0 / 264.0) * ((sh1[2][2] * sh5[0][10] - sh1[2][0] * sh5[0][0]) - (sh1[0][2] * sh5[10][10] - sh1[0][0] * sh5[10][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[2]);

        sh6[3][0] = sqrtf(9.0 / 44.0) * (sh1[1][2] * sh5[2][0] + sh1[1][0] * sh5[2][10]) + sqrtf(3.0 / 22.0) * ((sh1[2][2] * sh5[3][0] + sh1[2][0] * sh5[3][10]) + (sh1[0][2] * sh5[7][0] + sh1[0][0] * sh5[7][10])) - sqrtf(1.0 / 88.0) * ((sh1[2][2] * sh5[1][0] + sh1[2][0] * sh5[1][10]) - (sh1[0][2] * sh5[9][0] + sh1[0][0] * sh5[9][10]));
        sh6[3][1] = sqrtf(27.0 / 11.0) * sh1[1][1] * sh5[2][0] + sqrtf(18.0 / 11.0) * (sh1[2][1] * sh5[3][0] + sh1[0][1] * sh5[7][0]) - sqrtf(3.0 / 22.0) * (sh1[2][1] * sh5[1][0] - sh1[0][1] * sh5[9][0]);
        sh6[3][2] = sqrtf(27.0 / 20.0) * sh1[1][1] * sh5[2][1] + kSqrt09_10 * (sh1[2][1] * sh5[3][1] + sh1[0][1] * sh5[7][1]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh5[1][1] - sh1[0][1] * sh5[9][1]);
        sh6[3][3] = sh1[1][1] * sh5[2][2] + kSqrt02_03 * (sh1[2][1] * sh5[3][2] + sh1[0][1] * sh5[7][2]) - kSqrt01_18 * (sh1[2][1] * sh5[1][2] - sh1[0][1] * sh5[9][2]);
        sh6[3][4] = sqrtf(27.0 / 32.0) * sh1[1][1] * sh5[2][3] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh5[3][3] + sh1[0][1] * sh5[7][3]) - sqrtf(3.0 / 64.0) * (sh1[2][1] * sh5[1][3] - sh1[0][1] * sh5[9][3]);
        sh6[3][5] = sqrtf(27.0 / 35.0) * sh1[1][1] * sh5[2][4] + sqrtf(18.0 / 35.0) * (sh1[2][1] * sh5[3][4] + sh1[0][1] * sh5[7][4]) - sqrtf(3.0 / 70.0) * (sh1[2][1] * sh5[1][4] - sh1[0][1] * sh5[9][4]);
        sh6[3][6] = kSqrt03_04 * sh1[1][1] * sh5[2][5] + kSqrt01_02 * (sh1[2][1] * sh5[3][5] + sh1[0][1] * sh5[7][5]) - kSqrt01_24 * (sh1[2][1] * sh5[1][5] - sh1[0][1] * sh5[9][5]);
        sh6[3][7] = sqrtf(27.0 / 35.0) * sh1[1][1] * sh5[2][6] + sqrtf(18.0 / 35.0) * (sh1[2][1] * sh5[3][6] + sh1[0][1] * sh5[7][6]) - sqrtf(3.0 / 70.0) * (sh1[2][1] * sh5[1][6] - sh1[0][1] * sh5[9][6]);
        sh6[3][8] = sqrtf(27.0 / 32.0) * sh1[1][1] * sh5[2][7] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh5[3][7] + sh1[0][1] * sh5[7][7]) - sqrtf(3.0 / 64.0) * (sh1[2][1] * sh5[1][7] - sh1[0][1] * sh5[9][7]);
        sh6[3][9] = sh1[1][1] * sh5[2][8] + kSqrt02_03 * (sh1[2][1] * sh5[3][8] + sh1[0][1] * sh5[7][8]) - kSqrt01_18 * (sh1[2][1] * sh5[1][8] - sh1[0][1] * sh5[9][8]);
        sh6[3][10] = sqrtf(27.0 / 20.0) * sh1[1][1] * sh5[2][9] + kSqrt09_10 * (sh1[2][1] * sh5[3][9] + sh1[0][1] * sh5[7][9]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh5[1][9] - sh1[0][1] * sh5[9][9]);
        sh6[3][11] = sqrtf(27.0 / 11.0) * sh1[1][1] * sh5[2][10] + sqrtf(18.0 / 11.0) * (sh1[2][1] * sh5[3][10] + sh1[0][1] * sh5[7][10]) - sqrtf(3.0 / 22.0) * (sh1[2][1] * sh5[1][10] - sh1[0][1] * sh5[9][10]);
        sh6[3][12] = sqrtf(9.0 / 44.0) * (sh1[1][2] * sh5[2][10] - sh1[1][0] * sh5[2][0]) + sqrtf(3.0 / 22.0) * ((sh1[2][2] * sh5[3][10] - sh1[2][0] * sh5[3][0]) + (sh1[0][2] * sh5[7][10] - sh1[0][0] * sh5[7][0])) - sqrtf(1.0 / 88.0) * ((sh1[2][2] * sh5[1][10] - sh1[2][0] * sh5[1][0]) - (sh1[0][2] * sh5[9][10] - sh1[0][0] * sh5[9][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[3]);

        sh6[4][0] = sqrtf(8.0 / 33.0) * (sh1[1][2] * sh5[3][0] + sh1[1][0] * sh5[3][10]) + sqrtf(7.0 / 66.0) * ((sh1[2][2] * sh5[4][0] + sh1[2][0] * sh5[4][10]) + (sh1[0][2] * sh5[6][0] + sh1[0][0] * sh5[6][10])) - sqrtf(1.0 / 44.0) * ((sh1[2][2] * sh5[2][0] + sh1[2][0] * sh5[2][10]) - (sh1[0][2] * sh5[8][0] + sh1[0][0] * sh5[8][10]));
        sh6[4][1] = sqrtf(32.0 / 11.0) * sh1[1][1] * sh5[3][0] + sqrtf(14.0 / 11.0) * (sh1[2][1] * sh5[4][0] + sh1[0][1] * sh5[6][0]) - sqrtf(3.0 / 11.0) * (sh1[2][1] * sh5[2][0] - sh1[0][1] * sh5[8][0]);
        sh6[4][2] = kSqrt08_05 * sh1[1][1] * sh5[3][1] + kSqrt07_10 * (sh1[2][1] * sh5[4][1] + sh1[0][1] * sh5[6][1]) - sqrtf(3.0 / 20.0) * (sh1[2][1] * sh5[2][1] - sh1[0][1] * sh5[8][1]);
        sh6[4][3] = sqrtf(32.0 / 27.0) * sh1[1][1] * sh5[3][2] + sqrtf(14.0 / 27.0) * (sh1[2][1] * sh5[4][2] + sh1[0][1] * sh5[6][2]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh5[2][2] - sh1[0][1] * sh5[8][2]);
        sh6[4][4] = sh1[1][1] * sh5[3][3] + kSqrt07_16 * (sh1[2][1] * sh5[4][3] + sh1[0][1] * sh5[6][3]) - kSqrt03_32 * (sh1[2][1] * sh5[2][3] - sh1[0][1] * sh5[8][3]);
        sh6[4][5] = sqrtf(32.0 / 35.0) * sh1[1][1] * sh5[3][4] + kSqrt02_05 * (sh1[2][1] * sh5[4][4] + sh1[0][1] * sh5[6][4]) - sqrtf(3.0 / 35.0) * (sh1[2][1] * sh5[2][4] - sh1[0][1] * sh5[8][4]);
        sh6[4][6] = kSqrt08_09 * sh1[1][1] * sh5[3][5] + sqrtf(7.0 / 18.0) * (sh1[2][1] * sh5[4][5] + sh1[0][1] * sh5[6][5]) - kSqrt01_12 * (sh1[2][1] * sh5[2][5] - sh1[0][1] * sh5[8][5]);
        sh6[4][7] = sqrtf(32.0 / 35.0) * sh1[1][1] * sh5[3][6] + kSqrt02_05 * (sh1[2][1] * sh5[4][6] + sh1[0][1] * sh5[6][6]) - sqrtf(3.0 / 35.0) * (sh1[2][1] * sh5[2][6] - sh1[0][1] * sh5[8][6]);
        sh6[4][8] = sh1[1][1] * sh5[3][7] + kSqrt07_16 * (sh1[2][1] * sh5[4][7] + sh1[0][1] * sh5[6][7]) - kSqrt03_32 * (sh1[2][1] * sh5[2][7] - sh1[0][1] * sh5[8][7]);
        sh6[4][9] = sqrtf(32.0 / 27.0) * sh1[1][1] * sh5[3][8] + sqrtf(14.0 / 27.0) * (sh1[2][1] * sh5[4][8] + sh1[0][1] * sh5[6][8]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh5[2][8] - sh1[0][1] * sh5[8][8]);
        sh6[4][10] = kSqrt08_05 * sh1[1][1] * sh5[3][9] + kSqrt07_10 * (sh1[2][1] * sh5[4][9] + sh1[0][1] * sh5[6][9]) - sqrtf(3.0 / 20.0) * (sh1[2][1] * sh5[2][9] - sh1[0][1] * sh5[8][9]);
        sh6[4][11] = sqrtf(32.0 / 11.0) * sh1[1][1] * sh5[3][10] + sqrtf(14.0 / 11.0) * (sh1[2][1] * sh5[4][10] + sh1[0][1] * sh5[6][10]) - sqrtf(3.0 / 11.0) * (sh1[2][1] * sh5[2][10] - sh1[0][1] * sh5[8][10]);
        sh6[4][12] = sqrtf(8.0 / 33.0) * (sh1[1][2] * sh5[3][10] - sh1[1][0] * sh5[3][0]) + sqrtf(7.0 / 66.0) * ((sh1[2][2] * sh5[4][10] - sh1[2][0] * sh5[4][0]) + (sh1[0][2] * sh5[6][10] - sh1[0][0] * sh5[6][0])) - sqrtf(1.0 / 44.0) * ((sh1[2][2] * sh5[2][10] - sh1[2][0] * sh5[2][0]) - (sh1[0][2] * sh5[8][10] - sh1[0][0] * sh5[8][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[4]);

        sh6[5][0] = sqrtf(35.0 / 132.0) * (sh1[1][2] * sh5[4][0] + sh1[1][0] * sh5[4][10]) + sqrtf(7.0 / 44.0) * (sh1[0][2] * sh5[5][0] + sh1[0][0] * sh5[5][10]) - sqrtf(5.0 / 132.0) * ((sh1[2][2] * sh5[3][0] + sh1[2][0] * sh5[3][10]) - (sh1[0][2] * sh5[7][0] + sh1[0][0] * sh5[7][10]));
        sh6[5][1] = sqrtf(35.0 / 11.0) * sh1[1][1] * sh5[4][0] + sqrtf(21.0 / 11.0) * sh1[0][1] * sh5[5][0] - sqrtf(5.0 / 11.0) * (sh1[2][1] * sh5[3][0] - sh1[0][1] * sh5[7][0]);
        sh6[5][2] = sqrtf(7.0 / 4.0) * sh1[1][1] * sh5[4][1] + sqrtf(21.0 / 20.0) * sh1[0][1] * sh5[5][1] - kSqrt01_04 * (sh1[2][1] * sh5[3][1] - sh1[0][1] * sh5[7][1]);
        sh6[5][3] = sqrtf(35.0 / 27.0) * sh1[1][1] * sh5[4][2] + sqrtf(7.0 / 9.0) * sh1[0][1] * sh5[5][2] - sqrtf(5.0 / 27.0) * (sh1[2][1] * sh5[3][2] - sh1[0][1] * sh5[7][2]);
        sh6[5][4] = sqrtf(35.0 / 32.0) * sh1[1][1] * sh5[4][3] + kSqrt21_32 * sh1[0][1] * sh5[5][3] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh5[3][3] - sh1[0][1] * sh5[7][3]);
        sh6[5][5] = sh1[1][1] * sh5[4][4] + kSqrt03_05 * sh1[0][1] * sh5[5][4] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh5[3][4] - sh1[0][1] * sh5[7][4]);
        sh6[5][6] = sqrtf(35.0 / 36.0) * sh1[1][1] * sh5[4][5] + kSqrt07_12 * sh1[0][1] * sh5[5][5] - sqrtf(5.0 / 36.0) * (sh1[2][1] * sh5[3][5] - sh1[0][1] * sh5[7][5]);
        sh6[5][7] = sh1[1][1] * sh5[4][6] + kSqrt03_05 * sh1[0][1] * sh5[5][6] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh5[3][6] - sh1[0][1] * sh5[7][6]);
        sh6[5][8] = sqrtf(35.0 / 32.0) * sh1[1][1] * sh5[4][7] + kSqrt21_32 * sh1[0][1] * sh5[5][7] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh5[3][7] - sh1[0][1] * sh5[7][7]);
        sh6[5][9] = sqrtf(35.0 / 27.0) * sh1[1][1] * sh5[4][8] + sqrtf(7.0 / 9.0) * sh1[0][1] * sh5[5][8] - sqrtf(5.0 / 27.0) * (sh1[2][1] * sh5[3][8] - sh1[0][1] * sh5[7][8]);
        sh6[5][10] = sqrtf(7.0 / 4.0) * sh1[1][1] * sh5[4][9] + sqrtf(21.0 / 20.0) * sh1[0][1] * sh5[5][9] - kSqrt01_04 * (sh1[2][1] * sh5[3][9] - sh1[0][1] * sh5[7][9]);
        sh6[5][11] = sqrtf(35.0 / 11.0) * sh1[1][1] * sh5[4][10] + sqrtf(21.0 / 11.0) * sh1[0][1] * sh5[5][10] - sqrtf(5.0 / 11.0) * (sh1[2][1] * sh5[3][10] - sh1[0][1] * sh5[7][10]);
        sh6[5][12] = sqrtf(35.0 / 132.0) * (sh1[1][2] * sh5[4][10] - sh1[1][0] * sh5[4][0]) + sqrtf(7.0 / 44.0) * (sh1[0][2] * sh5[5][10] - sh1[0][0] * sh5[5][0]) - sqrtf(5.0 / 132.0) * ((sh1[2][2] * sh5[3][10] - sh1[2][0] * sh5[3][0]) - (sh1[0][2] * sh5[7][10] - sh1[0][0] * sh5[7][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[5]);

        sh6[6][0] = sqrtf(3.0 / 11.0) * (sh1[1][2] * sh5[5][0] + sh1[1][0] * sh5[5][10]) - sqrtf(5.0 / 44.0) * ((sh1[2][2] * sh5[6][0] + sh1[2][0] * sh5[6][10]) + (sh1[0][2] * sh5[4][0] + sh1[0][0] * sh5[4][10]));
        sh6[6][1] = sqrtf(36.0 / 11.0) * sh1[1][1] * sh5[5][0] - sqrtf(15.0 / 11.0) * (sh1[2][1] * sh5[6][0] + sh1[0][1] * sh5[4][0]);
        sh6[6][2] = kSqrt09_05 * sh1[1][1] * sh5[5][1] - kSqrt03_04 * (sh1[2][1] * sh5[6][1] + sh1[0][1] * sh5[4][1]);
        sh6[6][3] = kSqrt04_03 * sh1[1][1] * sh5[5][2] - kSqrt05_09 * (sh1[2][1] * sh5[6][2] + sh1[0][1] * sh5[4][2]);
        sh6[6][4] = kSqrt09_08 * sh1[1][1] * sh5[5][3] - kSqrt15_32 * (sh1[2][1] * sh5[6][3] + sh1[0][1] * sh5[4][3]);
        sh6[6][5] = sqrtf(36.0 / 35.0) * sh1[1][1] * sh5[5][4] - sqrtf(3.0 / 7.0) * (sh1[2][1] * sh5[6][4] + sh1[0][1] * sh5[4][4]);
        sh6[6][6] = sh1[1][1] * sh5[5][5] - sqrtf(5.0 / 12.0) * (sh1[2][1] * sh5[6][5] + sh1[0][1] * sh5[4][5]);
        sh6[6][7] = sqrtf(36.0 / 35.0) * sh1[1][1] * sh5[5][6] - sqrtf(3.0 / 7.0) * (sh1[2][1] * sh5[6][6] + sh1[0][1] * sh5[4][6]);
        sh6[6][8] = kSqrt09_08 * sh1[1][1] * sh5[5][7] - kSqrt15_32 * (sh1[2][1] * sh5[6][7] + sh1[0][1] * sh5[4][7]);
        sh6[6][9] = kSqrt04_03 * sh1[1][1] * sh5[5][8] - kSqrt05_09 * (sh1[2][1] * sh5[6][8] + sh1[0][1] * sh5[4][8]);
        sh6[6][10] = kSqrt09_05 * sh1[1][1] * sh5[5][9] - kSqrt03_04 * (sh1[2][1] * sh5[6][9] + sh1[0][1] * sh5[4][9]);
        sh6[6][11] = sqrtf(36.0 / 11.0) * sh1[1][1] * sh5[5][10] - sqrtf(15.0 / 11.0) * (sh1[2][1] * sh5[6][10] + sh1[0][1] * sh5[4][10]);
        sh6[6][12] = sqrtf(3.0 / 11.0) * (sh1[1][2] * sh5[5][10] - sh1[1][0] * sh5[5][0]) - sqrtf(5.0 / 44.0) * ((sh1[2][2] * sh5[6][10] - sh1[2][0] * sh5[6][0]) + (sh1[0][2] * sh5[4][10] - sh1[0][0] * sh5[4][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[6]);

        sh6[7][0] = sqrtf(35.0 / 132.0) * (sh1[1][2] * sh5[6][0] + sh1[1][0] * sh5[6][10]) + sqrtf(7.0 / 44.0) * (sh1[2][2] * sh5[5][0] + sh1[2][0] * sh5[5][10]) - sqrtf(5.0 / 132.0) * ((sh1[2][2] * sh5[7][0] + sh1[2][0] * sh5[7][10]) + (sh1[0][2] * sh5[3][0] + sh1[0][0] * sh5[3][10]));
        sh6[7][1] = sqrtf(35.0 / 11.0) * sh1[1][1] * sh5[6][0] + sqrtf(21.0 / 11.0) * sh1[2][1] * sh5[5][0] - sqrtf(5.0 / 11.0) * (sh1[2][1] * sh5[7][0] + sh1[0][1] * sh5[3][0]);
        sh6[7][2] = sqrtf(7.0 / 4.0) * sh1[1][1] * sh5[6][1] + sqrtf(21.0 / 20.0) * sh1[2][1] * sh5[5][1] - kSqrt01_04 * (sh1[2][1] * sh5[7][1] + sh1[0][1] * sh5[3][1]);
        sh6[7][3] = sqrtf(35.0 / 27.0) * sh1[1][1] * sh5[6][2] + sqrtf(7.0 / 9.0) * sh1[2][1] * sh5[5][2] - sqrtf(5.0 / 27.0) * (sh1[2][1] * sh5[7][2] + sh1[0][1] * sh5[3][2]);
        sh6[7][4] = sqrtf(35.0 / 32.0) * sh1[1][1] * sh5[6][3] + kSqrt21_32 * sh1[2][1] * sh5[5][3] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh5[7][3] + sh1[0][1] * sh5[3][3]);
        sh6[7][5] = sh1[1][1] * sh5[6][4] + kSqrt03_05 * sh1[2][1] * sh5[5][4] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh5[7][4] + sh1[0][1] * sh5[3][4]);
        sh6[7][6] = sqrtf(35.0 / 36.0) * sh1[1][1] * sh5[6][5] + kSqrt07_12 * sh1[2][1] * sh5[5][5] - sqrtf(5.0 / 36.0) * (sh1[2][1] * sh5[7][5] + sh1[0][1] * sh5[3][5]);
        sh6[7][7] = sh1[1][1] * sh5[6][6] + kSqrt03_05 * sh1[2][1] * sh5[5][6] - sqrtf(1.0 / 7.0) * (sh1[2][1] * sh5[7][6] + sh1[0][1] * sh5[3][6]);
        sh6[7][8] = sqrtf(35.0 / 32.0) * sh1[1][1] * sh5[6][7] + kSqrt21_32 * sh1[2][1] * sh5[5][7] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh5[7][7] + sh1[0][1] * sh5[3][7]);
        sh6[7][9] = sqrtf(35.0 / 27.0) * sh1[1][1] * sh5[6][8] + sqrtf(7.0 / 9.0) * sh1[2][1] * sh5[5][8] - sqrtf(5.0 / 27.0) * (sh1[2][1] * sh5[7][8] + sh1[0][1] * sh5[3][8]);
        sh6[7][10] = sqrtf(7.0 / 4.0) * sh1[1][1] * sh5[6][9] + sqrtf(21.0 / 20.0) * sh1[2][1] * sh5[5][9] - kSqrt01_04 * (sh1[2][1] * sh5[7][9] + sh1[0][1] * sh5[3][9]);
        sh6[7][11] = sqrtf(35.0 / 11.0) * sh1[1][1] * sh5[6][10] + sqrtf(21.0 / 11.0) * sh1[2][1] * sh5[5][10] - sqrtf(5.0 / 11.0) * (sh1[2][1] * sh5[7][10] + sh1[0][1] * sh5[3][10]);
        sh6[7][12] = sqrtf(35.0 / 132.0) * (sh1[1][2] * sh5[6][10] - sh1[1][0] * sh5[6][0]) + sqrtf(7.0 / 44.0) * (sh1[2][2] * sh5[5][10] - sh1[2][0] * sh5[5][0]) - sqrtf(5.0 / 132.0) * ((sh1[2][2] * sh5[7][10] - sh1[2][0] * sh5[7][0]) + (sh1[0][2] * sh5[3][10] - sh1[0][0] * sh5[3][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[7]);

        sh6[8][0] = sqrtf(8.0 / 33.0) * (sh1[1][2] * sh5[7][0] + sh1[1][0] * sh5[7][10]) + sqrtf(7.0 / 66.0) * ((sh1[2][2] * sh5[6][0] + sh1[2][0] * sh5[6][10]) - (sh1[0][2] * sh5[4][0] + sh1[0][0] * sh5[4][10])) - sqrtf(1.0 / 44.0) * ((sh1[2][2] * sh5[8][0] + sh1[2][0] * sh5[8][10]) + (sh1[0][2] * sh5[2][0] + sh1[0][0] * sh5[2][10]));
        sh6[8][1] = sqrtf(32.0 / 11.0) * sh1[1][1] * sh5[7][0] + sqrtf(14.0 / 11.0) * (sh1[2][1] * sh5[6][0] - sh1[0][1] * sh5[4][0]) - sqrtf(3.0 / 11.0) * (sh1[2][1] * sh5[8][0] + sh1[0][1] * sh5[2][0]);
        sh6[8][2] = kSqrt08_05 * sh1[1][1] * sh5[7][1] + kSqrt07_10 * (sh1[2][1] * sh5[6][1] - sh1[0][1] * sh5[4][1]) - sqrtf(3.0 / 20.0) * (sh1[2][1] * sh5[8][1] + sh1[0][1] * sh5[2][1]);
        sh6[8][3] = sqrtf(32.0 / 27.0) * sh1[1][1] * sh5[7][2] + sqrtf(14.0 / 27.0) * (sh1[2][1] * sh5[6][2] - sh1[0][1] * sh5[4][2]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh5[8][2] + sh1[0][1] * sh5[2][2]);
        sh6[8][4] = sh1[1][1] * sh5[7][3] + kSqrt07_16 * (sh1[2][1] * sh5[6][3] - sh1[0][1] * sh5[4][3]) - kSqrt03_32 * (sh1[2][1] * sh5[8][3] + sh1[0][1] * sh5[2][3]);
        sh6[8][5] = sqrtf(32.0 / 35.0) * sh1[1][1] * sh5[7][4] + kSqrt02_05 * (sh1[2][1] * sh5[6][4] - sh1[0][1] * sh5[4][4]) - sqrtf(3.0 / 35.0) * (sh1[2][1] * sh5[8][4] + sh1[0][1] * sh5[2][4]);
        sh6[8][6] = kSqrt08_09 * sh1[1][1] * sh5[7][5] + sqrtf(7.0 / 18.0) * (sh1[2][1] * sh5[6][5] - sh1[0][1] * sh5[4][5]) - kSqrt01_12 * (sh1[2][1] * sh5[8][5] + sh1[0][1] * sh5[2][5]);
        sh6[8][7] = sqrtf(32.0 / 35.0) * sh1[1][1] * sh5[7][6] + kSqrt02_05 * (sh1[2][1] * sh5[6][6] - sh1[0][1] * sh5[4][6]) - sqrtf(3.0 / 35.0) * (sh1[2][1] * sh5[8][6] + sh1[0][1] * sh5[2][6]);
        sh6[8][8] = sh1[1][1] * sh5[7][7] + kSqrt07_16 * (sh1[2][1] * sh5[6][7] - sh1[0][1] * sh5[4][7]) - kSqrt03_32 * (sh1[2][1] * sh5[8][7] + sh1[0][1] * sh5[2][7]);
        sh6[8][9] = sqrtf(32.0 / 27.0) * sh1[1][1] * sh5[7][8] + sqrtf(14.0 / 27.0) * (sh1[2][1] * sh5[6][8] - sh1[0][1] * sh5[4][8]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh5[8][8] + sh1[0][1] * sh5[2][8]);
        sh6[8][10] = kSqrt08_05 * sh1[1][1] * sh5[7][9] + kSqrt07_10 * (sh1[2][1] * sh5[6][9] - sh1[0][1] * sh5[4][9]) - sqrtf(3.0 / 20.0) * (sh1[2][1] * sh5[8][9] + sh1[0][1] * sh5[2][9]);
        sh6[8][11] = sqrtf(32.0 / 11.0) * sh1[1][1] * sh5[7][10] + sqrtf(14.0 / 11.0) * (sh1[2][1] * sh5[6][10] - sh1[0][1] * sh5[4][10]) - sqrtf(3.0 / 11.0) * (sh1[2][1] * sh5[8][10] + sh1[0][1] * sh5[2][10]);
        sh6[8][12] = sqrtf(8.0 / 33.0) * (sh1[1][2] * sh5[7][10] - sh1[1][0] * sh5[7][0]) + sqrtf(7.0 / 66.0) * ((sh1[2][2] * sh5[6][10] - sh1[2][0] * sh5[6][0]) - (sh1[0][2] * sh5[4][10] - sh1[0][0] * sh5[4][0])) - sqrtf(1.0 / 44.0) * ((sh1[2][2] * sh5[8][10] - sh1[2][0] * sh5[8][0]) + (sh1[0][2] * sh5[2][10] - sh1[0][0] * sh5[2][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[8]);

        sh6[9][0] = sqrtf(9.0 / 44.0) * (sh1[1][2] * sh5[8][0] + sh1[1][0] * sh5[8][10]) + sqrtf(3.0 / 22.0) * ((sh1[2][2] * sh5[7][0] + sh1[2][0] * sh5[7][10]) - (sh1[0][2] * sh5[3][0] + sh1[0][0] * sh5[3][10])) - sqrtf(1.0 / 88.0) * ((sh1[2][2] * sh5[9][0] + sh1[2][0] * sh5[9][10]) + (sh1[0][2] * sh5[1][0] + sh1[0][0] * sh5[1][10]));
        sh6[9][1] = sqrtf(27.0 / 11.0) * sh1[1][1] * sh5[8][0] + sqrtf(18.0 / 11.0) * (sh1[2][1] * sh5[7][0] - sh1[0][1] * sh5[3][0]) - sqrtf(3.0 / 22.0) * (sh1[2][1] * sh5[9][0] + sh1[0][1] * sh5[1][0]);
        sh6[9][2] = sqrtf(27.0 / 20.0) * sh1[1][1] * sh5[8][1] + sqrtf(9.0 / 10.0) * (sh1[2][1] * sh5[7][1] - sh1[0][1] * sh5[3][1]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh5[9][1] + sh1[0][1] * sh5[1][1]);
        sh6[9][3] = sh1[1][1] * sh5[8][2] + kSqrt02_03 * (sh1[2][1] * sh5[7][2] - sh1[0][1] * sh5[3][2]) - kSqrt01_18 * (sh1[2][1] * sh5[9][2] + sh1[0][1] * sh5[1][2]);
        sh6[9][4] = sqrtf(27.0 / 32.0) * sh1[1][1] * sh5[8][3] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh5[7][3] - sh1[0][1] * sh5[3][3]) - sqrtf(3.0 / 64.0) * (sh1[2][1] * sh5[9][3] + sh1[0][1] * sh5[1][3]);
        sh6[9][5] = sqrtf(27.0 / 35.0) * sh1[1][1] * sh5[8][4] + sqrtf(18.0 / 35.0) * (sh1[2][1] * sh5[7][4] - sh1[0][1] * sh5[3][4]) - sqrtf(3.0 / 70.0) * (sh1[2][1] * sh5[9][4] + sh1[0][1] * sh5[1][4]);
        sh6[9][6] = kSqrt03_04 * sh1[1][1] * sh5[8][5] + kSqrt01_02 * (sh1[2][1] * sh5[7][5] - sh1[0][1] * sh5[3][5]) - kSqrt01_24 * (sh1[2][1] * sh5[9][5] + sh1[0][1] * sh5[1][5]);
        sh6[9][7] = sqrtf(27.0 / 35.0) * sh1[1][1] * sh5[8][6] + sqrtf(18.0 / 35.0) * (sh1[2][1] * sh5[7][6] - sh1[0][1] * sh5[3][6]) - sqrtf(3.0 / 70.0) * (sh1[2][1] * sh5[9][6] + sh1[0][1] * sh5[1][6]);
        sh6[9][8] = sqrtf(27.0 / 32.0) * sh1[1][1] * sh5[8][7] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh5[7][7] - sh1[0][1] * sh5[3][7]) - sqrtf(3.0 / 64.0) * (sh1[2][1] * sh5[9][7] + sh1[0][1] * sh5[1][7]);
        sh6[9][9] = sh1[1][1] * sh5[8][8] + kSqrt02_03 * (sh1[2][1] * sh5[7][8] - sh1[0][1] * sh5[3][8]) - kSqrt01_18 * (sh1[2][1] * sh5[9][8] + sh1[0][1] * sh5[1][8]);
        sh6[9][10] = sqrtf(27.0 / 20.0) * sh1[1][1] * sh5[8][9] + sqrtf(9.0 / 10.0) * (sh1[2][1] * sh5[7][9] - sh1[0][1] * sh5[3][9]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh5[9][9] + sh1[0][1] * sh5[1][9]);
        sh6[9][11] = sqrtf(27.0 / 11.0) * sh1[1][1] * sh5[8][10] + sqrtf(18.0 / 11.0) * (sh1[2][1] * sh5[7][10] - sh1[0][1] * sh5[3][10]) - sqrtf(3.0 / 22.0) * (sh1[2][1] * sh5[9][10] + sh1[0][1] * sh5[1][10]);
        sh6[9][12] = sqrtf(9.0 / 44.0) * (sh1[1][2] * sh5[8][10] - sh1[1][0] * sh5[8][0]) + sqrtf(3.0 / 22.0) * ((sh1[2][2] * sh5[7][10] - sh1[2][0] * sh5[7][0]) - (sh1[0][2] * sh5[3][10] - sh1[0][0] * sh5[3][0])) - sqrtf(1.0 / 88.0) * ((sh1[2][2] * sh5[9][10] - sh1[2][0] * sh5[9][0]) + (sh1[0][2] * sh5[1][10] - sh1[0][0] * sh5[1][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[9]);

        sh6[10][0] = sqrtf(5.0 / 33.0) * (sh1[1][2] * sh5[9][0] + sh1[1][0] * sh5[9][10]) + sqrtf(15.0 / 88.0) * ((sh1[2][2] * sh5[8][0] + sh1[2][0] * sh5[8][10]) - (sh1[0][2] * sh5[2][0] + sh1[0][0] * sh5[2][10])) - sqrtf(1.0 / 264.0) * ((sh1[2][2] * sh5[10][0] + sh1[2][0] * sh5[10][10]) + (sh1[0][2] * sh5[0][0] + sh1[0][0] * sh5[0][10]));
        sh6[10][1] = sqrtf(20.0 / 11.0) * sh1[1][1] * sh5[9][0] + sqrtf(45.0 / 22.0) * (sh1[2][1] * sh5[8][0] - sh1[0][1] * sh5[2][0]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh5[10][0] + sh1[0][1] * sh5[0][0]);
        sh6[10][2] = sh1[1][1] * sh5[9][1] + kSqrt09_08 * (sh1[2][1] * sh5[8][1] - sh1[0][1] * sh5[2][1]) - sqrtf(1.0 / 40.0) * (sh1[2][1] * sh5[10][1] + sh1[0][1] * sh5[0][1]);
        sh6[10][3] = sqrtf(20.0 / 27.0) * sh1[1][1] * sh5[9][2] + kSqrt05_06 * (sh1[2][1] * sh5[8][2] - sh1[0][1] * sh5[2][2]) - sqrtf(1.0 / 54.0) * (sh1[2][1] * sh5[10][2] + sh1[0][1] * sh5[0][2]);
        sh6[10][4] = kSqrt05_08 * sh1[1][1] * sh5[9][3] + sqrtf(45.0 / 64.0) * (sh1[2][1] * sh5[8][3] - sh1[0][1] * sh5[2][3]) - sqrtf(1.0 / 64.0) * (sh1[2][1] * sh5[10][3] + sh1[0][1] * sh5[0][3]);
        sh6[10][5] = sqrtf(4.0 / 7.0) * sh1[1][1] * sh5[9][4] + sqrtf(9.0 / 14.0) * (sh1[2][1] * sh5[8][4] - sh1[0][1] * sh5[2][4]) - sqrtf(1.0 / 70.0) * (sh1[2][1] * sh5[10][4] + sh1[0][1] * sh5[0][4]);
        sh6[10][6] = kSqrt05_09 * sh1[1][1] * sh5[9][5] + kSqrt05_08 * (sh1[2][1] * sh5[8][5] - sh1[0][1] * sh5[2][5]) - sqrtf(1.0 / 72.0) * (sh1[2][1] * sh5[10][5] + sh1[0][1] * sh5[0][5]);
        sh6[10][7] = sqrtf(4.0 / 7.0) * sh1[1][1] * sh5[9][6] + sqrtf(9.0 / 14.0) * (sh1[2][1] * sh5[8][6] - sh1[0][1] * sh5[2][6]) - sqrtf(1.0 / 70.0) * (sh1[2][1] * sh5[10][6] + sh1[0][1] * sh5[0][6]);
        sh6[10][8] = kSqrt05_08 * sh1[1][1] * sh5[9][7] + sqrtf(45.0 / 64.0) * (sh1[2][1] * sh5[8][7] - sh1[0][1] * sh5[2][7]) - sqrtf(1.0 / 64.0) * (sh1[2][1] * sh5[10][7] + sh1[0][1] * sh5[0][7]);
        sh6[10][9] = sqrtf(20.0 / 27.0) * sh1[1][1] * sh5[9][8] + kSqrt05_06 * (sh1[2][1] * sh5[8][8] - sh1[0][1] * sh5[2][8]) - sqrtf(1.0 / 54.0) * (sh1[2][1] * sh5[10][8] + sh1[0][1] * sh5[0][8]);
        sh6[10][10] = sh1[1][1] * sh5[9][9] + kSqrt09_08 * (sh1[2][1] * sh5[8][9] - sh1[0][1] * sh5[2][9]) - sqrtf(1.0 / 40.0) * (sh1[2][1] * sh5[10][9] + sh1[0][1] * sh5[0][9]);
        sh6[10][11] = sqrtf(20.0 / 11.0) * sh1[1][1] * sh5[9][10] + sqrtf(45.0 / 22.0) * (sh1[2][1] * sh5[8][10] - sh1[0][1] * sh5[2][10]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh5[10][10] + sh1[0][1] * sh5[0][10]);
        sh6[10][12] = sqrtf(5.0 / 33.0) * (sh1[1][2] * sh5[9][10] - sh1[1][0] * sh5[9][0]) + sqrtf(15.0 / 88.0) * ((sh1[2][2] * sh5[8][10] - sh1[2][0] * sh5[8][0]) - (sh1[0][2] * sh5[2][10] - sh1[0][0] * sh5[2][0])) - sqrtf(1.0 / 264.0) * ((sh1[2][2] * sh5[10][10] - sh1[2][0] * sh5[10][0]) + (sh1[0][2] * sh5[0][10] - sh1[0][0] * sh5[0][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[10]);

        sh6[11][0] = kSqrt01_12 * (sh1[1][2] * sh5[10][0] + sh1[1][0] * sh5[10][10]) + sqrtf(5.0 / 24.0) * ((sh1[2][2] * sh5[9][0] + sh1[2][0] * sh5[9][10]) - (sh1[0][2] * sh5[1][0] + sh1[0][0] * sh5[1][10]));
        sh6[11][1] = sh1[1][1] * sh5[10][0] + sqrtf(5.0 / 2.0) * (sh1[2][1] * sh5[9][0] - sh1[0][1] * sh5[1][0]);
        sh6[11][2] = sqrtf(11.0 / 20.0) * sh1[1][1] * sh5[10][1] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh5[9][1] - sh1[0][1] * sh5[1][1]);
        sh6[11][3] = sqrtf(11.0 / 27.0) * sh1[1][1] * sh5[10][2] + sqrtf(55.0 / 54.0) * (sh1[2][1] * sh5[9][2] - sh1[0][1] * sh5[1][2]);
        sh6[11][4] = sqrtf(11.0 / 32.0) * sh1[1][1] * sh5[10][3] + sqrtf(55.0 / 64.0) * (sh1[2][1] * sh5[9][3] - sh1[0][1] * sh5[1][3]);
        sh6[11][5] = sqrtf(11.0 / 35.0) * sh1[1][1] * sh5[10][4] + sqrtf(11.0 / 14.0) * (sh1[2][1] * sh5[9][4] - sh1[0][1] * sh5[1][4]);
        sh6[11][6] = sqrtf(11.0 / 36.0) * sh1[1][1] * sh5[10][5] + sqrtf(55.0 / 72.0) * (sh1[2][1] * sh5[9][5] - sh1[0][1] * sh5[1][5]);
        sh6[11][7] = sqrtf(11.0 / 35.0) * sh1[1][1] * sh5[10][6] + sqrtf(11.0 / 14.0) * (sh1[2][1] * sh5[9][6] - sh1[0][1] * sh5[1][6]);
        sh6[11][8] = sqrtf(11.0 / 32.0) * sh1[1][1] * sh5[10][7] + sqrtf(55.0 / 64.0) * (sh1[2][1] * sh5[9][7] - sh1[0][1] * sh5[1][7]);
        sh6[11][9] = sqrtf(11.0 / 27.0) * sh1[1][1] * sh5[10][8] + sqrtf(55.0 / 54.0) * (sh1[2][1] * sh5[9][8] - sh1[0][1] * sh5[1][8]);
        sh6[11][10] = sqrtf(11.0 / 20.0) * sh1[1][1] * sh5[10][9] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh5[9][9] - sh1[0][1] * sh5[1][9]);
        sh6[11][11] = sh1[1][1] * sh5[10][10] + sqrtf(5.0 / 2.0) * (sh1[2][1] * sh5[9][10] - sh1[0][1] * sh5[1][10]);
        sh6[11][12] = kSqrt01_12 * (sh1[1][2] * sh5[10][10] - sh1[1][0] * sh5[10][0]) + sqrtf(5.0 / 24.0) * ((sh1[2][2] * sh5[9][10] - sh1[2][0] * sh5[9][0]) - (sh1[0][2] * sh5[1][10] - sh1[0][0] * sh5[1][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[11]);

        sh6[12][0] = kSqrt01_04 * ((sh1[2][2] * sh5[10][0] + sh1[2][0] * sh5[10][10]) - (sh1[0][2] * sh5[0][0] + sh1[0][0] * sh5[0][10]));
        sh6[12][1] = sqrtf(3.0 / 1.0) * (sh1[2][1] * sh5[10][0] - sh1[0][1] * sh5[0][0]);
        sh6[12][2] = sqrtf(33.0 / 20.0) * (sh1[2][1] * sh5[10][1] - sh1[0][1] * sh5[0][1]);
        sh6[12][3] = sqrtf(11.0 / 9.0) * (sh1[2][1] * sh5[10][2] - sh1[0][1] * sh5[0][2]);
        sh6[12][4] = sqrtf(33.0 / 32.0) * (sh1[2][1] * sh5[10][3] - sh1[0][1] * sh5[0][3]);
        sh6[12][5] = sqrtf(33.0 / 35.0) * (sh1[2][1] * sh5[10][4] - sh1[0][1] * sh5[0][4]);
        sh6[12][6] = sqrtf(11.0 / 12.0) * (sh1[2][1] * sh5[10][5] - sh1[0][1] * sh5[0][5]);
        sh6[12][7] = sqrtf(33.0 / 35.0) * (sh1[2][1] * sh5[10][6] - sh1[0][1] * sh5[0][6]);
        sh6[12][8] = sqrtf(33.0 / 32.0) * (sh1[2][1] * sh5[10][7] - sh1[0][1] * sh5[0][7]);
        sh6[12][9] = sqrtf(11.0 / 9.0) * (sh1[2][1] * sh5[10][8] - sh1[0][1] * sh5[0][8]);
        sh6[12][10] = sqrtf(33.0 / 20.0) * (sh1[2][1] * sh5[10][9] - sh1[0][1] * sh5[0][9]);
        sh6[12][11] = sqrtf(3.0 / 1.0) * (sh1[2][1] * sh5[10][10] - sh1[0][1] * sh5[0][10]);
        sh6[12][12] = kSqrt01_04 * ((sh1[2][2] * sh5[10][10] - sh1[2][0] * sh5[10][0]) - (sh1[0][2] * sh5[0][10] - sh1[0][0] * sh5[0][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[12]);

        if (n < 8)
            return;

        coeffsIn += 13;
        float sh7[15][15];

        sh7[0][0] = kSqrt01_04 * ((sh1[2][2] * sh6[0][0] + sh1[2][0] * sh6[0][12]) + (sh1[0][2] * sh6[12][0] + sh1[0][0] * sh6[12][12]));
        sh7[0][1] = sqrtf(7.0 / 2.0) * (sh1[2][1] * sh6[0][0] + sh1[0][1] * sh6[12][0]);
        sh7[0][2] = sqrtf(91.0 / 48.0) * (sh1[2][1] * sh6[0][1] + sh1[0][1] * sh6[12][1]);
        sh7[0][3] = sqrtf(91.0 / 66.0) * (sh1[2][1] * sh6[0][2] + sh1[0][1] * sh6[12][2]);
        sh7[0][4] = sqrtf(91.0 / 80.0) * (sh1[2][1] * sh6[0][3] + sh1[0][1] * sh6[12][3]);
        sh7[0][5] = sqrtf(91.0 / 90.0) * (sh1[2][1] * sh6[0][4] + sh1[0][1] * sh6[12][4]);
        sh7[0][6] = sqrtf(91.0 / 96.0) * (sh1[2][1] * sh6[0][5] + sh1[0][1] * sh6[12][5]);
        sh7[0][7] = sqrtf(13.0 / 14.0) * (sh1[2][1] * sh6[0][6] + sh1[0][1] * sh6[12][6]);
        sh7[0][8] = sqrtf(91.0 / 96.0) * (sh1[2][1] * sh6[0][7] + sh1[0][1] * sh6[12][7]);
        sh7[0][9] = sqrtf(91.0 / 90.0) * (sh1[2][1] * sh6[0][8] + sh1[0][1] * sh6[12][8]);
        sh7[0][10] = sqrtf(91.0 / 80.0) * (sh1[2][1] * sh6[0][9] + sh1[0][1] * sh6[12][9]);
        sh7[0][11] = sqrtf(91.0 / 66.0) * (sh1[2][1] * sh6[0][10] + sh1[0][1] * sh6[12][10]);
        sh7[0][12] = sqrtf(91.0 / 48.0) * (sh1[2][1] * sh6[0][11] + sh1[0][1] * sh6[12][11]);
        sh7[0][13] = sqrtf(7.0 / 2.0) * (sh1[2][1] * sh6[0][12] + sh1[0][1] * sh6[12][12]);
        sh7[0][14] = kSqrt01_04 * ((sh1[2][2] * sh6[0][12] - sh1[2][0] * sh6[0][0]) + (sh1[0][2] * sh6[12][12] - sh1[0][0] * sh6[12][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[0]);

        sh7[1][0] = kSqrt01_14 * (sh1[1][2] * sh6[0][0] + sh1[1][0] * sh6[0][12]) + kSqrt03_14 * ((sh1[2][2] * sh6[1][0] + sh1[2][0] * sh6[1][12]) + (sh1[0][2] * sh6[11][0] + sh1[0][0] * sh6[11][12]));
        sh7[1][1] = sh1[1][1] * sh6[0][0] + sqrtf(3.0 / 1.0) * (sh1[2][1] * sh6[1][0] + sh1[0][1] * sh6[11][0]);
        sh7[1][2] = sqrtf(13.0 / 24.0) * sh1[1][1] * sh6[0][1] + sqrtf(13.0 / 8.0) * (sh1[2][1] * sh6[1][1] + sh1[0][1] * sh6[11][1]);
        sh7[1][3] = sqrtf(13.0 / 33.0) * sh1[1][1] * sh6[0][2] + sqrtf(13.0 / 11.0) * (sh1[2][1] * sh6[1][2] + sh1[0][1] * sh6[11][2]);
        sh7[1][4] = sqrtf(13.0 / 40.0) * sh1[1][1] * sh6[0][3] + sqrtf(39.0 / 40.0) * (sh1[2][1] * sh6[1][3] + sh1[0][1] * sh6[11][3]);
        sh7[1][5] = sqrtf(13.0 / 45.0) * sh1[1][1] * sh6[0][4] + sqrtf(13.0 / 15.0) * (sh1[2][1] * sh6[1][4] + sh1[0][1] * sh6[11][4]);
        sh7[1][6] = sqrtf(13.0 / 48.0) * sh1[1][1] * sh6[0][5] + sqrtf(13.0 / 16.0) * (sh1[2][1] * sh6[1][5] + sh1[0][1] * sh6[11][5]);
        sh7[1][7] = sqrtf(13.0 / 49.0) * sh1[1][1] * sh6[0][6] + sqrtf(39.0 / 49.0) * (sh1[2][1] * sh6[1][6] + sh1[0][1] * sh6[11][6]);
        sh7[1][8] = sqrtf(13.0 / 48.0) * sh1[1][1] * sh6[0][7] + sqrtf(13.0 / 16.0) * (sh1[2][1] * sh6[1][7] + sh1[0][1] * sh6[11][7]);
        sh7[1][9] = sqrtf(13.0 / 45.0) * sh1[1][1] * sh6[0][8] + sqrtf(13.0 / 15.0) * (sh1[2][1] * sh6[1][8] + sh1[0][1] * sh6[11][8]);
        sh7[1][10] = sqrtf(13.0 / 40.0) * sh1[1][1] * sh6[0][9] + sqrtf(39.0 / 40.0) * (sh1[2][1] * sh6[1][9] + sh1[0][1] * sh6[11][9]);
        sh7[1][11] = sqrtf(13.0 / 33.0) * sh1[1][1] * sh6[0][10] + sqrtf(13.0 / 11.0) * (sh1[2][1] * sh6[1][10] + sh1[0][1] * sh6[11][10]);
        sh7[1][12] = sqrtf(13.0 / 24.0) * sh1[1][1] * sh6[0][11] + sqrtf(13.0 / 8.0) * (sh1[2][1] * sh6[1][11] + sh1[0][1] * sh6[11][11]);
        sh7[1][13] = sh1[1][1] * sh6[0][12] + sqrtf(3.0 / 1.0) * (sh1[2][1] * sh6[1][12] + sh1[0][1] * sh6[11][12]);
        sh7[1][14] = kSqrt01_14 * (sh1[1][2] * sh6[0][12] - sh1[1][0] * sh6[0][0]) + kSqrt03_14 * ((sh1[2][2] * sh6[1][12] - sh1[2][0] * sh6[1][0]) + (sh1[0][2] * sh6[11][12] - sh1[0][0] * sh6[11][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[1]);

        sh7[2][0] = sqrtf(12.0 / 91.0) * (sh1[1][2] * sh6[1][0] + sh1[1][0] * sh6[1][12]) + sqrtf(33.0 / 182.0) * ((sh1[2][2] * sh6[2][0] + sh1[2][0] * sh6[2][12]) + (sh1[0][2] * sh6[10][0] + sh1[0][0] * sh6[10][12])) - sqrtf(1.0 / 364.0) * ((sh1[2][2] * sh6[0][0] + sh1[2][0] * sh6[0][12]) - (sh1[0][2] * sh6[12][0] + sh1[0][0] * sh6[12][12]));
        sh7[2][1] = sqrtf(24.0 / 13.0) * sh1[1][1] * sh6[1][0] + sqrtf(33.0 / 13.0) * (sh1[2][1] * sh6[2][0] + sh1[0][1] * sh6[10][0]) - sqrtf(1.0 / 26.0) * (sh1[2][1] * sh6[0][0] - sh1[0][1] * sh6[12][0]);
        sh7[2][2] = sh1[1][1] * sh6[1][1] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh6[2][1] + sh1[0][1] * sh6[10][1]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh6[0][1] - sh1[0][1] * sh6[12][1]);
        sh7[2][3] = sqrtf(8.0 / 11.0) * sh1[1][1] * sh6[1][2] + (sh1[2][1] * sh6[2][2] + sh1[0][1] * sh6[10][2]) - sqrtf(1.0 / 66.0) * (sh1[2][1] * sh6[0][2] - sh1[0][1] * sh6[12][2]);
        sh7[2][4] = kSqrt03_05 * sh1[1][1] * sh6[1][3] + sqrtf(33.0 / 40.0) * (sh1[2][1] * sh6[2][3] + sh1[0][1] * sh6[10][3]) - sqrtf(1.0 / 80.0) * (sh1[2][1] * sh6[0][3] - sh1[0][1] * sh6[12][3]);
        sh7[2][5] = sqrtf(8.0 / 15.0) * sh1[1][1] * sh6[1][4] + sqrtf(11.0 / 15.0) * (sh1[2][1] * sh6[2][4] + sh1[0][1] * sh6[10][4]) - sqrtf(1.0 / 90.0) * (sh1[2][1] * sh6[0][4] - sh1[0][1] * sh6[12][4]);
        sh7[2][6] = kSqrt01_02 * sh1[1][1] * sh6[1][5] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[2][5] + sh1[0][1] * sh6[10][5]) - sqrtf(1.0 / 96.0) * (sh1[2][1] * sh6[0][5] - sh1[0][1] * sh6[12][5]);
        sh7[2][7] = sqrtf(24.0 / 49.0) * sh1[1][1] * sh6[1][6] + sqrtf(33.0 / 49.0) * (sh1[2][1] * sh6[2][6] + sh1[0][1] * sh6[10][6]) - sqrtf(1.0 / 98.0) * (sh1[2][1] * sh6[0][6] - sh1[0][1] * sh6[12][6]);
        sh7[2][8] = kSqrt01_02 * sh1[1][1] * sh6[1][7] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[2][7] + sh1[0][1] * sh6[10][7]) - sqrtf(1.0 / 96.0) * (sh1[2][1] * sh6[0][7] - sh1[0][1] * sh6[12][7]);
        sh7[2][9] = sqrtf(8.0 / 15.0) * sh1[1][1] * sh6[1][8] + sqrtf(11.0 / 15.0) * (sh1[2][1] * sh6[2][8] + sh1[0][1] * sh6[10][8]) - sqrtf(1.0 / 90.0) * (sh1[2][1] * sh6[0][8] - sh1[0][1] * sh6[12][8]);
        sh7[2][10] = kSqrt03_05 * sh1[1][1] * sh6[1][9] + sqrtf(33.0 / 40.0) * (sh1[2][1] * sh6[2][9] + sh1[0][1] * sh6[10][9]) - sqrtf(1.0 / 80.0) * (sh1[2][1] * sh6[0][9] - sh1[0][1] * sh6[12][9]);
        sh7[2][11] = sqrtf(8.0 / 11.0) * sh1[1][1] * sh6[1][10] + (sh1[2][1] * sh6[2][10] + sh1[0][1] * sh6[10][10]) - sqrtf(1.0 / 66.0) * (sh1[2][1] * sh6[0][10] - sh1[0][1] * sh6[12][10]);
        sh7[2][12] = sh1[1][1] * sh6[1][11] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh6[2][11] + sh1[0][1] * sh6[10][11]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh6[0][11] - sh1[0][1] * sh6[12][11]);
        sh7[2][13] = sqrtf(24.0 / 13.0) * sh1[1][1] * sh6[1][12] + sqrtf(33.0 / 13.0) * (sh1[2][1] * sh6[2][12] + sh1[0][1] * sh6[10][12]) - sqrtf(1.0 / 26.0) * (sh1[2][1] * sh6[0][12] - sh1[0][1] * sh6[12][12]);
        sh7[2][14] = sqrtf(12.0 / 91.0) * (sh1[1][2] * sh6[1][12] - sh1[1][0] * sh6[1][0]) + sqrtf(33.0 / 182.0) * ((sh1[2][2] * sh6[2][12] - sh1[2][0] * sh6[2][0]) + (sh1[0][2] * sh6[10][12] - sh1[0][0] * sh6[10][0])) - sqrtf(1.0 / 364.0) * ((sh1[2][2] * sh6[0][12] - sh1[2][0] * sh6[0][0]) - (sh1[0][2] * sh6[12][12] - sh1[0][0] * sh6[12][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[2]);

        sh7[3][0] = sqrtf(33.0 / 182.0) * (sh1[1][2] * sh6[2][0] + sh1[1][0] * sh6[2][12]) + sqrtf(55.0 / 364.0) * ((sh1[2][2] * sh6[3][0] + sh1[2][0] * sh6[3][12]) + (sh1[0][2] * sh6[9][0] + sh1[0][0] * sh6[9][12])) - sqrtf(3.0 / 364.0) * ((sh1[2][2] * sh6[1][0] + sh1[2][0] * sh6[1][12]) - (sh1[0][2] * sh6[11][0] + sh1[0][0] * sh6[11][12]));
        sh7[3][1] = sqrtf(33.0 / 13.0) * sh1[1][1] * sh6[2][0] + sqrtf(55.0 / 26.0) * (sh1[2][1] * sh6[3][0] + sh1[0][1] * sh6[9][0]) - sqrtf(3.0 / 26.0) * (sh1[2][1] * sh6[1][0] - sh1[0][1] * sh6[11][0]);
        sh7[3][2] = sqrtf(11.0 / 8.0) * sh1[1][1] * sh6[2][1] + sqrtf(55.0 / 48.0) * (sh1[2][1] * sh6[3][1] + sh1[0][1] * sh6[9][1]) - kSqrt01_16 * (sh1[2][1] * sh6[1][1] - sh1[0][1] * sh6[11][1]);
        sh7[3][3] = sh1[1][1] * sh6[2][2] + kSqrt05_06 * (sh1[2][1] * sh6[3][2] + sh1[0][1] * sh6[9][2]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh6[1][2] - sh1[0][1] * sh6[11][2]);
        sh7[3][4] = sqrtf(33.0 / 40.0) * sh1[1][1] * sh6[2][3] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[3][3] + sh1[0][1] * sh6[9][3]) - sqrtf(3.0 / 80.0) * (sh1[2][1] * sh6[1][3] - sh1[0][1] * sh6[11][3]);
        sh7[3][5] = sqrtf(11.0 / 15.0) * sh1[1][1] * sh6[2][4] + sqrtf(11.0 / 18.0) * (sh1[2][1] * sh6[3][4] + sh1[0][1] * sh6[9][4]) - kSqrt01_30 * (sh1[2][1] * sh6[1][4] - sh1[0][1] * sh6[11][4]);
        sh7[3][6] = sqrtf(11.0 / 16.0) * sh1[1][1] * sh6[2][5] + sqrtf(55.0 / 96.0) * (sh1[2][1] * sh6[3][5] + sh1[0][1] * sh6[9][5]) - kSqrt01_32 * (sh1[2][1] * sh6[1][5] - sh1[0][1] * sh6[11][5]);
        sh7[3][7] = sqrtf(33.0 / 49.0) * sh1[1][1] * sh6[2][6] + sqrtf(55.0 / 98.0) * (sh1[2][1] * sh6[3][6] + sh1[0][1] * sh6[9][6]) - sqrtf(3.0 / 98.0) * (sh1[2][1] * sh6[1][6] - sh1[0][1] * sh6[11][6]);
        sh7[3][8] = sqrtf(11.0 / 16.0) * sh1[1][1] * sh6[2][7] + sqrtf(55.0 / 96.0) * (sh1[2][1] * sh6[3][7] + sh1[0][1] * sh6[9][7]) - kSqrt01_32 * (sh1[2][1] * sh6[1][7] - sh1[0][1] * sh6[11][7]);
        sh7[3][9] = sqrtf(11.0 / 15.0) * sh1[1][1] * sh6[2][8] + sqrtf(11.0 / 18.0) * (sh1[2][1] * sh6[3][8] + sh1[0][1] * sh6[9][8]) - kSqrt01_30 * (sh1[2][1] * sh6[1][8] - sh1[0][1] * sh6[11][8]);
        sh7[3][10] = sqrtf(33.0 / 40.0) * sh1[1][1] * sh6[2][9] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[3][9] + sh1[0][1] * sh6[9][9]) - sqrtf(3.0 / 80.0) * (sh1[2][1] * sh6[1][9] - sh1[0][1] * sh6[11][9]);
        sh7[3][11] = sh1[1][1] * sh6[2][10] + kSqrt05_06 * (sh1[2][1] * sh6[3][10] + sh1[0][1] * sh6[9][10]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh6[1][10] - sh1[0][1] * sh6[11][10]);
        sh7[3][12] = sqrtf(11.0 / 8.0) * sh1[1][1] * sh6[2][11] + sqrtf(55.0 / 48.0) * (sh1[2][1] * sh6[3][11] + sh1[0][1] * sh6[9][11]) - kSqrt01_16 * (sh1[2][1] * sh6[1][11] - sh1[0][1] * sh6[11][11]);
        sh7[3][13] = sqrtf(33.0 / 13.0) * sh1[1][1] * sh6[2][12] + sqrtf(55.0 / 26.0) * (sh1[2][1] * sh6[3][12] + sh1[0][1] * sh6[9][12]) - sqrtf(3.0 / 26.0) * (sh1[2][1] * sh6[1][12] - sh1[0][1] * sh6[11][12]);
        sh7[3][14] = sqrtf(33.0 / 182.0) * (sh1[1][2] * sh6[2][12] - sh1[1][0] * sh6[2][0]) + sqrtf(55.0 / 364.0) * ((sh1[2][2] * sh6[3][12] - sh1[2][0] * sh6[3][0]) + (sh1[0][2] * sh6[9][12] - sh1[0][0] * sh6[9][0])) - sqrtf(3.0 / 364.0) * ((sh1[2][2] * sh6[1][12] - sh1[2][0] * sh6[1][0]) - (sh1[0][2] * sh6[11][12] - sh1[0][0] * sh6[11][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[3]);

        sh7[4][0] = sqrtf(20.0 / 91.0) * (sh1[1][2] * sh6[3][0] + sh1[1][0] * sh6[3][12]) + sqrtf(45.0 / 364.0) * ((sh1[2][2] * sh6[4][0] + sh1[2][0] * sh6[4][12]) + (sh1[0][2] * sh6[8][0] + sh1[0][0] * sh6[8][12])) - sqrtf(3.0 / 182.0) * ((sh1[2][2] * sh6[2][0] + sh1[2][0] * sh6[2][12]) - (sh1[0][2] * sh6[10][0] + sh1[0][0] * sh6[10][12]));
        sh7[4][1] = sqrtf(40.0 / 13.0) * sh1[1][1] * sh6[3][0] + sqrtf(45.0 / 26.0) * (sh1[2][1] * sh6[4][0] + sh1[0][1] * sh6[8][0]) - sqrtf(3.0 / 13.0) * (sh1[2][1] * sh6[2][0] - sh1[0][1] * sh6[10][0]);
        sh7[4][2] = sqrtf(5.0 / 3.0) * sh1[1][1] * sh6[3][1] + kSqrt15_16 * (sh1[2][1] * sh6[4][1] + sh1[0][1] * sh6[8][1]) - kSqrt01_08 * (sh1[2][1] * sh6[2][1] - sh1[0][1] * sh6[10][1]);
        sh7[4][3] = sqrtf(40.0 / 33.0) * sh1[1][1] * sh6[3][2] + sqrtf(15.0 / 22.0) * (sh1[2][1] * sh6[4][2] + sh1[0][1] * sh6[8][2]) - sqrtf(1.0 / 11.0) * (sh1[2][1] * sh6[2][2] - sh1[0][1] * sh6[10][2]);
        sh7[4][4] = sh1[1][1] * sh6[3][3] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh6[4][3] + sh1[0][1] * sh6[8][3]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh6[2][3] - sh1[0][1] * sh6[10][3]);
        sh7[4][5] = kSqrt08_09 * sh1[1][1] * sh6[3][4] + kSqrt01_02 * (sh1[2][1] * sh6[4][4] + sh1[0][1] * sh6[8][4]) - sqrtf(1.0 / 15.0) * (sh1[2][1] * sh6[2][4] - sh1[0][1] * sh6[10][4]);
        sh7[4][6] = kSqrt05_06 * sh1[1][1] * sh6[3][5] + kSqrt15_32 * (sh1[2][1] * sh6[4][5] + sh1[0][1] * sh6[8][5]) - kSqrt01_16 * (sh1[2][1] * sh6[2][5] - sh1[0][1] * sh6[10][5]);
        sh7[4][7] = sqrtf(40.0 / 49.0) * sh1[1][1] * sh6[3][6] + sqrtf(45.0 / 98.0) * (sh1[2][1] * sh6[4][6] + sh1[0][1] * sh6[8][6]) - sqrtf(3.0 / 49.0) * (sh1[2][1] * sh6[2][6] - sh1[0][1] * sh6[10][6]);
        sh7[4][8] = kSqrt05_06 * sh1[1][1] * sh6[3][7] + kSqrt15_32 * (sh1[2][1] * sh6[4][7] + sh1[0][1] * sh6[8][7]) - kSqrt01_16 * (sh1[2][1] * sh6[2][7] - sh1[0][1] * sh6[10][7]);
        sh7[4][9] = kSqrt08_09 * sh1[1][1] * sh6[3][8] + kSqrt01_02 * (sh1[2][1] * sh6[4][8] + sh1[0][1] * sh6[8][8]) - sqrtf(1.0 / 15.0) * (sh1[2][1] * sh6[2][8] - sh1[0][1] * sh6[10][8]);
        sh7[4][10] = sh1[1][1] * sh6[3][9] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh6[4][9] + sh1[0][1] * sh6[8][9]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh6[2][9] - sh1[0][1] * sh6[10][9]);
        sh7[4][11] = sqrtf(40.0 / 33.0) * sh1[1][1] * sh6[3][10] + sqrtf(15.0 / 22.0) * (sh1[2][1] * sh6[4][10] + sh1[0][1] * sh6[8][10]) - sqrtf(1.0 / 11.0) * (sh1[2][1] * sh6[2][10] - sh1[0][1] * sh6[10][10]);
        sh7[4][12] = sqrtf(5.0 / 3.0) * sh1[1][1] * sh6[3][11] + kSqrt15_16 * (sh1[2][1] * sh6[4][11] + sh1[0][1] * sh6[8][11]) - kSqrt01_08 * (sh1[2][1] * sh6[2][11] - sh1[0][1] * sh6[10][11]);
        sh7[4][13] = sqrtf(40.0 / 13.0) * sh1[1][1] * sh6[3][12] + sqrtf(45.0 / 26.0) * (sh1[2][1] * sh6[4][12] + sh1[0][1] * sh6[8][12]) - sqrtf(3.0 / 13.0) * (sh1[2][1] * sh6[2][12] - sh1[0][1] * sh6[10][12]);
        sh7[4][14] = sqrtf(20.0 / 91.0) * (sh1[1][2] * sh6[3][12] - sh1[1][0] * sh6[3][0]) + sqrtf(45.0 / 364.0) * ((sh1[2][2] * sh6[4][12] - sh1[2][0] * sh6[4][0]) + (sh1[0][2] * sh6[8][12] - sh1[0][0] * sh6[8][0])) - sqrtf(3.0 / 182.0) * ((sh1[2][2] * sh6[2][12] - sh1[2][0] * sh6[2][0]) - (sh1[0][2] * sh6[10][12] - sh1[0][0] * sh6[10][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[4]);

        sh7[5][0] = sqrtf(45.0 / 182.0) * (sh1[1][2] * sh6[4][0] + sh1[1][0] * sh6[4][12]) + sqrtf(9.0 / 91.0) * ((sh1[2][2] * sh6[5][0] + sh1[2][0] * sh6[5][12]) + (sh1[0][2] * sh6[7][0] + sh1[0][0] * sh6[7][12])) - sqrtf(5.0 / 182.0) * ((sh1[2][2] * sh6[3][0] + sh1[2][0] * sh6[3][12]) - (sh1[0][2] * sh6[9][0] + sh1[0][0] * sh6[9][12]));
        sh7[5][1] = sqrtf(45.0 / 13.0) * sh1[1][1] * sh6[4][0] + sqrtf(18.0 / 13.0) * (sh1[2][1] * sh6[5][0] + sh1[0][1] * sh6[7][0]) - sqrtf(5.0 / 13.0) * (sh1[2][1] * sh6[3][0] - sh1[0][1] * sh6[9][0]);
        sh7[5][2] = sqrtf(15.0 / 8.0) * sh1[1][1] * sh6[4][1] + kSqrt03_04 * (sh1[2][1] * sh6[5][1] + sh1[0][1] * sh6[7][1]) - sqrtf(5.0 / 24.0) * (sh1[2][1] * sh6[3][1] - sh1[0][1] * sh6[9][1]);
        sh7[5][3] = sqrtf(15.0 / 11.0) * sh1[1][1] * sh6[4][2] + sqrtf(6.0 / 11.0) * (sh1[2][1] * sh6[5][2] + sh1[0][1] * sh6[7][2]) - sqrtf(5.0 / 33.0) * (sh1[2][1] * sh6[3][2] - sh1[0][1] * sh6[9][2]);
        sh7[5][4] = kSqrt09_08 * sh1[1][1] * sh6[4][3] + sqrtf(9.0 / 20.0) * (sh1[2][1] * sh6[5][3] + sh1[0][1] * sh6[7][3]) - kSqrt01_08 * (sh1[2][1] * sh6[3][3] - sh1[0][1] * sh6[9][3]);
        sh7[5][5] = sh1[1][1] * sh6[4][4] + kSqrt02_05 * (sh1[2][1] * sh6[5][4] + sh1[0][1] * sh6[7][4]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh6[3][4] - sh1[0][1] * sh6[9][4]);
        sh7[5][6] = kSqrt15_16 * sh1[1][1] * sh6[4][5] + kSqrt03_08 * (sh1[2][1] * sh6[5][5] + sh1[0][1] * sh6[7][5]) - sqrtf(5.0 / 48.0) * (sh1[2][1] * sh6[3][5] - sh1[0][1] * sh6[9][5]);
        sh7[5][7] = sqrtf(45.0 / 49.0) * sh1[1][1] * sh6[4][6] + sqrtf(18.0 / 49.0) * (sh1[2][1] * sh6[5][6] + sh1[0][1] * sh6[7][6]) - sqrtf(5.0 / 49.0) * (sh1[2][1] * sh6[3][6] - sh1[0][1] * sh6[9][6]);
        sh7[5][8] = kSqrt15_16 * sh1[1][1] * sh6[4][7] + kSqrt03_08 * (sh1[2][1] * sh6[5][7] + sh1[0][1] * sh6[7][7]) - sqrtf(5.0 / 48.0) * (sh1[2][1] * sh6[3][7] - sh1[0][1] * sh6[9][7]);
        sh7[5][9] = sh1[1][1] * sh6[4][8] + kSqrt02_05 * (sh1[2][1] * sh6[5][8] + sh1[0][1] * sh6[7][8]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh6[3][8] - sh1[0][1] * sh6[9][8]);
        sh7[5][10] = kSqrt09_08 * sh1[1][1] * sh6[4][9] + sqrtf(9.0 / 20.0) * (sh1[2][1] * sh6[5][9] + sh1[0][1] * sh6[7][9]) - kSqrt01_08 * (sh1[2][1] * sh6[3][9] - sh1[0][1] * sh6[9][9]);
        sh7[5][11] = sqrtf(15.0 / 11.0) * sh1[1][1] * sh6[4][10] + sqrtf(6.0 / 11.0) * (sh1[2][1] * sh6[5][10] + sh1[0][1] * sh6[7][10]) - sqrtf(5.0 / 33.0) * (sh1[2][1] * sh6[3][10] - sh1[0][1] * sh6[9][10]);
        sh7[5][12] = sqrtf(15.0 / 8.0) * sh1[1][1] * sh6[4][11] + kSqrt03_04 * (sh1[2][1] * sh6[5][11] + sh1[0][1] * sh6[7][11]) - sqrtf(5.0 / 24.0) * (sh1[2][1] * sh6[3][11] - sh1[0][1] * sh6[9][11]);
        sh7[5][13] = sqrtf(45.0 / 13.0) * sh1[1][1] * sh6[4][12] + sqrtf(18.0 / 13.0) * (sh1[2][1] * sh6[5][12] + sh1[0][1] * sh6[7][12]) - sqrtf(5.0 / 13.0) * (sh1[2][1] * sh6[3][12] - sh1[0][1] * sh6[9][12]);
        sh7[5][14] = sqrtf(45.0 / 182.0) * (sh1[1][2] * sh6[4][12] - sh1[1][0] * sh6[4][0]) + sqrtf(9.0 / 91.0) * ((sh1[2][2] * sh6[5][12] - sh1[2][0] * sh6[5][0]) + (sh1[0][2] * sh6[7][12] - sh1[0][0] * sh6[7][0])) - sqrtf(5.0 / 182.0) * ((sh1[2][2] * sh6[3][12] - sh1[2][0] * sh6[3][0]) - (sh1[0][2] * sh6[9][12] - sh1[0][0] * sh6[9][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[5]);

        sh7[6][0] = sqrtf(24.0 / 91.0) * (sh1[1][2] * sh6[5][0] + sh1[1][0] * sh6[5][12]) + sqrtf(2.0 / 13.0) * (sh1[0][2] * sh6[6][0] + sh1[0][0] * sh6[6][12]) - sqrtf(15.0 / 364.0) * ((sh1[2][2] * sh6[4][0] + sh1[2][0] * sh6[4][12]) - (sh1[0][2] * sh6[8][0] + sh1[0][0] * sh6[8][12]));
        sh7[6][1] = sqrtf(48.0 / 13.0) * sh1[1][1] * sh6[5][0] + sqrtf(28.0 / 13.0) * sh1[0][1] * sh6[6][0] - sqrtf(15.0 / 26.0) * (sh1[2][1] * sh6[4][0] - sh1[0][1] * sh6[8][0]);
        sh7[6][2] = kSqrt02_01 * sh1[1][1] * sh6[5][1] + kSqrt07_06 * sh1[0][1] * sh6[6][1] - sqrtf(5.0 / 16.0) * (sh1[2][1] * sh6[4][1] - sh1[0][1] * sh6[8][1]);
        sh7[6][3] = sqrtf(16.0 / 11.0) * sh1[1][1] * sh6[5][2] + sqrtf(28.0 / 33.0) * sh1[0][1] * sh6[6][2] - sqrtf(5.0 / 22.0) * (sh1[2][1] * sh6[4][2] - sh1[0][1] * sh6[8][2]);
        sh7[6][4] = kSqrt06_05 * sh1[1][1] * sh6[5][3] + kSqrt07_10 * sh1[0][1] * sh6[6][3] - kSqrt03_16 * (sh1[2][1] * sh6[4][3] - sh1[0][1] * sh6[8][3]);
        sh7[6][5] = kSqrt16_15 * sh1[1][1] * sh6[5][4] + sqrtf(28.0 / 45.0) * sh1[0][1] * sh6[6][4] - kSqrt01_06 * (sh1[2][1] * sh6[4][4] - sh1[0][1] * sh6[8][4]);
        sh7[6][6] = sh1[1][1] * sh6[5][5] + kSqrt07_12 * sh1[0][1] * sh6[6][5] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh6[4][5] - sh1[0][1] * sh6[8][5]);
        sh7[6][7] = sqrtf(48.0 / 49.0) * sh1[1][1] * sh6[5][6] + sqrtf(4.0 / 7.0) * sh1[0][1] * sh6[6][6] - sqrtf(15.0 / 98.0) * (sh1[2][1] * sh6[4][6] - sh1[0][1] * sh6[8][6]);
        sh7[6][8] = sh1[1][1] * sh6[5][7] + kSqrt07_12 * sh1[0][1] * sh6[6][7] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh6[4][7] - sh1[0][1] * sh6[8][7]);
        sh7[6][9] = kSqrt16_15 * sh1[1][1] * sh6[5][8] + sqrtf(28.0 / 45.0) * sh1[0][1] * sh6[6][8] - kSqrt01_06 * (sh1[2][1] * sh6[4][8] - sh1[0][1] * sh6[8][8]);
        sh7[6][10] = kSqrt06_05 * sh1[1][1] * sh6[5][9] + kSqrt07_10 * sh1[0][1] * sh6[6][9] - kSqrt03_16 * (sh1[2][1] * sh6[4][9] - sh1[0][1] * sh6[8][9]);
        sh7[6][11] = sqrtf(16.0 / 11.0) * sh1[1][1] * sh6[5][10] + sqrtf(28.0 / 33.0) * sh1[0][1] * sh6[6][10] - sqrtf(5.0 / 22.0) * (sh1[2][1] * sh6[4][10] - sh1[0][1] * sh6[8][10]);
        sh7[6][12] = kSqrt02_01 * sh1[1][1] * sh6[5][11] + kSqrt07_06 * sh1[0][1] * sh6[6][11] - sqrtf(5.0 / 16.0) * (sh1[2][1] * sh6[4][11] - sh1[0][1] * sh6[8][11]);
        sh7[6][13] = sqrtf(48.0 / 13.0) * sh1[1][1] * sh6[5][12] + sqrtf(28.0 / 13.0) * sh1[0][1] * sh6[6][12] - sqrtf(15.0 / 26.0) * (sh1[2][1] * sh6[4][12] - sh1[0][1] * sh6[8][12]);
        sh7[6][14] = sqrtf(24.0 / 91.0) * (sh1[1][2] * sh6[5][12] - sh1[1][0] * sh6[5][0]) + sqrtf(2.0 / 13.0) * (sh1[0][2] * sh6[6][12] - sh1[0][0] * sh6[6][0]) - sqrtf(15.0 / 364.0) * ((sh1[2][2] * sh6[4][12] - sh1[2][0] * sh6[4][0]) - (sh1[0][2] * sh6[8][12] - sh1[0][0] * sh6[8][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[6]);

        sh7[7][0] = sqrtf(7.0 / 26.0) * (sh1[1][2] * sh6[6][0] + sh1[1][0] * sh6[6][12]) - sqrtf(3.0 / 26.0) * ((sh1[2][2] * sh6[7][0] + sh1[2][0] * sh6[7][12]) + (sh1[0][2] * sh6[5][0] + sh1[0][0] * sh6[5][12]));
        sh7[7][1] = sqrtf(49.0 / 13.0) * sh1[1][1] * sh6[6][0] - sqrtf(21.0 / 13.0) * (sh1[2][1] * sh6[7][0] + sh1[0][1] * sh6[5][0]);
        sh7[7][2] = sqrtf(49.0 / 24.0) * sh1[1][1] * sh6[6][1] - kSqrt07_08 * (sh1[2][1] * sh6[7][1] + sh1[0][1] * sh6[5][1]);
        sh7[7][3] = sqrtf(49.0 / 33.0) * sh1[1][1] * sh6[6][2] - sqrtf(7.0 / 11.0) * (sh1[2][1] * sh6[7][2] + sh1[0][1] * sh6[5][2]);
        sh7[7][4] = sqrtf(49.0 / 40.0) * sh1[1][1] * sh6[6][3] - sqrtf(21.0 / 40.0) * (sh1[2][1] * sh6[7][3] + sh1[0][1] * sh6[5][3]);
        sh7[7][5] = sqrtf(49.0 / 45.0) * sh1[1][1] * sh6[6][4] - kSqrt07_15 * (sh1[2][1] * sh6[7][4] + sh1[0][1] * sh6[5][4]);
        sh7[7][6] = sqrtf(49.0 / 48.0) * sh1[1][1] * sh6[6][5] - kSqrt07_16 * (sh1[2][1] * sh6[7][5] + sh1[0][1] * sh6[5][5]);
        sh7[7][7] = sh1[1][1] * sh6[6][6] - sqrtf(3.0 / 7.0) * (sh1[2][1] * sh6[7][6] + sh1[0][1] * sh6[5][6]);
        sh7[7][8] = sqrtf(49.0 / 48.0) * sh1[1][1] * sh6[6][7] - kSqrt07_16 * (sh1[2][1] * sh6[7][7] + sh1[0][1] * sh6[5][7]);
        sh7[7][9] = sqrtf(49.0 / 45.0) * sh1[1][1] * sh6[6][8] - kSqrt07_15 * (sh1[2][1] * sh6[7][8] + sh1[0][1] * sh6[5][8]);
        sh7[7][10] = sqrtf(49.0 / 40.0) * sh1[1][1] * sh6[6][9] - sqrtf(21.0 / 40.0) * (sh1[2][1] * sh6[7][9] + sh1[0][1] * sh6[5][9]);
        sh7[7][11] = sqrtf(49.0 / 33.0) * sh1[1][1] * sh6[6][10] - sqrtf(7.0 / 11.0) * (sh1[2][1] * sh6[7][10] + sh1[0][1] * sh6[5][10]);
        sh7[7][12] = sqrtf(49.0 / 24.0) * sh1[1][1] * sh6[6][11] - kSqrt07_08 * (sh1[2][1] * sh6[7][11] + sh1[0][1] * sh6[5][11]);
        sh7[7][13] = sqrtf(49.0 / 13.0) * sh1[1][1] * sh6[6][12] - sqrtf(21.0 / 13.0) * (sh1[2][1] * sh6[7][12] + sh1[0][1] * sh6[5][12]);
        sh7[7][14] = sqrtf(7.0 / 26.0) * (sh1[1][2] * sh6[6][12] - sh1[1][0] * sh6[6][0]) - sqrtf(3.0 / 26.0) * ((sh1[2][2] * sh6[7][12] - sh1[2][0] * sh6[7][0]) + (sh1[0][2] * sh6[5][12] - sh1[0][0] * sh6[5][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[7]);

        sh7[8][0] = sqrtf(24.0 / 91.0) * (sh1[1][2] * sh6[7][0] + sh1[1][0] * sh6[7][12]) + sqrtf(2.0 / 13.0) * (sh1[2][2] * sh6[6][0] + sh1[2][0] * sh6[6][12]) - sqrtf(15.0 / 364.0) * ((sh1[2][2] * sh6[8][0] + sh1[2][0] * sh6[8][12]) + (sh1[0][2] * sh6[4][0] + sh1[0][0] * sh6[4][12]));
        sh7[8][1] = sqrtf(48.0 / 13.0) * sh1[1][1] * sh6[7][0] + sqrtf(28.0 / 13.0) * sh1[2][1] * sh6[6][0] - sqrtf(15.0 / 26.0) * (sh1[2][1] * sh6[8][0] + sh1[0][1] * sh6[4][0]);
        sh7[8][2] = kSqrt02_01 * sh1[1][1] * sh6[7][1] + kSqrt07_06 * sh1[2][1] * sh6[6][1] - sqrtf(5.0 / 16.0) * (sh1[2][1] * sh6[8][1] + sh1[0][1] * sh6[4][1]);
        sh7[8][3] = sqrtf(16.0 / 11.0) * sh1[1][1] * sh6[7][2] + sqrtf(28.0 / 33.0) * sh1[2][1] * sh6[6][2] - sqrtf(5.0 / 22.0) * (sh1[2][1] * sh6[8][2] + sh1[0][1] * sh6[4][2]);
        sh7[8][4] = kSqrt06_05 * sh1[1][1] * sh6[7][3] + kSqrt07_10 * sh1[2][1] * sh6[6][3] - kSqrt03_16 * (sh1[2][1] * sh6[8][3] + sh1[0][1] * sh6[4][3]);
        sh7[8][5] = kSqrt16_15 * sh1[1][1] * sh6[7][4] + sqrtf(28.0 / 45.0) * sh1[2][1] * sh6[6][4] - kSqrt01_06 * (sh1[2][1] * sh6[8][4] + sh1[0][1] * sh6[4][4]);
        sh7[8][6] = sh1[1][1] * sh6[7][5] + kSqrt07_12 * sh1[2][1] * sh6[6][5] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh6[8][5] + sh1[0][1] * sh6[4][5]);
        sh7[8][7] = sqrtf(48.0 / 49.0) * sh1[1][1] * sh6[7][6] + sqrtf(4.0 / 7.0) * sh1[2][1] * sh6[6][6] - sqrtf(15.0 / 98.0) * (sh1[2][1] * sh6[8][6] + sh1[0][1] * sh6[4][6]);
        sh7[8][8] = sh1[1][1] * sh6[7][7] + kSqrt07_12 * sh1[2][1] * sh6[6][7] - sqrtf(5.0 / 32.0) * (sh1[2][1] * sh6[8][7] + sh1[0][1] * sh6[4][7]);
        sh7[8][9] = kSqrt16_15 * sh1[1][1] * sh6[7][8] + sqrtf(28.0 / 45.0) * sh1[2][1] * sh6[6][8] - kSqrt01_06 * (sh1[2][1] * sh6[8][8] + sh1[0][1] * sh6[4][8]);
        sh7[8][10] = kSqrt06_05 * sh1[1][1] * sh6[7][9] + kSqrt07_10 * sh1[2][1] * sh6[6][9] - kSqrt03_16 * (sh1[2][1] * sh6[8][9] + sh1[0][1] * sh6[4][9]);
        sh7[8][11] = sqrtf(16.0 / 11.0) * sh1[1][1] * sh6[7][10] + sqrtf(28.0 / 33.0) * sh1[2][1] * sh6[6][10] - sqrtf(5.0 / 22.0) * (sh1[2][1] * sh6[8][10] + sh1[0][1] * sh6[4][10]);
        sh7[8][12] = kSqrt02_01 * sh1[1][1] * sh6[7][11] + kSqrt07_06 * sh1[2][1] * sh6[6][11] - sqrtf(5.0 / 16.0) * (sh1[2][1] * sh6[8][11] + sh1[0][1] * sh6[4][11]);
        sh7[8][13] = sqrtf(48.0 / 13.0) * sh1[1][1] * sh6[7][12] + sqrtf(28.0 / 13.0) * sh1[2][1] * sh6[6][12] - sqrtf(15.0 / 26.0) * (sh1[2][1] * sh6[8][12] + sh1[0][1] * sh6[4][12]);
        sh7[8][14] = sqrtf(24.0 / 91.0) * (sh1[1][2] * sh6[7][12] - sh1[1][0] * sh6[7][0]) + sqrtf(2.0 / 13.0) * (sh1[2][2] * sh6[6][12] - sh1[2][0] * sh6[6][0]) - sqrtf(15.0 / 364.0) * ((sh1[2][2] * sh6[8][12] - sh1[2][0] * sh6[8][0]) + (sh1[0][2] * sh6[4][12] - sh1[0][0] * sh6[4][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[8]);

        sh7[9][0] = sqrtf(45.0 / 182.0) * (sh1[1][2] * sh6[8][0] + sh1[1][0] * sh6[8][12]) + sqrtf(9.0 / 91.0) * ((sh1[2][2] * sh6[7][0] + sh1[2][0] * sh6[7][12]) - (sh1[0][2] * sh6[5][0] + sh1[0][0] * sh6[5][12])) - sqrtf(5.0 / 182.0) * ((sh1[2][2] * sh6[9][0] + sh1[2][0] * sh6[9][12]) + (sh1[0][2] * sh6[3][0] + sh1[0][0] * sh6[3][12]));
        sh7[9][1] = sqrtf(45.0 / 13.0) * sh1[1][1] * sh6[8][0] + sqrtf(18.0 / 13.0) * (sh1[2][1] * sh6[7][0] - sh1[0][1] * sh6[5][0]) - sqrtf(5.0 / 13.0) * (sh1[2][1] * sh6[9][0] + sh1[0][1] * sh6[3][0]);
        sh7[9][2] = sqrtf(15.0 / 8.0) * sh1[1][1] * sh6[8][1] + kSqrt03_04 * (sh1[2][1] * sh6[7][1] - sh1[0][1] * sh6[5][1]) - sqrtf(5.0 / 24.0) * (sh1[2][1] * sh6[9][1] + sh1[0][1] * sh6[3][1]);
        sh7[9][3] = sqrtf(15.0 / 11.0) * sh1[1][1] * sh6[8][2] + sqrtf(6.0 / 11.0) * (sh1[2][1] * sh6[7][2] - sh1[0][1] * sh6[5][2]) - sqrtf(5.0 / 33.0) * (sh1[2][1] * sh6[9][2] + sh1[0][1] * sh6[3][2]);
        sh7[9][4] = kSqrt09_08 * sh1[1][1] * sh6[8][3] + sqrtf(9.0 / 20.0) * (sh1[2][1] * sh6[7][3] - sh1[0][1] * sh6[5][3]) - kSqrt01_08 * (sh1[2][1] * sh6[9][3] + sh1[0][1] * sh6[3][3]);
        sh7[9][5] = sh1[1][1] * sh6[8][4] + kSqrt02_05 * (sh1[2][1] * sh6[7][4] - sh1[0][1] * sh6[5][4]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh6[9][4] + sh1[0][1] * sh6[3][4]);
        sh7[9][6] = kSqrt15_16 * sh1[1][1] * sh6[8][5] + kSqrt03_08 * (sh1[2][1] * sh6[7][5] - sh1[0][1] * sh6[5][5]) - sqrtf(5.0 / 48.0) * (sh1[2][1] * sh6[9][5] + sh1[0][1] * sh6[3][5]);
        sh7[9][7] = sqrtf(45.0 / 49.0) * sh1[1][1] * sh6[8][6] + sqrtf(18.0 / 49.0) * (sh1[2][1] * sh6[7][6] - sh1[0][1] * sh6[5][6]) - sqrtf(5.0 / 49.0) * (sh1[2][1] * sh6[9][6] + sh1[0][1] * sh6[3][6]);
        sh7[9][8] = kSqrt15_16 * sh1[1][1] * sh6[8][7] + kSqrt03_08 * (sh1[2][1] * sh6[7][7] - sh1[0][1] * sh6[5][7]) - sqrtf(5.0 / 48.0) * (sh1[2][1] * sh6[9][7] + sh1[0][1] * sh6[3][7]);
        sh7[9][9] = sh1[1][1] * sh6[8][8] + kSqrt02_05 * (sh1[2][1] * sh6[7][8] - sh1[0][1] * sh6[5][8]) - sqrtf(1.0 / 9.0) * (sh1[2][1] * sh6[9][8] + sh1[0][1] * sh6[3][8]);
        sh7[9][10] = kSqrt09_08 * sh1[1][1] * sh6[8][9] + sqrtf(9.0 / 20.0) * (sh1[2][1] * sh6[7][9] - sh1[0][1] * sh6[5][9]) - kSqrt01_08 * (sh1[2][1] * sh6[9][9] + sh1[0][1] * sh6[3][9]);
        sh7[9][11] = sqrtf(15.0 / 11.0) * sh1[1][1] * sh6[8][10] + sqrtf(6.0 / 11.0) * (sh1[2][1] * sh6[7][10] - sh1[0][1] * sh6[5][10]) - sqrtf(5.0 / 33.0) * (sh1[2][1] * sh6[9][10] + sh1[0][1] * sh6[3][10]);
        sh7[9][12] = sqrtf(15.0 / 8.0) * sh1[1][1] * sh6[8][11] + kSqrt03_04 * (sh1[2][1] * sh6[7][11] - sh1[0][1] * sh6[5][11]) - sqrtf(5.0 / 24.0) * (sh1[2][1] * sh6[9][11] + sh1[0][1] * sh6[3][11]);
        sh7[9][13] = sqrtf(45.0 / 13.0) * sh1[1][1] * sh6[8][12] + sqrtf(18.0 / 13.0) * (sh1[2][1] * sh6[7][12] - sh1[0][1] * sh6[5][12]) - sqrtf(5.0 / 13.0) * (sh1[2][1] * sh6[9][12] + sh1[0][1] * sh6[3][12]);
        sh7[9][14] = sqrtf(45.0 / 182.0) * (sh1[1][2] * sh6[8][12] - sh1[1][0] * sh6[8][0]) + sqrtf(9.0 / 91.0) * ((sh1[2][2] * sh6[7][12] - sh1[2][0] * sh6[7][0]) - (sh1[0][2] * sh6[5][12] - sh1[0][0] * sh6[5][0])) - sqrtf(5.0 / 182.0) * ((sh1[2][2] * sh6[9][12] - sh1[2][0] * sh6[9][0]) + (sh1[0][2] * sh6[3][12] - sh1[0][0] * sh6[3][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[9]);

        sh7[10][0] = sqrtf(20.0 / 91.0) * (sh1[1][2] * sh6[9][0] + sh1[1][0] * sh6[9][12]) + sqrtf(45.0 / 364.0) * ((sh1[2][2] * sh6[8][0] + sh1[2][0] * sh6[8][12]) - (sh1[0][2] * sh6[4][0] + sh1[0][0] * sh6[4][12])) - sqrtf(3.0 / 182.0) * ((sh1[2][2] * sh6[10][0] + sh1[2][0] * sh6[10][12]) + (sh1[0][2] * sh6[2][0] + sh1[0][0] * sh6[2][12]));
        sh7[10][1] = sqrtf(40.0 / 13.0) * sh1[1][1] * sh6[9][0] + sqrtf(45.0 / 26.0) * (sh1[2][1] * sh6[8][0] - sh1[0][1] * sh6[4][0]) - sqrtf(3.0 / 13.0) * (sh1[2][1] * sh6[10][0] + sh1[0][1] * sh6[2][0]);
        sh7[10][2] = sqrtf(5.0 / 3.0) * sh1[1][1] * sh6[9][1] + kSqrt15_16 * (sh1[2][1] * sh6[8][1] - sh1[0][1] * sh6[4][1]) - kSqrt01_08 * (sh1[2][1] * sh6[10][1] + sh1[0][1] * sh6[2][1]);
        sh7[10][3] = sqrtf(40.0 / 33.0) * sh1[1][1] * sh6[9][2] + sqrtf(15.0 / 22.0) * (sh1[2][1] * sh6[8][2] - sh1[0][1] * sh6[4][2]) - sqrtf(1.0 / 11.0) * (sh1[2][1] * sh6[10][2] + sh1[0][1] * sh6[2][2]);
        sh7[10][4] = sh1[1][1] * sh6[9][3] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh6[8][3] - sh1[0][1] * sh6[4][3]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh6[10][3] + sh1[0][1] * sh6[2][3]);
        sh7[10][5] = kSqrt08_09 * sh1[1][1] * sh6[9][4] + kSqrt01_02 * (sh1[2][1] * sh6[8][4] - sh1[0][1] * sh6[4][4]) - sqrtf(1.0 / 15.0) * (sh1[2][1] * sh6[10][4] + sh1[0][1] * sh6[2][4]);
        sh7[10][6] = kSqrt05_06 * sh1[1][1] * sh6[9][5] + kSqrt15_32 * (sh1[2][1] * sh6[8][5] - sh1[0][1] * sh6[4][5]) - kSqrt01_16 * (sh1[2][1] * sh6[10][5] + sh1[0][1] * sh6[2][5]);
        sh7[10][7] = sqrtf(40.0 / 49.0) * sh1[1][1] * sh6[9][6] + sqrtf(45.0 / 98.0) * (sh1[2][1] * sh6[8][6] - sh1[0][1] * sh6[4][6]) - sqrtf(3.0 / 49.0) * (sh1[2][1] * sh6[10][6] + sh1[0][1] * sh6[2][6]);
        sh7[10][8] = kSqrt05_06 * sh1[1][1] * sh6[9][7] + kSqrt15_32 * (sh1[2][1] * sh6[8][7] - sh1[0][1] * sh6[4][7]) - kSqrt01_16 * (sh1[2][1] * sh6[10][7] + sh1[0][1] * sh6[2][7]);
        sh7[10][9] = kSqrt08_09 * sh1[1][1] * sh6[9][8] + kSqrt01_02 * (sh1[2][1] * sh6[8][8] - sh1[0][1] * sh6[4][8]) - sqrtf(1.0 / 15.0) * (sh1[2][1] * sh6[10][8] + sh1[0][1] * sh6[2][8]);
        sh7[10][10] = sh1[1][1] * sh6[9][9] + sqrtf(9.0 / 16.0) * (sh1[2][1] * sh6[8][9] - sh1[0][1] * sh6[4][9]) - sqrtf(3.0 / 40.0) * (sh1[2][1] * sh6[10][9] + sh1[0][1] * sh6[2][9]);
        sh7[10][11] = sqrtf(40.0 / 33.0) * sh1[1][1] * sh6[9][10] + sqrtf(15.0 / 22.0) * (sh1[2][1] * sh6[8][10] - sh1[0][1] * sh6[4][10]) - sqrtf(1.0 / 11.0) * (sh1[2][1] * sh6[10][10] + sh1[0][1] * sh6[2][10]);
        sh7[10][12] = sqrtf(5.0 / 3.0) * sh1[1][1] * sh6[9][11] + kSqrt15_16 * (sh1[2][1] * sh6[8][11] - sh1[0][1] * sh6[4][11]) - kSqrt01_08 * (sh1[2][1] * sh6[10][11] + sh1[0][1] * sh6[2][11]);
        sh7[10][13] = sqrtf(40.0 / 13.0) * sh1[1][1] * sh6[9][12] + sqrtf(45.0 / 26.0) * (sh1[2][1] * sh6[8][12] - sh1[0][1] * sh6[4][12]) - sqrtf(3.0 / 13.0) * (sh1[2][1] * sh6[10][12] + sh1[0][1] * sh6[2][12]);
        sh7[10][14] = sqrtf(20.0 / 91.0) * (sh1[1][2] * sh6[9][12] - sh1[1][0] * sh6[9][0]) + sqrtf(45.0 / 364.0) * ((sh1[2][2] * sh6[8][12] - sh1[2][0] * sh6[8][0]) - (sh1[0][2] * sh6[4][12] - sh1[0][0] * sh6[4][0])) - sqrtf(3.0 / 182.0) * ((sh1[2][2] * sh6[10][12] - sh1[2][0] * sh6[10][0]) + (sh1[0][2] * sh6[2][12] - sh1[0][0] * sh6[2][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[10]);

        sh7[11][0] = sqrtf(33.0 / 182.0) * (sh1[1][2] * sh6[10][0] + sh1[1][0] * sh6[10][12]) + sqrtf(55.0 / 364.0) * ((sh1[2][2] * sh6[9][0] + sh1[2][0] * sh6[9][12]) - (sh1[0][2] * sh6[3][0] + sh1[0][0] * sh6[3][12])) - sqrtf(3.0 / 364.0) * ((sh1[2][2] * sh6[11][0] + sh1[2][0] * sh6[11][12]) + (sh1[0][2] * sh6[1][0] + sh1[0][0] * sh6[1][12]));
        sh7[11][1] = sqrtf(33.0 / 13.0) * sh1[1][1] * sh6[10][0] + sqrtf(55.0 / 26.0) * (sh1[2][1] * sh6[9][0] - sh1[0][1] * sh6[3][0]) - sqrtf(3.0 / 26.0) * (sh1[2][1] * sh6[11][0] + sh1[0][1] * sh6[1][0]);
        sh7[11][2] = sqrtf(11.0 / 8.0) * sh1[1][1] * sh6[10][1] + sqrtf(55.0 / 48.0) * (sh1[2][1] * sh6[9][1] - sh1[0][1] * sh6[3][1]) - kSqrt01_16 * (sh1[2][1] * sh6[11][1] + sh1[0][1] * sh6[1][1]);
        sh7[11][3] = sh1[1][1] * sh6[10][2] + kSqrt05_06 * (sh1[2][1] * sh6[9][2] - sh1[0][1] * sh6[3][2]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh6[11][2] + sh1[0][1] * sh6[1][2]);
        sh7[11][4] = sqrtf(33.0 / 40.0) * sh1[1][1] * sh6[10][3] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[9][3] - sh1[0][1] * sh6[3][3]) - sqrtf(3.0 / 80.0) * (sh1[2][1] * sh6[11][3] + sh1[0][1] * sh6[1][3]);
        sh7[11][5] = sqrtf(11.0 / 15.0) * sh1[1][1] * sh6[10][4] + sqrtf(11.0 / 18.0) * (sh1[2][1] * sh6[9][4] - sh1[0][1] * sh6[3][4]) - kSqrt01_30 * (sh1[2][1] * sh6[11][4] + sh1[0][1] * sh6[1][4]);
        sh7[11][6] = sqrtf(11.0 / 16.0) * sh1[1][1] * sh6[10][5] + sqrtf(55.0 / 96.0) * (sh1[2][1] * sh6[9][5] - sh1[0][1] * sh6[3][5]) - kSqrt01_32 * (sh1[2][1] * sh6[11][5] + sh1[0][1] * sh6[1][5]);
        sh7[11][7] = sqrtf(33.0 / 49.0) * sh1[1][1] * sh6[10][6] + sqrtf(55.0 / 98.0) * (sh1[2][1] * sh6[9][6] - sh1[0][1] * sh6[3][6]) - sqrtf(3.0 / 98.0) * (sh1[2][1] * sh6[11][6] + sh1[0][1] * sh6[1][6]);
        sh7[11][8] = sqrtf(11.0 / 16.0) * sh1[1][1] * sh6[10][7] + sqrtf(55.0 / 96.0) * (sh1[2][1] * sh6[9][7] - sh1[0][1] * sh6[3][7]) - kSqrt01_32 * (sh1[2][1] * sh6[11][7] + sh1[0][1] * sh6[1][7]);
        sh7[11][9] = sqrtf(11.0 / 15.0) * sh1[1][1] * sh6[10][8] + sqrtf(11.0 / 18.0) * (sh1[2][1] * sh6[9][8] - sh1[0][1] * sh6[3][8]) - kSqrt01_30 * (sh1[2][1] * sh6[11][8] + sh1[0][1] * sh6[1][8]);
        sh7[11][10] = sqrtf(33.0 / 40.0) * sh1[1][1] * sh6[10][9] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[9][9] - sh1[0][1] * sh6[3][9]) - sqrtf(3.0 / 80.0) * (sh1[2][1] * sh6[11][9] + sh1[0][1] * sh6[1][9]);
        sh7[11][11] = sh1[1][1] * sh6[10][10] + kSqrt05_06 * (sh1[2][1] * sh6[9][10] - sh1[0][1] * sh6[3][10]) - sqrtf(1.0 / 22.0) * (sh1[2][1] * sh6[11][10] + sh1[0][1] * sh6[1][10]);
        sh7[11][12] = sqrtf(11.0 / 8.0) * sh1[1][1] * sh6[10][11] + sqrtf(55.0 / 48.0) * (sh1[2][1] * sh6[9][11] - sh1[0][1] * sh6[3][11]) - kSqrt01_16 * (sh1[2][1] * sh6[11][11] + sh1[0][1] * sh6[1][11]);
        sh7[11][13] = sqrtf(33.0 / 13.0) * sh1[1][1] * sh6[10][12] + sqrtf(55.0 / 26.0) * (sh1[2][1] * sh6[9][12] - sh1[0][1] * sh6[3][12]) - sqrtf(3.0 / 26.0) * (sh1[2][1] * sh6[11][12] + sh1[0][1] * sh6[1][12]);
        sh7[11][14] = sqrtf(33.0 / 182.0) * (sh1[1][2] * sh6[10][12] - sh1[1][0] * sh6[10][0]) + sqrtf(55.0 / 364.0) * ((sh1[2][2] * sh6[9][12] - sh1[2][0] * sh6[9][0]) - (sh1[0][2] * sh6[3][12] - sh1[0][0] * sh6[3][0])) - sqrtf(3.0 / 364.0) * ((sh1[2][2] * sh6[11][12] - sh1[2][0] * sh6[11][0]) + (sh1[0][2] * sh6[1][12] - sh1[0][0] * sh6[1][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[11]);

        sh7[12][0] = sqrtf(12.0 / 91.0) * (sh1[1][2] * sh6[11][0] + sh1[1][0] * sh6[11][12]) + sqrtf(33.0 / 182.0) * ((sh1[2][2] * sh6[10][0] + sh1[2][0] * sh6[10][12]) - (sh1[0][2] * sh6[2][0] + sh1[0][0] * sh6[2][12])) - sqrtf(1.0 / 364.0) * ((sh1[2][2] * sh6[12][0] + sh1[2][0] * sh6[12][12]) + (sh1[0][2] * sh6[0][0] + sh1[0][0] * sh6[0][12]));
        sh7[12][1] = sqrtf(24.0 / 13.0) * sh1[1][1] * sh6[11][0] + sqrtf(33.0 / 13.0) * (sh1[2][1] * sh6[10][0] - sh1[0][1] * sh6[2][0]) - sqrtf(1.0 / 26.0) * (sh1[2][1] * sh6[12][0] + sh1[0][1] * sh6[0][0]);
        sh7[12][2] = sh1[1][1] * sh6[11][1] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh6[10][1] - sh1[0][1] * sh6[2][1]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh6[12][1] + sh1[0][1] * sh6[0][1]);
        sh7[12][3] = sqrtf(8.0 / 11.0) * sh1[1][1] * sh6[11][2] + (sh1[2][1] * sh6[10][2] - sh1[0][1] * sh6[2][2]) - sqrtf(1.0 / 66.0) * (sh1[2][1] * sh6[12][2] + sh1[0][1] * sh6[0][2]);
        sh7[12][4] = kSqrt03_05 * sh1[1][1] * sh6[11][3] + sqrtf(33.0 / 40.0) * (sh1[2][1] * sh6[10][3] - sh1[0][1] * sh6[2][3]) - sqrtf(1.0 / 80.0) * (sh1[2][1] * sh6[12][3] + sh1[0][1] * sh6[0][3]);
        sh7[12][5] = sqrtf(8.0 / 15.0) * sh1[1][1] * sh6[11][4] + sqrtf(11.0 / 15.0) * (sh1[2][1] * sh6[10][4] - sh1[0][1] * sh6[2][4]) - sqrtf(1.0 / 90.0) * (sh1[2][1] * sh6[12][4] + sh1[0][1] * sh6[0][4]);
        sh7[12][6] = kSqrt01_02 * sh1[1][1] * sh6[11][5] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[10][5] - sh1[0][1] * sh6[2][5]) - sqrtf(1.0 / 96.0) * (sh1[2][1] * sh6[12][5] + sh1[0][1] * sh6[0][5]);
        sh7[12][7] = sqrtf(24.0 / 49.0) * sh1[1][1] * sh6[11][6] + sqrtf(33.0 / 49.0) * (sh1[2][1] * sh6[10][6] - sh1[0][1] * sh6[2][6]) - sqrtf(1.0 / 98.0) * (sh1[2][1] * sh6[12][6] + sh1[0][1] * sh6[0][6]);
        sh7[12][8] = kSqrt01_02 * sh1[1][1] * sh6[11][7] + sqrtf(11.0 / 16.0) * (sh1[2][1] * sh6[10][7] - sh1[0][1] * sh6[2][7]) - sqrtf(1.0 / 96.0) * (sh1[2][1] * sh6[12][7] + sh1[0][1] * sh6[0][7]);
        sh7[12][9] = sqrtf(8.0 / 15.0) * sh1[1][1] * sh6[11][8] + sqrtf(11.0 / 15.0) * (sh1[2][1] * sh6[10][8] - sh1[0][1] * sh6[2][8]) - sqrtf(1.0 / 90.0) * (sh1[2][1] * sh6[12][8] + sh1[0][1] * sh6[0][8]);
        sh7[12][10] = kSqrt03_05 * sh1[1][1] * sh6[11][9] + sqrtf(33.0 / 40.0) * (sh1[2][1] * sh6[10][9] - sh1[0][1] * sh6[2][9]) - sqrtf(1.0 / 80.0) * (sh1[2][1] * sh6[12][9] + sh1[0][1] * sh6[0][9]);
        sh7[12][11] = sqrtf(8.0 / 11.0) * sh1[1][1] * sh6[11][10] + (sh1[2][1] * sh6[10][10] - sh1[0][1] * sh6[2][10]) - sqrtf(1.0 / 66.0) * (sh1[2][1] * sh6[12][10] + sh1[0][1] * sh6[0][10]);
        sh7[12][12] = sh1[1][1] * sh6[11][11] + sqrtf(11.0 / 8.0) * (sh1[2][1] * sh6[10][11] - sh1[0][1] * sh6[2][11]) - sqrtf(1.0 / 48.0) * (sh1[2][1] * sh6[12][11] + sh1[0][1] * sh6[0][11]);
        sh7[12][13] = sqrtf(24.0 / 13.0) * sh1[1][1] * sh6[11][12] + sqrtf(33.0 / 13.0) * (sh1[2][1] * sh6[10][12] - sh1[0][1] * sh6[2][12]) - sqrtf(1.0 / 26.0) * (sh1[2][1] * sh6[12][12] + sh1[0][1] * sh6[0][12]);
        sh7[12][14] = sqrtf(12.0 / 91.0) * (sh1[1][2] * sh6[11][12] - sh1[1][0] * sh6[11][0]) + sqrtf(33.0 / 182.0) * ((sh1[2][2] * sh6[10][12] - sh1[2][0] * sh6[10][0]) - (sh1[0][2] * sh6[2][12] - sh1[0][0] * sh6[2][0])) - sqrtf(1.0 / 364.0) * ((sh1[2][2] * sh6[12][12] - sh1[2][0] * sh6[12][0]) + (sh1[0][2] * sh6[0][12] - sh1[0][0] * sh6[0][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[12]);

        sh7[13][0] = kSqrt01_14 * (sh1[1][2] * sh6[12][0] + sh1[1][0] * sh6[12][12]) + kSqrt03_14 * ((sh1[2][2] * sh6[11][0] + sh1[2][0] * sh6[11][12]) - (sh1[0][2] * sh6[1][0] + sh1[0][0] * sh6[1][12]));
        sh7[13][1] = sh1[1][1] * sh6[12][0] + sqrtf(3.0 / 1.0) * (sh1[2][1] * sh6[11][0] - sh1[0][1] * sh6[1][0]);
        sh7[13][2] = sqrtf(13.0 / 24.0) * sh1[1][1] * sh6[12][1] + sqrtf(13.0 / 8.0) * (sh1[2][1] * sh6[11][1] - sh1[0][1] * sh6[1][1]);
        sh7[13][3] = sqrtf(13.0 / 33.0) * sh1[1][1] * sh6[12][2] + sqrtf(13.0 / 11.0) * (sh1[2][1] * sh6[11][2] - sh1[0][1] * sh6[1][2]);
        sh7[13][4] = sqrtf(13.0 / 40.0) * sh1[1][1] * sh6[12][3] + sqrtf(39.0 / 40.0) * (sh1[2][1] * sh6[11][3] - sh1[0][1] * sh6[1][3]);
        sh7[13][5] = sqrtf(13.0 / 45.0) * sh1[1][1] * sh6[12][4] + sqrtf(13.0 / 15.0) * (sh1[2][1] * sh6[11][4] - sh1[0][1] * sh6[1][4]);
        sh7[13][6] = sqrtf(13.0 / 48.0) * sh1[1][1] * sh6[12][5] + sqrtf(13.0 / 16.0) * (sh1[2][1] * sh6[11][5] - sh1[0][1] * sh6[1][5]);
        sh7[13][7] = sqrtf(13.0 / 49.0) * sh1[1][1] * sh6[12][6] + sqrtf(39.0 / 49.0) * (sh1[2][1] * sh6[11][6] - sh1[0][1] * sh6[1][6]);
        sh7[13][8] = sqrtf(13.0 / 48.0) * sh1[1][1] * sh6[12][7] + sqrtf(13.0 / 16.0) * (sh1[2][1] * sh6[11][7] - sh1[0][1] * sh6[1][7]);
        sh7[13][9] = sqrtf(13.0 / 45.0) * sh1[1][1] * sh6[12][8] + sqrtf(13.0 / 15.0) * (sh1[2][1] * sh6[11][8] - sh1[0][1] * sh6[1][8]);
        sh7[13][10] = sqrtf(13.0 / 40.0) * sh1[1][1] * sh6[12][9] + sqrtf(39.0 / 40.0) * (sh1[2][1] * sh6[11][9] - sh1[0][1] * sh6[1][9]);
        sh7[13][11] = sqrtf(13.0 / 33.0) * sh1[1][1] * sh6[12][10] + sqrtf(13.0 / 11.0) * (sh1[2][1] * sh6[11][10] - sh1[0][1] * sh6[1][10]);
        sh7[13][12] = sqrtf(13.0 / 24.0) * sh1[1][1] * sh6[12][11] + sqrtf(13.0 / 8.0) * (sh1[2][1] * sh6[11][11] - sh1[0][1] * sh6[1][11]);
        sh7[13][13] = sh1[1][1] * sh6[12][12] + sqrtf(3.0 / 1.0) * (sh1[2][1] * sh6[11][12] - sh1[0][1] * sh6[1][12]);
        sh7[13][14] = kSqrt01_14 * (sh1[1][2] * sh6[12][12] - sh1[1][0] * sh6[12][0]) + kSqrt03_14 * ((sh1[2][2] * sh6[11][12] - sh1[2][0] * sh6[11][0]) - (sh1[0][2] * sh6[1][12] - sh1[0][0] * sh6[1][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[13]);

        sh7[14][0] = kSqrt01_04 * ((sh1[2][2] * sh6[12][0] + sh1[2][0] * sh6[12][12]) - (sh1[0][2] * sh6[0][0] + sh1[0][0] * sh6[0][12]));
        sh7[14][1] = sqrtf(7.0 / 2.0) * (sh1[2][1] * sh6[12][0] - sh1[0][1] * sh6[0][0]);
        sh7[14][2] = sqrtf(91.0 / 48.0) * (sh1[2][1] * sh6[12][1] - sh1[0][1] * sh6[0][1]);
        sh7[14][3] = sqrtf(91.0 / 66.0) * (sh1[2][1] * sh6[12][2] - sh1[0][1] * sh6[0][2]);
        sh7[14][4] = sqrtf(91.0 / 80.0) * (sh1[2][1] * sh6[12][3] - sh1[0][1] * sh6[0][3]);
        sh7[14][5] = sqrtf(91.0 / 90.0) * (sh1[2][1] * sh6[12][4] - sh1[0][1] * sh6[0][4]);
        sh7[14][6] = sqrtf(91.0 / 96.0) * (sh1[2][1] * sh6[12][5] - sh1[0][1] * sh6[0][5]);
        sh7[14][7] = sqrtf(13.0 / 14.0) * (sh1[2][1] * sh6[12][6] - sh1[0][1] * sh6[0][6]);
        sh7[14][8] = sqrtf(91.0 / 96.0) * (sh1[2][1] * sh6[12][7] - sh1[0][1] * sh6[0][7]);
        sh7[14][9] = sqrtf(91.0 / 90.0) * (sh1[2][1] * sh6[12][8] - sh1[0][1] * sh6[0][8]);
        sh7[14][10] = sqrtf(91.0 / 80.0) * (sh1[2][1] * sh6[12][9] - sh1[0][1] * sh6[0][9]);
        sh7[14][11] = sqrtf(91.0 / 66.0) * (sh1[2][1] * sh6[12][10] - sh1[0][1] * sh6[0][10]);
        sh7[14][12] = sqrtf(91.0 / 48.0) * (sh1[2][1] * sh6[12][11] - sh1[0][1] * sh6[0][11]);
        sh7[14][13] = sqrtf(7.0 / 2.0) * (sh1[2][1] * sh6[12][12] - sh1[0][1] * sh6[0][12]);
        sh7[14][14] = kSqrt01_04 * ((sh1[2][2] * sh6[12][12] - sh1[2][0] * sh6[12][0]) - (sh1[0][2] * sh6[0][12] - sh1[0][0] * sh6[0][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[14]);
    }
}

void SHL::RRotateSH(const Mat3f& rot_row, int n, const float* coeffsIn, float* coeffs)
{
    ::RotateSH<float>(rot_row, n, coeffsIn, coeffs);
}

void SHL::RRotateSH(const Mat3f& rot_row, int n, const Vec4f* coeffsIn, Vec4f* coeffs)
{
    ::RotateSH<Vec4f>(rot_row, n, coeffsIn, coeffs);
}

void SHL::CRotateSH(const Mat3f& rot_col, int n, const float* coeffsIn, float* coeffs)
{
    ::RotateSH<float>(trans(rot_col), n, coeffsIn, coeffs);
}

void SHL::CRotateSH(const Mat3f& rot_col, int n, const Vec4f* coeffsIn, Vec4f* coeffs)
{
    ::RotateSH<Vec4f>(trans(rot_col), n, coeffsIn, coeffs);
}

namespace
{
    template<typename T> inline void ApplyNormalizationConstants(int n, T* coeffs)
    {
        VL_ASSERT(n > 0);

        (*coeffs++) *= kSH_Y_00;

        if (n < 2)
            return;

        (*coeffs++) *= kSH_Y_10;
        (*coeffs++) *= kSH_Y_11;
        (*coeffs++) *= kSH_Y_12;

        if (n < 3)
            return;

        (*coeffs++) *= kSH_Y_20;
        (*coeffs++) *= kSH_Y_21;
        (*coeffs++) *= kSH_Y_22;
        (*coeffs++) *= kSH_Y_23;
        (*coeffs++) *= kSH_Y_24;

        if (n < 4)
            return;

        (*coeffs++) *= kSH_Y_30;
        (*coeffs++) *= kSH_Y_31;
        (*coeffs++) *= kSH_Y_32;
        (*coeffs++) *= kSH_Y_33;
        (*coeffs++) *= kSH_Y_34;
        (*coeffs++) *= kSH_Y_35;
        (*coeffs++) *= kSH_Y_36;

        if (n < 5)
            return;

        (*coeffs++) *= kSH_Y_40;
        (*coeffs++) *= kSH_Y_41;
        (*coeffs++) *= kSH_Y_42;
        (*coeffs++) *= kSH_Y_43;
        (*coeffs++) *= kSH_Y_44;
        (*coeffs++) *= kSH_Y_45;
        (*coeffs++) *= kSH_Y_46;
        (*coeffs++) *= kSH_Y_47;
        (*coeffs++) *= kSH_Y_48;

        VL_ASSERT(n < 6);
    }

    template<typename T> void RemoveNormalizationConstants(int n, T* coeffs)
    {
        VL_ASSERT(n > 0);

        (*coeffs++) *= 1.0f / kSH_Y_00;

        if (n < 2)
            return;

        (*coeffs++) *= 1.0f / kSH_Y_10;
        (*coeffs++) *= 1.0f / kSH_Y_11;
        (*coeffs++) *= 1.0f / kSH_Y_12;

        if (n < 3)
            return;

        (*coeffs++) *= 1.0f / kSH_Y_20;
        (*coeffs++) *= 1.0f / kSH_Y_21;
        (*coeffs++) *= 1.0f / kSH_Y_22;
        (*coeffs++) *= 1.0f / kSH_Y_23;
        (*coeffs++) *= 1.0f / kSH_Y_24;

        if (n < 4)
            return;

        (*coeffs++) *= 1.0f / kSH_Y_30;
        (*coeffs++) *= 1.0f / kSH_Y_31;
        (*coeffs++) *= 1.0f / kSH_Y_32;
        (*coeffs++) *= 1.0f / kSH_Y_33;
        (*coeffs++) *= 1.0f / kSH_Y_34;
        (*coeffs++) *= 1.0f / kSH_Y_35;
        (*coeffs++) *= 1.0f / kSH_Y_36;

        if (n < 5)
            return;

        (*coeffs++) *= 1.0f / kSH_Y_40;
        (*coeffs++) *= 1.0f / kSH_Y_41;
        (*coeffs++) *= 1.0f / kSH_Y_42;
        (*coeffs++) *= 1.0f / kSH_Y_43;
        (*coeffs++) *= 1.0f / kSH_Y_44;
        (*coeffs++) *= 1.0f / kSH_Y_45;
        (*coeffs++) *= 1.0f / kSH_Y_46;
        (*coeffs++) *= 1.0f / kSH_Y_47;
        (*coeffs++) *= 1.0f / kSH_Y_48;

        VL_ASSERT(n < 6);
    }

}

void SHL::ApplyNormalizationConstants(int n, float* coeffs)
{
    ::ApplyNormalizationConstants<float>(n, coeffs);
}

void SHL::ApplyNormalizationConstants(int n, Vec4f* coeffs)
{
    ::ApplyNormalizationConstants<Vec4f>(n, coeffs);
}

void SHL::RemoveNormalizationConstants(int n, float* coeffs)
{
    ::RemoveNormalizationConstants<float>(n, coeffs);
}

void SHL::RemoveNormalizationConstants(int n, Vec4f* coeffs)
{
    ::RemoveNormalizationConstants<Vec4f>(n, coeffs);
}

namespace
{
    const float kSH_Ym_00 = 1.0             / kSH_Y_00; // 1
    const float kSH_Ym_10 = 1.0             / kSH_Y_10; // y
    const float kSH_Ym_11 = 1.0             / kSH_Y_11; // z
    const float kSH_Ym_12 = 1.0             / kSH_Y_12; // x
    const float kSH_Ym_20 = 2.0             / kSH_Y_20; // xy
    const float kSH_Ym_21 = 2.0             / kSH_Y_21; // yz
    const float kSH_Ym_22 = 0.5             / kSH_Y_22; // (3z^2 - 1)
    const float kSH_Ym_23 = 2.0             / kSH_Y_23; // xz
    const float kSH_Ym_24 = 1.0             / kSH_Y_24; // (x^2 - y^2)
    const float kSH_Ym_30 = 1.0             / kSH_Y_30; // y (3x^2 - y^2)
    const float kSH_Ym_31 = sqrtf(27.0)     / kSH_Y_31; // z xy
    const float kSH_Ym_32 = 0.726184        / kSH_Y_32; // y (5z^2 - 1)
    const float kSH_Ym_33 = 0.5             / kSH_Y_33; // z (5z^2 - 3)
    const float kSH_Ym_34 = 0.726184        / kSH_Y_34; // x (5z^2 - 1)
    const float kSH_Ym_35 = sqrtf(27.0 / 4) / kSH_Y_35; // z (x^2 - y^2)
    const float kSH_Ym_36 = 1.0             / kSH_Y_36; // x (x^2 - 3y^2)
    const float kSH_Ym_40 = 4.0             / kSH_Y_40; // (x^2 - y^2)  xy
    const float kSH_Ym_41 = 3.079202        / kSH_Y_41; // (3x^2 - y^2) yz
    const float kSH_Ym_42 = 14.0 / 9.0      / kSH_Y_42; // (7z^2 - 1)   xy
    const float kSH_Ym_43 = 0.946946        / kSH_Y_43; // (7z^2 - 3)   yz
    const float kSH_Ym_44 = 0.125           / kSH_Y_44; // (3 - 30z^2 + 35z^4)
    const float kSH_Ym_45 = 0.946946        / kSH_Y_45; // (7z^2 - 3)   xz
    const float kSH_Ym_46 = 7.0 / 9.0       / kSH_Y_46; // (7z^2 - 1)   (x^2 -y^2)
    const float kSH_Ym_47 = 3.079202        / kSH_Y_47; // (x^2 - 3y^2) xz
    const float kSH_Ym_48 = 1.0             / kSH_Y_48; // (x^4 - 6 x^2 y^2 + y^4)
}

void SHL::ApplyMaxScale(int n, Vec4f* coeffs)
{
    VL_ASSERT(n > 0);

    (*coeffs++) *= kSH_Ym_00;

    if (n < 2)
        return;

    (*coeffs++) *= kSH_Ym_10;
    (*coeffs++) *= kSH_Ym_11;
    (*coeffs++) *= kSH_Ym_12;

    if (n < 3)
        return;

    (*coeffs++) *= kSH_Ym_20;
    (*coeffs++) *= kSH_Ym_21;
    (*coeffs++) *= kSH_Ym_22;
    (*coeffs++) *= kSH_Ym_23;
    (*coeffs++) *= kSH_Ym_24;

    if (n < 4)
        return;

    (*coeffs++) *= kSH_Ym_30;
    (*coeffs++) *= kSH_Ym_31;
    (*coeffs++) *= kSH_Ym_32;
    (*coeffs++) *= kSH_Ym_33;
    (*coeffs++) *= kSH_Ym_34;
    (*coeffs++) *= kSH_Ym_35;
    (*coeffs++) *= kSH_Ym_36;

    if (n < 5)
        return;

    (*coeffs++) *= kSH_Ym_40;
    (*coeffs++) *= kSH_Ym_41;
    (*coeffs++) *= kSH_Ym_42;
    (*coeffs++) *= kSH_Ym_43;
    (*coeffs++) *= kSH_Ym_44;
    (*coeffs++) *= kSH_Ym_45;
    (*coeffs++) *= kSH_Ym_46;
    (*coeffs++) *= kSH_Ym_47;
    (*coeffs++) *= kSH_Ym_48;

    VL_ASSERT(n < 6);
}

namespace
{
    template<typename T> inline void ApplyDiffuseBRDF(int n, T* coeffs)
    {
        // constant coeff stays the same: convolving constant with diffuse
        // leaves it untouched.
        if (n < 2)
            return;

        coeffs[1] *= kSH_A1;
        coeffs[2] *= kSH_A1;
        coeffs[3] *= kSH_A1;

        if (n < 3)
            return;

        for (int i = 4; i < 9; i++)
            coeffs[i] *= kSH_A2;

        if (n < 4)
            return;

        for (int i = 9; i < 16; i++)
            coeffs[i] = vl_0;

        if (n < 5)
            return;

        for (int i = 16; i < 25; i++)
            coeffs[i] *= kSH_A4;

        if (n < 6)
            return;

        int nc = n * n;
        for (int i = 25; i < nc; i++)
            coeffs[i] = vl_0;
    }
}

void SHL::ApplyDiffuseBRDF(int n, float* coeffs)
{
    ::ApplyDiffuseBRDF<float>(n, coeffs);
}

void SHL::ApplyDiffuseBRDF(int n, Vec4f* coeffs)
{
    ::ApplyDiffuseBRDF<Vec4f>(n, coeffs);
}

namespace
{
    template<typename T> inline void ApplyGlossyBRDF(float gloss, int n, T coeffs[])
    {
        if (n < 2)
            return;

        if (gloss >= 1.0f)
            return;

        float diff = 1.0f - gloss;

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
}

void SHL::ApplyGlossyBRDF(float gloss, int n, float coeffs[])
{
    ::ApplyGlossyBRDF<float>(gloss, n, coeffs);
}

void SHL::ApplyGlossyBRDF(float gloss, int n, Vec4f coeffs[])
{
    ::ApplyGlossyBRDF<Vec4f>(gloss, n, coeffs);
}


/////////////////////////////////////////////////////////////////////////////
// SH convolution with ZH -- super fast!
//
// Useful to convolve an environment map with any radially-symmetric BRDF. (The vector
// that the BRDF is symmetric about is then used to look up the reflection value in
// the convolved environment, e.g. the normal for diffuse, the reflection vector for
// specular BRDFs.)
//
// (Note: SH convolution with SH is not repesentable as a spherical function, i.e. you can't represent the result with SH.)

void SHL::ConvolveSHWithZH(int numBands, const float zcoeffs[], const float coeffsIn[], float coeffsOut[])
{
    /*
        Example: Calling this with CalcCosPowerSatZH7(1, ...) should be equivalent to the classic band-weighting from ApplyDiffuseBRDF()

        const float kZH_Y_0 = sqrtf( 1 / (   4 * pi));
        const float kZH_Y_1 = sqrtf( 3 / (   4 * pi));
        const float kZH_Y_2 = sqrtf( 5 / (  16 * pi));
        const float kZH_Y_3 = sqrtf( 7 / (  16 * pi));
        const float kZH_Y_4 = sqrtf( 9 / ( 256 * pi));

        zcoeffs[0] =   1 / 2 *  2 * pi * kZH_Y_0   * sqrtf( 4 * pi / 1) = pi
        zcoeffs[1] =   1 / 3 *  2 * pi * kZH_Y_1   * sqrtf( 4 * pi / 3) = 2 / 3 pi
        zcoeffs[2] =   1 / 4 *  2 * pi * kZH_Y_2   * sqrtf( 4 * pi / 5) = 2 / 4 / 2 = 1/4 pi
        zcoeffs[3] =   0
        zcoeffs[4] =  -1/6   *  2 * pi * kZH_Y_4   * sqrtf( 4 * pi / 9) = -2 / 6 / 8 = - 1 / 24  pi

        kSH_A0 = 1.0,
        kSH_A1 = (2.0 / 3.0),
        kSH_A2 = (1.0 / 4.0),
        kSH_A3 = 0.0,
        kSH_A4 = -(1.0 / 24.0);

        Checks, except for pi factor -- reflective of CalcCosPowerSatZH7 being sat(cos(theta)) -- integral of that
        over sphere is pi, and if we want our diffuse BRDF to be normalised, it should thus be CalcCosPowerSatZH7(1) / pi.
    */

    for (int i = 0; i < numBands; i++)
    {
        int n = (2 * i + 1);
        float alpha = sqrtf(kFourPi / n);
        float alphaZ = alpha * zcoeffs[i];

        for (int j = 0; j < n; j++)
            coeffsOut[j] = coeffsIn[j] * alphaZ;

        coeffsIn += n;
        coeffsOut += n;
    }
}

void SHL::ConvolveSHWithZH(int numBands, const float zcoeffs[], const Vec4f coeffsIn[], Vec4f coeffsOut[])
{
    for (int i = 0; i < numBands; i++)
    {
        int n = (2 * i + 1);
        float alpha = sqrtf(kFourPi / n);
        float alphaZ = alpha * zcoeffs[i];

        for (int j = 0; j < n; j++)
            coeffsOut[j] = coeffsIn[j] * alphaZ;

        coeffsIn += n;
        coeffsOut += n;
    }
}

// Windowing
void SHL::ApplySHWindowing(float gamma, int n, float* coeffs)
{
    for (int i = 0; i < n; i++)
    {
        float w = WindowScale(i, gamma);

        for (int j = 0, nj = 2 * i + 1; j < nj; j++)
            *coeffs++ *= w;
    }
}

void SHL::ApplySHWindowing(float gamma, int n, Vec4f* coeffs)
{
    for (int i = 0; i < n; i++)
    {
        float w = WindowScale(i, gamma);

        for (int j = 0, nj = 2 * i + 1; j < nj; j++)
            *coeffs++ *= w;
    }
}


/////////////////////////////////////////////////////////////////////////////
// SH Multiplication. Requires triple-product basis coefficients.
//

void SHL::MultiplySHByZH3(const float a[9], const float b[3], float c[9])
{
    VL_ASSERT(a != b && b != c);

    c[0]  = a[0] * b[0] * kSH_Y_00;
    c[0] += a[2] * b[1] * kSH_Y_00;
    c[0] += a[6] * b[2] * kSH_Y_00;

    c[1]  = a[1] * b[0] * kSH_Y_00;
    c[1] += a[5] * b[1] * kSH_Y_00 * kSqrt03_05;
    c[1] += a[1] * b[2] * kSH_Y_00 * -kSqrt01_05; // kSH_Y_00 * kSqrt01_05
    c[2]  = a[2] * b[0] * kSH_Y_00;
    c[2] += a[0] * b[1] * kSH_Y_00;
    c[2] += a[6] * b[1] * kSH_Y_00 * kSqrt04_05; // kSH_Y_00 * kSqrt04_05
    c[2] += a[2] * b[2] * kSH_Y_00 * kSqrt04_05;
    c[3]  = a[3] * b[0] * kSH_Y_00;
    c[3] += a[7] * b[1] * kSH_Y_00 * kSqrt03_05;
    c[3] += a[3] * b[2] * kSH_Y_00 * -kSqrt01_05;

    c[4]  = a[4] * b[0] * kSH_Y_00;
    c[4] += a[4] * b[2] * -0.180258f;
    c[5]  = a[5] * b[0] * kSH_Y_00;
    c[5] += a[1] * b[1] * kSH_Y_00 * kSqrt03_05;
    c[5] += a[5] * b[2] * 0.0901372f;
    c[6]  = a[6] * b[0] * kSH_Y_00;
    c[6] += a[2] * b[1] * kSH_Y_00 * kSqrt04_05;
    c[6] += a[0] * b[2] * kSH_Y_00;
    c[6] += a[6] * b[2] * 0.180203f;
    c[7]  = a[7] * b[0] * kSH_Y_00;
    c[7] += a[3] * b[1] * kSH_Y_00 * kSqrt03_05;
    c[7] += a[7] * b[2] * 0.0900888f;
    c[8]  = a[8] * b[0] * kSH_Y_00;
    c[8] += a[8] * b[2] * -0.180213f;
}

namespace
{
    struct BasisTriple
    {
        int8_t i;
        int8_t j;
        int8_t k;
        float  s;
    };

    const BasisTriple kBasisTriples2[] =
    {
        { 0, 0, 0,  kSH_Y_00 },
        { 0, 1, 1,  kSH_Y_00 },
        { 0, 2, 2,  kSH_Y_00 },
        { 0, 3, 3,  kSH_Y_00 },

        { -1, -1, -1, 0.0 }
    };

    const BasisTriple kBasisTriples3[] =
    {
        { 0, 0, 0,  kSH_Y_00 },
        { 0, 1, 1,  kSH_Y_00 },
        { 0, 2, 2,  kSH_Y_00 },
        { 0, 3, 3,  kSH_Y_00 },
        { 0, 4, 4,  kSH_Y_00 },
        { 0, 5, 5,  kSH_Y_00 },
        { 0, 6, 6,  kSH_Y_00 },
        { 0, 7, 7,  kSH_Y_00 },
        { 0, 8, 8,  kSH_Y_00 },

        { 1, 1, 6, -0.126147 },
        { 1, 1, 8, -0.218467 },
        { 1, 2, 5,  0.218513 },
        { 1, 3, 4,  0.218540 },

        { 2, 2, 6,  0.252306 },
        { 2, 3, 7,  0.218495 },

        { 3, 3, 6, -0.126184 },
        { 3, 3, 8,  0.218513 },

        { 4, 4, 6, -0.180258 },
        { 4, 5, 7,  0.156089 },

        { 5, 5, 6,  0.090137 },
        { 5, 5, 8, -0.156058 },

        { 6, 6, 6,  0.180203 },
        { 6, 7, 7,  0.090088 },
        { 6, 8, 8, -0.180213 },

        { 7, 7, 8,  0.156056 },

        { -1, -1, -1, 0.0 }
    };

    const BasisTriple kBasisTriples4[] =
    {
        { 0,  0,  0,  kSH_Y_00 },
        { 0,  1,  1,  kSH_Y_00 },
        { 0,  2,  2,  kSH_Y_00 },
        { 0,  3,  3,  kSH_Y_00 },
        { 0,  4,  4,  kSH_Y_00 },
        { 0,  5,  5,  kSH_Y_00 },
        { 0,  6,  6,  kSH_Y_00 },
        { 0,  7,  7,  kSH_Y_00 },
        { 0,  8,  8,  kSH_Y_00 },
        { 0,  9,  9,  kSH_Y_00 },
        { 0, 10, 10,  kSH_Y_00 },
        { 0, 11, 11,  kSH_Y_00 },
        { 0, 12, 12,  kSH_Y_00 },
        { 0, 13, 13,  kSH_Y_00 },
        { 0, 14, 14,  kSH_Y_00 },
        { 0, 15, 15,  kSH_Y_00 },

        { 1,  1,  6, -0.126147 },
        { 1,  1,  8, -0.218467 },
        { 1,  2,  5,  0.218513 },
        { 1,  3,  4,  0.218540 },
        { 1,  4, 13, -0.058417 },
        { 1,  4, 15, -0.226203 },
        { 1,  5, 12, -0.143029 },
        { 1,  5, 14, -0.184651 },
        { 1,  6, 11,  0.202316 },
        { 1,  7, 10,  0.184687 },
        { 1,  8,  9,  0.226142 },
        { 1,  8, 11,  0.058377 },

        { 2,  2,  6,  0.252306 },
        { 2,  3,  7,  0.218495 },
        { 2,  4, 10,  0.184687 },
        { 2,  5, 11,  0.233626 },
        { 2,  6, 12,  0.247754 },
        { 2,  7, 13,  0.233563 },
        { 2,  8, 14,  0.184650 },

        { 3,  3,  6, -0.126184 },
        { 3,  3,  8,  0.218513 },
        { 3,  4,  9,  0.226223 },
        { 3,  4, 11, -0.058417 },
        { 3,  5, 10,  0.184687 },
        { 3,  6, 13,  0.202305 },
        { 3,  7, 12, -0.143054 },
        { 3,  7, 14,  0.184648 },
        { 3,  8, 13, -0.058422 },
        { 3,  8, 15,  0.226178 },

        { 4,  4,  6, -0.180258 },
        { 4,  5,  7,  0.156089 },
        { 4,  9, 13, -0.094056 },
        { 4, 10, 12, -0.188072 },
        { 4, 11, 13,  0.145700 },
        { 4, 11, 15,  0.094055 },

        { 5,  5,  6,  0.090137 },
        { 5,  5,  8, -0.156058 },
        { 5,  9, 14,  0.148663 },
        { 5, 10, 13,  0.115178 },
        { 5, 10, 15, -0.148675 },
        { 5, 11, 12,  0.059489 },
        { 5, 11, 14, -0.115172 },

        { 6,  6,  6,  0.180203 },
        { 6,  7,  7,  0.090088 },
        { 6,  8,  8, -0.180213 },
        { 6,  9,  9, -0.210263 },
        { 6, 11, 11,  0.126194 },
        { 6, 12, 12,  0.168179 },
        { 6, 13, 13,  0.126122 },
        { 6, 15, 15, -0.210285 },

        { 7,  7,  8,  0.156056 },
        { 7,  9, 10,  0.148697 },
        { 7, 10, 11,  0.115178 },
        { 7, 12, 13,  0.059476 },
        { 7, 13, 14,  0.115107 },
        { 7, 14, 15,  0.148656 },

        { 8,  9, 11, -0.094007 },
        { 8, 11, 11, -0.145659 },
        { 8, 12, 14, -0.188046 },
        { 8, 13, 13,  0.145650 },
        { 8, 13, 15, -0.094047 },

        { -1, -1, -1, 0.0 }
    };

    const BasisTriple kBasisTriples5[] =
    {
        { 0,  0,  0, kSH_Y_00 },
        { 0,  1,  1, kSH_Y_00 },
        { 0,  2,  2, kSH_Y_00 },
        { 0,  3,  3, kSH_Y_00 },
        { 0,  4,  4, kSH_Y_00 },
        { 0,  5,  5, kSH_Y_00 },
        { 0,  6,  6, kSH_Y_00 },
        { 0,  7,  7, kSH_Y_00 },
        { 0,  8,  8, kSH_Y_00 },
        { 0,  9,  9, kSH_Y_00 },
        { 0, 10, 10, kSH_Y_00 },
        { 0, 11, 11, kSH_Y_00 },
        { 0, 12, 12, kSH_Y_00 },
        { 0, 13, 13, kSH_Y_00 },
        { 0, 14, 14, kSH_Y_00 },
        { 0, 15, 15, kSH_Y_00 },
        { 0, 16, 16, kSH_Y_00 },
        { 0, 17, 17, kSH_Y_00 },
        { 0, 18, 18, kSH_Y_00 },
        { 0, 19, 19, kSH_Y_00 },
        { 0, 20, 20, kSH_Y_00 },
        { 0, 21, 21, kSH_Y_00 },
        { 0, 22, 22, kSH_Y_00 },
        { 0, 23, 23, kSH_Y_00 },
        { 0, 24, 24, kSH_Y_00 },

        { 1,  1,  6, -0.126155 },
        { 1,  1,  8, -0.218511 },
        { 1,  2,  5,  0.218512 },
        { 1,  3,  4,  0.218507 },
        { 1,  4, 13, -0.058392 },
        { 1,  4, 15, -0.226175 },
        { 1,  5, 12, -0.14305  },
        { 1,  5, 14, -0.18467  },
        { 1,  6, 11,  0.202    },
        { 1,  7, 10,  0.18468  },
        { 1,  8,  9,  0.22618  },
        { 1,  8, 11,  0.058402 },
        { 1,  9, 22, -0.043528 },
        { 1,  9, 24, -0.23033  },
        { 1, 10, 21, -0.075395 },
        { 1, 10, 23, -0.19947  },
        { 1, 11, 20, -0.15078  },
        { 1, 11, 22, -0.16858  },
        { 1, 12, 19,  0.19466  },
        { 1, 13, 18,  0.16858  },
        { 1, 14, 17,  0.19947  },
        { 1, 14, 19,  0.075389 },
        { 1, 15, 16,  0.23032  },
        { 1, 15, 18,  0.043523 },

        { 2,  2,  6,  0.25230  },
        { 2,  3,  7,  0.21851  },
        { 2,  4, 10,  0.18468  },
        { 2,  5, 11,  0.23359  },
        { 2,  6, 12,  0.24776  },
        { 2,  7, 13,  0.23358  },
        { 2,  8, 14,  0.18467  },
        { 2,  9, 17,  0.16287  },
        { 2, 10, 18,  0.21325  },
        { 2, 11, 19,  0.23841  },
        { 2, 12, 20,  0.24623  },
        { 2, 13, 21,  0.23840  },
        { 2, 14, 22,  0.21323  },
        { 2, 15, 23,  0.16286  },

        { 3,  3,  6, -0.12615  },
        { 3,  3,  8,  0.21851  },
        { 3,  4,  9,  0.22617  },
        { 3,  4, 11, -0.058391 },
        { 3,  5, 10,  0.18468  },
        { 3,  6, 13,  0.20229  },
        { 3,  7, 12, -0.14305  },
        { 3,  7, 14,  0.18467  },
        { 3,  8, 13, -0.058400 },
        { 3,  8, 15,  0.22618  },
        { 3,  9, 16,  0.2303   },
        { 3,  9, 18, -0.043515 },
        { 3, 10, 17,  0.19948  },
        { 3, 10, 19, -0.075395 },
        { 3, 11, 18,  0.16858  },
        { 3, 12, 21,  0.19466  },
        { 3, 13, 20, -0.15078  },
        { 3, 13, 22,  0.16857  },
        { 3, 14, 21, -0.075403 },
        { 3, 14, 23,  0.19947  },
        { 3, 15, 22, -0.043530 },
        { 3, 15, 24,  0.23033  },

        { 4,  4,  6, -0.18021  },
        { 4,  4, 20,  0.040286 },
        { 4,  4, 24, -0.23841  },
        { 4,  5,  7,  0.15608  },
        { 4,  5, 21, -0.063720 },
        { 4,  5, 23, -0.16858  },
        { 4,  6, 18,  0.15607  },
        { 4,  7, 17,  0.16859  },
        { 4,  7, 19, -0.063720 },
        { 4,  8, 16,  0.23840  },
        { 4,  9, 13, -0.094020 },
        { 4, 10, 12, -0.1880   },
        { 4, 11, 13,  0.1456   },
        { 4, 11, 15,  0.094026 },
        { 4, 16, 22, -0.075073 },
        { 4, 17, 21, -0.11262  },
        { 4, 18, 20, -0.19036  },
        { 4, 18, 24,  0.075069 },
        { 4, 19, 21,  0.14189  },
        { 4, 19, 23,  0.11262  },

        { 5,  5,  6,  0.090112 },
        { 5,  5,  8, -0.15607  },
        { 5,  5, 20, -0.16120  },
        { 5,  5, 22, -0.18022  },
        { 5,  6, 19,  0.22072  },
        { 5,  7, 18,  0.18023  },
        { 5,  8, 17,  0.16858  },
        { 5,  8, 19,  0.06371  },
        { 5,  9, 14,  0.14867  },
        { 5, 10, 13,  0.11516  },
        { 5, 10, 15, -0.14867  },
        { 5, 11, 12,  0.05946  },
        { 5, 11, 14, -0.11516  },
        { 5, 16, 23,  0.1404   },
        { 5, 17, 22,  0.13272  },
        { 5, 17, 24, -0.14046  },
        { 5, 18, 21,  0.090299 },
        { 5, 18, 23, -0.13272  },
        { 5, 19, 20,  0.044860 },
        { 5, 19, 22, -0.090302 },

        { 6,  6,  6,  0.18022  },
        { 6,  6, 20,  0.24179  },
        { 6,  7,  7,  0.090104 },
        { 6,  7, 21,  0.22072  },
        { 6,  8,  8, -0.18022  },
        { 6,  8, 22,  0.15607  },
        { 6,  9,  9, -0.21025  },
        { 6, 11, 11,  0.12615  },
        { 6, 12, 12,  0.16820  },
        { 6, 13, 13,  0.12615  },
        { 6, 15, 15, -0.2102   },
        { 6, 16, 16, -0.22936  },
        { 6, 17, 17, -0.057344 },
        { 6, 18, 18,  0.06554  },
        { 6, 19, 19,  0.13925  },
        { 6, 20, 20,  0.1638   },
        { 6, 21, 21,  0.13925  },
        { 6, 22, 22,  0.065529 },
        { 6, 23, 23, -0.057348 },
        { 6, 24, 24, -0.22937  },

        { 7,  7,  8,  0.15607  },
        { 7,  7, 20, -0.16119  },
        { 7,  7, 22,  0.18021  },
        { 7,  8, 21, -0.063727 },
        { 7,  8, 23,  0.16858  },
        { 7,  9, 10,  0.14868  },
        { 7, 10, 11,  0.11516  },
        { 7, 12, 13,  0.059467 },
        { 7, 13, 14,  0.11515  },
        { 7, 14, 15,  0.14867  },
        { 7, 16, 17,  0.14046  },
        { 7, 17, 18,  0.13273  },
        { 7, 18, 19,  0.090299 },
        { 7, 20, 21,  0.04487  },
        { 7, 21, 22,  0.090285 },
        { 7, 22, 23,  0.13271  },
        { 7, 23, 24,  0.1404   },

        { 8,  8, 20,  0.040300 },
        { 8,  8, 24,  0.23842  },
        { 8,  9, 11, -0.094032 },
        { 8, 11, 11, -0.14567  },
        { 8, 12, 14, -0.18806  },
        { 8, 13, 13,  0.14566  },
        { 8, 13, 15, -0.094033 },
        { 8, 16, 18, -0.075073 },
        { 8, 17, 19, -0.11262  },
        { 8, 19, 19, -0.14188  },
        { 8, 20, 22, -0.19036  },
        { 8, 21, 21,  0.14189  },
        { 8, 21, 23, -0.11262  },
        { 8, 22, 24, -0.075090 },

        { 9,  9, 20,  0.07692  },
        { 9, 10, 21, -0.099327 },
        { 9, 11, 22,  0.13325  },
        { 9, 11, 24,  0.11752  },
        { 9, 12, 17, -0.20355  },
        { 9, 13, 16, -0.11750  },
        { 9, 13, 18,  0.1332   },
        { 9, 14, 19, -0.099322 },

        { 10, 10, 20, -0.17952  },
        { 10, 10, 24, -0.15172  },
        { 10, 11, 21,  0.10258  },
        { 10, 11, 23, -0.067851 },
        { 10, 12, 18, -0.04442  },
        { 10, 13, 17,  0.067855 },
        { 10, 13, 19,  0.10258  },
        { 10, 14, 16,  0.15171  },
        { 10, 15, 19,  0.099323 },

        { 11, 11, 20,  0.025631 },
        { 11, 11, 22, -0.11468  },
        { 11, 12, 19,  0.099312 },
        { 11, 13, 18,  0.11469  },
        { 11, 14, 17,  0.067852 },
        { 11, 14, 19, -0.10258  },
        { 11, 15, 16, -0.11751  },
        { 11, 15, 18, -0.13325  },

        { 12, 12, 20,  0.15387  },
        { 12, 13, 21,  0.09931  },
        { 12, 14, 22, -0.044418 },
        { 12, 15, 23, -0.20355  },

        { 13, 13, 20,  0.025640 },
        { 13, 13, 22,  0.11467  },
        { 13, 14, 21,  0.10257  },
        { 13, 14, 23,  0.067841 },
        { 13, 15, 22,  0.1332   },
        { 13, 15, 24, -0.11752  },

        { 14, 14, 20, -0.17951  },
        { 14, 14, 24,  0.15172  },
        { 14, 15, 21, -0.099328 },

        { 15, 15, 20,  0.076930 },

        { 16, 16, 20,  0.10651  },
        { 16, 17, 21, -0.119    },
        { 16, 18, 22,  0.13503  },
        { 16, 19, 23, -0.11909  },

        { 17, 17, 20, -0.15979  },
        { 17, 18, 21,  0.045014 },
        { 17, 19, 22,  0.045017 },
        { 17, 19, 24,  0.11910  },

        { 18, 18, 20, -0.083713 },
        { 18, 18, 24, -0.13504  },
        { 18, 19, 21,  0.10208  },
        { 18, 19, 23, -0.045016 },

        { 19, 19, 20,  0.068462 },
        { 19, 19, 22, -0.10208  },

        { 20, 20, 20,  0.13697  },
        { 20, 21, 21,  0.068467 },
        { 20, 22, 22, -0.083689 },
        { 20, 23, 23, -0.15978  },
        { 20, 24, 24,  0.10652  },

        { 21, 21, 22,  0.10208  },
        { 21, 22, 23,  0.045006 },
        { 21, 23, 24, -0.11910  },

        { 22, 22, 24,  0.13505  },

        { -1, -1, -1, 0.0 }
    };

    // 238 coeffs out of 2925 terms = 0.0813675 non-zero
    // density for 25 elts = 9.52


    template<typename T> inline void MultiplySH(int n, const T* a, const T* b, T* c, const BasisTriple* triples)
    /// Multiplies two sets of coeffs together using given triple table.
    {
        VL_ASSERT(a != c && b != c);

        n *= n;
        for (int i = 0; i < n; i++)
            c[i] = vl_0;

        while (triples->i >= 0)
        {
            VL_ASSERT(triples->i < n && triples->j < n && triples->k < n);

            if (triples->i == triples->j)
            {
                c[triples->k] += triples->s * a[triples->i] * b[triples->j];
                if (triples->j != triples->k)
                    c[triples->i] += triples->s * (a[triples->j] * b[triples->k] + a[triples->k] * b[triples->j]);
            }
            else if (triples->j == triples->k)
            {
                c[triples->i] += triples->s * a[triples->j] * b[triples->j];
                c[triples->j] += triples->s * (a[triples->k] * b[triples->i] + a[triples->i] * b[triples->k]);
            }
            else
            {
                // i != j != k
                c[triples->i] += triples->s * (a[triples->j] * b[triples->k] + a[triples->k] * b[triples->j]);
                c[triples->j] += triples->s * (a[triples->k] * b[triples->i] + a[triples->i] * b[triples->k]);
                c[triples->k] += triples->s * (a[triples->i] * b[triples->j] + a[triples->j] * b[triples->i]);
            }

            triples++;
        }
    }
}

void SHL::MultiplySH(int numBands, const float* a, const float* b, float* c)
{
    switch (numBands)
    {
    case 1: c[0] = a[0] * b[0] * kSH_Y_00; break;
    case 2: MultiplySH(2, a, b, c, kBasisTriples2); break;
    case 3: MultiplySH(3, a, b, c, kBasisTriples3); break;
    case 4: MultiplySH(4, a, b, c, kBasisTriples4); break;
    case 5: MultiplySH(5, a, b, c, kBasisTriples5); break;
    default:
        VL_ERROR("Unhandled band count\n");
    }
}


void SHL::MultiplySH(int numBands, const Vec4f* a, const Vec4f* b, Vec4f* c)
{
    switch (numBands)
    {
    case 1: c[0] = a[0] * b[0] * kSH_Y_00; break;
    case 2: MultiplySH(2, a, b, c, kBasisTriples2); break;
    case 3: MultiplySH(3, a, b, c, kBasisTriples3); break;
    case 4: MultiplySH(4, a, b, c, kBasisTriples4); break;
    case 5: MultiplySH(5, a, b, c, kBasisTriples5); break;
    default:
        VL_ERROR("Unhandled band count\n");
    }
}

// Environment map

namespace
{
    /// Given a particular axis direction, how do we map (x, y, z) to (u, v, n)?
    const uint8_t kSwizzleTable[3][4] =
    {
        { 0, 1, 2, 0 }, // +Z
        { 1, 2, 0, 0 }, // +X
        { 2, 0, 1, 0 }, // +Y
    };

    float Uint8ToUnitFloat(uint8_t x)
    {
        return x / 255.0f;
    }

    enum
    {
        kU32RShift = 16,
        kU32GShift = 8,
        kU32BShift = 0,
        kU32AShift = 24
    };

    inline bool IsPowerOfTwo(uint32_t a)
    {
        return (a & (a - 1)) == 0;
    }


    inline void FindHemiEnvThetaPhi(float domega0, float fx, float fy, float fr2, float* thetaOut, float* phiOut, float* domega)
    {
        // this mapping comes from devebec's probes -- entire range of theta [0, pi] is mapped to r.
        float& theta = *thetaOut;
        float& phi   = *phiOut;

        theta = sqrtf(fr2) * vl_pi;
        phi   = atan2f(fy, fx);
        *domega = domega0 * sinf(theta) / theta;
    }
}

void SHL::FindSHCoeffsFromHemiEnvMap(const ImageData32* image, int n, float* coeffsR, float* coeffsG, float* coeffsB)
{
    int height = image->mHeight;
    int width  = image->mWidth;

    const uint32_t* data = image->mData;

    float domega0 = sqr(vl_twoPi / float(width));

    for (int y = 0; y < height; y++)
    {
        float fy = 2.0f * (y + 0.5f) / float(height) - 1.0f;

        const uint32_t* span = data + width * y;

        for (int x = 0; x < width; x++)
        {
            float fx = 2.0f * (x + 0.5f) / float(width) - 1.0f;
            float fr2 = sqr(fx) + sqr(fy);

            if (fr2 > 1.0f)
                continue;

            float theta, phi, domega;
            FindHemiEnvThetaPhi(domega0, fx, fy, fr2, &theta, &phi, &domega);

            float r = Uint8ToUnitFloat((span[x] >> kU32RShift) & 0xFF);
            float g = Uint8ToUnitFloat((span[x] >> kU32GShift) & 0xFF);
            float b = Uint8ToUnitFloat((span[x] >> kU32BShift) & 0xFF);

            int l = 0;

            for (int i = 0; i < n; i++)
            {
                if (i >= sqr(l + 1))
                    l++;
                int m = i - l * (l + 1);

                float shScale = SH(l, m, theta, phi) * domega;
                coeffsR[i] += r * shScale;
                coeffsG[i] += g * shScale;
                coeffsB[i] += b * shScale;
            }
        }
    }
}

void SHL::FindSHCoeffsFromHemiEnvMap(const ImageData32* image, int n, Vec4f* coeffs)
{
    int height = image->mHeight;
    int width  = image->mWidth;

    const uint32_t* data  = image->mData;

    float domega0 = sqr(vl_twoPi / float(width));

    for (int y = 0; y < height; y++)
    {
        float fy = -2.0f * (y + 0.5f) / float(height) + 1.0f;

        const uint32_t* span = data + width * y;

        for (int x = 0; x < width; x++)
        {
            float fx = 2.0f * (x + 0.5f) / float(width) - 1.0f;
            float fr2 = sqr(fx) + sqr(fy);

            if (fr2 > 1.0f)
                continue;

            float theta, phi, domega;
            FindHemiEnvThetaPhi(domega0, fx, fy, fr2, &theta, &phi, &domega);

            float r = Uint8ToUnitFloat((span[x] >> kU32RShift) & 0xFF);
            float g = Uint8ToUnitFloat((span[x] >> kU32GShift) & 0xFF);
            float b = Uint8ToUnitFloat((span[x] >> kU32BShift) & 0xFF);
            float a = Uint8ToUnitFloat((span[x] >> kU32AShift) & 0xFF);

            int l = 0;

            for (int i = 0; i < n; i++)
            {
                if (i >= sqr(l + 1))
                    l++;
                int m = i - l * (l + 1);

                float shScale = SH(l, m, theta, phi) * domega;
                coeffs[i][0] += r * shScale;
                coeffs[i][1] += g * shScale;
                coeffs[i][2] += b * shScale;
                coeffs[i][3] += a * shScale;
            }
        }
    }
}

void SHL::FindSHCoeffsFromHemiEnvMap(const ImageData48* image, int n, float* coeffsR, float* coeffsG, float* coeffsB)
{
    int height  = image->mHeight;
    int width   = image->mWidth;
    int strideY = width * 3;

    const uint16_t* data  = image->mData;

    float domega0 = sqr(vl_twoPi / float(width)) * (1.0f / 65535.0f);

    for (int y = 0; y < height; y++)
    {
        float fy = 2.0f * (y + 0.5f) / float(height) - 1.0f;

        const uint16_t* span = data + strideY * y;

        for (int x = 0; x < width; x++)
        {
            float fx = 2.0f * (x + 0.5f) / float(width) - 1.0f;
            float fr2 = sqr(fx) + sqr(fy);

            if (fr2 > 1.0f)
                continue;

            float theta, phi, domega;
            FindHemiEnvThetaPhi(domega0, fx, fy, fr2, &theta, &phi, &domega);

            const uint16_t* pixel = span + 3 * x;
            float r = pixel[x + 0];
            float g = pixel[x + 1];
            float b = pixel[x + 2];

            int l = 0;

            for (int i = 0; i < n; i++)
            {
                if (i >= sqr(l + 1))
                    l++;
                int m = i - l * (l + 1);

                float shScale = SH(l, m, theta, phi) * domega;
                coeffsR[i] += r * shScale;
                coeffsG[i] += g * shScale;
                coeffsB[i] += b * shScale;
            }
        }
    }
}

void SHL::FindSHCoeffsFromHemiEnvMap(const ImageData48* image, int n, Vec4f* coeffs)
{
    int height      = image->mHeight;
    int width       = image->mWidth;
    int strideY     = width * 3;

    const uint16_t* data  = image->mData;

    float domega0 = sqr(vl_twoPi / float(width)) * (1.0f / 65535.0f);

    for (int y = 0; y < height; y++)
    {
        float fy = -2.0f * (y + 0.5f) / float(height) + 1.0f;

        const uint16_t* span = data + strideY * y;

        for (int x = 0; x < width; x++)
        {
            float fx = 2.0f * (x + 0.5f) / float(width) - 1.0f;
            float fr2 = sqr(fx) + sqr(fy);

            if (fr2 > 1.0f)
                continue;

            float theta, phi, domega;
            FindHemiEnvThetaPhi(domega0, fx, fy, fr2, &theta, &phi, &domega);

            const uint16_t* pixel = span + 3 * x;
            float r = pixel[x + 0];
            float g = pixel[x + 1];
            float b = pixel[x + 2];

            int l = 0;

            for (int i = 0; i < n; i++)
            {
                if (i >= sqr(l + 1))
                    l++;
                int m = i - l * (l + 1);

                float shScale = SH(l, m, theta, phi) * domega;
                coeffs[i][0] += r * shScale;
                coeffs[i][1] += g * shScale;
                coeffs[i][2] += b * shScale;
            }
        }
    }
}

namespace
{
    /*
        domega for a cube map is:
            d = 2/w
            domega = d^2 / (x^2 + y^2 + 1) ^ 3/2

        it may seem odd that this will give us a total factor of 4pi steradians, but:
            integral(-1, 1, x, integral(-1, 1, y, (x^2 + y^2 + 1) ^ -3/2)) = pi/6.

    */

    struct HDRFaceInfo
    {
        int   mTX, mTY;     // tile
        float mUX, mUY;     // u = mUX * x + mUY * y
        float mVX, mVY;
    };

    const HDRFaceInfo kFaceTable[6] =
    {
        // fx, fy, sx, sy
        { 3, 1,   -1,  0,   0,  1  },  // PZ
        { 1, 1,   -1,  0,   0,  1  },  // NZ
        { 2, 1,    0,  1,   1,  0  },  // PX
        { 0, 1,    0, -1,  -1,  0  },  // NX
        { 1, 0,    0,  1,   1,  0  },  // PY
        { 1, 2,    0,  1,   1,  0  }   // NY
    };
}

void SHL::FindSHCoeffsFromHDRCubeMap(const ImageData48* image, int numBands, Vec4f* coeffs)
{
    int height = image->mHeight;
    int width  = image->mWidth;

    int faceSize = image->mWidth / 4;

    if (!IsPowerOfTwo(faceSize) || faceSize * 4 != width || faceSize * 3 != height)
    {
        VL_ERROR("invalid cube map!\n");
        return;
    }

    Vec4f c(0, 0, 0, 0);
    const uint16_t* data  = image->mData;

    const int strideX = 3;
    int       strideY = width * strideX;

    float domega0 = (1.0f / 65535.0f) * 4.0f / sqr(faceSize);

    for (int face = 0; face < 6; face++)
    {
        const uint16_t* faceData = data + strideY * faceSize * kFaceTable[face].mTY + strideX * faceSize * kFaceTable[face].mTX;

        for (int y = 0; y < faceSize; y++)
        {
            float fy = -2 * (y + 0.5f) / float(faceSize) + 1;

            const uint16_t* pixel = faceData + strideY * y;

            for (int x = 0; x < faceSize; x++)
            {
                float fx = 2 * (x + 0.5f) / float(faceSize) - 1;

                float fr2 = sqr(fx) + sqr(fy) + 1.0f;
                float dn = fr2 * sqrtf(fr2);
                float domega = domega0 / dn;     // 4/w^2  *  (x^2 + y^2 + 1) ^ -3/2

                float u = kFaceTable[face].mUX * fx + kFaceTable[face].mUY * fy;
                float v = kFaceTable[face].mVX * fx + kFaceTable[face].mVY * fy;

                VL_ASSERT(pixel >= data && pixel < data + width * height * 3);

                c[0] = pixel[0] * domega;
                c[1] = pixel[1] * domega;
                c[2] = pixel[2] * domega;

                Vec3f dir;
                float invLen = 1.0f / sqrtf(sqr(u) + sqr(v) + 1);

                const uint8_t* swizzle = kSwizzleTable[face >> 1];
                float w = (face & 1) ? -invLen : invLen;

                dir[swizzle[0]] = u * w;
                dir[swizzle[1]] = v * invLen;
                dir[swizzle[2]] = w;

                AddSHSample(c, dir, numBands, coeffs);

                pixel += strideX;
            }
        }
    }
}

#ifdef REFERENCE
// This is the unoptimized version, to show more clearly what we're doing. Which
// is just using the GatedSpot light model to represent the light as a sphere,
// and accumulating everything into the destination.
void SHL::AddSphereLighting(Vec3f pos, int numLights, SphereLight* lights, Vec4f* coeffs)
{
    float zcoeffs[7];
    Vec4f zcoeffsColour[7];

    for (int i = 0; i < numLights; i++)
    {
        const Vec4f& colour = lights[i].mColourAndIntensity;

        Vec3f dir(pos - (Vec3f&) lights[i].mPositionAndSize);
        float r2 = sqrlen(dir);

        float strength = lights[i].mColourAndIntensity[3];

        float a = lights[i].mPositionAndSize[3];
        float a2 = sqr(a);

        // The zcoeffs below are for f(z) = 1 for z > t.
        // So t is just the z component of (a, 0, r) normalized.
        float t = sqrtf(r2 / (a2 + r2));
        dir /= sqrtf(r2);

        CalcGatedSpotZH7(t, zcoeffs);

        for (int i = 0; i < 5; i++)
            zcoeffsColour[i] = (zcoeffs[i] * strength) * colour;

        RotateZHToSHAdd(dir, 5, zcoeffsColour, coeffs);
    }
}
#endif

void SHL::AddSphereLighting(Vec3f pos, float scale, int numLights, const SphereLight* lights, Vec4f* coeffs)
{
    // we have no order dependency, and we might as well just do all lights, as it's almost
    // certainly quicker than calculating which are the N most important.
    // in particular, we tailor how many terms we add in to the strength of the light,
    // so distant lights effectively only add to the base ambient/x/y/z terms.

    Vec4f zcoeffs[4];

    for (int i = 0; i < numLights; i++)
    {
        const Vec4f& colour = lights[i].mColourAndIntensity;

        Vec3f dir((Vec3f&) lights[i].mPositionAndSize - pos);
        float r2 = sqrlen(dir);

        float strength = lights[i].mColourAndIntensity[3] * scale;

        float a = lights[i].mPositionAndSize[3];
        float a2 = sqr(a);

        // See CalcGatedSpotZH7
        // The zcoeffs below are for f(z) = 1 for z > t.
        // So t is just the z component of (a, 0, r) normalized.
        float t2 = r2 / (a2 + r2);
        float t = sqrtf(t2);
        dir /= sqrtf(r2) + 1e-6f;

        const int bands = 4;
        float z0 = strength * (vl_twoPi * kZH_Y[0]) * (1.0f - t);
        float z1 = strength * (vl_twoPi * kZH_Y[1]) * 0.5f * (1.0f - t2);
        float z2 = strength * (vl_twoPi * kZH_Y[2]) * t * (1.0f - t2);
        float z3 = strength * (vl_twoPi * kZH_Y[3]) * 0.25f * (t2 * (6.0f - 5.0f * t2) - 1.0f);

        zcoeffs[0] = colour * z0;
        zcoeffs[1] = colour * z1;
        zcoeffs[2] = colour * z2;
        zcoeffs[3] = colour * z3;

        RotateZHToSHAdd(dir, bands, zcoeffs, coeffs);
    }
}
