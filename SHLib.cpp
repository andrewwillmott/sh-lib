//------------------------------------------------------------------------------
// SH Utility Library
//------------------------------------------------------------------------------

#include "SHLib.h"

namespace
{
    /// Given a particular axis direction, how do we map (x, y, z) to (u, v, n)?
    const uint8_t kSwizzleTable[3][4] =
    {
        { 0, 1, 2, 0 },        // +Z
        { 1, 2, 0, 0 },        // +X
        { 2, 0, 1, 0 },        // +Y
    };

    int32_t FloorToInt32(float s)
    {
        return int(floorf(s));
    }

    inline float Saturate(float x)
    {
        if (x < 0.0f)
            return 0.0f;
        if (x > 1.0f)
            return 1.0f;
        return x;
    }

    float Uint8ToUnitFloat(uint8_t x)
    {
        return x / 255.0f;
    }

    uint8_t ClampUnitFloatToUint8(float x)
    {
        return uint8_t(Saturate(x) * 255.0f);
    }

    enum
    {
        kU32RShift = 16,
        kU32GShift = 8,
        kU32BShift = 0,
        kU32AShift = 24
    };
    
    uint32_t ColourToU32(Vec4f c)
    {
        return (ClampUnitFloatToUint8(c[0])  << kU32RShift)
             | (ClampUnitFloatToUint8(c[1])  << kU32GShift)
             | (ClampUnitFloatToUint8(c[2])  << kU32BShift)
             | (ClampUnitFloatToUint8(c[3])  << kU32AShift);
    }

    inline bool IsPowerOfTwo(uint32_t a)
    {
        return (a & (a - 1)) == 0;
    }
}

namespace
{
    const float pi = vl_pi;
    
    // The normalization constants themselves are given by
    //   sqrt(((2l + 1) (l - m)!)  / (4 pi (l + m)!))
    // However, we also roll in any overall multipliers from the
    // basis function itself. (E.g., the b.f. for Y_22 is (x^2 - y^2) / 2.)
    
    // normalization constants                            basis function:
    const float kSH_Y_00  = sqrtf( 1 / ( 4 * pi));          // 1
    
    // 3, 3 -> 1, 1
    const float kSH_Y_1_1 = sqrtf( 3 / ( 4 * pi));          // y
    const float kSH_Y_10  = sqrtf( 3 / ( 4 * pi));          // z
    const float kSH_Y_11  = sqrtf( 3 / ( 4 * pi));          // x

    // 15, 15, 5 -> 3, 3, 1
    const float kSH_Y_2_2 = sqrtf(15 / ( 4 * pi));          // xy
    const float kSH_Y_2_1 = sqrtf(15 / ( 4 * pi));          // yz
    const float kSH_Y_20  = sqrtf( 5 / (16 * pi));          // (3z^2 - 1)
    const float kSH_Y_21  = sqrtf(15 / ( 4 * pi));          // xz
    const float kSH_Y_22  = sqrtf(15 / (16 * pi));          // (x^2 - y^2)

    // 35/2  105 21/2  7 -> 5/2, 15, 3/2, 1
    const float kSH_Y_3_3 = sqrtf( 35 / (32 * pi));         // y (3x^2 - y^2)
    const float kSH_Y_3_2 = sqrtf(105 / ( 4 * pi));         // z xy
    const float kSH_Y_3_1 = sqrtf( 21 / (32 * pi));         // y (5z^2 - 1)
    const float kSH_Y_30  = sqrtf(  7 / (16 * pi));         // z (5z^2 - 3)
    const float kSH_Y_31  = sqrtf( 21 / (32 * pi));         // x (5z^2 - 1)
    const float kSH_Y_32  = sqrtf(105 / (16 * pi));         // z (x^2 - y^2)
    const float kSH_Y_33  = sqrtf( 35 / (32 * pi));         // x (x^2 - 3y^2)

    // 315, 315/2, 45, 45/2, 9 -> 35, 35/2, 5, 5/2, 1
    const float kSH_Y_4_4 = sqrtf(315 / ( 16 * pi));        // (x^2 - y^2)  xy
    const float kSH_Y_4_3 = sqrtf(315 / ( 32 * pi));        // (3x^2 - y^2) yz
    const float kSH_Y_4_2 = sqrtf( 45 / ( 16 * pi));        // (7z^2 - 1)   xy
    const float kSH_Y_4_1 = sqrtf( 45 / ( 32 * pi));        // (7z^2 - 3)   yz
    const float kSH_Y_40  = sqrtf(  9 / (256 * pi));        // (3 - 30z^2 + 35z^4)
    const float kSH_Y_41  = sqrtf( 45 / ( 32 * pi));        // (7z^2 - 3)   xz
    const float kSH_Y_42  = sqrtf( 45 / ( 64 * pi));        // (7z^2 - 1)   (x^2 -y^2)
    const float kSH_Y_43  = sqrtf(315 / ( 32 * pi));        // (x^2 - 3y^2) xz
    const float kSH_Y_44  = sqrtf(315 / (256 * pi));        // (x^4 - 6 x^2 y^2 + y^4)


    const float kSH_Y_1x  = sqrtf(3 / ( 4 * pi));          // x, y, z
        
    // The 'central' (z-symmetric, m=0) basis function set is a legendre poly in z:
    //   1, z, (3z^2 - 1), (5z^3 - 3z), 3 - 30z^2 + 35z^4 ...
    // Sometimes referred to as zonal harmonics.

    // The ZH normalization constants are just sqrtf( (2l + 1) / (4 pi) )
    // As above we also roll in any constant factor from the polynomial.
    // The constant factor for these is easy -- just evaluate with z = 1
    //                                                  N       P(z) / N
    const float kZH_Y_0 = sqrtf( 1 / (   4 * pi));  //         1
    const float kZH_Y_1 = sqrtf( 3 / (   4 * pi));  //         z
    const float kZH_Y_2 = sqrtf( 5 / (  16 * pi));  // 1/2     (3 z^2 - 1)
    const float kZH_Y_3 = sqrtf( 7 / (  16 * pi));  // 1/2     (5 z^3 - 3 z)
    const float kZH_Y_4 = sqrtf( 9 / ( 256 * pi));  // 1/8     (35 z^4 - 30 z^2 + 3)
    const float kZH_Y_5 = sqrtf(11 / ( 256 * pi));  // 1/8     (63 z^5 - 70 z^3 + 15 z)
    const float kZH_Y_6 = sqrtf(13 / (1024 * pi));  // 1/16    (231 z^6 - 315 z^4 + 105 z^2 - 5)
    

    // diffuse reflectance convolution weights.
    // Multiply bands by these constants to convert incident radiance to exit
    // radiance for a diffuse surface.
    const float kSH_A0 = 1.0f;
    const float kSH_A1 = (2.0f / 3.0f);
    const float kSH_A2 = (1.0f / 4.0f);
    const float kSH_A3 = 0.0f;
    const float kSH_A4 = -(1.0f / 24.0f);

    // Integrating f(z) over the sphere is:
    //   integral[-1;1]  f(z) . 2pi dz


    // The VS compiler doesn't convert sqrt to constants
    const float kSqrt03_02    = sqrtf( 3.0 /  2.0);
    const float kSqrt01_03    = sqrtf( 1.0 /  3.0);
    const float kSqrt02_03    = sqrtf( 2.0 /  3.0);
    const float kSqrt04_03    = sqrtf( 4.0 /  3.0);
    const float kSqrt01_04    = sqrtf( 1.0 /  4.0);
    const float kSqrt03_04    = sqrtf( 3.0 /  4.0);
    const float kSqrt01_05    = sqrtf( 1.0 /  5.0);
    const float kSqrt03_05    = sqrtf( 3.0 /  5.0);
    const float kSqrt06_05    = sqrtf( 6.0 /  5.0);
    const float kSqrt08_05    = sqrtf( 8.0 /  5.0);
    const float kSqrt09_05    = sqrtf( 9.0 /  5.0);
    const float kSqrt05_06    = sqrtf( 5.0 /  6.0);
    const float kSqrt01_06    = sqrtf( 1.0 /  6.0);
    const float kSqrt03_08    = sqrtf( 3.0 /  8.0);
    const float kSqrt05_08    = sqrtf( 5.0 /  8.0);
    const float kSqrt07_08    = sqrtf( 7.0 /  8.0);
    const float kSqrt09_08    = sqrtf( 9.0 /  8.0);
    const float kSqrt05_09    = sqrtf( 5.0 /  9.0);
    const float kSqrt08_09    = sqrtf( 8.0 /  9.0);

    const float kSqrt01_10    = sqrtf( 1.0 / 10.0);
    const float kSqrt03_10    = sqrtf( 3.0 / 10.0);
    const float kSqrt01_12    = sqrtf( 1.0 / 12.0);
    const float kSqrt04_15    = sqrtf( 4.0 / 15.0);
    const float kSqrt01_16    = sqrtf( 1.0 / 16.0);
    const float kSqrt07_16    = sqrtf( 7.0 / 16.0);
    const float kSqrt15_16    = sqrtf(15.0 / 16.0);
    const float kSqrt01_18    = sqrtf( 1.0 / 18.0);
    const float kSqrt03_25    = sqrtf( 3.0 / 25.0);
    const float kSqrt14_25    = sqrtf(14.0 / 25.0);
    const float kSqrt15_25    = sqrtf(15.0 / 25.0);
    const float kSqrt18_25    = sqrtf(18.0 / 25.0);
    const float kSqrt01_32    = sqrtf( 1.0 / 32.0);
    const float kSqrt03_32    = sqrtf( 3.0 / 32.0);
    const float kSqrt15_32    = sqrtf(15.0 / 32.0);
    const float kSqrt21_32    = sqrtf(21.0 / 32.0);
    const float kSqrt01_50    = sqrtf( 1.0 / 50.0);
    const float kSqrt03_50    = sqrtf( 3.0 / 50.0);
    const float kSqrt21_50    = sqrtf(21.0 / 50.0);
        
    uint32_t kFactorialTable[] =
    {
        1,
        1,
        2,
        6,
        24,
        120,
        720,
        5040,
        40320,
        362880,
        3628800,
        39916800,
        479001600
    };

    inline float fact(int n)
    {
        if (n < 13)
            return float(kFactorialTable[n]);
        
        float result = float(kFactorialTable[12]);
        while (n >= 13)
            result *= n--;

        return result;
    }

    inline float dp3(const float* a, const float* b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    inline Vec4f dp3(const Vec4f* a, const float* b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    inline float dp(int n, const float* a, const float* b)
    {
        float result = (*a) * (*b);
        
        while (--n > 0)
        {
            a++;
            b++;
            result += (*a) * (*b);
        }
        
        return result;
    }

    inline Vec4f dp(int n, const Vec4f* a, const float* b)
    {
        Vec4f result = (*a) * (*b);
        
        while (--n > 0)
        {
            a++;
            b++;
            result += (*a) * (*b);
        }
        
        return result;
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

    inline float SH_K(int l, int pm)
    {
        if (pm == 0)    // factorials cancel out
            return sqrtf(float(2 * l + 1) / (4.0f * vl_pi));

        float top    = float(2 * l + 1) * fact(l - pm);
        float bottom = (4.0f * vl_pi)   * fact(l + pm);

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

    float SH(int l, int m, float theta, float phi)
    {
        const float kRoot2 = sqrtf(2.0f);
        
        int pm = abs(m);
        float term1 = SH_K(l, pm) * SH_P(l, pm, cosf(theta));

        phi += vl_pi;
        
        if (m == 0)
            return term1;
        else if (m > 0)
            return term1 * kRoot2 * cosf(pm * phi);
        else
            return term1 * kRoot2 * sinf(pm * phi);
    }

    // Chebyshev recurrence relations for sin(nx)/cos(nx) given sin(x)/cos(x)
    inline float SinN(int n, float s, float c)
    {
        float sn0;
        float sn1 = 0;
        float sn2 = s;

        for (int i = 2; i < n; i++)
        {
            sn0 = sn1;
            sn1 = sn2;
            sn2 = 2 * c * sn1 - sn0;
        }

        return sn2;
    }

    inline float CosN(int n, float s, float c)
    {
        float cn0;
        float cn1 = 0;
        float cn2 = c;

        for (int i = 2; i < n; i++)
        {
            cn0 = cn1;
            cn1 = cn2;
            cn2 = 2 * c * cn1 - cn0;
        }

        return cn2;
    }

    float SH(int l, int m, const Vec3f& dir)  // dir must be normalised
    {
        const float kRoot2 = sqrtf(2.0f);
        
        float cosTheta = dir[2];
        float invSinTheta = 1.0f / (sqrtf(1.0f - cosTheta * cosTheta) + 1e-6f * cosTheta); // epsilon to avoid 0 divide
        
        // minus signs recapitulate += pi from original.
        float cosPhi = -dir[0] * invSinTheta;
        float sinPhi = -dir[1] * invSinTheta;

        int pm = abs(m);
        float term1 = SH_K(l, pm) * SH_P(l, pm, cosTheta);

        if (m == 0)
            return term1;
        else if (m > 0)
            return term1 * kRoot2 * CosN(pm, -sinPhi, -cosPhi);
        else
            return term1 * kRoot2 * SinN(pm, -sinPhi, -cosPhi);
    }
}

Vec4f SHL::SampleColour(const Vec4f* coeffs, int numBands, const Vec3f& v)
{
    CL_ASSERT(numBands >= 2 && numBands <= 5);
    
    Vec4f s(coeffs[0] * kSH_Y_00);
    
    const float& x  = v[0];
    const float& y  = v[1];
    const float& z  = v[2];

    float  x2 = sqr(x);
    float  y2 = sqr(y);
    float  z2 = sqr(z);
    
    s += coeffs[1] * kSH_Y_1x * y;
    s += coeffs[2] * kSH_Y_1x * z;
    s += coeffs[3] * kSH_Y_1x * x;

    if (numBands < 3)
        return s;
    
    float  xy = x * y;
    float  yz = y * z;
    float  xz = x * z;
    
    s += coeffs[4] * kSH_Y_2_2 * xy;
    s += coeffs[5] * kSH_Y_2_1 * yz;
    s += coeffs[6] * kSH_Y_20  * (3.0f * z2 - 1.0f);
    s += coeffs[7] * kSH_Y_21  * xz;
    s += coeffs[8] * kSH_Y_22  * (x2 - y2);

    if (numBands < 4)
        return s;
    
    s += coeffs[ 9] * kSH_Y_3_3 * (3.0f * x2 - y2)   * y;
    s += coeffs[10] * kSH_Y_3_2 * xy                 * z;
    s += coeffs[11] * kSH_Y_3_1 * (5.0f * z2 - 1.0f) * y;
    s += coeffs[12] * kSH_Y_30  * (5.0f * z2 - 3.0f) * z;
    s += coeffs[13] * kSH_Y_31  * (5.0f * z2 - 1.0f) * x;
    s += coeffs[14] * kSH_Y_32  * (x2 - y2)          * z;
    s += coeffs[15] * kSH_Y_33  * (x2 - 3.0f * y2)   * x;
            
    if (numBands < 5)
        return s;

    s += coeffs[16] * kSH_Y_4_4 * (x2 - y2)          * xy;
    s += coeffs[17] * kSH_Y_4_3 * (3.0f * x2 - y2)   * yz;
    s += coeffs[18] * kSH_Y_4_2 * (7.0f * z2 - 1.0f) * xy;
    s += coeffs[19] * kSH_Y_4_1 * (7.0f * z2 - 3.0f) * yz;
    s += coeffs[20] * kSH_Y_40  * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
    s += coeffs[21] * kSH_Y_41  * (7.0f * z2 - 3.0f) * xz;
    s += coeffs[22] * kSH_Y_42  * (7.0f * z2 - 1.0f) * (x2 - y2);
    s += coeffs[23] * kSH_Y_43  * (x2 - 3.0f * y2)   * xz;
    s += coeffs[24] * kSH_Y_44  * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);
    
    return s;
}


void SHL::AddColourSample(const Vec4f& s, const Vec3f& v, int numBands, Vec4f* coeffs)
{
    CL_ASSERT(numBands >= 3 && numBands <= 5);
    
    coeffs[0] += s * kSH_Y_00;
    
    const float& x = v[0];
    const float& y = v[1];
    const float& z = v[2];
    
    coeffs[1] += s * kSH_Y_1x * y;
    coeffs[2] += s * kSH_Y_1x * z;
    coeffs[3] += s * kSH_Y_1x * x;

    if (numBands < 3)
        return;
    
    float x2 = sqr(x);
    float y2 = sqr(y);
    float z2 = sqr(z);
    float xy = x * y;
    float yz = y * z;
    float xz = x * z;
    
    coeffs[4] += s * kSH_Y_2_2 * xy;
    coeffs[5] += s * kSH_Y_2_1 * yz;
    coeffs[6] += s * kSH_Y_20  * (3.0f * z2 - 1.0f);
    coeffs[7] += s * kSH_Y_21  * xz;
    coeffs[8] += s * kSH_Y_22  * (x2 - y2);

    if (numBands < 4)
        return;
    
    coeffs[ 9] += s * kSH_Y_3_3 * (3.0f * x2 - y2) * y;
    coeffs[10] += s * kSH_Y_3_2 * xy            * z;
    coeffs[11] += s * kSH_Y_3_1 * (5.0f * z2 - 1.0f)  * y;
    coeffs[12] += s * kSH_Y_30  * (5.0f * z2 - 3.0f)  * z;
    coeffs[13] += s * kSH_Y_31  * (5.0f * z2 - 1.0f)  * x;
    coeffs[14] += s * kSH_Y_32  * (x2 - y2)     * z;
    coeffs[15] += s * kSH_Y_33  * (x2 - 3.0f * y2) * x;
            
    if (numBands < 5)
        return;

    coeffs[16] += s * kSH_Y_4_4 * (x2 - y2)     * xy;
    coeffs[17] += s * kSH_Y_4_3 * (3.0f * x2 - y2) * yz;
    coeffs[18] += s * kSH_Y_4_2 * (7.0f * z2 - 1.0f)  * xy;
    coeffs[19] += s * kSH_Y_4_1 * (7.0f * z2 - 3.0f)  * yz;
    coeffs[20] += s * kSH_Y_40  * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
    coeffs[21] += s * kSH_Y_41  * (7.0f * z2 - 3.0f)  * xz;
    coeffs[22] += s * kSH_Y_42  * (7.0f * z2 - 1.0f)  * (x2 - y2);
    coeffs[23] += s * kSH_Y_43  * (x2 - 3.0f * y2) * xz;
    coeffs[24] += s * kSH_Y_44  * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);
}

void SHL::AddReflectedColourSample(const Vec4f& sIn, const Vec3f& v, int numBands, Vec4f* coeffs)
{
    CL_ASSERT(numBands >= 3 && numBands <= 5);
    
    Vec4f s(2.0f * sIn);

    coeffs[0] += s * kSH_Y_00;
    
    const float& x  = v[0];
    const float& y  = v[1];
    
    coeffs[1] += s * kSH_Y_1x * y;
    coeffs[3] += s * kSH_Y_1x * x;

    if (numBands < 3)
        return;
    
    const float& z  = v[2];
    float x2 = sqr(x);
    float y2 = sqr(y);
    float z2 = sqr(z);
    float xy = x * y;
    
    coeffs[4] += s * kSH_Y_2_2 * xy;
    coeffs[6] += s * kSH_Y_20  * (3.0f * z2 - 1.0f);
    coeffs[8] += s * kSH_Y_22  * (x2 - y2);

    if (numBands < 4)
        return;
    
    coeffs[ 9] += s * kSH_Y_3_3 * (3.0f * x2 -        y2) * y;
    coeffs[11] += s * kSH_Y_3_1 * (5.0f * z2 - 1.0f     )  * y;
    coeffs[13] += s * kSH_Y_31  * (5.0f * z2 - 1.0f     )  * x;
    coeffs[15] += s * kSH_Y_33  * (       x2 - 3.0f * y2) * x;
            
    if (numBands < 5)
        return;

    coeffs[16] += s * kSH_Y_4_4 * (x2 - y2)           * xy;
    coeffs[18] += s * kSH_Y_4_2 * (7.0f * z2 - 1.0f)  * xy;
    coeffs[20] += s * kSH_Y_40  * (3.0f - 30.0f * z2 + 35.0f * z2 * z2);
    coeffs[22] += s * kSH_Y_42  * (7.0f * z2 - 1.0f)  * (x2 - y2);
    coeffs[24] += s * kSH_Y_44  * (x2 * x2 - 6.0f * x2 * y2 + y2 * y2);
}


void SHL::CalcCosPowerZH7(int n, float zcoeffs[7])
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

void SHL::CalcCosPowerSatZH7(float n, float zcoeffs[7])
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

void SHL::GetHemisphereZH7(float zcoeffs[7])
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

void SHL::CalcGatedSpotZH7(float t, float zcoeffs[7])
{
    // Just the integral of P_i(z) 2pi dz over [t, 1].
    float t2 = sqr(t);
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;
    float t6 = t3 * t3;
    float t7 = t4 * t3;
    
    zcoeffs[0] = 1 - t;
    zcoeffs[1] = 0.5f * (1.0f - t2);
    zcoeffs[2] = t - t3;
    zcoeffs[3] = 0.25f * (-1.0f + t2 * 6.0f - 5.0f * t4);
    zcoeffs[4] = -3.0f * t + 10.0f * t3 - 7.0f * t5;
    zcoeffs[5] = 0.5f * (1.0f - 21.0f * t6 + 35.0f * t4 - 15.0f * t2);
    zcoeffs[6] = 5.0f * t - 35.0f * t3 + 63.0f * t5 - 33.0f * t7;

    // apply norm constants
    zcoeffs[0] *= vl_twoPi * kZH_Y_0;    
    zcoeffs[1] *= vl_twoPi * kZH_Y_1;
    zcoeffs[2] *= vl_twoPi * kZH_Y_2;
    zcoeffs[3] *= vl_twoPi * kZH_Y_3;
    zcoeffs[4] *= vl_twoPi * kZH_Y_4;
    zcoeffs[5] *= vl_twoPi * kZH_Y_5;
    zcoeffs[6] *= vl_twoPi * kZH_Y_6;
}

void SHL::CalcGatedCosZH7(float t, float zcoeffs[7])
{
    // Just the integral of P_i(z) 2pi dz over [t, 1].
    float t2 = sqr(t);
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;
    float t6 = t3 * t3;
    float t7 = t4 * t3;
    float t8 = t4 * t4;
    
    zcoeffs[0] = 0.5f * (1.0f - t2);
    zcoeffs[1] = (1.0f - t3) / 3.0f;
    zcoeffs[2] = 0.25f * (1.0f - 3.0f * t4 - 2.0f * t2);
    zcoeffs[3] = t3 - t5;
    zcoeffs[4] = -1.0f / 6.0f * (1.0f + 9.0f * t2 - 45.0f * t4 + 35.0f * t6);
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

void SHL::CalcGatedCosBands5(float t, float bandMult[5])
{
    // Just the integral of P_i(z) 2pi dz over [t, 1].
    float t2 = sqr(t);
    float t3 = t2 * t;
    float t4 = t2 * t2;
    float t5 = t3 * t2;
    float t6 = t3 * t3;
    
    bandMult[0] = 0.5f * (1.0f - t2);
    bandMult[1] = (1.0f - t3) / 3.0f;
    bandMult[2] = 0.25f * (1.0f - 3.0f * t4 - 2.0f * t2);
    bandMult[3] = t3 - t5;
    bandMult[4] = -1.0f / 6.0f * (1.0f + 9.0f * t2 - 45.0f * t4 + 35.0f * t6);

    float invB0 = 1.0f / bandMult[0];
    
    bandMult[0] = 1.0f;
    bandMult[1] *= invB0;
    bandMult[2] *= invB0 * 0.5f;
    bandMult[3] *= invB0 * 0.5f;
    bandMult[4] *= invB0 * 0.125f;
}

void SHL::CalcHGPhaseZH(float g, float weight, int n, float zcoeffs[])
{
    // The Heyney-Greenstein phase function turns out to be a
    // simple power series in g when decomposed into spherical harmonics.
    for (int i = 0; i < n; i++)
    {
        zcoeffs[i] = sqrtf(4.0f * vl_pi * float(i * 2 + 1)) * weight;
        weight *= g;
    }
}

void SHL::CalcNormalizedHGPhaseZH(float g, float weight, int n, float zcoeffs[])
{
    // normalize so that max of function is always 1.
    weight *= sqr(1 - fabsf(g)) / (1 + fabsf(g));

    for (int i = 0; i < n; i++)
    {
        zcoeffs[i] = sqrtf(4.0f * vl_pi * float(i * 2 + 1)) * weight;
        weight *= g;
    }
}

namespace
{
    template<class T_SRC, class T_DST> 
    void RotateZHToSH(const Vec3f& dir, int n, const T_SRC* zcoeffs, T_DST* coeffs)
    {
        CL_ASSERT(n >= 2 && n <= 8);

        coeffs[0] = zcoeffs[0];
        
        if (n < 2)
            return;
            
        coeffs++;
        
        float sh1[3] = { dir[1], dir[2], dir[0] };
        coeffs[0] = zcoeffs[1] * sh1[0];
        coeffs[1] = zcoeffs[1] * sh1[1];
        coeffs[2] = zcoeffs[1] * sh1[2];

    // band 3:

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

    // band 4:

        if (n < 4)
            return;

        coeffs += 5;
        float sh3[7];

        sh3[0] =  kSqrt05_06 * (sh1[2] * sh2[0] + sh1[0] * sh2[4]);
        sh3[1] =  kSqrt05_09 * (sh1[1] * sh2[0] + (sh1[2] * sh2[1] + sh1[0] * sh2[3]));
        sh3[2] =  kSqrt08_09 *  sh1[1] * sh2[1] + kSqrt02_03 * sh1[0] * sh2[2] + -kSqrt01_18 * (sh1[2] * sh2[0] - sh1[0] * sh2[4]);
        sh3[3] = -kSqrt01_03 * (sh1[2] * sh2[3] + sh1[0] * sh2[1]) + sh1[1] * sh2[2];
        sh3[4] =  kSqrt08_09 *  sh1[1] * sh2[3] + kSqrt02_03 * sh1[2] * sh2[2] + -kSqrt01_18 * (sh1[2] * sh2[4] + sh1[0] * sh2[0]);
        sh3[5] =  kSqrt05_09 * (sh1[1] * sh2[4] + (sh1[2] * sh2[3] - sh1[0] * sh2[1]));
        sh3[6] =  kSqrt05_06 * (sh1[2] * sh2[4] - sh1[0] * sh2[0]);

        for (int i = 0; i < 7; i++)
            coeffs[i] = zcoeffs[3] * sh3[i];
        
    // band 5:

        if (n < 5)
            return;

        coeffs += 7;
        float sh4[9];

        sh4[0] =  kSqrt07_08 * (sh1[2] * sh3[0] + sh1[0] * sh3[6]);
        sh4[1] =  kSqrt07_16 *  sh1[1] * sh3[0] + kSqrt21_32 * (sh1[2] * sh3[1] + sh1[0] * sh3[5]);
        sh4[2] =  kSqrt03_04 *  sh1[1] * sh3[1] + kSqrt15_32 * (sh1[2] * sh3[2] + sh1[0] * sh3[4]) + -kSqrt01_32 * (sh1[2] * sh3[0] - sh1[0] * sh3[6]);
        sh4[3] =  kSqrt15_16 *  sh1[1] * sh3[2] + kSqrt05_08 * sh1[0] * sh3[3] + -kSqrt03_32 * (sh1[2] * sh3[1] - sh1[0] * sh3[5]);
        sh4[4] = -kSqrt03_08 * (sh1[2] * sh3[4] + sh1[0] * sh3[2]) + sh1[1] * sh3[3];
        sh4[5] =  kSqrt15_16 *  sh1[1] * sh3[4] + kSqrt05_08 * sh1[2] * sh3[3] + -kSqrt03_32 * (sh1[2] * sh3[5] + sh1[0] * sh3[1]);
        sh4[6] =  kSqrt03_04 *  sh1[1] * sh3[5] + kSqrt15_32 * (sh1[2] * sh3[4] - sh1[0] * sh3[2]) + -kSqrt01_32 * (sh1[2] * sh3[6] + sh1[0] * sh3[0]);
        sh4[7] =  kSqrt07_16 *  sh1[1] * sh3[6] + kSqrt21_32 * (sh1[2] * sh3[5] - sh1[0] * sh3[1]);
        sh4[8] =  kSqrt07_08 * (sh1[2] * sh3[6] - sh1[0] * sh3[0]);

        for (int i = 0; i < 9; i++)
            coeffs[i] = zcoeffs[4] * sh4[i];

        // TODO: optimize the bands below (and move/delete this message)

    // band 6:

        if (n < 6)
            return;

        coeffs += 9;
        float sh5[11];

        sh5[ 0] = sqrt(9.0 / 10.0) * (sh1[2] * sh4[0] + sh1[0] * sh4[8]);
        sh5[ 1] = sqrt(9.0 / 25.0) * sh1[1] * sh4[0] + kSqrt18_25 * (sh1[2] * sh4[1] + sh1[0] * sh4[7]);
        sh5[ 2] = sqrt(16.0 / 25.0) * sh1[1] * sh4[1] + kSqrt14_25 * (sh1[2] * sh4[2] + sh1[0] * sh4[6]) + -kSqrt01_50 * (sh1[2] * sh4[0] - sh1[0] * sh4[8]);
        sh5[ 3] = sqrt(21.0 / 25.0) * sh1[1] * sh4[2] + kSqrt21_50 * (sh1[2] * sh4[3] + sh1[0] * sh4[5]) + -kSqrt03_50 * (sh1[2] * sh4[1] - sh1[0] * sh4[7]);
        sh5[ 4] = sqrt(24.0 / 25.0) * sh1[1] * sh4[3] + kSqrt03_05 * sh1[0] * sh4[4] + -kSqrt03_25 * (sh1[2] * sh4[2] - sh1[0] * sh4[6]);
        sh5[ 5] = sh1[1] * sh4[4] + -sqrt(2.0 / 5.0) * (sh1[2] * sh4[5] + sh1[0] * sh4[3]);
        sh5[ 6] = sqrt(24.0 / 25.0) * sh1[1] * sh4[5] + kSqrt03_05 * sh1[2] * sh4[4] + -kSqrt03_25 * (sh1[2] * sh4[6] + sh1[0] * sh4[2]);
        sh5[ 7] = sqrt(21.0 / 25.0) * sh1[1] * sh4[6] + kSqrt21_50 * (sh1[2] * sh4[5] - sh1[0] * sh4[3]) + -kSqrt03_50 * (sh1[2] * sh4[7] + sh1[0] * sh4[1]);
        sh5[ 8] = sqrt(16.0 / 25.0) * sh1[1] * sh4[7] + kSqrt14_25 * (sh1[2] * sh4[6] - sh1[0] * sh4[2]) + -kSqrt01_50 * (sh1[2] * sh4[8] + sh1[0] * sh4[0]);
        sh5[ 9] = sqrt(9.0 / 25.0) * sh1[1] * sh4[8] + kSqrt18_25 * (sh1[2] * sh4[7] - sh1[0] * sh4[1]);
        sh5[10] = sqrt(9.0 / 10.0) * (sh1[2] * sh4[8] - sh1[0] * sh4[0]);

        for (int i = 0; i < 11; i++)
            coeffs[i] = zcoeffs[5] * sh5[i];

    // band 7:

        if (n < 7)
            return;

        coeffs += 11;
        float sh6[13];

        sh6[ 0] = sqrt(11.0 / 12.0) * (sh1[2] * sh5[0] + sh1[0] * sh5[10]);
        sh6[ 1] = sqrt(11.0 / 36.0) * sh1[1] * sh5[0] + sqrt(55.0 / 72.0) * (sh1[2] * sh5[1] + sh1[0] * sh5[9]);
        sh6[ 2] = kSqrt05_09 * sh1[1] * sh5[1] + kSqrt05_08 * (sh1[2] * sh5[2] + sh1[0] * sh5[8]) + -sqrt(1.0 / 72.0) * (sh1[2] * sh5[0] - sh1[0] * sh5[10]);
        sh6[ 3] = kSqrt03_04 * sh1[1] * sh5[2] + sqrt(1.0 / 2.0) * (sh1[2] * sh5[3] + sh1[0] * sh5[7]) + -sqrt(1.0 / 24.0) * (sh1[2] * sh5[1] - sh1[0] * sh5[9]);
        sh6[ 4] = kSqrt08_09 * sh1[1] * sh5[3] + sqrt(7.0 / 18.0) * (sh1[2] * sh5[4] + sh1[0] * sh5[6]) + -kSqrt01_12 * (sh1[2] * sh5[2] - sh1[0] * sh5[8]);
        sh6[ 5] = sqrt(35.0 / 36.0) * sh1[1] * sh5[4] + sqrt(7.0 / 12.0) * sh1[0] * sh5[5] + -sqrt(5.0 / 36.0) * (sh1[2] * sh5[3] - sh1[0] * sh5[7]);
        sh6[ 6] = sh1[1] * sh5[5] + -sqrt(5.0 / 12.0) * (sh1[2] * sh5[6] + sh1[0] * sh5[4]);
        sh6[ 7] = sqrt(35.0 / 36.0) * sh1[1] * sh5[6] + sqrt(7.0 / 12.0) * sh1[2] * sh5[5] + -sqrt(5.0 / 36.0) * (sh1[2] * sh5[7] + sh1[0] * sh5[3]);
        sh6[ 8] = kSqrt08_09 * sh1[1] * sh5[7] + sqrt(7.0 / 18.0) * (sh1[2] * sh5[6] - sh1[0] * sh5[4]) + -kSqrt01_12 * (sh1[2] * sh5[8] + sh1[0] * sh5[2]);
        sh6[ 9] = kSqrt03_04 * sh1[1] * sh5[8] + sqrt(1.0 / 2.0) * (sh1[2] * sh5[7] - sh1[0] * sh5[3]) + -sqrt(1.0 / 24.0) * (sh1[2] * sh5[9] + sh1[0] * sh5[1]);
        sh6[10] = kSqrt05_09 * sh1[1] * sh5[9] + kSqrt05_08 * (sh1[2] * sh5[8] - sh1[0] * sh5[2]) + -sqrt(1.0 / 72.0) * (sh1[2] * sh5[10] + sh1[0] * sh5[0]);
        sh6[11] = sqrt(11.0 / 36.0) * sh1[1] * sh5[10] + sqrt(55.0 / 72.0) * (sh1[2] * sh5[9] - sh1[0] * sh5[1]);
        sh6[12] = sqrt(11.0 / 12.0) * (sh1[2] * sh5[10] - sh1[0] * sh5[0]);

        for (int i = 0; i < 13; i++)
            coeffs[i] = zcoeffs[6] * sh6[i];

    // band 8:

        if (n < 8)
            return;

        coeffs += 13;
        float sh7[15];

        sh7[ 0] = sqrt(13.0 / 14.0) * (sh1[2] * sh6[0] + sh1[0] * sh6[12]);
        sh7[ 1] = sqrt(13.0 / 49.0) * sh1[1] * sh6[0] + sqrt(39.0 / 49.0) * (sh1[2] * sh6[1] + sh1[0] * sh6[11]);
        sh7[ 2] = sqrt(24.0 / 49.0) * sh1[1] * sh6[1] + sqrt(33.0 / 49.0) * (sh1[2] * sh6[2] + sh1[0] * sh6[10]) + -sqrt(1.0 / 98.0) * (sh1[2] * sh6[0] - sh1[0] * sh6[12]);
        sh7[ 3] = sqrt(33.0 / 49.0) * sh1[1] * sh6[2] + sqrt(55.0 / 98.0) * (sh1[2] * sh6[3] + sh1[0] * sh6[9]) + -sqrt(3.0 / 98.0) * (sh1[2] * sh6[1] - sh1[0] * sh6[11]);
        sh7[ 4] = sqrt(40.0 / 49.0) * sh1[1] * sh6[3] + sqrt(45.0 / 98.0) * (sh1[2] * sh6[4] + sh1[0] * sh6[8]) + -sqrt(3.0 / 49.0) * (sh1[2] * sh6[2] - sh1[0] * sh6[10]);
        sh7[ 5] = sqrt(45.0 / 49.0) * sh1[1] * sh6[4] + sqrt(18.0 / 49.0) * (sh1[2] * sh6[5] + sh1[0] * sh6[7]) + -sqrt(5.0 / 49.0) * (sh1[2] * sh6[3] - sh1[0] * sh6[9]);
        sh7[ 6] = sqrt(48.0 / 49.0) * sh1[1] * sh6[5] + sqrt(4.0 / 7.0) * sh1[0] * sh6[6] + -sqrt(15.0 / 98.0) * (sh1[2] * sh6[4] - sh1[0] * sh6[8]);
        sh7[ 7] =                     sh1[1] * sh6[6] + -sqrt(3.0 / 7.0) * (sh1[2] * sh6[7] + sh1[0] * sh6[5]);
        sh7[ 8] = sqrt(48.0 / 49.0) * sh1[1] * sh6[7] + sqrt(4.0 / 7.0) * sh1[2] * sh6[6] + -sqrt(15.0 / 98.0) * (sh1[2] * sh6[8] + sh1[0] * sh6[4]);
        sh7[ 9] = sqrt(45.0 / 49.0) * sh1[1] * sh6[8] + sqrt(18.0 / 49.0) * (sh1[2] * sh6[7] - sh1[0] * sh6[5]) + -sqrt(5.0 / 49.0) * (sh1[2] * sh6[9] + sh1[0] * sh6[3]);
        sh7[10] = sqrt(40.0 / 49.0) * sh1[1] * sh6[9] + sqrt(45.0 / 98.0) * (sh1[2] * sh6[8] - sh1[0] * sh6[4]) + -sqrt(3.0 / 49.0) * (sh1[2] * sh6[10] + sh1[0] * sh6[2]);
        sh7[11] = sqrt(33.0 / 49.0) * sh1[1] * sh6[10] + sqrt(55.0 / 98.0) * (sh1[2] * sh6[9] - sh1[0] * sh6[3]) + -sqrt(3.0 / 98.0) * (sh1[2] * sh6[11] + sh1[0] * sh6[1]);
        sh7[12] = sqrt(24.0 / 49.0) * sh1[1] * sh6[11] + sqrt(33.0 / 49.0) * (sh1[2] * sh6[10] - sh1[0] * sh6[2]) + -sqrt(1.0 / 98.0) * (sh1[2] * sh6[12] + sh1[0] * sh6[0]);
        sh7[13] = sqrt(13.0 / 49.0) * sh1[1] * sh6[12] + sqrt(39.0 / 49.0) * (sh1[2] * sh6[11] - sh1[0] * sh6[1]);
        sh7[14] = sqrt(13.0 / 14.0) * (sh1[2] * sh6[12] - sh1[0] * sh6[0]);

        for (int i = 0; i < 15; i++)
            coeffs[i] = zcoeffs[7] * sh7[i];

    }
}

void SHL::RotateZHToSH(const Vec3f& dir, int n, const float* zcoeffs, float* coeffs)
{
    ::RotateZHToSH<float, float>(dir, n, zcoeffs, coeffs);
}

void SHL::RotateZHToSH(const Vec3f& dir, int n, const Vec4f* zcoeffs, Vec4f* coeffs)
{
    ::RotateZHToSH<Vec4f, Vec4f>(dir, n, zcoeffs, coeffs);
}

namespace
{
    template<class T_SRC, class T_DST>
    void RotateZHToSHAdd(const Vec3f& dir, int n, const T_SRC* zcoeffs, T_DST* coeffs)
    {
        CL_ASSERT(n >= 2 && n <= 8);

        coeffs[0] += zcoeffs[0];
        
        if (n < 2)
            return;
            
        coeffs++;
        
        float sh1[3] = { dir[1], dir[2], dir[0] };
        coeffs[0] += (zcoeffs[1] * sh1[0]);
        coeffs[1] += (zcoeffs[1] * sh1[1]);
        coeffs[2] += (zcoeffs[1] * sh1[2]);

    // band 3:

        if (n < 3)
            return;

        coeffs += 3;
        float sh2[5];

        sh2[0] = kSqrt03_04 * (sh1[2] * sh1[0] + sh1[0] * sh1[2]);
        sh2[1] = kSqrt03_04 * (sh1[1] * sh1[0] + sh1[0] * sh1[1]);
        sh2[2] = sh1[1] * sh1[1] - kSqrt01_04 * (sh1[2] * sh1[2] + sh1[0] * sh1[0]);
        sh2[3] = kSqrt03_04 * (sh1[1] * sh1[2] + sh1[2] * sh1[1]);
        sh2[4] = kSqrt03_04 * (sh1[2] * sh1[2] - sh1[0] * sh1[0]);

        for (int i = 0; i < 5; i++)
            coeffs[i] += (zcoeffs[2] * sh2[i]);

    // band 4:

        if (n < 4)
            return;

        coeffs += 5;
        float sh3[7];

        sh3[0] = kSqrt05_06 * (sh1[2] * sh2[0] + sh1[0] * sh2[4]);
        sh3[1] = kSqrt05_09 * (sh1[1] * sh2[0] + (sh1[2] * sh2[1] + sh1[0] * sh2[3]));
        sh3[2] = kSqrt08_09 * sh1[1] * sh2[1] + kSqrt02_03 * sh1[0] * sh2[2] + -kSqrt01_18 * (sh1[2] * sh2[0] - sh1[0] * sh2[4]);
        sh3[3] = sh1[1] * sh2[2] + -kSqrt01_03 * (sh1[2] * sh2[3] + sh1[0] * sh2[1]);
        sh3[4] = kSqrt08_09 * sh1[1] * sh2[3] + kSqrt02_03 * sh1[2] * sh2[2] + -kSqrt01_18 * (sh1[2] * sh2[4] + sh1[0] * sh2[0]);
        sh3[5] = kSqrt05_09 * (sh1[1] * sh2[4] + (sh1[2] * sh2[3] - sh1[0] * sh2[1]));
        sh3[6] = kSqrt05_06 * (sh1[2] * sh2[4] - sh1[0] * sh2[0]);

        for (int i = 0; i < 7; i++)
            coeffs[i] = (zcoeffs[3] * sh3[i]);

        
    // band 5:

        if (n < 5)
            return;

        coeffs += 7;
        float sh4[9];

        sh4[0] = kSqrt07_08 * (sh1[2] * sh3[0] + sh1[0] * sh3[6]);
        sh4[1] = kSqrt07_16 * sh1[1] * sh3[0] + kSqrt21_32 * (sh1[2] * sh3[1] + sh1[0] * sh3[5]);
        sh4[2] = kSqrt03_04 * sh1[1] * sh3[1] + kSqrt15_32 * (sh1[2] * sh3[2] + sh1[0] * sh3[4]) + -kSqrt01_32 * (sh1[2] * sh3[0] - sh1[0] * sh3[6]);
        sh4[3] = kSqrt15_16 * sh1[1] * sh3[2] + kSqrt05_08 * sh1[0] * sh3[3] + -kSqrt03_32 * (sh1[2] * sh3[1] - sh1[0] * sh3[5]);
        sh4[4] =              sh1[1] * sh3[3] + -kSqrt03_08 * (sh1[2] * sh3[4] + sh1[0] * sh3[2]);
        sh4[5] = kSqrt15_16 * sh1[1] * sh3[4] + kSqrt05_08 * sh1[2] * sh3[3] + -kSqrt03_32 * (sh1[2] * sh3[5] + sh1[0] * sh3[1]);
        sh4[6] = kSqrt03_04 * sh1[1] * sh3[5] + kSqrt15_32 * (sh1[2] * sh3[4] - sh1[0] * sh3[2]) + -kSqrt01_32 * (sh1[2] * sh3[6] + sh1[0] * sh3[0]);
        sh4[7] = kSqrt07_16 * sh1[1] * sh3[6] + kSqrt21_32 * (sh1[2] * sh3[5] - sh1[0] * sh3[1]);
        sh4[8] = kSqrt07_08 * (sh1[2] * sh3[6] - sh1[0] * sh3[0]);

        for (int i = 0; i < 9; i++)
            coeffs[i] = (zcoeffs[4] * sh4[i]);

        // TODO: This code section below needs sqrt replaced (move this message as bands are optimized)

    // band 6:

        if (n < 6)
            return;

        coeffs += 9;
        float sh5[11];

        sh5[ 0] = sqrtf(9.0f / 10.0f) * (sh1[2] * sh4[0] + sh1[0] * sh4[8]);
        sh5[ 1] = sqrt(9.0 / 25.0) * sh1[1] * sh4[0] + kSqrt18_25 * (sh1[2] * sh4[1] + sh1[0] * sh4[7]);
        sh5[ 2] = sqrt(16.0 / 25.0) * sh1[1] * sh4[1] + kSqrt14_25 * (sh1[2] * sh4[2] + sh1[0] * sh4[6]) + -kSqrt01_50 * (sh1[2] * sh4[0] - sh1[0] * sh4[8]);
        sh5[ 3] = sqrt(21.0 / 25.0) * sh1[1] * sh4[2] + kSqrt21_50 * (sh1[2] * sh4[3] + sh1[0] * sh4[5]) + -kSqrt03_50 * (sh1[2] * sh4[1] - sh1[0] * sh4[7]);
        sh5[ 4] = sqrt(24.0 / 25.0) * sh1[1] * sh4[3] + kSqrt03_05 * sh1[0] * sh4[4] + -kSqrt03_25 * (sh1[2] * sh4[2] - sh1[0] * sh4[6]);
        sh5[ 5] = sh1[1] * sh4[4] + -sqrt(2.0 / 5.0) * (sh1[2] * sh4[5] + sh1[0] * sh4[3]);
        sh5[ 6] = sqrt(24.0 / 25.0) * sh1[1] * sh4[5] + kSqrt03_05 * sh1[2] * sh4[4] + -kSqrt03_25 * (sh1[2] * sh4[6] + sh1[0] * sh4[2]);
        sh5[ 7] = sqrt(21.0 / 25.0) * sh1[1] * sh4[6] + kSqrt21_50 * (sh1[2] * sh4[5] - sh1[0] * sh4[3]) + -kSqrt03_50 * (sh1[2] * sh4[7] + sh1[0] * sh4[1]);
        sh5[ 8] = sqrt(16.0 / 25.0) * sh1[1] * sh4[7] + kSqrt14_25 * (sh1[2] * sh4[6] - sh1[0] * sh4[2]) + -kSqrt01_50 * (sh1[2] * sh4[8] + sh1[0] * sh4[0]);
        sh5[ 9] = sqrt(9.0 / 25.0) * sh1[1] * sh4[8] + kSqrt18_25 * (sh1[2] * sh4[7] - sh1[0] * sh4[1]);
        sh5[10] = sqrt(9.0 / 10.0) * (sh1[2] * sh4[8] - sh1[0] * sh4[0]);

        for (int i = 0; i < 11; i++)
            coeffs[i] = (zcoeffs[5] * sh5[i]);

    // band 7:

        if (n < 7)
            return;

        coeffs += 11;
        float sh6[13];

        sh6[ 0] = sqrt(11.0 / 12.0) * (sh1[2] * sh5[0] + sh1[0] * sh5[10]);
        sh6[ 1] = sqrt(11.0 / 36.0) * sh1[1] * sh5[0] + sqrt(55.0 / 72.0) * (sh1[2] * sh5[1] + sh1[0] * sh5[9]);
        sh6[ 2] = kSqrt05_09 * sh1[1] * sh5[1] + kSqrt05_08 * (sh1[2] * sh5[2] + sh1[0] * sh5[8]) + -sqrt(1.0 / 72.0) * (sh1[2] * sh5[0] - sh1[0] * sh5[10]);
        sh6[ 3] = kSqrt03_04 * sh1[1] * sh5[2] + sqrt(1.0 / 2.0) * (sh1[2] * sh5[3] + sh1[0] * sh5[7]) + -sqrt(1.0 / 24.0) * (sh1[2] * sh5[1] - sh1[0] * sh5[9]);
        sh6[ 4] = kSqrt08_09 * sh1[1] * sh5[3] + sqrt(7.0 / 18.0) * (sh1[2] * sh5[4] + sh1[0] * sh5[6]) + -kSqrt01_12 * (sh1[2] * sh5[2] - sh1[0] * sh5[8]);
        sh6[ 5] = sqrt(35.0 / 36.0) * sh1[1] * sh5[4] + sqrt(7.0 / 12.0) * sh1[0] * sh5[5] + -sqrt(5.0 / 36.0) * (sh1[2] * sh5[3] - sh1[0] * sh5[7]);
        sh6[ 6] = sh1[1] * sh5[5] + -sqrt(5.0 / 12.0) * (sh1[2] * sh5[6] + sh1[0] * sh5[4]);
        sh6[ 7] = sqrt(35.0 / 36.0) * sh1[1] * sh5[6] + sqrt(7.0 / 12.0) * sh1[2] * sh5[5] + -sqrt(5.0 / 36.0) * (sh1[2] * sh5[7] + sh1[0] * sh5[3]);
        sh6[ 8] = kSqrt08_09 * sh1[1] * sh5[7] + sqrt(7.0 / 18.0) * (sh1[2] * sh5[6] - sh1[0] * sh5[4]) + -kSqrt01_12 * (sh1[2] * sh5[8] + sh1[0] * sh5[2]);
        sh6[ 9] = kSqrt03_04 * sh1[1] * sh5[8] + sqrt(1.0 / 2.0) * (sh1[2] * sh5[7] - sh1[0] * sh5[3]) + -sqrt(1.0 / 24.0) * (sh1[2] * sh5[9] + sh1[0] * sh5[1]);
        sh6[10] = kSqrt05_09 * sh1[1] * sh5[9] + kSqrt05_08 * (sh1[2] * sh5[8] - sh1[0] * sh5[2]) + -sqrt(1.0 / 72.0) * (sh1[2] * sh5[10] + sh1[0] * sh5[0]);
        sh6[11] = sqrt(11.0 / 36.0) * sh1[1] * sh5[10] + sqrt(55.0 / 72.0) * (sh1[2] * sh5[9] - sh1[0] * sh5[1]);
        sh6[12] = sqrt(11.0 / 12.0) * (sh1[2] * sh5[10] - sh1[0] * sh5[0]);

        for (int i = 0; i < 13; i++)
            coeffs[i] = (zcoeffs[6] * sh6[i]);

    // band 8:

        if (n < 8)
            return;

        coeffs += 13;
        float sh7[15];

        sh7[ 0] = sqrt(13.0 / 14.0) * (sh1[2] * sh6[0] + sh1[0] * sh6[12]);
        sh7[ 1] = sqrt(13.0 / 49.0) * sh1[1] * sh6[0] + sqrt(39.0 / 49.0) * (sh1[2] * sh6[1] + sh1[0] * sh6[11]);
        sh7[ 2] = sqrt(24.0 / 49.0) * sh1[1] * sh6[1] + sqrt(33.0 / 49.0) * (sh1[2] * sh6[2] + sh1[0] * sh6[10]) + -sqrt(1.0 / 98.0) * (sh1[2] * sh6[0] - sh1[0] * sh6[12]);
        sh7[ 3] = sqrt(33.0 / 49.0) * sh1[1] * sh6[2] + sqrt(55.0 / 98.0) * (sh1[2] * sh6[3] + sh1[0] * sh6[9]) + -sqrt(3.0 / 98.0) * (sh1[2] * sh6[1] - sh1[0] * sh6[11]);
        sh7[ 4] = sqrt(40.0 / 49.0) * sh1[1] * sh6[3] + sqrt(45.0 / 98.0) * (sh1[2] * sh6[4] + sh1[0] * sh6[8]) + -sqrt(3.0 / 49.0) * (sh1[2] * sh6[2] - sh1[0] * sh6[10]);
        sh7[ 5] = sqrt(45.0 / 49.0) * sh1[1] * sh6[4] + sqrt(18.0 / 49.0) * (sh1[2] * sh6[5] + sh1[0] * sh6[7]) + -sqrt(5.0 / 49.0) * (sh1[2] * sh6[3] - sh1[0] * sh6[9]);
        sh7[ 6] = sqrt(48.0 / 49.0) * sh1[1] * sh6[5] + sqrt(4.0 / 7.0) * sh1[0] * sh6[6] + -sqrt(15.0 / 98.0) * (sh1[2] * sh6[4] - sh1[0] * sh6[8]);
        sh7[ 7] = sh1[1] * sh6[6] + -sqrt(3.0 / 7.0) * (sh1[2] * sh6[7] + sh1[0] * sh6[5]);
        sh7[ 8] = sqrt(48.0 / 49.0) * sh1[1] * sh6[7] + sqrt(4.0 / 7.0) * sh1[2] * sh6[6] + -sqrt(15.0 / 98.0) * (sh1[2] * sh6[8] + sh1[0] * sh6[4]);
        sh7[ 9] = sqrt(45.0 / 49.0) * sh1[1] * sh6[8] + sqrt(18.0 / 49.0) * (sh1[2] * sh6[7] - sh1[0] * sh6[5]) + -sqrt(5.0 / 49.0) * (sh1[2] * sh6[9] + sh1[0] * sh6[3]);
        sh7[10] = sqrt(40.0 / 49.0) * sh1[1] * sh6[9] + sqrt(45.0 / 98.0) * (sh1[2] * sh6[8] - sh1[0] * sh6[4]) + -sqrt(3.0 / 49.0) * (sh1[2] * sh6[10] + sh1[0] * sh6[2]);
        sh7[11] = sqrt(33.0 / 49.0) * sh1[1] * sh6[10] + sqrt(55.0 / 98.0) * (sh1[2] * sh6[9] - sh1[0] * sh6[3]) + -sqrt(3.0 / 98.0) * (sh1[2] * sh6[11] + sh1[0] * sh6[1]);
        sh7[12] = sqrt(24.0 / 49.0) * sh1[1] * sh6[11] + sqrt(33.0 / 49.0) * (sh1[2] * sh6[10] - sh1[0] * sh6[2]) + -sqrt(1.0 / 98.0) * (sh1[2] * sh6[12] + sh1[0] * sh6[0]);
        sh7[13] = sqrt(13.0 / 49.0) * sh1[1] * sh6[12] + sqrt(39.0 / 49.0) * (sh1[2] * sh6[11] - sh1[0] * sh6[1]);
        sh7[14] = sqrt(13.0 / 14.0) * (sh1[2] * sh6[12] - sh1[0] * sh6[0]);

        for (int i = 0; i < 15; i++)
            coeffs[i] = (zcoeffs[7] * sh7[i]);
    }
}

void SHL::RotateZHToSHAdd(const Vec3f& dir, int n, const float* zcoeffs, float* coeffs)
{
    ::RotateZHToSHAdd<float, float>(dir, n, zcoeffs, coeffs);
}

void SHL::RotateZHToSHAdd(const Vec3f& dir, int n, const Vec4f* zcoeffs, Vec4f* coeffs)
{
    ::RotateZHToSHAdd<Vec4f, Vec4f>(dir, n, zcoeffs, coeffs);
}

void SHL::RotateZHToSHAdd(const Vec3f& dir, int n, const Vec4f& c, const float* zcoeffs, Vec4f* coeffs)
{
    if (n > 8)
        n = 8;

    float scalarCoeffs[64];

    RotateZHToSH(dir, n, zcoeffs, scalarCoeffs);

    int n2 = n * n;

    for (int i = 0; i < n2; i++)
        coeffs[i] += c * scalarCoeffs[i];
}

void SHL::RotateSHAboutZ(float theta, int n, const float* coeffsIn, float* coeffs)
{
    // generic rotation about Z
    coeffs[0] = coeffsIn[0];
        
    if (n < 2)
        return;
    
    coeffs += 1;
    coeffsIn += 1;
    
    float ct = cosf(theta);
    float st = sinf(theta);
    
    coeffs[0] = coeffsIn[0] * ct + coeffsIn[2] * st;
    coeffs[1] = coeffsIn[1];
    coeffs[2] = coeffsIn[2] * ct - coeffsIn[0] * st;
    
    if (n < 3)
        return;
        
    CL_ASSERT(n <= 129); // 129 bands -> 2.128 coeffs 
    float rotCoeffs[256];
    rotCoeffs[0] = ct;
    rotCoeffs[1] = st;
    int nr = 2;
    
    for (int band = 2; band < n; band++)
    {
        // on to the next band
        coeffs += 2 * band - 1;
        coeffsIn += 2 * band - 1;
        
        // find rot coeffs for outermost band coeffs
        rotCoeffs[nr + 0] = ct * rotCoeffs[nr - 2] - st * rotCoeffs[nr - 1];
        rotCoeffs[nr + 1] = ct * rotCoeffs[nr - 1] + st * rotCoeffs[nr - 2];
        nr += 2;

        // rotate outer bands
        for (int i = 0; i < band; i++)
        {
            int j = 2 * band - i;
            int k = 2 * (band - 1) - 2 * i;
            
            coeffs[i] = coeffsIn[i] * rotCoeffs[k] + coeffsIn[j] * rotCoeffs[k + 1];
            coeffs[j] = coeffsIn[j] * rotCoeffs[k] - coeffsIn[i] * rotCoeffs[k + 1];
        }

        // centre coeff doesn't change
        coeffs[band] = coeffsIn[band]; 
    }
}

void SHL::RotateSHAboutZ(float theta, int n, float* coeffs)
{
    // generic in-place rotation about Z
    if (n < 2)
        return;
    
    coeffs++;
    float ct = cosf(theta);
    float st = sinf(theta);
    
    float tmp = coeffs[0] * ct + coeffs[2] * st;
    coeffs[2] = coeffs[2] * ct - coeffs[0] * st;
    coeffs[0] = tmp;
    
    if (n < 3)
        return;
        
    CL_ASSERT(n <= 129); // 129 bands -> 2.128 coeffs 
    float rotCoeffs[256];
    rotCoeffs[0] = ct;
    rotCoeffs[1] = st;
    int nr = 2;
    
    for (int band = 2; band < n; band++)
    {
        // on to the next band
        coeffs += 2 * band - 1;
        
        // find rot coeffs for outermost band coeffs
        rotCoeffs[nr + 0] = ct * rotCoeffs[nr - 2] - st * rotCoeffs[nr - 1];
        rotCoeffs[nr + 1] = ct * rotCoeffs[nr - 1] + st * rotCoeffs[nr - 2];
        nr += 2;

        // rotate outer bands
        for (int i = 0; i < band; i++)
        {
            int j = 2 * band - i;
            int k = 2 * (band - 1) - 2 * i;
            
            float tmp = coeffs[i] * rotCoeffs[k] + coeffs[j] * rotCoeffs[k + 1];
            coeffs[j] = coeffs[j] * rotCoeffs[k] - coeffs[i] * rotCoeffs[k + 1];
            coeffs[i] = tmp;
        }
    }
}

void SHL::RotateSHAboutZ(float theta, int n, Vec4f* coeffs)
{
    // generic in-place rotation about Z
    if (n < 2)
        return;
    
    coeffs++;
    float ct = cosf(theta);
    float st = sinf(theta);
    
    Vec4f tmp = coeffs[0] * ct + coeffs[2] * st;
    coeffs[2] = coeffs[2] * ct - coeffs[0] * st;
    coeffs[0] = tmp;
    
    if (n < 3)
        return;
        
    CL_ASSERT(n <= 129); // 129 bands -> 2.128 coeffs 
    float rotCoeffs[256];
    rotCoeffs[0] = ct;
    rotCoeffs[1] = st;
    int nr = 2;
    
    for (int band = 2; band < n; band++)
    {
        // on to the next band
        coeffs += 2 * band - 1;
        
        // find rot coeffs for outermost band coeffs
        rotCoeffs[nr + 0] = ct * rotCoeffs[nr - 2] - st * rotCoeffs[nr - 1];
        rotCoeffs[nr + 1] = ct * rotCoeffs[nr - 1] + st * rotCoeffs[nr - 2];
        nr += 2;

        // rotate outer bands
        for (int i = 0; i < band; i++)
        {
            int j = 2 * band - i;
            int k = 2 * (band - 1) - 2 * i;
            
            Vec4f tmp = coeffs[i] * rotCoeffs[k] + coeffs[j] * rotCoeffs[k + 1];
            coeffs[j] = coeffs[j] * rotCoeffs[k] - coeffs[i] * rotCoeffs[k + 1];
            coeffs[i] = tmp;
        }
    }
}

namespace
{
    template<class T> void RotateSH(const Mat3f& orient, int n, const T* coeffsIn, T* coeffs)
    /// We could make faster "specialized" versions, as for the last band we don't
    /// need the full sh matrix, just a row.
    {
        CL_ASSERT(n >= 1 && n <= 8);

        (*coeffs++) = coeffsIn[0];
        
        if (n < 2)
            return;
            
        coeffsIn += 1;
        
    #ifdef VL_ROW_ORIENT
        // for row vectors, v' = v M
        float sh1[3][3] =
        {
            orient(1, 1), orient(2, 1), orient(0, 1),
            orient(1, 2), orient(2, 2), orient(0, 2),
            orient(1, 0), orient(2, 0), orient(0, 0)
        };
    #else
        // for column vectors, v' = M v
        float sh1[3][3] =
        {
            orient(1, 1), orient(1, 2), orient(1, 0),
            orient(2, 1), orient(2, 2), orient(2, 0),
            orient(0, 1), orient(0, 2), orient(0, 0)
        };
    #endif

        (*coeffs++) = dp3(coeffsIn, sh1[0]);
        (*coeffs++) = dp3(coeffsIn, sh1[1]);
        (*coeffs++) = dp3(coeffsIn, sh1[2]);
        
    // band 3:

        if (n < 3)
            return;

        coeffsIn += 3;
        float sh2[5][5];

        sh2[0][0] = kSqrt01_04 * ((sh1[2][2] * sh1[0][0] + sh1[2][0] * sh1[0][2]) + (sh1[0][2] * sh1[2][0] + sh1[0][0] * sh1[2][2]));
        sh2[0][1] = (sh1[2][1] * sh1[0][0] + sh1[0][1] * sh1[2][0]);
        sh2[0][2] = kSqrt03_04 * (sh1[2][1] * sh1[0][1] + sh1[0][1] * sh1[2][1]);
        sh2[0][3] = (sh1[2][1] * sh1[0][2] + sh1[0][1] * sh1[2][2]);
        sh2[0][4] = kSqrt01_04 * ((sh1[2][2] * sh1[0][2] - sh1[2][0] * sh1[0][0]) + (sh1[0][2] * sh1[2][2] - sh1[0][0] * sh1[2][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[0]);

        sh2[1][0] = kSqrt01_04 * ((sh1[1][2] * sh1[0][0] + sh1[1][0] * sh1[0][2]) + (sh1[0][2] * sh1[1][0] + sh1[0][0] * sh1[1][2]));
        sh2[1][1] = sh1[1][1] * sh1[0][0] + sh1[0][1] * sh1[1][0];
        sh2[1][2] = kSqrt03_04 * (sh1[1][1] * sh1[0][1] + sh1[0][1] * sh1[1][1]);
        sh2[1][3] = sh1[1][1] * sh1[0][2] + sh1[0][1] * sh1[1][2];
        sh2[1][4] = kSqrt01_04 * ((sh1[1][2] * sh1[0][2] - sh1[1][0] * sh1[0][0]) + (sh1[0][2] * sh1[1][2] - sh1[0][0] * sh1[1][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[1]);

        sh2[2][0] = kSqrt01_03 * (sh1[1][2] * sh1[1][0] + sh1[1][0] * sh1[1][2]) + -kSqrt01_12 * ((sh1[2][2] * sh1[2][0] + sh1[2][0] * sh1[2][2]) + (sh1[0][2] * sh1[0][0] + sh1[0][0] * sh1[0][2]));
        sh2[2][1] = kSqrt04_03 * sh1[1][1] * sh1[1][0] + -kSqrt01_03 * (sh1[2][1] * sh1[2][0] + sh1[0][1] * sh1[0][0]);
        sh2[2][2] = sh1[1][1] * sh1[1][1] + -kSqrt01_04 * (sh1[2][1] * sh1[2][1] + sh1[0][1] * sh1[0][1]);
        sh2[2][3] = kSqrt04_03 * sh1[1][1] * sh1[1][2] + -kSqrt01_03 * (sh1[2][1] * sh1[2][2] + sh1[0][1] * sh1[0][2]);
        sh2[2][4] = kSqrt01_03 * (sh1[1][2] * sh1[1][2] - sh1[1][0] * sh1[1][0]) + -kSqrt01_12 * ((sh1[2][2] * sh1[2][2] - sh1[2][0] * sh1[2][0]) + (sh1[0][2] * sh1[0][2] - sh1[0][0] * sh1[0][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[2]);

        sh2[3][0] = kSqrt01_04 * ((sh1[1][2] * sh1[2][0] + sh1[1][0] * sh1[2][2]) + (sh1[2][2] * sh1[1][0] + sh1[2][0] * sh1[1][2]));
        sh2[3][1] = sh1[1][1] * sh1[2][0] + sh1[2][1] * sh1[1][0];
        sh2[3][2] = kSqrt03_04 * (sh1[1][1] * sh1[2][1] + sh1[2][1] * sh1[1][1]);
        sh2[3][3] = sh1[1][1] * sh1[2][2] + sh1[2][1] * sh1[1][2];
        sh2[3][4] = kSqrt01_04 * ((sh1[1][2] * sh1[2][2] - sh1[1][0] * sh1[2][0]) + (sh1[2][2] * sh1[1][2] - sh1[2][0] * sh1[1][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[3]);

        sh2[4][0] = kSqrt01_04 * ((sh1[2][2] * sh1[2][0] + sh1[2][0] * sh1[2][2]) - (sh1[0][2] * sh1[0][0] + sh1[0][0] * sh1[0][2]));
        sh2[4][1] = (sh1[2][1] * sh1[2][0] - sh1[0][1] * sh1[0][0]);
        sh2[4][2] = kSqrt03_04 * (sh1[2][1] * sh1[2][1] - sh1[0][1] * sh1[0][1]);
        sh2[4][3] = (sh1[2][1] * sh1[2][2] - sh1[0][1] * sh1[0][2]);
        sh2[4][4] = kSqrt01_04 * ((sh1[2][2] * sh1[2][2] - sh1[2][0] * sh1[2][0]) - (sh1[0][2] * sh1[0][2] - sh1[0][0] * sh1[0][0]));

        (*coeffs++) = dp(5, coeffsIn, sh2[4]);

    // band 4:

        if (n < 4)
            return;

        coeffsIn += 5;
        float sh3[7][7];

        sh3[0][0] = kSqrt01_04 * ((sh1[2][2] * sh2[0][0] + sh1[2][0] * sh2[0][4]) + (sh1[0][2] * sh2[4][0] + sh1[0][0] * sh2[4][4]));
        sh3[0][1] = kSqrt03_02 * (sh1[2][1] * sh2[0][0] + sh1[0][1] * sh2[4][0]);
        sh3[0][2] = kSqrt15_16 * (sh1[2][1] * sh2[0][1] + sh1[0][1] * sh2[4][1]);
        sh3[0][3] = kSqrt05_06 * (sh1[2][1] * sh2[0][2] + sh1[0][1] * sh2[4][2]);
        sh3[0][4] = kSqrt15_16 * (sh1[2][1] * sh2[0][3] + sh1[0][1] * sh2[4][3]);
        sh3[0][5] = kSqrt03_02 * (sh1[2][1] * sh2[0][4] + sh1[0][1] * sh2[4][4]);
        sh3[0][6] = kSqrt01_04 * ((sh1[2][2] * sh2[0][4] - sh1[2][0] * sh2[0][0]) + (sh1[0][2] * sh2[4][4] - sh1[0][0] * sh2[4][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[0]);

        sh3[1][0] = kSqrt01_06 * (sh1[1][2] * sh2[0][0] + sh1[1][0] * sh2[0][4]) + kSqrt01_06 * ((sh1[2][2] * sh2[1][0] + sh1[2][0] * sh2[1][4]) + (sh1[0][2] * sh2[3][0] + sh1[0][0] * sh2[3][4]));
        sh3[1][1] = sh1[1][1] * sh2[0][0] + (sh1[2][1] * sh2[1][0] + sh1[0][1] * sh2[3][0]);
        sh3[1][2] = kSqrt05_08 * sh1[1][1] * sh2[0][1] + kSqrt05_08 * (sh1[2][1] * sh2[1][1] + sh1[0][1] * sh2[3][1]);
        sh3[1][3] = kSqrt05_09 * sh1[1][1] * sh2[0][2] + kSqrt05_09 * (sh1[2][1] * sh2[1][2] + sh1[0][1] * sh2[3][2]);
        sh3[1][4] = kSqrt05_08 * sh1[1][1] * sh2[0][3] + kSqrt05_08 * (sh1[2][1] * sh2[1][3] + sh1[0][1] * sh2[3][3]);
        sh3[1][5] = sh1[1][1] * sh2[0][4] + (sh1[2][1] * sh2[1][4] + sh1[0][1] * sh2[3][4]);
        sh3[1][6] = kSqrt01_06 * (sh1[1][2] * sh2[0][4] - sh1[1][0] * sh2[0][0]) + kSqrt01_06 * ((sh1[2][2] * sh2[1][4] - sh1[2][0] * sh2[1][0]) + (sh1[0][2] * sh2[3][4] - sh1[0][0] * sh2[3][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[1]);

        sh3[2][0] = kSqrt04_15 * (sh1[1][2] * sh2[1][0] + sh1[1][0] * sh2[1][4]) + kSqrt01_05 * (sh1[0][2] * sh2[2][0] + sh1[0][0] * sh2[2][4]) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh2[0][0] + sh1[2][0] * sh2[0][4]) - (sh1[0][2] * sh2[4][0] + sh1[0][0] * sh2[4][4]));
        sh3[2][1] = kSqrt08_05 * sh1[1][1] * sh2[1][0] + kSqrt06_05 * sh1[0][1] * sh2[2][0] + -kSqrt01_10 * (sh1[2][1] * sh2[0][0] - sh1[0][1] * sh2[4][0]);
        sh3[2][2] = sh1[1][1] * sh2[1][1] + kSqrt03_04 * sh1[0][1] * sh2[2][1] + -kSqrt01_16 * (sh1[2][1] * sh2[0][1] - sh1[0][1] * sh2[4][1]);
        sh3[2][3] = kSqrt08_09 * sh1[1][1] * sh2[1][2] + kSqrt02_03 * sh1[0][1] * sh2[2][2] + -kSqrt01_18 * (sh1[2][1] * sh2[0][2] - sh1[0][1] * sh2[4][2]);
        sh3[2][4] = sh1[1][1] * sh2[1][3] + kSqrt03_04 * sh1[0][1] * sh2[2][3] + -kSqrt01_16 * (sh1[2][1] * sh2[0][3] - sh1[0][1] * sh2[4][3]);
        sh3[2][5] = kSqrt08_05 * sh1[1][1] * sh2[1][4] + kSqrt06_05 * sh1[0][1] * sh2[2][4] + -kSqrt01_10 * (sh1[2][1] * sh2[0][4] - sh1[0][1] * sh2[4][4]);
        sh3[2][6] = kSqrt04_15 * (sh1[1][2] * sh2[1][4] - sh1[1][0] * sh2[1][0]) + kSqrt01_05 * (sh1[0][2] * sh2[2][4] - sh1[0][0] * sh2[2][0]) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh2[0][4] - sh1[2][0] * sh2[0][0]) - (sh1[0][2] * sh2[4][4] - sh1[0][0] * sh2[4][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[2]);

        sh3[3][0] = kSqrt03_10 * (sh1[1][2] * sh2[2][0] + sh1[1][0] * sh2[2][4]) + -kSqrt01_10 * ((sh1[2][2] * sh2[3][0] + sh1[2][0] * sh2[3][4]) + (sh1[0][2] * sh2[1][0] + sh1[0][0] * sh2[1][4]));
        sh3[3][1] = kSqrt09_05 * sh1[1][1] * sh2[2][0] + -kSqrt03_05 * (sh1[2][1] * sh2[3][0] + sh1[0][1] * sh2[1][0]);
        sh3[3][2] = kSqrt09_08 * sh1[1][1] * sh2[2][1] + -kSqrt03_08 * (sh1[2][1] * sh2[3][1] + sh1[0][1] * sh2[1][1]);
        sh3[3][3] = sh1[1][1] * sh2[2][2] + -kSqrt01_03 * (sh1[2][1] * sh2[3][2] + sh1[0][1] * sh2[1][2]);
        sh3[3][4] = kSqrt09_08 * sh1[1][1] * sh2[2][3] + -kSqrt03_08 * (sh1[2][1] * sh2[3][3] + sh1[0][1] * sh2[1][3]);
        sh3[3][5] = kSqrt09_05 * sh1[1][1] * sh2[2][4] + -kSqrt03_05 * (sh1[2][1] * sh2[3][4] + sh1[0][1] * sh2[1][4]);
        sh3[3][6] = kSqrt03_10 * (sh1[1][2] * sh2[2][4] - sh1[1][0] * sh2[2][0]) + -kSqrt01_10 * ((sh1[2][2] * sh2[3][4] - sh1[2][0] * sh2[3][0]) + (sh1[0][2] * sh2[1][4] - sh1[0][0] * sh2[1][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[3]);

        sh3[4][0] = kSqrt04_15 * (sh1[1][2] * sh2[3][0] + sh1[1][0] * sh2[3][4]) + kSqrt01_05 * (sh1[2][2] * sh2[2][0] + sh1[2][0] * sh2[2][4]) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh2[4][0] + sh1[2][0] * sh2[4][4]) + (sh1[0][2] * sh2[0][0] + sh1[0][0] * sh2[0][4]));
        sh3[4][1] = kSqrt08_05 * sh1[1][1] * sh2[3][0] + kSqrt06_05 * sh1[2][1] * sh2[2][0] + -kSqrt01_10 * (sh1[2][1] * sh2[4][0] + sh1[0][1] * sh2[0][0]);
        sh3[4][2] = sh1[1][1] * sh2[3][1] + kSqrt03_04 * sh1[2][1] * sh2[2][1] + -kSqrt01_16 * (sh1[2][1] * sh2[4][1] + sh1[0][1] * sh2[0][1]);
        sh3[4][3] = kSqrt08_09 * sh1[1][1] * sh2[3][2] + kSqrt02_03 * sh1[2][1] * sh2[2][2] + -kSqrt01_18 * (sh1[2][1] * sh2[4][2] + sh1[0][1] * sh2[0][2]);
        sh3[4][4] = sh1[1][1] * sh2[3][3] + kSqrt03_04 * sh1[2][1] * sh2[2][3] + -kSqrt01_16 * (sh1[2][1] * sh2[4][3] + sh1[0][1] * sh2[0][3]);
        sh3[4][5] = kSqrt08_05 * sh1[1][1] * sh2[3][4] + kSqrt06_05 * sh1[2][1] * sh2[2][4] + -kSqrt01_10 * (sh1[2][1] * sh2[4][4] + sh1[0][1] * sh2[0][4]);
        sh3[4][6] = kSqrt04_15 * (sh1[1][2] * sh2[3][4] - sh1[1][0] * sh2[3][0]) + kSqrt01_05 * (sh1[2][2] * sh2[2][4] - sh1[2][0] * sh2[2][0]) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh2[4][4] - sh1[2][0] * sh2[4][0]) + (sh1[0][2] * sh2[0][4] - sh1[0][0] * sh2[0][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[4]);

        sh3[5][0] = kSqrt01_06 * (sh1[1][2] * sh2[4][0] + sh1[1][0] * sh2[4][4]) + kSqrt01_06 * ((sh1[2][2] * sh2[3][0] + sh1[2][0] * sh2[3][4]) - (sh1[0][2] * sh2[1][0] + sh1[0][0] * sh2[1][4]));
        sh3[5][1] = sh1[1][1] * sh2[4][0] + (sh1[2][1] * sh2[3][0] - sh1[0][1] * sh2[1][0]);
        sh3[5][2] = kSqrt05_08 * sh1[1][1] * sh2[4][1] + kSqrt05_08 * (sh1[2][1] * sh2[3][1] - sh1[0][1] * sh2[1][1]);
        sh3[5][3] = kSqrt05_09 * sh1[1][1] * sh2[4][2] + kSqrt05_09 * (sh1[2][1] * sh2[3][2] - sh1[0][1] * sh2[1][2]);
        sh3[5][4] = kSqrt05_08 * sh1[1][1] * sh2[4][3] + kSqrt05_08 * (sh1[2][1] * sh2[3][3] - sh1[0][1] * sh2[1][3]);
        sh3[5][5] = sh1[1][1] * sh2[4][4] + (sh1[2][1] * sh2[3][4] - sh1[0][1] * sh2[1][4]);
        sh3[5][6] = kSqrt01_06 * (sh1[1][2] * sh2[4][4] - sh1[1][0] * sh2[4][0]) + kSqrt01_06 * ((sh1[2][2] * sh2[3][4] - sh1[2][0] * sh2[3][0]) - (sh1[0][2] * sh2[1][4] - sh1[0][0] * sh2[1][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[5]);

        sh3[6][0] = kSqrt01_04 * ((sh1[2][2] * sh2[4][0] + sh1[2][0] * sh2[4][4]) - (sh1[0][2] * sh2[0][0] + sh1[0][0] * sh2[0][4]));
        sh3[6][1] = kSqrt03_02 * (sh1[2][1] * sh2[4][0] - sh1[0][1] * sh2[0][0]);
        sh3[6][2] = kSqrt15_16 * (sh1[2][1] * sh2[4][1] - sh1[0][1] * sh2[0][1]);
        sh3[6][3] = kSqrt05_06 * (sh1[2][1] * sh2[4][2] - sh1[0][1] * sh2[0][2]);
        sh3[6][4] = kSqrt15_16 * (sh1[2][1] * sh2[4][3] - sh1[0][1] * sh2[0][3]);
        sh3[6][5] = kSqrt03_02 * (sh1[2][1] * sh2[4][4] - sh1[0][1] * sh2[0][4]);
        sh3[6][6] = kSqrt01_04 * ((sh1[2][2] * sh2[4][4] - sh1[2][0] * sh2[4][0]) - (sh1[0][2] * sh2[0][4] - sh1[0][0] * sh2[0][0]));

        (*coeffs++) = dp(7, coeffsIn, sh3[6]);

        // TODO: This code section below needs sqrt replaced

    // band 5:

        if (n < 5)
            return;

        coeffsIn += 7;
        float sh4[9][9];

        sh4[0][0] = kSqrt01_04 * ((sh1[2][2] * sh3[0][0] + sh1[2][0] * sh3[0][6]) + (sh1[0][2] * sh3[6][0] + sh1[0][0] * sh3[6][6]));
        sh4[0][1] = sqrt(2.0 / 1.0) * (sh1[2][1] * sh3[0][0] + sh1[0][1] * sh3[6][0]);
        sh4[0][2] = sqrt(7.0 / 6.0) * (sh1[2][1] * sh3[0][1] + sh1[0][1] * sh3[6][1]);
        sh4[0][3] = sqrt(14.0 / 15.0) * (sh1[2][1] * sh3[0][2] + sh1[0][1] * sh3[6][2]);
        sh4[0][4] = kSqrt07_08 * (sh1[2][1] * sh3[0][3] + sh1[0][1] * sh3[6][3]);
        sh4[0][5] = sqrt(14.0 / 15.0) * (sh1[2][1] * sh3[0][4] + sh1[0][1] * sh3[6][4]);
        sh4[0][6] = sqrt(7.0 / 6.0) * (sh1[2][1] * sh3[0][5] + sh1[0][1] * sh3[6][5]);
        sh4[0][7] = sqrt(2.0 / 1.0) * (sh1[2][1] * sh3[0][6] + sh1[0][1] * sh3[6][6]);
        sh4[0][8] = kSqrt01_04 * ((sh1[2][2] * sh3[0][6] - sh1[2][0] * sh3[0][0]) + (sh1[0][2] * sh3[6][6] - sh1[0][0] * sh3[6][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[0]);

        sh4[1][0] = sqrt(1.0 / 8.0) * (sh1[1][2] * sh3[0][0] + sh1[1][0] * sh3[0][6]) + sqrt(3.0 / 16.0) * ((sh1[2][2] * sh3[1][0] + sh1[2][0] * sh3[1][6]) + (sh1[0][2] * sh3[5][0] + sh1[0][0] * sh3[5][6]));
        sh4[1][1] = sh1[1][1] * sh3[0][0] + kSqrt03_02 * (sh1[2][1] * sh3[1][0] + sh1[0][1] * sh3[5][0]);
        sh4[1][2] = sqrt(7.0 / 12.0) * sh1[1][1] * sh3[0][1] + kSqrt07_08 * (sh1[2][1] * sh3[1][1] + sh1[0][1] * sh3[5][1]);
        sh4[1][3] = sqrt(7.0 / 15.0) * sh1[1][1] * sh3[0][2] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh3[1][2] + sh1[0][1] * sh3[5][2]);
        sh4[1][4] = kSqrt07_16 * sh1[1][1] * sh3[0][3] + kSqrt21_32 * (sh1[2][1] * sh3[1][3] + sh1[0][1] * sh3[5][3]);
        sh4[1][5] = sqrt(7.0 / 15.0) * sh1[1][1] * sh3[0][4] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh3[1][4] + sh1[0][1] * sh3[5][4]);
        sh4[1][6] = sqrt(7.0 / 12.0) * sh1[1][1] * sh3[0][5] + kSqrt07_08 * (sh1[2][1] * sh3[1][5] + sh1[0][1] * sh3[5][5]);
        sh4[1][7] = sh1[1][1] * sh3[0][6] + kSqrt03_02 * (sh1[2][1] * sh3[1][6] + sh1[0][1] * sh3[5][6]);
        sh4[1][8] = sqrt(1.0 / 8.0) * (sh1[1][2] * sh3[0][6] - sh1[1][0] * sh3[0][0]) + sqrt(3.0 / 16.0) * ((sh1[2][2] * sh3[1][6] - sh1[2][0] * sh3[1][0]) + (sh1[0][2] * sh3[5][6] - sh1[0][0] * sh3[5][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[1]);

        sh4[2][0] = sqrt(3.0 / 14.0) * (sh1[1][2] * sh3[1][0] + sh1[1][0] * sh3[1][6]) + sqrt(15.0 / 112.0) * ((sh1[2][2] * sh3[2][0] + sh1[2][0] * sh3[2][6]) + (sh1[0][2] * sh3[4][0] + sh1[0][0] * sh3[4][6])) + -sqrt(1.0 / 112.0) * ((sh1[2][2] * sh3[0][0] + sh1[2][0] * sh3[0][6]) - (sh1[0][2] * sh3[6][0] + sh1[0][0] * sh3[6][6]));
        sh4[2][1] = sqrt(12.0 / 7.0) * sh1[1][1] * sh3[1][0] + sqrt(15.0 / 14.0) * (sh1[2][1] * sh3[2][0] + sh1[0][1] * sh3[4][0]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh3[0][0] - sh1[0][1] * sh3[6][0]);
        sh4[2][2] = sh1[1][1] * sh3[1][1] + kSqrt05_08 * (sh1[2][1] * sh3[2][1] + sh1[0][1] * sh3[4][1]) + -sqrt(1.0 / 24.0) * (sh1[2][1] * sh3[0][1] - sh1[0][1] * sh3[6][1]);
        sh4[2][3] = sqrt(4.0 / 5.0) * sh1[1][1] * sh3[1][2] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh3[2][2] + sh1[0][1] * sh3[4][2]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh3[0][2] - sh1[0][1] * sh3[6][2]);
        sh4[2][4] = kSqrt03_04 * sh1[1][1] * sh3[1][3] + kSqrt15_32 * (sh1[2][1] * sh3[2][3] + sh1[0][1] * sh3[4][3]) + -kSqrt01_32 * (sh1[2][1] * sh3[0][3] - sh1[0][1] * sh3[6][3]);
        sh4[2][5] = sqrt(4.0 / 5.0) * sh1[1][1] * sh3[1][4] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh3[2][4] + sh1[0][1] * sh3[4][4]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh3[0][4] - sh1[0][1] * sh3[6][4]);
        sh4[2][6] = sh1[1][1] * sh3[1][5] + kSqrt05_08 * (sh1[2][1] * sh3[2][5] + sh1[0][1] * sh3[4][5]) + -sqrt(1.0 / 24.0) * (sh1[2][1] * sh3[0][5] - sh1[0][1] * sh3[6][5]);
        sh4[2][7] = sqrt(12.0 / 7.0) * sh1[1][1] * sh3[1][6] + sqrt(15.0 / 14.0) * (sh1[2][1] * sh3[2][6] + sh1[0][1] * sh3[4][6]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh3[0][6] - sh1[0][1] * sh3[6][6]);
        sh4[2][8] = sqrt(3.0 / 14.0) * (sh1[1][2] * sh3[1][6] - sh1[1][0] * sh3[1][0]) + sqrt(15.0 / 112.0) * ((sh1[2][2] * sh3[2][6] - sh1[2][0] * sh3[2][0]) + (sh1[0][2] * sh3[4][6] - sh1[0][0] * sh3[4][0])) + -sqrt(1.0 / 112.0) * ((sh1[2][2] * sh3[0][6] - sh1[2][0] * sh3[0][0]) - (sh1[0][2] * sh3[6][6] - sh1[0][0] * sh3[6][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[2]);

        sh4[3][0] = sqrt(15.0 / 56.0) * (sh1[1][2] * sh3[2][0] + sh1[1][0] * sh3[2][6]) + sqrt(5.0 / 28.0) * (sh1[0][2] * sh3[3][0] + sh1[0][0] * sh3[3][6]) + -sqrt(3.0 / 112.0) * ((sh1[2][2] * sh3[1][0] + sh1[2][0] * sh3[1][6]) - (sh1[0][2] * sh3[5][0] + sh1[0][0] * sh3[5][6]));
        sh4[3][1] = sqrt(15.0 / 7.0) * sh1[1][1] * sh3[2][0] + sqrt(10.0 / 7.0) * sh1[0][1] * sh3[3][0] + -sqrt(3.0 / 14.0) * (sh1[2][1] * sh3[1][0] - sh1[0][1] * sh3[5][0]);
        sh4[3][2] = sqrt(5.0 / 4.0) * sh1[1][1] * sh3[2][1] + kSqrt05_06 * sh1[0][1] * sh3[3][1] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh3[1][1] - sh1[0][1] * sh3[5][1]);
        sh4[3][3] = sh1[1][1] * sh3[2][2] + kSqrt02_03 * sh1[0][1] * sh3[3][2] + -kSqrt01_10 * (sh1[2][1] * sh3[1][2] - sh1[0][1] * sh3[5][2]);
        sh4[3][4] = kSqrt15_16 * sh1[1][1] * sh3[2][3] + kSqrt05_08 * sh1[0][1] * sh3[3][3] + -kSqrt03_32 * (sh1[2][1] * sh3[1][3] - sh1[0][1] * sh3[5][3]);
        sh4[3][5] = sh1[1][1] * sh3[2][4] + kSqrt02_03 * sh1[0][1] * sh3[3][4] + -kSqrt01_10 * (sh1[2][1] * sh3[1][4] - sh1[0][1] * sh3[5][4]);
        sh4[3][6] = sqrt(5.0 / 4.0) * sh1[1][1] * sh3[2][5] + kSqrt05_06 * sh1[0][1] * sh3[3][5] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh3[1][5] - sh1[0][1] * sh3[5][5]);
        sh4[3][7] = sqrt(15.0 / 7.0) * sh1[1][1] * sh3[2][6] + sqrt(10.0 / 7.0) * sh1[0][1] * sh3[3][6] + -sqrt(3.0 / 14.0) * (sh1[2][1] * sh3[1][6] - sh1[0][1] * sh3[5][6]);
        sh4[3][8] = sqrt(15.0 / 56.0) * (sh1[1][2] * sh3[2][6] - sh1[1][0] * sh3[2][0]) + sqrt(5.0 / 28.0) * (sh1[0][2] * sh3[3][6] - sh1[0][0] * sh3[3][0]) + -sqrt(3.0 / 112.0) * ((sh1[2][2] * sh3[1][6] - sh1[2][0] * sh3[1][0]) - (sh1[0][2] * sh3[5][6] - sh1[0][0] * sh3[5][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[3]);

        sh4[4][0] = sqrt(2.0 / 7.0) * (sh1[1][2] * sh3[3][0] + sh1[1][0] * sh3[3][6]) + -sqrt(3.0 / 28.0) * ((sh1[2][2] * sh3[4][0] + sh1[2][0] * sh3[4][6]) + (sh1[0][2] * sh3[2][0] + sh1[0][0] * sh3[2][6]));
        sh4[4][1] = sqrt(16.0 / 7.0) * sh1[1][1] * sh3[3][0] + -sqrt(6.0 / 7.0) * (sh1[2][1] * sh3[4][0] + sh1[0][1] * sh3[2][0]);
        sh4[4][2] = kSqrt04_03 * sh1[1][1] * sh3[3][1] + -sqrt(1.0 / 2.0) * (sh1[2][1] * sh3[4][1] + sh1[0][1] * sh3[2][1]);
        sh4[4][3] = sqrt(16.0 / 15.0) * sh1[1][1] * sh3[3][2] + -sqrt(2.0 / 5.0) * (sh1[2][1] * sh3[4][2] + sh1[0][1] * sh3[2][2]);
        sh4[4][4] = sh1[1][1] * sh3[3][3] + -kSqrt03_08 * (sh1[2][1] * sh3[4][3] + sh1[0][1] * sh3[2][3]);
        sh4[4][5] = sqrt(16.0 / 15.0) * sh1[1][1] * sh3[3][4] + -sqrt(2.0 / 5.0) * (sh1[2][1] * sh3[4][4] + sh1[0][1] * sh3[2][4]);
        sh4[4][6] = kSqrt04_03 * sh1[1][1] * sh3[3][5] + -sqrt(1.0 / 2.0) * (sh1[2][1] * sh3[4][5] + sh1[0][1] * sh3[2][5]);
        sh4[4][7] = sqrt(16.0 / 7.0) * sh1[1][1] * sh3[3][6] + -sqrt(6.0 / 7.0) * (sh1[2][1] * sh3[4][6] + sh1[0][1] * sh3[2][6]);
        sh4[4][8] = sqrt(2.0 / 7.0) * (sh1[1][2] * sh3[3][6] - sh1[1][0] * sh3[3][0]) + -sqrt(3.0 / 28.0) * ((sh1[2][2] * sh3[4][6] - sh1[2][0] * sh3[4][0]) + (sh1[0][2] * sh3[2][6] - sh1[0][0] * sh3[2][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[4]);

        sh4[5][0] = sqrt(15.0 / 56.0) * (sh1[1][2] * sh3[4][0] + sh1[1][0] * sh3[4][6]) + sqrt(5.0 / 28.0) * (sh1[2][2] * sh3[3][0] + sh1[2][0] * sh3[3][6]) + -sqrt(3.0 / 112.0) * ((sh1[2][2] * sh3[5][0] + sh1[2][0] * sh3[5][6]) + (sh1[0][2] * sh3[1][0] + sh1[0][0] * sh3[1][6]));
        sh4[5][1] = sqrt(15.0 / 7.0) * sh1[1][1] * sh3[4][0] + sqrt(10.0 / 7.0) * sh1[2][1] * sh3[3][0] + -sqrt(3.0 / 14.0) * (sh1[2][1] * sh3[5][0] + sh1[0][1] * sh3[1][0]);
        sh4[5][2] = sqrt(5.0 / 4.0) * sh1[1][1] * sh3[4][1] + kSqrt05_06 * sh1[2][1] * sh3[3][1] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh3[5][1] + sh1[0][1] * sh3[1][1]);
        sh4[5][3] = sh1[1][1] * sh3[4][2] + kSqrt02_03 * sh1[2][1] * sh3[3][2] + -kSqrt01_10 * (sh1[2][1] * sh3[5][2] + sh1[0][1] * sh3[1][2]);
        sh4[5][4] = kSqrt15_16 * sh1[1][1] * sh3[4][3] + kSqrt05_08 * sh1[2][1] * sh3[3][3] + -kSqrt03_32 * (sh1[2][1] * sh3[5][3] + sh1[0][1] * sh3[1][3]);
        sh4[5][5] = sh1[1][1] * sh3[4][4] + kSqrt02_03 * sh1[2][1] * sh3[3][4] + -kSqrt01_10 * (sh1[2][1] * sh3[5][4] + sh1[0][1] * sh3[1][4]);
        sh4[5][6] = sqrt(5.0 / 4.0) * sh1[1][1] * sh3[4][5] + kSqrt05_06 * sh1[2][1] * sh3[3][5] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh3[5][5] + sh1[0][1] * sh3[1][5]);
        sh4[5][7] = sqrt(15.0 / 7.0) * sh1[1][1] * sh3[4][6] + sqrt(10.0 / 7.0) * sh1[2][1] * sh3[3][6] + -sqrt(3.0 / 14.0) * (sh1[2][1] * sh3[5][6] + sh1[0][1] * sh3[1][6]);
        sh4[5][8] = sqrt(15.0 / 56.0) * (sh1[1][2] * sh3[4][6] - sh1[1][0] * sh3[4][0]) + sqrt(5.0 / 28.0) * (sh1[2][2] * sh3[3][6] - sh1[2][0] * sh3[3][0]) + -sqrt(3.0 / 112.0) * ((sh1[2][2] * sh3[5][6] - sh1[2][0] * sh3[5][0]) + (sh1[0][2] * sh3[1][6] - sh1[0][0] * sh3[1][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[5]);

        sh4[6][0] = sqrt(3.0 / 14.0) * (sh1[1][2] * sh3[5][0] + sh1[1][0] * sh3[5][6]) + sqrt(15.0 / 112.0) * ((sh1[2][2] * sh3[4][0] + sh1[2][0] * sh3[4][6]) - (sh1[0][2] * sh3[2][0] + sh1[0][0] * sh3[2][6])) + -sqrt(1.0 / 112.0) * ((sh1[2][2] * sh3[6][0] + sh1[2][0] * sh3[6][6]) + (sh1[0][2] * sh3[0][0] + sh1[0][0] * sh3[0][6]));
        sh4[6][1] = sqrt(12.0 / 7.0) * sh1[1][1] * sh3[5][0] + sqrt(15.0 / 14.0) * (sh1[2][1] * sh3[4][0] - sh1[0][1] * sh3[2][0]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh3[6][0] + sh1[0][1] * sh3[0][0]);
        sh4[6][2] = sh1[1][1] * sh3[5][1] + kSqrt05_08 * (sh1[2][1] * sh3[4][1] - sh1[0][1] * sh3[2][1]) + -sqrt(1.0 / 24.0) * (sh1[2][1] * sh3[6][1] + sh1[0][1] * sh3[0][1]);
        sh4[6][3] = sqrt(4.0 / 5.0) * sh1[1][1] * sh3[5][2] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh3[4][2] - sh1[0][1] * sh3[2][2]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh3[6][2] + sh1[0][1] * sh3[0][2]);
        sh4[6][4] = kSqrt03_04 * sh1[1][1] * sh3[5][3] + kSqrt15_32 * (sh1[2][1] * sh3[4][3] - sh1[0][1] * sh3[2][3]) + -kSqrt01_32 * (sh1[2][1] * sh3[6][3] + sh1[0][1] * sh3[0][3]);
        sh4[6][5] = sqrt(4.0 / 5.0) * sh1[1][1] * sh3[5][4] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh3[4][4] - sh1[0][1] * sh3[2][4]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh3[6][4] + sh1[0][1] * sh3[0][4]);
        sh4[6][6] = sh1[1][1] * sh3[5][5] + kSqrt05_08 * (sh1[2][1] * sh3[4][5] - sh1[0][1] * sh3[2][5]) + -sqrt(1.0 / 24.0) * (sh1[2][1] * sh3[6][5] + sh1[0][1] * sh3[0][5]);
        sh4[6][7] = sqrt(12.0 / 7.0) * sh1[1][1] * sh3[5][6] + sqrt(15.0 / 14.0) * (sh1[2][1] * sh3[4][6] - sh1[0][1] * sh3[2][6]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh3[6][6] + sh1[0][1] * sh3[0][6]);
        sh4[6][8] = sqrt(3.0 / 14.0) * (sh1[1][2] * sh3[5][6] - sh1[1][0] * sh3[5][0]) + sqrt(15.0 / 112.0) * ((sh1[2][2] * sh3[4][6] - sh1[2][0] * sh3[4][0]) - (sh1[0][2] * sh3[2][6] - sh1[0][0] * sh3[2][0])) + -sqrt(1.0 / 112.0) * ((sh1[2][2] * sh3[6][6] - sh1[2][0] * sh3[6][0]) + (sh1[0][2] * sh3[0][6] - sh1[0][0] * sh3[0][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[6]);

        sh4[7][0] = sqrt(1.0 / 8.0) * (sh1[1][2] * sh3[6][0] + sh1[1][0] * sh3[6][6]) + sqrt(3.0 / 16.0) * ((sh1[2][2] * sh3[5][0] + sh1[2][0] * sh3[5][6]) - (sh1[0][2] * sh3[1][0] + sh1[0][0] * sh3[1][6]));
        sh4[7][1] = sh1[1][1] * sh3[6][0] + kSqrt03_02 * (sh1[2][1] * sh3[5][0] - sh1[0][1] * sh3[1][0]);
        sh4[7][2] = sqrt(7.0 / 12.0) * sh1[1][1] * sh3[6][1] + kSqrt07_08 * (sh1[2][1] * sh3[5][1] - sh1[0][1] * sh3[1][1]);
        sh4[7][3] = sqrt(7.0 / 15.0) * sh1[1][1] * sh3[6][2] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh3[5][2] - sh1[0][1] * sh3[1][2]);
        sh4[7][4] = kSqrt07_16 * sh1[1][1] * sh3[6][3] + kSqrt21_32 * (sh1[2][1] * sh3[5][3] - sh1[0][1] * sh3[1][3]);
        sh4[7][5] = sqrt(7.0 / 15.0) * sh1[1][1] * sh3[6][4] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh3[5][4] - sh1[0][1] * sh3[1][4]);
        sh4[7][6] = sqrt(7.0 / 12.0) * sh1[1][1] * sh3[6][5] + kSqrt07_08 * (sh1[2][1] * sh3[5][5] - sh1[0][1] * sh3[1][5]);
        sh4[7][7] = sh1[1][1] * sh3[6][6] + kSqrt03_02 * (sh1[2][1] * sh3[5][6] - sh1[0][1] * sh3[1][6]);
        sh4[7][8] = sqrt(1.0 / 8.0) * (sh1[1][2] * sh3[6][6] - sh1[1][0] * sh3[6][0]) + sqrt(3.0 / 16.0) * ((sh1[2][2] * sh3[5][6] - sh1[2][0] * sh3[5][0]) - (sh1[0][2] * sh3[1][6] - sh1[0][0] * sh3[1][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[7]);

        sh4[8][0] = kSqrt01_04 * ((sh1[2][2] * sh3[6][0] + sh1[2][0] * sh3[6][6]) - (sh1[0][2] * sh3[0][0] + sh1[0][0] * sh3[0][6]));
        sh4[8][1] = sqrt(2.0 / 1.0) * (sh1[2][1] * sh3[6][0] - sh1[0][1] * sh3[0][0]);
        sh4[8][2] = sqrt(7.0 / 6.0) * (sh1[2][1] * sh3[6][1] - sh1[0][1] * sh3[0][1]);
        sh4[8][3] = sqrt(14.0 / 15.0) * (sh1[2][1] * sh3[6][2] - sh1[0][1] * sh3[0][2]);
        sh4[8][4] = kSqrt07_08 * (sh1[2][1] * sh3[6][3] - sh1[0][1] * sh3[0][3]);
        sh4[8][5] = sqrt(14.0 / 15.0) * (sh1[2][1] * sh3[6][4] - sh1[0][1] * sh3[0][4]);
        sh4[8][6] = sqrt(7.0 / 6.0) * (sh1[2][1] * sh3[6][5] - sh1[0][1] * sh3[0][5]);
        sh4[8][7] = sqrt(2.0 / 1.0) * (sh1[2][1] * sh3[6][6] - sh1[0][1] * sh3[0][6]);
        sh4[8][8] = kSqrt01_04 * ((sh1[2][2] * sh3[6][6] - sh1[2][0] * sh3[6][0]) - (sh1[0][2] * sh3[0][6] - sh1[0][0] * sh3[0][0]));

        (*coeffs++) = dp(9, coeffsIn, sh4[8]);

    // band 6:

        if (n < 6)
            return;

        coeffsIn += 9;
        float sh5[11][11];

        sh5[0][0] = kSqrt01_04 * ((sh1[2][2] * sh4[0][0] + sh1[2][0] * sh4[0][8]) + (sh1[0][2] * sh4[8][0] + sh1[0][0] * sh4[8][8]));
        sh5[0][1] = sqrt(5.0 / 2.0) * (sh1[2][1] * sh4[0][0] + sh1[0][1] * sh4[8][0]);
        sh5[0][2] = sqrt(45.0 / 32.0) * (sh1[2][1] * sh4[0][1] + sh1[0][1] * sh4[8][1]);
        sh5[0][3] = sqrt(15.0 / 14.0) * (sh1[2][1] * sh4[0][2] + sh1[0][1] * sh4[8][2]);
        sh5[0][4] = kSqrt15_16 * (sh1[2][1] * sh4[0][3] + sh1[0][1] * sh4[8][3]);
        sh5[0][5] = sqrt(9.0 / 10.0) * (sh1[2][1] * sh4[0][4] + sh1[0][1] * sh4[8][4]);
        sh5[0][6] = kSqrt15_16 * (sh1[2][1] * sh4[0][5] + sh1[0][1] * sh4[8][5]);
        sh5[0][7] = sqrt(15.0 / 14.0) * (sh1[2][1] * sh4[0][6] + sh1[0][1] * sh4[8][6]);
        sh5[0][8] = sqrt(45.0 / 32.0) * (sh1[2][1] * sh4[0][7] + sh1[0][1] * sh4[8][7]);
        sh5[0][9] = sqrt(5.0 / 2.0) * (sh1[2][1] * sh4[0][8] + sh1[0][1] * sh4[8][8]);
        sh5[0][10] = kSqrt01_04 * ((sh1[2][2] * sh4[0][8] - sh1[2][0] * sh4[0][0]) + (sh1[0][2] * sh4[8][8] - sh1[0][0] * sh4[8][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[0]);

        sh5[1][0] = kSqrt01_10 * (sh1[1][2] * sh4[0][0] + sh1[1][0] * sh4[0][8]) + kSqrt01_05 * ((sh1[2][2] * sh4[1][0] + sh1[2][0] * sh4[1][8]) + (sh1[0][2] * sh4[7][0] + sh1[0][0] * sh4[7][8]));
        sh5[1][1] = sh1[1][1] * sh4[0][0] + sqrt(2.0 / 1.0) * (sh1[2][1] * sh4[1][0] + sh1[0][1] * sh4[7][0]);
        sh5[1][2] = sqrt(9.0 / 16.0) * sh1[1][1] * sh4[0][1] + kSqrt09_08 * (sh1[2][1] * sh4[1][1] + sh1[0][1] * sh4[7][1]);
        sh5[1][3] = sqrt(3.0 / 7.0) * sh1[1][1] * sh4[0][2] + sqrt(6.0 / 7.0) * (sh1[2][1] * sh4[1][2] + sh1[0][1] * sh4[7][2]);
        sh5[1][4] = kSqrt03_08 * sh1[1][1] * sh4[0][3] + kSqrt03_04 * (sh1[2][1] * sh4[1][3] + sh1[0][1] * sh4[7][3]);
        sh5[1][5] = sqrt(9.0 / 25.0) * sh1[1][1] * sh4[0][4] + kSqrt18_25 * (sh1[2][1] * sh4[1][4] + sh1[0][1] * sh4[7][4]);
        sh5[1][6] = kSqrt03_08 * sh1[1][1] * sh4[0][5] + kSqrt03_04 * (sh1[2][1] * sh4[1][5] + sh1[0][1] * sh4[7][5]);
        sh5[1][7] = sqrt(3.0 / 7.0) * sh1[1][1] * sh4[0][6] + sqrt(6.0 / 7.0) * (sh1[2][1] * sh4[1][6] + sh1[0][1] * sh4[7][6]);
        sh5[1][8] = sqrt(9.0 / 16.0) * sh1[1][1] * sh4[0][7] + kSqrt09_08 * (sh1[2][1] * sh4[1][7] + sh1[0][1] * sh4[7][7]);
        sh5[1][9] = sh1[1][1] * sh4[0][8] + sqrt(2.0 / 1.0) * (sh1[2][1] * sh4[1][8] + sh1[0][1] * sh4[7][8]);
        sh5[1][10] = kSqrt01_10 * (sh1[1][2] * sh4[0][8] - sh1[1][0] * sh4[0][0]) + kSqrt01_05 * ((sh1[2][2] * sh4[1][8] - sh1[2][0] * sh4[1][0]) + (sh1[0][2] * sh4[7][8] - sh1[0][0] * sh4[7][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[1]);

        sh5[2][0] = sqrt(8.0 / 45.0) * (sh1[1][2] * sh4[1][0] + sh1[1][0] * sh4[1][8]) + sqrt(7.0 / 45.0) * ((sh1[2][2] * sh4[2][0] + sh1[2][0] * sh4[2][8]) + (sh1[0][2] * sh4[6][0] + sh1[0][0] * sh4[6][8])) + -sqrt(1.0 / 180.0) * ((sh1[2][2] * sh4[0][0] + sh1[2][0] * sh4[0][8]) - (sh1[0][2] * sh4[8][0] + sh1[0][0] * sh4[8][8]));
        sh5[2][1] = sqrt(16.0 / 9.0) * sh1[1][1] * sh4[1][0] + sqrt(14.0 / 9.0) * (sh1[2][1] * sh4[2][0] + sh1[0][1] * sh4[6][0]) + -kSqrt01_18 * (sh1[2][1] * sh4[0][0] - sh1[0][1] * sh4[8][0]);
        sh5[2][2] = sh1[1][1] * sh4[1][1] + kSqrt07_08 * (sh1[2][1] * sh4[2][1] + sh1[0][1] * sh4[6][1]) + -kSqrt01_32 * (sh1[2][1] * sh4[0][1] - sh1[0][1] * sh4[8][1]);
        sh5[2][3] = sqrt(16.0 / 21.0) * sh1[1][1] * sh4[1][2] + kSqrt02_03 * (sh1[2][1] * sh4[2][2] + sh1[0][1] * sh4[6][2]) + -sqrt(1.0 / 42.0) * (sh1[2][1] * sh4[0][2] - sh1[0][1] * sh4[8][2]);
        sh5[2][4] = kSqrt02_03 * sh1[1][1] * sh4[1][3] + sqrt(7.0 / 12.0) * (sh1[2][1] * sh4[2][3] + sh1[0][1] * sh4[6][3]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh4[0][3] - sh1[0][1] * sh4[8][3]);
        sh5[2][5] = sqrt(16.0 / 25.0) * sh1[1][1] * sh4[1][4] + kSqrt14_25 * (sh1[2][1] * sh4[2][4] + sh1[0][1] * sh4[6][4]) + -kSqrt01_50 * (sh1[2][1] * sh4[0][4] - sh1[0][1] * sh4[8][4]);
        sh5[2][6] = kSqrt02_03 * sh1[1][1] * sh4[1][5] + sqrt(7.0 / 12.0) * (sh1[2][1] * sh4[2][5] + sh1[0][1] * sh4[6][5]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh4[0][5] - sh1[0][1] * sh4[8][5]);
        sh5[2][7] = sqrt(16.0 / 21.0) * sh1[1][1] * sh4[1][6] + kSqrt02_03 * (sh1[2][1] * sh4[2][6] + sh1[0][1] * sh4[6][6]) + -sqrt(1.0 / 42.0) * (sh1[2][1] * sh4[0][6] - sh1[0][1] * sh4[8][6]);
        sh5[2][8] = sh1[1][1] * sh4[1][7] + kSqrt07_08 * (sh1[2][1] * sh4[2][7] + sh1[0][1] * sh4[6][7]) + -kSqrt01_32 * (sh1[2][1] * sh4[0][7] - sh1[0][1] * sh4[8][7]);
        sh5[2][9] = sqrt(16.0 / 9.0) * sh1[1][1] * sh4[1][8] + sqrt(14.0 / 9.0) * (sh1[2][1] * sh4[2][8] + sh1[0][1] * sh4[6][8]) + -kSqrt01_18 * (sh1[2][1] * sh4[0][8] - sh1[0][1] * sh4[8][8]);
        sh5[2][10] = sqrt(8.0 / 45.0) * (sh1[1][2] * sh4[1][8] - sh1[1][0] * sh4[1][0]) + sqrt(7.0 / 45.0) * ((sh1[2][2] * sh4[2][8] - sh1[2][0] * sh4[2][0]) + (sh1[0][2] * sh4[6][8] - sh1[0][0] * sh4[6][0])) + -sqrt(1.0 / 180.0) * ((sh1[2][2] * sh4[0][8] - sh1[2][0] * sh4[0][0]) - (sh1[0][2] * sh4[8][8] - sh1[0][0] * sh4[8][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[2]);

        sh5[3][0] = sqrt(7.0 / 30.0) * (sh1[1][2] * sh4[2][0] + sh1[1][0] * sh4[2][8]) + sqrt(7.0 / 60.0) * ((sh1[2][2] * sh4[3][0] + sh1[2][0] * sh4[3][8]) + (sh1[0][2] * sh4[5][0] + sh1[0][0] * sh4[5][8])) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh4[1][0] + sh1[2][0] * sh4[1][8]) - (sh1[0][2] * sh4[7][0] + sh1[0][0] * sh4[7][8]));
        sh5[3][1] = sqrt(7.0 / 3.0) * sh1[1][1] * sh4[2][0] + sqrt(7.0 / 6.0) * (sh1[2][1] * sh4[3][0] + sh1[0][1] * sh4[5][0]) + -kSqrt01_06 * (sh1[2][1] * sh4[1][0] - sh1[0][1] * sh4[7][0]);
        sh5[3][2] = sqrt(21.0 / 16.0) * sh1[1][1] * sh4[2][1] + kSqrt21_32 * (sh1[2][1] * sh4[3][1] + sh1[0][1] * sh4[5][1]) + -kSqrt03_32 * (sh1[2][1] * sh4[1][1] - sh1[0][1] * sh4[7][1]);
        sh5[3][3] = sh1[1][1] * sh4[2][2] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh4[3][2] + sh1[0][1] * sh4[5][2]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh4[1][2] - sh1[0][1] * sh4[7][2]);
        sh5[3][4] = kSqrt07_08 * sh1[1][1] * sh4[2][3] + kSqrt07_16 * (sh1[2][1] * sh4[3][3] + sh1[0][1] * sh4[5][3]) + -kSqrt01_16 * (sh1[2][1] * sh4[1][3] - sh1[0][1] * sh4[7][3]);
        sh5[3][5] = sqrt(21.0 / 25.0) * sh1[1][1] * sh4[2][4] + kSqrt21_50 * (sh1[2][1] * sh4[3][4] + sh1[0][1] * sh4[5][4]) + -kSqrt03_50 * (sh1[2][1] * sh4[1][4] - sh1[0][1] * sh4[7][4]);
        sh5[3][6] = kSqrt07_08 * sh1[1][1] * sh4[2][5] + kSqrt07_16 * (sh1[2][1] * sh4[3][5] + sh1[0][1] * sh4[5][5]) + -kSqrt01_16 * (sh1[2][1] * sh4[1][5] - sh1[0][1] * sh4[7][5]);
        sh5[3][7] = sh1[1][1] * sh4[2][6] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh4[3][6] + sh1[0][1] * sh4[5][6]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh4[1][6] - sh1[0][1] * sh4[7][6]);
        sh5[3][8] = sqrt(21.0 / 16.0) * sh1[1][1] * sh4[2][7] + kSqrt21_32 * (sh1[2][1] * sh4[3][7] + sh1[0][1] * sh4[5][7]) + -kSqrt03_32 * (sh1[2][1] * sh4[1][7] - sh1[0][1] * sh4[7][7]);
        sh5[3][9] = sqrt(7.0 / 3.0) * sh1[1][1] * sh4[2][8] + sqrt(7.0 / 6.0) * (sh1[2][1] * sh4[3][8] + sh1[0][1] * sh4[5][8]) + -kSqrt01_06 * (sh1[2][1] * sh4[1][8] - sh1[0][1] * sh4[7][8]);
        sh5[3][10] = sqrt(7.0 / 30.0) * (sh1[1][2] * sh4[2][8] - sh1[1][0] * sh4[2][0]) + sqrt(7.0 / 60.0) * ((sh1[2][2] * sh4[3][8] - sh1[2][0] * sh4[3][0]) + (sh1[0][2] * sh4[5][8] - sh1[0][0] * sh4[5][0])) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh4[1][8] - sh1[2][0] * sh4[1][0]) - (sh1[0][2] * sh4[7][8] - sh1[0][0] * sh4[7][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[3]);

        sh5[4][0] = kSqrt04_15 * (sh1[1][2] * sh4[3][0] + sh1[1][0] * sh4[3][8]) + kSqrt01_06 * (sh1[0][2] * sh4[4][0] + sh1[0][0] * sh4[4][8]) + -sqrt(1.0 / 30.0) * ((sh1[2][2] * sh4[2][0] + sh1[2][0] * sh4[2][8]) - (sh1[0][2] * sh4[6][0] + sh1[0][0] * sh4[6][8]));
        sh5[4][1] = sqrt(8.0 / 3.0) * sh1[1][1] * sh4[3][0] + sqrt(5.0 / 3.0) * sh1[0][1] * sh4[4][0] + -kSqrt01_03 * (sh1[2][1] * sh4[2][0] - sh1[0][1] * sh4[6][0]);
        sh5[4][2] = kSqrt03_02 * sh1[1][1] * sh4[3][1] + kSqrt15_16 * sh1[0][1] * sh4[4][1] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh4[2][1] - sh1[0][1] * sh4[6][1]);
        sh5[4][3] = sqrt(8.0 / 7.0) * sh1[1][1] * sh4[3][2] + sqrt(5.0 / 7.0) * sh1[0][1] * sh4[4][2] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh4[2][2] - sh1[0][1] * sh4[6][2]);
        sh5[4][4] = sh1[1][1] * sh4[3][3] + kSqrt05_08 * sh1[0][1] * sh4[4][3] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh4[2][3] - sh1[0][1] * sh4[6][3]);
        sh5[4][5] = sqrt(24.0 / 25.0) * sh1[1][1] * sh4[3][4] + kSqrt03_05 * sh1[0][1] * sh4[4][4] + -kSqrt03_25 * (sh1[2][1] * sh4[2][4] - sh1[0][1] * sh4[6][4]);
        sh5[4][6] = sh1[1][1] * sh4[3][5] + kSqrt05_08 * sh1[0][1] * sh4[4][5] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh4[2][5] - sh1[0][1] * sh4[6][5]);
        sh5[4][7] = sqrt(8.0 / 7.0) * sh1[1][1] * sh4[3][6] + sqrt(5.0 / 7.0) * sh1[0][1] * sh4[4][6] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh4[2][6] - sh1[0][1] * sh4[6][6]);
        sh5[4][8] = kSqrt03_02 * sh1[1][1] * sh4[3][7] + kSqrt15_16 * sh1[0][1] * sh4[4][7] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh4[2][7] - sh1[0][1] * sh4[6][7]);
        sh5[4][9] = sqrt(8.0 / 3.0) * sh1[1][1] * sh4[3][8] + sqrt(5.0 / 3.0) * sh1[0][1] * sh4[4][8] + -kSqrt01_03 * (sh1[2][1] * sh4[2][8] - sh1[0][1] * sh4[6][8]);
        sh5[4][10] = kSqrt04_15 * (sh1[1][2] * sh4[3][8] - sh1[1][0] * sh4[3][0]) + kSqrt01_06 * (sh1[0][2] * sh4[4][8] - sh1[0][0] * sh4[4][0]) + -sqrt(1.0 / 30.0) * ((sh1[2][2] * sh4[2][8] - sh1[2][0] * sh4[2][0]) - (sh1[0][2] * sh4[6][8] - sh1[0][0] * sh4[6][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[4]);

        sh5[5][0] = sqrt(5.0 / 18.0) * (sh1[1][2] * sh4[4][0] + sh1[1][0] * sh4[4][8]) + -sqrt(1.0 / 9.0) * ((sh1[2][2] * sh4[5][0] + sh1[2][0] * sh4[5][8]) + (sh1[0][2] * sh4[3][0] + sh1[0][0] * sh4[3][8]));
        sh5[5][1] = sqrt(25.0 / 9.0) * sh1[1][1] * sh4[4][0] + -sqrt(10.0 / 9.0) * (sh1[2][1] * sh4[5][0] + sh1[0][1] * sh4[3][0]);
        sh5[5][2] = sqrt(25.0 / 16.0) * sh1[1][1] * sh4[4][1] + -kSqrt05_08 * (sh1[2][1] * sh4[5][1] + sh1[0][1] * sh4[3][1]);
        sh5[5][3] = sqrt(25.0 / 21.0) * sh1[1][1] * sh4[4][2] + -sqrt(10.0 / 21.0) * (sh1[2][1] * sh4[5][2] + sh1[0][1] * sh4[3][2]);
        sh5[5][4] = sqrt(25.0 / 24.0) * sh1[1][1] * sh4[4][3] + -sqrt(5.0 / 12.0) * (sh1[2][1] * sh4[5][3] + sh1[0][1] * sh4[3][3]);
        sh5[5][5] = sh1[1][1] * sh4[4][4] + -sqrt(2.0 / 5.0) * (sh1[2][1] * sh4[5][4] + sh1[0][1] * sh4[3][4]);
        sh5[5][6] = sqrt(25.0 / 24.0) * sh1[1][1] * sh4[4][5] + -sqrt(5.0 / 12.0) * (sh1[2][1] * sh4[5][5] + sh1[0][1] * sh4[3][5]);
        sh5[5][7] = sqrt(25.0 / 21.0) * sh1[1][1] * sh4[4][6] + -sqrt(10.0 / 21.0) * (sh1[2][1] * sh4[5][6] + sh1[0][1] * sh4[3][6]);
        sh5[5][8] = sqrt(25.0 / 16.0) * sh1[1][1] * sh4[4][7] + -kSqrt05_08 * (sh1[2][1] * sh4[5][7] + sh1[0][1] * sh4[3][7]);
        sh5[5][9] = sqrt(25.0 / 9.0) * sh1[1][1] * sh4[4][8] + -sqrt(10.0 / 9.0) * (sh1[2][1] * sh4[5][8] + sh1[0][1] * sh4[3][8]);
        sh5[5][10] = sqrt(5.0 / 18.0) * (sh1[1][2] * sh4[4][8] - sh1[1][0] * sh4[4][0]) + -sqrt(1.0 / 9.0) * ((sh1[2][2] * sh4[5][8] - sh1[2][0] * sh4[5][0]) + (sh1[0][2] * sh4[3][8] - sh1[0][0] * sh4[3][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[5]);

        sh5[6][0] = kSqrt04_15 * (sh1[1][2] * sh4[5][0] + sh1[1][0] * sh4[5][8]) + kSqrt01_06 * (sh1[2][2] * sh4[4][0] + sh1[2][0] * sh4[4][8]) + -sqrt(1.0 / 30.0) * ((sh1[2][2] * sh4[6][0] + sh1[2][0] * sh4[6][8]) + (sh1[0][2] * sh4[2][0] + sh1[0][0] * sh4[2][8]));
        sh5[6][1] = sqrt(8.0 / 3.0) * sh1[1][1] * sh4[5][0] + sqrt(5.0 / 3.0) * sh1[2][1] * sh4[4][0] + -kSqrt01_03 * (sh1[2][1] * sh4[6][0] + sh1[0][1] * sh4[2][0]);
        sh5[6][2] = kSqrt03_02 * sh1[1][1] * sh4[5][1] + kSqrt15_16 * sh1[2][1] * sh4[4][1] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh4[6][1] + sh1[0][1] * sh4[2][1]);
        sh5[6][3] = sqrt(8.0 / 7.0) * sh1[1][1] * sh4[5][2] + sqrt(5.0 / 7.0) * sh1[2][1] * sh4[4][2] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh4[6][2] + sh1[0][1] * sh4[2][2]);
        sh5[6][4] = sh1[1][1] * sh4[5][3] + kSqrt05_08 * sh1[2][1] * sh4[4][3] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh4[6][3] + sh1[0][1] * sh4[2][3]);
        sh5[6][5] = sqrt(24.0 / 25.0) * sh1[1][1] * sh4[5][4] + kSqrt03_05 * sh1[2][1] * sh4[4][4] + -kSqrt03_25 * (sh1[2][1] * sh4[6][4] + sh1[0][1] * sh4[2][4]);
        sh5[6][6] = sh1[1][1] * sh4[5][5] + kSqrt05_08 * sh1[2][1] * sh4[4][5] + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh4[6][5] + sh1[0][1] * sh4[2][5]);
        sh5[6][7] = sqrt(8.0 / 7.0) * sh1[1][1] * sh4[5][6] + sqrt(5.0 / 7.0) * sh1[2][1] * sh4[4][6] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh4[6][6] + sh1[0][1] * sh4[2][6]);
        sh5[6][8] = kSqrt03_02 * sh1[1][1] * sh4[5][7] + kSqrt15_16 * sh1[2][1] * sh4[4][7] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh4[6][7] + sh1[0][1] * sh4[2][7]);
        sh5[6][9] = sqrt(8.0 / 3.0) * sh1[1][1] * sh4[5][8] + sqrt(5.0 / 3.0) * sh1[2][1] * sh4[4][8] + -kSqrt01_03 * (sh1[2][1] * sh4[6][8] + sh1[0][1] * sh4[2][8]);
        sh5[6][10] = kSqrt04_15 * (sh1[1][2] * sh4[5][8] - sh1[1][0] * sh4[5][0]) + kSqrt01_06 * (sh1[2][2] * sh4[4][8] - sh1[2][0] * sh4[4][0]) + -sqrt(1.0 / 30.0) * ((sh1[2][2] * sh4[6][8] - sh1[2][0] * sh4[6][0]) + (sh1[0][2] * sh4[2][8] - sh1[0][0] * sh4[2][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[6]);

        sh5[7][0] = sqrt(7.0 / 30.0) * (sh1[1][2] * sh4[6][0] + sh1[1][0] * sh4[6][8]) + sqrt(7.0 / 60.0) * ((sh1[2][2] * sh4[5][0] + sh1[2][0] * sh4[5][8]) - (sh1[0][2] * sh4[3][0] + sh1[0][0] * sh4[3][8])) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh4[7][0] + sh1[2][0] * sh4[7][8]) + (sh1[0][2] * sh4[1][0] + sh1[0][0] * sh4[1][8]));
        sh5[7][1] = sqrt(7.0 / 3.0) * sh1[1][1] * sh4[6][0] + sqrt(7.0 / 6.0) * (sh1[2][1] * sh4[5][0] - sh1[0][1] * sh4[3][0]) + -kSqrt01_06 * (sh1[2][1] * sh4[7][0] + sh1[0][1] * sh4[1][0]);
        sh5[7][2] = sqrt(21.0 / 16.0) * sh1[1][1] * sh4[6][1] + kSqrt21_32 * (sh1[2][1] * sh4[5][1] - sh1[0][1] * sh4[3][1]) + -kSqrt03_32 * (sh1[2][1] * sh4[7][1] + sh1[0][1] * sh4[1][1]);
        sh5[7][3] = sh1[1][1] * sh4[6][2] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh4[5][2] - sh1[0][1] * sh4[3][2]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh4[7][2] + sh1[0][1] * sh4[1][2]);
        sh5[7][4] = kSqrt07_08 * sh1[1][1] * sh4[6][3] + kSqrt07_16 * (sh1[2][1] * sh4[5][3] - sh1[0][1] * sh4[3][3]) + -kSqrt01_16 * (sh1[2][1] * sh4[7][3] + sh1[0][1] * sh4[1][3]);
        sh5[7][5] = sqrt(21.0 / 25.0) * sh1[1][1] * sh4[6][4] + kSqrt21_50 * (sh1[2][1] * sh4[5][4] - sh1[0][1] * sh4[3][4]) + -kSqrt03_50 * (sh1[2][1] * sh4[7][4] + sh1[0][1] * sh4[1][4]);
        sh5[7][6] = kSqrt07_08 * sh1[1][1] * sh4[6][5] + kSqrt07_16 * (sh1[2][1] * sh4[5][5] - sh1[0][1] * sh4[3][5]) + -kSqrt01_16 * (sh1[2][1] * sh4[7][5] + sh1[0][1] * sh4[1][5]);
        sh5[7][7] = sh1[1][1] * sh4[6][6] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh4[5][6] - sh1[0][1] * sh4[3][6]) + -sqrt(1.0 / 14.0) * (sh1[2][1] * sh4[7][6] + sh1[0][1] * sh4[1][6]);
        sh5[7][8] = sqrt(21.0 / 16.0) * sh1[1][1] * sh4[6][7] + kSqrt21_32 * (sh1[2][1] * sh4[5][7] - sh1[0][1] * sh4[3][7]) + -kSqrt03_32 * (sh1[2][1] * sh4[7][7] + sh1[0][1] * sh4[1][7]);
        sh5[7][9] = sqrt(7.0 / 3.0) * sh1[1][1] * sh4[6][8] + sqrt(7.0 / 6.0) * (sh1[2][1] * sh4[5][8] - sh1[0][1] * sh4[3][8]) + -kSqrt01_06 * (sh1[2][1] * sh4[7][8] + sh1[0][1] * sh4[1][8]);
        sh5[7][10] = sqrt(7.0 / 30.0) * (sh1[1][2] * sh4[6][8] - sh1[1][0] * sh4[6][0]) + sqrt(7.0 / 60.0) * ((sh1[2][2] * sh4[5][8] - sh1[2][0] * sh4[5][0]) - (sh1[0][2] * sh4[3][8] - sh1[0][0] * sh4[3][0])) + -sqrt(1.0 / 60.0) * ((sh1[2][2] * sh4[7][8] - sh1[2][0] * sh4[7][0]) + (sh1[0][2] * sh4[1][8] - sh1[0][0] * sh4[1][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[7]);

        sh5[8][0] = sqrt(8.0 / 45.0) * (sh1[1][2] * sh4[7][0] + sh1[1][0] * sh4[7][8]) + sqrt(7.0 / 45.0) * ((sh1[2][2] * sh4[6][0] + sh1[2][0] * sh4[6][8]) - (sh1[0][2] * sh4[2][0] + sh1[0][0] * sh4[2][8])) + -sqrt(1.0 / 180.0) * ((sh1[2][2] * sh4[8][0] + sh1[2][0] * sh4[8][8]) + (sh1[0][2] * sh4[0][0] + sh1[0][0] * sh4[0][8]));
        sh5[8][1] = sqrt(16.0 / 9.0) * sh1[1][1] * sh4[7][0] + sqrt(14.0 / 9.0) * (sh1[2][1] * sh4[6][0] - sh1[0][1] * sh4[2][0]) + -kSqrt01_18 * (sh1[2][1] * sh4[8][0] + sh1[0][1] * sh4[0][0]);
        sh5[8][2] = sh1[1][1] * sh4[7][1] + kSqrt07_08 * (sh1[2][1] * sh4[6][1] - sh1[0][1] * sh4[2][1]) + -kSqrt01_32 * (sh1[2][1] * sh4[8][1] + sh1[0][1] * sh4[0][1]);
        sh5[8][3] = sqrt(16.0 / 21.0) * sh1[1][1] * sh4[7][2] + kSqrt02_03 * (sh1[2][1] * sh4[6][2] - sh1[0][1] * sh4[2][2]) + -sqrt(1.0 / 42.0) * (sh1[2][1] * sh4[8][2] + sh1[0][1] * sh4[0][2]);
        sh5[8][4] = kSqrt02_03 * sh1[1][1] * sh4[7][3] + sqrt(7.0 / 12.0) * (sh1[2][1] * sh4[6][3] - sh1[0][1] * sh4[2][3]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh4[8][3] + sh1[0][1] * sh4[0][3]);
        sh5[8][5] = sqrt(16.0 / 25.0) * sh1[1][1] * sh4[7][4] + kSqrt14_25 * (sh1[2][1] * sh4[6][4] - sh1[0][1] * sh4[2][4]) + -kSqrt01_50 * (sh1[2][1] * sh4[8][4] + sh1[0][1] * sh4[0][4]);
        sh5[8][6] = kSqrt02_03 * sh1[1][1] * sh4[7][5] + sqrt(7.0 / 12.0) * (sh1[2][1] * sh4[6][5] - sh1[0][1] * sh4[2][5]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh4[8][5] + sh1[0][1] * sh4[0][5]);
        sh5[8][7] = sqrt(16.0 / 21.0) * sh1[1][1] * sh4[7][6] + kSqrt02_03 * (sh1[2][1] * sh4[6][6] - sh1[0][1] * sh4[2][6]) + -sqrt(1.0 / 42.0) * (sh1[2][1] * sh4[8][6] + sh1[0][1] * sh4[0][6]);
        sh5[8][8] = sh1[1][1] * sh4[7][7] + kSqrt07_08 * (sh1[2][1] * sh4[6][7] - sh1[0][1] * sh4[2][7]) + -kSqrt01_32 * (sh1[2][1] * sh4[8][7] + sh1[0][1] * sh4[0][7]);
        sh5[8][9] = sqrt(16.0 / 9.0) * sh1[1][1] * sh4[7][8] + sqrt(14.0 / 9.0) * (sh1[2][1] * sh4[6][8] - sh1[0][1] * sh4[2][8]) + -kSqrt01_18 * (sh1[2][1] * sh4[8][8] + sh1[0][1] * sh4[0][8]);
        sh5[8][10] = sqrt(8.0 / 45.0) * (sh1[1][2] * sh4[7][8] - sh1[1][0] * sh4[7][0]) + sqrt(7.0 / 45.0) * ((sh1[2][2] * sh4[6][8] - sh1[2][0] * sh4[6][0]) - (sh1[0][2] * sh4[2][8] - sh1[0][0] * sh4[2][0])) + -sqrt(1.0 / 180.0) * ((sh1[2][2] * sh4[8][8] - sh1[2][0] * sh4[8][0]) + (sh1[0][2] * sh4[0][8] - sh1[0][0] * sh4[0][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[8]);

        sh5[9][0] = kSqrt01_10 * (sh1[1][2] * sh4[8][0] + sh1[1][0] * sh4[8][8]) + kSqrt01_05 * ((sh1[2][2] * sh4[7][0] + sh1[2][0] * sh4[7][8]) - (sh1[0][2] * sh4[1][0] + sh1[0][0] * sh4[1][8]));
        sh5[9][1] = sh1[1][1] * sh4[8][0] + sqrt(2.0 / 1.0) * (sh1[2][1] * sh4[7][0] - sh1[0][1] * sh4[1][0]);
        sh5[9][2] = sqrt(9.0 / 16.0) * sh1[1][1] * sh4[8][1] + kSqrt09_08 * (sh1[2][1] * sh4[7][1] - sh1[0][1] * sh4[1][1]);
        sh5[9][3] = sqrt(3.0 / 7.0) * sh1[1][1] * sh4[8][2] + sqrt(6.0 / 7.0) * (sh1[2][1] * sh4[7][2] - sh1[0][1] * sh4[1][2]);
        sh5[9][4] = kSqrt03_08 * sh1[1][1] * sh4[8][3] + kSqrt03_04 * (sh1[2][1] * sh4[7][3] - sh1[0][1] * sh4[1][3]);
        sh5[9][5] = sqrt(9.0 / 25.0) * sh1[1][1] * sh4[8][4] + kSqrt18_25 * (sh1[2][1] * sh4[7][4] - sh1[0][1] * sh4[1][4]);
        sh5[9][6] = kSqrt03_08 * sh1[1][1] * sh4[8][5] + kSqrt03_04 * (sh1[2][1] * sh4[7][5] - sh1[0][1] * sh4[1][5]);
        sh5[9][7] = sqrt(3.0 / 7.0) * sh1[1][1] * sh4[8][6] + sqrt(6.0 / 7.0) * (sh1[2][1] * sh4[7][6] - sh1[0][1] * sh4[1][6]);
        sh5[9][8] = sqrt(9.0 / 16.0) * sh1[1][1] * sh4[8][7] + kSqrt09_08 * (sh1[2][1] * sh4[7][7] - sh1[0][1] * sh4[1][7]);
        sh5[9][9] = sh1[1][1] * sh4[8][8] + sqrt(2.0 / 1.0) * (sh1[2][1] * sh4[7][8] - sh1[0][1] * sh4[1][8]);
        sh5[9][10] = kSqrt01_10 * (sh1[1][2] * sh4[8][8] - sh1[1][0] * sh4[8][0]) + kSqrt01_05 * ((sh1[2][2] * sh4[7][8] - sh1[2][0] * sh4[7][0]) - (sh1[0][2] * sh4[1][8] - sh1[0][0] * sh4[1][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[9]);

        sh5[10][0] = kSqrt01_04 * ((sh1[2][2] * sh4[8][0] + sh1[2][0] * sh4[8][8]) - (sh1[0][2] * sh4[0][0] + sh1[0][0] * sh4[0][8]));
        sh5[10][1] = sqrt(5.0 / 2.0) * (sh1[2][1] * sh4[8][0] - sh1[0][1] * sh4[0][0]);
        sh5[10][2] = sqrt(45.0 / 32.0) * (sh1[2][1] * sh4[8][1] - sh1[0][1] * sh4[0][1]);
        sh5[10][3] = sqrt(15.0 / 14.0) * (sh1[2][1] * sh4[8][2] - sh1[0][1] * sh4[0][2]);
        sh5[10][4] = kSqrt15_16 * (sh1[2][1] * sh4[8][3] - sh1[0][1] * sh4[0][3]);
        sh5[10][5] = sqrt(9.0 / 10.0) * (sh1[2][1] * sh4[8][4] - sh1[0][1] * sh4[0][4]);
        sh5[10][6] = kSqrt15_16 * (sh1[2][1] * sh4[8][5] - sh1[0][1] * sh4[0][5]);
        sh5[10][7] = sqrt(15.0 / 14.0) * (sh1[2][1] * sh4[8][6] - sh1[0][1] * sh4[0][6]);
        sh5[10][8] = sqrt(45.0 / 32.0) * (sh1[2][1] * sh4[8][7] - sh1[0][1] * sh4[0][7]);
        sh5[10][9] = sqrt(5.0 / 2.0) * (sh1[2][1] * sh4[8][8] - sh1[0][1] * sh4[0][8]);
        sh5[10][10] = kSqrt01_04 * ((sh1[2][2] * sh4[8][8] - sh1[2][0] * sh4[8][0]) - (sh1[0][2] * sh4[0][8] - sh1[0][0] * sh4[0][0]));

        (*coeffs++) = dp(11, coeffsIn, sh5[10]);

    // band 7:

        if (n < 7)
            return;

        coeffsIn += 11;
        float sh6[13][13];

        sh6[0][0] = kSqrt01_04 * ((sh1[2][2] * sh5[0][0] + sh1[2][0] * sh5[0][10]) + (sh1[0][2] * sh5[10][0] + sh1[0][0] * sh5[10][10]));
        sh6[0][1] = sqrt(3.0 / 1.0) * (sh1[2][1] * sh5[0][0] + sh1[0][1] * sh5[10][0]);
        sh6[0][2] = sqrt(33.0 / 20.0) * (sh1[2][1] * sh5[0][1] + sh1[0][1] * sh5[10][1]);
        sh6[0][3] = sqrt(11.0 / 9.0) * (sh1[2][1] * sh5[0][2] + sh1[0][1] * sh5[10][2]);
        sh6[0][4] = sqrt(33.0 / 32.0) * (sh1[2][1] * sh5[0][3] + sh1[0][1] * sh5[10][3]);
        sh6[0][5] = sqrt(33.0 / 35.0) * (sh1[2][1] * sh5[0][4] + sh1[0][1] * sh5[10][4]);
        sh6[0][6] = sqrt(11.0 / 12.0) * (sh1[2][1] * sh5[0][5] + sh1[0][1] * sh5[10][5]);
        sh6[0][7] = sqrt(33.0 / 35.0) * (sh1[2][1] * sh5[0][6] + sh1[0][1] * sh5[10][6]);
        sh6[0][8] = sqrt(33.0 / 32.0) * (sh1[2][1] * sh5[0][7] + sh1[0][1] * sh5[10][7]);
        sh6[0][9] = sqrt(11.0 / 9.0) * (sh1[2][1] * sh5[0][8] + sh1[0][1] * sh5[10][8]);
        sh6[0][10] = sqrt(33.0 / 20.0) * (sh1[2][1] * sh5[0][9] + sh1[0][1] * sh5[10][9]);
        sh6[0][11] = sqrt(3.0 / 1.0) * (sh1[2][1] * sh5[0][10] + sh1[0][1] * sh5[10][10]);
        sh6[0][12] = kSqrt01_04 * ((sh1[2][2] * sh5[0][10] - sh1[2][0] * sh5[0][0]) + (sh1[0][2] * sh5[10][10] - sh1[0][0] * sh5[10][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[0]);

        sh6[1][0] = kSqrt01_12 * (sh1[1][2] * sh5[0][0] + sh1[1][0] * sh5[0][10]) + sqrt(5.0 / 24.0) * ((sh1[2][2] * sh5[1][0] + sh1[2][0] * sh5[1][10]) + (sh1[0][2] * sh5[9][0] + sh1[0][0] * sh5[9][10]));
        sh6[1][1] = sh1[1][1] * sh5[0][0] + sqrt(5.0 / 2.0) * (sh1[2][1] * sh5[1][0] + sh1[0][1] * sh5[9][0]);
        sh6[1][2] = sqrt(11.0 / 20.0) * sh1[1][1] * sh5[0][1] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh5[1][1] + sh1[0][1] * sh5[9][1]);
        sh6[1][3] = sqrt(11.0 / 27.0) * sh1[1][1] * sh5[0][2] + sqrt(55.0 / 54.0) * (sh1[2][1] * sh5[1][2] + sh1[0][1] * sh5[9][2]);
        sh6[1][4] = sqrt(11.0 / 32.0) * sh1[1][1] * sh5[0][3] + sqrt(55.0 / 64.0) * (sh1[2][1] * sh5[1][3] + sh1[0][1] * sh5[9][3]);
        sh6[1][5] = sqrt(11.0 / 35.0) * sh1[1][1] * sh5[0][4] + sqrt(11.0 / 14.0) * (sh1[2][1] * sh5[1][4] + sh1[0][1] * sh5[9][4]);
        sh6[1][6] = sqrt(11.0 / 36.0) * sh1[1][1] * sh5[0][5] + sqrt(55.0 / 72.0) * (sh1[2][1] * sh5[1][5] + sh1[0][1] * sh5[9][5]);
        sh6[1][7] = sqrt(11.0 / 35.0) * sh1[1][1] * sh5[0][6] + sqrt(11.0 / 14.0) * (sh1[2][1] * sh5[1][6] + sh1[0][1] * sh5[9][6]);
        sh6[1][8] = sqrt(11.0 / 32.0) * sh1[1][1] * sh5[0][7] + sqrt(55.0 / 64.0) * (sh1[2][1] * sh5[1][7] + sh1[0][1] * sh5[9][7]);
        sh6[1][9] = sqrt(11.0 / 27.0) * sh1[1][1] * sh5[0][8] + sqrt(55.0 / 54.0) * (sh1[2][1] * sh5[1][8] + sh1[0][1] * sh5[9][8]);
        sh6[1][10] = sqrt(11.0 / 20.0) * sh1[1][1] * sh5[0][9] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh5[1][9] + sh1[0][1] * sh5[9][9]);
        sh6[1][11] = sh1[1][1] * sh5[0][10] + sqrt(5.0 / 2.0) * (sh1[2][1] * sh5[1][10] + sh1[0][1] * sh5[9][10]);
        sh6[1][12] = kSqrt01_12 * (sh1[1][2] * sh5[0][10] - sh1[1][0] * sh5[0][0]) + sqrt(5.0 / 24.0) * ((sh1[2][2] * sh5[1][10] - sh1[2][0] * sh5[1][0]) + (sh1[0][2] * sh5[9][10] - sh1[0][0] * sh5[9][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[1]);

        sh6[2][0] = sqrt(5.0 / 33.0) * (sh1[1][2] * sh5[1][0] + sh1[1][0] * sh5[1][10]) + sqrt(15.0 / 88.0) * ((sh1[2][2] * sh5[2][0] + sh1[2][0] * sh5[2][10]) + (sh1[0][2] * sh5[8][0] + sh1[0][0] * sh5[8][10])) + -sqrt(1.0 / 264.0) * ((sh1[2][2] * sh5[0][0] + sh1[2][0] * sh5[0][10]) - (sh1[0][2] * sh5[10][0] + sh1[0][0] * sh5[10][10]));
        sh6[2][1] = sqrt(20.0 / 11.0) * sh1[1][1] * sh5[1][0] + sqrt(45.0 / 22.0) * (sh1[2][1] * sh5[2][0] + sh1[0][1] * sh5[8][0]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh5[0][0] - sh1[0][1] * sh5[10][0]);
        sh6[2][2] = sh1[1][1] * sh5[1][1] + kSqrt09_08 * (sh1[2][1] * sh5[2][1] + sh1[0][1] * sh5[8][1]) + -sqrt(1.0 / 40.0) * (sh1[2][1] * sh5[0][1] - sh1[0][1] * sh5[10][1]);
        sh6[2][3] = sqrt(20.0 / 27.0) * sh1[1][1] * sh5[1][2] + kSqrt05_06 * (sh1[2][1] * sh5[2][2] + sh1[0][1] * sh5[8][2]) + -sqrt(1.0 / 54.0) * (sh1[2][1] * sh5[0][2] - sh1[0][1] * sh5[10][2]);
        sh6[2][4] = kSqrt05_08 * sh1[1][1] * sh5[1][3] + sqrt(45.0 / 64.0) * (sh1[2][1] * sh5[2][3] + sh1[0][1] * sh5[8][3]) + -sqrt(1.0 / 64.0) * (sh1[2][1] * sh5[0][3] - sh1[0][1] * sh5[10][3]);
        sh6[2][5] = sqrt(4.0 / 7.0) * sh1[1][1] * sh5[1][4] + sqrt(9.0 / 14.0) * (sh1[2][1] * sh5[2][4] + sh1[0][1] * sh5[8][4]) + -sqrt(1.0 / 70.0) * (sh1[2][1] * sh5[0][4] - sh1[0][1] * sh5[10][4]);
        sh6[2][6] = kSqrt05_09 * sh1[1][1] * sh5[1][5] + kSqrt05_08 * (sh1[2][1] * sh5[2][5] + sh1[0][1] * sh5[8][5]) + -sqrt(1.0 / 72.0) * (sh1[2][1] * sh5[0][5] - sh1[0][1] * sh5[10][5]);
        sh6[2][7] = sqrt(4.0 / 7.0) * sh1[1][1] * sh5[1][6] + sqrt(9.0 / 14.0) * (sh1[2][1] * sh5[2][6] + sh1[0][1] * sh5[8][6]) + -sqrt(1.0 / 70.0) * (sh1[2][1] * sh5[0][6] - sh1[0][1] * sh5[10][6]);
        sh6[2][8] = kSqrt05_08 * sh1[1][1] * sh5[1][7] + sqrt(45.0 / 64.0) * (sh1[2][1] * sh5[2][7] + sh1[0][1] * sh5[8][7]) + -sqrt(1.0 / 64.0) * (sh1[2][1] * sh5[0][7] - sh1[0][1] * sh5[10][7]);
        sh6[2][9] = sqrt(20.0 / 27.0) * sh1[1][1] * sh5[1][8] + kSqrt05_06 * (sh1[2][1] * sh5[2][8] + sh1[0][1] * sh5[8][8]) + -sqrt(1.0 / 54.0) * (sh1[2][1] * sh5[0][8] - sh1[0][1] * sh5[10][8]);
        sh6[2][10] = sh1[1][1] * sh5[1][9] + kSqrt09_08 * (sh1[2][1] * sh5[2][9] + sh1[0][1] * sh5[8][9]) + -sqrt(1.0 / 40.0) * (sh1[2][1] * sh5[0][9] - sh1[0][1] * sh5[10][9]);
        sh6[2][11] = sqrt(20.0 / 11.0) * sh1[1][1] * sh5[1][10] + sqrt(45.0 / 22.0) * (sh1[2][1] * sh5[2][10] + sh1[0][1] * sh5[8][10]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh5[0][10] - sh1[0][1] * sh5[10][10]);
        sh6[2][12] = sqrt(5.0 / 33.0) * (sh1[1][2] * sh5[1][10] - sh1[1][0] * sh5[1][0]) + sqrt(15.0 / 88.0) * ((sh1[2][2] * sh5[2][10] - sh1[2][0] * sh5[2][0]) + (sh1[0][2] * sh5[8][10] - sh1[0][0] * sh5[8][0])) + -sqrt(1.0 / 264.0) * ((sh1[2][2] * sh5[0][10] - sh1[2][0] * sh5[0][0]) - (sh1[0][2] * sh5[10][10] - sh1[0][0] * sh5[10][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[2]);

        sh6[3][0] = sqrt(9.0 / 44.0) * (sh1[1][2] * sh5[2][0] + sh1[1][0] * sh5[2][10]) + sqrt(3.0 / 22.0) * ((sh1[2][2] * sh5[3][0] + sh1[2][0] * sh5[3][10]) + (sh1[0][2] * sh5[7][0] + sh1[0][0] * sh5[7][10])) + -sqrt(1.0 / 88.0) * ((sh1[2][2] * sh5[1][0] + sh1[2][0] * sh5[1][10]) - (sh1[0][2] * sh5[9][0] + sh1[0][0] * sh5[9][10]));
        sh6[3][1] = sqrt(27.0 / 11.0) * sh1[1][1] * sh5[2][0] + sqrt(18.0 / 11.0) * (sh1[2][1] * sh5[3][0] + sh1[0][1] * sh5[7][0]) + -sqrt(3.0 / 22.0) * (sh1[2][1] * sh5[1][0] - sh1[0][1] * sh5[9][0]);
        sh6[3][2] = sqrt(27.0 / 20.0) * sh1[1][1] * sh5[2][1] + sqrt(9.0 / 10.0) * (sh1[2][1] * sh5[3][1] + sh1[0][1] * sh5[7][1]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh5[1][1] - sh1[0][1] * sh5[9][1]);
        sh6[3][3] = sh1[1][1] * sh5[2][2] + kSqrt02_03 * (sh1[2][1] * sh5[3][2] + sh1[0][1] * sh5[7][2]) + -kSqrt01_18 * (sh1[2][1] * sh5[1][2] - sh1[0][1] * sh5[9][2]);
        sh6[3][4] = sqrt(27.0 / 32.0) * sh1[1][1] * sh5[2][3] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh5[3][3] + sh1[0][1] * sh5[7][3]) + -sqrt(3.0 / 64.0) * (sh1[2][1] * sh5[1][3] - sh1[0][1] * sh5[9][3]);
        sh6[3][5] = sqrt(27.0 / 35.0) * sh1[1][1] * sh5[2][4] + sqrt(18.0 / 35.0) * (sh1[2][1] * sh5[3][4] + sh1[0][1] * sh5[7][4]) + -sqrt(3.0 / 70.0) * (sh1[2][1] * sh5[1][4] - sh1[0][1] * sh5[9][4]);
        sh6[3][6] = kSqrt03_04 * sh1[1][1] * sh5[2][5] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh5[3][5] + sh1[0][1] * sh5[7][5]) + -sqrt(1.0 / 24.0) * (sh1[2][1] * sh5[1][5] - sh1[0][1] * sh5[9][5]);
        sh6[3][7] = sqrt(27.0 / 35.0) * sh1[1][1] * sh5[2][6] + sqrt(18.0 / 35.0) * (sh1[2][1] * sh5[3][6] + sh1[0][1] * sh5[7][6]) + -sqrt(3.0 / 70.0) * (sh1[2][1] * sh5[1][6] - sh1[0][1] * sh5[9][6]);
        sh6[3][8] = sqrt(27.0 / 32.0) * sh1[1][1] * sh5[2][7] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh5[3][7] + sh1[0][1] * sh5[7][7]) + -sqrt(3.0 / 64.0) * (sh1[2][1] * sh5[1][7] - sh1[0][1] * sh5[9][7]);
        sh6[3][9] = sh1[1][1] * sh5[2][8] + kSqrt02_03 * (sh1[2][1] * sh5[3][8] + sh1[0][1] * sh5[7][8]) + -kSqrt01_18 * (sh1[2][1] * sh5[1][8] - sh1[0][1] * sh5[9][8]);
        sh6[3][10] = sqrt(27.0 / 20.0) * sh1[1][1] * sh5[2][9] + sqrt(9.0 / 10.0) * (sh1[2][1] * sh5[3][9] + sh1[0][1] * sh5[7][9]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh5[1][9] - sh1[0][1] * sh5[9][9]);
        sh6[3][11] = sqrt(27.0 / 11.0) * sh1[1][1] * sh5[2][10] + sqrt(18.0 / 11.0) * (sh1[2][1] * sh5[3][10] + sh1[0][1] * sh5[7][10]) + -sqrt(3.0 / 22.0) * (sh1[2][1] * sh5[1][10] - sh1[0][1] * sh5[9][10]);
        sh6[3][12] = sqrt(9.0 / 44.0) * (sh1[1][2] * sh5[2][10] - sh1[1][0] * sh5[2][0]) + sqrt(3.0 / 22.0) * ((sh1[2][2] * sh5[3][10] - sh1[2][0] * sh5[3][0]) + (sh1[0][2] * sh5[7][10] - sh1[0][0] * sh5[7][0])) + -sqrt(1.0 / 88.0) * ((sh1[2][2] * sh5[1][10] - sh1[2][0] * sh5[1][0]) - (sh1[0][2] * sh5[9][10] - sh1[0][0] * sh5[9][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[3]);

        sh6[4][0] = sqrt(8.0 / 33.0) * (sh1[1][2] * sh5[3][0] + sh1[1][0] * sh5[3][10]) + sqrt(7.0 / 66.0) * ((sh1[2][2] * sh5[4][0] + sh1[2][0] * sh5[4][10]) + (sh1[0][2] * sh5[6][0] + sh1[0][0] * sh5[6][10])) + -sqrt(1.0 / 44.0) * ((sh1[2][2] * sh5[2][0] + sh1[2][0] * sh5[2][10]) - (sh1[0][2] * sh5[8][0] + sh1[0][0] * sh5[8][10]));
        sh6[4][1] = sqrt(32.0 / 11.0) * sh1[1][1] * sh5[3][0] + sqrt(14.0 / 11.0) * (sh1[2][1] * sh5[4][0] + sh1[0][1] * sh5[6][0]) + -sqrt(3.0 / 11.0) * (sh1[2][1] * sh5[2][0] - sh1[0][1] * sh5[8][0]);
        sh6[4][2] = kSqrt08_05 * sh1[1][1] * sh5[3][1] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh5[4][1] + sh1[0][1] * sh5[6][1]) + -sqrt(3.0 / 20.0) * (sh1[2][1] * sh5[2][1] - sh1[0][1] * sh5[8][1]);
        sh6[4][3] = sqrt(32.0 / 27.0) * sh1[1][1] * sh5[3][2] + sqrt(14.0 / 27.0) * (sh1[2][1] * sh5[4][2] + sh1[0][1] * sh5[6][2]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh5[2][2] - sh1[0][1] * sh5[8][2]);
        sh6[4][4] = sh1[1][1] * sh5[3][3] + kSqrt07_16 * (sh1[2][1] * sh5[4][3] + sh1[0][1] * sh5[6][3]) + -kSqrt03_32 * (sh1[2][1] * sh5[2][3] - sh1[0][1] * sh5[8][3]);
        sh6[4][5] = sqrt(32.0 / 35.0) * sh1[1][1] * sh5[3][4] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh5[4][4] + sh1[0][1] * sh5[6][4]) + -sqrt(3.0 / 35.0) * (sh1[2][1] * sh5[2][4] - sh1[0][1] * sh5[8][4]);
        sh6[4][6] = kSqrt08_09 * sh1[1][1] * sh5[3][5] + sqrt(7.0 / 18.0) * (sh1[2][1] * sh5[4][5] + sh1[0][1] * sh5[6][5]) + -kSqrt01_12 * (sh1[2][1] * sh5[2][5] - sh1[0][1] * sh5[8][5]);
        sh6[4][7] = sqrt(32.0 / 35.0) * sh1[1][1] * sh5[3][6] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh5[4][6] + sh1[0][1] * sh5[6][6]) + -sqrt(3.0 / 35.0) * (sh1[2][1] * sh5[2][6] - sh1[0][1] * sh5[8][6]);
        sh6[4][8] = sh1[1][1] * sh5[3][7] + kSqrt07_16 * (sh1[2][1] * sh5[4][7] + sh1[0][1] * sh5[6][7]) + -kSqrt03_32 * (sh1[2][1] * sh5[2][7] - sh1[0][1] * sh5[8][7]);
        sh6[4][9] = sqrt(32.0 / 27.0) * sh1[1][1] * sh5[3][8] + sqrt(14.0 / 27.0) * (sh1[2][1] * sh5[4][8] + sh1[0][1] * sh5[6][8]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh5[2][8] - sh1[0][1] * sh5[8][8]);
        sh6[4][10] = kSqrt08_05 * sh1[1][1] * sh5[3][9] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh5[4][9] + sh1[0][1] * sh5[6][9]) + -sqrt(3.0 / 20.0) * (sh1[2][1] * sh5[2][9] - sh1[0][1] * sh5[8][9]);
        sh6[4][11] = sqrt(32.0 / 11.0) * sh1[1][1] * sh5[3][10] + sqrt(14.0 / 11.0) * (sh1[2][1] * sh5[4][10] + sh1[0][1] * sh5[6][10]) + -sqrt(3.0 / 11.0) * (sh1[2][1] * sh5[2][10] - sh1[0][1] * sh5[8][10]);
        sh6[4][12] = sqrt(8.0 / 33.0) * (sh1[1][2] * sh5[3][10] - sh1[1][0] * sh5[3][0]) + sqrt(7.0 / 66.0) * ((sh1[2][2] * sh5[4][10] - sh1[2][0] * sh5[4][0]) + (sh1[0][2] * sh5[6][10] - sh1[0][0] * sh5[6][0])) + -sqrt(1.0 / 44.0) * ((sh1[2][2] * sh5[2][10] - sh1[2][0] * sh5[2][0]) - (sh1[0][2] * sh5[8][10] - sh1[0][0] * sh5[8][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[4]);

        sh6[5][0] = sqrt(35.0 / 132.0) * (sh1[1][2] * sh5[4][0] + sh1[1][0] * sh5[4][10]) + sqrt(7.0 / 44.0) * (sh1[0][2] * sh5[5][0] + sh1[0][0] * sh5[5][10]) + -sqrt(5.0 / 132.0) * ((sh1[2][2] * sh5[3][0] + sh1[2][0] * sh5[3][10]) - (sh1[0][2] * sh5[7][0] + sh1[0][0] * sh5[7][10]));
        sh6[5][1] = sqrt(35.0 / 11.0) * sh1[1][1] * sh5[4][0] + sqrt(21.0 / 11.0) * sh1[0][1] * sh5[5][0] + -sqrt(5.0 / 11.0) * (sh1[2][1] * sh5[3][0] - sh1[0][1] * sh5[7][0]);
        sh6[5][2] = sqrt(7.0 / 4.0) * sh1[1][1] * sh5[4][1] + sqrt(21.0 / 20.0) * sh1[0][1] * sh5[5][1] + -kSqrt01_04 * (sh1[2][1] * sh5[3][1] - sh1[0][1] * sh5[7][1]);
        sh6[5][3] = sqrt(35.0 / 27.0) * sh1[1][1] * sh5[4][2] + sqrt(7.0 / 9.0) * sh1[0][1] * sh5[5][2] + -sqrt(5.0 / 27.0) * (sh1[2][1] * sh5[3][2] - sh1[0][1] * sh5[7][2]);
        sh6[5][4] = sqrt(35.0 / 32.0) * sh1[1][1] * sh5[4][3] + kSqrt21_32 * sh1[0][1] * sh5[5][3] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh5[3][3] - sh1[0][1] * sh5[7][3]);
        sh6[5][5] = sh1[1][1] * sh5[4][4] + kSqrt03_05 * sh1[0][1] * sh5[5][4] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh5[3][4] - sh1[0][1] * sh5[7][4]);
        sh6[5][6] = sqrt(35.0 / 36.0) * sh1[1][1] * sh5[4][5] + sqrt(7.0 / 12.0) * sh1[0][1] * sh5[5][5] + -sqrt(5.0 / 36.0) * (sh1[2][1] * sh5[3][5] - sh1[0][1] * sh5[7][5]);
        sh6[5][7] = sh1[1][1] * sh5[4][6] + kSqrt03_05 * sh1[0][1] * sh5[5][6] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh5[3][6] - sh1[0][1] * sh5[7][6]);
        sh6[5][8] = sqrt(35.0 / 32.0) * sh1[1][1] * sh5[4][7] + kSqrt21_32 * sh1[0][1] * sh5[5][7] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh5[3][7] - sh1[0][1] * sh5[7][7]);
        sh6[5][9] = sqrt(35.0 / 27.0) * sh1[1][1] * sh5[4][8] + sqrt(7.0 / 9.0) * sh1[0][1] * sh5[5][8] + -sqrt(5.0 / 27.0) * (sh1[2][1] * sh5[3][8] - sh1[0][1] * sh5[7][8]);
        sh6[5][10] = sqrt(7.0 / 4.0) * sh1[1][1] * sh5[4][9] + sqrt(21.0 / 20.0) * sh1[0][1] * sh5[5][9] + -kSqrt01_04 * (sh1[2][1] * sh5[3][9] - sh1[0][1] * sh5[7][9]);
        sh6[5][11] = sqrt(35.0 / 11.0) * sh1[1][1] * sh5[4][10] + sqrt(21.0 / 11.0) * sh1[0][1] * sh5[5][10] + -sqrt(5.0 / 11.0) * (sh1[2][1] * sh5[3][10] - sh1[0][1] * sh5[7][10]);
        sh6[5][12] = sqrt(35.0 / 132.0) * (sh1[1][2] * sh5[4][10] - sh1[1][0] * sh5[4][0]) + sqrt(7.0 / 44.0) * (sh1[0][2] * sh5[5][10] - sh1[0][0] * sh5[5][0]) + -sqrt(5.0 / 132.0) * ((sh1[2][2] * sh5[3][10] - sh1[2][0] * sh5[3][0]) - (sh1[0][2] * sh5[7][10] - sh1[0][0] * sh5[7][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[5]);

        sh6[6][0] = sqrt(3.0 / 11.0) * (sh1[1][2] * sh5[5][0] + sh1[1][0] * sh5[5][10]) + -sqrt(5.0 / 44.0) * ((sh1[2][2] * sh5[6][0] + sh1[2][0] * sh5[6][10]) + (sh1[0][2] * sh5[4][0] + sh1[0][0] * sh5[4][10]));
        sh6[6][1] = sqrt(36.0 / 11.0) * sh1[1][1] * sh5[5][0] + -sqrt(15.0 / 11.0) * (sh1[2][1] * sh5[6][0] + sh1[0][1] * sh5[4][0]);
        sh6[6][2] = kSqrt09_05 * sh1[1][1] * sh5[5][1] + -kSqrt03_04 * (sh1[2][1] * sh5[6][1] + sh1[0][1] * sh5[4][1]);
        sh6[6][3] = kSqrt04_03 * sh1[1][1] * sh5[5][2] + -kSqrt05_09 * (sh1[2][1] * sh5[6][2] + sh1[0][1] * sh5[4][2]);
        sh6[6][4] = kSqrt09_08 * sh1[1][1] * sh5[5][3] + -kSqrt15_32 * (sh1[2][1] * sh5[6][3] + sh1[0][1] * sh5[4][3]);
        sh6[6][5] = sqrt(36.0 / 35.0) * sh1[1][1] * sh5[5][4] + -sqrt(3.0 / 7.0) * (sh1[2][1] * sh5[6][4] + sh1[0][1] * sh5[4][4]);
        sh6[6][6] = sh1[1][1] * sh5[5][5] + -sqrt(5.0 / 12.0) * (sh1[2][1] * sh5[6][5] + sh1[0][1] * sh5[4][5]);
        sh6[6][7] = sqrt(36.0 / 35.0) * sh1[1][1] * sh5[5][6] + -sqrt(3.0 / 7.0) * (sh1[2][1] * sh5[6][6] + sh1[0][1] * sh5[4][6]);
        sh6[6][8] = kSqrt09_08 * sh1[1][1] * sh5[5][7] + -kSqrt15_32 * (sh1[2][1] * sh5[6][7] + sh1[0][1] * sh5[4][7]);
        sh6[6][9] = kSqrt04_03 * sh1[1][1] * sh5[5][8] + -kSqrt05_09 * (sh1[2][1] * sh5[6][8] + sh1[0][1] * sh5[4][8]);
        sh6[6][10] = kSqrt09_05 * sh1[1][1] * sh5[5][9] + -kSqrt03_04 * (sh1[2][1] * sh5[6][9] + sh1[0][1] * sh5[4][9]);
        sh6[6][11] = sqrt(36.0 / 11.0) * sh1[1][1] * sh5[5][10] + -sqrt(15.0 / 11.0) * (sh1[2][1] * sh5[6][10] + sh1[0][1] * sh5[4][10]);
        sh6[6][12] = sqrt(3.0 / 11.0) * (sh1[1][2] * sh5[5][10] - sh1[1][0] * sh5[5][0]) + -sqrt(5.0 / 44.0) * ((sh1[2][2] * sh5[6][10] - sh1[2][0] * sh5[6][0]) + (sh1[0][2] * sh5[4][10] - sh1[0][0] * sh5[4][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[6]);

        sh6[7][0] = sqrt(35.0 / 132.0) * (sh1[1][2] * sh5[6][0] + sh1[1][0] * sh5[6][10]) + sqrt(7.0 / 44.0) * (sh1[2][2] * sh5[5][0] + sh1[2][0] * sh5[5][10]) + -sqrt(5.0 / 132.0) * ((sh1[2][2] * sh5[7][0] + sh1[2][0] * sh5[7][10]) + (sh1[0][2] * sh5[3][0] + sh1[0][0] * sh5[3][10]));
        sh6[7][1] = sqrt(35.0 / 11.0) * sh1[1][1] * sh5[6][0] + sqrt(21.0 / 11.0) * sh1[2][1] * sh5[5][0] + -sqrt(5.0 / 11.0) * (sh1[2][1] * sh5[7][0] + sh1[0][1] * sh5[3][0]);
        sh6[7][2] = sqrt(7.0 / 4.0) * sh1[1][1] * sh5[6][1] + sqrt(21.0 / 20.0) * sh1[2][1] * sh5[5][1] + -kSqrt01_04 * (sh1[2][1] * sh5[7][1] + sh1[0][1] * sh5[3][1]);
        sh6[7][3] = sqrt(35.0 / 27.0) * sh1[1][1] * sh5[6][2] + sqrt(7.0 / 9.0) * sh1[2][1] * sh5[5][2] + -sqrt(5.0 / 27.0) * (sh1[2][1] * sh5[7][2] + sh1[0][1] * sh5[3][2]);
        sh6[7][4] = sqrt(35.0 / 32.0) * sh1[1][1] * sh5[6][3] + kSqrt21_32 * sh1[2][1] * sh5[5][3] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh5[7][3] + sh1[0][1] * sh5[3][3]);
        sh6[7][5] = sh1[1][1] * sh5[6][4] + kSqrt03_05 * sh1[2][1] * sh5[5][4] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh5[7][4] + sh1[0][1] * sh5[3][4]);
        sh6[7][6] = sqrt(35.0 / 36.0) * sh1[1][1] * sh5[6][5] + sqrt(7.0 / 12.0) * sh1[2][1] * sh5[5][5] + -sqrt(5.0 / 36.0) * (sh1[2][1] * sh5[7][5] + sh1[0][1] * sh5[3][5]);
        sh6[7][7] = sh1[1][1] * sh5[6][6] + kSqrt03_05 * sh1[2][1] * sh5[5][6] + -sqrt(1.0 / 7.0) * (sh1[2][1] * sh5[7][6] + sh1[0][1] * sh5[3][6]);
        sh6[7][8] = sqrt(35.0 / 32.0) * sh1[1][1] * sh5[6][7] + kSqrt21_32 * sh1[2][1] * sh5[5][7] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh5[7][7] + sh1[0][1] * sh5[3][7]);
        sh6[7][9] = sqrt(35.0 / 27.0) * sh1[1][1] * sh5[6][8] + sqrt(7.0 / 9.0) * sh1[2][1] * sh5[5][8] + -sqrt(5.0 / 27.0) * (sh1[2][1] * sh5[7][8] + sh1[0][1] * sh5[3][8]);
        sh6[7][10] = sqrt(7.0 / 4.0) * sh1[1][1] * sh5[6][9] + sqrt(21.0 / 20.0) * sh1[2][1] * sh5[5][9] + -kSqrt01_04 * (sh1[2][1] * sh5[7][9] + sh1[0][1] * sh5[3][9]);
        sh6[7][11] = sqrt(35.0 / 11.0) * sh1[1][1] * sh5[6][10] + sqrt(21.0 / 11.0) * sh1[2][1] * sh5[5][10] + -sqrt(5.0 / 11.0) * (sh1[2][1] * sh5[7][10] + sh1[0][1] * sh5[3][10]);
        sh6[7][12] = sqrt(35.0 / 132.0) * (sh1[1][2] * sh5[6][10] - sh1[1][0] * sh5[6][0]) + sqrt(7.0 / 44.0) * (sh1[2][2] * sh5[5][10] - sh1[2][0] * sh5[5][0]) + -sqrt(5.0 / 132.0) * ((sh1[2][2] * sh5[7][10] - sh1[2][0] * sh5[7][0]) + (sh1[0][2] * sh5[3][10] - sh1[0][0] * sh5[3][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[7]);

        sh6[8][0] = sqrt(8.0 / 33.0) * (sh1[1][2] * sh5[7][0] + sh1[1][0] * sh5[7][10]) + sqrt(7.0 / 66.0) * ((sh1[2][2] * sh5[6][0] + sh1[2][0] * sh5[6][10]) - (sh1[0][2] * sh5[4][0] + sh1[0][0] * sh5[4][10])) + -sqrt(1.0 / 44.0) * ((sh1[2][2] * sh5[8][0] + sh1[2][0] * sh5[8][10]) + (sh1[0][2] * sh5[2][0] + sh1[0][0] * sh5[2][10]));
        sh6[8][1] = sqrt(32.0 / 11.0) * sh1[1][1] * sh5[7][0] + sqrt(14.0 / 11.0) * (sh1[2][1] * sh5[6][0] - sh1[0][1] * sh5[4][0]) + -sqrt(3.0 / 11.0) * (sh1[2][1] * sh5[8][0] + sh1[0][1] * sh5[2][0]);
        sh6[8][2] = kSqrt08_05 * sh1[1][1] * sh5[7][1] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh5[6][1] - sh1[0][1] * sh5[4][1]) + -sqrt(3.0 / 20.0) * (sh1[2][1] * sh5[8][1] + sh1[0][1] * sh5[2][1]);
        sh6[8][3] = sqrt(32.0 / 27.0) * sh1[1][1] * sh5[7][2] + sqrt(14.0 / 27.0) * (sh1[2][1] * sh5[6][2] - sh1[0][1] * sh5[4][2]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh5[8][2] + sh1[0][1] * sh5[2][2]);
        sh6[8][4] = sh1[1][1] * sh5[7][3] + kSqrt07_16 * (sh1[2][1] * sh5[6][3] - sh1[0][1] * sh5[4][3]) + -kSqrt03_32 * (sh1[2][1] * sh5[8][3] + sh1[0][1] * sh5[2][3]);
        sh6[8][5] = sqrt(32.0 / 35.0) * sh1[1][1] * sh5[7][4] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh5[6][4] - sh1[0][1] * sh5[4][4]) + -sqrt(3.0 / 35.0) * (sh1[2][1] * sh5[8][4] + sh1[0][1] * sh5[2][4]);
        sh6[8][6] = kSqrt08_09 * sh1[1][1] * sh5[7][5] + sqrt(7.0 / 18.0) * (sh1[2][1] * sh5[6][5] - sh1[0][1] * sh5[4][5]) + -kSqrt01_12 * (sh1[2][1] * sh5[8][5] + sh1[0][1] * sh5[2][5]);
        sh6[8][7] = sqrt(32.0 / 35.0) * sh1[1][1] * sh5[7][6] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh5[6][6] - sh1[0][1] * sh5[4][6]) + -sqrt(3.0 / 35.0) * (sh1[2][1] * sh5[8][6] + sh1[0][1] * sh5[2][6]);
        sh6[8][8] = sh1[1][1] * sh5[7][7] + kSqrt07_16 * (sh1[2][1] * sh5[6][7] - sh1[0][1] * sh5[4][7]) + -kSqrt03_32 * (sh1[2][1] * sh5[8][7] + sh1[0][1] * sh5[2][7]);
        sh6[8][9] = sqrt(32.0 / 27.0) * sh1[1][1] * sh5[7][8] + sqrt(14.0 / 27.0) * (sh1[2][1] * sh5[6][8] - sh1[0][1] * sh5[4][8]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh5[8][8] + sh1[0][1] * sh5[2][8]);
        sh6[8][10] = kSqrt08_05 * sh1[1][1] * sh5[7][9] + sqrt(7.0 / 10.0) * (sh1[2][1] * sh5[6][9] - sh1[0][1] * sh5[4][9]) + -sqrt(3.0 / 20.0) * (sh1[2][1] * sh5[8][9] + sh1[0][1] * sh5[2][9]);
        sh6[8][11] = sqrt(32.0 / 11.0) * sh1[1][1] * sh5[7][10] + sqrt(14.0 / 11.0) * (sh1[2][1] * sh5[6][10] - sh1[0][1] * sh5[4][10]) + -sqrt(3.0 / 11.0) * (sh1[2][1] * sh5[8][10] + sh1[0][1] * sh5[2][10]);
        sh6[8][12] = sqrt(8.0 / 33.0) * (sh1[1][2] * sh5[7][10] - sh1[1][0] * sh5[7][0]) + sqrt(7.0 / 66.0) * ((sh1[2][2] * sh5[6][10] - sh1[2][0] * sh5[6][0]) - (sh1[0][2] * sh5[4][10] - sh1[0][0] * sh5[4][0])) + -sqrt(1.0 / 44.0) * ((sh1[2][2] * sh5[8][10] - sh1[2][0] * sh5[8][0]) + (sh1[0][2] * sh5[2][10] - sh1[0][0] * sh5[2][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[8]);

        sh6[9][0] = sqrt(9.0 / 44.0) * (sh1[1][2] * sh5[8][0] + sh1[1][0] * sh5[8][10]) + sqrt(3.0 / 22.0) * ((sh1[2][2] * sh5[7][0] + sh1[2][0] * sh5[7][10]) - (sh1[0][2] * sh5[3][0] + sh1[0][0] * sh5[3][10])) + -sqrt(1.0 / 88.0) * ((sh1[2][2] * sh5[9][0] + sh1[2][0] * sh5[9][10]) + (sh1[0][2] * sh5[1][0] + sh1[0][0] * sh5[1][10]));
        sh6[9][1] = sqrt(27.0 / 11.0) * sh1[1][1] * sh5[8][0] + sqrt(18.0 / 11.0) * (sh1[2][1] * sh5[7][0] - sh1[0][1] * sh5[3][0]) + -sqrt(3.0 / 22.0) * (sh1[2][1] * sh5[9][0] + sh1[0][1] * sh5[1][0]);
        sh6[9][2] = sqrt(27.0 / 20.0) * sh1[1][1] * sh5[8][1] + sqrt(9.0 / 10.0) * (sh1[2][1] * sh5[7][1] - sh1[0][1] * sh5[3][1]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh5[9][1] + sh1[0][1] * sh5[1][1]);
        sh6[9][3] = sh1[1][1] * sh5[8][2] + kSqrt02_03 * (sh1[2][1] * sh5[7][2] - sh1[0][1] * sh5[3][2]) + -kSqrt01_18 * (sh1[2][1] * sh5[9][2] + sh1[0][1] * sh5[1][2]);
        sh6[9][4] = sqrt(27.0 / 32.0) * sh1[1][1] * sh5[8][3] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh5[7][3] - sh1[0][1] * sh5[3][3]) + -sqrt(3.0 / 64.0) * (sh1[2][1] * sh5[9][3] + sh1[0][1] * sh5[1][3]);
        sh6[9][5] = sqrt(27.0 / 35.0) * sh1[1][1] * sh5[8][4] + sqrt(18.0 / 35.0) * (sh1[2][1] * sh5[7][4] - sh1[0][1] * sh5[3][4]) + -sqrt(3.0 / 70.0) * (sh1[2][1] * sh5[9][4] + sh1[0][1] * sh5[1][4]);
        sh6[9][6] = kSqrt03_04 * sh1[1][1] * sh5[8][5] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh5[7][5] - sh1[0][1] * sh5[3][5]) + -sqrt(1.0 / 24.0) * (sh1[2][1] * sh5[9][5] + sh1[0][1] * sh5[1][5]);
        sh6[9][7] = sqrt(27.0 / 35.0) * sh1[1][1] * sh5[8][6] + sqrt(18.0 / 35.0) * (sh1[2][1] * sh5[7][6] - sh1[0][1] * sh5[3][6]) + -sqrt(3.0 / 70.0) * (sh1[2][1] * sh5[9][6] + sh1[0][1] * sh5[1][6]);
        sh6[9][8] = sqrt(27.0 / 32.0) * sh1[1][1] * sh5[8][7] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh5[7][7] - sh1[0][1] * sh5[3][7]) + -sqrt(3.0 / 64.0) * (sh1[2][1] * sh5[9][7] + sh1[0][1] * sh5[1][7]);
        sh6[9][9] = sh1[1][1] * sh5[8][8] + kSqrt02_03 * (sh1[2][1] * sh5[7][8] - sh1[0][1] * sh5[3][8]) + -kSqrt01_18 * (sh1[2][1] * sh5[9][8] + sh1[0][1] * sh5[1][8]);
        sh6[9][10] = sqrt(27.0 / 20.0) * sh1[1][1] * sh5[8][9] + sqrt(9.0 / 10.0) * (sh1[2][1] * sh5[7][9] - sh1[0][1] * sh5[3][9]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh5[9][9] + sh1[0][1] * sh5[1][9]);
        sh6[9][11] = sqrt(27.0 / 11.0) * sh1[1][1] * sh5[8][10] + sqrt(18.0 / 11.0) * (sh1[2][1] * sh5[7][10] - sh1[0][1] * sh5[3][10]) + -sqrt(3.0 / 22.0) * (sh1[2][1] * sh5[9][10] + sh1[0][1] * sh5[1][10]);
        sh6[9][12] = sqrt(9.0 / 44.0) * (sh1[1][2] * sh5[8][10] - sh1[1][0] * sh5[8][0]) + sqrt(3.0 / 22.0) * ((sh1[2][2] * sh5[7][10] - sh1[2][0] * sh5[7][0]) - (sh1[0][2] * sh5[3][10] - sh1[0][0] * sh5[3][0])) + -sqrt(1.0 / 88.0) * ((sh1[2][2] * sh5[9][10] - sh1[2][0] * sh5[9][0]) + (sh1[0][2] * sh5[1][10] - sh1[0][0] * sh5[1][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[9]);

        sh6[10][0] = sqrt(5.0 / 33.0) * (sh1[1][2] * sh5[9][0] + sh1[1][0] * sh5[9][10]) + sqrt(15.0 / 88.0) * ((sh1[2][2] * sh5[8][0] + sh1[2][0] * sh5[8][10]) - (sh1[0][2] * sh5[2][0] + sh1[0][0] * sh5[2][10])) + -sqrt(1.0 / 264.0) * ((sh1[2][2] * sh5[10][0] + sh1[2][0] * sh5[10][10]) + (sh1[0][2] * sh5[0][0] + sh1[0][0] * sh5[0][10]));
        sh6[10][1] = sqrt(20.0 / 11.0) * sh1[1][1] * sh5[9][0] + sqrt(45.0 / 22.0) * (sh1[2][1] * sh5[8][0] - sh1[0][1] * sh5[2][0]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh5[10][0] + sh1[0][1] * sh5[0][0]);
        sh6[10][2] = sh1[1][1] * sh5[9][1] + kSqrt09_08 * (sh1[2][1] * sh5[8][1] - sh1[0][1] * sh5[2][1]) + -sqrt(1.0 / 40.0) * (sh1[2][1] * sh5[10][1] + sh1[0][1] * sh5[0][1]);
        sh6[10][3] = sqrt(20.0 / 27.0) * sh1[1][1] * sh5[9][2] + kSqrt05_06 * (sh1[2][1] * sh5[8][2] - sh1[0][1] * sh5[2][2]) + -sqrt(1.0 / 54.0) * (sh1[2][1] * sh5[10][2] + sh1[0][1] * sh5[0][2]);
        sh6[10][4] = kSqrt05_08 * sh1[1][1] * sh5[9][3] + sqrt(45.0 / 64.0) * (sh1[2][1] * sh5[8][3] - sh1[0][1] * sh5[2][3]) + -sqrt(1.0 / 64.0) * (sh1[2][1] * sh5[10][3] + sh1[0][1] * sh5[0][3]);
        sh6[10][5] = sqrt(4.0 / 7.0) * sh1[1][1] * sh5[9][4] + sqrt(9.0 / 14.0) * (sh1[2][1] * sh5[8][4] - sh1[0][1] * sh5[2][4]) + -sqrt(1.0 / 70.0) * (sh1[2][1] * sh5[10][4] + sh1[0][1] * sh5[0][4]);
        sh6[10][6] = kSqrt05_09 * sh1[1][1] * sh5[9][5] + kSqrt05_08 * (sh1[2][1] * sh5[8][5] - sh1[0][1] * sh5[2][5]) + -sqrt(1.0 / 72.0) * (sh1[2][1] * sh5[10][5] + sh1[0][1] * sh5[0][5]);
        sh6[10][7] = sqrt(4.0 / 7.0) * sh1[1][1] * sh5[9][6] + sqrt(9.0 / 14.0) * (sh1[2][1] * sh5[8][6] - sh1[0][1] * sh5[2][6]) + -sqrt(1.0 / 70.0) * (sh1[2][1] * sh5[10][6] + sh1[0][1] * sh5[0][6]);
        sh6[10][8] = kSqrt05_08 * sh1[1][1] * sh5[9][7] + sqrt(45.0 / 64.0) * (sh1[2][1] * sh5[8][7] - sh1[0][1] * sh5[2][7]) + -sqrt(1.0 / 64.0) * (sh1[2][1] * sh5[10][7] + sh1[0][1] * sh5[0][7]);
        sh6[10][9] = sqrt(20.0 / 27.0) * sh1[1][1] * sh5[9][8] + kSqrt05_06 * (sh1[2][1] * sh5[8][8] - sh1[0][1] * sh5[2][8]) + -sqrt(1.0 / 54.0) * (sh1[2][1] * sh5[10][8] + sh1[0][1] * sh5[0][8]);
        sh6[10][10] = sh1[1][1] * sh5[9][9] + kSqrt09_08 * (sh1[2][1] * sh5[8][9] - sh1[0][1] * sh5[2][9]) + -sqrt(1.0 / 40.0) * (sh1[2][1] * sh5[10][9] + sh1[0][1] * sh5[0][9]);
        sh6[10][11] = sqrt(20.0 / 11.0) * sh1[1][1] * sh5[9][10] + sqrt(45.0 / 22.0) * (sh1[2][1] * sh5[8][10] - sh1[0][1] * sh5[2][10]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh5[10][10] + sh1[0][1] * sh5[0][10]);
        sh6[10][12] = sqrt(5.0 / 33.0) * (sh1[1][2] * sh5[9][10] - sh1[1][0] * sh5[9][0]) + sqrt(15.0 / 88.0) * ((sh1[2][2] * sh5[8][10] - sh1[2][0] * sh5[8][0]) - (sh1[0][2] * sh5[2][10] - sh1[0][0] * sh5[2][0])) + -sqrt(1.0 / 264.0) * ((sh1[2][2] * sh5[10][10] - sh1[2][0] * sh5[10][0]) + (sh1[0][2] * sh5[0][10] - sh1[0][0] * sh5[0][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[10]);

        sh6[11][0] = kSqrt01_12 * (sh1[1][2] * sh5[10][0] + sh1[1][0] * sh5[10][10]) + sqrt(5.0 / 24.0) * ((sh1[2][2] * sh5[9][0] + sh1[2][0] * sh5[9][10]) - (sh1[0][2] * sh5[1][0] + sh1[0][0] * sh5[1][10]));
        sh6[11][1] = sh1[1][1] * sh5[10][0] + sqrt(5.0 / 2.0) * (sh1[2][1] * sh5[9][0] - sh1[0][1] * sh5[1][0]);
        sh6[11][2] = sqrt(11.0 / 20.0) * sh1[1][1] * sh5[10][1] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh5[9][1] - sh1[0][1] * sh5[1][1]);
        sh6[11][3] = sqrt(11.0 / 27.0) * sh1[1][1] * sh5[10][2] + sqrt(55.0 / 54.0) * (sh1[2][1] * sh5[9][2] - sh1[0][1] * sh5[1][2]);
        sh6[11][4] = sqrt(11.0 / 32.0) * sh1[1][1] * sh5[10][3] + sqrt(55.0 / 64.0) * (sh1[2][1] * sh5[9][3] - sh1[0][1] * sh5[1][3]);
        sh6[11][5] = sqrt(11.0 / 35.0) * sh1[1][1] * sh5[10][4] + sqrt(11.0 / 14.0) * (sh1[2][1] * sh5[9][4] - sh1[0][1] * sh5[1][4]);
        sh6[11][6] = sqrt(11.0 / 36.0) * sh1[1][1] * sh5[10][5] + sqrt(55.0 / 72.0) * (sh1[2][1] * sh5[9][5] - sh1[0][1] * sh5[1][5]);
        sh6[11][7] = sqrt(11.0 / 35.0) * sh1[1][1] * sh5[10][6] + sqrt(11.0 / 14.0) * (sh1[2][1] * sh5[9][6] - sh1[0][1] * sh5[1][6]);
        sh6[11][8] = sqrt(11.0 / 32.0) * sh1[1][1] * sh5[10][7] + sqrt(55.0 / 64.0) * (sh1[2][1] * sh5[9][7] - sh1[0][1] * sh5[1][7]);
        sh6[11][9] = sqrt(11.0 / 27.0) * sh1[1][1] * sh5[10][8] + sqrt(55.0 / 54.0) * (sh1[2][1] * sh5[9][8] - sh1[0][1] * sh5[1][8]);
        sh6[11][10] = sqrt(11.0 / 20.0) * sh1[1][1] * sh5[10][9] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh5[9][9] - sh1[0][1] * sh5[1][9]);
        sh6[11][11] = sh1[1][1] * sh5[10][10] + sqrt(5.0 / 2.0) * (sh1[2][1] * sh5[9][10] - sh1[0][1] * sh5[1][10]);
        sh6[11][12] = kSqrt01_12 * (sh1[1][2] * sh5[10][10] - sh1[1][0] * sh5[10][0]) + sqrt(5.0 / 24.0) * ((sh1[2][2] * sh5[9][10] - sh1[2][0] * sh5[9][0]) - (sh1[0][2] * sh5[1][10] - sh1[0][0] * sh5[1][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[11]);

        sh6[12][0] = kSqrt01_04 * ((sh1[2][2] * sh5[10][0] + sh1[2][0] * sh5[10][10]) - (sh1[0][2] * sh5[0][0] + sh1[0][0] * sh5[0][10]));
        sh6[12][1] = sqrt(3.0 / 1.0) * (sh1[2][1] * sh5[10][0] - sh1[0][1] * sh5[0][0]);
        sh6[12][2] = sqrt(33.0 / 20.0) * (sh1[2][1] * sh5[10][1] - sh1[0][1] * sh5[0][1]);
        sh6[12][3] = sqrt(11.0 / 9.0) * (sh1[2][1] * sh5[10][2] - sh1[0][1] * sh5[0][2]);
        sh6[12][4] = sqrt(33.0 / 32.0) * (sh1[2][1] * sh5[10][3] - sh1[0][1] * sh5[0][3]);
        sh6[12][5] = sqrt(33.0 / 35.0) * (sh1[2][1] * sh5[10][4] - sh1[0][1] * sh5[0][4]);
        sh6[12][6] = sqrt(11.0 / 12.0) * (sh1[2][1] * sh5[10][5] - sh1[0][1] * sh5[0][5]);
        sh6[12][7] = sqrt(33.0 / 35.0) * (sh1[2][1] * sh5[10][6] - sh1[0][1] * sh5[0][6]);
        sh6[12][8] = sqrt(33.0 / 32.0) * (sh1[2][1] * sh5[10][7] - sh1[0][1] * sh5[0][7]);
        sh6[12][9] = sqrt(11.0 / 9.0) * (sh1[2][1] * sh5[10][8] - sh1[0][1] * sh5[0][8]);
        sh6[12][10] = sqrt(33.0 / 20.0) * (sh1[2][1] * sh5[10][9] - sh1[0][1] * sh5[0][9]);
        sh6[12][11] = sqrt(3.0 / 1.0) * (sh1[2][1] * sh5[10][10] - sh1[0][1] * sh5[0][10]);
        sh6[12][12] = kSqrt01_04 * ((sh1[2][2] * sh5[10][10] - sh1[2][0] * sh5[10][0]) - (sh1[0][2] * sh5[0][10] - sh1[0][0] * sh5[0][0]));

        (*coeffs++) = dp(13, coeffsIn, sh6[12]);

    // band 8:

        if (n < 8)
            return;

        coeffsIn += 13;
        float sh7[15][15];

        sh7[0][0] = kSqrt01_04 * ((sh1[2][2] * sh6[0][0] + sh1[2][0] * sh6[0][12]) + (sh1[0][2] * sh6[12][0] + sh1[0][0] * sh6[12][12]));
        sh7[0][1] = sqrt(7.0 / 2.0) * (sh1[2][1] * sh6[0][0] + sh1[0][1] * sh6[12][0]);
        sh7[0][2] = sqrt(91.0 / 48.0) * (sh1[2][1] * sh6[0][1] + sh1[0][1] * sh6[12][1]);
        sh7[0][3] = sqrt(91.0 / 66.0) * (sh1[2][1] * sh6[0][2] + sh1[0][1] * sh6[12][2]);
        sh7[0][4] = sqrt(91.0 / 80.0) * (sh1[2][1] * sh6[0][3] + sh1[0][1] * sh6[12][3]);
        sh7[0][5] = sqrt(91.0 / 90.0) * (sh1[2][1] * sh6[0][4] + sh1[0][1] * sh6[12][4]);
        sh7[0][6] = sqrt(91.0 / 96.0) * (sh1[2][1] * sh6[0][5] + sh1[0][1] * sh6[12][5]);
        sh7[0][7] = sqrt(13.0 / 14.0) * (sh1[2][1] * sh6[0][6] + sh1[0][1] * sh6[12][6]);
        sh7[0][8] = sqrt(91.0 / 96.0) * (sh1[2][1] * sh6[0][7] + sh1[0][1] * sh6[12][7]);
        sh7[0][9] = sqrt(91.0 / 90.0) * (sh1[2][1] * sh6[0][8] + sh1[0][1] * sh6[12][8]);
        sh7[0][10] = sqrt(91.0 / 80.0) * (sh1[2][1] * sh6[0][9] + sh1[0][1] * sh6[12][9]);
        sh7[0][11] = sqrt(91.0 / 66.0) * (sh1[2][1] * sh6[0][10] + sh1[0][1] * sh6[12][10]);
        sh7[0][12] = sqrt(91.0 / 48.0) * (sh1[2][1] * sh6[0][11] + sh1[0][1] * sh6[12][11]);
        sh7[0][13] = sqrt(7.0 / 2.0) * (sh1[2][1] * sh6[0][12] + sh1[0][1] * sh6[12][12]);
        sh7[0][14] = kSqrt01_04 * ((sh1[2][2] * sh6[0][12] - sh1[2][0] * sh6[0][0]) + (sh1[0][2] * sh6[12][12] - sh1[0][0] * sh6[12][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[0]);

        sh7[1][0] = sqrt(1.0 / 14.0) * (sh1[1][2] * sh6[0][0] + sh1[1][0] * sh6[0][12]) + sqrt(3.0 / 14.0) * ((sh1[2][2] * sh6[1][0] + sh1[2][0] * sh6[1][12]) + (sh1[0][2] * sh6[11][0] + sh1[0][0] * sh6[11][12]));
        sh7[1][1] = sh1[1][1] * sh6[0][0] + sqrt(3.0 / 1.0) * (sh1[2][1] * sh6[1][0] + sh1[0][1] * sh6[11][0]);
        sh7[1][2] = sqrt(13.0 / 24.0) * sh1[1][1] * sh6[0][1] + sqrt(13.0 / 8.0) * (sh1[2][1] * sh6[1][1] + sh1[0][1] * sh6[11][1]);
        sh7[1][3] = sqrt(13.0 / 33.0) * sh1[1][1] * sh6[0][2] + sqrt(13.0 / 11.0) * (sh1[2][1] * sh6[1][2] + sh1[0][1] * sh6[11][2]);
        sh7[1][4] = sqrt(13.0 / 40.0) * sh1[1][1] * sh6[0][3] + sqrt(39.0 / 40.0) * (sh1[2][1] * sh6[1][3] + sh1[0][1] * sh6[11][3]);
        sh7[1][5] = sqrt(13.0 / 45.0) * sh1[1][1] * sh6[0][4] + sqrt(13.0 / 15.0) * (sh1[2][1] * sh6[1][4] + sh1[0][1] * sh6[11][4]);
        sh7[1][6] = sqrt(13.0 / 48.0) * sh1[1][1] * sh6[0][5] + sqrt(13.0 / 16.0) * (sh1[2][1] * sh6[1][5] + sh1[0][1] * sh6[11][5]);
        sh7[1][7] = sqrt(13.0 / 49.0) * sh1[1][1] * sh6[0][6] + sqrt(39.0 / 49.0) * (sh1[2][1] * sh6[1][6] + sh1[0][1] * sh6[11][6]);
        sh7[1][8] = sqrt(13.0 / 48.0) * sh1[1][1] * sh6[0][7] + sqrt(13.0 / 16.0) * (sh1[2][1] * sh6[1][7] + sh1[0][1] * sh6[11][7]);
        sh7[1][9] = sqrt(13.0 / 45.0) * sh1[1][1] * sh6[0][8] + sqrt(13.0 / 15.0) * (sh1[2][1] * sh6[1][8] + sh1[0][1] * sh6[11][8]);
        sh7[1][10] = sqrt(13.0 / 40.0) * sh1[1][1] * sh6[0][9] + sqrt(39.0 / 40.0) * (sh1[2][1] * sh6[1][9] + sh1[0][1] * sh6[11][9]);
        sh7[1][11] = sqrt(13.0 / 33.0) * sh1[1][1] * sh6[0][10] + sqrt(13.0 / 11.0) * (sh1[2][1] * sh6[1][10] + sh1[0][1] * sh6[11][10]);
        sh7[1][12] = sqrt(13.0 / 24.0) * sh1[1][1] * sh6[0][11] + sqrt(13.0 / 8.0) * (sh1[2][1] * sh6[1][11] + sh1[0][1] * sh6[11][11]);
        sh7[1][13] = sh1[1][1] * sh6[0][12] + sqrt(3.0 / 1.0) * (sh1[2][1] * sh6[1][12] + sh1[0][1] * sh6[11][12]);
        sh7[1][14] = sqrt(1.0 / 14.0) * (sh1[1][2] * sh6[0][12] - sh1[1][0] * sh6[0][0]) + sqrt(3.0 / 14.0) * ((sh1[2][2] * sh6[1][12] - sh1[2][0] * sh6[1][0]) + (sh1[0][2] * sh6[11][12] - sh1[0][0] * sh6[11][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[1]);

        sh7[2][0] = sqrt(12.0 / 91.0) * (sh1[1][2] * sh6[1][0] + sh1[1][0] * sh6[1][12]) + sqrt(33.0 / 182.0) * ((sh1[2][2] * sh6[2][0] + sh1[2][0] * sh6[2][12]) + (sh1[0][2] * sh6[10][0] + sh1[0][0] * sh6[10][12])) + -sqrt(1.0 / 364.0) * ((sh1[2][2] * sh6[0][0] + sh1[2][0] * sh6[0][12]) - (sh1[0][2] * sh6[12][0] + sh1[0][0] * sh6[12][12]));
        sh7[2][1] = sqrt(24.0 / 13.0) * sh1[1][1] * sh6[1][0] + sqrt(33.0 / 13.0) * (sh1[2][1] * sh6[2][0] + sh1[0][1] * sh6[10][0]) + -sqrt(1.0 / 26.0) * (sh1[2][1] * sh6[0][0] - sh1[0][1] * sh6[12][0]);
        sh7[2][2] = sh1[1][1] * sh6[1][1] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh6[2][1] + sh1[0][1] * sh6[10][1]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh6[0][1] - sh1[0][1] * sh6[12][1]);
        sh7[2][3] = sqrt(8.0 / 11.0) * sh1[1][1] * sh6[1][2] + (sh1[2][1] * sh6[2][2] + sh1[0][1] * sh6[10][2]) + -sqrt(1.0 / 66.0) * (sh1[2][1] * sh6[0][2] - sh1[0][1] * sh6[12][2]);
        sh7[2][4] = kSqrt03_05 * sh1[1][1] * sh6[1][3] + sqrt(33.0 / 40.0) * (sh1[2][1] * sh6[2][3] + sh1[0][1] * sh6[10][3]) + -sqrt(1.0 / 80.0) * (sh1[2][1] * sh6[0][3] - sh1[0][1] * sh6[12][3]);
        sh7[2][5] = sqrt(8.0 / 15.0) * sh1[1][1] * sh6[1][4] + sqrt(11.0 / 15.0) * (sh1[2][1] * sh6[2][4] + sh1[0][1] * sh6[10][4]) + -sqrt(1.0 / 90.0) * (sh1[2][1] * sh6[0][4] - sh1[0][1] * sh6[12][4]);
        sh7[2][6] = sqrt(1.0 / 2.0) * sh1[1][1] * sh6[1][5] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[2][5] + sh1[0][1] * sh6[10][5]) + -sqrt(1.0 / 96.0) * (sh1[2][1] * sh6[0][5] - sh1[0][1] * sh6[12][5]);
        sh7[2][7] = sqrt(24.0 / 49.0) * sh1[1][1] * sh6[1][6] + sqrt(33.0 / 49.0) * (sh1[2][1] * sh6[2][6] + sh1[0][1] * sh6[10][6]) + -sqrt(1.0 / 98.0) * (sh1[2][1] * sh6[0][6] - sh1[0][1] * sh6[12][6]);
        sh7[2][8] = sqrt(1.0 / 2.0) * sh1[1][1] * sh6[1][7] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[2][7] + sh1[0][1] * sh6[10][7]) + -sqrt(1.0 / 96.0) * (sh1[2][1] * sh6[0][7] - sh1[0][1] * sh6[12][7]);
        sh7[2][9] = sqrt(8.0 / 15.0) * sh1[1][1] * sh6[1][8] + sqrt(11.0 / 15.0) * (sh1[2][1] * sh6[2][8] + sh1[0][1] * sh6[10][8]) + -sqrt(1.0 / 90.0) * (sh1[2][1] * sh6[0][8] - sh1[0][1] * sh6[12][8]);
        sh7[2][10] = kSqrt03_05 * sh1[1][1] * sh6[1][9] + sqrt(33.0 / 40.0) * (sh1[2][1] * sh6[2][9] + sh1[0][1] * sh6[10][9]) + -sqrt(1.0 / 80.0) * (sh1[2][1] * sh6[0][9] - sh1[0][1] * sh6[12][9]);
        sh7[2][11] = sqrt(8.0 / 11.0) * sh1[1][1] * sh6[1][10] + (sh1[2][1] * sh6[2][10] + sh1[0][1] * sh6[10][10]) + -sqrt(1.0 / 66.0) * (sh1[2][1] * sh6[0][10] - sh1[0][1] * sh6[12][10]);
        sh7[2][12] = sh1[1][1] * sh6[1][11] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh6[2][11] + sh1[0][1] * sh6[10][11]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh6[0][11] - sh1[0][1] * sh6[12][11]);
        sh7[2][13] = sqrt(24.0 / 13.0) * sh1[1][1] * sh6[1][12] + sqrt(33.0 / 13.0) * (sh1[2][1] * sh6[2][12] + sh1[0][1] * sh6[10][12]) + -sqrt(1.0 / 26.0) * (sh1[2][1] * sh6[0][12] - sh1[0][1] * sh6[12][12]);
        sh7[2][14] = sqrt(12.0 / 91.0) * (sh1[1][2] * sh6[1][12] - sh1[1][0] * sh6[1][0]) + sqrt(33.0 / 182.0) * ((sh1[2][2] * sh6[2][12] - sh1[2][0] * sh6[2][0]) + (sh1[0][2] * sh6[10][12] - sh1[0][0] * sh6[10][0])) + -sqrt(1.0 / 364.0) * ((sh1[2][2] * sh6[0][12] - sh1[2][0] * sh6[0][0]) - (sh1[0][2] * sh6[12][12] - sh1[0][0] * sh6[12][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[2]);

        sh7[3][0] = sqrt(33.0 / 182.0) * (sh1[1][2] * sh6[2][0] + sh1[1][0] * sh6[2][12]) + sqrt(55.0 / 364.0) * ((sh1[2][2] * sh6[3][0] + sh1[2][0] * sh6[3][12]) + (sh1[0][2] * sh6[9][0] + sh1[0][0] * sh6[9][12])) + -sqrt(3.0 / 364.0) * ((sh1[2][2] * sh6[1][0] + sh1[2][0] * sh6[1][12]) - (sh1[0][2] * sh6[11][0] + sh1[0][0] * sh6[11][12]));
        sh7[3][1] = sqrt(33.0 / 13.0) * sh1[1][1] * sh6[2][0] + sqrt(55.0 / 26.0) * (sh1[2][1] * sh6[3][0] + sh1[0][1] * sh6[9][0]) + -sqrt(3.0 / 26.0) * (sh1[2][1] * sh6[1][0] - sh1[0][1] * sh6[11][0]);
        sh7[3][2] = sqrt(11.0 / 8.0) * sh1[1][1] * sh6[2][1] + sqrt(55.0 / 48.0) * (sh1[2][1] * sh6[3][1] + sh1[0][1] * sh6[9][1]) + -kSqrt01_16 * (sh1[2][1] * sh6[1][1] - sh1[0][1] * sh6[11][1]);
        sh7[3][3] = sh1[1][1] * sh6[2][2] + kSqrt05_06 * (sh1[2][1] * sh6[3][2] + sh1[0][1] * sh6[9][2]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh6[1][2] - sh1[0][1] * sh6[11][2]);
        sh7[3][4] = sqrt(33.0 / 40.0) * sh1[1][1] * sh6[2][3] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[3][3] + sh1[0][1] * sh6[9][3]) + -sqrt(3.0 / 80.0) * (sh1[2][1] * sh6[1][3] - sh1[0][1] * sh6[11][3]);
        sh7[3][5] = sqrt(11.0 / 15.0) * sh1[1][1] * sh6[2][4] + sqrt(11.0 / 18.0) * (sh1[2][1] * sh6[3][4] + sh1[0][1] * sh6[9][4]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh6[1][4] - sh1[0][1] * sh6[11][4]);
        sh7[3][6] = sqrt(11.0 / 16.0) * sh1[1][1] * sh6[2][5] + sqrt(55.0 / 96.0) * (sh1[2][1] * sh6[3][5] + sh1[0][1] * sh6[9][5]) + -kSqrt01_32 * (sh1[2][1] * sh6[1][5] - sh1[0][1] * sh6[11][5]);
        sh7[3][7] = sqrt(33.0 / 49.0) * sh1[1][1] * sh6[2][6] + sqrt(55.0 / 98.0) * (sh1[2][1] * sh6[3][6] + sh1[0][1] * sh6[9][6]) + -sqrt(3.0 / 98.0) * (sh1[2][1] * sh6[1][6] - sh1[0][1] * sh6[11][6]);
        sh7[3][8] = sqrt(11.0 / 16.0) * sh1[1][1] * sh6[2][7] + sqrt(55.0 / 96.0) * (sh1[2][1] * sh6[3][7] + sh1[0][1] * sh6[9][7]) + -kSqrt01_32 * (sh1[2][1] * sh6[1][7] - sh1[0][1] * sh6[11][7]);
        sh7[3][9] = sqrt(11.0 / 15.0) * sh1[1][1] * sh6[2][8] + sqrt(11.0 / 18.0) * (sh1[2][1] * sh6[3][8] + sh1[0][1] * sh6[9][8]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh6[1][8] - sh1[0][1] * sh6[11][8]);
        sh7[3][10] = sqrt(33.0 / 40.0) * sh1[1][1] * sh6[2][9] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[3][9] + sh1[0][1] * sh6[9][9]) + -sqrt(3.0 / 80.0) * (sh1[2][1] * sh6[1][9] - sh1[0][1] * sh6[11][9]);
        sh7[3][11] = sh1[1][1] * sh6[2][10] + kSqrt05_06 * (sh1[2][1] * sh6[3][10] + sh1[0][1] * sh6[9][10]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh6[1][10] - sh1[0][1] * sh6[11][10]);
        sh7[3][12] = sqrt(11.0 / 8.0) * sh1[1][1] * sh6[2][11] + sqrt(55.0 / 48.0) * (sh1[2][1] * sh6[3][11] + sh1[0][1] * sh6[9][11]) + -kSqrt01_16 * (sh1[2][1] * sh6[1][11] - sh1[0][1] * sh6[11][11]);
        sh7[3][13] = sqrt(33.0 / 13.0) * sh1[1][1] * sh6[2][12] + sqrt(55.0 / 26.0) * (sh1[2][1] * sh6[3][12] + sh1[0][1] * sh6[9][12]) + -sqrt(3.0 / 26.0) * (sh1[2][1] * sh6[1][12] - sh1[0][1] * sh6[11][12]);
        sh7[3][14] = sqrt(33.0 / 182.0) * (sh1[1][2] * sh6[2][12] - sh1[1][0] * sh6[2][0]) + sqrt(55.0 / 364.0) * ((sh1[2][2] * sh6[3][12] - sh1[2][0] * sh6[3][0]) + (sh1[0][2] * sh6[9][12] - sh1[0][0] * sh6[9][0])) + -sqrt(3.0 / 364.0) * ((sh1[2][2] * sh6[1][12] - sh1[2][0] * sh6[1][0]) - (sh1[0][2] * sh6[11][12] - sh1[0][0] * sh6[11][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[3]);

        sh7[4][0] = sqrt(20.0 / 91.0) * (sh1[1][2] * sh6[3][0] + sh1[1][0] * sh6[3][12]) + sqrt(45.0 / 364.0) * ((sh1[2][2] * sh6[4][0] + sh1[2][0] * sh6[4][12]) + (sh1[0][2] * sh6[8][0] + sh1[0][0] * sh6[8][12])) + -sqrt(3.0 / 182.0) * ((sh1[2][2] * sh6[2][0] + sh1[2][0] * sh6[2][12]) - (sh1[0][2] * sh6[10][0] + sh1[0][0] * sh6[10][12]));
        sh7[4][1] = sqrt(40.0 / 13.0) * sh1[1][1] * sh6[3][0] + sqrt(45.0 / 26.0) * (sh1[2][1] * sh6[4][0] + sh1[0][1] * sh6[8][0]) + -sqrt(3.0 / 13.0) * (sh1[2][1] * sh6[2][0] - sh1[0][1] * sh6[10][0]);
        sh7[4][2] = sqrt(5.0 / 3.0) * sh1[1][1] * sh6[3][1] + kSqrt15_16 * (sh1[2][1] * sh6[4][1] + sh1[0][1] * sh6[8][1]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[2][1] - sh1[0][1] * sh6[10][1]);
        sh7[4][3] = sqrt(40.0 / 33.0) * sh1[1][1] * sh6[3][2] + sqrt(15.0 / 22.0) * (sh1[2][1] * sh6[4][2] + sh1[0][1] * sh6[8][2]) + -sqrt(1.0 / 11.0) * (sh1[2][1] * sh6[2][2] - sh1[0][1] * sh6[10][2]);
        sh7[4][4] = sh1[1][1] * sh6[3][3] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh6[4][3] + sh1[0][1] * sh6[8][3]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh6[2][3] - sh1[0][1] * sh6[10][3]);
        sh7[4][5] = kSqrt08_09 * sh1[1][1] * sh6[3][4] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh6[4][4] + sh1[0][1] * sh6[8][4]) + -sqrt(1.0 / 15.0) * (sh1[2][1] * sh6[2][4] - sh1[0][1] * sh6[10][4]);
        sh7[4][6] = kSqrt05_06 * sh1[1][1] * sh6[3][5] + kSqrt15_32 * (sh1[2][1] * sh6[4][5] + sh1[0][1] * sh6[8][5]) + -kSqrt01_16 * (sh1[2][1] * sh6[2][5] - sh1[0][1] * sh6[10][5]);
        sh7[4][7] = sqrt(40.0 / 49.0) * sh1[1][1] * sh6[3][6] + sqrt(45.0 / 98.0) * (sh1[2][1] * sh6[4][6] + sh1[0][1] * sh6[8][6]) + -sqrt(3.0 / 49.0) * (sh1[2][1] * sh6[2][6] - sh1[0][1] * sh6[10][6]);
        sh7[4][8] = kSqrt05_06 * sh1[1][1] * sh6[3][7] + kSqrt15_32 * (sh1[2][1] * sh6[4][7] + sh1[0][1] * sh6[8][7]) + -kSqrt01_16 * (sh1[2][1] * sh6[2][7] - sh1[0][1] * sh6[10][7]);
        sh7[4][9] = kSqrt08_09 * sh1[1][1] * sh6[3][8] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh6[4][8] + sh1[0][1] * sh6[8][8]) + -sqrt(1.0 / 15.0) * (sh1[2][1] * sh6[2][8] - sh1[0][1] * sh6[10][8]);
        sh7[4][10] = sh1[1][1] * sh6[3][9] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh6[4][9] + sh1[0][1] * sh6[8][9]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh6[2][9] - sh1[0][1] * sh6[10][9]);
        sh7[4][11] = sqrt(40.0 / 33.0) * sh1[1][1] * sh6[3][10] + sqrt(15.0 / 22.0) * (sh1[2][1] * sh6[4][10] + sh1[0][1] * sh6[8][10]) + -sqrt(1.0 / 11.0) * (sh1[2][1] * sh6[2][10] - sh1[0][1] * sh6[10][10]);
        sh7[4][12] = sqrt(5.0 / 3.0) * sh1[1][1] * sh6[3][11] + kSqrt15_16 * (sh1[2][1] * sh6[4][11] + sh1[0][1] * sh6[8][11]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[2][11] - sh1[0][1] * sh6[10][11]);
        sh7[4][13] = sqrt(40.0 / 13.0) * sh1[1][1] * sh6[3][12] + sqrt(45.0 / 26.0) * (sh1[2][1] * sh6[4][12] + sh1[0][1] * sh6[8][12]) + -sqrt(3.0 / 13.0) * (sh1[2][1] * sh6[2][12] - sh1[0][1] * sh6[10][12]);
        sh7[4][14] = sqrt(20.0 / 91.0) * (sh1[1][2] * sh6[3][12] - sh1[1][0] * sh6[3][0]) + sqrt(45.0 / 364.0) * ((sh1[2][2] * sh6[4][12] - sh1[2][0] * sh6[4][0]) + (sh1[0][2] * sh6[8][12] - sh1[0][0] * sh6[8][0])) + -sqrt(3.0 / 182.0) * ((sh1[2][2] * sh6[2][12] - sh1[2][0] * sh6[2][0]) - (sh1[0][2] * sh6[10][12] - sh1[0][0] * sh6[10][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[4]);

        sh7[5][0] = sqrt(45.0 / 182.0) * (sh1[1][2] * sh6[4][0] + sh1[1][0] * sh6[4][12]) + sqrt(9.0 / 91.0) * ((sh1[2][2] * sh6[5][0] + sh1[2][0] * sh6[5][12]) + (sh1[0][2] * sh6[7][0] + sh1[0][0] * sh6[7][12])) + -sqrt(5.0 / 182.0) * ((sh1[2][2] * sh6[3][0] + sh1[2][0] * sh6[3][12]) - (sh1[0][2] * sh6[9][0] + sh1[0][0] * sh6[9][12]));
        sh7[5][1] = sqrt(45.0 / 13.0) * sh1[1][1] * sh6[4][0] + sqrt(18.0 / 13.0) * (sh1[2][1] * sh6[5][0] + sh1[0][1] * sh6[7][0]) + -sqrt(5.0 / 13.0) * (sh1[2][1] * sh6[3][0] - sh1[0][1] * sh6[9][0]);
        sh7[5][2] = sqrt(15.0 / 8.0) * sh1[1][1] * sh6[4][1] + kSqrt03_04 * (sh1[2][1] * sh6[5][1] + sh1[0][1] * sh6[7][1]) + -sqrt(5.0 / 24.0) * (sh1[2][1] * sh6[3][1] - sh1[0][1] * sh6[9][1]);
        sh7[5][3] = sqrt(15.0 / 11.0) * sh1[1][1] * sh6[4][2] + sqrt(6.0 / 11.0) * (sh1[2][1] * sh6[5][2] + sh1[0][1] * sh6[7][2]) + -sqrt(5.0 / 33.0) * (sh1[2][1] * sh6[3][2] - sh1[0][1] * sh6[9][2]);
        sh7[5][4] = kSqrt09_08 * sh1[1][1] * sh6[4][3] + sqrt(9.0 / 20.0) * (sh1[2][1] * sh6[5][3] + sh1[0][1] * sh6[7][3]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[3][3] - sh1[0][1] * sh6[9][3]);
        sh7[5][5] = sh1[1][1] * sh6[4][4] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh6[5][4] + sh1[0][1] * sh6[7][4]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh6[3][4] - sh1[0][1] * sh6[9][4]);
        sh7[5][6] = kSqrt15_16 * sh1[1][1] * sh6[4][5] + kSqrt03_08 * (sh1[2][1] * sh6[5][5] + sh1[0][1] * sh6[7][5]) + -sqrt(5.0 / 48.0) * (sh1[2][1] * sh6[3][5] - sh1[0][1] * sh6[9][5]);
        sh7[5][7] = sqrt(45.0 / 49.0) * sh1[1][1] * sh6[4][6] + sqrt(18.0 / 49.0) * (sh1[2][1] * sh6[5][6] + sh1[0][1] * sh6[7][6]) + -sqrt(5.0 / 49.0) * (sh1[2][1] * sh6[3][6] - sh1[0][1] * sh6[9][6]);
        sh7[5][8] = kSqrt15_16 * sh1[1][1] * sh6[4][7] + kSqrt03_08 * (sh1[2][1] * sh6[5][7] + sh1[0][1] * sh6[7][7]) + -sqrt(5.0 / 48.0) * (sh1[2][1] * sh6[3][7] - sh1[0][1] * sh6[9][7]);
        sh7[5][9] = sh1[1][1] * sh6[4][8] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh6[5][8] + sh1[0][1] * sh6[7][8]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh6[3][8] - sh1[0][1] * sh6[9][8]);
        sh7[5][10] = kSqrt09_08 * sh1[1][1] * sh6[4][9] + sqrt(9.0 / 20.0) * (sh1[2][1] * sh6[5][9] + sh1[0][1] * sh6[7][9]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[3][9] - sh1[0][1] * sh6[9][9]);
        sh7[5][11] = sqrt(15.0 / 11.0) * sh1[1][1] * sh6[4][10] + sqrt(6.0 / 11.0) * (sh1[2][1] * sh6[5][10] + sh1[0][1] * sh6[7][10]) + -sqrt(5.0 / 33.0) * (sh1[2][1] * sh6[3][10] - sh1[0][1] * sh6[9][10]);
        sh7[5][12] = sqrt(15.0 / 8.0) * sh1[1][1] * sh6[4][11] + kSqrt03_04 * (sh1[2][1] * sh6[5][11] + sh1[0][1] * sh6[7][11]) + -sqrt(5.0 / 24.0) * (sh1[2][1] * sh6[3][11] - sh1[0][1] * sh6[9][11]);
        sh7[5][13] = sqrt(45.0 / 13.0) * sh1[1][1] * sh6[4][12] + sqrt(18.0 / 13.0) * (sh1[2][1] * sh6[5][12] + sh1[0][1] * sh6[7][12]) + -sqrt(5.0 / 13.0) * (sh1[2][1] * sh6[3][12] - sh1[0][1] * sh6[9][12]);
        sh7[5][14] = sqrt(45.0 / 182.0) * (sh1[1][2] * sh6[4][12] - sh1[1][0] * sh6[4][0]) + sqrt(9.0 / 91.0) * ((sh1[2][2] * sh6[5][12] - sh1[2][0] * sh6[5][0]) + (sh1[0][2] * sh6[7][12] - sh1[0][0] * sh6[7][0])) + -sqrt(5.0 / 182.0) * ((sh1[2][2] * sh6[3][12] - sh1[2][0] * sh6[3][0]) - (sh1[0][2] * sh6[9][12] - sh1[0][0] * sh6[9][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[5]);

        sh7[6][0] = sqrt(24.0 / 91.0) * (sh1[1][2] * sh6[5][0] + sh1[1][0] * sh6[5][12]) + sqrt(2.0 / 13.0) * (sh1[0][2] * sh6[6][0] + sh1[0][0] * sh6[6][12]) + -sqrt(15.0 / 364.0) * ((sh1[2][2] * sh6[4][0] + sh1[2][0] * sh6[4][12]) - (sh1[0][2] * sh6[8][0] + sh1[0][0] * sh6[8][12]));
        sh7[6][1] = sqrt(48.0 / 13.0) * sh1[1][1] * sh6[5][0] + sqrt(28.0 / 13.0) * sh1[0][1] * sh6[6][0] + -sqrt(15.0 / 26.0) * (sh1[2][1] * sh6[4][0] - sh1[0][1] * sh6[8][0]);
        sh7[6][2] = sqrt(2.0 / 1.0) * sh1[1][1] * sh6[5][1] + sqrt(7.0 / 6.0) * sh1[0][1] * sh6[6][1] + -sqrt(5.0 / 16.0) * (sh1[2][1] * sh6[4][1] - sh1[0][1] * sh6[8][1]);
        sh7[6][3] = sqrt(16.0 / 11.0) * sh1[1][1] * sh6[5][2] + sqrt(28.0 / 33.0) * sh1[0][1] * sh6[6][2] + -sqrt(5.0 / 22.0) * (sh1[2][1] * sh6[4][2] - sh1[0][1] * sh6[8][2]);
        sh7[6][4] = kSqrt06_05 * sh1[1][1] * sh6[5][3] + sqrt(7.0 / 10.0) * sh1[0][1] * sh6[6][3] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh6[4][3] - sh1[0][1] * sh6[8][3]);
        sh7[6][5] = sqrt(16.0 / 15.0) * sh1[1][1] * sh6[5][4] + sqrt(28.0 / 45.0) * sh1[0][1] * sh6[6][4] + -kSqrt01_06 * (sh1[2][1] * sh6[4][4] - sh1[0][1] * sh6[8][4]);
        sh7[6][6] = sh1[1][1] * sh6[5][5] + sqrt(7.0 / 12.0) * sh1[0][1] * sh6[6][5] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh6[4][5] - sh1[0][1] * sh6[8][5]);
        sh7[6][7] = sqrt(48.0 / 49.0) * sh1[1][1] * sh6[5][6] + sqrt(4.0 / 7.0) * sh1[0][1] * sh6[6][6] + -sqrt(15.0 / 98.0) * (sh1[2][1] * sh6[4][6] - sh1[0][1] * sh6[8][6]);
        sh7[6][8] = sh1[1][1] * sh6[5][7] + sqrt(7.0 / 12.0) * sh1[0][1] * sh6[6][7] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh6[4][7] - sh1[0][1] * sh6[8][7]);
        sh7[6][9] = sqrt(16.0 / 15.0) * sh1[1][1] * sh6[5][8] + sqrt(28.0 / 45.0) * sh1[0][1] * sh6[6][8] + -kSqrt01_06 * (sh1[2][1] * sh6[4][8] - sh1[0][1] * sh6[8][8]);
        sh7[6][10] = kSqrt06_05 * sh1[1][1] * sh6[5][9] + sqrt(7.0 / 10.0) * sh1[0][1] * sh6[6][9] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh6[4][9] - sh1[0][1] * sh6[8][9]);
        sh7[6][11] = sqrt(16.0 / 11.0) * sh1[1][1] * sh6[5][10] + sqrt(28.0 / 33.0) * sh1[0][1] * sh6[6][10] + -sqrt(5.0 / 22.0) * (sh1[2][1] * sh6[4][10] - sh1[0][1] * sh6[8][10]);
        sh7[6][12] = sqrt(2.0 / 1.0) * sh1[1][1] * sh6[5][11] + sqrt(7.0 / 6.0) * sh1[0][1] * sh6[6][11] + -sqrt(5.0 / 16.0) * (sh1[2][1] * sh6[4][11] - sh1[0][1] * sh6[8][11]);
        sh7[6][13] = sqrt(48.0 / 13.0) * sh1[1][1] * sh6[5][12] + sqrt(28.0 / 13.0) * sh1[0][1] * sh6[6][12] + -sqrt(15.0 / 26.0) * (sh1[2][1] * sh6[4][12] - sh1[0][1] * sh6[8][12]);
        sh7[6][14] = sqrt(24.0 / 91.0) * (sh1[1][2] * sh6[5][12] - sh1[1][0] * sh6[5][0]) + sqrt(2.0 / 13.0) * (sh1[0][2] * sh6[6][12] - sh1[0][0] * sh6[6][0]) + -sqrt(15.0 / 364.0) * ((sh1[2][2] * sh6[4][12] - sh1[2][0] * sh6[4][0]) - (sh1[0][2] * sh6[8][12] - sh1[0][0] * sh6[8][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[6]);

        sh7[7][0] = sqrt(7.0 / 26.0) * (sh1[1][2] * sh6[6][0] + sh1[1][0] * sh6[6][12]) + -sqrt(3.0 / 26.0) * ((sh1[2][2] * sh6[7][0] + sh1[2][0] * sh6[7][12]) + (sh1[0][2] * sh6[5][0] + sh1[0][0] * sh6[5][12]));
        sh7[7][1] = sqrt(49.0 / 13.0) * sh1[1][1] * sh6[6][0] + -sqrt(21.0 / 13.0) * (sh1[2][1] * sh6[7][0] + sh1[0][1] * sh6[5][0]);
        sh7[7][2] = sqrt(49.0 / 24.0) * sh1[1][1] * sh6[6][1] + -kSqrt07_08 * (sh1[2][1] * sh6[7][1] + sh1[0][1] * sh6[5][1]);
        sh7[7][3] = sqrt(49.0 / 33.0) * sh1[1][1] * sh6[6][2] + -sqrt(7.0 / 11.0) * (sh1[2][1] * sh6[7][2] + sh1[0][1] * sh6[5][2]);
        sh7[7][4] = sqrt(49.0 / 40.0) * sh1[1][1] * sh6[6][3] + -sqrt(21.0 / 40.0) * (sh1[2][1] * sh6[7][3] + sh1[0][1] * sh6[5][3]);
        sh7[7][5] = sqrt(49.0 / 45.0) * sh1[1][1] * sh6[6][4] + -sqrt(7.0 / 15.0) * (sh1[2][1] * sh6[7][4] + sh1[0][1] * sh6[5][4]);
        sh7[7][6] = sqrt(49.0 / 48.0) * sh1[1][1] * sh6[6][5] + -kSqrt07_16 * (sh1[2][1] * sh6[7][5] + sh1[0][1] * sh6[5][5]);
        sh7[7][7] = sh1[1][1] * sh6[6][6] + -sqrt(3.0 / 7.0) * (sh1[2][1] * sh6[7][6] + sh1[0][1] * sh6[5][6]);
        sh7[7][8] = sqrt(49.0 / 48.0) * sh1[1][1] * sh6[6][7] + -kSqrt07_16 * (sh1[2][1] * sh6[7][7] + sh1[0][1] * sh6[5][7]);
        sh7[7][9] = sqrt(49.0 / 45.0) * sh1[1][1] * sh6[6][8] + -sqrt(7.0 / 15.0) * (sh1[2][1] * sh6[7][8] + sh1[0][1] * sh6[5][8]);
        sh7[7][10] = sqrt(49.0 / 40.0) * sh1[1][1] * sh6[6][9] + -sqrt(21.0 / 40.0) * (sh1[2][1] * sh6[7][9] + sh1[0][1] * sh6[5][9]);
        sh7[7][11] = sqrt(49.0 / 33.0) * sh1[1][1] * sh6[6][10] + -sqrt(7.0 / 11.0) * (sh1[2][1] * sh6[7][10] + sh1[0][1] * sh6[5][10]);
        sh7[7][12] = sqrt(49.0 / 24.0) * sh1[1][1] * sh6[6][11] + -kSqrt07_08 * (sh1[2][1] * sh6[7][11] + sh1[0][1] * sh6[5][11]);
        sh7[7][13] = sqrt(49.0 / 13.0) * sh1[1][1] * sh6[6][12] + -sqrt(21.0 / 13.0) * (sh1[2][1] * sh6[7][12] + sh1[0][1] * sh6[5][12]);
        sh7[7][14] = sqrt(7.0 / 26.0) * (sh1[1][2] * sh6[6][12] - sh1[1][0] * sh6[6][0]) + -sqrt(3.0 / 26.0) * ((sh1[2][2] * sh6[7][12] - sh1[2][0] * sh6[7][0]) + (sh1[0][2] * sh6[5][12] - sh1[0][0] * sh6[5][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[7]);

        sh7[8][0] = sqrt(24.0 / 91.0) * (sh1[1][2] * sh6[7][0] + sh1[1][0] * sh6[7][12]) + sqrt(2.0 / 13.0) * (sh1[2][2] * sh6[6][0] + sh1[2][0] * sh6[6][12]) + -sqrt(15.0 / 364.0) * ((sh1[2][2] * sh6[8][0] + sh1[2][0] * sh6[8][12]) + (sh1[0][2] * sh6[4][0] + sh1[0][0] * sh6[4][12]));
        sh7[8][1] = sqrt(48.0 / 13.0) * sh1[1][1] * sh6[7][0] + sqrt(28.0 / 13.0) * sh1[2][1] * sh6[6][0] + -sqrt(15.0 / 26.0) * (sh1[2][1] * sh6[8][0] + sh1[0][1] * sh6[4][0]);
        sh7[8][2] = sqrt(2.0 / 1.0) * sh1[1][1] * sh6[7][1] + sqrt(7.0 / 6.0) * sh1[2][1] * sh6[6][1] + -sqrt(5.0 / 16.0) * (sh1[2][1] * sh6[8][1] + sh1[0][1] * sh6[4][1]);
        sh7[8][3] = sqrt(16.0 / 11.0) * sh1[1][1] * sh6[7][2] + sqrt(28.0 / 33.0) * sh1[2][1] * sh6[6][2] + -sqrt(5.0 / 22.0) * (sh1[2][1] * sh6[8][2] + sh1[0][1] * sh6[4][2]);
        sh7[8][4] = kSqrt06_05 * sh1[1][1] * sh6[7][3] + sqrt(7.0 / 10.0) * sh1[2][1] * sh6[6][3] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh6[8][3] + sh1[0][1] * sh6[4][3]);
        sh7[8][5] = sqrt(16.0 / 15.0) * sh1[1][1] * sh6[7][4] + sqrt(28.0 / 45.0) * sh1[2][1] * sh6[6][4] + -kSqrt01_06 * (sh1[2][1] * sh6[8][4] + sh1[0][1] * sh6[4][4]);
        sh7[8][6] = sh1[1][1] * sh6[7][5] + sqrt(7.0 / 12.0) * sh1[2][1] * sh6[6][5] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh6[8][5] + sh1[0][1] * sh6[4][5]);
        sh7[8][7] = sqrt(48.0 / 49.0) * sh1[1][1] * sh6[7][6] + sqrt(4.0 / 7.0) * sh1[2][1] * sh6[6][6] + -sqrt(15.0 / 98.0) * (sh1[2][1] * sh6[8][6] + sh1[0][1] * sh6[4][6]);
        sh7[8][8] = sh1[1][1] * sh6[7][7] + sqrt(7.0 / 12.0) * sh1[2][1] * sh6[6][7] + -sqrt(5.0 / 32.0) * (sh1[2][1] * sh6[8][7] + sh1[0][1] * sh6[4][7]);
        sh7[8][9] = sqrt(16.0 / 15.0) * sh1[1][1] * sh6[7][8] + sqrt(28.0 / 45.0) * sh1[2][1] * sh6[6][8] + -kSqrt01_06 * (sh1[2][1] * sh6[8][8] + sh1[0][1] * sh6[4][8]);
        sh7[8][10] = kSqrt06_05 * sh1[1][1] * sh6[7][9] + sqrt(7.0 / 10.0) * sh1[2][1] * sh6[6][9] + -sqrt(3.0 / 16.0) * (sh1[2][1] * sh6[8][9] + sh1[0][1] * sh6[4][9]);
        sh7[8][11] = sqrt(16.0 / 11.0) * sh1[1][1] * sh6[7][10] + sqrt(28.0 / 33.0) * sh1[2][1] * sh6[6][10] + -sqrt(5.0 / 22.0) * (sh1[2][1] * sh6[8][10] + sh1[0][1] * sh6[4][10]);
        sh7[8][12] = sqrt(2.0 / 1.0) * sh1[1][1] * sh6[7][11] + sqrt(7.0 / 6.0) * sh1[2][1] * sh6[6][11] + -sqrt(5.0 / 16.0) * (sh1[2][1] * sh6[8][11] + sh1[0][1] * sh6[4][11]);
        sh7[8][13] = sqrt(48.0 / 13.0) * sh1[1][1] * sh6[7][12] + sqrt(28.0 / 13.0) * sh1[2][1] * sh6[6][12] + -sqrt(15.0 / 26.0) * (sh1[2][1] * sh6[8][12] + sh1[0][1] * sh6[4][12]);
        sh7[8][14] = sqrt(24.0 / 91.0) * (sh1[1][2] * sh6[7][12] - sh1[1][0] * sh6[7][0]) + sqrt(2.0 / 13.0) * (sh1[2][2] * sh6[6][12] - sh1[2][0] * sh6[6][0]) + -sqrt(15.0 / 364.0) * ((sh1[2][2] * sh6[8][12] - sh1[2][0] * sh6[8][0]) + (sh1[0][2] * sh6[4][12] - sh1[0][0] * sh6[4][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[8]);

        sh7[9][0] = sqrt(45.0 / 182.0) * (sh1[1][2] * sh6[8][0] + sh1[1][0] * sh6[8][12]) + sqrt(9.0 / 91.0) * ((sh1[2][2] * sh6[7][0] + sh1[2][0] * sh6[7][12]) - (sh1[0][2] * sh6[5][0] + sh1[0][0] * sh6[5][12])) + -sqrt(5.0 / 182.0) * ((sh1[2][2] * sh6[9][0] + sh1[2][0] * sh6[9][12]) + (sh1[0][2] * sh6[3][0] + sh1[0][0] * sh6[3][12]));
        sh7[9][1] = sqrt(45.0 / 13.0) * sh1[1][1] * sh6[8][0] + sqrt(18.0 / 13.0) * (sh1[2][1] * sh6[7][0] - sh1[0][1] * sh6[5][0]) + -sqrt(5.0 / 13.0) * (sh1[2][1] * sh6[9][0] + sh1[0][1] * sh6[3][0]);
        sh7[9][2] = sqrt(15.0 / 8.0) * sh1[1][1] * sh6[8][1] + kSqrt03_04 * (sh1[2][1] * sh6[7][1] - sh1[0][1] * sh6[5][1]) + -sqrt(5.0 / 24.0) * (sh1[2][1] * sh6[9][1] + sh1[0][1] * sh6[3][1]);
        sh7[9][3] = sqrt(15.0 / 11.0) * sh1[1][1] * sh6[8][2] + sqrt(6.0 / 11.0) * (sh1[2][1] * sh6[7][2] - sh1[0][1] * sh6[5][2]) + -sqrt(5.0 / 33.0) * (sh1[2][1] * sh6[9][2] + sh1[0][1] * sh6[3][2]);
        sh7[9][4] = kSqrt09_08 * sh1[1][1] * sh6[8][3] + sqrt(9.0 / 20.0) * (sh1[2][1] * sh6[7][3] - sh1[0][1] * sh6[5][3]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[9][3] + sh1[0][1] * sh6[3][3]);
        sh7[9][5] = sh1[1][1] * sh6[8][4] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh6[7][4] - sh1[0][1] * sh6[5][4]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh6[9][4] + sh1[0][1] * sh6[3][4]);
        sh7[9][6] = kSqrt15_16 * sh1[1][1] * sh6[8][5] + kSqrt03_08 * (sh1[2][1] * sh6[7][5] - sh1[0][1] * sh6[5][5]) + -sqrt(5.0 / 48.0) * (sh1[2][1] * sh6[9][5] + sh1[0][1] * sh6[3][5]);
        sh7[9][7] = sqrt(45.0 / 49.0) * sh1[1][1] * sh6[8][6] + sqrt(18.0 / 49.0) * (sh1[2][1] * sh6[7][6] - sh1[0][1] * sh6[5][6]) + -sqrt(5.0 / 49.0) * (sh1[2][1] * sh6[9][6] + sh1[0][1] * sh6[3][6]);
        sh7[9][8] = kSqrt15_16 * sh1[1][1] * sh6[8][7] + kSqrt03_08 * (sh1[2][1] * sh6[7][7] - sh1[0][1] * sh6[5][7]) + -sqrt(5.0 / 48.0) * (sh1[2][1] * sh6[9][7] + sh1[0][1] * sh6[3][7]);
        sh7[9][9] = sh1[1][1] * sh6[8][8] + sqrt(2.0 / 5.0) * (sh1[2][1] * sh6[7][8] - sh1[0][1] * sh6[5][8]) + -sqrt(1.0 / 9.0) * (sh1[2][1] * sh6[9][8] + sh1[0][1] * sh6[3][8]);
        sh7[9][10] = kSqrt09_08 * sh1[1][1] * sh6[8][9] + sqrt(9.0 / 20.0) * (sh1[2][1] * sh6[7][9] - sh1[0][1] * sh6[5][9]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[9][9] + sh1[0][1] * sh6[3][9]);
        sh7[9][11] = sqrt(15.0 / 11.0) * sh1[1][1] * sh6[8][10] + sqrt(6.0 / 11.0) * (sh1[2][1] * sh6[7][10] - sh1[0][1] * sh6[5][10]) + -sqrt(5.0 / 33.0) * (sh1[2][1] * sh6[9][10] + sh1[0][1] * sh6[3][10]);
        sh7[9][12] = sqrt(15.0 / 8.0) * sh1[1][1] * sh6[8][11] + kSqrt03_04 * (sh1[2][1] * sh6[7][11] - sh1[0][1] * sh6[5][11]) + -sqrt(5.0 / 24.0) * (sh1[2][1] * sh6[9][11] + sh1[0][1] * sh6[3][11]);
        sh7[9][13] = sqrt(45.0 / 13.0) * sh1[1][1] * sh6[8][12] + sqrt(18.0 / 13.0) * (sh1[2][1] * sh6[7][12] - sh1[0][1] * sh6[5][12]) + -sqrt(5.0 / 13.0) * (sh1[2][1] * sh6[9][12] + sh1[0][1] * sh6[3][12]);
        sh7[9][14] = sqrt(45.0 / 182.0) * (sh1[1][2] * sh6[8][12] - sh1[1][0] * sh6[8][0]) + sqrt(9.0 / 91.0) * ((sh1[2][2] * sh6[7][12] - sh1[2][0] * sh6[7][0]) - (sh1[0][2] * sh6[5][12] - sh1[0][0] * sh6[5][0])) + -sqrt(5.0 / 182.0) * ((sh1[2][2] * sh6[9][12] - sh1[2][0] * sh6[9][0]) + (sh1[0][2] * sh6[3][12] - sh1[0][0] * sh6[3][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[9]);

        sh7[10][0] = sqrt(20.0 / 91.0) * (sh1[1][2] * sh6[9][0] + sh1[1][0] * sh6[9][12]) + sqrt(45.0 / 364.0) * ((sh1[2][2] * sh6[8][0] + sh1[2][0] * sh6[8][12]) - (sh1[0][2] * sh6[4][0] + sh1[0][0] * sh6[4][12])) + -sqrt(3.0 / 182.0) * ((sh1[2][2] * sh6[10][0] + sh1[2][0] * sh6[10][12]) + (sh1[0][2] * sh6[2][0] + sh1[0][0] * sh6[2][12]));
        sh7[10][1] = sqrt(40.0 / 13.0) * sh1[1][1] * sh6[9][0] + sqrt(45.0 / 26.0) * (sh1[2][1] * sh6[8][0] - sh1[0][1] * sh6[4][0]) + -sqrt(3.0 / 13.0) * (sh1[2][1] * sh6[10][0] + sh1[0][1] * sh6[2][0]);
        sh7[10][2] = sqrt(5.0 / 3.0) * sh1[1][1] * sh6[9][1] + kSqrt15_16 * (sh1[2][1] * sh6[8][1] - sh1[0][1] * sh6[4][1]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[10][1] + sh1[0][1] * sh6[2][1]);
        sh7[10][3] = sqrt(40.0 / 33.0) * sh1[1][1] * sh6[9][2] + sqrt(15.0 / 22.0) * (sh1[2][1] * sh6[8][2] - sh1[0][1] * sh6[4][2]) + -sqrt(1.0 / 11.0) * (sh1[2][1] * sh6[10][2] + sh1[0][1] * sh6[2][2]);
        sh7[10][4] = sh1[1][1] * sh6[9][3] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh6[8][3] - sh1[0][1] * sh6[4][3]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh6[10][3] + sh1[0][1] * sh6[2][3]);
        sh7[10][5] = kSqrt08_09 * sh1[1][1] * sh6[9][4] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh6[8][4] - sh1[0][1] * sh6[4][4]) + -sqrt(1.0 / 15.0) * (sh1[2][1] * sh6[10][4] + sh1[0][1] * sh6[2][4]);
        sh7[10][6] = kSqrt05_06 * sh1[1][1] * sh6[9][5] + kSqrt15_32 * (sh1[2][1] * sh6[8][5] - sh1[0][1] * sh6[4][5]) + -kSqrt01_16 * (sh1[2][1] * sh6[10][5] + sh1[0][1] * sh6[2][5]);
        sh7[10][7] = sqrt(40.0 / 49.0) * sh1[1][1] * sh6[9][6] + sqrt(45.0 / 98.0) * (sh1[2][1] * sh6[8][6] - sh1[0][1] * sh6[4][6]) + -sqrt(3.0 / 49.0) * (sh1[2][1] * sh6[10][6] + sh1[0][1] * sh6[2][6]);
        sh7[10][8] = kSqrt05_06 * sh1[1][1] * sh6[9][7] + kSqrt15_32 * (sh1[2][1] * sh6[8][7] - sh1[0][1] * sh6[4][7]) + -kSqrt01_16 * (sh1[2][1] * sh6[10][7] + sh1[0][1] * sh6[2][7]);
        sh7[10][9] = kSqrt08_09 * sh1[1][1] * sh6[9][8] + sqrt(1.0 / 2.0) * (sh1[2][1] * sh6[8][8] - sh1[0][1] * sh6[4][8]) + -sqrt(1.0 / 15.0) * (sh1[2][1] * sh6[10][8] + sh1[0][1] * sh6[2][8]);
        sh7[10][10] = sh1[1][1] * sh6[9][9] + sqrt(9.0 / 16.0) * (sh1[2][1] * sh6[8][9] - sh1[0][1] * sh6[4][9]) + -sqrt(3.0 / 40.0) * (sh1[2][1] * sh6[10][9] + sh1[0][1] * sh6[2][9]);
        sh7[10][11] = sqrt(40.0 / 33.0) * sh1[1][1] * sh6[9][10] + sqrt(15.0 / 22.0) * (sh1[2][1] * sh6[8][10] - sh1[0][1] * sh6[4][10]) + -sqrt(1.0 / 11.0) * (sh1[2][1] * sh6[10][10] + sh1[0][1] * sh6[2][10]);
        sh7[10][12] = sqrt(5.0 / 3.0) * sh1[1][1] * sh6[9][11] + kSqrt15_16 * (sh1[2][1] * sh6[8][11] - sh1[0][1] * sh6[4][11]) + -sqrt(1.0 / 8.0) * (sh1[2][1] * sh6[10][11] + sh1[0][1] * sh6[2][11]);
        sh7[10][13] = sqrt(40.0 / 13.0) * sh1[1][1] * sh6[9][12] + sqrt(45.0 / 26.0) * (sh1[2][1] * sh6[8][12] - sh1[0][1] * sh6[4][12]) + -sqrt(3.0 / 13.0) * (sh1[2][1] * sh6[10][12] + sh1[0][1] * sh6[2][12]);
        sh7[10][14] = sqrt(20.0 / 91.0) * (sh1[1][2] * sh6[9][12] - sh1[1][0] * sh6[9][0]) + sqrt(45.0 / 364.0) * ((sh1[2][2] * sh6[8][12] - sh1[2][0] * sh6[8][0]) - (sh1[0][2] * sh6[4][12] - sh1[0][0] * sh6[4][0])) + -sqrt(3.0 / 182.0) * ((sh1[2][2] * sh6[10][12] - sh1[2][0] * sh6[10][0]) + (sh1[0][2] * sh6[2][12] - sh1[0][0] * sh6[2][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[10]);

        sh7[11][0] = sqrt(33.0 / 182.0) * (sh1[1][2] * sh6[10][0] + sh1[1][0] * sh6[10][12]) + sqrt(55.0 / 364.0) * ((sh1[2][2] * sh6[9][0] + sh1[2][0] * sh6[9][12]) - (sh1[0][2] * sh6[3][0] + sh1[0][0] * sh6[3][12])) + -sqrt(3.0 / 364.0) * ((sh1[2][2] * sh6[11][0] + sh1[2][0] * sh6[11][12]) + (sh1[0][2] * sh6[1][0] + sh1[0][0] * sh6[1][12]));
        sh7[11][1] = sqrt(33.0 / 13.0) * sh1[1][1] * sh6[10][0] + sqrt(55.0 / 26.0) * (sh1[2][1] * sh6[9][0] - sh1[0][1] * sh6[3][0]) + -sqrt(3.0 / 26.0) * (sh1[2][1] * sh6[11][0] + sh1[0][1] * sh6[1][0]);
        sh7[11][2] = sqrt(11.0 / 8.0) * sh1[1][1] * sh6[10][1] + sqrt(55.0 / 48.0) * (sh1[2][1] * sh6[9][1] - sh1[0][1] * sh6[3][1]) + -kSqrt01_16 * (sh1[2][1] * sh6[11][1] + sh1[0][1] * sh6[1][1]);
        sh7[11][3] = sh1[1][1] * sh6[10][2] + kSqrt05_06 * (sh1[2][1] * sh6[9][2] - sh1[0][1] * sh6[3][2]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh6[11][2] + sh1[0][1] * sh6[1][2]);
        sh7[11][4] = sqrt(33.0 / 40.0) * sh1[1][1] * sh6[10][3] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[9][3] - sh1[0][1] * sh6[3][3]) + -sqrt(3.0 / 80.0) * (sh1[2][1] * sh6[11][3] + sh1[0][1] * sh6[1][3]);
        sh7[11][5] = sqrt(11.0 / 15.0) * sh1[1][1] * sh6[10][4] + sqrt(11.0 / 18.0) * (sh1[2][1] * sh6[9][4] - sh1[0][1] * sh6[3][4]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh6[11][4] + sh1[0][1] * sh6[1][4]);
        sh7[11][6] = sqrt(11.0 / 16.0) * sh1[1][1] * sh6[10][5] + sqrt(55.0 / 96.0) * (sh1[2][1] * sh6[9][5] - sh1[0][1] * sh6[3][5]) + -kSqrt01_32 * (sh1[2][1] * sh6[11][5] + sh1[0][1] * sh6[1][5]);
        sh7[11][7] = sqrt(33.0 / 49.0) * sh1[1][1] * sh6[10][6] + sqrt(55.0 / 98.0) * (sh1[2][1] * sh6[9][6] - sh1[0][1] * sh6[3][6]) + -sqrt(3.0 / 98.0) * (sh1[2][1] * sh6[11][6] + sh1[0][1] * sh6[1][6]);
        sh7[11][8] = sqrt(11.0 / 16.0) * sh1[1][1] * sh6[10][7] + sqrt(55.0 / 96.0) * (sh1[2][1] * sh6[9][7] - sh1[0][1] * sh6[3][7]) + -kSqrt01_32 * (sh1[2][1] * sh6[11][7] + sh1[0][1] * sh6[1][7]);
        sh7[11][9] = sqrt(11.0 / 15.0) * sh1[1][1] * sh6[10][8] + sqrt(11.0 / 18.0) * (sh1[2][1] * sh6[9][8] - sh1[0][1] * sh6[3][8]) + -sqrt(1.0 / 30.0) * (sh1[2][1] * sh6[11][8] + sh1[0][1] * sh6[1][8]);
        sh7[11][10] = sqrt(33.0 / 40.0) * sh1[1][1] * sh6[10][9] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[9][9] - sh1[0][1] * sh6[3][9]) + -sqrt(3.0 / 80.0) * (sh1[2][1] * sh6[11][9] + sh1[0][1] * sh6[1][9]);
        sh7[11][11] = sh1[1][1] * sh6[10][10] + kSqrt05_06 * (sh1[2][1] * sh6[9][10] - sh1[0][1] * sh6[3][10]) + -sqrt(1.0 / 22.0) * (sh1[2][1] * sh6[11][10] + sh1[0][1] * sh6[1][10]);
        sh7[11][12] = sqrt(11.0 / 8.0) * sh1[1][1] * sh6[10][11] + sqrt(55.0 / 48.0) * (sh1[2][1] * sh6[9][11] - sh1[0][1] * sh6[3][11]) + -kSqrt01_16 * (sh1[2][1] * sh6[11][11] + sh1[0][1] * sh6[1][11]);
        sh7[11][13] = sqrt(33.0 / 13.0) * sh1[1][1] * sh6[10][12] + sqrt(55.0 / 26.0) * (sh1[2][1] * sh6[9][12] - sh1[0][1] * sh6[3][12]) + -sqrt(3.0 / 26.0) * (sh1[2][1] * sh6[11][12] + sh1[0][1] * sh6[1][12]);
        sh7[11][14] = sqrt(33.0 / 182.0) * (sh1[1][2] * sh6[10][12] - sh1[1][0] * sh6[10][0]) + sqrt(55.0 / 364.0) * ((sh1[2][2] * sh6[9][12] - sh1[2][0] * sh6[9][0]) - (sh1[0][2] * sh6[3][12] - sh1[0][0] * sh6[3][0])) + -sqrt(3.0 / 364.0) * ((sh1[2][2] * sh6[11][12] - sh1[2][0] * sh6[11][0]) + (sh1[0][2] * sh6[1][12] - sh1[0][0] * sh6[1][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[11]);

        sh7[12][0] = sqrt(12.0 / 91.0) * (sh1[1][2] * sh6[11][0] + sh1[1][0] * sh6[11][12]) + sqrt(33.0 / 182.0) * ((sh1[2][2] * sh6[10][0] + sh1[2][0] * sh6[10][12]) - (sh1[0][2] * sh6[2][0] + sh1[0][0] * sh6[2][12])) + -sqrt(1.0 / 364.0) * ((sh1[2][2] * sh6[12][0] + sh1[2][0] * sh6[12][12]) + (sh1[0][2] * sh6[0][0] + sh1[0][0] * sh6[0][12]));
        sh7[12][1] = sqrt(24.0 / 13.0) * sh1[1][1] * sh6[11][0] + sqrt(33.0 / 13.0) * (sh1[2][1] * sh6[10][0] - sh1[0][1] * sh6[2][0]) + -sqrt(1.0 / 26.0) * (sh1[2][1] * sh6[12][0] + sh1[0][1] * sh6[0][0]);
        sh7[12][2] = sh1[1][1] * sh6[11][1] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh6[10][1] - sh1[0][1] * sh6[2][1]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh6[12][1] + sh1[0][1] * sh6[0][1]);
        sh7[12][3] = sqrt(8.0 / 11.0) * sh1[1][1] * sh6[11][2] + (sh1[2][1] * sh6[10][2] - sh1[0][1] * sh6[2][2]) + -sqrt(1.0 / 66.0) * (sh1[2][1] * sh6[12][2] + sh1[0][1] * sh6[0][2]);
        sh7[12][4] = kSqrt03_05 * sh1[1][1] * sh6[11][3] + sqrt(33.0 / 40.0) * (sh1[2][1] * sh6[10][3] - sh1[0][1] * sh6[2][3]) + -sqrt(1.0 / 80.0) * (sh1[2][1] * sh6[12][3] + sh1[0][1] * sh6[0][3]);
        sh7[12][5] = sqrt(8.0 / 15.0) * sh1[1][1] * sh6[11][4] + sqrt(11.0 / 15.0) * (sh1[2][1] * sh6[10][4] - sh1[0][1] * sh6[2][4]) + -sqrt(1.0 / 90.0) * (sh1[2][1] * sh6[12][4] + sh1[0][1] * sh6[0][4]);
        sh7[12][6] = sqrt(1.0 / 2.0) * sh1[1][1] * sh6[11][5] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[10][5] - sh1[0][1] * sh6[2][5]) + -sqrt(1.0 / 96.0) * (sh1[2][1] * sh6[12][5] + sh1[0][1] * sh6[0][5]);
        sh7[12][7] = sqrt(24.0 / 49.0) * sh1[1][1] * sh6[11][6] + sqrt(33.0 / 49.0) * (sh1[2][1] * sh6[10][6] - sh1[0][1] * sh6[2][6]) + -sqrt(1.0 / 98.0) * (sh1[2][1] * sh6[12][6] + sh1[0][1] * sh6[0][6]);
        sh7[12][8] = sqrt(1.0 / 2.0) * sh1[1][1] * sh6[11][7] + sqrt(11.0 / 16.0) * (sh1[2][1] * sh6[10][7] - sh1[0][1] * sh6[2][7]) + -sqrt(1.0 / 96.0) * (sh1[2][1] * sh6[12][7] + sh1[0][1] * sh6[0][7]);
        sh7[12][9] = sqrt(8.0 / 15.0) * sh1[1][1] * sh6[11][8] + sqrt(11.0 / 15.0) * (sh1[2][1] * sh6[10][8] - sh1[0][1] * sh6[2][8]) + -sqrt(1.0 / 90.0) * (sh1[2][1] * sh6[12][8] + sh1[0][1] * sh6[0][8]);
        sh7[12][10] = kSqrt03_05 * sh1[1][1] * sh6[11][9] + sqrt(33.0 / 40.0) * (sh1[2][1] * sh6[10][9] - sh1[0][1] * sh6[2][9]) + -sqrt(1.0 / 80.0) * (sh1[2][1] * sh6[12][9] + sh1[0][1] * sh6[0][9]);
        sh7[12][11] = sqrt(8.0 / 11.0) * sh1[1][1] * sh6[11][10] + (sh1[2][1] * sh6[10][10] - sh1[0][1] * sh6[2][10]) + -sqrt(1.0 / 66.0) * (sh1[2][1] * sh6[12][10] + sh1[0][1] * sh6[0][10]);
        sh7[12][12] = sh1[1][1] * sh6[11][11] + sqrt(11.0 / 8.0) * (sh1[2][1] * sh6[10][11] - sh1[0][1] * sh6[2][11]) + -sqrt(1.0 / 48.0) * (sh1[2][1] * sh6[12][11] + sh1[0][1] * sh6[0][11]);
        sh7[12][13] = sqrt(24.0 / 13.0) * sh1[1][1] * sh6[11][12] + sqrt(33.0 / 13.0) * (sh1[2][1] * sh6[10][12] - sh1[0][1] * sh6[2][12]) + -sqrt(1.0 / 26.0) * (sh1[2][1] * sh6[12][12] + sh1[0][1] * sh6[0][12]);
        sh7[12][14] = sqrt(12.0 / 91.0) * (sh1[1][2] * sh6[11][12] - sh1[1][0] * sh6[11][0]) + sqrt(33.0 / 182.0) * ((sh1[2][2] * sh6[10][12] - sh1[2][0] * sh6[10][0]) - (sh1[0][2] * sh6[2][12] - sh1[0][0] * sh6[2][0])) + -sqrt(1.0 / 364.0) * ((sh1[2][2] * sh6[12][12] - sh1[2][0] * sh6[12][0]) + (sh1[0][2] * sh6[0][12] - sh1[0][0] * sh6[0][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[12]);

        sh7[13][0] = sqrt(1.0 / 14.0) * (sh1[1][2] * sh6[12][0] + sh1[1][0] * sh6[12][12]) + sqrt(3.0 / 14.0) * ((sh1[2][2] * sh6[11][0] + sh1[2][0] * sh6[11][12]) - (sh1[0][2] * sh6[1][0] + sh1[0][0] * sh6[1][12]));
        sh7[13][1] = sh1[1][1] * sh6[12][0] + sqrt(3.0 / 1.0) * (sh1[2][1] * sh6[11][0] - sh1[0][1] * sh6[1][0]);
        sh7[13][2] = sqrt(13.0 / 24.0) * sh1[1][1] * sh6[12][1] + sqrt(13.0 / 8.0) * (sh1[2][1] * sh6[11][1] - sh1[0][1] * sh6[1][1]);
        sh7[13][3] = sqrt(13.0 / 33.0) * sh1[1][1] * sh6[12][2] + sqrt(13.0 / 11.0) * (sh1[2][1] * sh6[11][2] - sh1[0][1] * sh6[1][2]);
        sh7[13][4] = sqrt(13.0 / 40.0) * sh1[1][1] * sh6[12][3] + sqrt(39.0 / 40.0) * (sh1[2][1] * sh6[11][3] - sh1[0][1] * sh6[1][3]);
        sh7[13][5] = sqrt(13.0 / 45.0) * sh1[1][1] * sh6[12][4] + sqrt(13.0 / 15.0) * (sh1[2][1] * sh6[11][4] - sh1[0][1] * sh6[1][4]);
        sh7[13][6] = sqrt(13.0 / 48.0) * sh1[1][1] * sh6[12][5] + sqrt(13.0 / 16.0) * (sh1[2][1] * sh6[11][5] - sh1[0][1] * sh6[1][5]);
        sh7[13][7] = sqrt(13.0 / 49.0) * sh1[1][1] * sh6[12][6] + sqrt(39.0 / 49.0) * (sh1[2][1] * sh6[11][6] - sh1[0][1] * sh6[1][6]);
        sh7[13][8] = sqrt(13.0 / 48.0) * sh1[1][1] * sh6[12][7] + sqrt(13.0 / 16.0) * (sh1[2][1] * sh6[11][7] - sh1[0][1] * sh6[1][7]);
        sh7[13][9] = sqrt(13.0 / 45.0) * sh1[1][1] * sh6[12][8] + sqrt(13.0 / 15.0) * (sh1[2][1] * sh6[11][8] - sh1[0][1] * sh6[1][8]);
        sh7[13][10] = sqrt(13.0 / 40.0) * sh1[1][1] * sh6[12][9] + sqrt(39.0 / 40.0) * (sh1[2][1] * sh6[11][9] - sh1[0][1] * sh6[1][9]);
        sh7[13][11] = sqrt(13.0 / 33.0) * sh1[1][1] * sh6[12][10] + sqrt(13.0 / 11.0) * (sh1[2][1] * sh6[11][10] - sh1[0][1] * sh6[1][10]);
        sh7[13][12] = sqrt(13.0 / 24.0) * sh1[1][1] * sh6[12][11] + sqrt(13.0 / 8.0) * (sh1[2][1] * sh6[11][11] - sh1[0][1] * sh6[1][11]);
        sh7[13][13] = sh1[1][1] * sh6[12][12] + sqrt(3.0 / 1.0) * (sh1[2][1] * sh6[11][12] - sh1[0][1] * sh6[1][12]);
        sh7[13][14] = sqrt(1.0 / 14.0) * (sh1[1][2] * sh6[12][12] - sh1[1][0] * sh6[12][0]) + sqrt(3.0 / 14.0) * ((sh1[2][2] * sh6[11][12] - sh1[2][0] * sh6[11][0]) - (sh1[0][2] * sh6[1][12] - sh1[0][0] * sh6[1][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[13]);

        sh7[14][0] = kSqrt01_04 * ((sh1[2][2] * sh6[12][0] + sh1[2][0] * sh6[12][12]) - (sh1[0][2] * sh6[0][0] + sh1[0][0] * sh6[0][12]));
        sh7[14][1] = sqrt(7.0 / 2.0) * (sh1[2][1] * sh6[12][0] - sh1[0][1] * sh6[0][0]);
        sh7[14][2] = sqrt(91.0 / 48.0) * (sh1[2][1] * sh6[12][1] - sh1[0][1] * sh6[0][1]);
        sh7[14][3] = sqrt(91.0 / 66.0) * (sh1[2][1] * sh6[12][2] - sh1[0][1] * sh6[0][2]);
        sh7[14][4] = sqrt(91.0 / 80.0) * (sh1[2][1] * sh6[12][3] - sh1[0][1] * sh6[0][3]);
        sh7[14][5] = sqrt(91.0 / 90.0) * (sh1[2][1] * sh6[12][4] - sh1[0][1] * sh6[0][4]);
        sh7[14][6] = sqrt(91.0 / 96.0) * (sh1[2][1] * sh6[12][5] - sh1[0][1] * sh6[0][5]);
        sh7[14][7] = sqrt(13.0 / 14.0) * (sh1[2][1] * sh6[12][6] - sh1[0][1] * sh6[0][6]);
        sh7[14][8] = sqrt(91.0 / 96.0) * (sh1[2][1] * sh6[12][7] - sh1[0][1] * sh6[0][7]);
        sh7[14][9] = sqrt(91.0 / 90.0) * (sh1[2][1] * sh6[12][8] - sh1[0][1] * sh6[0][8]);
        sh7[14][10] = sqrt(91.0 / 80.0) * (sh1[2][1] * sh6[12][9] - sh1[0][1] * sh6[0][9]);
        sh7[14][11] = sqrt(91.0 / 66.0) * (sh1[2][1] * sh6[12][10] - sh1[0][1] * sh6[0][10]);
        sh7[14][12] = sqrt(91.0 / 48.0) * (sh1[2][1] * sh6[12][11] - sh1[0][1] * sh6[0][11]);
        sh7[14][13] = sqrt(7.0 / 2.0) * (sh1[2][1] * sh6[12][12] - sh1[0][1] * sh6[0][12]);
        sh7[14][14] = kSqrt01_04 * ((sh1[2][2] * sh6[12][12] - sh1[2][0] * sh6[12][0]) - (sh1[0][2] * sh6[0][12] - sh1[0][0] * sh6[0][0]));

        (*coeffs++) = dp(15, coeffsIn, sh7[14]);
    }
}

void SHL::RotateSH(const Mat3f& orient, int n, const float* coeffsIn, float* coeffs)
{
    ::RotateSH<float>(orient, n, coeffsIn, coeffs);
}

void SHL::RotateSH(const Mat3f& orient, int n, const Vec4f* coeffsIn, Vec4f* coeffs)
{
    ::RotateSH<Vec4f>(orient, n, coeffsIn, coeffs);
}


void SHL::ApplyNormalizationConstants(int n, Vec4f* coeffs)
{
    CL_ASSERT(n > 0);

    (*coeffs++) *= kSH_Y_00;

    if (n < 2)
        return;

    (*coeffs++) *= kSH_Y_1x;
    (*coeffs++) *= kSH_Y_1x;
    (*coeffs++) *= kSH_Y_1x;

    if (n < 3)
        return;

    (*coeffs++) *= kSH_Y_2_2;
    (*coeffs++) *= kSH_Y_2_1;
    (*coeffs++) *= kSH_Y_20;
    (*coeffs++) *= kSH_Y_21;
    (*coeffs++) *= kSH_Y_22;

    if (n < 4)
        return;

    (*coeffs++) *= kSH_Y_3_3;
    (*coeffs++) *= kSH_Y_3_2;
    (*coeffs++) *= kSH_Y_3_1;
    (*coeffs++) *= kSH_Y_30;
    (*coeffs++) *= kSH_Y_31;
    (*coeffs++) *= kSH_Y_32;
    (*coeffs++) *= kSH_Y_33;

    if (n < 5)
        return;

    (*coeffs++) *= kSH_Y_4_4;
    (*coeffs++) *= kSH_Y_4_3;
    (*coeffs++) *= kSH_Y_4_2;
    (*coeffs++) *= kSH_Y_4_1;
    (*coeffs++) *= kSH_Y_40;
    (*coeffs++) *= kSH_Y_41;
    (*coeffs++) *= kSH_Y_42;
    (*coeffs++) *= kSH_Y_43;
    (*coeffs++) *= kSH_Y_44;

    CL_ASSERT(n < 6);
}

void SHL::RemoveNormalizationConstants(int n, Vec4f* coeffs)
{
    CL_ASSERT(n > 0);

    (*coeffs++) *= 1.0f / kSH_Y_00;

    if (n < 2)
        return;

    (*coeffs++) *= 1.0f / kSH_Y_1x;
    (*coeffs++) *= 1.0f / kSH_Y_1x;
    (*coeffs++) *= 1.0f / kSH_Y_1x;

    if (n < 3)
        return;

    (*coeffs++) *= 1.0f / kSH_Y_2_2;
    (*coeffs++) *= 1.0f / kSH_Y_2_1;
    (*coeffs++) *= 1.0f / kSH_Y_20;
    (*coeffs++) *= 1.0f / kSH_Y_21;
    (*coeffs++) *= 1.0f / kSH_Y_22;

    if (n < 4)
        return;

    (*coeffs++) *= 1.0f / kSH_Y_3_3;
    (*coeffs++) *= 1.0f / kSH_Y_3_2;
    (*coeffs++) *= 1.0f / kSH_Y_3_1;
    (*coeffs++) *= 1.0f / kSH_Y_30;
    (*coeffs++) *= 1.0f / kSH_Y_31;
    (*coeffs++) *= 1.0f / kSH_Y_32;
    (*coeffs++) *= 1.0f / kSH_Y_33;

    if (n < 5)
        return;

    (*coeffs++) *= 1.0f / kSH_Y_4_4;
    (*coeffs++) *= 1.0f / kSH_Y_4_3;
    (*coeffs++) *= 1.0f / kSH_Y_4_2;
    (*coeffs++) *= 1.0f / kSH_Y_4_1;
    (*coeffs++) *= 1.0f / kSH_Y_40;
    (*coeffs++) *= 1.0f / kSH_Y_41;
    (*coeffs++) *= 1.0f / kSH_Y_42;
    (*coeffs++) *= 1.0f / kSH_Y_43;
    (*coeffs++) *= 1.0f / kSH_Y_44;

    CL_ASSERT(n < 6);
}

void SHL::ApplyDiffuseReflection(int n, float* coeffs)
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
        
    coeffs[4] *= kSH_A2;
    coeffs[5] *= kSH_A2;
    coeffs[6] *= kSH_A2;
    coeffs[7] *= kSH_A2;
    coeffs[8] *= kSH_A2;

    if (n < 4)
        return;

    for (int i = 9; i < 16; i++)
        coeffs[i] = 0;
    
    if (n < 5)
        return;
    
    for (int i = 16; i < 25; i++)
        coeffs[i] *= kSH_A4;
        
    if (n < 6)
        return;
    
    int nc = n * n;
    for (int i = 25; i < nc; i++)
        coeffs[i] = 0;
}

void SHL::ApplyGlossyReflection(float diff, float spec, int n, float coeffs[])
{
    coeffs[0] *= (spec + diff);
    
    if (n < 2)
        return;
    
    float w = spec + diff * kSH_A1;
    coeffs[1] *= w;
    coeffs[2] *= w;
    coeffs[3] *= w;

    if (n < 3)
        return;
        
    w = spec + diff * kSH_A2;
    coeffs[4] *= w;
    coeffs[5] *= w;
    coeffs[6] *= w;
    coeffs[7] *= w;
    coeffs[8] *= w;

    if (n < 4)
        return;

    for (int i = 9; i < 16; i++)
        coeffs[i] *= spec;
    
    if (n < 5)
        return;

    w = spec + diff * kSH_A4;
    for (int i = 16; i < 25; i++)
        coeffs[i] *= w;
        
    if (n < 6)
        return;
    
    int nc = n * n;
    for (int i = 25; i < nc; i++)
        coeffs[i] *= spec;
}

void SHL::ApplyGlossyReflection(float diff, float spec, int n, Vec4f coeffs[])
{
    coeffs[0] *= (spec + diff);
    
    if (n < 2)
        return;
    
    float w = spec + diff * kSH_A1;
    coeffs[1] *= w;
    coeffs[2] *= w;
    coeffs[3] *= w;

    if (n < 3)
        return;
        
    w = spec + diff * kSH_A2;
    coeffs[4] *= w;
    coeffs[5] *= w;
    coeffs[6] *= w;
    coeffs[7] *= w;
    coeffs[8] *= w;

    if (n < 4)
        return;

    for (int i = 9; i < 16; i++)
        coeffs[i] *= spec;
    
    if (n < 5)
        return;

    w = spec + diff * kSH_A4;
    for (int i = 16; i < 25; i++)
        coeffs[i] *= w;
        
    if (n < 6)
        return;
    
    int nc = n * n;
    for (int i = 25; i < nc; i++)
        coeffs[i] *= spec;
}

void SHL::ApplyGlossyReflection(const Vec4f& diff, const Vec4f& spec, int n, Vec4f coeffs[])
{
    coeffs[0] *= (spec + diff);
    
    if (n < 2)
        return;
    
    Vec4f w = spec + diff * kSH_A1;
    coeffs[1] *= w;
    coeffs[2] *= w;
    coeffs[3] *= w;

    if (n < 3)
        return;
        
    w = spec + diff * kSH_A2;
    coeffs[4] *= w;
    coeffs[5] *= w;
    coeffs[6] *= w;
    coeffs[7] *= w;
    coeffs[8] *= w;

    if (n < 4)
        return;

    for (int i = 9; i < 16; i++)
        coeffs[i] *= spec;
    
    if (n < 5)
        return;

    w = spec + diff * kSH_A4;
    for (int i = 16; i < 25; i++)
        coeffs[i] *= w;
        
    if (n < 6)
        return;
    
    int nc = n * n;
    for (int i = 25; i < nc; i++)
        coeffs[i] *= spec;
}

namespace
{
    inline float WindowScale(int n, float gamma)    // From PPS's SSHT
    {
        float nt = n * (n + 1.0f);
        return 1.0f / (1.0f + gamma * nt * nt);
    }
}

void SHL::ApplyZHWindowing(float gamma, int n, float* coeffs)
{
    for (int i = 0; i < n; i++)
        coeffs[i] *= WindowScale(i, gamma);
}

void SHL::ApplyZHWindowing(float gamma, int n, Vec4f* coeffs)
{
    for (int i = 0; i < n; i++)
        coeffs[i] *= WindowScale(i, gamma);
}

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

void SHL::ApplyDiffuseReflZH7(float zcoeffs[7])
{
    zcoeffs[1] *= kSH_A1;
    zcoeffs[2] *= kSH_A2;
    zcoeffs[3] = 0;
    zcoeffs[4] *= kSH_A4;
    zcoeffs[5] = 0;
    zcoeffs[6] = 0;
}

void SHL::ApplyDiffuseReflZH7(const float* zcoeffsIn, float* zcoeffsOut)
{
    zcoeffsOut[0] = zcoeffsIn[0];
    zcoeffsOut[1] = kSH_A1 * zcoeffsIn[1];
    zcoeffsOut[2] = kSH_A2 * zcoeffsIn[2];
    zcoeffsOut[3] = 0;
    zcoeffsOut[4] = kSH_A4 * zcoeffsIn[4];
    zcoeffsOut[5] = 0;
    zcoeffsOut[6] = 0;
}


float SHL::SampleZH7(float z, const float zcoeffs[7])
{
    float result = kZH_Y_0 * zcoeffs[0];

    result += kZH_Y_1 * z * zcoeffs[1];

    float z2 = z * z;
    result += kZH_Y_2 * (3 * z2 - 1) * zcoeffs[2];

    float z3 = z2 * z;
    result += kZH_Y_3 * (5 * z3 - 3 * z) * zcoeffs[3];

    float z4 = z2 * z2;
    result += kZH_Y_4 * (35  * z4 -  30 * z2 +  3) * zcoeffs[4];

    float z5 = z3 * z2;
    result += kZH_Y_5 * (63  * z5 -  70 * z3 +  15 * z) * zcoeffs[5];

    float z6 = z3 * z3;
    result += kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5) * zcoeffs[6];

    return result;
}

void SHL::AddZH7Sample(float x, float z, float zcoeffs[7])
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

float SHL::SampleZH(float z, int numBands, const float zcoeffs[])
{
    float result = kZH_Y_0 * zcoeffs[0];

    float z2 =  z *  z;
    float z3 =  z * z2;
    float z4 = z2 * z2;
    float z5 = z2 * z3;
    float z6 = z3 * z3;

    switch (numBands)
    {
    case 7:
        result += kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5) * zcoeffs[6];
    case 6:
        result += kZH_Y_5 * (63  * z5 -  70 * z3 +  15 * z) * zcoeffs[5];
    case 5:
        result += kZH_Y_4 * (35  * z4 -  30 * z2 +  3) * zcoeffs[4];
    case 4:
        result += kZH_Y_3 * (5 * z3 - 3 * z) * zcoeffs[3];
    case 3:
        result += kZH_Y_2 * (3 * z2 - 1) * zcoeffs[2];
    case 2:
        result += kZH_Y_1 * z * zcoeffs[1];
    }

    return result;
}

void SHL::AddZHSample(float x, float z, int numBands, float zcoeffs[])
{
    zcoeffs[0] += kZH_Y_0 * x;

    if (numBands < 2)
        return;

    zcoeffs[1] += kZH_Y_1 * z * x;

    if (numBands < 3)
        return;

    float z2 = z * z;
    zcoeffs[2] += kZH_Y_2 * (3 * z2 - 1) * x;

    if (numBands < 4)
        return;

    float z3 = z2 * z;
    zcoeffs[3] += kZH_Y_3 * (5 * z3 - 3 * z) * x;

    if (numBands < 5)
        return;

    float z4 = z2 * z2;
    zcoeffs[4] += kZH_Y_4 * (35  * z4 -  30 * z2 +  3) * x;

    if (numBands < 6)
        return;

    float z5 = z3 * z2;
    zcoeffs[5] += kZH_Y_5 * (63  * z5 -  70 * z3 +  15 * z) * x;

    if (numBands < 7)
        return;

    float z6 = z3 * z3;
    zcoeffs[6] += kZH_Y_6 * (231 * z6 - 315 * z4 + 105 * z2 - 5) * x;
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

void SHL::ConvolveSHWithZH(int numBands, const float* zcoeffs, float* coeffs)
{
    /*
        Example: Calling this with CalcCosPowerSatZH7(1, numBands) should be equivalent to the classic band-weighting from ApplyDiffuseReflection()

        const float kZH_Y_0 = sqrt( 1 / (   4 * pi)); 
        const float kZH_Y_1 = sqrt( 3 / (   4 * pi));
        const float kZH_Y_2 = sqrt( 5 / (  16 * pi));
        const float kZH_Y_3 = sqrt( 7 / (  16 * pi));
        const float kZH_Y_4 = sqrt( 9 / ( 256 * pi));

        zcoeffs[0] =   1 / 2 *  2 * pi * kZH_Y_0   * sqrt( 4 * pi / 1) = pi
        zcoeffs[1] =   1 / 3 *  2 * pi * kZH_Y_1   * sqrt( 4 * pi / 3) = 2 / 3 pi
        zcoeffs[2] =   1 / 4 *  2 * pi * kZH_Y_2   * sqrt( 4 * pi / 5) = 2 / 4 / 2 = 1/4 pi
        zcoeffs[3] =   0     
        zcoeffs[4] =  -1/6   *  2 * pi * kZH_Y_4   * sqrt( 4 * pi / 9) = -2 / 6 / 8 = - 1 / 24  pi

        kSH_A0 = 1.0,
        kSH_A1 = (2.0 / 3.0),
        kSH_A2 = (1.0 / 4.0),
        kSH_A3 = 0.0,
        kSH_A4 = -(1.0 / 24.0);

        Checks, except for pi factor -- reflective of CalcCosPowerSatZH7 being sat(cos(theta)) -- integral of that 
        // over sphere is pi, and if we want our diffuse BRDF to be normalised, it should thus be CalcCosPowerSatZH7(1) / pi.
    */

    for (int i = 0; i < numBands; i++)
    {
        int n = (2 * i + 1);
        float alpha = sqrtf(4.0f * vl_pi / n);
        float alphaZ = alpha * zcoeffs[i];

        for (int j = 0; j < n; j++)
            coeffs[j] *= alphaZ;

        coeffs += n;
    }
}

void SHL::ConvolveSHWithZH(int numBands, const float* zcoeffs, Vec4f* coeffs)
{
    for (int i = 0; i < numBands; i++)
    {
        int n = (2 * i + 1);
        float alpha = sqrtf(4.0f * vl_pi / n);
        float alphaZ = alpha * zcoeffs[i];

        for (int j = 0; j < n; j++)
            coeffs[j] *= alphaZ;

        coeffs += n;
    }
}


/////////////////////////////////////////////////////////////////////////////
// SH Multiplication. Requires triple-product basis coefficients.
//

void SHL::MultiplySHByZH3(const float* a, const float* b, float* c)
{
    c[0]  = a[0] * b[0] * kSH_Y_00;
    c[0] += a[1] * b[2] * kSH_Y_00;
    c[0] += a[2] * b[6] * kSH_Y_00;
    
    c[1]  = a[0] * b[1] * kSH_Y_00;
    c[1] += a[1] * b[5] * 0.218513f;
    c[1] += a[2] * b[1] * -0.126147f;
    c[2]  = a[0] * b[2] * kSH_Y_00;
    c[2] += a[1] * b[0] * kSH_Y_00;
    c[2] += a[1] * b[6] * 0.252306f;
    c[2] += a[2] * b[2] * 0.252306f;
    c[3]  = a[0] * b[3] * kSH_Y_00;
    c[3] += a[1] * b[7] * 0.218495f;
    c[3] += a[2] * b[3] * -0.126184f;
    
    c[4]  = a[0] * b[4] * 0.282134f;
    c[4] += a[2] * b[4] * -0.180258f;
    c[5]  = a[0] * b[5] * kSH_Y_00;
    c[5] += a[1] * b[1] * 0.218513f;
    c[5] += a[2] * b[5] * 0.0901372f;
    c[6]  = a[0] * b[6] * kSH_Y_00;
    c[6] += a[1] * b[2] * 0.252306f;
    c[6] += a[2] * b[0] * kSH_Y_00;
    c[6] += a[2] * b[6] * 0.180203f;
    c[7]  = a[0] * b[7] * kSH_Y_00;
    c[7] += a[1] * b[3] * 0.218495f;
    c[7] += a[2] * b[7] * 0.0900888f;
    c[8]  = a[8] * b[6] * -0.180213f;
    c[8] += a[0] * b[8] * kSH_Y_00;
    c[8] += a[2] * b[8] * -0.180213f;
}

namespace
{
    struct cBasisTriple
    {
        int8_t i;
        int8_t j;
        int8_t k;
        float  s;
    };

    cBasisTriple kBasisTriples3[] =
    {    
        0, 0, 0,  0.282095f,
        0, 1, 1,  0.282089f,
        0, 2, 2,  0.282087f,
        0, 3, 3,  0.282108f,
        0, 4, 4,  0.282134f,
        0, 5, 5,  0.282099f,
        0, 6, 6,  0.282096f,
        0, 7, 7,  0.282076f,
        0, 8, 8,  0.282069f,

        1, 1, 6, -0.126147f,
        1, 1, 8, -0.218467f,
        1, 2, 5,  0.218513f,
        1, 3, 4,  0.218540f,

        2, 2, 6,  0.252306f,
        2, 3, 7,  0.218495f,

        3, 3, 6, -0.126184f,
        3, 3, 8,  0.218513f,

        4, 4, 6, -0.180258f,
        4, 5, 7,  0.156089f,

        5, 5, 6,  0.0901372f,
        5, 5, 8, -0.156058f,

        6, 6, 6,  0.180203f,
        6, 7, 7,  0.0900888f,
        6, 8, 8, -0.180213f,

        7, 7, 8,  0.156056f,

        -1, -1, -1, 0.0f
    };

    cBasisTriple kBasisTriples4[] =
    {
        0,  0,  0,  0.282095f,
        0,  1,  1,  0.282089f,
        0,  2,  2,  0.282087f,
        0,  3,  3,  0.282108f,
        0,  4,  4,  0.282134f,
        0,  5,  5,  0.282099f,
        0,  6,  6,  0.282096f,
        0,  7,  7,  0.282076f,
        0,  8,  8,  0.282069f,
        0,  9,  9,  0.282099f,
        0, 10, 10,  0.282115f,
        0, 11, 11,  0.282121f,
        0, 12, 12,  0.282081f,
        0, 13, 13,  0.282082f,
        0, 14, 14,  0.282057f,
        0, 15, 15,  0.282109f,

        1,  1,  6, -0.126147f,
        1,  1,  8, -0.218467f,
        1,  2,  5,  0.218513f,
        1,  3,  4,  0.218540f,
        1,  4, 13, -0.0584171f,
        1,  4, 15, -0.226203f,
        1,  5, 12, -0.143029f,
        1,  5, 14, -0.184651f,
        1,  6, 11,  0.202316f,
        1,  7, 10,  0.184687f,
        1,  8,  9,  0.226142f,
        1,  8, 11,  0.0583776f,

        2,  2,  6,  0.252306f,
        2,  3,  7,  0.218495f,
        2,  4, 10,  0.184687f,
        2,  5, 11,  0.233626f,
        2,  6, 12,  0.247754f,
        2,  7, 13,  0.233563f,
        2,  8, 14,  0.184650f,

        3,  3,  6, -0.126184f,
        3,  3,  8,  0.218513f,
        3,  4,  9,  0.226223f,
        3,  4, 11, -0.0584171f,
        3,  5, 10,  0.184687f,
        3,  6, 13,  0.202305f,
        3,  7, 12, -0.143054f,
        3,  7, 14,  0.184648f,
        3,  8, 13, -0.0584227f,
        3,  8, 15,  0.226178f,

        4,  4,  6, -0.180258f,
        4,  5,  7,  0.156089f,
        4,  9, 13, -0.0940566f,
        4, 10, 12, -0.188072f,
        4, 11, 13,  0.145700f,
        4, 11, 15,  0.0940553f,

        5,  5,  6,  0.0901372f,
        5,  5,  8, -0.156058f,
        5,  9, 14,  0.148663f,
        5, 10, 13,  0.115178f,
        5, 10, 15, -0.148675f,
        5, 11, 12,  0.0594896f,
        5, 11, 14, -0.115172f,

        6,  6,  6,  0.180203f,
        6,  7,  7,  0.0900888f,
        6,  8,  8, -0.180213f,
        6,  9,  9, -0.210263f,
        6, 11, 11,  0.126194f,
        6, 12, 12,  0.168179f,
        6, 13, 13,  0.126122f,
        6, 15, 15, -0.210285f,

        7,  7,  8,  0.156056f,
        7,  9, 10,  0.148697f,
        7, 10, 11,  0.115178f,
        7, 12, 13,  0.0594762f,
        7, 13, 14,  0.115107f,
        7, 14, 15,  0.148656f,

        8,  9, 11, -0.0940076f,
        8, 11, 11, -0.145659f,
        8, 12, 14, -0.188046f,
        8, 13, 13,  0.145650f,
        8, 13, 15, -0.0940476f,

        -1, -1, -1, 0.0f
    };

    cBasisTriple kBasisTriples5[] =
    {    
         0,  0,  0, 0.282094791774f,
         0,  1,  1, 0.282088928036f,
         0,  2,  2, 0.282087283767f,
         0,  3,  3, 0.282108163520f,
         0,  4,  4, 0.282133668515f,
         0,  5,  5, 0.282099059533f,
         0,  6,  6, 0.282096391102f,
         0,  7,  7, 0.282075878233f,
         0,  8,  8, 0.282068961488f,
         0,  9,  9, 0.282098740287f,
         0, 10, 10, 0.282114588610f,
         0, 11, 11, 0.282121411308f,
         0, 12, 12, 0.282081078666f,
         0, 13, 13, 0.282081848443f,
         0, 14, 14, 0.282056774188f,
         0, 15, 15, 0.282109100917f,
         0, 16, 16, 0.282101344947f,
         0, 17, 17, 0.282100692635f,
         0, 18, 18, 0.282137443079f,
         0, 19, 19, 0.282115356411f,
         0, 20, 20, 0.282071599413f,
         0, 21, 21, 0.282093675305f,
         0, 22, 22, 0.282049260566f,
         0, 23, 23, 0.282072987179f,
         0, 24, 24, 0.282110766430f,

         1,  1,  6, -0.126155f,
         1,  1,  8, -0.218511f,
         1,  2,  5,  0.218512f,
         1,  3,  4,  0.218507f,
         1,  4, 13, -0.0583919f,
         1,  4, 15, -0.226175f,
         1,  5, 12, -0.14305f,
         1,  5, 14, -0.18467f,
         1,  6, 11,  0.202f,
         1,  7, 10,  0.18468f,
         1,  8,  9,  0.22618f,
         1,  8, 11,  0.058402f,
         1,  9, 22, -0.043528f,
         1,  9, 24, -0.23033f,
         1, 10, 21, -0.075395f,
         1, 10, 23, -0.19947f,
         1, 11, 20, -0.15078f,
         1, 11, 22, -0.16858f,
         1, 12, 19,  0.19466f,
         1, 13, 18,  0.16858f,
         1, 14, 17,  0.19947f,
         1, 14, 19,  0.075389f,
         1, 15, 16,  0.23032f,
         1, 15, 18,  0.043523f,

         2,  2,  6,  0.25230f,
         2,  3,  7,  0.21851f,
         2,  4, 10,  0.18468f,
         2,  5, 11,  0.23359f,
         2,  6, 12,  0.24776f,
         2,  7, 13,  0.23358f,
         2,  8, 14,  0.18467f,
         2,  9, 17,  0.16287f,
         2, 10, 18,  0.21325f,
         2, 11, 19,  0.23841f,
         2, 12, 20,  0.24623f,
         2, 13, 21,  0.23840f,
         2, 14, 22,  0.21323f,
         2, 15, 23,  0.16286f,

         3,  3,  6, -0.12615f,
         3,  3,  8,  0.21851f,
         3,  4,  9,  0.22617f,
         3,  4, 11, -0.058391f,
         3,  5, 10,  0.18468f,
         3,  6, 13,  0.20229f,
         3,  7, 12, -0.14305f,
         3,  7, 14,  0.18467f,
         3,  8, 13, -0.058400f,
         3,  8, 15,  0.22618f,
         3,  9, 16,  0.2303f,
         3,  9, 18, -0.043515f,
         3, 10, 17,  0.19948f,
         3, 10, 19, -0.075395f,
         3, 11, 18,  0.16858f,
         3, 12, 21,  0.19466f,
         3, 13, 20, -0.15078f,
         3, 13, 22,  0.16857f,
         3, 14, 21, -0.075403f,
         3, 14, 23,  0.19947f,
         3, 15, 22, -0.043530f,
         3, 15, 24,  0.23033f,

         4,  4,  6, -0.18021f,
         4,  4, 20,  0.040286f,
         4,  4, 24, -0.23841f,
         4,  5,  7,  0.15608f,
         4,  5, 21, -0.063720f,
         4,  5, 23, -0.16858f,
         4,  6, 18,  0.15607f,
         4,  7, 17,  0.16859f,
         4,  7, 19, -0.063720f,
         4,  8, 16,  0.23840f,
         4,  9, 13, -0.094020f,
         4, 10, 12, -0.1880f,
         4, 11, 13,  0.1456f,
         4, 11, 15,  0.094026f,
         4, 16, 22, -0.075073f,
         4, 17, 21, -0.11262f,
         4, 18, 20, -0.19036f,
         4, 18, 24,  0.075069f,
         4, 19, 21,  0.14189f,
         4, 19, 23,  0.11262f,

         5,  5,  6,  0.090112f,
         5,  5,  8, -0.15607f,
         5,  5, 20, -0.16120f,
         5,  5, 22, -0.18022f,
         5,  6, 19,  0.22072f,
         5,  7, 18,  0.18023f,
         5,  8, 17,  0.16858f,
         5,  8, 19,  0.06371f,
         5,  9, 14,  0.14867f,
         5, 10, 13,  0.11516f,
         5, 10, 15, -0.14867f,
         5, 11, 12,  0.05946f,
         5, 11, 14, -0.11516f,
         5, 16, 23,  0.1404f,
         5, 17, 22,  0.13272f,
         5, 17, 24, -0.14046f,
         5, 18, 21,  0.090299f,
         5, 18, 23, -0.13272f,
         5, 19, 20,  0.044860f,
         5, 19, 22, -0.090302f,

         6,  6,  6,  0.18022f,
         6,  6, 20,  0.24179f,
         6,  7,  7,  0.090104f,
         6,  7, 21,  0.22072f,
         6,  8,  8, -0.18022f,
         6,  8, 22,  0.15607f,
         6,  9,  9, -0.21025f,
         6, 11, 11,  0.12615f,
         6, 12, 12,  0.16820f,
         6, 13, 13,  0.12615f,
         6, 15, 15, -0.2102f,
         6, 16, 16, -0.22936f,
         6, 17, 17, -0.057344f,
         6, 18, 18,  0.06554f,
         6, 19, 19,  0.13925f,
         6, 20, 20,  0.1638f,
         6, 21, 21,  0.13925f,
         6, 22, 22,  0.065529f,
         6, 23, 23, -0.057348f,
         6, 24, 24, -0.22937f,

         7,  7,  8,  0.15607f,
         7,  7, 20, -0.16119f,
         7,  7, 22,  0.18021f,
         7,  8, 21, -0.063727f,
         7,  8, 23,  0.16858f,
         7,  9, 10,  0.14868f,
         7, 10, 11,  0.11516f,
         7, 12, 13,  0.059467f,
         7, 13, 14,  0.11515f,
         7, 14, 15,  0.14867f,
         7, 16, 17,  0.14046f,
         7, 17, 18,  0.13273f,
         7, 18, 19,  0.090299f,
         7, 20, 21,  0.04487f,
         7, 21, 22,  0.090285f,
         7, 22, 23,  0.13271f,
         7, 23, 24,  0.1404f,

         8,  8, 20,  0.040300f,
         8,  8, 24,  0.23842f,
         8,  9, 11, -0.094032f,
         8, 11, 11, -0.14567f,
         8, 12, 14, -0.18806f,
         8, 13, 13,  0.14566f,
         8, 13, 15, -0.094033f,
         8, 16, 18, -0.075073f,
         8, 17, 19, -0.11262f,
         8, 19, 19, -0.14188f,
         8, 20, 22, -0.19036f,
         8, 21, 21,  0.14189f,
         8, 21, 23, -0.11262f,
         8, 22, 24, -0.075090f,

         9,  9, 20,  0.07692f,
         9, 10, 21, -0.099327f,
         9, 11, 22,  0.13325f,
         9, 11, 24,  0.11752f,
         9, 12, 17, -0.20355f,
         9, 13, 16, -0.11750f,
         9, 13, 18,  0.1332f,
         9, 14, 19, -0.099322f,

        10, 10, 20, -0.17952f,
        10, 10, 24, -0.15172f,
        10, 11, 21,  0.10258f,
        10, 11, 23, -0.067851f,
        10, 12, 18, -0.04442f,
        10, 13, 17,  0.067855f,
        10, 13, 19,  0.10258f,
        10, 14, 16,  0.15171f,
        10, 15, 19,  0.099323f,

        11, 11, 20,  0.025631f,
        11, 11, 22, -0.11468f,
        11, 12, 19,  0.099312f,
        11, 13, 18,  0.11469f,
        11, 14, 17,  0.067852f,
        11, 14, 19, -0.10258f,
        11, 15, 16, -0.11751f,
        11, 15, 18, -0.13325f,

        12, 12, 20,  0.15387f,
        12, 13, 21,  0.09931f,
        12, 14, 22, -0.044418f,
        12, 15, 23, -0.20355f,

        13, 13, 20,  0.025640f,
        13, 13, 22,  0.11467f,
        13, 14, 21,  0.10257f,
        13, 14, 23,  0.067841f,
        13, 15, 22,  0.1332f,
        13, 15, 24, -0.11752f,

        14, 14, 20, -0.17951f,
        14, 14, 24,  0.15172f,
        14, 15, 21, -0.099328f,

        15, 15, 20,  0.076930f,

        16, 16, 20,  0.10651f,
        16, 17, 21, -0.119f,
        16, 18, 22,  0.13503f,
        16, 19, 23, -0.11909f,

        17, 17, 20, -0.15979f,
        17, 18, 21,  0.045014f,
        17, 19, 22,  0.045017f,
        17, 19, 24,  0.11910f,

        18, 18, 20, -0.083713f,
        18, 18, 24, -0.13504f,
        18, 19, 21,  0.10208f,
        18, 19, 23, -0.045016f,

        19, 19, 20,  0.068462f,
        19, 19, 22, -0.10208f,

        20, 20, 20,  0.13697f,
        20, 21, 21,  0.068467f,
        20, 22, 22, -0.083689f,
        20, 23, 23, -0.15978f,
        20, 24, 24,  0.10652f,

        21, 21, 22,  0.10208f,
        21, 22, 23,  0.045006f,
        21, 23, 24, -0.11910f,

        22, 22, 24,  0.13505f,
    };

    // 238 coeffs out of 2925 terms = 0.0813675 non-zero
    // density for 25 elts = 9.52


    void MultiplySH(int n, const float* a, const float* b, float* c, cBasisTriple* coeff)
    /// Multiplies two sets of coeffs together using given triple table.
    {
        for (int i = 0; i < n; i++)
            c[i] = 0;
            
        // tells: b[i > 0] = 0 -> copy a[i] * b[0]
        // i <= j <= k
        while (coeff->i >= 0)
        {
            if (coeff->i == coeff->j)
            {
                c[coeff->k] += coeff->s * a[coeff->i] * b[coeff->j];
                if (coeff->j != coeff->k)
                    c[coeff->i] += coeff->s * (a[coeff->j] * b[coeff->k] + a[coeff->k] * b[coeff->j]);            
            }
            else if (coeff->j == coeff->k)
            {
                c[coeff->i] += coeff->s * a[coeff->j] * b[coeff->j];
                c[coeff->j] += coeff->s * (a[coeff->i] * b[coeff->k] + a[coeff->k] * b[coeff->i]);        
            }
            else
            {
                // i != j != k
                c[coeff->i] += coeff->s * (a[coeff->j] * b[coeff->k] + a[coeff->k] * b[coeff->j]);
                c[coeff->j] += coeff->s * (a[coeff->i] * b[coeff->k] + a[coeff->k] * b[coeff->i]);
                c[coeff->k] += coeff->s * (a[coeff->i] * b[coeff->j] + a[coeff->j] * b[coeff->i]);
            }
            
            coeff++;
        }
    }

    void MultiplySH(int n, const Vec4f* a, const Vec4f* b, Vec4f* c, cBasisTriple* coeff)
    /// Multiplies two sets of coeffs together using given triple table.
    {
        for (int i = 0; i < n; i++)
            c[i] = vl_0;

        // tells: b[i > 0] = 0 -> copy a[i] * b[0]
        // i <= j <= k
        while (coeff->i >= 0)
        {
            if (coeff->i == coeff->j)
            {
                c[coeff->k] += coeff->s * a[coeff->i] * b[coeff->j];
                if (coeff->j != coeff->k)
                    c[coeff->i] += coeff->s * (a[coeff->j] * b[coeff->k] + a[coeff->k] * b[coeff->j]);            
            }
            else if (coeff->j == coeff->k)
            {
                c[coeff->i] += coeff->s * a[coeff->j] * b[coeff->j];
                c[coeff->j] += coeff->s * (a[coeff->i] * b[coeff->k] + a[coeff->k] * b[coeff->i]);        
            }
            else
            {
                // i != j != k
                c[coeff->i] += coeff->s * (a[coeff->j] * b[coeff->k] + a[coeff->k] * b[coeff->j]);
                c[coeff->j] += coeff->s * (a[coeff->i] * b[coeff->k] + a[coeff->k] * b[coeff->i]);
                c[coeff->k] += coeff->s * (a[coeff->i] * b[coeff->j] + a[coeff->j] * b[coeff->i]);
            }
            
            coeff++;
        }
    }
}

void SHL::MultiplySH3(const float* a, const float* b, float* c)
{
    MultiplySH(9, a, b, c, kBasisTriples3);
}

void SHL::MultiplySH4(const float* a, const float* b, float* c)
{
    MultiplySH(16, a, b, c, kBasisTriples4);
}

void SHL::MultiplySH5(const float* a, const float* b, float* c)
{
    MultiplySH(25, a, b, c, kBasisTriples5);
}

void SHL::MultiplySH(int numBands, const float* a, const float* b, float* c)
{
    switch (numBands)
    {
    case 3:
        MultiplySH(9,  a, b, c, kBasisTriples3);
        break;
    case 4:
        MultiplySH(16, a, b, c, kBasisTriples4);
        break;
    case 5:
        MultiplySH(25, a, b, c, kBasisTriples5);
        break;
    default:
        CL_ERROR("Unhandled band count\n");
    }
}


void SHL::MultiplySH(int numBands, const Vec4f* a, const Vec4f* b, Vec4f* c)
{
    switch (numBands)
    {
    case 3:
        MultiplySH(9,  a, b, c, kBasisTriples3);
        break;
    case 4:
        MultiplySH(16, a, b, c, kBasisTriples4);
        break;
    case 5:
        MultiplySH(25, a, b, c, kBasisTriples5);
        break;
    default:
        CL_ERROR("Unhandled band count\n");
    }
}

void SHL::FindSHCoeffsFromHemiEnvMap(const cImageData32* image, int n, float* coeffsR, float* coeffsG, float* coeffsB)
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

void SHL::FindSHCoeffsFromHemiEnvMap(const cImageData32* image, int n, Vec4f* coeffs)
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


void SHL::FindSHCoeffsFromHemiEnvMap(const cImageData48* image, int n, float* coeffsR, float* coeffsG, float* coeffsB)
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

void SHL::FindSHCoeffsFromHemiEnvMap(const cImageData48* image, int n, Vec4f* coeffs)
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

    struct cHDRFaceInfo
    {
        int   mTX, mTY;     // tile
        float mUX, mUY;     // u = mUX * x + mUY * y
        float mVX, mVY;
    };

    const cHDRFaceInfo kFaceTable[6] = 
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

void SHL::FindSHCoeffsFromHDRCubeMap(const cImageData48* image, int numBands, Vec4f* coeffs)
{
    int height = image->mHeight;
    int width  = image->mWidth;

    int faceSize = image->mWidth / 4;

    if (!IsPowerOfTwo(faceSize) || faceSize * 4 != image->mWidth || faceSize * 3 != image->mHeight)
    {
        CL_ERROR("invalid cube map!\n");
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
                float domega = domega0 / fr2;     // 4/w^2  *  (x^2 + y^2 + 1) ^ -3/2

                float u = kFaceTable[face].mUX * fx + kFaceTable[face].mUY * fy;
                float v = kFaceTable[face].mVX * fx + kFaceTable[face].mVY * fy;

                CL_ASSERT(pixel >= data && pixel < data + width * height * 3);

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

                AddColourSample(c, dir, numBands, coeffs);

                pixel += strideX;
            }
        }
    }
}

namespace
{
    void CalcHGPhaseZH(float g, const Vec4f& weightIn, int n, Vec4f* zcoeffs)
    {
        Vec4f weight(weightIn);

        for (int i = 0; i < n; i++)
        {
            zcoeffs[i] = sqrtf(4.0f * vl_pi * (i * 2 + 1)) * weight;
            weight *= g;
        }
    }

    void CalcHGPhaseZHAdd(float g, const Vec4f& weightIn, int n, Vec4f* zcoeffs)
    {
        Vec4f weight(weightIn);

        for (int i = 0; i < n; i++)
        {
            zcoeffs[i] += sqrtf(4.0f * vl_pi * (i * 2 + 1)) * weight;
            weight *= g;
        }
    }

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

    template<class M_SRC1, class M_SRC2> void MultiplyZH(int n, const M_SRC1* a, const M_SRC2* b, M_SRC2* c, cBasisTriple* coeff)
    /// Multiplies two sets of z coeffs together using given triple table.
    {
        // TODO: codegen this properly

        for (int i = 0; i < n; i++)
           c[i] = vl_0;

        // tells: b[i > 0] = 0 -> copy a[i] * b[0]
        // i <= j <= k
        while (coeff->i >= 0)
        {
            int bi = FloorToInt32(sqrtf(float(coeff->i)));
            int bj = FloorToInt32(sqrtf(float(coeff->j)));
            int bk = FloorToInt32(sqrtf(float(coeff->k)));

            if (   bi * (bi + 1) == coeff->i
                && bj * (bj + 1) == coeff->j
                && bk * (bk + 1) == coeff->k)
            {
                float bs = coeff->s;
                if (bi == bj)
                {
                    c[bk] += bs * a[bi] * b[bi];
                    if (bi != bk)
                        c[bi] += bs * (a[bk] * b[bi] + b[bk] * a[bi]);
                }
                else if (bj == bk)
                {
                    c[bi] += bs * a[bj] * b[bj];
                    c[bj] += bs * (a[bi] * b[bj] + b[bi] * a[bj]);
                }
                else
                {
                    // i != j != k
                    c[bi] += bs * (a[bj] * b[bk] + b[bj] * a[bk]);
                    c[bj] += bs * (a[bk] * b[bi] + b[bk] * a[bi]);
                    c[bk] += bs * (a[bi] * b[bj] + b[bi] * a[bj]);
                }
            }
            
            coeff++;
        }
    }
}

void SHL::CalcAtmosphereZH
(
    int          numPhases,
    const Vec4f* colourPhases,
    int          bands,
    Vec4f*       zcoeffs
)
{
    Vec4f phaseCoeffs[7];
 
    if (numPhases > 0)
        ::CalcHGPhaseZH(colourPhases[0][3], Vec4f(colourPhases[0]), 7, phaseCoeffs);

    for (int i = 1; i < numPhases; i++)
        ::CalcHGPhaseZHAdd(colourPhases[i][3], Vec4f(colourPhases[i]), 7, phaseCoeffs);

    float termZHCoeffs[7];
    GetTerminatorZH7(termZHCoeffs);

    switch (bands)
    {
    case 3:
        MultiplyZH<float, Vec4f>(bands, termZHCoeffs, phaseCoeffs, zcoeffs, kBasisTriples3);
        break;
    case 4:
        MultiplyZH<float, Vec4f>(bands, termZHCoeffs, phaseCoeffs, zcoeffs, kBasisTriples4);
        break;
    case 5:
        MultiplyZH<float, Vec4f>(bands, termZHCoeffs, phaseCoeffs, zcoeffs, kBasisTriples5);
        break;
    default:
        CL_ERROR("unhandled bands\n");
    }
}



void SHL::AddGroundBounce
(
    int          numBands,
    Vec4f*       coeffs,
    const Vec4f& diffColour,
    const Vec4f& specColour
)
{   
    CL_ASSERT(numBands <= 5);

    int nc = sqr(numBands);
    
    // get half mask
    float maskCoeffs[7];
    GetHemisphereZH7(maskCoeffs);

    // TODO: codegen a MultiplyZHSH(float*, Vec4f*, Vec4f*);
    Vec4f mask[25];
    for (int i = 0; i < nc; i++)
        mask[i] = vl_0;

    for (int i = 0; i < numBands; i++)
        mask[i * (i + 1)] = Vec4f(maskCoeffs[i], maskCoeffs[i], maskCoeffs[i], maskCoeffs[i]);

    Vec4f sky[25];
    switch (numBands)
    {
    case 3:
        MultiplySH(9,  coeffs, mask, sky, kBasisTriples3);
        break;
    case 4:
        MultiplySH(16, coeffs, mask, sky, kBasisTriples4);
        break;
    case 5:
        MultiplySH(25, coeffs, mask, sky, kBasisTriples5);
        break;
    default:
        CL_ASSERT(false);
    }
    for (int i = 0; i < nc; i++)
        coeffs[i] = sky[i];

    MirrorSHInZ(numBands, coeffs);    
    ApplyGlossyReflection(diffColour, specColour, numBands, coeffs);
    
    for (int i = 0; i < nc; i++)
        coeffs[i] += sky[i];
}

#ifdef EXAMPLE
// This is the unoptimized version, to show more clearly what we're doing. Which
// is just using the GatedSpot light model to represent the light as a sphere,
// and accumulating everything into the destination.
void SHL::AddSphereLighting(Vec3f pos, int numLights, cSphereLight* lights, Vec4f* coeffs)
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

        RotateZHToSHAdd(dir, 5, zcoeffs, coeffs);
    }
}
#endif

void SHL::AddSphereLighting(Vec3f pos, float scale, int numLights, const cSphereLight* lights, Vec4f* coeffs)
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
        float z0 = strength * (vl_twoPi * kZH_Y_0) * (1.0f - t);
        float z1 = strength * (vl_twoPi * kZH_Y_1) * 0.5f * (1.0f - t2);
        float z2 = strength * (vl_twoPi * kZH_Y_2) * t * (1.0f - t2);
        float z3 = strength * (vl_twoPi * kZH_Y_3) * 0.25f * (t2 * (6.0f - 5.0f * t2) - 1.0f); 

        zcoeffs[0] = colour * z0;
        zcoeffs[1] = colour * z1;
        zcoeffs[2] = colour * z2;
        zcoeffs[3] = colour * z3;

        // TODO: break this out and insert above. Also, time -- is it faster to do scalar and then fold in colour at the end, or to keep to an SSE path?
        // actually, TODO: write version that does 4 lights at a time.
        RotateZHToSHAdd(dir, bands, zcoeffs, coeffs); 
    }
}

void SHL::CreateSphereLightBRDF(int width, int height, uint32_t* pixels)
{
    float zcoeffsSpec[7];
    float zcoeffsDiff[7];
    float t;

    for (int i = 0; i < height; i++)
    {
        if (i == height - 1)
            t = 1.0f;   // ensure last row is 0, or we'll contribute light off to infinity.
        else
            t = (i + 0.5f) / height;

        CalcGatedSpotZH7(t, zcoeffsSpec);
        ApplyDiffuseReflZH7(zcoeffsSpec, zcoeffsDiff);

        for (int j = 0; j < width; j++)
        {
            float z = (j + 0.5f) / width;
            z = 2 * z - 1;

            float diff = SampleZH7(z, zcoeffsDiff);
            float spec = SampleZH7(z, zcoeffsSpec);

            diff = Saturate(diff);
            spec = Saturate(spec);

            *pixels++ = ColourToU32(Vec4f(diff, spec, 0, 1));  // TODO: extra stuff we can cram in?
        }
    }
}

