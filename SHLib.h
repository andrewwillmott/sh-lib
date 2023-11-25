//------------------------------------------------------------------------------
// Spherical/Zonal Harmonics Utility Library
//------------------------------------------------------------------------------

#ifndef SH_LIB_H
#define SH_LIB_H

#include "VL234f.h"

namespace SHL
{
    // Basic SH utilities
    float SampleSH     (         Vec3f d, int n, const float coeffs[]);  // Sample SH coeffs in direction 'd'.
    Vec4f SampleSH     (         Vec3f d, int n, const Vec4f coeffs[]);  // Sample SH coeffs in direction 'd'.
    void  CalcSHWeights(         Vec3f d, int n, float weights[]);       // Generate SH weights for direction 'd'
    void  AddSHSample  (float s, Vec3f d, int n, float coeffs[]);        // Add given sample 's' in dir 'd'  into  SH coeffs.
    void  AddSHSample  (Vec4f s, Vec3f d, int n, Vec4f coeffs[]);        // Add given sample 's' in dir 'd'  into  SH coeffs.

    void  AddSHSampleHemi(float s, Vec3f d, int n, float coeffs[]);      // Add sample 's' in dir 'd'  and mirror_z(d). Useful for functions only defined on the hemisphere.
    void  AddSHSampleHemi(Vec4f s, Vec3f d, int n, Vec4f coeffs[]);      // Add given sample 's' in dir 'd'  into  SH coeffs.

    float SH(int l, int m, Vec3f d);                // Evaluate SH basis function. 'l' = band, 'm' = -l .. +l. 'd' must be normalised
    float SH(int l, int m, float theta, float phi); // Evaluate SH basis function. 'l' = band, 'm' = -l .. +l. 'theta' and 'phi' in radians.

    extern const float  kSH_Y[25];  // First 25 SH coeffients (5 bands). Note that these include any constant factor from the corresponding polynomial.
    extern const float* kSH_Y0;     // kSH_Y0[1]
    extern const float* kSH_Y1;     // kSH_Y1[3]
    extern const float* kSH_Y2;     // kSH_Y2[5]
    extern const float* kSH_Y3;     // kSH_Y3[7]
    extern const float* kSH_Y4;     // kSH_Y4[9]

    // Rotation/mirror utilities. Note: the coordinate system is assumed to be right-handed!
    void RotateZHToSH(const Vec3f& dir, int n, const float zCoeffs[], float shCoeffs[]);  // Apply rotate zCoeffs to align with 'dir', and store the result in shCoeffs. Handles n <= 8.
    void RotateZHToSH(const Vec3f& dir, int n, const Vec4f zCoeffs[], Vec4f shCoeffs[]);  // Apply rotate zCoeffs to align with 'dir', and store the result in shCoeffs. Handles n <= 8.

    void RotateZHToSHAdd(const Vec3f& dir, int n, const float zCoeffs[], float shCoeffs[]);  // Accumulating version of RotateZHToSH
    void RotateZHToSHAdd(const Vec3f& dir, int n, const Vec4f zCoeffs[], Vec4f shCoeffs[]);  // Accumulating version of RotateZHToSH

    void RotateSHAboutZ(float theta, int n, float coeffs[]);  // Rotate given SH coefficients by 'theta' about Z. Fast -- O(n).
    void RotateSHAboutZ(float theta, int n, Vec4f coeffs[]);  // Rotate given SH coefficients by 'theta' about Z. Fast -- O(n).

    void MirrorSHInZ(int n, float coeffs[]);  // Mirror given SH coeffs in Z.
    void MirrorSHInZ(int n, Vec4f coeffs[]);  // Mirror given SH coeffs in Z.

    void RRotateSH(const Mat3f& rot_row, int n, const float coeffsIn[], float coeffsOut[]);  // Generic SH rotation by given rotation on row vectors. Pretty fast, but O(n^2). Handles n <= 8.
    void RRotateSH(const Mat3f& rot_row, int n, const Vec4f coeffsIn[], Vec4f coeffsOut[]);  // Generic SH rotation by given rotation on row vectors.
    void CRotateSH(const Mat3f& rot_col, int n, const float coeffsIn[], float coeffsOut[]);  // Generic SH rotation by given rotation on column vectors.
    void CRotateSH(const Mat3f& rot_col, int n, const Vec4f coeffsIn[], Vec4f coeffsOut[]);  // Generic SH rotation by given rotation on column vectors.

    // Normalisation and convolution
    void ApplyNormalizationConstants (int n, float coeffs[]);  // Premultiplies in normalization constants to coeffs. Only handles n <= 5.
    void ApplyNormalizationConstants (int n, Vec4f coeffs[]);  // Premultiplies in normalization constants to coeffs. Only handles n <= 5.
    void RemoveNormalizationConstants(int n, float coeffs[]);  // Undoes the effects of ApplyNormalizationConstants. Only handles n <= 5.
    void RemoveNormalizationConstants(int n, Vec4f coeffs[]);  // Undoes the effects of ApplyNormalizationConstants. Only handles n <= 5.

    void ApplyMaxScale(int n, float coeffs[]);  // Scales coeffs so +1 = maximum value of 1 in corresponding basis function. Useful for basis function display. Handles up to 5 bands.
    void ApplyMaxScale(int n, Vec4f coeffs[]);  // Scales coeffs so +1 = maximum value of 1 in corresponding basis function.

    void ApplyDiffuseBRDF(int n, float coeffs[]);  // Convolves lighting environment 'coeffs' with the diffuse BRDF
    void ApplyDiffuseBRDF(int n, Vec4f coeffs[]);

    void ApplyGlossyBRDF(float gloss, int n, float coeffs[]);
    void ApplyGlossyBRDF(float gloss, int n, Vec4f coeffs[]);
    // Convolves the lighting environment with a glossy BRDF, defined as a mix between a diffuse
    // BRDF and the band-limited 'specular' BRDF from band limiting to n. gloss=0 -> diffuse, gloss=1 -> max gloss.

    void ConvolveSHWithZH(int n, const float zCoeffs[], const float shCoeffsIn[], float shCoeffsOut[]);
    void ConvolveSHWithZH(int n, const float zCoeffs[], const Vec4f shCoeffsIn[], Vec4f shCoeffsOut[]);
    // Convolve n-band SH coeffs with n-band ZH coeffs. shCoeffsIn == shCoeffsOut is valid.

    // Windowing
    void ApplySHWindowing(float gamma, int n, float coeffs[]);  // Apply windowing to coefficients -- trades ringing for blur. Larger gamma = more, but start with very small epsilon
    void ApplySHWindowing(float gamma, int n, Vec4f coeffs[]);  // Apply windowing to coefficients -- trades ringing for blur. Larger gamma = more, but start with very small epsilon

    // Multiplication
    void MultiplySH(int n, const float a[], const float b[], float c[]); // Multiplies two sets of coeffs together, n = 1..5
    void MultiplySH(int n, const Vec4f a[], const Vec4f b[], Vec4f c[]); // Multiplies two sets of coeffs together, n = 1..5

    void MultiplySHByZH3(const float a[9], const float b[3], float c[9]); // Multiply 3-coeff ZH 'a' by 9-coeff SH 'b', store in 9-coeff SH 'c'.

    // Environment map support
    struct ImageData32
    {
        const uint32_t* mData;  // 32-bit BGRA8
        int             mWidth;
        int             mHeight;
    };

    struct ImageData48
    {
        const uint16_t* mData;  // 48-bit RGB16
        int             mWidth;
        int             mHeight;
    };

    void FindSHCoeffsFromHemiEnvMap(const ImageData32* image, int n, float* coeffsR, float* coeffsG, float* coeffsB);
    void FindSHCoeffsFromHemiEnvMap(const ImageData48* image, int n, float* coeffsR, float* coeffsG, float* coeffsB);
    // Extract coeffs from Devebec-style hemi env map.

    void FindSHCoeffsFromHemiEnvMap(const ImageData32* image, int n, Vec4f coeffs[]);
    void FindSHCoeffsFromHemiEnvMap(const ImageData48* image, int n, Vec4f coeffs[]);
    // Extract coeffs from Devebec-style hemi env map.

    void FindSHCoeffsFromHDRCubeMap(const ImageData48* image, int n, Vec4f coeffs[]);
    // Extract coeffs from 16-bit cube map.

    // Misc lighting
    struct SphereLight
    {
        Vec4f   mColourAndIntensity;    // r, g, b, intensity
        Vec4f   mPositionAndSize;       // location, size
    };

    void AddSphereLighting(Vec3f pos, float scale, int numLights, const SphereLight lights[], Vec4f coeffs[]);
    // Add in given sphere lights. 'coeffs' must be at least four bands (16 coefficients.)
    // Does not expect them to be sorted in importance order -- will choose the appropriate approximation depending on light size and distance.
}

#endif
