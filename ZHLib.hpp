//
// ZHLib.hpp
//
// Spherical Harmonics utility library
//
// Andrew Willmott
//

#ifndef ZH_LIB_H
#define ZH_LIB_H

#include "VL234f.hpp"

namespace ZHL
{
    // ZH utilities
    void  CalcZHWeights(float z, int n, float weights[]);           // Generate 1-7 band ZH weights for 'z'
    float SampleZH     (float z, int n, const float zCoeffs[]);     // Sample 1-7 band ZH.
    float SampleZH_p1  (         int n, const float zCoeffs[]);     // Sample zh(1.0), useful for finding max
    void  AddZHSample  (float x, float z, int n, float zCoeffs[]);  // Add given sample f(z) = x to 1-7 band zCoeffs.

    float ZH(int l, float z);     // Evaluate ZH basis function. 'l' = band, z is cos(theta), -1 to 1

    extern const float kZH_Y[7];  // First 7 ZH coeffients. Note that these include any constant factor from the corresponding polynomial.

    // ZH7 variants. These always calculate 7 bands.
    void  CalcZH7Weights(float z, float weights[7]);           // Generate ZH coefficient weights for 'z'
    float SampleZH7     (float z, const float zCoeffs[7]);     // Sample given set of ZH coeffs according to z.
    float SampleZH7_p1  (         const float zCoeffs[7]);     // Sample zh(1.0), useful for finding max
    void  AddZH7Sample  (float x, float z, float zCoeffs[7]);  // Add given sample f(z) = x to zCoeffs. Remember to weight x by dz

    void  CalcCosPowerZH7   (int   n, float zCoeffs[7]);  // Find first 7 ZH coeffs for cos(theta)^n function. Note: this is signed! You likely want CalcCosPowerSatZH7() instead.
    void  CalcCosPowerSatZH7(float n, float zCoeffs[7]);  // Find first 7 ZH coeffs for max(0, cos(theta)^n) function.
    void  CalcHemisphereZH7 (         float zCoeffs[7]);  // Fill in first 7 ZH coeffs for upper hemisphere: f(z) = delta(z).
    void  CalcGatedSpotZH7  (float t, float zCoeffs[7]);  // Find first 7 ZH coeffs for given gated spot illumination: f(z) = 1 for z > t, 0 otherwise.
    void  CalcGatedCosZH7   (float t, float zCoeffs[7]);  // Find first 7 ZH coeffs for given gated cos(theta) illumination: f(z) = z for z > t, 0 otherwise.

    void  CalcGatedCosBands5(float t, float bandScale[5]); // Find band weights for energy-preserving convolution with given gated cos function. This can be used to vary between diffuse (t=0) and glossy (t=1) reflection.

    void  MirrorZH(int n, float coeffs[]);  // Mirror given ZH coeffs in Z.
    void  MirrorZH(int n, Vec4f coeffs[]);  // Mirror given ZH coeffs in Z.

    void  ApplyDiffuseBRDFZH7(float zCoeffs[7]);                             // Apply diffuse BRDF convolution to zCoeffs.
    void  ApplyDiffuseBRDFZH7(const float zCoeffsIn[], float zCoeffsOut[]);  // Apply diffuse BRDF convolution to zCoeffs.

    void  CalcHGPhaseZH          (float g, float weight, int n, float zCoeffs[]);  // Find ZH coeffs for given Henyey-Greenstein phase function.
    void  CalcNormalizedHGPhaseZH(float g, float weight, int n, float zCoeffs[]);  // CalcHGPhaseZH, but normalized so that the max of the function is always 1.

    void  ConvolveZH(int n, const float brdfCoeffs[], const float zCoeffsIn[], float zCoeffsOut[]);

    void  ApplyZHMaxScale(int n, float zCoeffs[]);

    void  ApplyZHWindowing(float gamma, int n, float coeffs[]);  // Apply windowing to coefficients -- trades ringing for blur. Larger gamma = more, but start with very small epsilon
    void  ApplyZHWindowing(float gamma, int n, Vec4f coeffs[]);  // Apply windowing to coefficients -- trades ringing for blur. Larger gamma = more, but start with very small epsilon

    float WindowScale(int n, float gamma); // returns windowing scaling factor for band n

    // Multiplication
    void MultiplyZH(int n, const float a[], const float b[], float c[]); // Multiplies two sets of coeffs together, n = 1..7
    void MultiplyZH(int n, const Vec4f a[], const Vec4f b[], Vec4f c[]); // Multiplies two sets of coeffs together, n = 1..7

    // Atmosphere
    void CalcAtmosphereZH(int numPhases, const Vec4f colourPhases[], int n, Vec4f zCoeffs[]); // Find ZH coeffs for given atmosphere model, which is just a sum of colour * HG(phase).


    //
    // Inlines
    //

    inline float WindowScale(int n, float gamma)
    {
        float nt = n * (n + 1.0f);
        return 1.0f / (1.0f + gamma * nt * nt);
    }
}

#endif
