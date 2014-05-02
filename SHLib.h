//------------------------------------------------------------------------------
// Spherical/Zonal Harmonics Utility Library
//------------------------------------------------------------------------------

#ifndef SH_LIB_H
#define SH_LIB_H

#include "VL234f.h"

namespace SHL
{
    // Basic SH utilities

    void AddColourSample(const Vec4f& s, const Vec3f& dir, int numBands, Vec4f* coeffs);
    ///< Accumulate colour s sampled in direction dir into the given coeffs. Handles 2-5 bands.
    void AddReflectedColourSample(const Vec4f& s, const Vec3f& dir, int numBands, Vec4f* coeffs);
    ///< Cheaper version of AddColourSample that assumes the destination is reflected about the xy plane,
    ///< useful for functions defined only over one hemisphere.

    Vec4f SampleColour(const Vec4f* coeffs, int numBands, const Vec3f& v);
    ///< Reconstruct colour in given direction. Handles between 2 and 5 bands. (4 to 25 coeffs.)

    // 7-band ZH utilities, all assume 7 coefficients.

    void CalcCosPowerZH7(int n, float zcoeffs[7]);
    ///< Find first 7 ZH coeffs for cos(theta)^n function. Note: this is signed! You likely want CalcCosPowerSatZH7() instead.
    void CalcCosPowerSatZH7(float n, float zcoeffs[7]);
    ///< Find first 7 ZH coeffs for max(0, cos(theta)^n) function.
    void GetHemisphereZH7(float zcoeffs[7]);
    ///< Fill in first 7 ZH coeffs for upper hemisphere: f(z) = delta(z).
    void CalcGatedSpotZH7(float t, float zcoeffs[7]);
    ///< Find first 7 ZH coeffs for given gated spot illumination: f(z) = 1 for z > t, 0 otherwise.
    void CalcGatedCosZH7(float t, float zcoeffs[7]);
    ///< Find first 7 ZH coeffs for given gated cos(theta) illumination: f(z) = z for z > t, 0 otherwise.
    
    void CalcGatedCosBands5(float t, float bandMult[5]);
    ///< Find band weights for energy-preserving convolution with given gated cos function.
    ///< This can be used to vary between diffuse (t=0) and glossy (t=1) reflection.
    
    void CalcHGPhaseZH(float g, float weight, int n, float zcoeffs[]);
    ///< Find ZH coeffs for given Henyey-Greenstein phase function.
    void CalcNormalizedHGPhaseZH(float g, float weight, int n, float zcoeffs[]);
    // CalcHGPhaseZH, but normalized so that the max of the function is always 1.

    void ApplyDiffuseReflZH7(float zcoeffs[7]);
    void ApplyDiffuseReflZH7(const float* zcoeffsIn, float* zcoeffsOut);
    ///< Apply diffuse convolution to zcoeffs.

    float SampleZH7(float z, const float zcoeffs[7]);
    ///< Sample given set of ZH coeffs according to z.
    void  AddZH7Sample(float x, float z, float zcoeffs[7]);
    ///< Add given sample f(z) = x to zcoeffs. Remember to weight x by dz

    float SampleZH(float z, int numBands, const float zcoeffs[]);
    ///< Sample 1-7 band ZH.
    void  AddZHSample(float x, float z, int numBands, float zcoeffs[]);
    ///< Add given sample f(z) = x to 1-7 band zcoeffs.


    // ZH/SH rotation utilities.

    void RotateZHToSH(const Vec3f& dir, int numBands, const float* zcoeffs, float* coeffs);
    void RotateZHToSH(const Vec3f& dir, int numBands, const Vec4f* zcoeffs, Vec4f* coeffs);
    ///< Apply orientation to given ZHs, store result in SHs. Because the source
    ///< function is rotationally symmetric about z, we only need an orientation
    ///< direction, not the full 3x3 matrix. Handles n <= 8.

    void RotateZHToSHAdd(const Vec3f& dir, int numBands, const float* zcoeffs, float* coeffs);
    void RotateZHToSHAdd(const Vec3f& dir, int numBands, const Vec4f* zcoeffs, Vec4f* coeffs);
    void RotateZHToSHAdd(const Vec3f& dir, int numBands, const Vec4f& c, const float* zcoeffs, Vec4f* coeffs);
    ///< Apply orientation to given ZHs, multiply by c, and accumulate the result in coeffs.
    ///< Handles n <= 8.

    void RotateSHAboutZ(float theta, int numBands, float* coeffs);
    ///< Rotate given SH coefficients by 'theta' about Z. Fast -- O(n).
    void RotateSHAboutZ(float theta, int numBands, Vec4f* coeffs);
    ///< Rotate given SH coefficients by 'theta' about Z. Fast -- O(n).
    void RotateSHAboutZ(float theta, int n, const float* coeffsIn, float* coeffs);
    ///< Rotate given SH coefficients by 'theta' about Z. Fast -- O(n).

    void RotateSH(const Mat3f& orient, int n, const float* coeffsIn, float* coeffs);
    ///< Generic SH coefficient rotation by 'orient'. Pretty fast, but O(n^2). Handles n <= 8.
    void RotateSH(const Mat3f& orient, int n, const Vec4f* coeffsIn, Vec4f* coeffs);
    ///< Generic SH coefficient rotation by 'orient'. Pretty fast, but O(n^2). Handles n <= 8.

    // Normalisation and convolution
    
    void ApplyNormalizationConstants(int n, Vec4f* coeffs);
    ///< Premultiplies in normalization constants to coeffs. Only handles n <= 5.
    void RemoveNormalizationConstants(int n, Vec4f* coeffs);
    ///< Undoes the effects of ApplyNormalizationConstants. Only handles n <= 5.

    void ApplyDiffuseReflection(int n, float coeffs[]);
    ///< Simulates convolution of lighting environment with a diffuse BRDF
    void ApplyGlossyReflection(float diff, float spec, int n, float coeffs[]);
    ///< Simulates convolution of lighting environment with a glossy BRDF, defined as a mix between a diffuse 
    ///< BRDF and the band-limited 'specular' BRDF from band limiting to n. Normally diff + spec = 1 to maintain energy.
    void ApplyGlossyReflection(float diff, float spec, int n, Vec4f coeffs[]);
    void ApplyGlossyReflection(const Vec4f& diff, const Vec4f& spec, int n, Vec4f coeffs[]);
    ///< Glossy reflection for full RGBA.

    void ApplyZHWindowing(float gamma, int n, float*   coeffs);
    void ApplyZHWindowing(float gamma, int n, Vec4f* coeffs);
    void ApplySHWindowing(float gamma, int n, float*   coeffs);
    void ApplySHWindowing(float gamma, int n, Vec4f* coeffs);
    ///< Apply windowing to coefficients -- trades ringing for blur. Larger gamma = more, but start with very small epsilon

    void MirrorSHInZ(int n, float* coeffs);
    ///< Mirror given SH coeffs about the XY plane.
    void MirrorSHInZ(int n, Vec4f* coeffs);
    ///< Mirror given SH coeffs about the XY plane.

    void ConvolveSHWithZH(int numBands, const float* zcoeffs, float* coeffs);
    void ConvolveSHWithZH(int numBands, const float* zcoeffs, Vec4f* coeffs);
    ///< Convolve n-band SH coeffs with n-band ZH coeffs.

    void MultiplySHByZH3(const float* a, const float* b, float* c);
    ///< Multiply 3-coeff ZH 'a' by 9-coeff SH 'b', store in 9-coeff SH 'c'.
    void MultiplySH3(const float* a, const float* b, float* c);
    ///< Multiply two 3-band SH coeff sets together.
    void MultiplySH4(const float* a, const float* b, float* c);
    ///< Multiply two 4-band SH coeff sets together.
    void MultiplySH5(const float* a, const float* b, float* c);
    ///< Multiply two 5-band SH coeff sets together.
    void MultiplySH(int numBands, const float* a, const float* b, float* c);
    /// Multiplies two sets of coeffs together -- currently only handles numBands = 3, 4, 5.
    void MultiplySH(int numBands, const Vec4f* a, const Vec4f* b, Vec4f* c);
    /// Multiplies two sets of coeffs together -- currently only handles numBands = 3, 4, 5.

    // Environment map support
    struct cImageData32
    {
        const uint32_t* mData;  // 32-bit BGRA8
        int mWidth;
        int mHeight;
    };
    struct cImageData48
    {
        const uint16_t* mData;  // 48-bit RGB16
        int mWidth;
        int mHeight;
    };

    void FindSHCoeffsFromHemiEnvMap(const cImageData32* image, int n, float* coeffsR, float* coeffsG, float* coeffsB);
    void FindSHCoeffsFromHemiEnvMap(const cImageData48* image, int n, float* coeffsR, float* coeffsG, float* coeffsB);
    ///< Extract coeffs from Devebec-style hemi env map.

    void FindSHCoeffsFromHemiEnvMap(const cImageData32* image, int n, Vec4f* coeffs);
    void FindSHCoeffsFromHemiEnvMap(const cImageData48* image, int n, Vec4f* coeffs);
    ///< Extract coeffs from Devebec-style hemi env map.

    void FindSHCoeffsFromHDRCubeMap(const cImageData48* image, int n, Vec4f* coeffs);
    ///< Extract coeffs from 16-bit cube map.

    void CalcAtmosphereZH
    (
        int          numPhases,
        const Vec4f* colourPhases,   // (r, g, b, phase)
        int          n,
        Vec4f*       zcoeffs
    );
    ///< Find ZH coeffs for given atmosphere model, which is just a sum of colour * HG(phase).

    void AddGroundBounce
    (
        int             numBands,
        Vec4f*          coeffs, 
        const Vec4f&    diffColour,
        const Vec4f&    specColour
    );
    ///< Add in ground bounce controlled by spec/diff. Z is assumed to be local up.


    struct cSphereLight
    {
        Vec4f   mColourAndIntensity;    ///< r, g, b, intensity
        Vec4f   mPositionAndSize;       ///< location, size
    };

    void AddSphereLighting(Vec3f pos, float scale, int numLights, const cSphereLight* lights, Vec4f* coeffs);
    ///< Add in given sphere lights. 'coeffs' must be at least four bands (16 coefficients.)
    ///< Does not expect them to be sorted in importance order -- will choose the appropriate
    ///< approximation depending on light size and distance.

    void CreateSphereLightBRDF(int width, int height, uint32_t* pixels);
    ///< Creates a texture lookup table F(z, t) that returns diff/spec BRDF values given
    ///< z = dot(F_n, L_n) and t = r^2 / (r^2 + a^2). (F = fragment, L = light.)
}

#endif
