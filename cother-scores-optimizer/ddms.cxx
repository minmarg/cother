/***************************************************************************
 *   Copyright (C) 2020-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include <stdio.h>
#include <xmmintrin.h>
#include "ddms.h"

// arithmetics on bytes packed into a 32-bit integer [QB QE PB PE]; QB and 
// QE, the query beginning and end positions; PB and PE, the db profile 
// beginning and end positions
//
#define CSO_DDMS_DP_COORDS__Q_BE_xff__P_BE_xff 0xffffffff
#define CSO_DDMS_DP_COORDS__MAKE(QB, QE, PB, PE) (((QB&0xff)<<24)|((QE&0xff)<<16)|((PB&0xff)<<8)|(PE&0xff))

#define DO_PRAGMA(x) _Pragma (#x)
#define PRAGMA_UNROLL_IMPL(x) DO_PRAGMA("GCC unroll " #x)
#define PRAGMA_UNROLL(x) PRAGMA_UNROLL_IMPL(x)

#define FP0 0.0f

// #define CSO_DDMS_DP_TEST_DP2
// #define CSO_DDMS_DP_TEST_DP2_Reg

// =========================================================================
//
CSO_DDMS::CSO_DDMS(float theta, float absdiffexp, float avgdistexp)
:   cso_ddms_dscore_theta_(theta),
    cso_ddms_dscore_absdiffexp_(absdiffexp),
    cso_ddms_dscore_avgdistexp_(avgdistexp),
    dstable(
        CSO_DDMS_DScore_MaxDst,
        CSO_DDMS_DScore_MaxDst_inv,
        CSO_DDMS_DScore_Approx_Score00,
        CSO_DDMS_DScore_Approx_Score200,
        theta, absdiffexp, avgdistexp)
{
}
CSO_DDMS::~CSO_DDMS()
{
}



// =========================================================================
// Helper methods for registering max ddm score
//
inline
void AssignMaxDDMScore(
    float* __restrict__ maxscores,
    __m64* __restrict__ coords_maxscores,
    const int ndxto, const uint ndxfrom )
{
    maxscores[ndxto] = maxscores[ndxfrom];
    coords_maxscores[ndxto] = coords_maxscores[ndxfrom];
}

inline
void AssignMaxDDMScore(
    float* __restrict__ maxscores,
    __m64* __restrict__ coords_maxscores,
    const int ndxto, const float segsco, const __m64 segcoords )
{
    maxscores[ndxto] = segsco;
    coords_maxscores[ndxto] = segcoords;
}

// CoordsIntersect: return true if coordinates coords1 and coords2 
// intersect; coordinates are packed into a 64-bit integer (upper 32 unused)
// [QB QE PB PE]; QB and QE, the query beginning and end positions; 
// PB and PE, the db profile beginning and end positions
inline
bool CoordsIntersect(__m64 coords1, __m64 coords2)
{
    //coords1 = [QB QE PB PE]; coords2 = [QB' QE' PB' PE'];
    //target condition for intersection is:
    // min(QE,QE')>=max(QB,QB') OR min(PE,PE')>=max(PB,PB')
    __m64 r1 = _mm_max_pu8(coords1, coords2);
    __m64 r2 = _mm_min_pu8(coords1, coords2);
    return 
        (_mm_extract_pi16(r2, 1) & 0xff) >= (((_mm_extract_pi16(r1, 1)>>8) & 0xff)) ||
        (_mm_extract_pi16(r2, 0) & 0xff) >= (((_mm_extract_pi16(r1, 0)>>8) & 0xff));
//         ((r2>>16) & 0xff) >= ((r1>>24) & 0xff) ||
//         (r2 & 0xff) >= ((r1>>8) & 0xff);
}

// -------------------------------------------------------------------------
// RegisterMaxDDMScore: register distance distribution match score of an 
// ungapped region from two distance distributions;
// maxscores, array of maximum-scoring segments recorded so far
// coords_maxscores, coordinates of the maximum-scoring segments
// segsco, score of the segment under consideration
// segcoords, coordinates of the segment under consideration
// NOTE: arrays are indexed by constant indexing (loops manually unrolled)
//
// #define APPROX_RegisterMaxDDMScore

inline
void RegisterMaxDDMScore(
    float* __restrict__ maxscores,
    __m64* __restrict__ coords_maxscores,
    const float segsco, const __m64 segcoords )
{
    if(!segsco)
        return;

#ifdef CSO_DDMS_DP_TEST_DP2_Reg
    printf(" b New score: %.3f %04x %04x    Buffer:", segsco,
            _mm_extract_pi16(segcoords, 1) & 0xffff,
            _mm_extract_pi16(segcoords, 0) & 0xffff);
    for(int i = 0; i < CSO_DDMS_DP_NSEGM; i++)
        printf("  %.3f %04x %04x", maxscores[i],
            _mm_extract_pi16(coords_maxscores[i], 1) & 0xffff,
            _mm_extract_pi16(coords_maxscores[i], 0) & 0xffff);
    printf("\n");
#endif

    //register score
    if(maxscores[0] <= segsco) {
        AssignMaxDDMScore(maxscores, coords_maxscores, 2, 1);
        AssignMaxDDMScore(maxscores, coords_maxscores, 1, 0);
        AssignMaxDDMScore(maxscores, coords_maxscores, 0, segsco, segcoords);
#ifndef APPROX_RegisterMaxDDMScore
        //invalidate intersecting segments for the case of segsco being the highest
        if(maxscores[1] && CoordsIntersect(coords_maxscores[0], coords_maxscores[1]))
            maxscores[1] = FP0;
        if(maxscores[2] && CoordsIntersect(coords_maxscores[0], coords_maxscores[2]))
            maxscores[2] = FP0;
        if(maxscores[1] == FP0) {
            AssignMaxDDMScore(maxscores, coords_maxscores, 1, 2);
            maxscores[2] = FP0;
        }
#endif
    } else if(maxscores[1] <= segsco) {
#ifndef APPROX_RegisterMaxDDMScore
        //if intersects with the regions of a higher score, exit
        if(CoordsIntersect(segcoords, coords_maxscores[0]))
            return;
#endif
        AssignMaxDDMScore(maxscores, coords_maxscores, 2, 1);
        AssignMaxDDMScore(maxscores, coords_maxscores, 1, segsco, segcoords);
#ifndef APPROX_RegisterMaxDDMScore
        //invalidate intersecting segments
        if(CoordsIntersect(coords_maxscores[1], coords_maxscores[2]))
            maxscores[2] = FP0;
#endif
    } else if(maxscores[2] <= segsco) {
#ifndef APPROX_RegisterMaxDDMScore
        //if intersects with the regions of a higher score, exit
        if(CoordsIntersect(segcoords, coords_maxscores[0]) ||
            CoordsIntersect(segcoords, coords_maxscores[1]))
            return;
#endif
        AssignMaxDDMScore(maxscores, coords_maxscores, 2, segsco, segcoords);
    } else
        return;

#ifdef APPROX_RegisterMaxDDMScore
    //invalidate intersecting segments of lower score
    if(maxscores[1] && CoordsIntersect(coords_maxscores[0], coords_maxscores[1]))
        maxscores[1] = FP0;
    if(maxscores[2] && CoordsIntersect(coords_maxscores[0], coords_maxscores[2]))
        maxscores[2] = FP0;
    if(maxscores[1] && CoordsIntersect(coords_maxscores[1], coords_maxscores[2]))
        maxscores[2] = FP0;
    if(!maxscores[1]) {
        AssignMaxDDMScore(maxscores, coords_maxscores, 1, 2);
        maxscores[2] = FP0;
    }
#endif

#ifdef CSO_DDMS_DP_TEST_DP2_Reg
    printf(" A New score: %.3f %08x    Buffer:",segsco,segcoords);
    for(int i = 0; i < CSO_DDMS_DP_NSEGM; i++)
        printf("  %.3f %04x %04x", maxscores[i],
            _mm_extract_pi16(coords_maxscores[i], 1) & 0xffff,
            _mm_extract_pi16(coords_maxscores[i], 0) & 0xffff);
    printf("\n");
#endif
}

// =========================================================================
// UpdateSegScoreUpperh2: update working variables while processing one DP 
// diagonal in the upper part of the DP matrix
//
template<int scndx>
inline
void CSO_DDMS::UpdateSegScoreUpperh2(
    float scij,
    int i, int j,
    int& segqb, int& segpb,
    int& im, int& jim, int& ib, int& jib,
    float& segsco, float& scom)
{
    if(!segsco) {
        segqb = i+scndx;
        segpb = j+i+scndx;
    }
    segsco = mymax(segsco+scij, FP0);
    if(segsco) {
        if(scom < segsco && GetOptionDSTSEGM() <= i+scndx - segqb + 1) {
            scom = segsco;
            im = i+scndx; jim = j + i+scndx;
            ib = segqb; jib = segpb;
        }
    }
}

// -------------------------------------------------------------------------
// UpdateSegScoreLowerh2: update working variables while processing one DP 
// diagonal in the lower part of the DP matrix
//
template<int scndx>
inline
void CSO_DDMS::UpdateSegScoreLowerh2(
    float scij,
    int i, int j,
    int& segqb, int& segpb,
    int& ijm, int& jm, int& ijb, int& jb,
    float& segsco, float& scom)
{
    if(!segsco) {
        segqb = i+j+scndx;
        segpb = j+scndx;
    }
    segsco = mymax(segsco+scij, FP0);
    if(segsco) {
        if(scom < segsco && GetOptionDSTSEGM() <= j+scndx - segpb + 1) {
            scom = segsco;
            ijm = i+j+scndx; jm = j + scndx;
            ijb = segqb; jb = segpb;
        }
    }
}

// =========================================================================
// ReadDstMatchScoresh2Unrolled: compiler-based unrolling (through 
// recursion) for reading data and calculating relative scores
//
template<int N>
inline
void CSO_DDMS::ReadDstMatchScoresh2Unrolled(
    const int* __restrict__ vdst1,
    const int * __restrict__ vdst2,
    float* __restrict__ scij, int i, int j)
{
    scij[N+2] = GetPosDstMatchScoreh2(vdst1[i+N+2+1], vdst2[j+N+2+1]);
    ReadDstMatchScoresh2Unrolled<N+2>(vdst1, vdst2, scij, i, j);
#if CSO_DDMS_DP_HEURISTIC == 1
    scij[N+1] = 0.5f*(scij[N] + scij[N+2]);
#else
    scij[N+1] = GetPosDstMatchScoreh2(vdst1[i+N+1+1], vdst2[j+N+1+1]);
#endif
}
template<>
inline
void CSO_DDMS::ReadDstMatchScoresh2Unrolled<CSO_DDMS_UNROLLING_DEGREE>(
    const int* __restrict__,
    const int* __restrict__,
    float* __restrict__, int, int)
{
}
// UpdateSegScoreUpperh2Unrolled: compiler-based unrolling for updating a
// maximal score and related indices when processing diagonals in the upper 
// part of the DP matrix
//
template<int N>
inline
void CSO_DDMS::UpdateSegScoreUpperh2Unrolled(
    float* __restrict__ scij,
    int i, int j,
    int& segqb, int& segpb,
    int& im, int& jim, int& ib, int& jib,
    float& segsco, float& scom)
{
    UpdateSegScoreUpperh2<N+1>(scij[N+1], i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
    UpdateSegScoreUpperh2Unrolled<N+1>(scij, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
}
template<>
inline
void CSO_DDMS::UpdateSegScoreUpperh2Unrolled<CSO_DDMS_UNROLLING_DEGREE>(
    float* __restrict__, int, int, int&, int&,
    int&, int&, int&, int&, float&, float&)
{
}
// UpdateSegScoreLowerh2Unrolled: compiler-based unrolling for updating a
// maximal score and related indices when processing diagonals in the lower 
// part of the DP matrix
//
template<int N>
inline
void CSO_DDMS::UpdateSegScoreLowerh2Unrolled(
    float* __restrict__ scij,
    int i, int j,
    int& segqb, int& segpb,
    int& ijm, int& jm, int& ijb, int& jb,
    float& segsco, float& scom)
{
    UpdateSegScoreLowerh2<N+1>(scij[N+1], i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
    UpdateSegScoreLowerh2Unrolled<N+1>(scij, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
}
template<>
inline
void CSO_DDMS::UpdateSegScoreLowerh2Unrolled<CSO_DDMS_UNROLLING_DEGREE>(
    float* __restrict__, int, int, int&, int&,
    int&, int&, int&, int&, float&, float&)
{
}
// -------------------------------------------------------------------------
// ProcessDiagonalUpperh2: process (find a maximum score in) one diagonal in 
// the upper part of the DP matrix
//
inline
void CSO_DDMS::ProcessDiagonalUpperh2(
    const int* vdst1, const int* vdst2,
    const int len1, const int len2,
    const int j,
    float* __restrict__ maxscores,
    __m64* __restrict__ coords_maxscores)
{
    float segsco;//score of the current segment
    int segqb, segpb;//query and db profile beginning positions of a segment

    segqb = 0;
    segpb = j;
    float scij[CSO_DDMS_UNROLLING_DEGREE+1];
    scij[0] = GetPosDstMatchScoreh2(vdst1[1], vdst2[j+1]);
    segsco = mymax(scij[0], FP0);
    for(int i = 0; i < len1 && j+i < len2; i += CSO_DDMS_UNROLLING_DEGREE) {
        float scom = FP0;
        int im, jim, ib, jib;

        ReadDstMatchScoresh2Unrolled<0>(vdst1, vdst2, scij, i, j+i);

        UpdateSegScoreUpperh2Unrolled<0>(scij, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);

        RegisterMaxDDMScore(
            maxscores, coords_maxscores,
            scom,
            _mm_set_pi8(0,0,0,0,ib&0xff,im&0xff,jib&0xff,jim&0xff));//CSO_DDMS_DP_COORDS__MAKE(ib, im, jib, jim));

        scij[0] = scij[CSO_DDMS_UNROLLING_DEGREE];
    }
}

// -------------------------------------------------------------------------
// ProcessDiagonalLowerh2: process (find a maximum score in) one diagonal in 
// the lower part of the DP matrix
//
inline
void CSO_DDMS::ProcessDiagonalLowerh2(
    const int* vdst1, const int* vdst2,
    const int len1, const int len2,
    const int i,
    float* __restrict__ maxscores,
    __m64* __restrict__ coords_maxscores)
{
    float segsco;//score of the current segment
    int segqb, segpb;//query and db profile beginning positions of a segment

    segqb = i;
    segpb = 0;
    float scij[CSO_DDMS_UNROLLING_DEGREE+1];
    scij[0] = GetPosDstMatchScoreh2(vdst1[i+1], vdst2[1]);
    segsco = mymax(scij[0], FP0);
    for(int j = 0; j < len2 && i+j < len1; j += CSO_DDMS_UNROLLING_DEGREE) {
        float scom = FP0;
        int ijm, jm, ijb, jb;

        ReadDstMatchScoresh2Unrolled<0>(vdst1, vdst2, scij, i+j, j);

        UpdateSegScoreLowerh2Unrolled<0>(scij, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);

        RegisterMaxDDMScore(
            maxscores, coords_maxscores,
            scom,
            _mm_set_pi8(0,0,0,0,ijb&0xff,ijm&0xff,jb&0xff,jm&0xff));//DDMS_COORDS__MAKE(ijb, ijm, jb, jm));

        scij[0] = scij[CSO_DDMS_UNROLLING_DEGREE];
    }
}

// =========================================================================
// CalculateDDMS: calculate distance distribution match score given two 
// distributions by vdst1 and vdst2; write the result to ddms;
//
float CSO_DDMS::CalculateDDMS(int* vdst1, int* vdst2)
{
    //no check of addresses
    int len1 = *vdst1, len2 = *vdst2;
    float score = FP0;

    if((len1 >> GetDPBandFraction()) < GetOptionDSTSEGM() || 
        (len2 >> GetDPBandFraction()) < GetOptionDSTSEGM())
        return GetIgnoreCode();

    float maxscores[CSO_DDMS_DP_NSEGM];//max-scoring segments
    //coordinate positions of the segments:
    // beginning (QB) and end (QE) positions for position 1;
    // beginning (PB) and end (PE) positions for position 2: [QB QE PB PE]
    __m64 coords_maxscores[CSO_DDMS_DP_NSEGM];

    //initialize
#ifdef __GNUC__
    PRAGMA_UNROLL(CSO_DDMS_DP_NSEGM)
#endif
    for(int i = 0; i < CSO_DDMS_DP_NSEGM; i++)
        maxscores[i] = FP0;
#ifdef __GNUC__
    PRAGMA_UNROLL(CSO_DDMS_DP_NSEGM)
#endif
    for(int i = 0; i < CSO_DDMS_DP_NSEGM; i++)
        coords_maxscores[i] = _mm_set_pi8(0,0,0,0,0xff,0xff,0xff,0xff);//CSO_DDMS_DP_COORDS__Q_BE_xff__P_BE_xff


    //Local alignment
    //NOTE: split dynamic programming into two parts;
    // (for recently identified high-scoring fragments do not 
    // remove form the queue non-overlapping other high-scoring fragments):
    for(int j = len2 >> GetDPBandFraction(); 0 <= j; j -= GetDPStep()) {
        ProcessDiagonalUpperh2(vdst1, vdst2, len1, len2, j, maxscores, coords_maxscores);
    }

    // the second part of the DP corresponds to traversing from the 
    // bottom-left corner up to the main diagonal:
    for(int i = len1 >> GetDPBandFraction(); 0 < i; i -= GetDPStep()) {
        ProcessDiagonalLowerh2(vdst1, vdst2, len1, len2, i, maxscores, coords_maxscores);
    }


#ifdef __GNUC__
    PRAGMA_UNROLL(CSO_DDMS_DP_NSEGM)
#endif
    for(int i = 0; i < CSO_DDMS_DP_NSEGM; i++)
        score += maxscores[i];

#ifdef CSO_DDMS_DP_TEST_DP2
    printf("\n           ");
    for(int i = 0; i < len1; i++)
        printf(" %2d",i);
    printf("\n\n D[qry %2d]:", len1);
    for(int i = 0; i < len1; i++)
        printf(" %2d",vdst1[i+1]);
    printf("\n\n D[dbp %2d]:", len2);
    for(int i = 0; i < len2; i++)
        printf(" %2d",vdst2[i+1]);
    printf("\n\n Top 3 scores:");
    for(int i = 0; i < CSO_DDMS_DP_NSEGM; i++)
        if(maxscores[i])
            printf("  %.3f (%2u %2u %2u %2u)",
                maxscores[i],
                (_mm_extract_pi16(coords_maxscores[i], 1)>>8) & 0xff,_mm_extract_pi16(coords_maxscores[i], 1) & 0xff,
                (_mm_extract_pi16(coords_maxscores[i], 0)>>8) & 0xff,_mm_extract_pi16(coords_maxscores[i], 0) & 0xff);
    printf("\n");
    for(int i = 0; i < CSO_DDMS_DP_NSEGM; i++)
        if(maxscores[i]) {
            for(int j = (_mm_extract_pi16(coords_maxscores[i], 1)>>8) & 0xff;
                j <= (_mm_extract_pi16(coords_maxscores[i], 1) & 0xff); j++)
                printf(" %2d",vdst1[j+1]);
            printf("\n");
            for(int j = (_mm_extract_pi16(coords_maxscores[i], 0)>>8) & 0xff;
                j <= (_mm_extract_pi16(coords_maxscores[i], 0) & 0xff); j++)
                printf(" %2d",vdst2[j+1]);
            printf("\n\n");
        }
    printf("\n\n");
#endif

    return rintf(score * GetDScoreGranularity());
}
