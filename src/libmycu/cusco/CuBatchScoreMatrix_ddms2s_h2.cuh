/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_ddms2s_h2_h__
#define __CuBatchScoreMatrix_ddms2s_h2_h__

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/warpscan.cuh"
#include "libmycu/cupro/SerializedDstMatchScores.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_ddms2s.cuh"
#include "CuBatchScoreMatrix_com.h"

//constant shift for the DDMS scoring function during realignment
#define CONSTDDMSSHIFT  0.2f //when DDMSWGT=.3; 0.1f if less

// #define CUSCO_DDMS2S_SMEMUnroll2x_TEST_DP2_h2

//device functions for computing distance distribution match scores

//kernel declarations

template<int DSTEP, bool HEURISTIC, int DSTBNDF>
__global__ void CalcSM_DDMS2S_SMEMUnroll2x_h2(
    DMS_TYPE* __restrict__ ddms2scores,
    SerializedDstMatchScoresAttr dmsattr,
    const float ddmswgt,
    const int dstsegm,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores );

// -------------------------------------------------------------------------
// unrolling degree used in DP heuristic
#define DDMS_h2_UNROLL_DEG 4

// =========================================================================
// CacheIRDistanceValuesh2: cache IR distance values for each of the 
// SMINIT_2DCACHE_DIM positions of either query or db profile;
// ndxlens, index for 2D array dstLens;
// dstCache, address of the SMEM block for distances themselves;
// offsetDsts, offset for dstCache along the first axis;
// sectionoffset, offset of a section which can represent query of db 
// profile data;
// profileposition, profile position to be read (thread idx included);
// validposition, flag of whether position profileposition is valid and is 
// not outside allowed boundaries;
//
__device__ __forceinline__
void CacheIRDistanceValuesh2( 
    INTYPE (* __restrict__ dstCache)[SMINIT_2DCACHE_DIM],
    int sectionoffset,
    uint profileposition,
    bool validposition )
{
    __shared__ int dstLens[SMINIT_2DCACHE_DIM];//lengths of distance values for SMINIT_2DCACHE_DIM positions
    auto lfSaturate = [](INTYPE dst){return (dst < DMS_MaxDst || dst == 0xff)? dst: DMS_MaxDst;};

    if( validposition ) 
    {
        //number of distances is variable; 
        //first, cache the first 32 for each position
        if(threadIdx.y < ptr2DNoDstsPerPos) {
            dstCache[threadIdx.y][threadIdx.x] = 
                threadIdx.y? lfSaturate((INTYPE)
                    ((FPTYPE*)(dc_pm2dvfields_[sectionoffset+threadIdx.y]))[
                        profileposition
                    ]): ((INTYPE)
                    ((FPTYPE*)(dc_pm2dvfields_[sectionoffset+threadIdx.y]))[
                        profileposition
                    ]);
        }
    }

    __syncthreads();//block sync

    INTYPE nmaxdstvals = 0;
    if(threadIdx.y == 0) {
        if( validposition ) 
            nmaxdstvals = dstCache[0][threadIdx.x];
        nmaxdstvals = mywarpreducemax(nmaxdstvals);//max across the warp
        nmaxdstvals = __shfl_sync(0xffffffff, nmaxdstvals, 0);//populate across the warp
        dstLens[threadIdx.x] = (int)nmaxdstvals;//write to SMEM for the block threads' availability
    }

    __syncthreads();//block sync

    if( validposition ) 
    {
        int ndxy = threadIdx.y+blockDim.y; 
        //continue caching distance values
        for(; ndxy < dstCache[0]/*dstLens*/[threadIdx.x]+1; ndxy += blockDim.y)
        {
            dstCache[ndxy][threadIdx.x] = lfSaturate((INTYPE)
                ((FPTYPE*)(dc_pm2dvfields_[sectionoffset+ndxy]))[
                    profileposition
                ]);
        }

        if(ndxy < dstLens[threadIdx.x]+1+DDMS_h2_UNROLL_DEG)
            dstCache[ndxy][threadIdx.x] = 0xff;
    }
}



// =========================================================================
// GetPosDstMatchScoreh2: get the distance matching scores for a pair of 
// distance values packed in d1 and d2, obtained at the query and db profile 
// positions
//
__device__ __forceinline__
FPTYPE GetPosDstMatchScoreh2(
    const DMS_TYPE* __restrict__ ddms2sCache, const int card,
    INTYPE d1, INTYPE d2)
{
    return SerializedDstMatchScoresSM<DMS_TYPE,DMS_APPROXIMATE_SCORES>::GetDstScore(ddms2sCache, card, d1, d2);
}



// =========================================================================
// UpdateSegScoreUpperh2: update working variables while processing one DP 
// diagonal in the upper part of the DP matrix
//
template<int scndx>
__device__ __forceinline__
void UpdateSegScoreUpperh2(
    const int dstsegm,
    FPTYPE scij,
    int i, int j,
    int& segqb, int& segpb,
    int& im, int& jim, int& ib, int& jib,
    FPTYPE& segsco, FPTYPE& scom)
{
    if(!segsco) {
        segqb = i+scndx;
        segpb = j+i+scndx;
    }
    segsco = myhdmax(segsco+scij, FP0);
    if(segsco) {
        if(scom < segsco && dstsegm <= i+scndx - segqb + 1) {
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
__device__ __forceinline__
void UpdateSegScoreLowerh2(
    const int dstsegm,
    FPTYPE scij,
    int i, int j,
    int& segqb, int& segpb,
    int& ijm, int& jm, int& ijb, int& jb,
    FPTYPE& segsco, FPTYPE& scom)
{
    if(!segsco) {
        segqb = i+j+scndx;
        segpb = j+scndx;
    }
    segsco = myhdmax(segsco+scij, FP0);
    if(segsco) {
        if(scom < segsco && dstsegm <= j+scndx - segpb + 1) {
            scom = segsco;
            ijm = i+j+scndx; jm = j + scndx;
            ijb = segqb; jb = segpb;
        }
    }
}



// =========================================================================
// ProcessDiagonalUpperh2: process (find a maximum score) one diagonal in 
// the upper part of the DP matrix;
// Helper functions follow:
//
// ReadDstMatchScoresh2Unrolled: compiler-based unrolling (through 
// recursion) for reading data from SMEM and calculating relative scores
//
template<int N, bool HEURISTIC>
__device__ __forceinline__
void ReadDstMatchScoresh2Unrolled(
    const DMS_TYPE* __restrict__ ddms2sCache, const int card,
    const INTYPE (* __restrict__ qrdstCache)[SMINIT_2DCACHE_DIM],
    const INTYPE (* __restrict__ dbdstCache)[SMINIT_2DCACHE_DIM],
    FPTYPE* __restrict__ scij, int i, int j)
{
    scij[N+2] = GetPosDstMatchScoreh2(ddms2sCache, card,
        qrdstCache[i+N+2+1][threadIdx.y], dbdstCache[j+N+2+1][threadIdx.x]);
    ReadDstMatchScoresh2Unrolled<N+2,HEURISTIC>(ddms2sCache, card, qrdstCache, dbdstCache, scij, i, j);
    if(HEURISTIC)
        scij[N+1] = 0.5f*(scij[N] + scij[N+2]);
    else
        scij[N+1] = GetPosDstMatchScoreh2(ddms2sCache, card,
            qrdstCache[i+N+1+1][threadIdx.y], dbdstCache[j+N+1+1][threadIdx.x]);
}
template<>
__device__ __forceinline__
void ReadDstMatchScoresh2Unrolled<DDMS_h2_UNROLL_DEG,true>(
    const DMS_TYPE* __restrict__, const int,
    const INTYPE (*__restrict__)[SMINIT_2DCACHE_DIM],
    const INTYPE (*__restrict__)[SMINIT_2DCACHE_DIM],
    FPTYPE* __restrict__, int, int)
{
}
template<>
__device__ __forceinline__
void ReadDstMatchScoresh2Unrolled<DDMS_h2_UNROLL_DEG,false>(
    const DMS_TYPE* __restrict__, const int,
    const INTYPE (*__restrict__)[SMINIT_2DCACHE_DIM],
    const INTYPE (*__restrict__)[SMINIT_2DCACHE_DIM],
    FPTYPE* __restrict__, int, int)
{
}
// UpdateSegScoreUpperh2Unrolled: compiler-based unrolling for updating a
// maximal score and related indices when processing diagonals in the upper 
// part of the DP matrix
//
template<int N>
__device__ __forceinline__
void UpdateSegScoreUpperh2Unrolled(
    const int dstsegm,
    FPTYPE* __restrict__ scij,
    int i, int j,
    int& segqb, int& segpb,
    int& im, int& jim, int& ib, int& jib,
    FPTYPE& segsco, FPTYPE& scom)
{
    UpdateSegScoreUpperh2<N+1>(dstsegm, scij[N+1], i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
    UpdateSegScoreUpperh2Unrolled<N+1>(dstsegm, scij, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
}
template<>
__device__ __forceinline__
void UpdateSegScoreUpperh2Unrolled<DDMS_h2_UNROLL_DEG>(
    const int, FPTYPE* __restrict__, int, int, int&, int&,
    int&, int&, int&, int&, FPTYPE&, FPTYPE&)
{
}
// UpdateSegScoreLowerh2Unrolled: compiler-based unrolling for updating a
// maximal score and related indices when processing diagonals in the lower 
// part of the DP matrix
//
template<int N>
__device__ __forceinline__
void UpdateSegScoreLowerh2Unrolled(
    const int dstsegm,
    FPTYPE* __restrict__ scij,
    int i, int j,
    int& segqb, int& segpb,
    int& ijm, int& jm, int& ijb, int& jb,
    FPTYPE& segsco, FPTYPE& scom)
{
    UpdateSegScoreLowerh2<N+1>(dstsegm, scij[N+1], i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
    UpdateSegScoreLowerh2Unrolled<N+1>(dstsegm, scij, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
}
template<>
__device__ __forceinline__
void UpdateSegScoreLowerh2Unrolled<DDMS_h2_UNROLL_DEG>(
    const int, FPTYPE* __restrict__, int, int, int&, int&,
    int&, int&, int&, int&, FPTYPE&, FPTYPE&)
{
}
// -------------------------------------------------------------------------
// ProcessDiagonalUpperh2: process (find a maximum score) one diagonal in 
// the upper part of the DP matrix
//
template<bool HEURISTIC>
__device__ __forceinline__
void ProcessDiagonalUpperh2(
    const DMS_TYPE* __restrict__ ddms2sCache,
    const int card,
    const INTYPE (* __restrict__ qrdstCache)[SMINIT_2DCACHE_DIM],
    const INTYPE (* __restrict__ dbdstCache)[SMINIT_2DCACHE_DIM],
    const int dstsegm,
    const int len1, const int len2,
    const int j,
    FPTYPE* __restrict__ maxscores,
    uint* __restrict__ coords_maxscores)
{
    FPTYPE segsco;//score of the current segment
    int segqb, segpb;//query and db profile beginning positions of a segment

    segqb = 0;
    segpb = j;
    FPTYPE scij[DDMS_h2_UNROLL_DEG+1];
//     FPTYPE scij0 = 
    scij[0] = GetPosDstMatchScoreh2(ddms2sCache, card,
        qrdstCache[1][threadIdx.y], dbdstCache[j+1][threadIdx.x]);
    segsco = myhdmax(scij[0], FP0);
    for(int i = 0; i < len1 && j+i < len2; i += DDMS_h2_UNROLL_DEG) {
        FPTYPE scom = FP0;
        int im, jim, ib, jib;

        ReadDstMatchScoresh2Unrolled<0,HEURISTIC>(
            ddms2sCache, card, qrdstCache, dbdstCache,
            scij, i, j+i);

//         //manual unrolling for the case DDMS_h2_UNROLL_DEG==4
//         FPTYPE scij2 = GetPosDstMatchScoreh2(ddms2sCache, card,
//             qrdstCache[i+2+1][threadIdx.y], dbdstCache[j+i+2+1][threadIdx.x]);
//         FPTYPE scij4 = GetPosDstMatchScoreh2(ddms2sCache, card,
//             qrdstCache[i+4+1][threadIdx.y], dbdstCache[j+i+4+1][threadIdx.x]);
// 
//         FPTYPE scij1, scij3;
//         if(HEURISTIC) {
//             scij1 = 0.5f*(scij0 + scij2);
//             scij3 = 0.5f*(scij2 + scij4);
//         } else {
//             scij1 = GetPosDstMatchScoreh2(ddms2sCache, card,
//                 qrdstCache[i+1+1][threadIdx.y], dbdstCache[j+i+1+1][threadIdx.x]);
//             scij3 = GetPosDstMatchScoreh2(ddms2sCache, card,
//                 qrdstCache[i+3+1][threadIdx.y], dbdstCache[j+i+3+1][threadIdx.x]);
//         }

        UpdateSegScoreUpperh2Unrolled<0>(dstsegm, scij, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);

//         UpdateSegScoreUpperh2<1>(dstsegm, scij1, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
//         UpdateSegScoreUpperh2<2>(dstsegm, scij2, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
//         UpdateSegScoreUpperh2<3>(dstsegm, scij3, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
//         UpdateSegScoreUpperh2<4>(dstsegm, scij4, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);

        RegisterMaxDDMScore(maxscores, coords_maxscores,
            scom, DDMS_COORDS__MAKE(ib, im, jib, jim));

        scij[0] = scij[DDMS_h2_UNROLL_DEG];
    }
}

// -------------------------------------------------------------------------
// ProcessDiagonalLowerh2: process (find a maximum score) one diagonal in 
// the lower part of the DP matrix
//
template<bool HEURISTIC>
__device__ __forceinline__
void ProcessDiagonalLowerh2(
    const DMS_TYPE* __restrict__ ddms2sCache,
    const int card,
    const INTYPE (* __restrict__ qrdstCache)[SMINIT_2DCACHE_DIM],
    const INTYPE (* __restrict__ dbdstCache)[SMINIT_2DCACHE_DIM],
    const int dstsegm,
    const int len1, const int len2,
    const int i,
    FPTYPE* __restrict__ maxscores,
    uint* __restrict__ coords_maxscores)
{
    FPTYPE segsco;//score of the current segment
    int segqb, segpb;//query and db profile beginning positions of a segment

    segqb = i;
    segpb = 0;
    FPTYPE scij[DDMS_h2_UNROLL_DEG+1];
//     FPTYPE scij0 = 
    scij[0] = GetPosDstMatchScoreh2(ddms2sCache, card,
        qrdstCache[i+1][threadIdx.y], dbdstCache[1][threadIdx.x]);
    segsco = myhdmax(scij[0], FP0);
    for(int j = 0; j < len2 && i+j < len1; j += DDMS_h2_UNROLL_DEG) {
        FPTYPE scom = FP0;
        int ijm, jm, ijb, jb;

        ReadDstMatchScoresh2Unrolled<0,HEURISTIC>(
            ddms2sCache, card, qrdstCache, dbdstCache,
            scij, i+j, j);

//         //manual unrolling for the case DDMS_h2_UNROLL_DEG==4
//         FPTYPE scij2 = GetPosDstMatchScoreh2(ddms2sCache, card,
//             qrdstCache[i+j+2+1][threadIdx.y], dbdstCache[j+2+1][threadIdx.x]);
//         FPTYPE scij4 = GetPosDstMatchScoreh2(ddms2sCache, card,
//             qrdstCache[i+j+4+1][threadIdx.y], dbdstCache[j+4+1][threadIdx.x]);
// 
//         FPTYPE scij1, scij3;
//         if(HEURISTIC) {
//             scij1 = 0.5f*(scij0 + scij2);
//             scij3 = 0.5f*(scij2 + scij4);
//         } else {
//             scij1 = GetPosDstMatchScoreh2(ddms2sCache, card,
//                 qrdstCache[i+j+1+1][threadIdx.y], dbdstCache[j+1+1][threadIdx.x]);
//             scij3 = GetPosDstMatchScoreh2(ddms2sCache, card,
//                 qrdstCache[i+j+3+1][threadIdx.y], dbdstCache[j+3+1][threadIdx.x]);
//         }

        UpdateSegScoreLowerh2Unrolled<0>(dstsegm, scij, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);

//         UpdateSegScoreLowerh2<1>(dstsegm, scij1, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
//         UpdateSegScoreLowerh2<2>(dstsegm, scij2, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
//         UpdateSegScoreLowerh2<3>(dstsegm, scij3, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
//         UpdateSegScoreLowerh2<4>(dstsegm, scij4, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);

        RegisterMaxDDMScore(maxscores, coords_maxscores,
            scom, DDMS_COORDS__MAKE(ijb, ijm, jb, jm));

        scij[0] = scij[DDMS_h2_UNROLL_DEG];
    }
}



// =========================================================================
// CalcDDMS2ScoreSMEMh2: calculate the IR distance distribution match score 
// between one query position and one db profile position, using the 
// heuristic of linear inpterpolation for half of DP cells; 
// DSTEP, template argument, which determines the sparsity of processing the 
// DSTBNDF, template argument, fraction of distances size the band 
// occupies in the DP matrix; NOTE: expressed as x to 1/2^x
// DP matrix for a given position (step over diagonals);
// NOTE: The score corresponds to up to DDMS2S_DP_NSEGM(three) 
// maximum-scoring non-intersecting ungapped segments found by dynamic 
// programming; to secure a valid number of registers, 2x unrolling is used 
// outside the function;
// NOTE: sync-free score calculation for a pair of positions;
// ddms2sCache, SMEM cache of distance distribution match scores;
// card, cardinality of the distance match score table;
// ddmswgt, weight of DD match scores;
// dstsegm, parameter of the minimum length of an indivisible alignment 
// segment;
// qrenoCache, query ENO;
// qrdstCache, cache of query distances for each position considered;
// dbdstCache, cache of db profile distances for each position considered;
// score1, DDMS2S score for the position under consideration;
//
template<int DSTEP, bool HEURISTIC, int DSTBNDF>
__device__ __forceinline__
void CalcDDMS2ScoreSMEMh2( 
    const DMS_TYPE* __restrict__ ddms2sCache,
    const int card,
    const float ddmswgt,
    const int dstsegm,
    const FPTYPE qrenoCache,
    const INTYPE (* __restrict__ qrdstCache)[SMINIT_2DCACHE_DIM],
    const INTYPE (* __restrict__ dbdstCache)[SMINIT_2DCACHE_DIM],
    float& score1)
{
    score1 = 0.0f;

    int len1 = (int)qrdstCache[0][threadIdx.y];
    int len2 = (int)dbdstCache[0][threadIdx.x];

    if(len1 < dstsegm || len2 < dstsegm)
        //NOTE: some threads exit, DO NOT sync!
        return;

    //NOTE: in most cases, len1<len2, so do not swap lengths and pointers

    FPTYPE maxscores[DDMS2S_DP_NSEGM];//max-scoring segments
    //coordinate positions of the segments:
    // the meaning of the bytes are described below in the order from the 
    // most-significant to the least significant byte:
    // beginning (QB) and end (QE) positions for the query;
    // beginning (PB) and end (PE) positions for the db profile: [QB QE PB PE]
    uint coords_maxscores[DDMS2S_DP_NSEGM];

    //initialize local variables
    #pragma unroll
    for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
        maxscores[i] = FP0;//DDMS2S_INIT_SEGSCO;
    #pragma unroll
    for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
        coords_maxscores[i] = DDMS_COORDS__Q_BE_xff__P_BE_xff;

    //Local alignment
    //NOTE: split dynamic programming into two parts;
    // the first part corresponds to traversing from the upper-right corner
    // down to the main diagonal (this is done so that recently identified 
    // high-scoring fragments do not remove form the queue non-overlapping 
    // other high-scoring fragments):
    for(int j = len2 >> DSTBNDF; 0 <= j; j -= DSTEP) {
        ProcessDiagonalUpperh2<HEURISTIC>(ddms2sCache, card,
            qrdstCache, dbdstCache, dstsegm,
            len1, len2, j, maxscores, coords_maxscores);
    }

    // the second part of the DP corresponds to traversing from the 
    // bottom-left corner up to the main diagonal:
    for(int i = len1 >> DSTBNDF; 0 < i; i -= DSTEP) {
        ProcessDiagonalLowerh2<HEURISTIC>(ddms2sCache, card,
            qrdstCache, dbdstCache, dstsegm,
            len1, len2, i, maxscores, coords_maxscores);
    }

    #pragma unroll
    for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
        score1 += maxscores[i];

    //scale the ddms score first; then, apply regression
    // (allow the compiler to multiply the constants)
    score1 = CUSCO_DDMS2S_TRANS_SLOPE * CUSCO_DDMS2S_TRANS_SCALE * score1 + 
            CUSCO_DDMS2S_TRANS_INTERCEPT;
    score1 *= ddmswgt;

#ifdef CUSCO_DDMS2S_SMEMUnroll2x_TEST_DP2_h2
    if(blockIdx.x==14 && blockIdx.y==1 && threadIdx.x==3 && threadIdx.y==4) {
        printf("\n           ");
        for(int i = 0; i < len1; i++)
            printf(" %2d",i);
        printf("\n\n D[qry %2d]:", len1);
        for(int i = 0; i < len1; i++)
            printf(" %2d",qrdstCache[i+1][threadIdx.y]);
        printf("\n\n D[dbp %2d]:", len2);
        for(int i = 0; i < len2; i++)
            printf(" %2d",dbdstCache[i+1][threadIdx.x]);
        printf("\n\n Top 3 scores:");
        for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
            if(maxscores[i])
                printf("  %.3f (%2u %2u %2u %2u)",
                    maxscores[i],
                    (coords_maxscores[i]>>24)&0xff,(coords_maxscores[i]>>16)&0xff,
                    (coords_maxscores[i]>>8)&0xff,coords_maxscores[i]&0xff);
        printf("\n");
        for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
            if(maxscores[i]) {
                for(int j = (coords_maxscores[i]>>24)&0xff; j <= ((coords_maxscores[i]>>16)&0xff); j++)
                    printf(" %2d",qrdstCache[j+1][threadIdx.y]);
                printf("\n");
                for(int j = (coords_maxscores[i]>>8)&0xff; j <= (coords_maxscores[i]&0xff); j++)
                    printf(" %2d",dbdstCache[j+1][threadIdx.x]);
                printf("\n\n");
            }
        printf("\n\n");
    }
#endif
}

// #include "CuBatchScoreMatrix_ddms2s_h2.cu"

#endif//__CuBatchScoreMatrix_ddms2s_h2_h__
