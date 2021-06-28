/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_ddms2s_h__
#define __CuBatchScoreMatrix_ddms2s_h__

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/warpscan.cuh"
#include "libmycu/cupro/SerializedDstMatchScores.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_com.h"

#if MAX_N_IR_DISTANCE_VALUES_PLUS_LENGTH > 127
#   error   "MAX_N_IR_DISTANCE_VALUES_PLUS_LENGTH cannot be greater than 127."
#endif

//number of ungapped maximum-scoring segments to consider in dynamic 
// programming when calculating the distance distribution match score 
// between two positions
#define DDMS2S_DP_NSEGM 3
//initial score for a segment
//#define DDMS2S_INIT_SEGSCO -1.f

#if DDMS2S_DP_NSEGM != 3
#   error   "DDMS2S_DP_NSEGM <> 3: Modify RegisterMaxDDMScore appropriately wrt the new value."
#endif

// #define CUSCO_DDMS2S_SMEMUnroll2x_TEST_DP2
// #define CUSCO_DDMS2S_SMEMUnroll2x_TEST_DP2_Reg

//device functions for computing distance distribution match scores

//kernel declarations

__global__ void CalcSM_DDMS2S_SMEMUnroll2x(
    CUBSM_TYPE* __restrict__ ddms2scores,
    SerializedDstMatchScoresAttr dmsattr,
    const float ddmswgt,
    const int dstsegm,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores );

// -------------------------------------------------------------------------
// CacheIRDistanceValues: cache IR distance values for each of the 
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
void CacheIRDistanceValues( 
    FPTYPE (* __restrict__ dstCache)[SMINIT_2DCACHE_DIM],
    int sectionoffset,
    uint profileposition,
    bool validposition )
{
    __shared__ int dstLens[SMINIT_2DCACHE_DIM];//lengths of distance values for SMINIT_2DCACHE_DIM positions

    if( validposition ) 
    {
        //number of distances is variable; 
        //first, cache the first 32 for each position
        if(threadIdx.y < ptr2DNoDstsPerPos) {
            dstCache[threadIdx.y][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[sectionoffset+threadIdx.y]))[
                    profileposition
                ];
        }
    }

    __syncthreads();//block sync

    FPTYPE nmaxdstvals = 0.0f;
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
        //continue caching distance values
        for(int ndxy = threadIdx.y+blockDim.y; 
                ndxy < dstLens[threadIdx.x]+1; ndxy += blockDim.y)
        {
            dstCache[ndxy][threadIdx.x] = 
                ((FPTYPE*)(dc_pm2dvfields_[sectionoffset+ndxy]))[
                    profileposition
                ];
        }
    }
}

// =========================================================================
// arithmetics on bytes packed into a 32-bit integer [QB QE PB PE]; QB and 
// QE, the query beginning and end positions; PB and PE, the db profile 
// beginning and end positions
//
#define DDMS_COORDS__Q_BE_xff__P_BE_xff 0xffffffff
#define DDMS_COORDS__MAKE(QB, QE, PB, PE) (((QB&0xff)<<24)|((QE&0xff)<<16)|((PB&0xff)<<8)|(PE&0xff))
// //swap beg. with end: [QB QE PB PE] -> [QE QB PE PB]:
// #define DDMS_COORDS__SWAP(COORDS) (((COORDS<<8)&0xff00ff00) | ((COORDS>>8)&0x00ff00ff))

// -------------------------------------------------------------------------
// GetPosDstMatchScore: get the distance matching score for a pair of 
// distance values d1 and d2 obtained at the query and db profile positions
//
__device__ __forceinline__
FPTYPE GetPosDstMatchScore(
    const CUBSM_TYPE* __restrict__ ddms2sCache, const int card,
    FPTYPE d1, FPTYPE d2)
{
    return SerializedDstMatchScoresSM<CUBSM_TYPE,DMS_APPROXIMATE_SCORES>::GetDstScore(ddms2sCache, card, d1, d2);
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// Helper methods for registering max ddm score
//
__device__ __forceinline__
void AssignMaxDDMScore(
    FPTYPE* __restrict__ maxscores, uint* __restrict__ coords_maxscores,
    const int ndxto, const uint ndxfrom )
{
    maxscores[ndxto] = maxscores[ndxfrom];
    coords_maxscores[ndxto] = coords_maxscores[ndxfrom];
}

__device__ __forceinline__
void AssignMaxDDMScore(
    FPTYPE* __restrict__ maxscores, uint* __restrict__ coords_maxscores,
    const int ndxto, const FPTYPE segsco, const uint segcoords )
{
    maxscores[ndxto] = segsco;
    coords_maxscores[ndxto] = segcoords;
}

// CoordsIntersect: return true if coordinates coords1 and coords2 
// intersect; coordinates are packed into a 32-bit integer 
// [QB QE PB PE]; QB and QE, the query beginning and end positions; 
// PB and PE, the db profile beginning and end positions
__device__ __forceinline__
bool CoordsIntersect(uint coords1, uint coords2)
{
    //coords1 = [QB QE PB PE]; coords2 = [QB' QE' PB' PE'];
    //target condition for intersection is:
    // min(QE,QE')>=max(QB,QB') OR min(PE,PE')>=max(PB,PB')
    uint r1 = __vmaxu4(coords1, coords2);
    uint r2 = __vminu4(coords1, coords2);
    return 
        ((r2>>16) & 0xff) >= ((r1>>24) & 0xff) ||
        (r2 & 0xff) >= ((r1>>8) & 0xff);
}

// -------------------------------------------------------------------------
// RegisterMaxDDMScore: register distance distribution match score of an 
// ungapped region from two distance distributions;
// maxscores, array of maximum-scoring segments recorded so far
// coords_maxscores, coordinates of the maximum-scoring segments
// segsco, score of the segment under consideration
// segcoords, coordinates of the segment under consideration
// NOTE: arrays are indexed by constant indexing to avoid using local memory 
// (loops manually unrolled)
//
// #define APPROX_RegisterMaxDDMScore

__device__ __forceinline__
void RegisterMaxDDMScore(
    FPTYPE* __restrict__ maxscores, uint* __restrict__ coords_maxscores,
    const FPTYPE segsco, const uint segcoords )
{
    if(!segsco)
        return;

#ifdef CUSCO_DDMS2S_SMEMUnroll2x_TEST_DP2_Reg
    if(blockIdx.x==14 && blockIdx.y==1 && threadIdx.x==3 && threadIdx.y==4) {
        printf(" b New score: %.3f %08x    Buffer:",segsco,segcoords);
        for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
            printf("  %.3f %08x",maxscores[i],coords_maxscores[i]);
        printf("\n");
    }
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

#ifdef CUSCO_DDMS2S_SMEMUnroll2x_TEST_DP2_Reg
    if(blockIdx.x==14 && blockIdx.y==1 && threadIdx.x==3 && threadIdx.y==4) {
        printf(" A New score: %.3f %08x    Buffer:",segsco,segcoords);
        for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
            printf("  %.3f %08x",maxscores[i],coords_maxscores[i]);
        printf("\n");
    }
#endif
}

// -------------------------------------------------------------------------
// CalcDDMS2ScoreSMEM: calculate the IR distance distribution match score 
// between one query position and one db profile position; 
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
__device__ __forceinline__
void CalcDDMS2ScoreSMEM( 
    const CUBSM_TYPE* __restrict__ ddms2sCache,
    const int card,
    const float ddmswgt,
    const int dstsegm,
    const FPTYPE qrenoCache,
    const FPTYPE (* __restrict__ qrdstCache)[SMINIT_2DCACHE_DIM],
    const FPTYPE (* __restrict__ dbdstCache)[SMINIT_2DCACHE_DIM],
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

    FPTYPE segsco;//score of the current segment
    int segqb, segpb;//query and db profile beginning positions of a segment

    //Local alignment
    //NOTE: split dynamic programming into two parts;
    // the first part corresponds to traversing from the upper-right corner
    // down to the main diagonal (this is done so that recently identified 
    // high-scoring fragments do not remove form the queue non-overlapping 
    // other high-scoring fragments):
    for(int j = len2-dstsegm; 0 <= j; j--) {
        segsco = FP0;
        segqb = 0;
        segpb = j;
        for(int i = 0; i < len1 && j+i < len2; i++) {
            FPTYPE scij = 
            GetPosDstMatchScore(
                ddms2sCache, card,
                qrdstCache[i+1][threadIdx.y], dbdstCache[j+i+1][threadIdx.x]);
            if(!segsco) {
                segqb = i;
                segpb = j+i;
            }
            segsco = myhdmax(segsco+scij, FP0);
            if(i - segqb + 1 < dstsegm)
                continue;
            RegisterMaxDDMScore(maxscores, coords_maxscores,
                segsco, DDMS_COORDS__MAKE(segqb, i, segpb, j+i));
        }
    }

    // the second part of the DP corresponds to traversing from the 
    // bottom-left corner up to the main diagonal:
    for(int i = len1-dstsegm; 0 < i; i--) {
        segsco = FP0;
        segqb = i;
        segpb = 0;
        for(int j = 0; j < len2 && i+j < len1; j++) {
            FPTYPE scij = 
            GetPosDstMatchScore(
                ddms2sCache, card,
                qrdstCache[i+j+1][threadIdx.y], dbdstCache[j+1][threadIdx.x]);
            if(!segsco) {
                segqb = i+j;
                segpb = j;
            }
            segsco = myhdmax(segsco+scij, FP0);
            if(j - segpb + 1 < dstsegm)
                continue;
            RegisterMaxDDMScore(maxscores, coords_maxscores,
                segsco, DDMS_COORDS__MAKE(segqb, i+j, segpb, j));
        }
    }

    #pragma unroll
    for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
        score1 += maxscores[i];

    score1 = CUSCO_DDMS2S_TRANS_SLOPE * score1 + CUSCO_DDMS2S_TRANS_INTERCEPT;
    score1 *= ddmswgt;

#ifdef CUSCO_DDMS2S_SMEMUnroll2x_TEST_DP2
    if(blockIdx.x==14 && blockIdx.y==1 && threadIdx.x==3 && threadIdx.y==4) {
        printf("\n           ", len1);
        for(int i = 0; i < len1; i++)
            printf(" %2d",i);
        printf("\n\n D[qry %2d]:", len1);
        for(int i = 0; i < len1; i++)
            printf(" %2.0f",qrdstCache[i+1][threadIdx.y]);
        printf("\n\n D[dbp %2d]:", len2);
        for(int i = 0; i < len2; i++)
            printf(" %2.0f",dbdstCache[i+1][threadIdx.x]);
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
                    printf(" %2.0f",qrdstCache[j+1][threadIdx.y]);
                printf("\n");
                for(int j = (coords_maxscores[i]>>8)&0xff; j <= (coords_maxscores[i]&0xff); j++)
                    printf(" %2.0f",dbdstCache[j+1][threadIdx.x]);
                printf("\n\n");
            }
        printf("\n\n");
    }
#endif
}

#endif//__CuBatchScoreMatrix_ddms2s_h__
