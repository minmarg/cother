/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __CuBatchScoreMatrix_ddms2s_v2_h__
#define __CuBatchScoreMatrix_ddms2s_v2_h__

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/warpscan.cuh"
#include "libmycu/cupro/SerializedDstMatchScores.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_ddms2s.cuh"
#include "CuBatchScoreMatrix_com.h"

// #define CUSCO_DDMS2S_SMEMUnroll2xv2_TEST_DP2v2

//device functions for computing distance distribution match scores

//kernel declarations

__global__ void CalcSM_DDMS2S_SMEMUnroll2xv2(
    CUBSM_TYPE* __restrict__ ddms2scores,
    SerializedDstMatchScoresAttr dmsattr,
    const float ddmswgt,
    const int dstsegm,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores );

// =========================================================================
// degree of packing distance values
#define DDMS_DSTVAL_PACK_DEG 4
#define DDMS_DSTVAL_PACK_DEG_LOG 2

// uninformative distance values packed into one 32-bit variuable
#define DDMS_DSTVAL_UNINF_PACK 0xffffffff

// NOTE: pack distance values in reverse order so that unpacking is more efficient
#define DDMS_DSTVAL__PACK(V1, V2, V3, V4) DDMS_COORDS__MAKE(V4, V3, V2, V1)

// -------------------------------------------------------------------------
// CacheIRDistanceValuesv2: cache IR distance values for each of the 
// SMINIT_2DCACHE_DIM positions of either query or db profile and pack 
// them in compressed format;
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
void CacheIRDistanceValuesv2( 
    INTYPE (* __restrict__ dstCache)[SMINIT_2DCACHE_DIM],
    int sectionoffset,
    uint profileposition,
    bool validposition )
{
    //maximum length of distance values for SMINIT_2DCACHE_DIM positions:
    __shared__ int dstLens[SMINIT_2DCACHE_DIM];
    __shared__ INTYPE //temporary cache for a tile of IR Dsts at profile positions
            tmpdstCache[ptr2DNoDstsPerPos][SMINIT_2DCACHE_DIM];

    if( validposition ) 
    {
        //number of distances is variable; 
        //first, cache the first 32 for each position
        if(threadIdx.y < ptr2DNoDstsPerPos) {
            tmpdstCache[threadIdx.y][threadIdx.x] = (INTYPE)
                ((FPTYPE*)(dc_pm2dvfields_[sectionoffset+threadIdx.y]))[
                    profileposition
                ];
        }
    }

    __syncthreads();//block sync

    INTYPE nmaxdstvals = 0.0f;
    if(threadIdx.y == 0) {
        if( validposition ) 
            nmaxdstvals = tmpdstCache[0][threadIdx.x];
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
            tmpdstCache[ndxy][threadIdx.x] = (INTYPE)
                ((FPTYPE*)(dc_pm2dvfields_[sectionoffset+ndxy]))[
                    profileposition
                ];
        }
    }

    __syncthreads();//block sync

    //finally, pack/compress distance values for each position
    //NOTE: Important: it is assumed that the size of packed data does not 
    // exceed SMINIT_2DCACHE_DIM
    if( validposition && 
        threadIdx.y < blockDim.y && 
        threadIdx.y * DDMS_DSTVAL_PACK_DEG < tmpdstCache[0][threadIdx.x] + DDMS_DSTVAL_PACK_DEG) 
    {
        if(threadIdx.y == 0)
            dstCache[threadIdx.y][threadIdx.x] = tmpdstCache[threadIdx.y][threadIdx.x];
        else {
            if((threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 3 < tmpdstCache[0][threadIdx.x])
                dstCache[threadIdx.y][threadIdx.x] = DDMS_DSTVAL__PACK(
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 1][threadIdx.x],
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 2][threadIdx.x],
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 3][threadIdx.x],
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 4][threadIdx.x]
                );
            else if((threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 2 < tmpdstCache[0][threadIdx.x])
                dstCache[threadIdx.y][threadIdx.x] = DDMS_DSTVAL__PACK(
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 1][threadIdx.x],
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 2][threadIdx.x],
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 3][threadIdx.x],
                    0xff
                );
            else if((threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 1 < tmpdstCache[0][threadIdx.x])
                dstCache[threadIdx.y][threadIdx.x] = DDMS_DSTVAL__PACK(
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 1][threadIdx.x],
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 2][threadIdx.x],
                    0xff, 0xff
                );
            else if((threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG < tmpdstCache[0][threadIdx.x])
                dstCache[threadIdx.y][threadIdx.x] = DDMS_DSTVAL__PACK(
                    tmpdstCache[(threadIdx.y-1) * DDMS_DSTVAL_PACK_DEG + 1][threadIdx.x],
                    0xff, 0xff, 0xff
                );
            else
                dstCache[threadIdx.y][threadIdx.x] = DDMS_DSTVAL_UNINF_PACK;
        }
    }
}

// -------------------------------------------------------------------------
// GetPosDstMatchScorev2: get the distance matching scores for four pairs of 
// distance values packed in d1 and d2, obtained at the query and db profile 
// positions
//
__device__ __forceinline__
void GetPosDstMatchScorev2(
    const CUBSM_TYPE* __restrict__ ddms2sCache, const int card,
    INTYPE d1, INTYPE d2,
    FPTYPE& sc1, FPTYPE& sc2, FPTYPE& sc3, FPTYPE& sc4)
{
    sc1 = SerializedDstMatchScoresSM<CUBSM_TYPE>::GetDstScore(ddms2sCache, card, d1 & 0xff, d2 & 0xff);
    d1 >>= 8; d2 >>= 8;
    sc2 = SerializedDstMatchScoresSM<CUBSM_TYPE>::GetDstScore(ddms2sCache, card, d1 & 0xff, d2 & 0xff);
    d1 >>= 8; d2 >>= 8;
    sc3 = SerializedDstMatchScoresSM<CUBSM_TYPE>::GetDstScore(ddms2sCache, card, d1 & 0xff, d2 & 0xff);
    d1 >>= 8; d2 >>= 8;
    sc4 = SerializedDstMatchScoresSM<CUBSM_TYPE>::GetDstScore(ddms2sCache, card, d1 & 0xff, d2 & 0xff);
}

// -------------------------------------------------------------------------
// UpdateSegScoreUpperv2: 
//
template<int scndx>
__device__ __forceinline__
void UpdateSegScoreUpperv2(
    const int dstsegm,
    FPTYPE scij,
    int i, int j,
    int& segqb, int& segpb,
    int& im, int& jim, int& ib, int& jib,
    FPTYPE& segsco, FPTYPE& scom)
{
    if(segsco + scij > FP0) {
        segsco += scij;
        if(scom < segsco && i+scndx - segqb + 1 < dstsegm) {
            scom = segsco;
            im = i+scndx; jim = j + i+scndx;
            ib = segqb; jib = segpb;
        }
    } else {
        segsco = FP0; 
        segqb = i+scndx+1; segpb = j+i+scndx+1;
    }
}

// -------------------------------------------------------------------------
// UpdateSegScoreLowerv2: 
//
template<int scndx>
__device__ __forceinline__
void UpdateSegScoreLowerv2(
    const int dstsegm,
    FPTYPE scij,
    int i, int j,
    int& segqb, int& segpb,
    int& ijm, int& jm, int& ijb, int& jb,
    FPTYPE& segsco, FPTYPE& scom)
{
    if(segsco + scij > FP0) {
        segsco += scij;
        if(scom < segsco && j+scndx - segpb + 1 < dstsegm) {
            scom = segsco;
            ijm = i+j+scndx; jm = j + scndx;
            ijb = segqb; jb = segpb;
        }
    } else {
        segsco = FP0;
        segqb = i+j+scndx+1; segpb = j+scndx+1;
    }
}

// -------------------------------------------------------------------------
// CalcDDMS2ScoreSMEMv2: calculate the IR distance distribution match score 
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
void CalcDDMS2ScoreSMEMv2( 
    const CUBSM_TYPE* __restrict__ ddms2sCache,
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
        for(int i = 0; i < len1 && j+i < len2; i += DDMS_DSTVAL_PACK_DEG) {
            FPTYPE scij1, scij2, scij3, scij4, scom = FP0;
            int im = i, jim = j + i, ib = segqb, jib = segpb;
            GetPosDstMatchScorev2(
                ddms2sCache, card,
                //the first cell is the length (+1)
                qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y],
                dbdstCache[(j+i)/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x],
                scij1, scij2, scij3, scij4);
            //if(!segsco) {segqb = i; segpb = j+i;}
            UpdateSegScoreUpperv2<0>(dstsegm, scij1, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
            UpdateSegScoreUpperv2<1>(dstsegm, scij2, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
            UpdateSegScoreUpperv2<2>(dstsegm, scij3, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);
            UpdateSegScoreUpperv2<3>(dstsegm, scij4, i, j, segqb, segpb, im, jim, ib, jib, segsco, scom);

            RegisterMaxDDMScore(maxscores, coords_maxscores,
                scom, DDMS_COORDS__MAKE(ib, im, jib, jim));
        }
    }

    // the second part of the DP corresponds to traversing from the 
    // bottom-left corner up to the main diagonal:
    for(int i = len1-dstsegm; 0 < i; i--) {
        segsco = FP0;
        segqb = i;
        segpb = 0;
        for(int j = 0; j < len2 && i+j < len1; j += DDMS_DSTVAL_PACK_DEG) {
            FPTYPE scij1, scij2, scij3, scij4, scom = FP0; 
            int ijm = i + j, jm = j, ijb = segqb, jb = segpb;
            GetPosDstMatchScorev2(
                ddms2sCache, card,
                //the first cell is the length (+1)
                qrdstCache[(i+j)/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y],
                dbdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x],
                scij1, scij2, scij3, scij4);
            //if(!segsco) {segqb = i+j; segpb = j;}
            UpdateSegScoreLowerv2<0>(dstsegm, scij1, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
            UpdateSegScoreLowerv2<1>(dstsegm, scij2, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
            UpdateSegScoreLowerv2<2>(dstsegm, scij3, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
            UpdateSegScoreLowerv2<3>(dstsegm, scij4, i, j, segqb, segpb, ijm, jm, ijb, jb, segsco, scom);
            RegisterMaxDDMScore(maxscores, coords_maxscores,
                scom, DDMS_COORDS__MAKE(ijb, ijm, jb, jm));
        }
    }

    #pragma unroll
    for(int i = 0; i < DDMS2S_DP_NSEGM; i++)
        score1 += maxscores[i];

    score1 *= ddmswgt;

#ifdef CUSCO_DDMS2S_SMEMUnroll2xv2_TEST_DP2v2
    if(blockIdx.x==14 && blockIdx.y==1 && threadIdx.x==3 && threadIdx.y==4) {
        printf("\n           ", len1);
        for(int i = 0; i < len1; i++)
            printf(" %2d",i);
        printf("\n\n D[qry %2d]:", len1);
        for(int i = 0; i < len1; i += DDMS_DSTVAL_PACK_DEG)
            printf(" %2d %2d %2d %2d",
                qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]&0xff,
                (qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]>>8)&0xff,
                (qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]>>16)&0xff,
                (qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]>>24)&0xff);
        printf("\n\n D[dbp %2d]:", len2);
        for(int i = 0; i < len2; i += DDMS_DSTVAL_PACK_DEG)
            printf(" %2d %2d %2d %2d",
                dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]&0xff,
                (dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]>>8)&0xff,
                (dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]>>16)&0xff,
                (dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]>>24)&0xff
            );
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
                for(int j = (coords_maxscores[i]>>24)&0xff; j <= ((coords_maxscores[i]>>16)&0xff); j+=DDMS_DSTVAL_PACK_DEG)
                    printf(" %2d %2d %2d %2d",
                        qrdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]&0xff,
                        (qrdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]>>8)&0xff,
                        (qrdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]>>16)&0xff,
                        (qrdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.y]>>24)&0xff);
                printf("\n");
                for(int j = (coords_maxscores[i]>>8)&0xff; j <= (coords_maxscores[i]&0xff); j+=DDMS_DSTVAL_PACK_DEG)
                    printf(" %2d %2d %2d %2d",
                           dbdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]&0xff,
                           (dbdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]>>8)&0xff,
                           (dbdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]>>16)&0xff,
                           (dbdstCache[j/DDMS_DSTVAL_PACK_DEG+1][threadIdx.x]>>24)&0xff);
                printf("\n\n");
            }
        printf("\n\n");
    }
#endif
}

#endif//__CuBatchScoreMatrix_ddms2s_v2_h__
