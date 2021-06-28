/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#include "liblib/mybase.h"
#include "libmycu/cucom/cucommon.h"
#include "libmycu/cucom/warpscan.cuh"
#include "libmycu/cupro/SerializedScoresSM.cuh"
#include "libmycu/cupro/SerializedScoresAttr.h"
#include "libmycu/cupro/SerializedDstMatchScores.cuh"
#include "libmycu/cupro/PM2DVectorFields.h"
#include "libmycu/cupro/CuBatchProcessing_com.cuh"
#include "CuBatchScoreMatrix_com.h"
#include "CuBatchScoreMatrix_ddms2s_v2.cuh"

// =========================================================================

// NOTE: parameters are passed to the device via constant memory and are 
// limited to 4 KB
// 
// device functions for computing batch score matrices;
// ddms2scores, serialized distance match scores and the map between 
// accumulated match scores and translated scores;
// dmsattr, attributes of ddms2scores
// ddmswgt, weight of scoring DDMS2S scores
// dstsegm, indivisible fragment of consecutive distance values
// nqyposs, number of query positions to process
// ndb1poss, number of cached db profile positions to process
// ndbCposs, number of new db profile positions to process
// dbxpad, number of padded positions for memory alignment
// querposoffset, offset from the origin of the device buffers allocated for 
// queries;
// bdb1posoffset, offset from the origin of the device buffers allocated for 
// cached db profile data;
// bdbCposoffset, offset from the origin of the device buffers allocated for 
// new (read) db profile data;
//

// #define CUSCO_DDMS2S_SMEMUnroll2xv2_TESTPRINT 1
// -------------------------------------------------------------------------
// CalcSM_DDMS2S_SMEMUnroll2xv2: device code for calculating pairwise 
// distance distribution match scores;
// matrix using shared memory and twofold unrolling along the x axis;
// each thread block processes two matrix blocks actually; 
// query profile data remains the same for these two spatially prallel 
// blocks;
// NOTE: it uses more shared memory but allows for greater occupancy;
// although bank conflicts cannot be avoided (due to random acces to the 
// SMEM), using SMEM garantees much greater efficiency in comparison to 
// other types of memory;
// NOTE: output pointers should be aligned!
// NOTE: results add to outscores and rewrite outmodscores;
// 
__global__ void CalcSM_DDMS2S_SMEMUnroll2xv2(
    CUBSM_TYPE* __restrict__ ddms2scores,
    SerializedDstMatchScoresAttr dmsattr,
    const float ddmswgt,
    const int dstsegm,
    uint nqyposs, uint ndb1poss, uint ndbCposs, uint dbxpad,
    uint querposoffset, uint bdb1posoffset, uint bdbCposoffset,
    CUBSM_TYPE* __restrict__ outscores,
    CUBSM_TYPE* __restrict__ outmodscores )
{
    extern __shared__ CUBSM_TYPE ddms2sCache[];//cached serialized scores
    //
    __shared__ FPTYPE 
            qrenoCache;//cache for query ENO
//             dbenoCache[2][SMINIT_2DCACHE_DIM];//cache for the ENOs of db profiles
    __shared__ INTYPE //cache for a tile of IR Dsts at query positions
            //NOTE: DDMS_DSTVAL_PACK_DEG-degree packing with +2 for one 32-bit 
            // cell for the length and one to account for the rounding effect
            qrdstCache[ptr2DNoDstsPerPos/DDMS_DSTVAL_PACK_DEG+2][SMINIT_2DCACHE_DIM];
    __shared__ INTYPE //cache for a tile of IR Dsts at db profile positions
            dbdstCache[ptr2DNoDstsPerPos/DDMS_DSTVAL_PACK_DEG+2][SMINIT_2DCACHE_DIM];
    __shared__ INTYPE //cache for a tile of IR Dsts at sec. db profile positions
            db2dstCache[ptr2DNoDstsPerPos/DDMS_DSTVAL_PACK_DEG+2][SMINIT_2DCACHE_DIM];
    //
    uint blockbeg_y = blockIdx.y * blockDim.y;
    uint row = blockbeg_y + threadIdx.y;
    const uint col = blockIdx.x * blockDim.x * 2 + threadIdx.x;//logical column
    const uint col2 = col + blockDim.x;//logical column
    //physical indices:
    uint db1pos;
    uint db1pos2;
    uint dbfldsndx;
    uint dbfldsndx2;
    if( col < ndb1poss ) {  db1pos = col + bdb1posoffset;
                            dbfldsndx = pmv2DTotFlds;
    } else {                db1pos = col - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx = TIMES2(pmv2DTotFlds);
    }
    if( col2 < ndb1poss ) { db1pos2 = col2 + bdb1posoffset;
                            dbfldsndx2 = pmv2DTotFlds;
    } else {                db1pos2 = col2 - ndb1poss + bdbCposoffset;//jump to section ndbCposs
                            dbfldsndx2 = TIMES2(pmv2DTotFlds);
    }
    //
    uint dbpronr;
    //
    //{{CACHE SERIALIZED SCORES using coalescent access (reuse registers)
    //NOTE: IMPORTANT: total number of entries (attr.ntotents_) is 
    // assumed to be less than the block size;
    //if this is not the case, uncomment the for loop for n-fold caching!
    dbpronr = threadIdx.y * blockDim.x + threadIdx.x;
    if( dbpronr < dmsattr.ntotents_ )
        ddms2sCache[dbpronr] = ddms2scores[dbpronr];
    //for( dbpronr += blockDim.x; dbpronr < dmsattr.ntotents_; dbpronr += blockDim.x )
    //    ddms2sCache[dbpronr] = ddms2scores[dbpronr];
    //for( ; dbpronr < ssattr.ntotents_ + cvattr.ntotents_; dbpronr += blockDim.x )
    //    ddms2sCache[dbpronr] = ddms2scores[dbpronr-dmsattr.ntotents_];
    //}}
    //



    //for comments, see the CalcSMInit... kernels
    //
    //read query ENO
    //NOTE: valid when processing one query at a time
    if( threadIdx.y < 1 && threadIdx.x < 1 ) {
        uint qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[blockbeg_y+querposoffset];
        //read only one element per block (blockDim.y x blockDim.x)
        qrenoCache = ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr];
    }

    //cache query IR distance values
    CacheIRDistanceValuesv2(
        qrdstCache, pmv2DDDvalues,
        blockbeg_y + querposoffset + threadIdx.x,
        blockbeg_y + threadIdx.x < nqyposs );

    //cache db profile IR distance values
    CacheIRDistanceValuesv2(
        dbdstCache, dbfldsndx + pmv2DDDvalues,
        db1pos,
        col < (ndb1poss + ndbCposs) );

    //cache second db profile's IR distance values
    CacheIRDistanceValuesv2(
        db2dstCache, dbfldsndx2 + pmv2DDDvalues,
        db1pos2,
        col2 < (ndb1poss + ndbCposs) );



    //the warp reads data written by other warps, sync
    __syncthreads();

    if( nqyposs <= row || (ndb1poss + ndbCposs) <= col )
        //NOTE: NO sync after the exit of some of the threads
        return;

//     //reuse registers
//     dbpronr = (col2 < (ndb1poss + ndbCposs));

    float score1;

    CalcDDMS2ScoreSMEMv2(
        ddms2sCache,
        dmsattr.card_,
        ddmswgt,
        dstsegm,
        qrenoCache,
        qrdstCache,
        dbdstCache,
        score1);

    row = row * (ndb1poss + ndbCposs + dbxpad);

    //perform coalescent write of scores
    //atomicAdd is faster than += when we need coalescent write performed once
    if(score1) {
        atomicAdd( &outscores[row + col], score1 );
        //outmodscores[row + col] = score1 - CONSTCVSSHIFT * ddmswgt;
    }

// //     if( dbpronr/*col2 < (ndb1poss + ndbCposs)*/)
    if(col2 < (ndb1poss + ndbCposs)) {
        CalcDDMS2ScoreSMEMv2(
            ddms2sCache,
            dmsattr.card_,
            ddmswgt,
            dstsegm,
            qrenoCache,
            qrdstCache,
            db2dstCache,
            score1);

// //     if( dbpronr/*col2 < (ndb1poss + ndbCposs)*/) {
        if(score1) {
            atomicAdd( &outscores[row + col2], score1 );
            //outmodscores[row + col2] = score1 - CONSTCVSSHIFT * cvswgt;
        }
// //     }
    }

#ifdef CUSCO_DDMS2S_SMEMUnroll2xv2_TESTPRINT
    if(!blockIdx.x && !blockIdx.y && !threadIdx.x && !threadIdx.y)
        ;//SerializedDstMatchScoresSM<CUBSM_TYPE>::PrintScores(ddms2sCache, dmsattr.card_);
    //
    int querypos = blockIdx.y * blockDim.y + threadIdx.y + querposoffset;
    dbpronr = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx+pmv2DAddrPro]))[db1pos];
    uint qpronr = ((INTYPE*)(dc_pm2dvfields_[pmv2DAddrPro]))[querypos];
    uint dbpronr2 = 99999;
    MYASSERT( qrenoCache == ((FPTYPE*)(dc_pm2dvfields_[pps2DENO]))[qpronr], "Inconsistency.");
    int nqrdstvals = (int)qrdstCache[0][threadIdx.y];
    int ndbdstvals = (int)dbdstCache[0][threadIdx.x];
    int ndb2dstvals = -1;
    if( col2 < (ndb1poss + ndbCposs)) {
        dbpronr2 = ((INTYPE*)(dc_pm2dvfields_[dbfldsndx2+pmv2DAddrPro]))[db1pos2];
        ndb2dstvals = (int)db2dstCache[0][threadIdx.x];
    }
    if((dbpronr == CUSCO_DDMS2S_SMEMUnroll2xv2_TESTPRINT || dbpronr2 == CUSCO_DDMS2S_SMEMUnroll2xv2_TESTPRINT) &&
        blockIdx.x==14 && blockIdx.y==1 && threadIdx.x==0 && threadIdx.y==0) {
        printf("\n Query %d (ENO %.1f):\n",qpronr,qrenoCache);
        for(int q=0; q < blockDim.y; q++) {
            nqrdstvals = (int)qrdstCache[0][q];
            printf(" %4d D: %d", querypos+q+1, nqrdstvals);
            for( int i = 0; i < nqrdstvals; i += DDMS_DSTVAL_PACK_DEG ) 
                printf(" %d %d %d %d",
                    qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]&0xff, (qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>8)&0xff,
                    (qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>16)&0xff, (qrdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>24)&0xff);
            printf("\n");
        }
        if(dbpronr == CUSCO_DDMS2S_SMEMUnroll2xv2_TESTPRINT) {
            printf(" Db pro %d:\n", dbpronr);
            for(int q=0; q < blockDim.x; q++) {
                ndbdstvals = (int)dbdstCache[0][q];
                printf(" %4d D: %d", db1pos+q+1, ndbdstvals);
                for( int i = 0; i < ndbdstvals; i += DDMS_DSTVAL_PACK_DEG ) 
                    printf(" %d %d %d %d",
                        dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]&0xff, (dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>8)&0xff,
                        (dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>16)&0xff, (dbdstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>24)&0xff);
                printf("\n");
            }
        } else {
            printf(" Db2 pro %d:\n", dbpronr2);
            for(int q=0; q < blockDim.x; q++) {
                ndb2dstvals = (int)db2dstCache[0][q];
                printf(" %4d D: %d", db1pos2+q+1, ndb2dstvals);
                for( int i = 0; i < ndb2dstvals; i += DDMS_DSTVAL_PACK_DEG ) 
                    printf(" %d %d %d %d",
                        db2dstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]&0xff, (db2dstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>8)&0xff,
                        (db2dstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>16)&0xff, (db2dstCache[i/DDMS_DSTVAL_PACK_DEG+1][q]>>24)&0xff);
                printf("\n");
            }
        }
    }
#endif
}
