/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedDstMatchScores_h__
#define __SerializedDstMatchScores_h__

#include <stdlib.h>
#include <cuda_fp16.h>

#include "liblib/msg.h"
#include "liblib/mybase.h"
#include "libmycu/cucom/myassert.h"
#include "extsp/psl.h"

#include "SerializedDstMatchScoresAttr.h"

#define DMS_MaxDst 20 //max distance value considered; cardinality is 1 greater (includes 0)
#define DMS_MaxDst_inv 0.05f //inverse of DMS_MaxDst
#define DMS_Theta 0.2f //default theta parameter
#define DMS_AbsDiffExp 0.2f //default parameter for the exponent of the absolute difference
#define DMS_AvgDistExp 1.8f //default parameter for the exponent of the average distance

// #define DMS_Alpha 20.f //default alpha parameter
// #if DMS_Alpha <= 0
// #   error   "Invalid DMS_Alpha macro value."
// #endif



// =========================================================================
// class SerializedDstMatchScores
//
// Implementation of serialized distance match scores for parallel 
// processing; serialization transforms a 2D score table to a 1D 
// array
//
template <typename TScore, int APPROX>
class SerializedDstMatchScores
{
public:
    __host__ __device__ SerializedDstMatchScores(
        int szalloc, int card, int nelems )
    :   szalloc_(szalloc), card_(card), nelems_(nelems) 
    {};

    __host__ __device__ SerializedDstMatchScores();
    virtual __host__ __device__ ~SerializedDstMatchScores();

    virtual __host__ __device__ const TScore* GetScores() const = 0;
    __host__ __device__ int GetSizeAlloc() const { return szalloc_; }
    __host__ __device__ int GetCardinality() const { return card_; }
    __host__ __device__ int GetNEntries() const { return nelems_; }
    __host__ __device__ float GetScore( int row, int col) const;

protected:
    virtual __host__ __device__ TScore AccessScore(int ndx) const = 0;

protected:
    int szalloc_;//size (bytes) allocated for scores
    int card_;//cardinality determinator, # rows of a square score table
    int nelems_;//number of entries in one table
};

// -------------------------------------------------------------------------
// INLINES
//
// Default constructor
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
SerializedDstMatchScores<TScore,APPROX>::SerializedDstMatchScores()
:   szalloc_( 0 ),
    card_( 0 ),
    nelems_( 0 )
{
    CUMSG("SerializedDstMatchScores::SerializedDstMatchScores", 3);
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
SerializedDstMatchScores<TScore,APPROX>::~SerializedDstMatchScores()
{
    CUMSG("SerializedDstMatchScores::~SerializedDstMatchScores", 3);
}

// -------------------------------------------------------------------------
// GetScore: get score at table position indexed by row, col (0-based)
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
float SerializedDstMatchScores<TScore,APPROX>::GetScore(int row, int col) const
{
    CUMSG("SerializedDstMatchScores::GetScore", 5);
    MYASSERT(GetScores(), "SerializedDstMatchScores::GetScore: Memory access error.");
    MYASSERT(0 <= row && 0 <= col /*&& row < card_ && col < card_*/, 
        "SerializedDstMatchScores::GetScore: Invalid indices.");

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + ((row * (row-1)) >> 1);//index of the element now

    MYASSERT(col*(int)sizeof(TScore) < szalloc_,
        "SerializedDstMatchScores::GetScore: Index out of range.");

    return AccessScore(col);
}

// -------------------------------------------------------------------------
// SerializedDstMatchScores specializations --------------------------------

// -------------------------------------------------------------------------
// GetScore: get score at table position indexed by row, col (0-based)
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScores<float,DMS_AS_2>::GetScore(int row, int col) const
{
    CUMSG("SerializedDstMatchScores<float,2>::GetScore", 5);
    MYASSERT(0 <= row && 0 <= col /*&& row < card_ && col < card_*/, 
        "SerializedDstMatchScores<float,2>::GetScore: Invalid indices.");

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );

    MYASSERT(col < card_,
        "SerializedDstMatchScores<float,2>::GetScore: Index out of range.");

    float intpt = DMS_Approx2_Score00 - (float)(col) * DMS_Approx2_Score200 * DMS_MaxDst_inv;
    constexpr float slope = (DMS_Approx2_Score200 - DMS_Approx2_Score00) * DMS_MaxDst_inv;
    return (float)(row) * slope + intpt;
}

// -------------------------------------------------------------------------
// GetScore: get score at table position indexed by row, col (0-based)
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScores<__half2,DMS_AS_1>::GetScore(int row, int col) const
{
    CUMSG("SerializedDstMatchScores<__half2,1>::GetScore", 5);
    MYASSERT(GetScores(), "SerializedDstMatchScores<__half2,1>::GetScore: Memory access error.");
    MYASSERT(0 <= row && 0 <= col /*&& row < card_ && col < card_*/, 
        "SerializedDstMatchScores<__half2,1>::GetScore: Invalid indices.");

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );

    MYASSERT(col < card_,
        "SerializedDstMatchScores<__half2,1>::GetScore: Index out of range.");

    //__half2 si = scores[col];
    float2 sif = __half22float2(AccessScore(col));
    return (float)(row) * sif.x + sif.y;
}





// =========================================================================
// class SerializedDstMatchScoresSM
//
// Implementation of serialized distance match scores for parallel 
// processing (shared/general memory version); serialization transforms a 2D 
// score table to a 1D array
//
template <typename TScore, int APPROX>
class SerializedDstMatchScoresSM: public SerializedDstMatchScores<TScore,APPROX>
{
public:
    __host__ __device__ SerializedDstMatchScoresSM(
        TScore* scores, int szalloc, int card, int nelems )
    :   SerializedDstMatchScores<TScore,APPROX>(szalloc, card, nelems),
        h_scores_(scores)
    {};

    __host__ __device__ SerializedDstMatchScoresSM();
    virtual __host__ __device__ ~SerializedDstMatchScoresSM();

    virtual __host__ __device__ const TScore* GetScores() const { return h_scores_; }

    __host__ __device__ float GetScore(int row, int col) const
    {
        return SerializedDstMatchScores<TScore,APPROX>::GetScore(row, col);
    }

    static __host__ __device__ float GetDstScore(const TScore* scores, int card, TScore d1, TScore d2);
    static __host__ __device__ float GetDstScore(const TScore* scores, int card, int d1, int d2);
    static __host__ __device__ float GetScore(const TScore* scores, int row, int col);

    static __host__ __device__ void PrintScores(const TScore* scores, int card);

protected:
    virtual __host__ __device__ TScore AccessScore(int ndx) const { return h_scores_[ndx]; }

protected:
    //{{host/device data
    TScore* h_scores_;//multi-dimensional scores
    //}}
};

// -------------------------------------------------------------------------
// INLINES
//
// Default constructor
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
SerializedDstMatchScoresSM<TScore,APPROX>::SerializedDstMatchScoresSM()
:   SerializedDstMatchScores<TScore,APPROX>(),
    h_scores_( NULL )
{
    CUMSG("SerializedDstMatchScoresSM::SerializedDstMatchScoresSM", 3);
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
SerializedDstMatchScoresSM<TScore,APPROX>::~SerializedDstMatchScoresSM()
{
    CUMSG("SerializedDstMatchScoresSM::~SerializedDstMatchScoresSM", 3);
}

// -------------------------------------------------------------------------
// GetScore: get score at table position identified by row, col; 
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
float SerializedDstMatchScoresSM<TScore,APPROX>::GetScore( 
    const TScore* scores, int row, int col )
{
    CUMSG("SerializedDstMatchScoresSM::GetScore [static]", 5);
    //MYASSERT( scores, "SerializedDstMatchScoresSM::GetScore: Memory access error.");
    //MYASSERT( 0 <= row && 0 <= col, 
    //    "SerializedDstMatchScoresSM::GetScore: Invalid indices.");

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    if( row )
        col += row + ((row * (row-1)) >> 1);//index of the element now

    return scores[col];
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
float SerializedDstMatchScoresSM<TScore,APPROX>::GetDstScore(
    const TScore* scores, int card, TScore d1, TScore d2)
{
    CUMSG("SerializedDstMatchScoresSM::GetDstScore [static]", 5);
    int row = (int)d1;//rintf(d1);
    int col = (int)d2;//rintf(d2);
    row = SLC_MIN(row, card-1);
    col = SLC_MIN(col, card-1);
    return SerializedDstMatchScoresSM<TScore,APPROX>::GetScore(scores, row, col);
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
float SerializedDstMatchScoresSM<TScore,APPROX>::GetDstScore(
    const TScore* scores, int /*card*/, int d1, int d2)
{
    CUMSG("SerializedDstMatchScoresSM::GetDstScore [static int]", 5);
    if(d1 == 0xff || d2 == 0xff)
        return TScore(-1.f);
//     d1 = SLC_MIN(d1, card-1);
//     d2 = SLC_MIN(d2, card-1);
    return SerializedDstMatchScoresSM<TScore,APPROX>::GetScore(scores, d1, d2);
}

// -------------------------------------------------------------------------
// PrintScores: print all scores for testing
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
void SerializedDstMatchScoresSM<TScore,APPROX>::PrintScores(
    const TScore* scores, int card)
{
    CUMSG("SerializedDstMatchScoresSM::PrintScores [static]", 5);
    printf("\nSerializedDstMatchScoresSM:\n");
    for(int n=0, i=0; i < card; i++) {
        printf(" %2d:", i);
        for(int j=0; j <= i; j++)
            printf(" %.3f", scores[n++]);
        printf("\n");
    }
    printf("\n");
}


// -------------------------------------------------------------------------
// SerializedDstMatchScoresSM specializations ------------------------------

// -------------------------------------------------------------------------
// GetScore: get score by two-level linear interpolation given distance 
// values row, col. 
// The scores are calculated on the fly with the minimum number of 
// instructions; let n=DMS_MaxDst, then score at (n,n) is always 0; the 
// scores at the score matrix positions k1=s_(0,0)>0 and k2=s_(n,0)<0 are 
// given. 
// Other scores are calculated by two level interpolation: the first phase 
// of linear interpolation applies to scores (i,i) and (n,i), where i=1..n. 
// The second phase of linear interpolation applies to every column of 
// scores (i,j) separately for given j and i=j..n, with scores (j,j) and 
// (n,j) calculated in the first phase.
// The slope a and intercept b for a column have closed-form expressions:
// 
// a = (s_(j,j) - s_(n,j)) / (j-n) = ((k1-k1*j/n) - (k2-k2*j/n)) / (j-n) 
//   = (j-n)/n 1/(j-n) (k2-k1); 
// a = (k2-k1) / n (does not depend on j!);
// 
// b = (j s_(n,j) - n s_(j,j)) / (j-n) 
//   = (j (k2-k2*j/n) - n (k1-k1*j/n)) / (j-n)
//   = (j (k2+k1) - k2 j^2/n - n k1) / (j-n)
// The roots of the numerator is n and n k1/ k2; given the roots, the 
// multiplicative factor is found to be -k2 / n. Hence,
// b = (-k2/n) (j-n) (j - n k1 /k2) / (j-n) = (-k2/n) (j - n k1 /k2), or
// b = k1 - j k2 / n
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScoresSM<float,DMS_AS_2>::GetScore( 
    const float* /*scores*/, int row, int col )
{
    CUMSG("SerializedDstMatchScoresSM<float,2>::GetScore [static]", 5);
    //get the smaller distance value (select between row and col)
    CNDSWAP( int, row<col, row, col );
    float intpt = DMS_Approx2_Score00 - (float)(col) * DMS_Approx2_Score200 * DMS_MaxDst_inv;
    constexpr float slope = (DMS_Approx2_Score200 - DMS_Approx2_Score00) * DMS_MaxDst_inv;
    return (float)(row) * slope + intpt;
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScoresSM<float,DMS_AS_2>::GetDstScore(
    const float* /*scores*/, int /*card*/, float d1, float d2)
{
    CUMSG("SerializedDstMatchScoresSM<float,2>::GetDstScore [static]", 5);
    if(d1 == 255.f || d2 == 255.f)
        return (-1.f);
    CNDSWAP(float, d1<d2, d1, d2);
    float intpt = DMS_Approx2_Score00 - d2 * DMS_Approx2_Score200 * DMS_MaxDst_inv;
    constexpr float slope = (DMS_Approx2_Score200 - DMS_Approx2_Score00) * DMS_MaxDst_inv;
    return d1 * slope + intpt;
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScoresSM<float,DMS_AS_2>::GetDstScore(
    const float* scores, int /*card*/, int d1, int d2)
{
    CUMSG("SerializedDstMatchScoresSM<float,2>::GetDstScore [static int]", 5);
    if(d1 == 0xff || d2 == 0xff)
        return (-1.f);
//     d1 = (d1 < card)? d1: card-1;//SLC_MIN(d1, card-1);
//     d2 = (d2 < card)? d2: card-1;//SLC_MIN(d2, card-1);
    return SerializedDstMatchScoresSM<float,DMS_AS_2>::GetScore(scores, d1, d2);
}

// -------------------------------------------------------------------------
// PrintScores: print all scores for testing
// 
template <>
__host__ __device__
inline
void SerializedDstMatchScoresSM<float,DMS_AS_2>::PrintScores(
    const float* /*scores*/, int card)
{
    CUMSG("SerializedDstMatchScoresSM<float,2>::PrintScores [static]", 5);
    printf("\nSerializedDstMatchScoresSM<float,2>:\n");
    for(int i=0; i < card; i++) {
        printf(" %2d:", i);
        for(int j=0; j <= i; j++) {
            float intpt = DMS_Approx2_Score00 - (float)(j) * DMS_Approx2_Score200 * DMS_MaxDst_inv;
            constexpr float slope = (DMS_Approx2_Score200 - DMS_Approx2_Score00) * DMS_MaxDst_inv;
            printf(" %.3f", (float)(i) * slope + intpt);
        }
        printf("\n");
    }
    printf("\n");
}


// =========================================================================
// GetScore: get score by linear interpolation given distance 
// values row, col; 
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScoresSM<__half2,DMS_AS_1>::GetScore( 
    const __half2* scores, int row, int col )
{
    CUMSG("SerializedDstMatchScoresSM<__half2,1>::GetScore [static]", 5);

    //calculate index within a score table
    CNDSWAP( int, row<col, row, col );
    //__half2 si = scores[col];
    float2 sif = __half22float2(scores[col]);
    return (float)(row) * sif.x + sif.y;
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScoresSM<__half2,DMS_AS_1>::GetDstScore(
    const __half2* scores, int card, __half2 /*d1*/, __half2 /*d2*/)
{
    CUMSG("SerializedDstMatchScoresSM<__half2,1>::GetDstScore [static]", 5);
    MYASSERT(0, "SerializedDstMatchScoresSM<__half2,1>::GetDstScore: Cannot be called.");
    return SerializedDstMatchScoresSM<__half2,DMS_AS_1>::GetScore(scores, card-1, card-1);
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <>
__host__ __device__
inline
float SerializedDstMatchScoresSM<__half2,DMS_AS_1>::GetDstScore(
    const __half2* scores, int /*card*/, int d1, int d2)
{
    CUMSG("SerializedDstMatchScoresSM<__half2,1>::GetDstScore [static int]", 5);
    if(d1 == 0xff || d2 == 0xff)
        return (-1.f);
//     d1 = (d1 < card)? d1: card-1;//SLC_MIN(d1, card-1);
//     d2 = (d2 < card)? d2: card-1;//SLC_MIN(d2, card-1);
    return SerializedDstMatchScoresSM<__half2,DMS_AS_1>::GetScore(scores, d1, d2);
}

// -------------------------------------------------------------------------
// PrintScores: print all scores for testing
// 
template <>
__host__ __device__
inline
void SerializedDstMatchScoresSM<__half2,DMS_AS_1>::PrintScores(
    const __half2* scores, int card)
{
    CUMSG("SerializedDstMatchScoresSM<__half2,1>::PrintScores [static]", 5);
    printf("\nSerializedDstMatchScoresSM<__half2,1>:\n");
    for(int i=0; i < card; i++) {
        printf(" %2d:", i);
        for(int j=0; j <= i; j++) {
            __half2 si = scores[j];
            float2 sif = __half22float2(si);
            printf(" %.3f", (float)(i) * sif.x + sif.y);
        }
        printf("\n");
    }
    printf("\n");
}





// =========================================================================
// class SerializedScoresCtor
//
// Implementation of construction of serialized distance match scores for 
// parallel processing; serialization transforms a 2D score table to a 1D 
// array
//
template <typename TScore, int APPROX>
class SerializedDstMatchScoresCtor: public SerializedDstMatchScoresSM<TScore,APPROX>
{
    using SerializedDstMatchScoresSM<TScore,APPROX>::h_scores_;
    using SerializedDstMatchScoresSM<TScore,APPROX>::szalloc_;
    using SerializedDstMatchScoresSM<TScore,APPROX>::card_;
    using SerializedDstMatchScoresSM<TScore,APPROX>::nelems_;

public:
    SerializedDstMatchScoresCtor(int card = DMS_MaxDst+1);

    virtual __host__ __device__ ~SerializedDstMatchScoresCtor();

protected:
    void NewScores(int sztotal);
    float CalculateScore(int d1, int d2);
    __host__ __device__ void DestroyScores();
};

// -------------------------------------------------------------------------
// INLINES
//
// Destructor
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
SerializedDstMatchScoresCtor<TScore,APPROX>::~SerializedDstMatchScoresCtor()
{
    CUMSG("SerializedDstMatchScoresCtor::~SerializedDstMatchScoresCtor", 3);
    DestroyScores();
}

// -------------------------------------------------------------------------
// DestroyScores: destroy scores
// 
template <typename TScore, int APPROX>
__host__ __device__
inline
void SerializedDstMatchScoresCtor<TScore,APPROX>::DestroyScores()
{
    CUMSG("SerializedDstMatchScoresCtor::DestroyScores", 3);
    if( h_scores_ ) {
        free(h_scores_);
        h_scores_ = NULL;
    }
    szalloc_ = 0;
}

// -------------------------------------------------------------------------
// NewScores: allocate memory for scores
// 
template <typename TScore,int APPROX>
inline
void SerializedDstMatchScoresCtor<TScore,APPROX>::NewScores(int sztotal)
{
    CUMSG("SerializedDstMatchScoresCtor::NewScores", 3);
    DestroyScores();
    h_scores_ = (TScore*)malloc(sztotal);
    if( !h_scores_ )
        throw MYRUNTIME_ERROR("SerializedDstMatchScoresCtor::NewScores: Not enough memory.");
    szalloc_ = sztotal;
}

// -------------------------------------------------------------------------
// CalculateScore: calculate score given a pair of distance values d1 and d2
//
template <typename TScore,int APPROX>
inline
float SerializedDstMatchScoresCtor<TScore,APPROX>::CalculateScore(int d1, int d2)
{
    float davg = 0.5f * (float)(d1+d2);
    float dabs = (float)abs(d1-d2);
    float sco = DMS_Theta;
    if(davg)
        sco = 
            ////(DMS_Theta - (float)abs(d1-d2)/davg) * expf(-davg*davg / DMS_Alpha));
            (DMS_Theta - dabs/davg) * 
                //the below represents the score weight; the first factor 
                //implies scaling a (negative) score increasingly with 
                //increasing absolute difference; the second factor 
                //downscales the score as the average distance increases:
                (1.f + powf(dabs/(float)DMS_MaxDst, DMS_AbsDiffExp)) *
                (1.f - powf(davg/(float)DMS_MaxDst, DMS_AvgDistExp));
    return sco;
}

// -------------------------------------------------------------------------
// constructor: calculate and serialize 2D distance match scores
//
template <typename TScore, int APPROX>
inline
SerializedDstMatchScoresCtor<TScore,APPROX>::SerializedDstMatchScoresCtor(int card)
:   SerializedDstMatchScoresSM<TScore,APPROX>()
{
    CUMSG("SerializedDstMatchScoresCtor::SerializedDstMatchScoresCtor", 3);

    card_ = card;

    nelems_ = (card_ * (card_ + 1)) >> 1;//number of entries in the table

    //size calculated for scores:
    int szscores = nelems_ * sizeof(TScore);

    NewScores(szscores);

    int n = 0;

    //write all score tables
    for(int d1 = 0; d1 < card_; d1++)
        for(int d2 = 0; d2 <= d1; d2++) {
            *(h_scores_ + n++) = (TScore)CalculateScore(d1, d2);
        }
}

// -------------------------------------------------------------------------
// SerializedDstMatchScoresCtor specializations ----------------------------

// -------------------------------------------------------------------------
// constructor: specialization for two-level score interpolation: no 
// calculations
//
template <>
inline
SerializedDstMatchScoresCtor<float,DMS_AS_2>::SerializedDstMatchScoresCtor(int card)
:   SerializedDstMatchScoresSM<float,DMS_AS_2>()
{
    CUMSG("SerializedDstMatchScoresCtor<float,2>::SerializedDstMatchScoresCtor<float,2>", 3);

    card_ = card;

    //number of entries in the table is 0; no precalculation
    nelems_ = 0;
}

// -------------------------------------------------------------------------
// constructor: specialization for type half2: calculate and serialize 2D 
// distance match scores
//
template <>
inline
SerializedDstMatchScoresCtor<__half2,DMS_AS_1>::SerializedDstMatchScoresCtor(int card)
:   SerializedDstMatchScoresSM<__half2,DMS_AS_1>()
{
    CUMSG("SerializedDstMatchScoresCtor<__half2,1>::SerializedDstMatchScoresCtor<__half2,1>", 3);

    card_ = card;

    //number of entries in the table is card, i.e., packed slopes and intercepts
    nelems_ = card_;

    //size calculated for scores:
    int szscores = nelems_ * sizeof(__half2);

    NewScores(szscores);

    int n = 0;

    //write all score tables
    for(int d1 = 0; d1 < card_; d1++) {
        int d21 = d1, d22 = card_-1;
        float sco1 = CalculateScore(d1, d21);
        float sco2 = CalculateScore(d1, d22);
        float slope = 0.0f;
        float intpt = 0.0f;
        if(d1 - d22) {
            //slope for the linear apporoximation of column d1: y1-y2/(x1-x2)
            //intercept: x1y2-x2y1/(x1-x2)
            slope = (sco1 - sco2) / (float)(d1-d22);
            intpt = ((float)(d1)*sco2 - (float)(d22)*sco1) / (float)(d1-d22);
        }
        __half2 si = __floats2half2_rn(slope, intpt);
        *(h_scores_ + n++) = si;
    }
}

#endif//__SerializedScoresCtor_h__
