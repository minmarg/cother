/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __dstable_h__
#define __dstable_h__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>

#define MYSWAP(T,CND,a,b) if(CND) {T _my_tmp_=a; a=b; b=_my_tmp_;}
#define MYMIN( X, Y ) ((( X ) < ( Y ))? ( X ): ( Y ))

// =========================================================================
// class DstScoreTable
//
// Implementation of serialized distance match scores; serialization 
// transforms a 2D score table to a 1D array
//
template <typename TScore>
class DstScoreTable
{
public:
    DstScoreTable(
        int szalloc, int card, int nelems )
    :   szalloc_(szalloc), card_(card), nelems_(nelems) 
    {};

    DstScoreTable();
    virtual ~DstScoreTable();

    virtual const TScore* GetScores() const = 0;
    int GetSizeAlloc() const { return szalloc_; }
    int GetCardinality() const { return card_; }
    int GetNEntries() const { return nelems_; }
    TScore GetScore( int row, int col) const;

protected:
    virtual TScore AccessScore(int ndx) const = 0;

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
template <typename TScore>
inline
DstScoreTable<TScore>::DstScoreTable()
:   szalloc_( 0 ),
    card_( 0 ),
    nelems_( 0 )
{
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore>
inline
DstScoreTable<TScore>::~DstScoreTable()
{
}

// -------------------------------------------------------------------------
// GetScore: get score at table position indexed by row, col (0-based)
// 
template <typename TScore>
inline
TScore DstScoreTable<TScore>::GetScore(int row, int col) const
{
    //calculate index within a score table
    MYSWAP( int, row<col, row, col );

    if( row )
        col += row + ((row * (row-1)) >> 1);//index of the element now

    return AccessScore(col);
}



// =========================================================================
// class DstScoreTableSM
//
// Implementation of serialized distance match scores; serialization 
// transforms a 2D score table to a 1D array
//
template <typename TScore>
class DstScoreTableSM: public DstScoreTable<TScore>
{
protected:
    using DstScoreTable<TScore>::card_;

public:
    DstScoreTableSM(
        TScore* scores, int szalloc, int card, int nelems )
    :   DstScoreTable<TScore>(szalloc, card, nelems),
        h_scores_(scores)
    {};

    DstScoreTableSM();
    virtual ~DstScoreTableSM();

    virtual const TScore* GetScores() const { return h_scores_; }

    TScore GetScore(int row, int col) const
    {
        return DstScoreTable<TScore>::GetScore(row, col);
    }

    TScore GetDstScore(TScore d1, TScore d2);
    TScore GetDstScore(int d1, int d2);

    void PrintScores();

protected:
    virtual TScore AccessScore(int ndx) const { return h_scores_[ndx]; }

protected:
    TScore* h_scores_;//multi-dimensional scores
};

// -------------------------------------------------------------------------
// INLINES
//
// Default constructor
// 
template <typename TScore>
inline
DstScoreTableSM<TScore>::DstScoreTableSM()
:   DstScoreTable<TScore>(),
    h_scores_( NULL )
{
}

// -------------------------------------------------------------------------
// Destructor
// 
template <typename TScore>
inline
DstScoreTableSM<TScore>::~DstScoreTableSM()
{
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <typename TScore>
inline
TScore DstScoreTableSM<TScore>::GetDstScore(TScore d1, TScore d2)
{
    int row = (int)d1;//rintf(d1);
    int col = (int)d2;//rintf(d2);
    row = MYMIN(row, card_-1);
    col = MYMIN(col, card_-1);
    return DstScoreTableSM<TScore>::GetScore(row, col);
}

// -------------------------------------------------------------------------
// GetDstScore: get the score at the position identified by distance values 
// d1 and d2
// 
template <typename TScore>
inline
TScore DstScoreTableSM<TScore>::GetDstScore(int d1, int d2)
{
    if(d1 == 0xff || d2 == 0xff)
        return TScore(-1.f);
    d1 = MYMIN(d1, card_-1);
    d2 = MYMIN(d2, card_-1);
    return DstScoreTableSM<TScore>::GetScore(d1, d2);
}

// -------------------------------------------------------------------------
// PrintScores: print all scores for testing
// 
template <typename TScore>
inline
void DstScoreTableSM<TScore>::PrintScores()
{
    fprintf(stderr, "\n## DstScoreTableSM:\n");
    for(int n=0, i=0; i < card_; i++) {
        fprintf(stderr, "## %2d:", i);
        for(int j=0; j <= i; j++)
            fprintf(stderr, " %.3f", h_scores_[n++]);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}



// =========================================================================
// class SerializedScoresCtor
//
// Implementation of construction of serialized distance match scores; 
// serialization transforms a 2D score table to a 1D array
//
template <typename TScore, int APPROX>
class DstScoreTableCtor: public DstScoreTableSM<TScore>
{
    using DstScoreTableSM<TScore>::h_scores_;
    using DstScoreTableSM<TScore>::szalloc_;
    using DstScoreTableSM<TScore>::card_;
    using DstScoreTableSM<TScore>::nelems_;

public:
    DstScoreTableCtor(
        int maxdst, 
        float maxdstinv,
        float approx_score00,
        float approx_Score200,
        float parTheta, float parAbsDiffExp, float parAvgDistExp);

    virtual ~DstScoreTableCtor();

protected:
    void NewScores(int sztotal);
    float CalculateScore(int d1, int d2, int maxdst,
        float parTheta, float parAbsDiffExp, float parAvgDistExp);
    void DestroyScores();
};

// -------------------------------------------------------------------------
// INLINES
//
// Destructor
// 
template <typename TScore, int APPROX>
inline
DstScoreTableCtor<TScore, APPROX>::~DstScoreTableCtor()
{
    DestroyScores();
}

// -------------------------------------------------------------------------
// DestroyScores: destroy scores
// 
template <typename TScore, int APPROX>
inline
void DstScoreTableCtor<TScore, APPROX>::DestroyScores()
{
    if( h_scores_ ) {
        free(h_scores_);
        h_scores_ = NULL;
    }
    szalloc_ = 0;
}

// -------------------------------------------------------------------------
// NewScores: allocate memory for scores
// 
template <typename TScore, int APPROX>
inline
void DstScoreTableCtor<TScore, APPROX>::NewScores(int sztotal)
{
    DestroyScores();
    h_scores_ = (TScore*)malloc(sztotal);
    if( !h_scores_ ) {
        fprintf(stderr, "ERROR: DstScoreTableCtor::NewScores: Not enough memory.\n");
        abort();
    }
    szalloc_ = sztotal;
}

// -------------------------------------------------------------------------
// CalculateScore: calculate score given a pair of distance values d1 and d2
//
template <typename TScore, int APPROX>
inline
float DstScoreTableCtor<TScore, APPROX>::CalculateScore(int d1, int d2,
    int maxdst,
    float parTheta, float parAbsDiffExp, float parAvgDistExp)
{
    float davg = 0.5f * (float)(d1+d2);
    float dabs = (float)abs(d1-d2);
    float sco = parTheta;
    if(davg)
        sco = 
            ////(parTheta - dabs/davg) * expf(-davg*davg / 20.f));//DMS_Alpha));
            (parTheta - dabs/davg) * 
                //the below represents the score weight; the first factor 
                //implies scaling a (negative) score increasingly with 
                //increasing absolute difference; the second factor 
                //downscales the score as the average distance increases:
                (1.f + powf(dabs/(float)maxdst, parAbsDiffExp)) *
                (1.f - powf(davg/(float)maxdst, parAvgDistExp));
    return sco;
}

// -------------------------------------------------------------------------
// constructor: calculate and serialize 2D distance match scores
// maxdstinv, inverse of maxdst;
// approx_score00, score at (0,0) for two-level approximation;
// approx_Score200, score at (maxdst,0) for two-level approximation;
//
template <typename TScore, int APPROX>
inline
DstScoreTableCtor<TScore, APPROX>::DstScoreTableCtor(
    int maxdst, 
    float maxdstinv,
    float approx_score00,
    float approx_Score200,
    float parTheta, float parAbsDiffExp, float parAvgDistExp)
:   DstScoreTableSM<TScore>()
{
    if(maxdst <= 0 || 255 <= maxdst) {
        fprintf(stderr, 
        "ERROR: DstScoreTableCtor::DstScoreTableCtor: Invalid max distance value.\n");
        abort();
    }

    card_ = maxdst + 1;

    nelems_ = (card_ * (card_ + 1)) >> 1;//number of entries in the table

    //size calculated for scores:
    int szscores = nelems_ * sizeof(TScore);

    NewScores(szscores);

    std::vector<float> slopes, intpts;

    if(APPROX==2) {
        intpts.reserve(card_);
        for(int d1 = 0; d1 < card_; d1++) {
            float intpt = approx_score00 - (float)(d1) * approx_Score200 * maxdstinv;
            intpts.push_back(intpt);
        }
    }
    else if(APPROX==1) {
        slopes.reserve(card_);
        intpts.reserve(card_);

        for(int d1 = 0; d1 < card_; d1++)
        {
            int d21 = d1, d22 = card_-1;
            float sco1 = CalculateScore(d1, d21, maxdst, parTheta, parAbsDiffExp, parAvgDistExp);
            float sco2 = CalculateScore(d1, d22, maxdst, parTheta, parAbsDiffExp, parAvgDistExp);
            float slope = 0.0f;
            float intpt = 0.0f;
            if(d1 - d22) {
                //slope for the linear apporoximation of column d1: y1-y2/(x1-x2)
                //intercept: x1y2-x2y1/(x1-x2)
                slope = (sco1 - sco2) / (float)(d1-d22);
                intpt = ((float)(d1)*sco2 - (float)(d22)*sco1) / (float)(d1-d22);
            }
            slopes.push_back(slope);
            intpts.push_back(intpt);
        }
    }

    int n = 0;
    //slope for two-level approximation; does not depend on d1 and d2
    float slopeapp2 = (approx_Score200 - approx_score00) * maxdstinv;

    //write all score tables
    for(int d1 = 0; d1 < card_; d1++)
    {
        for(int d2 = 0; d2 <= d1; d2++) {
            float sco = 0.0f;
            if(APPROX==2)
                sco = (float)(d1) * slopeapp2 + intpts[d2];
            else if(APPROX==1)
                sco = (float)(d1) * slopes[d2] + intpts[d2];
            else
                sco = CalculateScore(d1, d2, maxdst, parTheta, parAbsDiffExp, parAvgDistExp);
            *(h_scores_ + n++) = (TScore)sco;
        }
    }
}

#endif//__dstable_h__
