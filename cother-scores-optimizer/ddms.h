/***************************************************************************
 *   Copyright (C) 2020-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __ddms_h__
#define __ddms_h__

#include <xmmintrin.h>
#include "dstable.h"

// -------------------------------------------------------------------------
// unrolling degree, which determines extra space required for arrays
#define CSO_DDMS_UNROLLING_DEGREE 4

// cother option DSTSEGM
#define COTHER_OPTION_DSTSEGM 6


//flag of whether the scores are approximated by linear interpolation
//1, one-level approximation;
//2, two-level approximation given two limit values for highest and lowest 
// starting scores along the diagonal and bottom line of the score matrix
#define CSO_DDMS_DScore_Approx 2
#define CSO_DDMS_DScore_Approx_Score00 0.35f //score at (0,0)
#define CSO_DDMS_DScore_Approx_Score200 -0.4f ////score at (20,0)

// Parameters for the score between two distance values:
//max distance value considered; cardinality is 1 greater (includes 0)
#define CSO_DDMS_DScore_MaxDst 20
//inverse of CSO_DDMS_DScore_MaxDst
#define CSO_DDMS_DScore_MaxDst_inv 0.05f
//default theta parameter
#define CSO_DDMS_DScore_Theta 0.2f
//default alpha parameter
// #define CSO_DDMS_DScore_Alpha 20.f
//default parameter for the exponent of the absolute difference
#define CSO_DDMS_DScore_AbsDiffExp 0.2f
//default parameter for the exponent of the average distance
#define CSO_DDMS_DScore_AvgDistExp 1.8f
//default parameter for the exponent of the average distance
#define CSO_DDMS_DScore_Granularity 16.f
//number of ungapped maximum-scoring segments to consider in dynamic 
// programming when calculating the distance distribution match score 
// between two positions
#define CSO_DDMS_DP_NSEGM 3

//using heuristics of linear inpterpolation for half of DP cells
#define CSO_DDMS_DP_HEURISTIC 1
//sparsity of processing the DP matrix (step over diagonals)
#define CSO_DDMS_DP_STEP 2
//fraction of the size of distances (number of distance values) the 
//(processing) band occupies, expressed as an exponent of the 
//denominator to 1/2^x
#define CSO_DDMS_DP_BAND_FRACTION 2

// code to ignore score
#define CSO_DDMS_IGNORE_CODE (-777.f)
// code to specify an error 
#define CSO_DDMS_ERR_CODE (-999.f)

// -------------------------------------------------------------------------
// class CSO_DDMS for calculating distance distribution match scores
//
class CSO_DDMS {
public:
    CSO_DDMS(
        float theta = CSO_DDMS_DScore_Theta,
        float absdiffexp = CSO_DDMS_DScore_AbsDiffExp,
        float avgdistexp = CSO_DDMS_DScore_AvgDistExp);
    ~CSO_DDMS();
    int GetUnrollingDegree() const {return CSO_DDMS_UNROLLING_DEGREE;}
    int GetOptionDSTSEGM() const {return COTHER_OPTION_DSTSEGM;}


    int GetDScoreApprox() const {return CSO_DDMS_DScore_Approx;}
    float GetDScoreApproxScore00() const {return CSO_DDMS_DScore_Approx_Score00;}
    float GetDScoreApproxScore200() const {return CSO_DDMS_DScore_Approx_Score200;}
    int GetDScoreMaxDst() const {return CSO_DDMS_DScore_MaxDst;}
    float GetDScoreTheta() const {return cso_ddms_dscore_theta_;}
    float GetDScoreAbsDiffExp() const {return cso_ddms_dscore_absdiffexp_;}
    float GetDScoreAvgDistExp() const {return cso_ddms_dscore_avgdistexp_;}
    float GetDScoreGranularity() const {return CSO_DDMS_DScore_Granularity;}
    int GetDPNSegments() const {return CSO_DDMS_DP_NSEGM;}

    int GetDPHeuristic() const {return CSO_DDMS_DP_HEURISTIC;}
    int GetDPStep() const {return CSO_DDMS_DP_STEP;}
    int GetDPBandFraction() const {return CSO_DDMS_DP_BAND_FRACTION;}


    float GetIgnoreCode() const {return CSO_DDMS_IGNORE_CODE;}
    float GetErrorCode() const {return CSO_DDMS_ERR_CODE;}

    void PrintTable() {dstable.PrintScores();}
    float CalculateDDMS(int*, int*);

protected:
    //distance matching score for a pair of distance values d1 and d2
    float GetPosDstMatchScoreh2(int d1, int d2) {
        return dstable.GetDstScore(d1, d2);
    }

    template<int scndx>
    void UpdateSegScoreUpperh2(
        float scij,
        int i, int j,
        int& segqb, int& segpb,
        int& im, int& jim, int& ib, int& jib,
        float& segsco, float& scom);

    template<int scndx>
    void UpdateSegScoreLowerh2(
        float scij,
        int i, int j,
        int& segqb, int& segpb,
        int& ijm, int& jm, int& ijb, int& jb,
        float& segsco, float& scom);

    template<int N>
    void ReadDstMatchScoresh2Unrolled(
        const int* vdst1,
        const int * vdst2,
        float* scij, int i, int j);

    template<int N>
    void UpdateSegScoreUpperh2Unrolled(
        float* scij,
        int i, int j,
        int& segqb, int& segpb,
        int& im, int& jim, int& ib, int& jib,
        float& segsco, float& scom);

    template<int N>
    void UpdateSegScoreLowerh2Unrolled(
        float* scij,
        int i, int j,
        int& segqb, int& segpb,
        int& ijm, int& jm, int& ijb, int& jb,
        float& segsco, float& scom);

    void ProcessDiagonalUpperh2(
        const int* vdst1, const int* vdst2,
        const int len1, const int len2,
        const int j,
        float* maxscores,
        __m64* coords_maxscores);

    void ProcessDiagonalLowerh2(
        const int* vdst1, const int* vdst2,
        const int len1, const int len2,
        const int i,
        float* maxscores,
        __m64* coords_maxscores);

private:
    float cso_ddms_dscore_theta_;
    float cso_ddms_dscore_absdiffexp_;
    float cso_ddms_dscore_avgdistexp_;
    DstScoreTableCtor<float, CSO_DDMS_DScore_Approx> dstable;
};

// INLINES -----------------------------------------------------------------
//
template <typename T> 
inline
T mymax( T a, T b ){ return a<b?b:a; }


#endif//__ddms_h__
