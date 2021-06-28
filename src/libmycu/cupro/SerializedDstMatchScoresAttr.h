/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __SerializedDstMatchScoresAttr_h__
#define __SerializedDstMatchScoresAttr_h__

#include <cuda_fp16.h>

// #include "liblib/msg.h"
// #include "liblib/mybase.h"
#include "libpro/srcpro/MOptions.h"

//the scores are approximated by linear interpolation
//NOTE: 2, two-level linear interpolation given two limit values for 
// highest and lowest starting scores along the diagonal and bottom 
// line of the score matrix;
#define DMS_AS_2 2
//NOTE: 1, linear interpolation resolved through template 
// specialization with __half2;
#define DMS_AS_1 1
//NOTE: 0, no approximation
#define DMS_AS_0 0
#define DMS_APPROXIMATE_SCORES DMS_AS_2
// #define DMS_Approx2_Scores_028_m08
// #define DMS_Approx2_Score00 0.28f //score at (0,0) for DMS_APPROXIMATE_SCORES==2
// #define DMS_Approx2_Score200 -0.8f //score at (20,0) for DMS_APPROXIMATE_SCORES==2
// #define DMS_Approx2_Scores_035_m08
// #define DMS_Approx2_Score00 0.35f //score at (0,0) for DMS_APPROXIMATE_SCORES==2
// #define DMS_Approx2_Score200 -0.8f //score at (20,0) for DMS_APPROXIMATE_SCORES==2
#define DMS_Approx2_Scores_035_m04
#define DMS_Approx2_Score00 0.35f //score at (0,0) for DMS_APPROXIMATE_SCORES==2
#define DMS_Approx2_Score200 -0.4f //score at (20,0) for DMS_APPROXIMATE_SCORES==2

//linear interpolation for half of the distance score map 
//between two profile positions
#define DMS_LINEAR_INTERP_HALF true
//step for traversing the distance score map between two profile 
//positions; defines processing scarcity
#define DMS_STEP 2
// Fraction of the size of distances (number of
// distance values) the (processing) band occupies 
// in the DP matrix [0-5]; it is expressed as the 
// denominator's exponent x to 1/2^x; the larger 
// this value, the narrower band (0 == full matrix)
// (programs: cother)
#define DMS_DSTBNDF 2

//{{DDMS2S scores
//scale for DDMS scores before regression
#define CUSCO_DDMS2S_TRANS_SCALE 16.f
//}}

//attributes of serialized distance match scores
struct SerializedDstMatchScoresAttr {
    SerializedDstMatchScoresAttr( 
            int ntotents, int szalloc,
            int card, int nelems)
    :   ntotents_(ntotents), szalloc_(szalloc),
        card_(card), nelems_(nelems)
    {}
    SerializedDstMatchScoresAttr( const SerializedDstMatchScoresAttr& attr )
    :   ntotents_(attr.ntotents_), szalloc_(attr.szalloc_),
        card_(attr.card_), nelems_(attr.nelems_)
    {}
    SerializedDstMatchScoresAttr()
    :   ntotents_(0), szalloc_(0),
        card_(0), nelems_(0)
    {}
    int ntotents_;//total number of entries in the buffer, included cardinalities, etc.
    int szalloc_;//size (bytes) allocated for scores
    int card_;//cardinality determinator, # rows of a square score table
    int nelems_;//number of entries in one table
};

#if DMS_APPROXIMATE_SCORES == DMS_AS_2
#   define DMS_TYPE float
#   if MAX_N_IR_DISTANCE_VALUES == 99
#       if DMS_LINEAR_INTERP_HALF == true
#           if DMS_STEP == 2
#               if DMS_DSTBNDF == 2
#                   ifdef DMS_Approx2_Scores_028_m08
#                       define CUSCO_DDMS2S_TRANS_SLOPE 0.082f
#                       define CUSCO_DDMS2S_TRANS_INTERCEPT -0.889f
#                   elif defined(DMS_Approx2_Scores_035_m08)
#                       define CUSCO_DDMS2S_TRANS_SLOPE 0.062f
#                       define CUSCO_DDMS2S_TRANS_INTERCEPT -0.996f
#                   elif defined(DMS_Approx2_Scores_035_m04)
#                       define CUSCO_DDMS2S_TRANS_SLOPE 0.053f
#                       define CUSCO_DDMS2S_TRANS_INTERCEPT -1.235f
#                   else
#                       error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#                   endif
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#           endif
#       else
#           error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#       endif
#   else
#       error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#   endif
#elif DMS_APPROXIMATE_SCORES == DMS_AS_1
#   define DMS_TYPE __half2
#   if MAX_N_IR_DISTANCE_VALUES == 99
#       if DMS_LINEAR_INTERP_HALF == true
#           if DMS_STEP == 2
#               if DMS_DSTBNDF == 2
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.075f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.413f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#           endif
#       else
#           error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#       endif
#   else
#       error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES!=99)."
#   endif
#elif DMS_APPROXIMATE_SCORES == DMS_AS_0
#   define DMS_TYPE float
#   if MAX_N_IR_DISTANCE_VALUES == 99
#       if DMS_LINEAR_INTERP_HALF == true
#           if DMS_STEP == 2
#               if DMS_DSTBNDF == 3
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.068f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.192f
#               elif DMS_DSTBNDF == 2
                    //NOTE: changing slope and intercept so that the line passes 
                    //      through the point x0 at y=0 is equivalent to increasing or 
                    //      decreasing application weight DDMSWGT
                    //slope of a linear fit to the DDMS2S scores
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.069f
                    //intercept of a linear fit to the DDMS2S scores
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.39f
#               elif DMS_DSTBNDF == 1
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.069f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.548f
#               elif DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.069f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.594f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==99)."
#               endif
#           elif DMS_STEP == 1
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.066f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.689f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==99)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==99)."
#           endif
#       else//DMS_LINEAR_INTERP_HALF == false
#           if DMS_STEP == 1
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.084f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.796f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==99)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==99)."
#           endif
#       endif
#   elif MAX_N_IR_DISTANCE_VALUES == 63
#       if DMS_LINEAR_INTERP_HALF == true
#           if DMS_STEP == 2
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.073f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.31f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==63)."
#               endif
#           elif DMS_STEP == 1
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.072f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.438f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==63)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==63)."
#           endif
#       else//DMS_LINEAR_INTERP_HALF == false
#           if DMS_STEP == 1
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.09f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -1.532f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==63)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==63)."
#           endif
#       endif
#   elif MAX_N_IR_DISTANCE_VALUES == 31
#       if DMS_LINEAR_INTERP_HALF == true
#           if DMS_STEP == 2
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.053f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -0.51f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==31)."
#               endif
#           elif DMS_STEP == 1
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.057f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -0.621f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==31)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==31)."
#           endif
#       else//DMS_LINEAR_INTERP_HALF == false
#           if DMS_STEP == 1
#               if DMS_DSTBNDF == 0
#                   define CUSCO_DDMS2S_TRANS_SLOPE 0.065f
#                   define CUSCO_DDMS2S_TRANS_INTERCEPT -0.63f
#               else
#                   error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==31)."
#               endif
#           else
#               error "No DDMS parameters for given configuration (MAX_N_IR_DISTANCE_VALUES==31)."
#           endif
#       endif
#   endif
#else//DMS_APPROXIMATE_SCORES
#   error "No DDMS parameters for given configuration (DMS_APPROXIMATE_SCORES!={0,1,2})."
#endif//DMS_APPROXIMATE_SCORES

#endif//__SerializedDstMatchScoresAttr_h__
