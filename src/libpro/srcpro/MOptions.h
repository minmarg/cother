/***************************************************************************
 *   Copyright (C) 2013-2021 Mindaugas Margelevicius                       *
 *   Institute of Biotechnology, Vilnius University                        *
 ***************************************************************************/

#ifndef __MOptions__
#define __MOptions__

#include "liblib/preproc.h"

// =========================================================================
// CONSTANTS
//
// #define MAX_VALUE_FOR_DSTNMAX 1000
#define MAX_N_IR_DISTANCE_VALUES 99
#define MAX_N_IR_DISTANCE_VALUES_PLUS_LENGTH (MAX_N_IR_DISTANCE_VALUES+1)
//
// =========================================================================

class mystring;

// the OPER argument in the following macro represents operator along with 
// operand, e.g., `*0.01'.
#define DECLAREOPTION( NAME, TYPE, TYPECAST, OPER ) \
  extern TYPE val##NAME##_; \
  static inline TYPECAST Get##NAME() { return (TYPECAST)val##NAME##_ OPER; } \
  ; // void Read##NAME();

namespace MOptions {

extern const char* SECOPTIONS;
void Read( const char* filename );

DECLAREOPTION( EVAL, float, float, );
DECLAREOPTION( NOHITS, int, int, );
DECLAREOPTION( NOALNS, int, int, );
DECLAREOPTION( SHOW, int, int, );
DECLAREOPTION( DSCLEN, int, int, );
DECLAREOPTION( DSCWIDTH, int, int, );
DECLAREOPTION( ALNWIDTH, int, int, );

DECLAREOPTION( X_UNINF, int, int, );
DECLAREOPTION( IDENTITY, int, float, *0.01f );
DECLAREOPTION( SHOWCMD, int, int, );

DECLAREOPTION( DELSTATE, int, int, );
DECLAREOPTION( PCFWEIGHT, float, float, );
DECLAREOPTION( MINALNFRN, int, float, *0.01f );
DECLAREOPTION( MINALNPOS, int, int, );

DECLAREOPTION( SUBMAT, mystring, mystring, );


// DECLAREOPTION( DSTBNDF, int, int, );
// DECLAREOPTION( DSTNMAX, int, int, );
DECLAREOPTION( DSTSEGM, int, int, );
DECLAREOPTION( DSTNPOS, int, int, );
DECLAREOPTION( DSTSEP6, float, float, );
DECLAREOPTION( DSTGAPL, float, float, );
DECLAREOPTION( DSTFREE, int, int, );


DECLAREOPTION( TFRMIX, mystring, mystring, );
DECLAREOPTION( MIXWGT, float, float, );
DECLAREOPTION( SCOADJ, mystring, mystring, );
DECLAREOPTION( SUPCLT, int, int, );
DECLAREOPTION( ADJWGT, float, float, );
DECLAREOPTION( cADJWGT, float, float, );

DECLAREOPTION( DDMSWGT, float, float, );

DECLAREOPTION( CVSWGT, float, float, );

DECLAREOPTION( SSSWGT, float, float, );

DECLAREOPTION( SSEMODEL, int, int, );

DECLAREOPTION( HCFILTER, int, int, );
DECLAREOPTION( HCWINDOW, int, int, );
DECLAREOPTION( HCLOWENT, float, float, );
DECLAREOPTION( HCHIGHENT, float, float, );

DECLAREOPTION( INVLCFILTER, int, int, );
DECLAREOPTION( LCFILTEREACH, int, int, );
DECLAREOPTION( LCWINDOW, int, int, );
DECLAREOPTION( LCLOWENT, float, float, );
DECLAREOPTION( LCHIGHENT, float, float, );
DECLAREOPTION( DISTANCE, float, float, );

DECLAREOPTION( SCHEME, mystring, mystring, );
DECLAREOPTION( MAPALN, int, int, );
DECLAREOPTION( MINPP, float, float, );

DECLAREOPTION( HSPLEN, int, int, );
DECLAREOPTION( HSPMINSCORE, int, int, );
DECLAREOPTION( HSPMAXDIST, int, int, );
DECLAREOPTION( NOHSPS, int, int, );

}//namespace MOptions

#endif//__MOptions__
