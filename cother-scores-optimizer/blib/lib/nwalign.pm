package nwalign;

##
## 2007 (C) Mindaugas Margelevicius
## Institute of Biotechnology
## Vilnius, Lithuania
##

use strict;


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

my  $stAlign = 0; ## state aligned
my  $stGapUp = 1; ## state gap up
my  $stGapLf = 2; ## state gap left
my  $noStats = 3; ## number of states

my  $drNo        = 0; ## 000, no valid direction
my  $drDiag      = 1; ## 001, diagonal
my  $drUp        = 2; ## 010, up
my  $drDiag_Up   = 3; ## 011, diagonal or up
my  $drLeft      = 4; ## 100, left
my  $drDiag_Left = 5; ## 101, diagonal or left
my  $drUp_Left   = 6; ## 110, up or left
my  $drAll       = 7; ## 111, diagonal, up, or left

my  $defOpen   = -15;  ## default gap open cost
my  $defExtend =  0;  ## default gap extend cost
my  $SCMIN = -32767;

## -------------------------------------------------------------------

BEGIN {
    my $X =   0;    ## score to align X
    my $G =  -6;    ## score to align -
    my $GG = 0;     ## score to align -vs.-

    # Hash for conversion of amino acid to number
    my %CODESAA = ( 'A' =>  0, 'R' =>  1, 'N' =>  2, 'D' =>  3, 'C' =>  4, 'Q' =>  5, 
                    'E' =>  6, 'G' =>  7, 'H' =>  8, 'I' =>  9, 'L' => 10, 'K' => 11, 
                    'M' => 12, 'F' => 13, 'P' => 14, 'S' => 15, 'T' => 16, 'W' => 17, 
                    'Y' => 18, 'V' => 19, 'B' => 20, 'Z' => 21, 'X' => 22, '*' => 23,
                    '-' => 24
    );
    my @BLOSUM62 = (
        ## A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *   -
        [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, $X, -4, $G ], # A
        [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, $X, -4, $G ], # R
        [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, $X, -4, $G ], # N
        [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, $X, -4, $G ], # D
        [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, $X, -4, $G ], # C
        [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, $X, -4, $G ], # Q
        [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, $X, -4, $G ], # E
        [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, $X, -4, $G ], # G
        [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, $X, -4, $G ], # H
        [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, $X, -4, $G ], # I
        [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, $X, -4, $G ], # L
        [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, $X, -4, $G ], # K
        [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, $X, -4, $G ], # M
        [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, $X, -4, $G ], # F
        [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, $X, -4, $G ], # P
        [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0, $X, -4, $G ], # S
        [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, $X, -4, $G ], # T
        [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, $X, -4, $G ], # W
        [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, $X, -4, $G ], # Y
        [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, $X, -4, $G ], # V
        [ -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, $X, -4, $G ], # B
        [ -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, $X, -4, $G ], # Z
        [ $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, $X, -4, $G ], # X
        [ -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1, -4 ], # *
        [ $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, $G, -4,$GG ]  # -
    );

    $stAlign = 0; ## state aligned
    $stGapUp = 1; ## state gap up
    $stGapLf = 2; ## state gap left
    $noStats = 3; ## number of states

    $drNo        = 0; ## 000, no valid direction
    $drDiag      = 1; ## 001, diagonal
    $drUp        = 2; ## 010, up
    $drDiag_Up   = 3; ## 011, diagonal or up
    $drLeft      = 4; ## 100, left
    $drDiag_Left = 5; ## 101, diagonal or left
    $drUp_Left   = 6; ## 110, up or left
    $drAll       = 7; ## 111, diagonal, up, or left

    $defOpen   = -15;  ## default gap open cost
    $defExtend =   0;  ## default gap extend cost

    $SCMIN = -32767;

    sub GetBLOSUM62   { return \@BLOSUM62; }
    sub GetCODESAA    { return \%CODESAA; }
}

## ===================================================================

sub new {
    my $that  = shift;
    my $class = ref( $that ) || $that;
    my $self;

    $self->{QUERY}  = '';   ## query sequence
    $self->{SBJCT}  = '';   ## subject sequence
    $self->{OTHER}  = [];   ## ref. to other subjects to align simultaneously
    $self->{GAP} = '-';     ## gap symbol

    while( scalar( @_ )) {
        $self->{uc( $_[0] )} = $_[1];
        shift, shift;
    }

    return bless  $self, $class;
}

## -------------------------------------------------------------------
## read/write member methods
##

sub Query { my $self = shift; if (@_) { $self->{QUERY} = shift } return $self->{QUERY}; }
sub Sbjct { my $self = shift; if (@_) { $self->{SBJCT} = shift } return $self->{SBJCT}; }
sub Other { my $self = shift; if (@_) { $self->{OTHER} = shift } return $self->{OTHER}; }
sub Gapsm { my $self = shift; if (@_) { $self->{GAP} = shift } return $self->{GAP}; }




## ===================================================================
## prefer match states to gaps
##
sub GetState {
    my  $dirct = shift;
    my  $state = ( $dirct & 1 )? $stAlign: (( $dirct & 2 )? $stGapUp: (( $dirct & 4 )? $stGapLf: $noStats ));
    return $state;
}

## prefer gaps to match states
##
sub GetState2 {
    my  $dirct = shift;
    my  $state = ( $dirct & 4 )? $stGapLf: (( $dirct & 2 )? $stGapUp: (( $dirct & 1 )? $stAlign: $noStats ));
    return $state;
}

## -------------------------------------------------------------------
## align pairwise locally sequences
##

sub Align
{
    my  $self = shift;
    my  $class = ref( $self ) || die( "ERROR: Align: Should be called by object." );

    my  $algnquery = shift; ## reference to aligned query sequence
    my  $algnsbjct = shift; ## reference to aligned subject sequence
    my  $algnother = shift; ## reference to sequences to be aligned in turn

    my  $queryseq = \$self->{QUERY};  ## reference to query sequence
    my  $sbjctseq = \$self->{SBJCT};  ## reference to subject sequence
    my  $rothersq =  $self->{OTHER};  ## ref. to other subject sequences to align simultaneously

    my  @DPM;         ## dynamkic programming matrix
    my  @BTpnt;       ## back-tracing pointer
    my  @APath;       ## alignment path

    $self->AlignWAffineGap( $queryseq, $sbjctseq, \@DPM, \@BTpnt );
    $self->MakeAlignmentPath( \@DPM, \@BTpnt, \@APath );
    $self->DrawAlignment( $queryseq, $sbjctseq, $rothersq, \@APath, $algnquery, $algnsbjct, $algnother );
}

## -------------------------------------------------------------------
## align sequences by Needleman-Wunsch using affine gap cost scheme
##

sub AlignWAffineGap
{
    my  $self = shift;
    my  $class = ref( $self ) || die( "ERROR: AlignWAffineGap: Should be called by object." );

    my  $queryseq = shift;    ## reference to query sequence
    my  $sbjctseq = shift;    ## reference to subject sequence
    my  $DPM = shift;         ## reference to dynamkic programming matrix
    my  $BTpnt = shift;       ## reference to back-tracing pointer
    my  $open = shift;        ## gap open cost
    my  $extend = shift;      ## gap extend cost

    my  $matrix = GetBLOSUM62();

    $open = $defOpen unless $open;
    $extend = $defExtend unless $extend;

    ## convert sequences into arrays
    my  @queryaas = split( //, uc( $$queryseq ));
    my  @sbjctaas = split( //, uc( $$sbjctseq ));

    ## change symbolic notation into numbers
    $_ = ${GetCODESAA()}{$_} foreach( @queryaas );
    $_ = ${GetCODESAA()}{$_} foreach( @sbjctaas );

    my ($n, $m, $s );

    my ($bestA, $currentA, $A );     ## the best, current, and operating scores for state stAlign
    my ($bestU, $currentU, $U );     ## the best, current, and operating scores for state stGapUp
    my ($bestL, $currentL, $L );     ## the best, current, and operating scores for state stGapLf

    my  $upOpen;                 ## score after the gap opening cost is evaluated for a query position
    my  $leftOpen;               ## score after the gap opening cost is evaluated for a subject position
    my  $upExtend;               ## score after the gap extension cost is evaluated for a query position
    my  $leftExtend;             ## score after the gap extension cost is evaluated for a subject position

    my  $scoring;                ## score to align the query and subject at a position
    my  $ptr;                    ## value of the pointer for one cell of the dynamic programming matrix


    ## initialize matrices
    ##
    $#{$DPM} = $#sbjctaas + 1;
    $#{$BTpnt} = $#sbjctaas + 1;

    for( $m = 0; $m <= scalar( @sbjctaas ); $m++ ) {
        $#{$DPM->[$m]} = $#queryaas + 1;
        $#{$BTpnt->[$m]} = $#queryaas + 1;

        for( $n = 0; $n <= scalar( @queryaas ); $n++ ) {
            $#{$DPM->[$m][$n]} = $noStats - 1;
            $#{$BTpnt->[$m][$n]} = $noStats - 1;

            if( $n == 0 && $m ) {
                $DPM->[$m][$n][$stGapUp] = $open + $extend * ( $m - 1 );
                $DPM->[$m][$n][$stGapLf] = $SCMIN;
                $DPM->[$m][$n][$stAlign] = $SCMIN;
                $BTpnt->[$m][$n][$stGapUp] = $drUp;
                $BTpnt->[$m][$n][$stGapLf] = 0;
                $BTpnt->[$m][$n][$stAlign] = 0;
            }
            elsif( $m == 0 && $n ) {
                $DPM->[$m][$n][$stGapLf] = $open + $extend * ( $n - 1 );
                $DPM->[$m][$n][$stGapUp] = $SCMIN;
                $DPM->[$m][$n][$stAlign] = $SCMIN;
                $BTpnt->[$m][$n][$stGapLf] = $drLeft;
                $BTpnt->[$m][$n][$stGapUp] = 0;
                $BTpnt->[$m][$n][$stAlign] = 0;
            }
            else {
                for( $s = 0; $s < $noStats; $s++ ) {
                    $DPM->[$m][$n][$s] = 0;
                    $BTpnt->[$m][$n][$s] = 0;
                }
            }
        }
    }

    $DPM->[0][0][$stAlign] = 0;
    $DPM->[0][0][$stGapUp] = $SCMIN;
    $DPM->[0][0][$stGapLf] = $SCMIN;


    for( $n = 1; $n <= scalar( @queryaas ); $n++ )
    {
        $bestA     = $DPM->[0][$n][$stAlign];
        $currentA  = $DPM->[0][$n-1][$stAlign];
        $bestU     = $DPM->[0][$n][$stGapUp];
        $currentL  = $DPM->[0][$n-1][$stGapLf];

        for( $m = 1; $m <= scalar( @sbjctaas ); $m++ )
        {                                       ## score to align at positions m and n
            $scoring  = $matrix->[ $sbjctaas[$m-1] ][ $queryaas[$n-1] ];
            $upOpen   = $bestA + $open;         ## gap open cost for query position
            $A        = $currentA;
            $currentA = $DPM->[$m][$n-1][$stAlign];
            $leftOpen = $currentA + $open;      ## gap open cost for subject position


            $L        = $currentL;
            $currentL = $DPM->[$m][$n-1][$stGapLf];
            $leftExtend = $currentL + $extend;  ## gap extend cost for subject position

            $upExtend = $bestU + $extend;       ## gap extend cost for query position
            $U = $DPM->[$m-1][$n-1][$stGapUp];

            ## process state stGapUp
            ##
            if( $upExtend < $upOpen ) {
                $bestU = $upOpen;
                $ptr = $drDiag;                 ## diagonal direction for backtracing
            } elsif( $upOpen < $upExtend ) {
                        $bestU = $upExtend;
                        $ptr = $drUp;           ## direction up
                    } else {                    ## $upOpen == $upExtend
                        $bestU = $upOpen;
                        $ptr = $drDiag_Up;      ## diagonal or up
                    }
            $DPM->[$m][$n][$stGapUp] = $bestU;
            $BTpnt->[$m][$n][$stGapUp] = $ptr;

            ## process state stGapLf
            ##
            if( $leftExtend < $leftOpen ) {
                $bestL = $leftOpen;
                $ptr = $drDiag;                 ## diagonal direction for backtracing
            } elsif( $leftOpen < $leftExtend ) {
                        $bestL = $leftExtend;
                        $ptr = $drLeft;         ## direction left
                    } else {                    ## $leftOpen == $leftExtend
                        $bestL = $leftOpen;
                        $ptr = $drDiag_Left;    ## diagonal or left
                    }
            $DPM->[$m][$n][$stGapLf] = $bestL;
            $BTpnt->[$m][$n][$stGapLf] = $ptr;

            ## process state stAlign
            ## check previous scores to correctly set ptr value
            ##
            if( $U < $A ) {
                if( $L < $A ) {
                    $bestA = $A;
                    $ptr = $drDiag;
                } elsif( $A < $L ) {
                            $bestA = $L;
                            $ptr = $drLeft;
                        } else {                ## $A == $L
                            $bestA = $A;
                            $ptr = $drDiag_Left;
                        }
            } elsif( $A < $U ) {
                        if( $L < $U ) {
                            $bestA = $U;
                            $ptr = $drUp;
                        } elsif( $U < $L ){
                                    $bestA = $L;
                                    $ptr = $drLeft;
                                } else {        ## $U == $L
                                    $bestA = $U;
                                    $ptr = $drUp_Left;
                                }
                    } else {                    ## $A == $U
                        if( $L < $A ){
                            $bestA = $A;
                            $ptr = $drDiag_Up;
                        } elsif( $A < $L ) {
                                    $bestA = $L;
                                    $ptr = $drLeft;
                                } else {        ## $A == $L
                                    $bestA = $A;
                                    $ptr = $drAll;
                                }
                    }
            $bestA += $scoring;
            $DPM->[$m][$n][$stAlign] = $bestA;
            $BTpnt->[$m][$n][$stAlign] = $ptr;
        }
    }
}

## -------------------------------------------------------------------
## make alignment path according to dynamic programming matrix
## computed
##

sub MakeAlignmentPath
{
    my  $self = shift;
    my  $class = ref( $self ) || die( "ERROR: MakeAlignmentPath: Should be called by object." );

    my  $DPM = shift;         ## reference to dynamkic programming matrix
    my  $BTpnt = shift;       ## reference to back-tracing pointer
    my  $APath = shift;       ## reference to alignment path

    my  $score = 0;           ## maximum score from the dynamic programming matrix
    my  $row = 0;             ## row index of the maximum score
    my  $column = 0;          ## column index of the maximum score
    my ($row2, $column2 );    ## half row and column
    my  $laststate = $stAlign;## last state back-tracing started with
    my  $state = $laststate;  ## state of back-tracing
    my  $step = 0;            ## alignment step (index)

    my ($m, $n );

    ## find the most distant maximum value
    ## (for local alignments)
    ##
##    for( $m = $#{$DPM}; 0 < $m; $m-- ) {
##        for( $n = $#{$DPM->[$m]}; 0 < $n; $n-- ) {
##            if( $score < $DPM->[$m][$n][$stAlign] ) {    ## check only stAlign, i.e. aligned states
##                $score = $DPM->[$m][$n][$stAlign];
##                $row = $m;
##                $column = $n;
##            }
##    }   }

##    return if $score <= 0;

    $row = $#{$DPM};
    $column = $#{$DPM->[$row]};
    ## prefer gaps to match states in the end of global alignment
    ##
    do { $score = $DPM->[$row][$column][$stGapLf]; $laststate = $stGapLf };
    do { $score = $DPM->[$row][$column][$stGapUp]; $laststate = $stGapUp } if $score < $DPM->[$row][$column][$stGapUp];
    do { $score = $DPM->[$row][$column][$stAlign]; $laststate = $stAlign } if $score < $DPM->[$row][$column][$stAlign];

    $row2 = int( $row / 2 ) + ( $row % 2 != 0 );
    $column2 = int( $column / 2 ) + ( $column % 2 != 0 );

    $APath->[$step][0] = $column;     ## index for query
    $APath->[$step][1] = $row;        ## index for subject
    $step++;

    while( 0 < $row || 0 < $column ) {
        $state = $laststate;
        if( $row2 < $row || $column2 < $column ) { 
                 $laststate = GetState( $BTpnt->[$row][$column][$state] ); ## GetState2
        } else { $laststate = GetState(  $BTpnt->[$row][$column][$state] ); }

        SWITCH: {
            if( $state == $stAlign ) { $row--; $column--; last SWITCH }
            if( $state == $stGapUp ) { $row--;            last SWITCH }
            if( $state == $stGapLf ) { $column--;         last SWITCH }
            die "ERROR: Unkown alignment state.";
        }

        $APath->[$step][0] = $column;
        $APath->[$step][1] = $row;
        $step++;
    }
}

## -------------------------------------------------------------------
## make query and subject alignment sequences according to alignment
## path generated
##

sub DrawAlignment
{
    my  $self = shift;
    my  $class = ref( $self ) || die( "ERROR: DrawAlignment: Should be called by object." );

    my  $queryseq = shift;    ## reference to query sequence
    my  $sbjctseq = shift;    ## reference to subject sequence
    my  $rothersq = shift;    ## reference to other sequences to align them simult.
    my  $APath = shift;       ## reference to alignment path
    my  $algnquery = shift;   ## reference to aligned query sequence
    my  $algnsbjct = shift;   ## reference to aligned subject sequence
    my  $algnother = shift;   ## reference to sequences aligned in turn with subject

    my  $step;
    my  $gap = $self->{GAP};
    my ($n, $ind );

    $$algnquery = '';
    $$algnsbjct = '';

    if( ref( $rothersq ) eq 'ARRAY' && ref( $algnother ) eq 'ARRAY' ) {
        for( $n = 0; $n <= $#{$rothersq}; $n++ ) {
            $$algnother[$n] = '';
        }
    }

    ## -1: the very beginning points to cell of score 0
    ##
    for( $step = $#{$APath} - 1; 0 <= $step; $step-- )
    {
        ## if previous index matches current index
        if( $step < $#{$APath}  &&  $APath->[$step][0] == $APath->[$step+1][0] ) {
            $$algnquery .= $gap;
        } else {
            ## -1 to compensate one extra index in dynamic prog. matrix...
            $$algnquery .= substr( $$queryseq, $APath->[$step][0] - 1, 1 );
        }

        if( $step < $#{$APath}  &&  $APath->[$step][1] == $APath->[$step+1][1] ) {
            $$algnsbjct .= $gap;
            if( ref( $rothersq ) eq 'ARRAY' && ref( $algnother ) eq 'ARRAY' ) {
                for( $n = 0; $n <= $#{$rothersq}; $n++ ) {
                    $$algnother[$n] .= $gap;
                }
            }
        } else {
            $$algnsbjct .= substr( $$sbjctseq, $APath->[$step][1] - 1, 1 );
            if( ref( $rothersq ) eq 'ARRAY' && ref( $algnother ) eq 'ARRAY' ) {
                $ind = $APath->[$step][1] - 1;
                for( $n = 0; $n <= $#{$rothersq}; $n++ ) {
                    $$algnother[$n] .= substr( $$rothersq[$n], $ind, 1 ) if $ind < length( $$rothersq[$n] );
                }
            }
        }
    }
}

## -------------------------------------------------------------------

1;

