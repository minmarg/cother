#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2021 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Script for making a map between distance distribution match scores and scores 
## translated according to their observed frequency along structural alignment.

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use Getopt::Long;
use List::Util qw(min max);
use POSIX qw(strftime);

use ddms;

use threads;
use threads::shared;

my $USINGTHREADS = eval 'use threads; 1';
my $NUMTHREADS = 1;

$|++;

my  $CMDL = join(' ', $0, @ARGV);
my  $MYPROGNAME = basename( $0 );
my  $EFFNOS = "1";## THRESHOLDS of EFFECTIVE No. SEQUENCES
my  $ZSCORE = 3.0;
my ($THETA, $ABSDIFFEXP, $AVGDISTEXP) = (0.2, 0.2, 1.8);

my  $usage = <<EOIN;

Program for making a map between distance distribution match scores and scores 
translated according to their observed frequency along structural alignment.
It also produces score distributions for structurally equivalent and non-
equivalent residues.
2021(C)Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$MYPROGNAME <Parameters>

Parameters:

--out <filename>   Output filename.
--aln <directory>  Directory of structural alignments.
--pro <directory>  Directory of COTHER profiles with predicted distances.
--prs <directory>  Directory of COTHER profiles with distances from structure.
--eff <numbers>    Make score maps regarding these levels of effective 
                   number of observations.
           Default=$EFFNOS
--zsc <score>      Consider alignments of at least this Z-score.
           Default=$ZSCORE
--theta <parameter> Theta parameter for scores.
           Default=$THETA
--adexp <parameter> Exponent for absolute difference (scores parameter).
           Default=$ABSDIFFEXP
--avexp <parameter> Exponent for average distance (scores parameter).
           Default=$AVGDISTEXP
--pws              Derive scores from pairwise alignments.
--ovl              Produce overlapping score tables.
--pos              Calculate positional sequence weights.
--nosw             Do not calculate sequence weights.
--ts  <number>     Number of concurrent threads to use if in use
                   (Threads can be used: $USINGTHREADS).
           Default=$NUMTHREADS
--printandexit     Print the table of scores between (pairwise) distance 
                   values and exit.
--help             This text.

EOIN

my  @EFF;##levels of efective number of sequences
my  $PRINTANDEXIT;
my  $OUTFILE;
my  $ALNDIR;
my  $PRODIR;
my  $PRSDIR;
my  $PWS = 0;
my  $OVL = 0;
my  $POS = 0;
my  $NOSW = 0;

my  $result = GetOptions(
               'out=s'     => \$OUTFILE,
               'aln=s'     => \$ALNDIR,
               'pro=s'     => \$PRODIR,
               'prs=s'     => \$PRSDIR,
               'eff=s'     => \$EFFNOS,
               'zsc=f'     => \$ZSCORE,
               'theta=f'     => \$THETA,
               'adexp=f'     => \$ABSDIFFEXP,
               'avexp=f'     => \$AVGDISTEXP,
               'pws'       => sub { $PWS = 1; },
               'ovl'       => sub { $OVL = 1; },
               'pos'       => sub { $POS = 1; },
               'nosw'      => sub { $NOSW = 1; },
               'ts=i'      => \$NUMTHREADS,
               'printandexit' => sub { $PRINTANDEXIT = 1; },
               'help|h'    => sub { print $usage; exit( 0 ); }
);


do { print $usage; exit( 1 ); }  unless $result;

die "ERROR: No output filename given." unless($PRINTANDEXIT || $OUTFILE);
die "ERROR: Alingment directory does not exist." unless($PRINTANDEXIT || ($ALNDIR && -d $ALNDIR));
die "ERROR: Profile directory does not exist." unless($PRINTANDEXIT || ($PRODIR && -d $PRODIR));
die "ERROR: Profile directory not found." unless($PRINTANDEXIT || ($PRSDIR && -d $PRSDIR));
die "ERROR: Invalid Z-score." if( $ZSCORE <= 0.0 || 1000.0 < $ZSCORE );

die "ERROR: Threads cannot be used." if( !$USINGTHREADS && $NUMTHREADS != 1 );
die "ERROR: Invalid number of threads." if( $USINGTHREADS && ( $NUMTHREADS <= 0 || 100 < $NUMTHREADS ));
$USINGTHREADS = $USINGTHREADS && 1<$NUMTHREADS;

@EFF = sort { $a <=> $b } split(',',$EFFNOS);
die "ERROR: No effective number of observations given." if( $#EFF < 0 );
for( my $i = 0; $i <= $#EFF; ){
    die "ERROR: Invalid effective number of observations." if( $EFF[$i] <= 0.0 || 20.0 < $EFF[$i] );
    ##remove duplicates
    do { splice( @EFF, $i, 1 ); next; } if( $i && $EFF[$i] == $EFF[$i-1]);
    $i++;
}

## ===================================================================
##
my  @alnfiles;
my (@profiles, @prsfiles);
my (@mythrds, @targs, $tret, $tref, $ti, $tii );
my  %THRES;

my ($f, $d, $ii, $t );
my  %PROS;##vector data at each position of each profile
my  %PRSS;##vector data for profiles with structural information
my  %ALNS :shared;##alignment query-subject names
my  %PAIRS;##profile pairs containing log odds scores for aligned positions
my  @FRQS;##frequencies
my (@BFS );##background frequencies
my  %SCOS;##derived score table
my ($noalns, @nopss) = (0,());
my  $NASC = -999999;
my  $INF = 999999;
my  $MINSC;
my  $MAXSC;
my ($ENTROPY, $EXPECTD, $IDV);
my @EINDDESC;

my $DDMS = new ddms::CSO_DDMS($THETA, $ABSDIFFEXP, $AVGDISTEXP);

## -------------------------------------------------------------------

sub ReadFiles;
sub GetDstDistrs;
sub UpdateFreqs;
sub ProcessTopHits;

sub MakeMultipleAlignment;
sub VerifyMultipleAlignment;
sub CalculateSequenceWeights;
sub ApplyWeights;
sub ApplyWeightsPw;##apply weights pairwise

sub DeriveScores;
sub PrintScores;

## -----------------------------------------------------------------------------
## GetTime: get time string
##
sub GetTime
{
    return strftime("%H:%M:%S ",localtime());
}

sub GetDatetime
{
    return strftime("%a %b %d %H:%M:%S %Z %Y ",localtime());
}

## -------------------------------------------------------------------
##

print( STDERR "$CMDL\n\n");
$DDMS->PrintTable();

exit(0) if $PRINTANDEXIT;

print( STDERR GetDatetime()."Reading directories...\n");

ReadFiles( $PRODIR, \@profiles );
ReadFiles( $PRSDIR, \@prsfiles );
ReadFiles( $ALNDIR, \@alnfiles );

printf( STDERR GetDatetime()."Reading profiles...\n");

foreach $f(@profiles) { 
    exit(1) unless GetDstDistrs("$PRODIR/$f", \%PROS); 
}

printf( STDERR GetDatetime()."Reading profiles w/ structural information...\n");

foreach $f(@prsfiles) { 
    exit(1) unless GetDstDistrs("$PRSDIR/$f", \%PRSS); 
}

##{{THREADS SECTION
my $nperthread = $#alnfiles+1;
if( $USINGTHREADS ) {
    $nperthread = int($nperthread/$NUMTHREADS);
    $nperthread++ if( ($#alnfiles+1) % $NUMTHREADS );
    $NUMTHREADS = $#alnfiles+1 if $nperthread < 2;
}
if( $USINGTHREADS ) {
    print( STDERR GetDatetime()."Processing alignments... ($nperthread/thread)\n");
    ##launch threads
    for( $t = 0; $t < $NUMTHREADS; $t++ ) {
        @targs = ( $t, $NUMTHREADS, \@alnfiles, \$noalns, \@nopss );
        $mythrds[$t] = threads->create({'context'=>'list'}, \&ThreadUpdateFreqs, @targs );
    }
    ##print( STDERR "    $NUMTHREADS threads launched.\n");
    ##join threads
    for( $t = 0; $t < $NUMTHREADS; $t++ ) {
        (%THRES) = $mythrds[$t]->join();
        print( STDERR GetDatetime()."  thread $t joined.\n");
        $tret = $THRES{RET};
        last unless $tret;
        ##infuse thread data
        $noalns += ${$THRES{ALNS}};
        $tref = $THRES{POSS};
        for( $ti = 0; $ti <= $#$tref; $ti++ ) {
            $nopss[$ti] += $$tref[$ti];
        }
        $tref = $THRES{FRQS};
        for( $ti = 0; $ti <= $#$tref; $ti++ ) {
            for( $tii = 0; $tii <= $#{$$tref[$ti]}; $tii++ ) {
                $FRQS[$ti][$tii]{$_} += $$tref[$ti][$tii]{$_} foreach( keys %{$$tref[$ti][$tii]});
            }
        }
        $tref = $THRES{BFS};
        for( $ti = 0; $ti <= $#$tref; $ti++ ) {
            for( $tii = 0; $tii <= $#{$$tref[$ti]}; $tii++ ) {
                $BFS[$ti][$tii] += $$tref[$ti][$tii];
            }
        }
        $tref = $THRES{EDSC};
        for( $ti = 0; $ti <= $#$tref; $ti++ ) {
            for( $tii = 0; $tii <= $#{$$tref[$ti]}; $tii++ ) {
                $EINDDESC[$ti][$tii] = $$tref[$ti][$tii] unless defined $EINDDESC[$ti][$tii];
            }
        }
    }
    unless( $tret ) {
        print( STDERR GetDatetime()."Errors detected. Detaching remaining threads...\n");
        for( $t++; $t < $NUMTHREADS; $t++ ) {
            $mythrds[$t]->detach();
        }
    }
} else {
    printf( STDERR GetDatetime()."Processing alignments...\n");
    my ($pointc) = (0);
    foreach $f( @alnfiles ) { 
        exit( 1 ) unless UpdateFreqs("$ALNDIR/$f", $ZSCORE, $PWS, \@EFF, \%PROS, \%PRSS, \%ALNS, \@FRQS, \@BFS, \$noalns, \@nopss, \%PAIRS);
        print( STDERR ".");
        if( ++$pointc%100 == 0 ) {
            print( STDERR "\n");
            $pointc = 0;
        }
    }
    print( STDERR "\n") unless( $pointc%100 == 0 );
}
##}}

printf( STDERR GetDatetime()." %d alignment(s) processed\n", $noalns );

my  $fd = *STDOUT;
if( open( OF, ">$OUTFILE" )) {
    $fd = *OF;
} else {
    print( STDERR "ERROR: Unable to open file for writing $OUTFILE; printing to stdout.\n" );
}

printf( $fd "ddms2s= %s\n", join(' ',@EFF));

for( $ii=0; $ii <= $#FRQS; $ii++ )
{
    $d = $EINDDESC[$ii][0];
    printf( STDERR GetDatetime()."  Deriving scores ($d)...\n");
    printf( STDERR GetDatetime()."    %d position(s) processed\n", $nopss[$ii]);

    next unless $nopss[$ii];
    undef %SCOS; $MINSC = $MAXSC = $ENTROPY = $EXPECTD = $IDV = 0;
    exit( 1 ) unless DeriveScores($NASC, $INF, $FRQS[$ii], $BFS[$ii], \%SCOS, \$MINSC, \$MAXSC, \$ENTROPY, \$EXPECTD, \$IDV);

    printf( STDERR GetDatetime()."    Printing scores...\n");

    printf( $fd "%s\n", $EINDDESC[$ii][1]);
    exit( 1 ) unless PrintScores($fd, $d, $noalns, $nopss[$ii], $NASC, $INF, $FRQS[$ii], \%SCOS, $MINSC, $MAXSC, $ENTROPY, $EXPECTD, $IDV);
}
printf( STDERR GetDatetime()."Finished.\n");
exit( 0 );

## -------------------------------------------------------------------
## -------------------------------------------------------------------

sub ThreadUpdateFreqs
{
    my $mytid = shift;
    my $mynthreads = shift;
    my $rafiles = shift;##ref. to alignment files
    my $rnoalns = shift;##ref. to no. alignments processed
    my $rnopss = shift;##ref. to array of no. positions
    my ($f, $i, $ret);
    my (%thres);
    $ret = 1;
    for($i = $mytid; $i <= $#$rafiles; $i+=$mynthreads) {
        $f = $$rafiles[$i];
        $ret = UpdateFreqs("$ALNDIR/$f", $ZSCORE, $PWS, \@EFF, \%PROS, \%PRSS, \%ALNS, \@FRQS, \@BFS, $rnoalns, $rnopss, \%PAIRS);
        last unless $ret;
    }
    $thres{RET} = $ret;
    $thres{FRQS} = \@FRQS;
    $thres{BFS} = \@BFS;
    $thres{ALNS} = $rnoalns;
    $thres{POSS} = $rnopss;
    $thres{EDSC} = \@EINDDESC;
    return( %thres );
}

## -------------------------------------------------------------------
## -------------------------------------------------------------------
## Calculate distance distribution match score for a pair of profile 
## positions
##
sub CalcPairDstDistrMatchScore
{
    my  $rqvec = shift;
    my  $rsvec = shift;
    my  $rpairlos = shift;
    ##die "ERROR: CalcPairDstDistrMatchScore: Invalid vector sizes." 
    ##    unless(0 <= $#$rqvec && 0 <= $#$rsvec);
    $$rpairlos = $DDMS->CalculateDDMS($rqvec, $rsvec);
    if($$rpairlos == $DDMS->GetErrorCode()) {
        die "ERROR: CalcPairDstDistrMatchScore: External method call failed.";
    }
}

## -------------------------------------------------------------------
## Get distance distribution for each profile position
##
sub GetDstDistrs
{
    my  $pfile = shift;
    my  $rpros = shift;
    my  $maxndsts = 127;
    my ($l, $c, $p, $r, $str, $name, $effnos, $ens );
    my ($vsz );

    open( PF, "<$pfile" ) or die "ERROR: GetDstDistrs: Unable to open file $pfile";

    $l = <PF>;
    return if $l !~ /^COTHER\s+profile/;

    $c = 0;
    while( <PF> ) {
        chomp;
        if( /^DESC:\s+(\S+)/) {
            $name = $1;
            next;
        }
        if( /^EFFNOS:\s+(\S+)/) {
            $effnos = $1;
            $$rpros{$name}[0][0] = $effnos;##data starts from position 1
            next;
        }
        if( /^(\d+)\s+([A-Z])\s+/) {
            die "ERROR: GetDstDistrs: No name defined in profile." unless $name;
            $c = 1;
            $p = $1;
            $r = $2;
            $$rpros{$name}[$p][0] = $r;
            next;
        }
        next unless $c;
        $c++;
        if( $c == 3 && /^\s+(\d+)/) {
            die "ERROR: GetDstDistrs: Undefined name in profile." unless $name;
            die "ERROR: GetDstDistrs: Undefined position." unless( $p && $r );
            $ens = sprintf("%.1f", $1/1000.0);##effective number of observations at the position
            $$rpros{$name}[$p][2] = $ens;
        }
        if( /^\s+D:\s+(.*)$/) {
            $str = $1;
            my @vals;
            die "ERROR: GetDstDistrs: Undefined name in profile." unless $name;
            die "ERROR: GetDstDistrs: Undefined position." unless( $p && $r );
            die "ERROR: GetDstDistrs: Undefined eff. no. observations at a position." unless($ens);
            @vals = map{int($_)} split(/\s+/,$str);
            push @vals, 0xff foreach(1..$DDMS->GetUnrollingDegree());
            $vsz = $vals[0];
            die "ERROR: GetDstDistrs: Inconsistent distribution size ($vsz) at line $. of $pfile" 
                if($vsz < 0 || $maxndsts < $vsz);
            die "ERROR: GetDstDistrs: Inconsistent distribution size at line $. of $pfile" 
                if($vsz + $DDMS->GetUnrollingDegree() != $#vals);

            ##skip this position's record if distances contain too small number of values
            next if $vsz < $DDMS->GetOptionDSTSEGM();

            $$rpros{$name}[$p][1] = ddms::new_intArray($#vals+1);#create an array for external routines
            ddms::intArray_setitem($$rpros{$name}[$p][1], $_, $vals[$_]) foreach(0..$#vals);
            ##delete_doubleArray($$rpros{$name}[$p][1]);#left for automatic destruction on exit
            $c = 0; undef $p; undef $r;
            undef $ens;
        }
    }

    close( PF );
    return 1;
}

## -------------------------------------------------------------------
## Update frequencies by analysing structural alignment files
##
sub UpdateFreqs
{
    my  $afile = shift;
    my  $zscore = shift;
    my  $pws = shift;
    my  $reff = shift;
    my  $rpros = shift;
    my  $rprss = shift;
    my  $ralns = shift;
    my  $rfrqs = shift;
    my  $rbfs = shift;
    my  $rnoalns = shift;
    my  $rnopss = shift;
    my  $rPAIRS = shift;
    my  $rec;
    my ($qsname, $qsscop, $qsbeg, $qs );
    my ($qSname, $qSscop, $qSbeg, $qSend, $qSzsc, $qS );
    my ($sSname, $sSscop, $sSbeg, $sSend, $sS );
    my ($ssname, $ssscop, $ssbeg, $ss );
    my  @tophits;

    unless( open( AF, "<$afile" )) {
        print STDERR GetDatetime()."ERROR: UpdateFreqs: Unable to open file $afile\n";
        return 0;
    }
    {
        local $/ = "\/\/\n";
        while( $rec = <AF> ) {
            if( $rec =~
                /^>(\S+)\s+(\S+)\s+(\d+)\s+SEQUENCE\s*\n
                   ([^>\/]+)\n
                  >(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+([\d\.eE\-\+]+)\s*\n
                   ([^>\/]+)\n
                  >(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s*\n
                   ([^>\/]+)\n
                  >(\S+)\s+(\S+)\s+(\d+)\s+SEQUENCE\s*\n
                   ([^>\/]+)\n
                  \/\/\n$
                /x )
            {
                $qsname = $1;  $qsscop = $2;  $qsbeg = $3;  $qs = $4;
                $qSname = $5;  $qSscop = $6;  $qSbeg = $7;  $qSend = $8;  $qSzsc = $9;  $qS = $10;
                $sSname = $11; $sSscop = $12; $sSbeg = $13; $sSend = $14; $sS = $15;
                $ssname = $16; $ssscop = $17; $ssbeg = $18; $ss = $19;

                $qs =~ s/\n//g; $qS =~ s/\n//g;
                $ss =~ s/\n//g; $sS =~ s/\n//g;

                unless(( $qsname eq $qSname )&&( $qsscop eq $qSscop )) { 
                    print STDERR "ERROR: UpdateFreqs: Inconsistent sequence-structure query names ($afile).\n"; 
                    return 0;
                }
                unless(( $ssname eq $sSname )&&( $ssscop eq $sSscop )) {
                    print STDERR "ERROR: UpdateFreqs: Inconsistent sequence-structure sbjct names ($afile).\n";
                    return 0;
                }
                unless(( length($qs) == length($qS))&&( length($ss) == length($sS))&&( length($qS) == length($ss))) {
                    print STDERR "ERROR: UpdateFreqs: Inconsistent sequence lengths ($afile).\n";
                    return 0;
                }
                unless( exists $$rpros{$qsname}) {
                    print STDERR "ERROR: UpdateFreqs: Name not found in profile library ($afile).\n";
                    return 0;
                }
                unless( exists $$rpros{$ssname}) {
                    print STDERR "ERROR: UpdateFreqs: Name not found in profile library ($afile).\n";
                    return 0;
                }
                unless( exists $$rprss{$qsname}) {
                    print STDERR "ERROR: UpdateFreqs: Name not found in struct. profile library ($afile).\n";
                    return 0;
                }
                unless( exists $$rprss{$ssname}) {
                    print STDERR "ERROR: UpdateFreqs: Name not found in struct. profile library ($afile).\n";
                    return 0;
                }

                ##next if $qSzsc < $zscore;
                last if $qSzsc < $zscore;
                if( $USINGTHREADS ) {
                    lock( %$ralns );
                    ##profile with predicted vs observed distances are aligned: no symmetry
                    #next if( exists( $$ralns{"$qSname$sSname"}) || exists( $$ralns{"$sSname$qSname"}));
                    next if( exists($$ralns{"$qSname$sSname"}));
                    $$ralns{"$qSname$sSname"} = 1;
                } else {
                    ##profile with predicted vs observed distances are aligned: no symmetry
                    #next if( exists( $$ralns{"$qSname$sSname"}) || exists( $$ralns{"$sSname$qSname"}));
                    next if( exists($$ralns{"$qSname$sSname"}));
                    $$ralns{"$qSname$sSname"} = 1;
                }

                push @tophits, [ $qsname, $qsscop, $qsbeg, $qs, 
                                 $qSname, $qSscop, $qSbeg, $qSend, $qSzsc, $qS,
                                 $sSname, $sSscop, $sSbeg, $sSend, $sS, 
                                 $ssname, $ssscop, $ssbeg, $ss ];
                $$rnoalns++;
            }
        }
    }
    close( AF );
    return 0 unless ProcessTopHits( $afile, $pws, $reff, $rpros, $rprss, $rfrqs, $rbfs, $rnopss, \@tophits, $rPAIRS );
    undef @tophits;
    return 1;
}

## -------------------------------------------------------------------
## process top hits of one query structure
##      
sub ProcessTopHits
{           
    my  $afile = shift;
    my  $pws = shift;
    my  $reff = shift;
    my  $rpros = shift;
    my  $rprss = shift;
    my  $rfrqs = shift;
    my  $rbfs = shift; 
    my  $rnopss = shift;
    my  $rtophits = shift;
    my  $rPAIRS = shift;
    my  @mypsqnwgs;

    return 0 unless MakeMultipleAlignment( $afile, $rtophits );
    return 0 unless VerifyMultipleAlignment( $afile, $rtophits );
    unless( $NOSW ) {
        return 0 unless CalculateSequenceWeights( $afile, $rtophits );
        if( $POS ) {
            return 0 unless DerivePositionalSequenceWeights( $afile, $rtophits, \@mypsqnwgs );
        }
    }
    if( $pws ) {
        return 0 unless ApplyWeightsPw( $afile, $reff, $rpros, $rprss, $rfrqs, $rbfs, $rnopss, $rtophits, $rPAIRS, \@mypsqnwgs );
    } else {
        return 0 unless ApplyWeights( $afile, $reff, $rpros, $rprss, $rfrqs, $rbfs, $rnopss, $rtophits, $rPAIRS, \@mypsqnwgs );
    }
    return 1;
}

## -------------------------------------------------------------------
## make multiple alignment given a set of pairwise 
## alignments with the same query
##
sub MakeMultipleAlignment
{
    my  $afile = shift;
    my  $rtophits = shift;
    my ($qsname, $qsscop, $qsbeg, $qs );
    my ($qSname, $qSscop, $qSbeg, $qSend, $qSzsc, $qS );
    my ($sSname, $sSscop, $sSbeg, $sSend, $sS );
    my ($ssname, $ssscop, $ssbeg, $ss );
    my ($pos, $qgap, $allmatch );
    my ($r, $c, $l, $len, $qsym, $sym );
    my ($minbeg, $maxend, $diff, @posits );

    $minbeg = 99999; ##large enough
    $maxend = 0;
    for( $r = 0; $r <= $#$rtophits; $r++ )
    {
        $qsname = $$rtophits[$r][0]; $qsscop = $$rtophits[$r][1]; $qsbeg = $$rtophits[$r][2]; $qs = \$$rtophits[$r][3];
        $qSname = $$rtophits[$r][4]; $qSscop = $$rtophits[$r][5]; $qSbeg = $$rtophits[$r][6]; $qSend = $$rtophits[$r][7]; 
                                     $qSzsc = $$rtophits[$r][8]; $qS = \$$rtophits[$r][9];
        $sSname = $$rtophits[$r][10]; $sSscop = $$rtophits[$r][11]; $sSbeg = $$rtophits[$r][12]; $sSend = $$rtophits[$r][13]; 
                                      $sS = \$$rtophits[$r][14];
        $ssname = $$rtophits[$r][15]; $ssscop = $$rtophits[$r][16]; $ssbeg = $$rtophits[$r][17]; $ss = \$$rtophits[$r][18];

        ##{{NOTE: to account for multiple chains that may implicitly appear in an alignment file
        ####TODO: MSA should be then built with respect to sequences (e.g. lacking a structure fragment)
        if( $qsbeg < $qSbeg ) {
            ##non-trivial instances when portions of structure lack
            ;#print(STDERR "WARNING: MakeMultipleAlignment: Structure sequence ".
             #            "beginning is greater than that of the sequence: $qsbeg<$qSbeg (sbjct: $sSname; $afile).\n");
        }
        #$diff = $qsbeg - $qSbeg;
        #$qSbeg = $qSbeg + $diff;
        #$qSend = $qSend + $diff;
        #$$rtophits[$r][6] = $qSbeg;
        #$$rtophits[$r][7] = $qSend;
        ##}}
        $posits[$r] = $qSbeg;
        $minbeg = $qSbeg if( $qSbeg < $minbeg );
        $maxend = $qSend if( $maxend < $qSend );
    }

    for( $c = $minbeg, $l = 0; $c <= $maxend; $l++ ) {
        undef $qsym;
        $qgap = 0;
        $allmatch = 1;
        for( $r = 0; $r <= $#$rtophits; $r++ )
        {
            $qsname = $$rtophits[$r][0]; $qsscop = $$rtophits[$r][1]; $qsbeg = $$rtophits[$r][2]; $qs = \$$rtophits[$r][3];
            $qSname = $$rtophits[$r][4]; $qSscop = $$rtophits[$r][5]; $qSbeg = $$rtophits[$r][6]; $qSend = $$rtophits[$r][7]; 
                                         $qSzsc = $$rtophits[$r][8]; $qS = \$$rtophits[$r][9];
            $sSname = $$rtophits[$r][10]; $sSscop = $$rtophits[$r][11]; $sSbeg = $$rtophits[$r][12]; $sSend = $$rtophits[$r][13]; 
                                          $sS = \$$rtophits[$r][14];
            $ssname = $$rtophits[$r][15]; $ssscop = $$rtophits[$r][16]; $ssbeg = $$rtophits[$r][17]; $ss = \$$rtophits[$r][18];

            $pos = $posits[$r];

            next if( $pos != $c );
            unless( $qsym ) {
                $qsym = substr( $$qS, $l, 1 );
                $qgap = 1 if $qsym eq '-';
                ##$qgap = 1 if( $qsym eq '-' && substr($$qs,$l,1) eq '-');##multiple chains; checking with respect to sequence
            }
            $sym = substr( $$qS, $l, 1 );
            ## added protection from livelock
            if( $sym ne '-' && $qsym ne '-' && uc($sym) ne uc($qsym)) {
                ##print( STDERR "ERROR: MakeMultipleAlignment: Inconsistent query sequences in alignment $qsym<=>$sym ".
                ##              "(pos: $l; sbjct: $sSname; $afile).\n" );
                ##return 0;
                print( STDERR "WARNING: MakeMultipleAlignment: Inconsistent sequences in alignment $qsym<=>$sym ".
                              "(pos: $l; sbjct: $sSname; $afile). Different chains likely. Skipped.\n" );
                splice(@$rtophits,$r,1);
                $r--;
                next;
            }
            do { $allmatch = 0; last; } if( uc($sym) ne uc($qsym));
        }

        for( $r = 0; $r <= $#$rtophits; $r++ )
        {
            $qsname = $$rtophits[$r][0]; $qsscop = $$rtophits[$r][1]; $qsbeg = $$rtophits[$r][2]; $qs = \$$rtophits[$r][3];
            $qSname = $$rtophits[$r][4]; $qSscop = $$rtophits[$r][5]; $qSbeg = $$rtophits[$r][6]; $qSend = $$rtophits[$r][7]; 
                                         $qSzsc = $$rtophits[$r][8]; $qS = \$$rtophits[$r][9];
            $sSname = $$rtophits[$r][10]; $sSscop = $$rtophits[$r][11]; $sSbeg = $$rtophits[$r][12]; $sSend = $$rtophits[$r][13]; 
                                          $sS = \$$rtophits[$r][14];
            $ssname = $$rtophits[$r][15]; $ssscop = $$rtophits[$r][16]; $ssbeg = $$rtophits[$r][17]; $ss = \$$rtophits[$r][18];

            $pos = $posits[$r];

            if( $pos != $c ) {
                if( $c < $pos ) {
                    substr( $$qS, $l, 1 )  = '-' . substr( $$qS, $l, 1 );
                    substr( $$sS, $l, 1 )  = '-' . substr( $$sS, $l, 1 );
                    substr( $$qs, $l, 1 )  = '-' . substr( $$qs, $l, 1 );
                    substr( $$ss, $l, 1 )  = '-' . substr( $$ss, $l, 1 );
                }
                else {
                    substr( $$qS, $l, 1 ) .= '-';
                    substr( $$qs, $l, 1 ) .= '-';
                    ##
                    if( substr( $$sS, $l, 1 ) ne '-' ) {
                            substr( $$sS, $l, 1 )  = '-' . substr( $$sS, $l, 1 );
                            substr( $$ss, $l, 1 )  = '-' . substr( $$ss, $l, 1 );
                    } else{
                            substr( $$sS, $l, 1 ) .= '-';
                            substr( $$ss, $l, 1 ) .= '-';
                    }
                }
                next;
            }

            next if $allmatch;

            $qsym = substr( $$qS, $l, 1 );
            $qgap = 1 if $qsym eq '-';
            ##$qgap = 1 if( $qsym eq '-' && substr($$qs,$l,1) eq '-');##multiple chains; checking with respect to sequence

            if( $qsym ne '-' ) {
                substr( $$qS, $l, 1 )  = '-' . substr( $$qS, $l, 1 );
                substr( $$sS, $l, 1 )  = '-' . substr( $$sS, $l, 1 );
                substr( $$qs, $l, 1 )  = '-' . substr( $$qs, $l, 1 );
                substr( $$ss, $l, 1 )  = '-' . substr( $$ss, $l, 1 );
            }
        }
        if( $allmatch ) {
            unless( $qgap ) {
                for( $r = 0; $r <= $#$rtophits; $r++ ) {
                    $qSend = $$rtophits[$r][7];
                    $pos = $posits[$r];
                    $posits[$r]++ if( $pos == $c && $c < $qSend );
                }
                $c++;
            }
        }
    }
##Note that gaps alingned to gaps indicate a possible lack of a 
##  structure fragment with respect to the sequence (e.g. astral)
##for( $r = 0; $r <= $#$rtophits; $r++ ){print STDERR "$$rtophits[$r][9]\n$$rtophits[$r][14]\n\n"}
    ## add terminal gaps if needed
    for( ;; $l++ ) {
        for( $r = 0; $r <= $#$rtophits; $r++ ) {
            $qS = \$$rtophits[$r][9];
            $sS = \$$rtophits[$r][14];
            last if $l < length( $$sS );
        }
        last if $#$rtophits < $r;
        for( $r = 0; $r <= $#$rtophits; $r++ ) {
            $qs = \$$rtophits[$r][3];
            $qS = \$$rtophits[$r][9];
            $sS = \$$rtophits[$r][14];
            $ss = \$$rtophits[$r][18];
            if( length( $$sS ) <= $l ) {
                substr( $$sS, $l, 1 ) .= '-';
                substr( $$qS, $l, 1 ) .= '-';
                substr( $$ss, $l, 1 ) .= '-';
                substr( $$qs, $l, 1 ) .= '-';
            }
        }
    }
    return 1;
}

## -------------------------------------------------------------------
## verify multiple alignment
##
sub VerifyMultipleAlignment
{
    my  $afile = shift;
    my  $rtophits = shift;
    my ($qS1, $qS1beg, $qS1end );
    my ($qs, $qS, $sS, $ss, $qSbeg, $qSend, $sSname );
    my ($r, $l, $qsym, $qsm1 );

    for( $r = 0; $r <= $#$rtophits; $r++ )
    {
        unless( $qS1 ) {
            $qS1 = \$$rtophits[$r][9];
            $qS1beg = $$rtophits[$r][6];
            $qS1end = $$rtophits[$r][7];
        }
        $qs = \$$rtophits[$r][3];
        $qS = \$$rtophits[$r][9];
        $sS = \$$rtophits[$r][14];
        $ss = \$$rtophits[$r][18];
        $qSbeg = $$rtophits[$r][6];
        $qSend = $$rtophits[$r][7];
        $sSname = $$rtophits[$r][10];

##print(STDOUT "$$qs\n$$qS\n$$sS\n$$ss\n//\n");

        if( length( $$qS ) != length( $$qs ) || length( $$qS ) != length( $$sS ) || length( $$qS ) != length( $$ss )) {
            printf( STDERR "ERROR: Lengths of query and subject are not equal ($afile).\n" );
            return 0;
        }
        if( length( $$qS ) != length( $$qS1 )) {
            printf( STDERR "ERROR: Lengths of query sequences are not equal ($afile).\n" );
            return 0;
        }

        for( $l = 0; $l < length( $$qS ); $l++ ) {
            $qsym = substr( $$qS, $l, 1 );
            $qsm1 = substr( $$qS1, $l, 1 );
            next if( $qsym eq '-' && $qsm1 eq '-' );
            if( $qsym ne '-' && $qsm1 ne '-' && uc($qsym) ne uc($qsm1)) {
                printf( STDERR "ERROR: VerifyMultipleAlignment: Inconsistent query sequences in alignment $qsm1<=>$qsym ".
                               "(pos: $l; sbjct: $sSname; $afile).\n$$qS1\n$$qS\n" );
                return 0;
            }
            if( $qsm1 eq '-' && $qsym ne '-' ) {
                substr( $$qS1, $l, 1 ) = uc($qsym);##uc will require subject's uc to (consider position)
            }
        }
        $$rtophits[0][6] = $qS1beg = $qSbeg if $qSbeg < $qS1beg;
        $$rtophits[0][7] = $qS1end = $qSend if $qS1end < $qSend;
    }
    return 1;
}

## -------------------------------------------------------------------
## Calculate sequence weights
##
sub CalculateSequenceWeights
{
    my  $afile = shift;
    my  $rtophits = shift;
    my ($qS1 );
    my ($qs, $qS, $sS, $ss, $sSname );
    my ($nadiffsyms, @column, @usdlen );
    my ($r, $p, $sym, $scd, $gwtsum );
    my  $MINUSDLEN = 30;
    my  $noeffress = 20;
    # Hash for conversion of amino acid to number
    my  %CODESAA = ( 'A' =>  0, 'R' =>  1, 'N' =>  2, 'D' =>  3, 'C' =>  4, 'Q' =>  5,
                     'E' =>  6, 'G' =>  7, 'H' =>  8, 'I' =>  9, 'L' => 10, 'K' => 11,
                     'M' => 12, 'F' => 13, 'P' => 14, 'S' => 15, 'T' => 16, 'W' => 17,
                     'Y' => 18, 'V' => 19, 'B' => 20, 'Z' => 21, 'J' => 22, 'O' => 23, 
                     'X' => 24
    );

    return 1 if $#$rtophits < 0;
    $qS1 = \$$rtophits[0][9];
    $usdlen[0] = 0;

    ##set bounds for each sequence
    for( $r = 0; $r <= $#$rtophits; $r++ )
    {
        $qs = \$$rtophits[$r][3];
        $qS = \$$rtophits[$r][9];
        $sS = \$$rtophits[$r][14];
        $ss = \$$rtophits[$r][18];

        if( length($$qS) < 1 || length($$qS) != length($$qs) || length($$qS) != length($$sS) || length($$qS) != length($$ss)) {
            printf( STDERR "ERROR: CalculateSequenceWeights: Invalid query/subject lengths ($afile).\n" );
            return 0;
        }
        $usdlen[$r+1] = 0;

        $$rtophits[$r][19] = $$rtophits[$r][20] = -1;
        for( $p = 0; $p < length($$qs); $p++ ) { do { $$rtophits[$r][19] = $p; last } if substr($$qs,$p,1) ne '-'; }
        for( $p = length($$qs)-1; 0 <= $p; $p-- ) { do { $$rtophits[$r][20] = $p; last } if substr($$qs,$p,1) ne '-'; }

        $$rtophits[$r][21] = $$rtophits[$r][22] = -1;
        for( $p = 0; $p < length($$qS); $p++ ) { do { $$rtophits[$r][21] = $p; last } if substr($$qS,$p,1) ne '-'; }
        for( $p = length($$qS)-1; 0 <= $p; $p-- ) { do { $$rtophits[$r][22] = $p; last } if substr($$qS,$p,1) ne '-'; }

        $$rtophits[$r][23] = $$rtophits[$r][24] = -1;
        for( $p = 0; $p < length($$sS); $p++ ) { do { $$rtophits[$r][23] = $p; last } if substr($$sS,$p,1) ne '-'; }
        for( $p = length($$sS)-1; 0 <= $p; $p-- ) { do { $$rtophits[$r][24] = $p; last } if substr($$sS,$p,1) ne '-'; }

        $$rtophits[$r][25] = $$rtophits[$r][26] = -1;
        for( $p = 0; $p < length($$ss); $p++ ) { do { $$rtophits[$r][25] = $p; last } if substr($$ss,$p,1) ne '-'; }
        for( $p = length($$ss)-1; 0 <= $p; $p-- ) { do { $$rtophits[$r][26] = $p; last } if substr($$ss,$p,1) ne '-'; }

        $$rtophits[$r][27] = $$rtophits[$r][28] = 0;##sequence weights
    }

    ##calculate weights
    for( $p = 0; $p < length($$qS1); $p++ ) {
        next if substr($$qS1,$p,1) eq '-';#insertions not aligned!
        $nadiffsyms = 0;
        undef @column;
        ##process query
        if( $$rtophits[0][21] <= $p && $p <= $$rtophits[0][22]) {
            $sym = substr($$qS1,$p,1);
            if( $sym ne '-') {
                $usdlen[0]++;
                unless( exists $CODESAA{uc($sym)}) {
                    printf( STDERR "ERROR: CalculateSequenceWeights: Unrecognized symbol in query: $sym ($afile).\n" );
                    return 0;
                }
                $scd = $CODESAA{uc($sym)};
                if( !exists $column[$scd] || $column[$scd] < 1 ) {
                    $column[$scd]++;
                    $nadiffsyms++;
                }
            }
        }
        ##process subjects
        for( $r = 0; $r <= $#$rtophits; $r++ )
        {
            $qs = \$$rtophits[$r][3];
            $qS = \$$rtophits[$r][9];
            $sS = \$$rtophits[$r][14];
            $ss = \$$rtophits[$r][18];
            $sSname = $$rtophits[$r][10];
            ##analyze subject's structure sequence sS
            next if( $p < $$rtophits[$r][23] || $$rtophits[$r][24] < $p );
            $sym = substr($$sS,$p,1);
            next if $sym eq '-';

            $usdlen[$r+1]++;
            unless( exists $CODESAA{uc($sym)}) {
                printf( STDERR "ERROR: CalculateSequenceWeights: Unrecognized symbol in sbjct $sSname: $sym ($afile).\n" );
                return 0;
            }
            $scd = $CODESAA{uc($sym)};
            if( !exists $column[$scd] || $column[$scd] < 1 ) {
                $column[$scd]++;
                $nadiffsyms++;
            }
        }
        $nadiffsyms = $noeffress if( $noeffress < $nadiffsyms );
        next if $nadiffsyms < 1;

        ##calculate weight for query
        if( $$rtophits[0][21] <= $p && $p <= $$rtophits[0][22]) {
            $sym = substr($$qS1,$p,1);
            if( $sym ne '-') {
                unless( exists $CODESAA{uc($sym)}) {
                    printf( STDERR "ERROR: CalculateSequenceWeights: Unrecognized symbol in query: $sym ($afile).\n" );
                    return 0;
                }
                $scd = $CODESAA{uc($sym)};
                if( !exists $column[$scd] || $column[$scd] < 1 ) {
                    printf( STDERR "ERROR: CalculateSequenceWeights: Invalid residue $sym distribution in query ($afile).\n" );
                    return 0;
                }
                $$rtophits[0][27] += 1.0 / ( $column[$scd] * $nadiffsyms );
            }
        }
        ##calculate weights for subjects
        for( $r = 0; $r <= $#$rtophits; $r++ )
        {
            $qs = \$$rtophits[$r][3];
            $qS = \$$rtophits[$r][9];
            $sS = \$$rtophits[$r][14];
            $ss = \$$rtophits[$r][18];
            $sSname = $$rtophits[$r][10];
            ##subject's structure sequence sS
            next if( $p < $$rtophits[$r][23] || $$rtophits[$r][24] < $p );
            $sym = substr($$sS,$p,1);
            next if $sym eq '-';

            unless( exists $CODESAA{uc($sym)}) {
                printf( STDERR "ERROR: CalculateSequenceWeights: Unrecognized symbol in sbjct $sSname: $sym ($afile).\n" );
                return 0;
            }
            $scd = $CODESAA{uc($sym)};
            if( !exists $column[$scd] || $column[$scd] < 1 ) {
                printf( STDERR "ERROR: CalculateSequenceWeights: Invalid residue $sym distribution in sbjct $sSname ($afile).\n" );
                return 0;
            }
            $$rtophits[$r][28] += 1.0 / ( $column[$scd] * $nadiffsyms );
        }
    }

    ##adjust weights by sequence lengths;
    ## each column has total weight =1
    $gwtsum = 0.0;
    unless( $$rtophits[0][27]) {
        printf( STDERR "ERROR: CalculateSequenceWeights: Null weight of query ($afile).\n" );
        return 0;
    }
    $$rtophits[0][27] /= max( $usdlen[0] + 1, $MINUSDLEN );
    $gwtsum += $$rtophits[0][27];
    for( $r = 0; $r <= $#$rtophits; $r++ ) {
        next unless( $$rtophits[$r][28]);
        $$rtophits[$r][28] /= max( $usdlen[$r+1] + 1, $MINUSDLEN );
        $gwtsum += $$rtophits[$r][28];
    }
    unless( $gwtsum ) {
        printf( STDERR "ERROR: CalculateSequenceWeights: Null weights ($afile).\n" );
        return 0;
    }
    ##normalize weights
    $$rtophits[0][27] /= $gwtsum;
    for( $r = 0; $r <= $#$rtophits; $r++ ) {
        $$rtophits[$r][28] /= $gwtsum if $$rtophits[$r][28];
    }
    return 1;
}

## -------------------------------------------------------------------
## Calculate sequence weights at each position of multiple alignment
##      
sub DerivePositionalSequenceWeights
{       
    my  $afile = shift;
    my  $rtophits = shift;
    my  $rpsqnwgs = shift;##positional sequence weights
    my ($qS1 );
    my ($qs, $qS, $sS, $ss, $sSname );
    my ($r, $p, $sym, $pwtsum );

    return 1 if $#$rtophits < 0;
    $qS1 = \$$rtophits[0][9];
    for( $p = 0; $p < length($$qS1); $p++ )
    {
        $pwtsum = 0.0;
        ## check query sequence
        if( $$rtophits[0][21] <= $p && $p <= $$rtophits[0][22]) {
            $sym = substr($$qS1,$p,1);
            $pwtsum = $$rtophits[0][27] if $sym ne '-';
        }
        ## positional weights for subject sequences
        for( $r = 0; $r <= $#$rtophits; $r++ ) {
            $$rpsqnwgs[$p][$r] = 0;##reset weights
            ##
            $qs = \$$rtophits[$r][3];
            $qS = \$$rtophits[$r][9];
            $sS = \$$rtophits[$r][14];
            $ss = \$$rtophits[$r][18];
            $sSname = $$rtophits[$r][10];
            ##subject's structure sequence sS
            next if( $p < $$rtophits[$r][23] || $$rtophits[$r][24] < $p );
            $sym = substr($$sS,$p,1);
            next if $sym eq '-';

            $pwtsum += $$rtophits[$r][28];
        }
        $$rpsqnwgs[$p][$r] = 0;##reset weight for query
        next unless $pwtsum;
        ## appropriately normalize positional weights
        for( $r = 0; $r <= $#$rtophits; $r++ ) {
            $sS = \$$rtophits[$r][14];
            $sSname = $$rtophits[$r][10];
            next if( $p < $$rtophits[$r][23] || $$rtophits[$r][24] < $p );
            $sym = substr($$sS,$p,1);
            next if $sym eq '-';
            if( $$rtophits[$r][28] <= 0.0 ) {
                printf( STDERR "ERROR: DerivePositionalSequenceWeights: Invalid sequence weight for sbjct $sSname: $sym ($afile).\n" );
                return 0;
            }
            $$rpsqnwgs[$p][$r] = $$rtophits[$r][28] / $pwtsum;
        }
        ## reserve last sequence position for query
        if( $$rtophits[0][21] <= $p && $p <= $$rtophits[0][22]) {
            $sym = substr($$qS1,$p,1);
            if( $sym ne '-') {
                if( $$rtophits[0][27] <= 0.0 ) {
                    printf( STDERR "ERROR: DerivePositionalSequenceWeights: Invalid sequence weight for query: $sym ($afile).\n" );
                    return 0;
                }
                $$rpsqnwgs[$p][$r] = $$rtophits[0][27] / $pwtsum;
            }
        }
    }
    return 1;
}

## -------------------------------------------------------------------
## apply sequence weights for each pair in the derived MSA when 
## calculating target frequencies
##
sub ApplyWeights
{
    my  $afile = shift;
    my  $reff = shift;
    my  $rpros = shift;
    my  $rprss = shift;
    my  $rfrqs = shift;
    my  $rbfs = shift;
    my  $rnopss = shift;
    my  $rtophits = shift;
    my  $rPAIRS = shift;
    my  $rpsqnwgs = shift;##positional sequence weights
    my ($rwght, $swght ) = (1,1);
    my ($rntmp, $sntmp );
    my ($rpostmp, $spostmp );
    my ($pairlos );
    my ($qsname_r, $qsbeg_r, $qs_r );
    my ($qSname_r, $qS_r );
    my ($sSname_r, $sS_r );
    my ($ssname_r, $ssbeg_r, $ss_r );
    my ($qsname_s, $qsbeg_s, $qs_s );
    my ($qSname_s, $qS_s );
    my ($sSname_s, $sS_s );
    my ($ssname_s, $ssbeg_s, $ss_s );
    my ($r, $s, $i, $len, $qsa, $qSa, $ssa, $sSa );
    my ($qpos, $spos);
    my ($rqvec, $rsvec );
    my ($tmp1, $tmp2, $pind );
    my ($ens1, $ens2 );
    my ($e, $e1, $e2, $ee, $ee1, $ee2, $ea, $eind );

    return 1 if $#$rtophits < 0;

    for( $r = -1; $r <= $#$rtophits; $r++ )
    {
        if( $r < 0 ) {
            $qsname_r = $$rtophits[0][0]; $qsbeg_r = $$rtophits[0][2]; $qs_r = \$$rtophits[0][3];
            $qSname_r = $$rtophits[0][4]; $qS_r = \$$rtophits[0][9];
            $sSname_r = $qSname_r; $sS_r = $qS_r;
            $ssname_r = $qsname_r; $ssbeg_r = $qsbeg_r; $ss_r = $qs_r;
            $rwght = $$rtophits[0][27] unless $NOSW;
        }
        else {
            $qsname_r = $$rtophits[$r][0]; $qsbeg_r = $$rtophits[$r][2]; $qs_r = \$$rtophits[$r][3];
            $qSname_r = $$rtophits[$r][4]; $qS_r = \$$rtophits[$r][9];
            $sSname_r = $$rtophits[$r][10]; $sS_r = \$$rtophits[$r][14];
            $ssname_r = $$rtophits[$r][15]; $ssbeg_r = $$rtophits[$r][17]; $ss_r = \$$rtophits[$r][18];
            $rwght = $$rtophits[$r][28] unless $NOSW;
        }
        unless( exists $$rpros{$qsname_r} && exists $$rprss{$qsname_r}) {
            print STDERR "ERROR: ApplyWeights: Name not found in profile library ($afile).\n";
            return 0;
        }
        unless( exists $$rpros{$ssname_r} && exists $$rprss{$ssname_r}) {
            print STDERR "ERROR: ApplyWeights: Name not found in profile library ($afile).\n";
            return 0;
        }
        if( $rwght <= 0.0 ) {
            print STDERR "ERROR: ApplyWeights: Invalid sequence $r weight ($afile).\n";
            return 0;
        }

        for( $s = $r+1; $s <= $#$rtophits; $s++ ) {

            $qsname_s = $$rtophits[$s][0]; $qsbeg_s = $$rtophits[$s][2]; $qs_s = \$$rtophits[$s][3];
            $qSname_s = $$rtophits[$s][4]; $qS_s = \$$rtophits[$s][9];
            $sSname_s = $$rtophits[$s][10]; $sS_s = \$$rtophits[$s][14];
            $ssname_s = $$rtophits[$s][15]; $ssbeg_s = $$rtophits[$s][17]; $ss_s = \$$rtophits[$s][18];
            $swght = $$rtophits[$s][28] unless $NOSW;

            unless( exists $$rpros{$qsname_s} && exists $$rprss{$qsname_s}) {
                print STDERR "ERROR: ApplyWeights: Name not found in profile library ($afile).\n";
                return 0;
            }
            unless( exists $$rpros{$ssname_s} && exists $$rprss{$ssname_s}) {
                print STDERR "ERROR: ApplyWeights: Name not found in profile library ($afile).\n";
                return 0;
            }
            if( $swght <= 0.0 ) {
                print STDERR "ERROR: ApplyWeights: Invalid sequence $s weight ($afile).\n";
                return 0;
            }

            $qpos = $ssbeg_r - 1;
            $spos = $ssbeg_s - 1;
            $len = length( $$ss_r );
            for( $i = 0; $i < $len; $i++ ) {
                $qsa = substr( $$ss_r, $i, 1);
                $qSa = substr( $$sS_r, $i, 1);
                $sSa = substr( $$sS_s, $i, 1);
                $ssa = substr( $$ss_s, $i, 1);
                $qpos++ if $qsa ne '-';
                $spos++ if $ssa ne '-';
                next if ( $qsa eq '-' || $qSa eq '-' || $sSa eq '-' || $ssa eq '-');
                ##if(( $qSa =~ /[a-z]/ && $sSa =~ /[A-Z]/ )||( $qSa =~ /[A-Z]/ && $sSa =~ /[a-z]/ )) {
                ##    CAN BE THE CASE BECAUSE OF INDIRECT ALIGNMENTS IN CONSTRUCTED MULTIPLE SEQUENCE ALIGNMENT
                ##    printf( STDERR "WARNING: ApplyWeights: Unmatched character case in alignment %s vs. %s (%s).\n", $qSname, $sSname, $afile );
                ##    next;
                ##}

                ## get positional weights
                if( !$NOSW && $POS ) {
                    if( $r < 0 ) { $rwght = $$rpsqnwgs[$i][$#$rtophits+1]; }## query positional weights given at last sequence position
                        else {     $rwght = $$rpsqnwgs[$i][$r]; }
                    $swght = $$rpsqnwgs[$i][$s];
                    if( $rwght <= 0.0 ) {
                        print STDERR "ERROR: ApplyWeights: Invalid positional sequence weights at $i: $rwght ($afile).\n";
                        return 0;
                    }
                    if( $swght <= 0.0 ) {
                        print STDERR "ERROR: ApplyWeights: Invalid positional sequence weights at $i: $swght ($afile).\n";
                        return 0;
                    }
                }

                unless( exists $$rpros{$ssname_r}[$qpos]) {
                    print STDERR "ERROR: ApplyWeights: Undefined position $qpos in profile $ssname_r ($afile).\n";
                    return 0;
                }
                #unless( exists $$rpros{$ssname_s}[$spos]) {
                unless( exists $$rprss{$ssname_s}[$spos]) {
                    print STDERR "ERROR: ApplyWeights: Undefined position $spos in profile $ssname_s ($afile).\n";
                    return 0;
                }
                unless( $$rpros{$ssname_r}[$qpos][0] eq $qsa ) {
                    print STDERR "ERROR: ApplyWeights: Unmatched residue character at $qpos in query $ssname_r sequence-structure data ($afile).\n";
                    return 0;
                }
                #unless( $$rpros{$ssname_s}[$spos][0] eq $ssa ) {
                unless( $$rprss{$ssname_s}[$spos][0] eq $ssa ) {
                    print STDERR "ERROR: ApplyWeights: Unmatched residue character at $spos in sbjct $ssname_s sequence-structure data ($afile).\n";
                    return 0;
                }

                ##omit position if profile's eff. no. sequences is not sufficient
                if( $$rpros{$ssname_r}[0][0] < 0 || 20 < $$rpros{$ssname_r}[0][0] ||
                    $$rprss{$ssname_s}[0][0] < 0 || 20 < $$rprss{$ssname_s}[0][0] ) {
                    print STDERR "ERROR: ApplyWeights: Invalid profile's effnos: $$rpros{$ssname_r}[0][0], $$rprss{$ssname_s}[0][0] ".
                                 "($ssname_r vs. $ssname_s; $afile).\n";
                    return 0;
                }
                next if( $$rpros{$ssname_r}[0][0] < $$reff[0] ||
                         $$rprss{$ssname_s}[0][0] < $$reff[0] );

                next unless exists $$rpros{$ssname_r}[$qpos][1];
                next unless exists $$rprss{$ssname_s}[$spos][1];

                $rqvec = $$rpros{$ssname_r}[$qpos][1];
                #$rsvec = $$rpros{$ssname_s}[$spos][1];
                $rsvec = $$rprss{$ssname_s}[$spos][1];

                ##unless(0 <= $#$rqvec && 0 <= $#$rsvec) {
                ##    print STDERR "ERROR: ApplyWeights: Invalid vectors: sizes: $#$rqvec, $#$rsvec ($ssname_r vs. $ssname_s; $afile).\n";
                ##    return 0;
                ##}

                ## check for calculated score
                ##$rntmp = $ssname_r; $rpostmp = $qpos;
                ##$sntmp = $ssname_s; $spostmp = $spos;
                ##if( $sntmp lt $rntmp ) { 
                ##    $rntmp = $ssname_s; $rpostmp = $spos; 
                ##    $sntmp = $ssname_r; $spostmp = $qpos; 
                ##}
                ##if( exists $$rPAIRS{$rntmp}{$sntmp} && exists $$rPAIRS{$rntmp}{$sntmp}{$rpostmp}{$spostmp}) {
                ##    $pairlos = $$rPAIRS{$rntmp}{$sntmp}{$rpostmp}{$spostmp};
                ##}
                ##else {
                    CalcPairDstDistrMatchScore($rqvec,$rsvec,\$pairlos);
                    ##$$rPAIRS{$rntmp}{$sntmp}{$rpostmp}{$spostmp} = $pairlos;##too much memory required
                ##}
                next if $pairlos == $DDMS->GetIgnoreCode();
                ##NOTE: rounding score moved to CalcPairDstDistrMatchScore: 
                ##$tmp1 = $pairlos-int($pairlos);
                ##if( 0 < $tmp1 ) {
                ##    if( $tmp1 < 0.25 ) { $pairlos = int($pairlos); }
                ##    elsif( $tmp1 < 0.75 ) { $pairlos = int($pairlos)+0.5; }
                ##    else { $pairlos = int($pairlos)+1.; }
                ##}
                ##elsif( $tmp1 < 0 ) {
                ##    if( -0.25 < $tmp1 ) { $pairlos = int($pairlos); }
                ##    elsif( -0.75 < $tmp1 ) { $pairlos = int($pairlos)-0.5; }
                ##    else { $pairlos = int($pairlos)-1.; }
                ##}



                ## table index regarding thresholds of eff. no. observations
                ##
                $ens1 = $$rpros{$ssname_r}[$qpos][2];
                #$ens2 = $$rpros{$ssname_s}[$spos][2];
                $ens2 = $$rprss{$ssname_s}[$spos][2];

                ##$ens1 = $$rpros{$ssname_r}[0][0]; ##
                ####$ens2 = $$rpros{$ssname_s}[0][0]; ##
                ##$ens2 = $$rprss{$ssname_s}[0][0]; ##

                if( $ens1 < 0 || 20 < $ens1  ||  $ens2 < 0 || 20 < $ens2 ) {
                    print STDERR "ERROR: ApplyWeights: Invalid effnos at profile position: $ens1, $ens2 ($ssname_r vs. $ssname_s; $afile).\n";
                    return 0;
                }

                next if( $ens1 < $$reff[0] ||  $ens2 < $$reff[0] );

                do{ $e = $ens1; $ens1 = $ens2; $ens2 = $e; } if( $ens2 < $ens1 );

                for( $e = 0; $e <= $#$reff; $e++ ) {
                    $e1 = $$reff[$e];
                    $e2 = 20.01;
                    $e2 = $$reff[$e+1] if $e < $#$reff;
                    next unless( $e1 <= $ens1 && $ens1 < $e2 );
                    last;
                }
                for( $ee = $e; $ee <= $#$reff; $ee++ ) {
                    $ee1 = $$reff[$ee];
                    $ee2 = 20.01;
                    $ee2 = $$reff[$ee+1] if $ee < $#$reff;
                    next unless( $ee1 <= $ens2 && $ens2 < $ee2 );
                    last;
                }
                if( $#$reff < $e || $#$reff < $ee || $ee < $e ) {
                    print STDERR "ERROR: ApplyWeights: Segment covering ENO not found: $ens1, $ens2 ($ssname_r vs. $ssname_s; $afile).\n";
                    return 0;
                }
                $eind = ($e*(2*($#$reff+1)-$e+1))/2 + $ee-$e;##index calculated as sum of arithm. series

                $e2 = '-' if 20 < $e2;
                $ee2 = '-' if 20 < $ee2;
                unless( defined( $EINDDESC[$eind])) {
                    $EINDDESC[$eind][0] = "[$e1,$e2) vs. [$ee1,$ee2)";
                    $EINDDESC[$eind][1] = "$e1$ee1";
                }

                for( ; 0<=$e; $e-- ) {
                    for( $ea=$ee; 0<=$ea && $e<=$ea; $ea-- ) {
                        $eind = ($e*(2*($#$reff+1)-$e+1))/2 + $ea-$e;
                        ##
                        $pind = 0;##positive
                        $pind = 1 if( $qSa =~ /[a-z]/ || $sSa =~ /[a-z]/);
                        $$rfrqs[$eind][$pind]{$pairlos} += $rwght * $swght;
                        $$rbfs[$eind][$pind] += $rwght * $swght;
                        $$rnopss[$eind]++;
                        last unless $OVL;
                    }
                    last unless $OVL;
                }
            }
        }
    }
    return 1;
}

## -------------------------------------------------------------------
## apply sequence weights pairwise in calculating target frequencies
##
sub ApplyWeightsPw
{
    my  $afile = shift;
    my  $reff = shift;
    my  $rpros = shift;
    my  $rprss = shift;
    my  $rfrqs = shift;
    my  $rbfs = shift;
    my  $rnopss = shift;
    my  $rtophits = shift;
    my  $rPAIRS = shift;
    my  $rpsqnwgs = shift;##positional sequence weights
    my ($noths, $qwght, $swght ) = (0,1,1);
    my ($rntmp, $sntmp );
    my ($rpostmp, $spostmp );
    my ($pairlos );
    my ($qsname, $qsscop, $qsbeg, $qs );
    my ($qSname, $qSscop, $qSbeg, $qSend, $qSzsc, $qS );
    my ($sSname, $sSscop, $sSbeg, $sSend, $sS );
    my ($ssname, $ssscop, $ssbeg, $ss );
    my ($r, $i, $len, $qsa, $qSa, $ssa, $sSa );
    my ($qpos, $spos );
    my ($rqvec, $rsvec );
    my ($tmp1, $tmp2, $pind );
    my ($ens1, $ens2 );
    my ($e, $e1, $e2, $ee, $ee1, $ee2, $ea, $eind );

    $noths = $#$rtophits + 1;
    $qwght = 1.0 / $noths if( 0 < $noths && !$NOSW );

    for( $r = 0; $r <= $#$rtophits; $r++ )
    {
        $qsname = $$rtophits[$r][0]; $qsscop = $$rtophits[$r][1]; $qsbeg = $$rtophits[$r][2]; $qs = \$$rtophits[$r][3];
        $qSname = $$rtophits[$r][4]; $qSscop = $$rtophits[$r][5]; $qSbeg = $$rtophits[$r][6]; $qSend = $$rtophits[$r][7]; 
                                     $qSzsc = $$rtophits[$r][8]; $qS = \$$rtophits[$r][9];
        $sSname = $$rtophits[$r][10]; $sSscop = $$rtophits[$r][11]; $sSbeg = $$rtophits[$r][12]; $sSend = $$rtophits[$r][13]; 
                                      $sS = \$$rtophits[$r][14];
        $ssname = $$rtophits[$r][15]; $ssscop = $$rtophits[$r][16]; $ssbeg = $$rtophits[$r][17]; $ss = \$$rtophits[$r][18];

        unless( exists $$rpros{$qsname} && exists $$rprss{$qsname}) {
            print STDERR "ERROR: ApplyWeightsPw: Name not found in profile library ($afile).\n";
            return 0;
        }
        unless( exists $$rpros{$ssname} && exists $$rprss{$ssname}) {
            print STDERR "ERROR: ApplyWeightsPw: Name not found in profile library ($afile).\n";
            return 0;
        }

        $qwght = $$rtophits[0][27] unless $NOSW;
        $swght = $$rtophits[$r][28] unless $NOSW;

        $qpos = $qsbeg - 1;
        $spos = $ssbeg - 1;
        $len = length( $$qs );
        for( $i = 0; $i < $len; $i++ ) {
            $qsa = substr( $$qs, $i, 1);
            $qSa = substr( $$qS, $i, 1);
            $sSa = substr( $$sS, $i, 1);
            $ssa = substr( $$ss, $i, 1);
            $qpos++ if $qsa ne '-';
            $spos++ if $ssa ne '-';
            next if ( $qsa eq '-' || $qSa eq '-' || $sSa eq '-' || $ssa eq '-');
            if(( $qSa =~ /[a-z]/ && $sSa =~ /[A-Z]/ )||( $qSa =~ /[A-Z]/ && $sSa =~ /[a-z]/ )) {
                printf( STDERR "WARNING: ApplyWeightsPw: Unmatched character case in alignment %s vs. %s (%s).\n", $qSname, $sSname, $afile );
                next;
            }

            ## get positional weights
            if( !$NOSW && $POS ) {
                $qwght = $$rpsqnwgs[$i][$#$rtophits+1]; ## query positional weights given at last sequence position
                $swght = $$rpsqnwgs[$i][$r];
                if( $qwght <= 0.0 ) {
                    print STDERR "ERROR: ApplyWeightsPw: Invalid positional sequence weights at $i: $qwght ($afile).\n";
                    return 0;
                }
                if( $swght <= 0.0 ) {
                    print STDERR "ERROR: ApplyWeightsPw: Invalid positional sequence weights at $i: $swght ($afile).\n";
                    return 0;
                }
            }

            unless( exists $$rpros{$qsname}[$qpos]) {
                print STDERR "ERROR: ApplyWeightsPw: Undefined position $qpos in profile $qsname ($afile).\n";
                return 0;
            }
            #unless( exists $$rpros{$ssname}[$spos]) {
            unless( exists $$rprss{$ssname}[$spos]) {
                print STDERR "ERROR: ApplyWeightsPw: Undefined position $spos in profile $ssname ($afile).\n";
                return 0;
            }
            unless( $$rpros{$qsname}[$qpos][0] eq $qsa ) {
                print STDERR "ERROR: ApplyWeightsPw: Unmatched residue character at $qpos in query $qsname sequence-structure data ($afile).\n";
                return 0;
            }
            #unless( $$rpros{$ssname}[$spos][0] eq $ssa ) {
            unless( $$rprss{$ssname}[$spos][0] eq $ssa ) {
                print STDERR "ERROR: ApplyWeightsPw: Unmatched residue character at $spos in sbjct $ssname sequence-structure data ($afile).\n";
                return 0;
            }

            ##omit position if profile's eff. no. sequences is not sufficient
            if( $$rpros{$qSname}[0][0] < 0 || 20 < $$rpros{$qSname}[0][0] ||
                $$rprss{$sSname}[0][0] < 0 || 20 < $$rprss{$sSname}[0][0] ) {
                print STDERR "ERROR: ApplyWeights: Invalid profile's effnos: $$rpros{$qSname}[0][0], $$rprss{$sSname}[0][0] ".
                             "($qsname vs. $ssname; $afile).\n";
                return 0;
            }
            next if( $$rpros{$qSname}[0][0] < $$reff[0] ||
                     $$rprss{$sSname}[0][0] < $$reff[0] );

            next unless exists $$rpros{$qSname}[$qpos][1];
            next unless exists $$rprss{$sSname}[$spos][1];

            $rqvec = $$rpros{$qSname}[$qpos][1];
            #$rsvec = $$rpros{$sSname}[$spos][1];
            $rsvec = $$rprss{$sSname}[$spos][1];

            ##unless(0 <= $#$rqvec && 0 <= $#$rsvec) {
            ##    print STDERR "ERROR: ApplyWeightsPw: Invalid vectors: sizes: $#$rqvec, $#$rsvec ($qsname vs. $ssname; $afile).\n";
            ##    return 0;
           ## }

            ## check for calculated score
            ##$rntmp = $qSname; $rpostmp = $qpos;
            ##$sntmp = $sSname; $spostmp = $spos;
            ##if( $sntmp lt $rntmp ) {
            ##    $rntmp = $sSname; $rpostmp = $spos;
            ##    $sntmp = $qSname; $spostmp = $qpos;
            ##}
            ##if( exists $$rPAIRS{$rntmp}{$sntmp} && exists $$rPAIRS{$rntmp}{$sntmp}{$rpostmp}{$spostmp}) {
            ##    $pairlos = $$rPAIRS{$rntmp}{$sntmp}{$rpostmp}{$spostmp};
            ##}
            ##else {
                CalcPairDstDistrMatchScore($rqvec,$rsvec,\$pairlos);
                ##$$rPAIRS{$rntmp}{$sntmp}{$rpostmp}{$spostmp} = $pairlos;
            ##}
            next if $pairlos == $DDMS->GetIgnoreCode();
            ##NOTE: rounding score moved to CalcPairDstDistrMatchScore:
            ##$tmp1 = $pairlos-int($pairlos);
            ##if( 0 < $tmp1 ) {
            ##    if( $tmp1 < 0.25 ) { $pairlos = int($pairlos); }
            ##    elsif( $tmp1 < 0.75 ) { $pairlos = int($pairlos)+0.5; }
            ##    else { $pairlos = int($pairlos)+1.; }
            ##}
            ##elsif( $tmp1 < 0 ) {
            ##    if( -0.25 < $tmp1 ) { $pairlos = int($pairlos); }
            ##    elsif( -0.75 < $tmp1 ) { $pairlos = int($pairlos)-0.5; }
            ##    else { $pairlos = int($pairlos)-1.; }
            ##}



            ## table index regarding thresholds of eff. no. sequences
            ##
            $ens1 = $$rpros{$qSname}[$qpos][2];
            #$ens2 = $$rpros{$sSname}[$spos][2];
            $ens2 = $$rprss{$sSname}[$spos][2];

            ##$ens1 = $$rpros{$qSname}[0][0]; ##
            ####$ens2 = $$rpros{$sSname}[0][0]; ##
            ##$ens2 = $$rprss{$sSname}[0][0]; ##

            if( $ens1 < 0 || 20 < $ens1  ||  $ens2 < 0 || 20 < $ens2 ) {
                print STDERR "ERROR: ApplyWeightsPw: Invalid effnos at profile position: $ens1, $ens2 ($qsname vs. $ssname; $afile).\n";
                return 0;
            }

            next if( $ens1 < $$reff[0] ||  $ens2 < $$reff[0] );

            do{ $e = $ens1; $ens1 = $ens2; $ens2 = $e; } if( $ens2 < $ens1 );

            for( $e = 0; $e <= $#$reff; $e++ ) {
                $e1 = $$reff[$e];
                $e2 = 20.01;
                $e2 = $$reff[$e+1] if $e < $#$reff;
                next unless( $e1 <= $ens1 && $ens1 < $e2 );
                last;
            }
            for( $ee = $e; $ee <= $#$reff; $ee++ ) {
                $ee1 = $$reff[$ee];
                $ee2 = 20.01;
                $ee2 = $$reff[$ee+1] if $ee < $#$reff;
                next unless( $ee1 <= $ens2 && $ens2 < $ee2 );
                last;
            }
            if( $#$reff < $e || $#$reff < $ee || $ee < $e ) {
                print STDERR "ERROR: ApplyWeightsPw: Segment covering eff. no. sequences not found: $ens1, $ens2 ($qsname vs. $ssname; $afile).\n";
                return 0;
            }
            $eind = ($e*(2*($#$reff+1)-$e+1))/2 + $ee-$e;##index calculated as sum of arithm. series

            $e2 = '-' if 20 < $e2;
            $ee2 = '-' if 20 < $ee2;
            unless( defined( $EINDDESC[$eind])) {
                $EINDDESC[$eind][0] = "[$e1,$e2) vs. [$ee1,$ee2)";
                $EINDDESC[$eind][1] = "$e1$ee1";
            }

            for( ; 0<=$e; $e-- ) {
                for( $ea=$ee; 0<=$ea && $e<=$ea; $ea-- ) {
                    $eind = ($e*(2*($#$reff+1)-$e+1))/2 + $ea-$e;
                    ##
                    $pind = 0;##positive
                    $pind = 1 if( $qSa =~ /[a-z]/ || $sSa =~ /[a-z]/);
                    $$rfrqs[$eind][$pind]{$pairlos} += $qwght * $swght;
                    $$rbfs[$eind][$pind] += $qwght * $swght;
                    $$rnopss[$eind]++;
                    last unless $OVL;
                }
                last unless $OVL;
            }
        }
    }
    return 1;
}

## -------------------------------------------------------------------
## derive scores given target frequencies
##
sub DeriveScores
{
    my  $nasc = shift;
    my  $inf = shift;
    my  $rfrqs = shift;
    my  $rbfs = shift;
    my  $rscos = shift;
    my  $rminsc = shift;
    my  $rmaxsc = shift;
    my  $rentropy = shift;
    my  $rexpectd = shift;
    my  $ridv = shift;
    my ($possum, $negsum ) = (0,0);
    my ($v1, $v2 ) = (0,0);
    my ($odds, $score );
    my ($i, $j );

    foreach $i( keys %{$$rfrqs[0]}) {
        $possum += $$rfrqs[0]{$i};
    }
    foreach $i( keys %{$$rfrqs[1]}) {
        $negsum += $$rfrqs[1]{$i};
    }

    if( $possum <= 0.0 || $negsum <= 0.0 ) {
        printf( STDERR "ERROR: DeriveScores: Invalid sums of target frequencies: $possum, $negsum.\n");
        return 0;
    }

    ##normalize frequencies
    foreach $i( keys %{$$rfrqs[0]}) {
        $$rfrqs[0]{$i} /= $possum;
    }
    foreach $i( keys %{$$rfrqs[1]}) {
        $$rfrqs[1]{$i} /= $negsum;
    }

    ##calculate scores
    $$rminsc = 0.0;
    $$rmaxsc = 0.0;
    $$rentropy = 0.0;
    $$rexpectd = 0.0;
    $$ridv = 0.0;
    foreach $i( keys %{$$rfrqs[0]}) {
        $odds = 0.0;
        $score = $inf;

        if( exists $$rfrqs[1]{$i}) {
            $odds = $$rfrqs[0]{$i} / $$rfrqs[1]{$i};
        }
        if( 0.0 < $odds ) {
            $score = log( $odds );
            $$rminsc = $score if $score < $$rminsc;
            $$rmaxsc = $score if $$rmaxsc < $score;

            $v1 = $score * $$rfrqs[0]{$i};
            $v2 = $score * $$rfrqs[1]{$i};

            $$rentropy += $v1;
            $$rexpectd += $v2;

            $$ridv += (($$rfrqs[0]{$i} < $$rfrqs[1]{$i})? $$rfrqs[0]{$i}: $$rfrqs[1]{$i});
        }
        $$rscos{$i} = $score;
    }
    return 1;
}

## -------------------------------------------------------------------
## print scores and related information
##
sub PrintScores
{
    my  $fd = shift;
    my  $desc = shift;##description
    my  $noalns = shift;
    my  $nopss = shift;
    my  $nasc = shift;
    my  $inf = shift;
    my  $rfrqs = shift;
    my  $rscos = shift;
    my  $minsc = shift;
    my  $maxsc = shift;
    my  $entropy = shift;
    my  $expectd = shift;
    my  $idv = shift;
    my ($i, @skeys );

    printf( $fd "## Scores from processing %d alignments, %d positions\n", $noalns, $nopss );
    printf( $fd "##  (Profile pairs processed with respect to ENO: %s)\n", $desc );
    printf( $fd "## Entropy = %.4f, Expected = %.4f\n", $entropy, $expectd );
    printf( $fd "## Minimum = %.4f, Maximum = %.4f\n", $minsc, $maxsc );
    printf( $fd "## Integrated distribution overlap: %.4f\n", $idv);
    printf( $fd "## Parameters: Z-score = %.1f (argument: DALI Z-score threshold)\n", $ZSCORE);
    printf( $fd "## Parameters: DSTSEGM = %d (option: indivisible fragment of consecutive values in DP)\n", 
            $DDMS->GetOptionDSTSEGM());
    if($DDMS->GetDScoreApprox()==2) {
        printf( $fd "## Parameters: Approximation by two-level linear interpolation: k1= %.2f k2= %.2f\n",
                $DDMS->GetDScoreApproxScore00(), $DDMS->GetDScoreApproxScore200());
    } elsif($DDMS->GetDScoreApprox()==1) {
        printf( $fd "## Parameters: Approximation by linear interpolation = yes\n"); 
    } else {
        printf( $fd "## Parameters: Approximation by linear interpolation = no\n"); 
    }
    printf( $fd "## Parameters: MaxDst = %d (maximum distance value)\n", 
            $DDMS->GetDScoreMaxDst());
    printf( $fd "## Parameters: Theta = %.2f ".
            "(allowed percentage variation from average for the score to be positive)\n", 
            $DDMS->GetDScoreTheta());
    printf( $fd "## Parameters: AbsDiffExp = %.2f ".
            "(exponent for the absolute difference between distances)\n", 
            $DDMS->GetDScoreAbsDiffExp());
    printf( $fd "## Parameters: AvgDistExp = %.2f (exponent for the average distance)\n", 
            $DDMS->GetDScoreAvgDistExp());
    printf( $fd "## Parameters: Granularity = %.0f (multiplication factor for scores)\n", 
            $DDMS->GetDScoreGranularity());
    printf( $fd "## Parameters: NSEGM = %d (#ungapped maximum-scoring DP segments)\n",
            $DDMS->GetDPNSegments());
    printf( $fd "## Parameters: HEURISTIC = %d (linear inpterpolation for half of DP cells)\n",
            $DDMS->GetDPHeuristic());
    printf( $fd "## Parameters: STEP = %d (processing sparsity, step over DP diagonals)\n",
            $DDMS->GetDPStep());
    printf( $fd "## Parameters: BAND_FRAC = %d (band fraction, expressed as x to 1/2^x)\n",
            $DDMS->GetDPBandFraction());
    printf( $fd "## Parameters: UNROLL = %d (unrolling degree used for DP diagonals)\n",
            $DDMS->GetUnrollingDegree());
    printf( $fd "Cardinality = %d\n", scalar(keys %$rscos));

    @skeys = sort { $a <=> $b } keys %$rscos;
    foreach $i( @skeys ) {
        printf( $fd " %8.1f", $i);
    }
    printf( $fd "\n");
    foreach $i( @skeys ) {
        if( exists $$rscos{$i} && $$rscos{$i} < $inf ) {
            printf( $fd " %8.4f", $$rscos{$i});
        }
        else {
            printf( $fd " %8s", '+');
        }
    }
    printf( $fd "\n\n##");
    foreach $i( @skeys ) {
        if( exists $$rfrqs[0]{$i}) {
            printf( $fd " %8.4f", $$rfrqs[0]{$i});
        }
        else {
            printf( $fd " %8s", 'NA');
        }
    }
    printf( $fd "\n##");
    foreach $i( @skeys ) {
        if( exists $$rfrqs[1]{$i}) {
            printf( $fd " %8.4f", $$rfrqs[1]{$i});
        }
        else {
            printf( $fd " %8s", 'NA');
        }
    }
    printf( $fd "\n\n");
    return 1;
}

## -------------------------------------------------------------------
## read file list from directory
##
sub ReadFiles {
    my  $dirname = shift;
    my  $refiles = shift;

    opendir( DIR, $dirname ) || die "ERROR: Cannot open directory $dirname.";

    @{$refiles} = grep { -f "$dirname/$_" } readdir( DIR );

    ##sort files by size
    @$refiles =
        map  { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map  { [$_, -s "$dirname/$_"] }
        @$refiles;

    closedir( DIR );
}

## -------------------------------------------------------------------
