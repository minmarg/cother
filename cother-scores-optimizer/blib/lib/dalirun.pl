#!/usr/bin/perl

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use nwalign;


##my  $DALIPROG = '/usr/local/bin/DaliLite';
##my  $DALIPROG = '/usr/local/install/DaliLite/DaliLite';
my  $DALIPROG = '/home/mindaugas/installed/DaliLite/DaliLite';

my  $SUPOSFILE = 'CA_1.pdb';
my  $ALIGNFILE = 'index.html';
my  $MATRIXFILE = 'matrix.txt';
my  $RANGESFILE = 'ranges.txt';

my  $LOGFILE = 'dali.log';
my  $ERRFILE = 'dali.err';
my  $LOCKFILE= 'dali.lock';

my  $ALNSTART = '<PRE>';
my  $ALNSTOP = qr{<\/PRE>};

my  $usage = "
Align the protein structure with a given set of protein structures 
by DALI. Sequences and structures are mutually aligned too.
2013(C)Mindaugas Margelevicius,VU IBT,Vilnius


Usage:
$0 <Options> 

Options:

<pdbfile>  - name of query pdb (domain) file
<pdbdir>   - name of directory of pdb files to compare
             <pdbfile> against
<fastadir> - name of directory of corresponding fasta
             files of sequences
<outdir>   - name of output directory to put results in
             (directory must exist)
<outfile>  - name of output file to write fasta 
             alignments in
             default=<name_determined_by_input>

";

die "ERROR: Too few parameters provided.\n$usage" if $#ARGV < 3;

my  $INPUTPDB = $ARGV[0];
my  $PDBREPOS = $ARGV[1];
my  $FASTADIR = $ARGV[2];
my  $OUTPUTDIR = $ARGV[3];
my  $outputname; $outputname = $ARGV[4] if 3 < $#ARGV;
my  $OUTPUT = \*STDOUT;

my  $NAMEPAT = 'top_hits_dali';
my  $EXT = 'out';

-f $DALIPROG || die "ERROR: Dali executable $DALIPROG does not exist.";

-f $INPUTPDB || die "ERROR: Input pdb file $INPUTPDB does not exist.";
-d $PDBREPOS || die "ERROR: Directory of pdb files $PDBREPOS does not exist.";
-d $FASTADIR || die "ERROR: Directory of fasta files $FASTADIR does not exist.";
-d $OUTPUTDIR || die "ERROR: Output directory $OUTPUTDIR does not exist.";


## ===================================================================
## main globals
##

my  @Files;
my  @Fastas;

my  @HITS;
#my  %PDBHASH;

my  $QUERY;
my  $SBJCT;
my  $Queryscop;
my  %Queryfa;

## -------------------------------------------------------------------
## functional interface
##

sub ReadPDBFiles;
sub ReadFastaFiles;
sub WriteAllHits;
sub WriteHit;
sub RunDaliLite;
sub ReadAlignment;
sub ReadPDB;
sub WrapFasta;
sub RemoveUnecFiles;
sub RemoveDir;

## -------------------------------------------------------------------
##

ReadFiles( $PDBREPOS, \@Files );
ReadFiles( $FASTADIR, \@Fastas );

if( !GetScopClass( $INPUTPDB, \@Fastas, \$QUERY, \$Queryscop )) {
    die "ERROR: No scop info obtained for query pdb file."
}

ReadFastaFile("$FASTADIR/${QUERY}-${Queryscop}.fa", \%Queryfa );

## -------------------------------------------------------------------
##

my  $UPDIR = '../';
my  $UP2DIR = $UPDIR . $UPDIR;

my  $curdir = `pwd`;
chomp( $curdir );

$INPUTPDB = "$curdir/$INPUTPDB" unless( substr( $INPUTPDB, 0, 1 ) eq '/' );
$PDBREPOS = "$curdir/$PDBREPOS" unless( substr( $PDBREPOS, 0, 1 ) eq '/' );
$FASTADIR = "$curdir/$FASTADIR" unless( substr( $FASTADIR, 0, 1 ) eq '/' );
$outputname = "$curdir/$outputname" if( $outputname && substr( $outputname, 0, 1 ) ne '/' );

chdir $OUTPUTDIR;

## -------------------------------------------------------------------


-d $QUERY && die "ERROR: Directory $QUERY already exists.";
mkdir( $QUERY );
chdir( $QUERY );

$outputname = "$NAMEPAT\_$QUERY.$EXT" unless( $outputname );

if( $outputname ) {
    open( $$OUTPUT, ">$outputname" ) or die "ERROR: Cannot open file $outputname for writing.";
}

my  $c = 0;
foreach my $pf ( @Files ) {
    my  $sbjctscop;
    my  $queryaln;
    my  $sbjctaln;
    my (%sbjctfa );
    my ($querynwfas, $querynwpdb );
    my ($sbjctnwfas, $sbjctnwpdb );
    my ($qfabeg, $sfabeg, $tmpqs, $tmpss );

    my ($querystart, $sbjctstart );
    my ($queryend,   $sbjctend );

    my  $nwalign = nwalign->new();
    my  $zscore = 0.0;

    if( !GetScopClass( $pf, \@Fastas, \$SBJCT, \$sbjctscop )) {
        next
    }

    next if $QUERY eq $SBJCT;

    my  $subdirname = "$QUERY-$SBJCT";

    -d $subdirname && die "ERROR: Directory $subdirname already exists.";
    mkdir( $subdirname );
    chdir( $subdirname );

##    ReadPDB( "$INPUTPDB",     \@{$PDBHASH{$QUERY}} ) unless exists $PDBHASH{$QUERY};
##    ReadPDB( "$PDBREPOS/$pf", \@{$PDBHASH{$SBJCT}} ) unless exists $PDBHASH{$SBJCT};

    ReadFastaFile("$FASTADIR/${SBJCT}-${sbjctscop}.fa", \%sbjctfa );

    printf( STDERR "Processing $QUERY $SBJCT ...\n" );

    if( RunDaliLite( $INPUTPDB, "$PDBREPOS/$pf" )) {
        unless( -f $LOCKFILE ) {
            if( ReadAlignment( $ALIGNFILE, 
                               \$zscore, 
                               \$queryaln, \$sbjctaln, 
                               \$querystart, \$sbjctstart, \$queryend, \$sbjctend ))
            {
                if( 0 < length( $queryaln ) && 0 < length( $sbjctaln )) {
                    if( 0.0 <= $zscore ) 
                    {
                        chdir( $UPDIR );
                        RemoveUnecFiles( $subdirname );

                        ## Align query fasta sequence against query alignment (against subject) obtained by DaliLite
                        $nwalign->Query( $Queryfa{SEQN});
                        $nwalign->Sbjct( $queryaln );
                        $nwalign->Other( [$sbjctaln] );
                        $nwalign->Align( \$sbjctfa{QFALN}, \$sbjctfa{QSALN}, \@{$sbjctfa{QOALN}});

                        ## Align subject fasta sequence against subject alignment (against query) obtained by DaliLite
                        $nwalign->Query( $sbjctfa{SEQN});
                        $nwalign->Sbjct( ${$sbjctfa{QOALN}}[0] );
                        $nwalign->Other( [$sbjctfa{QFALN}, $sbjctfa{QSALN}]);
                        $nwalign->Align( \$sbjctfa{SFALN}, \$sbjctfa{SSALN}, \@{$sbjctfa{SOALN}});

                        $querynwfas = ${$sbjctfa{SOALN}}[0];
                        $querynwpdb = ${$sbjctfa{SOALN}}[1];
                        $sbjctnwfas = $sbjctfa{SFALN};
                        $sbjctnwpdb = $sbjctfa{SSALN};

                        undef $nwalign;

                        $qfabeg = $sfabeg = 0;
                        $qfabeg = $+[0] if $querynwpdb =~ /^\-+/;
                        $qfabeg = $+[0] if $sbjctnwpdb =~ /^\-+/ && $qfabeg < $+[0];
                        if( 0 < $qfabeg ) {
                            $tmpqs = substr( $querynwfas, 0, $qfabeg );
                            $tmpss = substr( $sbjctnwfas, 0, $qfabeg );
                            $querynwfas = substr( $querynwfas, $qfabeg );
                            $querynwpdb = substr( $querynwpdb, $qfabeg );
                            $sbjctnwfas = substr( $sbjctnwfas, $qfabeg );
                            $sbjctnwpdb = substr( $sbjctnwpdb, $qfabeg );
                            $qfabeg = 0;
                            $qfabeg++ while $tmpqs =~ /[^\-]/g;
                            $sfabeg++ while $tmpss =~ /[^\-]/g;
                        }
                        $qfabeg++;
                        $sfabeg++;

                        $tmpqs = 0;
                        $tmpqs = $-[0] if $querynwpdb =~ /\-+$/;
                        $tmpqs = $-[0] if $sbjctnwpdb =~ /\-+$/ && $-[0] < $tmpqs;
                        if( 0 < $tmpqs ) {
                            $querynwfas = substr( $querynwfas, 0, $tmpqs );
                            $querynwpdb = substr( $querynwpdb, 0, $tmpqs );
                            $sbjctnwfas = substr( $sbjctnwfas, 0, $tmpqs );
                            $sbjctnwpdb = substr( $sbjctnwpdb, 0, $tmpqs );
                        }

##printf( STDERR ">Qseq $qfabeg\n%s\n>Qstruct\n%s\n>Sstruct\n%s\n>Sseq $sfabeg\n%s\n\n",$querynwfas,$querynwpdb,$sbjctnwpdb,$sbjctnwfas);

                        ##NOTE: Obsolete
                        ##push @HITS, [ $zscore, $QUERY, $SBJCT,  $Queryscop, $sbjctscop,
                        ##              $querystart, $sbjctstart, $queryend,  $sbjctend,
                        ##              $queryaln,   $sbjctaln ];
                        push @HITS, [ $zscore, $QUERY, $SBJCT,  $Queryscop, $sbjctscop,
                                      $querystart, $sbjctstart, $queryend,  $sbjctend,
                                      $querynwpdb, $sbjctnwpdb, $querynwfas,$sbjctnwfas, $qfabeg, $sfabeg ];
                        ##==> last if 2 < $c++;## <===
                        next
                    }
                    else {
                        printf( STDERR "WARNING: Z-score not extracted. Skipping that.\n" );
                    }
                }
                else {
                    printf( STDERR "    (Dissimilar)\n" );
                }
            }
        }
        else {
            printf( STDERR "WARNING: Dali FAILED; lock file left.\n" );
            ##next
        }
    }
    chdir( $UPDIR );
    RemoveDir( $subdirname );
}

WriteAllHits( $OUTPUT, \@HITS );

close( $OUTPUT ) if $outputname;

chdir $curdir;

## -------------------------------------------------------------------
## write processed hit to file
##

sub WriteAllHits
{
    my  $reffile = shift;
    my  $refhits = shift;

    return unless $refhits;
    @{$refhits} = sort { $b->[0] <=> $a->[0] } @{$refhits};

    foreach my $rec ( @{$refhits} ) {
        WriteHit( $reffile, 
                  $rec->[1], $rec->[2],
                  $rec->[3], $rec->[4],
                  $rec->[5], $rec->[6], $rec->[7], $rec->[8],
                  $rec->[0],
                  \$rec->[9], \$rec->[10], \$rec->[11], \$rec->[12], $rec->[13], $rec->[14],
        );
    }
}

## -------------------------------------------------------------------
## write processed hit to file
##

sub WriteHit
{
    my  $reffile = shift;
    my  $llquery = shift;
    my  $llsbjct = shift;
    my  $llqueryscop = shift;
    my  $llsbjctscop = shift;
    my  $llquerystart = shift;
    my  $llsbjctstart = shift;
    my  $llqueryend = shift;
    my  $llsbjctend = shift;
    my  $llzscore = shift;
    my  $rqudalialn = shift;
    my  $rsbdalialn = shift;
    my  $rqufastaln = shift;
    my  $rsbfastaln = shift;
    my  $rqubeg = shift;
    my  $rsbbeg = shift;

    printf( $reffile ">%-9s %-13s %4d SEQUENCE\n",
                $llquery, $llqueryscop, $rqubeg );
    WrapFasta( $reffile, $rqufastaln );
    printf( $reffile ">%-9s %-13s %4d %-4d %10g\n",
                $llquery, $llqueryscop, $llquerystart, $llqueryend, $llzscore );
    WrapFasta( $reffile, $rqudalialn );
    printf( $reffile ">%-9s %-13s %4d %-4d\n",
                $llsbjct, $llsbjctscop, $llsbjctstart, $llsbjctend );
    WrapFasta( $reffile, $rsbdalialn );
    printf( $reffile ">%-9s %-13s %4d SEQUENCE\n",
                $llsbjct, $llsbjctscop, $rsbbeg );
    WrapFasta( $reffile, $rsbfastaln );
    printf( $reffile "//\n" );
}


## -------------------------------------------------------------------
## read in alignment generated by DaliLite structural superposition
##

sub ReadAlignment
{
    my  $alnfile = shift;
    my  $refzscore = shift;
    my  $refqueryaln = shift;
    my  $refsbjctaln = shift;
    my  $refquerystartpos = shift; ## query beginning position
    my  $refsbjctstartpos = shift; ## sbjct beginning position
    my  $refqueryendpos = shift; ## query end position
    my  $refsbjctendpos = shift; ## sbjct end position

    my  $queryseq;
    my  $sbjctseq;
    my  $line;
    my  $start = 0;
    my ($n, $qr, $sr );

    $$refzscore = -1.0;
    $$refqueryaln = '';
    $$refsbjctaln = '';
    $$refquerystartpos = 1;
    $$refsbjctstartpos = 1;
    $$refqueryendpos = 0;
    $$refsbjctendpos = 0;

    unless( open( ALN, "<$alnfile" )) {
         printf( STDERR "ERROR: Cannot open alignment file %s.", $alnfile );
         return 0;
    }

    while( <ALN> ) {
        chomp;
        next if /^$/;
        $line = $_;
        last if( $start && $line =~ /^$ALNSTOP/);##FIX

        if( $line =~ /\s+Sbjct=\S+\s+Z\-score=([\d\.eE\-\+]+)</ ) {
            $$refzscore = $1;
            $start = 1;
        }

        ##do { $start = 1; next } if $line =~ /^$ALNSTART/;
        next unless $start;

        if( $line =~ /^Query\s+(\S+)\s+\d+\s*$/ ) {
            $$refqueryaln = $$refqueryaln . $1;#uc( $1 );
            next
        }
        if( $line =~ /^Sbjct\s+(\S+)\s+\d+\s*$/ ) {
            $$refsbjctaln = $$refsbjctaln . $1;#uc( $1 );
            next
        }
    }

    close( ALN );

    $$refqueryaln =~ s/\./\-/g;
    $$refsbjctaln =~ s/\./\-/g;

##print "$$refqueryaln\n\n$$refsbjctaln\n\n\n\n";

    if( length( $$refqueryaln ) != length( $$refsbjctaln )) {
        printf( STDERR "WARNING: Alignment lengths of query and subject are not equal.\n" );
        return 0
    }

    ## remove preceding and trailing insertions
    ##
    for( $n = 0; $n < length( $$refqueryaln ); ) {
        $qr = substr( $$refqueryaln, $n, 1 );
        $sr = substr( $$refsbjctaln, $n, 1 );
        last if( $qr ne '-' && $sr ne '-' );
        substr( $$refqueryaln, $n, 1 ) = '';
        substr( $$refsbjctaln, $n, 1 ) = '';
        $$refquerystartpos++ if( $qr ne '-' );
        $$refsbjctstartpos++ if( $sr ne '-' );
    }
    for( $n = length( $$refqueryaln ) - 1; 0 <= $n && $n < length( $$refqueryaln ); $n-- ) {
        $qr = substr( $$refqueryaln, $n, 1 );
        $sr = substr( $$refsbjctaln, $n, 1 );
        last if( $qr ne '-' && $sr ne '-' );
        substr( $$refqueryaln, $n, 1 ) = '';
        substr( $$refsbjctaln, $n, 1 ) = '';
    }

    $queryseq = $$refqueryaln; $queryseq =~ s/\-//g;
    $sbjctseq = $$refsbjctaln; $sbjctseq =~ s/\-//g;

    $$refqueryendpos = $$refquerystartpos + length( $queryseq ) - 1;
    $$refsbjctendpos = $$refsbjctstartpos + length( $sbjctseq ) - 1;

    return 1;
}

## -------------------------------------------------------------------
## run DaliLite to superpose two structures: query and subject
##

sub RunDaliLite
{
    my  $queryfile = shift;
    my  $sbjctfile = shift;

    -f $queryfile || die "ERROR: Missing query file $sbjctfile.";
    -f $sbjctfile || die "ERROR: Missing subject file $sbjctfile.";


    system( "$DALIPROG -pairwise $queryfile $sbjctfile >$LOGFILE 2>$ERRFILE" );

    die "ERROR: Failed to execute Dali: $!" if $? == -1;
    if( $? & 127 ) {
        printf( STDERR "ERROR: Dali died with signal %d, %s coredump.\n",
                       ( $? & 127 ), ( $? & 128 )? 'with' : 'without' );
        exit 1
    }
    else {
        if( $? >> 8 != 0 ) {
            printf( STDERR "WARNING: Dali failed and exited with status %d; See above ( %s %s ).\n",
                $? >> 8, $QUERY, $SBJCT );
            return 0
        }
        else {
            printf( STDERR "Dali succeeded.\n" );
        }
    }

    return 1;
}


## -------------------------------------------------------------------
## extract scop classification from fasta filename
##

sub GetScopClass {
    my  $filename = shift;
    my  $refastas = shift;
    my  $refdomn = shift;
    my  $refscop = shift;

    my  $pfname = basename( $filename );
    $$refdomn = $pfname; $$refdomn =~ s/^(.+)\.[^\.]+$/$1/;
    my  @fafiles = grep { /^$$refdomn/ } @{$refastas};

    if( $#fafiles < 0 ) {
        printf( STDERR "WARNING: No fasta file found for $$refdomn.\n" ); 
        return 0 
    }
    if( 0 < $#fafiles ) {
        printf( STDERR "WARNING: >1 fasta files found for $$refdomn.\n" ); 
    }

    $$refscop = $fafiles[0];
    $$refscop =~ s/^[^\-]+\-(.+)\.[^\.]+$/$1/;
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

    closedir( DIR );
}

## -------------------------------------------------------------------
## read PDB file into list
##

sub ReadPDB {
    my  $filename = shift;
    my  $reflines = shift;

    open( PFD, "<$filename" ) or die "ERROR: Cannot open PDB file $filename";

    while( <PFD> ) {
        chomp;
        push @{$reflines}, $_;
    }

    close( PFD );
}

## -------------------------------------------------------------------
## read query description and sequence data
##

sub ReadFastaFile
{
    my  $filename = shift;
    my  $rqueryin = shift; ## ref. to hash
    
    unless( open( FA, $filename )) {
        printf( STDERR "ERROR: Failed to open $filename: $!\n" );
        return 0; 
    }
    while( <FA> ) {
        chomp;
        if( /^>(.+)\s*$/ ) {
            $$rqueryin{DESC} = $1;
            $$rqueryin{SEQN} = '';
            next;
        }
        s/\s//g;
        $$rqueryin{SEQN} .= uc( "$_" );
    }
    close( FA );
    return 1;
}



## -------------------------------------------------------------------
## wrap Fasta sequence and write it to file
##

sub WrapFasta {
    my  $reffile = shift;
    my  $reffasta = shift;
    my  $width = 80;

    for( my $n = 0; $n < length( $$reffasta ); $n += $width ) {
        printf( $reffile "%s\n", substr( $$reffasta, $n, $width ));
    }
}

## -------------------------------------------------------------------
## remove unnecessary files
##

sub RemoveUnecFiles {
    my  $dirname = shift;
    my  @required = ( $SUPOSFILE, $ALIGNFILE, $MATRIXFILE, $RANGESFILE );

    RemoveDir( $dirname, \@required );
}

## -------------------------------------------------------------------
## remove directory
##

sub RemoveDir {
    my  $dirname = shift;
    my  $refrequired = shift;
    my  @entries;
    my  @files;
    my  @dirs;
    my  $a;

    opendir( DIR, $dirname ) || die "ERROR: Directory $dirname does not exist.";
    @entries = readdir( DIR );
    closedir( DIR );

    @files = grep { -f "$dirname/$_" } @entries;
    @dirs = grep { !/^\./ && -d "$dirname/$_" } @entries;

    if( $refrequired ) {
        for( $a = 0; $a <= $#files; ) {
            if( grep {/$files[$a]/} @{$refrequired} ) {
                splice( @files, $a, 1 );
                next
            }
            $a++
	}
    }

    my  $locurdir = `pwd`;
    chomp( $locurdir );
    chdir( $dirname );
    unlink( @files );

    foreach my $dd ( @dirs ) {
        opendir( DR, $dd ) || die "ERROR: Directory $dirname/$dd does not exist.";
        @files = grep { -f "$dd/$_" } readdir( DR );
        closedir( DR );
        if( scalar( @files )) {
            chdir( $dd );
            unlink( @files );
            chdir( $UPDIR );
	}
        rmdir( $dd );
    }

    chdir( $locurdir );
    rmdir( $dirname ) unless $refrequired;
}

