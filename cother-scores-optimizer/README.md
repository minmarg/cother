COTHER scores optimizer is developed for finding optimal parameters for 
distance distribution match scores (DDMS) employed by COTHER

(C)2021 Mindaugas Margelevicius,
Institute of Biotechnology, Vilnius University

# Description

   DDMS scores are obtained by dynamic programming (DP) using two distance
   distributions at different protein positions. It is the sum of several
   ungapped maximum-scoring segments identified during DP. The DDMS score
   depends on a number of parameters, including those of the score between
   two distance values from the respective distributions. These parameters
   are optimized with respect to the DDMS capability to discriminate 
   between structurally equivalent and non-equivalent positions. As DDMS 
   scores are always positive, the discrimination is expressed in positive
   and negative scores by a translation of the DDMS score to the 
   log-odds space. This information is printed on output.

# Building the optimizer

   Calculating DDMS scores is computationally intensive. Therefore, the 
   module for calculating DDMS scores is written in C++ and integrated 
   into the (Perl) optimizer `makeddms2s++.pl`.

   To build the optimizer, SWIG is assumed to be installed on the system.
   Make and C/C++ compiler toolkit are also required. The build steps 
   are as follows. First, an interface for the Perl program is generated 
   by the command:

```
swig -c++ -perl ddms.i
```

   Next, a makefile is generated using the following command:

```
perl Makefile.PL
```

   Finally, the target dynamic library required for the Perl optimizer 
   (`makeddms2s++.pl`) is obtained by typing

```
make
```

   For the optimizer to find the library in its directory, a link is 
   created, e,g. (linux), using

```
ln -s blib/arch/auto/ddms/ddms.so ddms.so
```

   Alternatively, and more conveniently, all these steps can be executed
   using a single command:

```
make -f ddms-optimizer.make
```

# Running the optimizer

   The available options of the makeddms2s++.pl program are listed below:

```
Usage:
makeddms2s++.pl <Parameters>

Parameters:

--out <filename>   Output filename.
--aln <directory>  Directory of structural alignments.
--pro <directory>  Directory of COTHER profiles with predicted distances.
--prs <directory>  Directory of COTHER profiles with distances from structure.
--eff <numbers>    Make score maps regarding these levels of effective 
                   number of observations.
           Default=1
--zsc <score>      Consider alignments of at least this Z-score.
           Default=3
--theta <parameter> Theta parameter for scores.
           Default=0.2
--adexp <parameter> Exponent for absolute difference (scores parameter).
           Default=0.2
--avexp <parameter> Exponent for average distance (scores parameter).
           Default=1.8
--pws              Derive scores from pairwise alignments.
--ovl              Produce overlapping score tables.
--pos              Calculate positional sequence weights.
--nosw             Do not calculate sequence weights.
--ts  <number>     Number of concurrent threads to use if in use
                   (Threads can be used: 1).
           Default=1
--printandexit     Print the table of scores between (pairwise) distance 
                   values and exit.
--help             This text.

```

   The user provides the output filename (`--out`), directory of 
   structural alignments (`--aln`) made by `dalirun.pl` (they are used 
   for the identification of structurally equivalent residues), directory
   of corresponding COTHER profiles (`--pro`) with predicted distances 
   (incorporated by `adddist` from the COTHER software package), and 
   directory of corresponding COTHER profiles (`--prs`) with the distances
   extracted from the structures (incorporated by `adddist` from the 
   COTHER software package).

   A subset of the DDMS score parameters can be specified too: `--theta`,
   `--adexp`, and `--avexp`. Other more advanced parameters should be
   chosen at compile time and can be modified in the upper part of the 
   ddms.h header file.

