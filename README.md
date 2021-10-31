COTHER2, a protein remote homology search and threading tool

(C)2020-2021 Mindaugas Margelevicius,
Institute of Biotechnology, Vilnius University

# Description

   COTHER is a fast deep learning-based threading method for protein homology 
   search and alignment generation.
   COTHER employs the [COMER2](https://github.com/minmarg/comer2) 
   GPU-accelerated search engine but produces alignments by integrating 
   inter-residue distance map comparison with profile-profile comparison. 
   For distance map comparison, COTHER evaluates how closely predicted 
   inter-residue distances for a query match the distances observed in a 
   protein structure.
   Hence, COTHER integrates sensitive homology detection by COMER2 and 
   identification of common structural features predicted by a deep learning 
   framework. 
   This combination makes homology search more sensitive and accurate.

   COTHER is licensed under GNU General Public License version 3. Please find
   the LICENSE and COPYING files.

<!--
# Web service

   COTHER is available as a web service hosted on 
   [The COMER web server](https://bioinformatics.lt/comer), 
   where the intermediate steps to the service functionality, including 
   profile construction and profile-profile search, are transparent to the 
   user. 

   For local use, an up-to-date 
   [PDB70 COTHER profile database](#cother-profile-database-availability) is 
   available.
-->

# Available Platforms

   The COTHER source code should compile and run on Linux and macOS. Please 
   note, however, that COTHER was tested on and the binaries 
   are provided for the following platforms:

  *  Linux x64

   COTHER was also compiled with the `clang` version 6 compiler and is
   expected to install and run on macOS.

# Hardware requirements

  *  CUDA-enabled GPU(s) with compute capability 3.5 (released in 2012) or
     above
  *  2 GB of RAM or more

# Structure of the Package

   The main directories are described below:

  *  build -- an empty directory to contain built files.

  *  Linux\_installer -- this directory contains the necessary files to 
   install the prebuilt COTHER software on Linux.

  *  src -- the directory of the source files.

  *  [cother-scores-optimizer](https://github.com/minmarg/cother/tree/master/cother-scores-optimizer) -- 
   COTHER scores optimizer developed for finding optimal parameters for 
   distance distribution match scores (DDMS) employed by COTHER.
   COTHER is configured with the optimal DDMS scores found using this 
   optimizer.

# Installation of pre-compiled binaries

   On Linux, run the shell script and follow the instructions:

     Linux_installer/COTHER-installer1.sh

   NOTE: system requirements for the COTHER software installed on Linux are 
   NVIDIA driver version >=418.87 and CUDA version >=10.1.

# Installation from source code

   ## Installation on Linux and macOS

   ### Software requirements

   To successfully build and install the COTHER software from the source code
   on Linux or macOS, these tools are required to be installed:

  *  CMake version 3.8 or greater

  *  GNU Make version 3.81 or greater

  *  GNU GCC compiler version 4.8 or greater, or LLVM clang compiler
     version 6 or greater (or another C++ compiler that supports C++11)

  *  [the NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-downloads) version 10.0 or greater

   ### Installation

   Run the shell script:

     BUILD_and_INSTALL_unix.sh

# Main executables

   The COTHER software (`cother`) runs on CUDA-capable GPU devices. 
   An appropriate NVIDIA driver must have been installed. (It can also be 
   installed during the installation of the CUDA toolkit; see 
   "Installation from source code.")

   The software package contains these main programs in the bin directory in
   the installation path:

  *  `makepro` and `makepro.sh`, developed for making COMER2(!) profiles. 
   These programs can be found in the 
   [COMER2](https://github.com/minmarg/comer2) software package but are 
   added here for convenience.
   It is recommended to use makepro.sh for enriching profiles with 
   secondary structure predictions.
   `makepro` and `makepro.sh` make COMER2 profiles in text format.

  *  `makecov`, calculate weighted cross-covariance matrix between the 
   positions of multiple sequence alignment (MSA).
   This matrix and the COMER2 profile are used for predicting an 
   inter-residue distance map using the 
   [ROPIUS0](https://github.com/minmarg/ropius0) deep learning framework.

  *  `adddist`, this program is used to add inter-residue distances 
   predicted by [ROPIUS0](https://github.com/minmarg/ropius0), or another 
   tool, to a COMER2 profile and produce a COTHER profile.

  *  `batchadddist.py`, make COTHER profiles in bulk given distance files and 
   COMER2 profiles. This utility is useful when generating COTHER profiles
   for database construction.

  *  `makedb`, make a COTHER profile database to be searched.
   `makedb` makes output profile databases in text format. They are
   cross-platform portable.

  *  db2bin, convert a COTHER profile database to binary format. 
   Only this format (i.e., `makedb` + `db2bin`) is valid for searching 
   using the `cother` program.
   Please note that the output of `db2bin` is platform-dependent, and 
   `db2bin` should be invoked on every platform.

  *  `cother`, the main program for homology search/threading using one or 
   more GPUs.

# Getting Started

   Homology search using COTHER corresponds to searching a database of 
   COTHER profiles with one or more query COTHER profiles. 

   Assume that a query profile `myprofile.tpro` and a profile database 
   `mydb[.bin]` have been obtained. Then the simplest way to run 
   `cother` is to type:

     cother -i myprofile.tpro -d mydb -o my_output_directory

   where `my_output_directory` is an output directory to store output
   alignments files for each query profile present in the input file
   `myprofile.tpro`.

   `cother` allows for multiple queries in the input file. In that case,
   profiles should be stacked one on top of the other. It is also possible to 
   search profile(s) against multiple profile databases:

     cother -i myprofile.tpro -d mydb1,mydb2,mydb3 -o my_output_directory

   or perform an all-against-all comparison:

     cother -i mydb -d mydb -o my_output_directory

   `cother` search, as well as the process of making profiles, can be
   controlled with options read from the options file options.txt in the `var`
   directory in the installation path:

     cother -i myprofile.pro -d mydb -o my_output_directory -p options.txt

   The user can copy the options file `options.txt` to a preferred location
   and modify option values.

  ## Profile construction

   A standard way to construct a COTHER profile from an MSA `mymsa.afa` 
   (e.g., in 
   [aligned FASTA format](https://github.com/minmarg/comer2#input-multiple-alignment)) 
   and using 
   [ROPIUS0](https://github.com/minmarg/ropius0) for inter-residue distance 
   map prediction includes the following steps:

  *  Make a COMER2 profile:

    makepro.sh -i mymsa.afa -o myprofile.pro

  *  Predict distances using ROPIUS0 (e.g., Docker container), assuming 
   `path_to_ropius0` is its installation directory:

    path_to_ropius0/srvs/distopred.sh -i mymsa.afa -o resdir

  *  Make a COTHER profile based on the constructed COMER2 profile 
   `myprofile.pro` and predicted distances `mymsa__pred__nonavg.prb` 
   found in the ROPIUS0 output directory `resdir`:

    adddist -v -i resdir/mymsa__pred__nonavg.prb -j myprofile.pro -o myprofile.tpro --dst=3,5,7 --prb=0.05

  ## Distance map prediction

   The [ROPIUS0](https://github.com/minmarg/ropius0) deep learning 
   framework can be used to predict an inter-residue distance map for a query 
   protein, as demonstrated above.

   Distances can also be predicted by any other suitable tool. 
   The format that `adddist` recognizes is the following. 

```
#R1 R2       Dst  Prob
1 7         11.0 0.094
1 8         15.0 0.087
1 9         14.0 0.080
1 10        11.0 0.075
...
```

   The first two columns give residue (one-based) indices followed by 
   distance information (distance(s), optional probabilities). 
   The file is required to provide the ordered upper triangle of the 
   distance matrix.

   The distances saved as a file (e.g., mymsa\_\_pred.dst) can then be 
   combined with a COMER2 profile myprofile.pro to make a COTHER profile:

    adddist -v -i mymsa__pred.dst -j myprofile.pro -o myprofile.tpro --dst=3

# COTHER Profile database construction

   Profile-profile search using COTHER employs a comparison of a 2D 
   inter-residue distance map predicted for a query with the map 
   intrinsic to a protein structure. The comparison of a profile with a 
   structure is more sensitive than comparing only profiles but requires the 
   structure of a target protein to be known.

   The COTHER profiles of a database are expected to represent protein 
   structures. 
   The database profiles can be constructed by first

  *  making distance files in the format specified in 
   [Distance map prediction](#distance-map-prediction) for each protein 
   structure using the structure and corresponding COMER2 profile as 
   input to the programs `infer/promage4cother_519.py` and `infer/msk2dst.py` 
   from the [ROPIUS0](https://github.com/minmarg/ropius0) software package

  *  and combining the inter-residue distances with the corresponding COMER2
   profiles using `batchadddist.py` or `adddist` as specified in 
   [Distance map prediction](#distance-map-prediction).

# COTHER Profile database availability

   For convenience, an up-to-date COTHER profile database for PDB70 is 
   available for download at:

<!--   [https://sourceforge.net/projects/cother/files/cother-profile-databases](https://sourceforge.net/projects/cother/files/cother-profile-database) -->

# Final Notes

   All executables in the COTHER software package invoked with the "-h"
   option print a list of valid command-line options.

<!-- # References -->

# Funding

The work was supported by the European Regional Development Fund 
[grant number 01.2.2-LMT-K-718-01-0028]

---

Contact: <mindaugas.margelevicius@bti.vu.lt>
