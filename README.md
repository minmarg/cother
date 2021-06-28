COTHER2, a cross-platform software package for
protein remote homology search and threading

(C)2020-2021 Mindaugas Margelevicius,
Institute of Biotechnology, Vilnius University

# Description

   The COTHER method based on sequence profile-profile comparison is one of
   the most sensitive and accurate computational tools developed for protein
   alignment and homology search. COTHER version 2 represents one
   of the fastest implementations of calculations for protein threading. 
   High COTHER performance is achieved by harnessing the
   power of the Graphics processing unit (GPU). 

   COTHER is licensed under GNU General Public License version 3. Please find
   the LICENSE and COPYING files included in the software package.

# Available Platforms

   The COTHER source code should compile and run on Linux, MS Windows, and 
   macOS. Please note, however, that COTHER was tested on and the binaries 
   are provided for the following platforms:

  *  MS Windows 64 (64 bit)
  *  Linux x64 (64 bit)

   COTHER was also compiled with the `clang` version 6 compiler and is
   expected to install and run on macOS.

# Hardware requirements

  *  CUDA-enabled GPU(s) with compute capability 3.5 (released in 2012) or
     above
  *  2 GB of RAM or more

# Getting the COTHER Software

   The package is available at:

   [https://github.com/minmarg/cother](https://github.com/minmarg/cother)

   The Docker image will be available at:

   https:\//hub.docker.com/r/minmar/cother

# Structure of the Package

   The main directories are described below:

  *  build -- an empty directory to contain built files.

  *  MS_Windows_installer -- contains an installer file for MS Windows.

  *  Linux_installer -- this directory contains the necessary files to 
     install the prebuilt COTHER software on Linux.

  *  src -- the directory of the source files.

# Installation of pre-compiled binaries

   On MS Windows, run the installer:

     MS_Windows_installer\COTHER-installer1.msi

   NOTE: system requirements for the COTHER software installed on Windows are 
   NVIDIA driver version >=425.25 and CUDA version >=10.1.

   On Linux, run the shell script and follow the instructions:

     Linux_installer/COTHER-installer1.sh

   NOTE: system requirements for the COTHER software installed on Linux are 
   NVIDIA driver version >=418.87 and CUDA version >=10.1.

# Installation from source code

   ## Installation on MS Windows

   ### Software requirements

   To successfully build and install the COTHER software from the source code
   on MS Windows, these tools are required to be installed:

  *  CMake version 3.8 or greater (free software)

  *  Visual C++ compiler, e.g., Visual Studio Community (free for open 
     source projects; COTHER is an open source project)

  *  [the NVIDIA CUDA toolkit](https://developer.nvidia.com/cuda-downloads) version 10.0 or greater 
     (free software)

   ### Installation

   Run the command (batch) file:

     BUILD_and_INSTALL_win64.cmd

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

# Getting Started

   The COTHER software runs on CUDA-capable GPU devices. An appropriate
   NVIDIA driver must have been installed. (It can also be installed during
   the installation of the CUDA toolkit; see "Installation from source 
   code.")

   The software package contains four main programs in the bin directory in
   the installation path:

  *  makepro and makepro.sh (makepro.cmd on MS Windows), developed for
     making COTHER profiles. It is recommended to use makepro.sh for enriching
     profiles with secondary structure predictions. makepro.sh, however,
     requires the external packages PSIPRED and PSI-BLAST to be installed.
     makepro and makepro.sh (makepro.cmd) make profiles in text format.
     Profiles, therefore, can be transferred between different platforms.

  *  makedb is developed for making a COTHER profile database to be searched.
     makedb makes output profile databases in text format. They are also 
     cross-platform portable.

  *  db2bin, developed for converting a COTHER profile database to binary
     format. For an n-fold read speedup, it is highly RECOMMENDED to 
     convert a profile database using db2bin before conducting homology
     search with the `cother` program. Please note that the output of db2bin
     is platform-dependent, and db2bin should be invoked on every platform.

  *  cother, the main program for homology search/threading using one or more GPUs.

   Assuming that a query profile myprofile.pro and a profile database mydb
   have been obtained, the simplest way to run `cother` is to type:

     cother -i myprofile.pro -d mydb -o my_output_directory

   where my_output_directory is an output directory to store output
   alignments files for each query profile present in the input file
   myprofile.pro.

   `cother` allows for multiple queries in the input file. In that case,
   profiles made using makepro or makepro.sh (makepro.cmd) should be stacked
   one on top of the other. It is also possible to search all profiles in
   one profile database against the profiles of another one:

     cother -i mydb1 -d mydb2 -o my_output_directory

   or perform an all-against-all comparison:

     cother -i mydb -d mydb -o my_output_directory

   Mutually aligning two profiles requires making a database of one of the
   two profiles:

     cother -i myprofile1.pro -d myprofile2_db -o my_output_directory

   `cother` search, as well as the process of making profiles, can be
   controlled with options read from the options file options.txt in the var
   directory in the installation path:

     cother -i myprofile.pro -d mydb -o my_output_directory -p options.txt

   The user can copy the options file options.txt to a preferred location
   and modify option values.

# Input Multiple Alignment

   The program makepro accepts input multiple alignment files in FASTA or
   STOCKHOLM format.

   The FASTA format can be described as follows. The section of each 
   sequence begins with a description line, whose first character is a ">"
   delimiter. Sequence data begins on the next line and can occupy multiple
   lines. An example of a multiple alignment in FASTA is provided below:

```
>d1qhka_ d.100.1.2 (A:) N-terminal domain of RNase HI...
GNFYAVRKGRE--T---G--------IYNTW---NECKNQVDGYG---GAIYKKFNSYEQAKSFLG
>gi|28379120|ref|NP_786012.1|:(2-47) ribonuclease H (putative)...
-KYYAVRKGRQ--P---G--------IYRTW---PETQKQVSGYP---QAQYKSFTSEKDAQDFMA
>gi|84386727|ref|ZP_00989753.1|:(2-47) hypothetical ribonuclease HI...
-KYYVVWKGRT--P---G--------IFTTW---NECKSQVDGFA---GARYKSFPTLGEAESAFG
>gi|116492108|ref|YP_803843.1|:(2-47) RNase H with double-stranded...
-KFYAVKKGRK--P---G--------LYLTW---DAAKQQVDGFA---GAVYKSFLTKAEAEEWMA
>gi|6323890|ref|NP_013961.1|:(1-47) Ribonuclease H1...
GNFYAVRKGRE--T---G--------IYNTW---NECKNQVDGYG---GAIYKKFNSYEQAKSFLG
```

   The package also contains the perl script blast2fa.pl to convert
   (PSI-)BLAST output to FASTA format. Please type `blast2fa.pl -h` for more
   information.

# Final Notes

   All executables in the COTHER software package invoked with the "-h"
   option print a list of valid command-line options.

# References


# Funding

The work was supported by the European Regional Development Fund 
[grant number 01.2.2-LMT-K-718-01-0028]

---

Contact: <mindaugas.margelevicius@bti.vu.lt>
