#!/bin/bash

MYDIR=$(dirname $0)

MAKEDDMS2S="${MYDIR}/makeddms2s++.pl"

DATDIR="/data/CASP14_datasets"

OUTDIR="${DATDIR}/pdb70_from_mmcif_200205__selection__makeddms2s__training"
STRALNDIR="${DATDIR}/pdb70_from_mmcif_200205__selection__Dali__training"
PREDPRODIR="${DATDIR}/pdb70_from_mmcif_200205__selection__cother_pro_pred01__training"
PRSTPRODIR="${DATDIR}/pdb70_from_mmcif_200205__selection__cother_pro_pdbdst__training"

NTHREADS=20

THETAbeg=0.2; THETAend=0.2; THETAstep=0.05
ABSDIFFEXPbeg=0.2; ABSDIFFEXPend=0.2; ABSDIFFEXPstep=0.1
AVGDISTEXPbeg=1.8; AVGDISTEXPend=1.8; AVGDISTEXPstep=0.1

thetavals=($(perl -e "for(\$t=${THETAbeg};\$t<=${THETAend};\$t+=${THETAstep}){print \" \$t\"}"))
adexpvals=($(perl -e "for(\$t=${ABSDIFFEXPbeg};\$t<=${ABSDIFFEXPend};\$t+=${ABSDIFFEXPstep}){print \" \$t\"}"))
avexpvals=($(perl -e "for(\$t=${AVGDISTEXPbeg};\$t<=${AVGDISTEXPend};\$t+=${AVGDISTEXPstep}){print \" \$t\"}"))

for theta in ${thetavals[@]}; do
  for adexp in ${adexpvals[@]}; do
    for avexp in ${avexpvals[@]}; do
      datestr="$(date +%F_%X)"
      "${MAKEDDMS2S}" --out "${OUTDIR}/makeddms2s_${datestr}.out" --aln "${STRALNDIR}" --pro "${PREDPRODIR}" --prs "${PRSTPRODIR}" --theta ${theta} --adexp ${adexp} --avexp ${avexp} --ts ${NTHREADS} >"${OUTDIR}/makeddms2s_${datestr}.log" 2>&1
    done
  done
done

echo Done.

