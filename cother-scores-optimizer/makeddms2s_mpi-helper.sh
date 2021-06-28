#!/bin/bash

MYDIR=$(dirname $0)

MAKEDDMS2S="/home/mindaugas/projects/cother/cother-scores-optimizer/makeddms2s++_mpi.pl"

SLURM_PARTITION=P16
#use half of available cores, physical cores, and instruct slurm to allocate 2 cores per task (768)
NTASKS=384

DATDIR="/share/write/mindaugas/cother-scop2.03-20-uniref50-folds-train-hhb"
RESDIR="/share/write/mindaugas/cother-scop2.03-20-uniref50-folds-train-hhb"

OUTDIR="${RESDIR}/makeddms2s.2.02.01.scaled.SS3.hdps.cvs2s"
STRALNDIR="${DATDIR}/top_hits_dali_scop_2.03_20_folds_train"
PREDPRODIR="${DATDIR}/tpro_pred005.2.02.01.scaled.SS3.hdps.cvs2s"
PRSTPRODIR="${DATDIR}/tpro_pdbdst.2.02.01.scaled.SS3.hdps.cvs2s"

THETAbeg=0.2; THETAend=0.2; THETAstep=0.05
ABSDIFFEXPbeg=0.2; ABSDIFFEXPend=0.2; ABSDIFFEXPstep=0.1
AVGDISTEXPbeg=1.8; AVGDISTEXPend=1.8; AVGDISTEXPstep=0.1

thetavals=($(perl -e "for(\$t=${THETAbeg};\$t<=${THETAend};\$t+=${THETAstep}){print \" \$t\"}"))
adexpvals=($(perl -e "for(\$t=${ABSDIFFEXPbeg};\$t<=${ABSDIFFEXPend};\$t+=${ABSDIFFEXPstep}){print \" \$t\"}"))
avexpvals=($(perl -e "for(\$t=${AVGDISTEXPbeg};\$t<=${AVGDISTEXPend};\$t+=${AVGDISTEXPstep}){print \" \$t\"}"))

datestr="$(date +%F_%X)"

theta=${THETAbeg}
adexp=${ABSDIFFEXPbeg}
avexp=${AVGDISTEXPbeg}

sbatch -p ${SLURM_PARTITION} -n ${NTASKS} -c 2 -t 4320 -J ddms2s \
    <<EOF
#!/bin/bash
echo "Job started: \$(date); Running on: \$(hostname)"
mpirun perl "${MAKEDDMS2S}" --out "${OUTDIR}/makeddms2s_${datestr}.out" --aln "${STRALNDIR}" --pro "${PREDPRODIR}" --prs "${PRSTPRODIR}" --theta ${theta} --adexp ${adexp} --avexp ${avexp}
EOF

