#!/bin/bash
ports="1025 1028"
#ports="1025"
for P in ${ports}; do
  cmd="rsync -av -e \"ssh -p${P}\" --exclude '*~' --exclude '.*' --exclude 'build*' --exclude '*.obs' * localhost:/home/mindaugas/projects/cother/";
  echo "${cmd}"; eval ${cmd}
done
cmd="rsync -av --exclude '*~' --exclude '.*' --exclude 'build*' --exclude '*.obs' * /media/mindaugas/TSB1T/Linux/mindaugas/projects/cother/"
echo "${cmd}"; eval ${cmd}
#cmd="rsync -av --exclude '*~' --exclude '*.obs' * /media/mindaugas/ADATA\\ SD700/Linux/mindaugas/projects/cother/"
#echo "${cmd}"; eval ${cmd}

