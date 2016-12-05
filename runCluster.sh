#!/bin/bash

numEntries=540800
entriesPerCluster=$((numEntries/4))
entriesPerCore=$((numEntries/256))
clusterNum=0
if   [ `hostname` == "anitai.phys.hawaii.edu" ];   then clusterNum=0;
elif [ `hostname` == "anitaii.phys.hawaii.edu" ];  then clusterNum=1;
elif [ `hostname` == "anitaiii.phys.hawaii.edu" ]; then clusterNum=2;
elif [ `hostname` == "anitaiv.phys.hawaii.edu" ];  then clusterNum=3;
else echo "You aren't on a server I recognize!"; fi

for core in `seq 0 64`; do
    startEv=$((entriesPerCluster*clusterNum + core*entriesPerCore))
    stopEv=$((entriesPerCluster*clusterNum + (core+1)*entriesPerCore))
    ./dTOffsetFinder rootFiles/dTOffsetFinder_stEv${startEv}.root ${startEv} ${stopEv}  >> logs/${i}.log &
done
