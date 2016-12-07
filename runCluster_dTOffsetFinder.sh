#!/bin/bash

outputDir=/home/brotter/nfsShared/results/dTSineFitting/

numEntries=540800
entriesPerCluster=$((numEntries/4))
entriesPerCore=$((numEntries/256))
clusterNum=0
if   [ `hostname` == "anitaI.phys.hawaii.edu" ];   then clusterNum=0;
elif [ `hostname` == "anitaII.phys.hawaii.edu" ];  then clusterNum=1;
elif [ `hostname` == "anitaIII.phys.hawaii.edu" ]; then clusterNum=2;
elif [ `hostname` == "anitaIV.phys.hawaii.edu" ];  then clusterNum=3;
else echo "You aren't on a server I recognize!"; fi

echo "You are on server "${clusterNum}

for core in `seq 0 64`; do
    startEv=$((entriesPerCluster*clusterNum + core*entriesPerCore))
    stopEv=$((entriesPerCluster*clusterNum + (core+1)*entriesPerCore))
    ./dTOffsetFinder ${outputDir}/dTOffsetFinder_stEv${startEv}.root ${startEv} ${stopEv} &> logs/${startEv}.log  &
done
