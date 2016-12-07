#!/bin/bash

outputDir=/home/brotter/nfsShared/results/findClockEdges/

numEntries=124
entriesPerCore=$((numEntries/64))
clusterNum=0
if  ! [ `hostname` == "anitaI.phys.hawaii.edu" ]; then echo "You aren't on a server I recognize!"; exit; fi

for core in `seq 0 64`; do
    startRun=$((core*2+1));
    stopRun=$(((core+1)*2+1));
    echo $startRun" "$stopRun
    root -b findClockEdges.C\(${startRun},${stopRun},\"${outputDir}/findClockEdges_${core}.root\"\) &>  ${outputDir}/logs/${core}.log &
done
