#!/bin/bash

#Directory containing bed files to be analyzed. Include ending "/"
WORKINGDIRECTORY="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly"

#Name of input bed files to analyze
PICLUSTERFILE="Aag2_piClusters.bed"

#Name of initial output SORTED bedfiles
TEFILESORTED="Aag2_Contigs_TEs_sorted.bed" #This should already exist from running FindClosestTE.sh
EVEFILESORTED="Aag_Contigs_EVEs_NEW_cut_sorted.bed" #This should already exist from running FindClosestTE.sh
PICLUSTERFILESORTED="Aag2_piClusters_sorted.bed"

#Sort bedfiles (if necessary)
bedtools sort -i ${WORKINGDIRECTORY}/${PICLUSTERFILE}  > ${WORKINGDIRECTORY}/${PICLUSTERFILESORTED}

#Get TEs and EVEs which overlap piClusters only. Only want overlap, but bedtools doesn't allow to ignore upstream and downstream. Ignore only downstream for now, but will need to filter for overlap only in actual scripts (ie only keep those where distance is 0)
bedtools closest -a ${WORKINGDIRECTORY}/${TEFILESORTED} -b ${WORKINGDIRECTORY}/${PICLUSTERFILESORTED} -id -D ref > ${WORKINGDIRECTORY}/${TEFILESORTED}'_closestTopiClusters.txt'

bedtools closest -a ${WORKINGDIRECTORY}/${EVEFILESORTED} -b ${WORKINGDIRECTORY}/${PICLUSTERFILESORTED} -id -D ref > ${WORKINGDIRECTORY}/${EVEFILESORTED}'_closestTopiClusters.txt'
