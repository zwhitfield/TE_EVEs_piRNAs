#!/bin/bash

INPUTDIR="/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/piRNAmapping"
FILENAME="Ty3_Ele104_TF000385"
FOLDERNAME="Ty3_Ele104"

cd ${INPUTDIR}/${FOLDERNAME}

bowtie-build ${FILENAME}.fasta ${FILENAME}

#Map to piRNAs from Aag2 infected with SINV
bowtie --suppress 6,7,8 -f -v 3 ${INPUTDIR}/${FOLDERNAME}/${FILENAME} ${INPUTDIR}/dsFluc-B_S2_c.fa ${INPUTDIR}/${FOLDERNAME}/'dsFluc-B_S2_c_'${FOLDERNAME}'_v3.map'

#Map to piRNAs from Aag2 infected with DENV
#bowtie --suppress 6,7,8 -f -v 3 ${INPUTDIR}/${FOLDERNAME}/${FILENAME} ${INPUTDIR}/I251_dsFluc-B_c.fa ${INPUTDIR}/${FOLDERNAME}/'I251_dsFluc-B_c_'${FOLDERNAME}'_v3.map'
