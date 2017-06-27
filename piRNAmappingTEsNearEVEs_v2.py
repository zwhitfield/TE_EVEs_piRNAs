# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 12:41:09 2016

@author: zwhitfield
"""

import sys
import pandas as pd
import numpy as np
import scipy.stats as sp
#scipy.__version__
#np.__version__
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('ggplot')

outputdir = "/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/"

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#FUNCTIONS-----------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#only uses 5' position of piRNA (technically 3' position if piRNA is antisense to genome)
def map_piRNAs (teFile, piFile,distanceStatus, kimuraStatus, allOrNearEVEs, sampleOrNot, needToCount):
    totalpiStats = pd.DataFrame()
    previousContig = ""
    if sampleOrNot == 'sample':
        teFile = teFile.sample(n=10000) #sample 1000 of all TEs for comparison. Entire file would be too big. Can try later.
    for index, row in teFile.iterrows():
        currentContig = row["ContigTE"]  # contig current EVE is located
        print (currentContig)

        if allOrNearEVEs == 'nearEVEs':
            EVEname = row["EVEdescription"]
            EVElength = row["EVEend"] - row["EVEstart"] + 1
        else:
            EVEname = 'None'
            EVElength = 1.0

        TEname = row["TEdescription"]
        TElength = row["TEend"] - row["TEstart"] + 1
        TEclass = row["TEclass"]
        TEfamily = row["TEfamily"]
        if distanceStatus == 'yes':
            EVE_TEdistance = row["Distance"]
        else:
            EVE_TEdistance = 'None'

        if kimuraStatus == 'yes':
            kimuraScore = row["KimuraDivergence"]
        else:
            kimuraScore = 'None'

        if currentContig != previousContig:
            if needToCount == 'yes':
                piRNAstartPositions = piFile.loc[piFile.ContigpiRNA == currentContig, "piRNAstart"]  # Start positions of all piRNAs mapping to currentContig
                piRNAstart, piRNAindexLocation, piRNAcount = np.unique(piRNAstartPositions, return_index=True, return_counts=True)  # Each unique start position, counter, counts of piRNAs with that start positon

                piRNA_DF = pd.DataFrame({"Contig": currentContig,
                                         "piRNAstart": piRNAstart,
                                         "Counts": piRNAcount})
                piRNA_DF["Counts"] = piRNA_DF["Counts"].astype('float')
            else:
                piRNA_DF = piFile[piFile['ContigpiRNA'] == currentContig]
                piRNA_DF = piRNA_DF.groupby(['ContigpiRNA', 'piRNAstart'], as_index=False)['Counts'].sum()
                piRNA_DF = pd.DataFrame(piRNA_DF)

        if allOrNearEVEs == 'nearEVEs':
            currentEVErange = range(row["EVEstart"], row["EVEend"] + 1)  # nucleotide positions of current EVE
            piEVEreads = piRNA_DF.loc[piRNA_DF['piRNAstart'].isin(currentEVErange).values, 'Counts']  # counts of each piRNA start positions in current EVE
        else:
            piEVEreads = [0.0,0.0]

        currentTErange = range(row["TEstart"], row["TEend"] + 1)  # nucleotide positions of current TE
        piTEreads = piRNA_DF.loc[piRNA_DF['piRNAstart'].isin(currentTErange).values, 'Counts']  # counts of each piRNA start positions in current TE

        # piStats = pd.DataFrame({"EVEdescription": [EVEname],
        #                         "piRNAmappingToEVE": [sum(piEVEreads)],
        #                         "piRNAsPerEVEnuc": [sum(piEVEreads) / EVElength],
        #                         "TEdescription": [TEname],
        #                         "piRNAmappingToTE": [sum(piTEreads)],
        #                         "piRNAsPerTEnuc": [sum(piTEreads) / TElength],
        #                         "TE_EVEdistance": [EVE_TEdistance],
        #                         "KimuraDS": [kimuraScore]}
        #                        )
        piStats = pd.DataFrame({"EVEdescription": [EVEname],
                                "piRNAmappingToEVE": [sum(piEVEreads)],
                                "piRNAsPerEVEnuc": [sum(piEVEreads) / EVElength],
                                "TEdescription": [TEname],
                                "TEclass": [TEclass],
                                "TEfamily":[TEfamily],
                                "piRNAmappingToTE": [sum(piTEreads)],
                                "piRNAsPerTEnuc": [sum(piTEreads) / TElength],
                                "TE_EVEdistance": [EVE_TEdistance],
                                "KimuraDS": [kimuraScore]}
                               )
        totalpiStats = pd.concat([totalpiStats, piStats])

        previousContig = currentContig

    return totalpiStats

def LoadDataAllTEs(filePathToEVEs, filterBypiCluster):
    if filterBypiCluster == 'no':
        allTEs = pd.read_csv(filePathToEVEs,
                             names=["ContigTE", "TEstart", "TEend", "TEdescription", "TEscore", "TEstrand", "TEfamily"],
                             sep="\t")

    if filterBypiCluster == 'yes':
        #Use following if loading in Aag2 with piClusters
        allTEs = pd.read_csv(filePathToEVEs,
                             names=["ContigTE", "TEstart", "TEend", "TEdescription", "TEscore", "TEstrand", "TEfamily", "piClusterContig", "piClusterStart", "piClusterEnd", "piClusterInfo1", "piClusterStrand", "piClusterInfo2", "piClusterInfo3", "Distance"],
                             sep="\t")
        allTEs = allTEs[allTEs['Distance'] == 0]

    # First, rename all LTR/Gypsy as LTR/Ty3_gyspsy. THere are both in the file, so need a consensus
    # First, rename all LTR/Copia as LTR/Ty1_copia. THere are both in the file, so need a consensus
    # First, rename all LTR/Pao as LTR/Pao_Bel. THere are both in the file, so need a consensus
    allTEs.loc[allTEs.TEfamily == 'LTR/Gypsy', 'TEfamily'] = "LTR/Ty3_gypsy"
    allTEs.loc[allTEs.TEfamily == 'LTR/Copia', 'TEfamily'] = "LTR/Ty1_copia"
    allTEs.loc[allTEs.TEfamily == 'LTR/Pao', 'TEfamily'] = "LTR/Pao_Bel"

    allTEs["TEclass"] = allTEs["TEfamily"].str.split('/').str[0]
    allTEs["TEfamily"] = allTEs["TEfamily"].str.split('/').str[1]
    return allTEs

def readInpiRNAdata (pathToFile,needToCount):
    if needToCount == 'yes':
        # Read in all piRNAs which map to the Aag2 PacBio genome (this is the non-compact form) so need to count each piRNA within the script
        piRNAreads = pd.read_csv(pathToFile,
            sep=",",
            names=["piRNAseq", "piRNAstrand", "ContigpiRNA", "piRNAstart", "piRNAquerySeq"])

    if needToCount == 'no':
        # Read in UNIQUELY MAPPING piRNAs ONLY which map to the Aag2 PacBio genome (this is the compact form, so already counted each piRNA)
        # *****Note: genomeSeq in following file means the genome sequence the piRNA mapped to. THIS INCLUDES MISMATCHES, so this is the rev. complement of the piRNA in the antisense case. It is the exact sequence in the sense case.
        # piRNAstart is the 5' position IN THE REFERENCE the piRNA maps to. In the case of antisense, this is the position mapping to the 3' end of the piRNA
        piRNAreads = pd.read_csv(pathToFile,
            sep="\t",
            names=["piRNAfastaID", "piRNAstrand", "ContigpiRNA", "piRNAstart", "genomeSeq"])
        # Need to reformat some of the columns
        piRNAreads["Counts"] = piRNAreads["piRNAfastaID"].str.split('-').str[1]
        piRNAreads["piRNAfastaID"] = piRNAreads["piRNAfastaID"].str.split('-').str[0]
        piRNAreads['Counts'] = piRNAreads['Counts'].astype('float')
        piRNAreads = piRNAreads[(piRNAreads["genomeSeq"].str.len() >= 24) & (piRNAreads["genomeSeq"].str.len() <= 30)]

    return piRNAreads

def readInTEnearestEVEsdata (pathToFile,filterBypiCluster, pathTopiClusterFile):
    # Read in file of all EVEs and nearest TEs
    TEsOverlapOrNearEVEs = pd.read_csv(pathToFile,
        sep="\t",
        )
    if filterBypiCluster == 'yes':
        # If desired, filter by which EVEs are in piClusters
        # Read in EVEs which were run with bedtools closest function to check for overlap with piClusters
        EVEsWithpiClusterDistance = pd.read_csv(pathTopiClusterFile,
            names=["ContigEVE", "EVEstart", "EVEend", "EVEdescription", "EVEscore", "EVEstrand",
                   "EVEpercIdent", "piClusterContig", "piClusterStart", "piClusterEnd", "piClusterInfo1",
                   "piClusterStrand", "piClusterInfo2", "piClusterInfo3", "Distance"],

            sep="\t")
        EVEsWithpiClusterDistance = EVEsWithpiClusterDistance[EVEsWithpiClusterDistance['Distance'] == 0]
        EVEsInpiClusters = EVEsWithpiClusterDistance.EVEdescription.unique()

        TEsOverlapOrNearEVEs = TEsOverlapOrNearEVEs[(TEsOverlapOrNearEVEs['EVEdescription'].isin(EVEsInpiClusters))]

    return TEsOverlapOrNearEVEs

def piClustersWithEVEs(pathTopiCluster_EVEs):
    EVEsWithpiClusterDistance = pd.read_csv(pathTopiCluster_EVEs,
                                            names=["ContigEVE", "EVEstart", "EVEend", "EVEdescription", "EVEscore",
                                                   "EVEstrand",
                                                   "EVEpercIdent", "piClusterContig", "piClusterStart", "piClusterEnd",
                                                   "piClusterInfo1",
                                                   "piClusterStrand", "piClusterInfo2", "piClusterInfo3", "Distance"],

                                            sep="\t")
    EVEsWithpiClusterDistance = EVEsWithpiClusterDistance[EVEsWithpiClusterDistance['Distance'] == 0]

    piClustersContainingEVEs = EVEsWithpiClusterDistance['piClusterStart'].map(str) + '_' \
                         + EVEsWithpiClusterDistance['piClusterEnd'].map(str) + '_' \
                         + EVEsWithpiClusterDistance['piClusterContig']

    return piClustersContainingEVEs


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Read in piRNAs mapping to genome of interest------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#Read in all piRNAs which map to the Aag2 PacBio genome (this is the non-compact form) so need to count each piRNA within the script
# piRNAreads = readInpiRNAdata("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/piRNAmapping/Aag2_PacBio/I234_dsFluc-B_Aag2-PB_Frozen_v1.csv", 'yes')

#OR--------------------------------------------------------

#Read in UNIQUELY MAPPING piRNAs ONLY which map to the Aag2 PacBio genome (this is the compact form, so already counted each piRNA)
piRNAreads = readInpiRNAdata("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/piRNAmapping/Aag2_PacBio_uniqueMappers/dsFluc-B_S2_c_Aag2_m1v1.map", 'no')

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Read in TEs nearest to EVEs-----------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#Read in file of all EVEs and nearest TEs

#Read in TEs nearest EVEs filtered by being in a piCluster
TEsOverlapOrNearEVEs = readInTEnearestEVEsdata("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/TEsClosestToEVEs_nearestOnly_withEVEtaxonomy.txt",
                                    'yes',
                                    '/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/Aag_Contigs_EVEs_NEW_cut_sorted.bed_closestTopiClusters.txt')

#Read in TEs nearest EVEs which are NOT filtered by being in a piCluster
TEsOverlapOrNearEVEs_ALL = readInTEnearestEVEsdata("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/TEsClosestToEVEs_nearestOnly_withEVEtaxonomy.txt",
                                    'no',
                                    '/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/Aag_Contigs_EVEs_NEW_cut_sorted.bed_closestTopiClusters.txt')

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Read in all TEs and narrow down to those on contigs WITHOUT EVEs OR in piClusters without EVEs
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#Load in original TE sorted bed file
#allAag2TEs = LoadDataAllTEs ("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/Aag2_Contigs_TEs_sorted.bed", 'no')

#OR---------------

#Load in TEs which were run with bedtools closest function to check for overlap with piClusters
allAag2TEs = LoadDataAllTEs ("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/Aag2_Contigs_TEs_sorted.bed_closestTopiClusters.txt",'yes')
#----------------------------------------------------

#Get list of Contigs with EVEs on them
contigsWithEVEs = TEsOverlapOrNearEVEs_ALL.ContigEVE.unique()

#Filter allTE list based on whether TEs are on contigs with or without an EVE(s)
allAag2TEs_withoutEVEs = allAag2TEs[~(allAag2TEs['ContigTE'].isin(contigsWithEVEs))]

allAag2TEs_withEVEs = allAag2TEs[(allAag2TEs['ContigTE'].isin(contigsWithEVEs))]

# len(allAag2TEs)
# len(allAag2TEs_withEVEs)
# len(allAag2TEs_withoutEVEs)


#OR--------------------------------------
#Look at TEs in piClusters with and without EVEs (a bit less arbitrary than comparing TEs nearest EVEs vs those on contigs without any EVEs)

piClustersContainingEVEs = piClustersWithEVEs('/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/Aag_Contigs_EVEs_NEW_cut_sorted.bed_closestTopiClusters.txt')

allAag2TEs['comboPiClusterInfo'] = allAag2TEs['piClusterStart'].map(str) \
                                                              + '_' + allAag2TEs['piClusterEnd'].map(str) \
                                                              + '_' + allAag2TEs['piClusterContig']

allAag2TEs_inPiClustersWithoutEVEs = allAag2TEs[~(allAag2TEs['comboPiClusterInfo'].isin(piClustersContainingEVEs))]
allAag2TEs_inPiClustersWithEVEs = allAag2TEs[(allAag2TEs['comboPiClusterInfo'].isin(piClustersContainingEVEs))]

# len(allAag2TEs)
# len(allAag2TEs_inPiClustersWithEVEs)
# len(allAag2TEs_inPiClustersWitouthEVEs)




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Determine which piRNAs 'map' to EVEs or TEs using Patrick's much quicker method-
#Then, save the datasets to file-------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#piRNAs mapping to nearest-EVE TEs
piRNAsMappingToEVEsAndNearestTEs = map_piRNAs(TEsOverlapOrNearEVEs,piRNAreads,'yes','no','nearEVEs','noSample', 'no')

piRNAsMappingToEVEsAndNearestTEs.to_csv(outputdir + "totalpiStats_EVEsAndNearestTEs_inPiClusters_uniqueMappersOnly.csv",
                                      header = True,
                                      index = False,
                                      sep = "\t")

#piRNAs mapping to TEs on contigs without EVEs
piRNAsMappingToTEs_contigsWithoutEVEs = map_piRNAs(allAag2TEs_withoutEVEs,piRNAreads,'no','no','all','noSample', 'no')

piRNAsMappingToTEs_contigsWithoutEVEs.to_csv(outputdir + "totalpiStats_TEsOnContigsWithoutEVEs_inPiClusters_uniqueMappersOnly.csv",
                                      header = True,
                                      index = False,
                                      sep = "\t")

#piRNAs mapping to TEs in piClusters without EVEs
piRNAsMappingToTEs_inPiClustersWithoutEVEs = map_piRNAs(allAag2TEs_inPiClustersWithoutEVEs,piRNAreads,'no','no','all','noSample', 'no')

piRNAsMappingToTEs_inPiClustersWithoutEVEs.to_csv(outputdir + "totalpiStats_TEsInPiClustersWithoutEVEs_uniqueMappersOnly.csv",
                                      header = True,
                                      index = False,
                                      sep = "\t")

#piRNAs mapping to TEs in piClusters with EVEs
piRNAsMappingToTEs_inPiClustersWithEVEs = map_piRNAs(allAag2TEs_inPiClustersWithEVEs,piRNAreads,'no','no','all','noSample', 'no')

piRNAsMappingToTEs_inPiClustersWithEVEs.to_csv(outputdir + "totalpiStats_TEsInPiClustersWithEVEs_uniqueMappersOnly.csv",
                                      header = True,
                                      index = False,
                                      sep = "\t")

#piRNAs mapping to all TEs in piClusters
piRNAsMappingToTEs_contigsAll = map_piRNAs(allAag2TEs,piRNAreads,'no','no','all','noSample', 'no')

piRNAsMappingToTEs_contigsAll.to_csv(outputdir + "totalpiStats_TEsInPiClustersAll_uniqueMappersOnly.csv",
                                      header = True,
                                      index = False,
                                      sep = "\t")

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Load above datasets without having to remake them----------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

#Read in one of the above datasets if desired
# piRNAsMappingToEVEsAndNearestTEs = pd.read_csv(outputdir + "totalpiStats_EVEandTE_" + elementToCompare + "_inPiClusters_uniqueMappersOnly.csv",
#                       sep = "\t",
#                       )
#
# piRNAsMappingToTEs_contigsWithoutEVEs = pd.read_csv(outputdir + "totalpiStats_TEsonContigsWithoutEVEs_" + elementToCompare + "_inPiClusters_uniqueMappersOnly.csv",
#                       sep = "\t",
#                       )




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Visualize piRNA counts mapping to TEs nearest EVEs vs. on contigs without EVEs--
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#
# fig = plt.figure(figsize=(8.5, 11), facecolor='white')
# ax1 = fig.add_subplot(2, 1, 1)  # says use 1 row,1 column for plotting area, and insert current graph into position 1
# ax1.hist(piRNAsMappingToEVEsAndNearestTEs['piRNAsPerTEnuc'], bins = np.linspace(0,10,40), bottom=0.1)
# ax1.set_xlim((0.0, 10.0))
# ax1.set_yscale('log')
# ax1.set_xlabel('#piRNAs/TE length')
# ax1.set_ylabel('Count (log10)')
# ax1.set_title("TEs nearest EVEs")
#
# ax2 = fig.add_subplot(2, 1, 2)  # says use 1 row,1 column for plotting area, and insert current graph into position 1
# ax2.hist(piRNAsMappingToTEs_contigsWithoutEVEs['piRNAsPerTEnuc'], bins = np.linspace(0,10,40), bottom=0.1)
# ax2.set_xlim((0.0, 10.0))
# ax2.set_yscale('log')
# ax2.set_xlabel('#piRNAs/TE length')
# ax2.set_ylabel('Count (log10)')
# ax2.set_title("TEs on contigs without EVEs")

# a=sns.kdeplot(piRNAsMappingToEVEsAndNearestTEs_specificElement['piRNAsPerTEnuc'], bw=0.2,shade=True,label="TEs nearest EVEs")
# sns.kdeplot(piRNAsMappingToTEs_contigsWithoutEVEs_specificElement['piRNAsPerTEnuc'], bw=0.2,shade=True,label="TEs on contigs without EVEs")
# a.set(xlabel='#piRNAs/TE length', ylabel='KDE (density)')

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#Visualize piRNA counts mapping to EVEs vs. their nearest TE---------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# fig = plt.figure(figsize=(8.5, 11), facecolor='white')
# ax1 = fig.add_subplot(2, 1, 1)  # says use 1 row,1 column for plotting area, and insert current graph into position 1
# ax1.hist(piRNAsMappingToEVEsAndNearestTEs['piRNAsPerTEnuc'], bins = np.linspace(0,50,50), bottom=0.1)
# ax1.set_xlim((0.0, 50.0))
# ax1.set_yscale('log')
# ax1.set_xlabel('#piRNAs/TE nucleotie')
# ax1.set_ylabel('Count (log10)')
# ax1.set_title("TEs nearest EVEs")
#
# ax2 = fig.add_subplot(2, 1, 2)  # says use 1 row,1 column for plotting area, and insert current graph into position 1
# ax2.hist(piRNAsMappingToEVEsAndNearestTEs['piRNAsPerEVEnuc'], bins = np.linspace(0,50,50), bottom=0.1)
# ax2.set_xlim((0.0, 50.0))
# ax2.set_yscale('log')
# ax2.set_xlabel('#piRNAs/EVE nucleotide')
# ax2.set_ylabel('Count (log10)')
# ax2.set_title("EVEs")
#
# a=sns.kdeplot(piRNAsMappingToEVEsAndNearestTEs_specificElement['piRNAsPerTEnuc'], bw=0.2,shade=True,label="TEs nearest EVEs")
# sns.kdeplot(piRNAsMappingToEVEsAndNearestTEs_specificElement['piRNAsPerEVEnuc'], bw=0.2,shade=True,label="EVEs")
# a.set(xlabel='#piRNAs/TE length', ylabel='KDE (density)')
#
# fig = plt.figure(figsize=(8.5, 11), facecolor='white')
# ax1 = fig.add_subplot(1, 1, 1)  # says use 1 row,1 column for plotting area, and insert current graph into position 1
# ax1.plot(piRNAsMappingToEVEsAndNearestTEs_specificElement['piRNAsPerEVEnuc'],
#          piRNAsMappingToEVEsAndNearestTEs_specificElement['piRNAsPerTEnuc'],
#          'ro')
# ax1.set_xlim((0.0001, 100.0))
# ax1.set_ylim((0.0001, 100.0))
# ax1.set_yscale('log')
# ax1.set_xscale('log')
# ax1.set_xlabel('#piRNAs/EVE nucleotide')
# ax1.set_ylabel('#piRNAs/TE nucleotide')
# ax1.set_title("EVEs and Nearest TEs in piClusters Only")
