library(ggplot2)
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##------------------------------------------FUNCTIONS-------------------------------------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------

getData = function(pathToFile, EVEorTE)
{
  piRNAmappings = read.table(pathToFile,
                             sep = "\t",
                             header = TRUE,
                             stringsAsFactors = FALSE)
  if (EVEorTE == 'EVE')
  {
  piRNAmappings$EVEproteinID = unlist(lapply(strsplit(piRNAmappings$EVEdescription,split = "|",fixed = TRUE),'[[',4)) #Only take second part of piRNAfastaID collumn (everything afte the '-')
  }
  return(piRNAmappings)
}

##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##----Distribution of piRNAs mapping to EVEs of the same protein--------------------------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
{
piRNAmappings_TEsNearEVEs = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_EVEsAndNearestTEs_inPiClusters_uniqueMappersOnly.csv", 'EVE')

# g = ggplot(piRNAmappings_TEsNearEVEs, aes(x = EVEproteinID, y=piRNAsPerEVEnuc))
g = ggplot(piRNAmappings_TEsNearEVEs, aes(x = TEdescription, y=piRNAsPerTEnuc))

#Organize by EVEproteinID of each EVE
# g = g + geom_jitter(aes(color = EVEproteinID), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
# g = g + geom_jitter(aes(color = TEdescription), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
g = g + geom_point() +
  scale_fill_brewer(palette = 'Spectral') +
  scale_y_log10() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  # ggtitle("Individual TE Kimura Divergence Scores") + 
  xlab("EVE protein")+ 
  ylab("piRNAs per bp")
g
}
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##----Distribution of piRNAs mapping to TEs near EVEs vs. on contigs without EVEs---------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
{
  
#piRNAs mapping to TEs in piClusters near EVEs vs. those on contigs without EVEs
piRNAmappings_WithEVEs = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_EVEsAndNearestTEs_inPiClusters_uniqueMappersOnly.csv", 'EVE')
piRNAmappings_WithoutEVEs = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_TEsOnContigsWithoutEVEs_inPiClusters_uniqueMappersOnly.csv", 'TE')
#OR------
#piRNAs mapping to TEs in piClusters WITH EVEs vs. piClusters WITHOUT EVEs
piRNAmappings_WithEVEs = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_TEsInPiClustersWithEVEs_uniqueMappersOnly.csv", 'TE')
piRNAmappings_WithoutEVEs = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_TEsInPiClustersWithoutEVEs_uniqueMappersOnly.csv", 'TE')
#-----------------------------

#Filter if desired
piRNAmappings_WithEVEs = piRNAmappings_WithEVEs[piRNAmappings_WithEVEs$TEfamily == 'Pao_Bel',]
piRNAmappings_WithoutEVEs = piRNAmappings_WithoutEVEs[piRNAmappings_WithoutEVEs$TEfamily == 'Pao_Bel',]

#Combine datasets for visualization with ggplot
toCombine_TEsNearEVEs = data.frame(piRNAsPerTEnuc = piRNAmappings_WithEVEs$piRNAsPerTEnuc, sample = 'TEsNearEVEs')
toCombine_TEsNOTnearEVEs = data.frame(piRNAsPerTEnuc = piRNAmappings_WithoutEVEs$piRNAsPerTEnuc, sample = 'TEsNOTnearEVEs')

compareTEs = rbind(toCombine_TEsNearEVEs, toCombine_TEsNOTnearEVEs)

# g = ggplot(compareTEs, aes(x = piRNAsPerTEnuc, y = ..density..))
# g = g + geom_histogram(aes(fill = sample), position = 'dodge', bins = 10) +
g = ggplot(compareTEs, aes(x = piRNAsPerTEnuc))
g = g + geom_density(aes(fill = sample), alpha = .75, color = "black", bw = 0.1) +
  # scale_fill_brewer(palette = 'Spectral') +
  # scale_y_sqrt() +
  scale_x_continuous(limits = c(0, 0.6)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20,face="bold")) +
  ggtitle("piRNAs mapping to Pao Bel elements in piClusters") +
  xlab("piRNAs Per TE (normalized to TE length)") + 
  ylab("Density")
g

}
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##----Distribution of piRNAs mapping to TEs of different classes--------------------------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
{
#Load data of piRNAs mapping to all TEs in piClusters
piRNAmappings_TEsInPiClustersAll = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_TEsInPiClustersAll_uniqueMappersOnly.csv", 'TE')
#OR
#Load data of piRNAs mapping to TEs nearest EVEs in piClusters
piRNAmappings_TEsInPiClustersAll = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_EVEsAndNearestTEs_inPiClusters_uniqueMappersOnly.csv", 'TE')
#OR
#Load data of piRNAs mapping to TEs in piClusters on contigs without EVEs
piRNAmappings_TEsInPiClustersAll = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_TEsonContigsWithoutEVEs_inPiClusters_uniqueMappersOnly.csv", 'TE')
#OR
#Load data of piRNAs mapping to TEs in piClusters without EVEs
piRNAmappings_TEsInPiClustersAll = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_TEsInPiClustersWithoutEVEs_uniqueMappersOnly.csv", 'TE')
#OR
#Load data of piRNAs mapping to TEs in piClusters with EVEs
piRNAmappings_TEsInPiClustersAll = getData("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/totalpiStats_TEsInPiClustersWithEVEs_uniqueMappersOnly.csv", 'TE')

#Filter if desired
piRNAmappings_TEsInPiClustersAll = piRNAmappings_TEsInPiClustersAll[piRNAmappings_TEsInPiClustersAll$piRNAsPerTEnuc <=0.75,]
piRNAmappings_TEsInPiClustersAll = piRNAmappings_TEsInPiClustersAll[piRNAmappings_TEsInPiClustersAll$TEclass == 'Unknown',]
piRNAmappings_TEsInPiClustersAll = piRNAmappings_TEsInPiClustersAll[piRNAmappings_TEsInPiClustersAll$TEfamily == 'Ty3_gypsy',]

# mean(piRNAmappings_TEsInPiClustersTest$piRNAsPerTEnuc)

##----------------------------------------------------------------------------------------
#Plot histogram of distribution of number of piRNAs mapping to TEs(per TE nucleotide)
##----------------------------------------------------------------------------------------

g = ggplot(piRNAmappings_TEsInPiClustersAll, aes(x = piRNAsPerTEnuc))
g = g + geom_histogram(aes(fill = TEclass), color = "black",position = 'dodge', bins = 10) +
# g = g + geom_violin(aes(fill = TEclass)) +
# g = g + geom_density(aes(fill = TEclass), alpha = .75, color = "black") +
  # scale_fill_brewer(palette = 'Spectral') +
  scale_y_sqrt() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20,face="bold")) +
  # ggtitle("Individual TE Kimura Divergence Scores") + 
  xlab("piRNAsPerTEperNuc") + 
  ylab("Count")
g

##----------------------------------------------------------------------------------------
#Plot average number piRNAs per TE nuc vs TEclass (scale size of point to counts of each TE type)
##----------------------------------------------------------------------------------------

#Make a 'counts' column for each type of TE using table function
piRNAmappings_TEsInPiClustersAll$counts = table(piRNAmappings_TEsInPiClustersAll$TEclass)[piRNAmappings_TEsInPiClustersAll$TEclass]
piRNAmappings_TEsInPiClustersAll$counts = as.numeric(piRNAmappings_TEsInPiClustersAll$counts)
# sum(unique(piRNAmappings_TEsInPiClustersAll$counts))

# g = ggplot(piRNAmappings_TEsInPiClustersAll, aes(x = TEclass, y = piRNAsPerTEnuc)) +
g = ggplot(piRNAmappings_TEsInPiClustersAll, aes(x = reorder(TEclass, piRNAsPerTEnuc, function(x){ -mean(x) }),y = piRNAsPerTEnuc)) + #reorder code from:https://stackoverflow.com/questions/33613385/sort-bar-chart-by-sum-of-values-in-ggplot
  stat_summary(aes(size = counts), fun.y="mean", geom = "point") +
  # geom_point(aes(size = counts))
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20,face="bold")) +
  # ggtitle("Individual TE Kimura Divergence Scores") + 
  xlab("Transposon Class") + 
  ylab("Avg. number of piRNAs mapping TE (normalized by TE length)")
g
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#----Statistics of TE piRNA distributions-------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
{
  #Check for difference between distributions based on both Mann-Whitney U test or Kolmogorov-Smirnov test.
  #Compare Kimura scores of TEs closest to EVEs vs. Kimura scores of all TEs in genome.
  
  KSresult = ks.test(as.numeric(piRNAmappings_TEsNearEVEs$piRNAsPerTEnuc), as.numeric(piRNAmappings_TEsNOTnearEVEs$piRNAsPerTEnuc))
  MUresult = wilcox.test(as.numeric(piRNAmappings_TEsNearEVEs$piRNAsPerTEnuc), as.numeric(piRNAmappings_TEsNOTnearEVEs$piRNAsPerTEnuc))
  
  statResults = data.frame(pVal_KS = KSresult$p.value, 
                           testStatistic_KS = KSresult$statistic, 
                           pVal_MUW = MUresult$p.value, 
                           testStatistic_MUW = MUresult$statistic, 
                           row.names = 1)
  write.table(statResults, 
              file = paste("/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/comparisonOfpiRNADistributionsTEs/TEsNearEVEsVsTEsNotOnContigsWithEVEs_uniqueMappes_inPiClusters_Ty3_gypsyOnly.txt",sep=""),
              quote = FALSE,
              row.names = FALSE,
              sep = "\t")
}

