##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##------------------------------------------FUNCTIONS-------------------------------------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------

getData = function(pathToFile)
{
  piRNAmappings = read.table(pathToFile,
                             sep = "\t",
                             header = FALSE,
                             col.names = c("piRNAfastaID","piRNAstrand","Contig","piRNAstart","genomeSeq"),
                             stringsAsFactors = FALSE)
  
  piRNAmappings = piRNAmappings[nchar(piRNAmappings$genomeSeq) >= 24 & nchar(piRNAmappings$genomeSeq) <= 30,] #Filter by size
  piRNAmappings$PhasedPair = "NONE" #Will eventually label sense-antisense piRNAs whose 5' ends are offset by 10 nucleotides
  piRNAmappings$Counts = as.numeric(unlist(lapply(strsplit(piRNAmappings$piRNAfastaID,split = "-"),'[[',2))) #Only take second part of piRNAfastaID collumn (everything afte the '-')
  piRNAmappings$piRNAfastaIDonly = as.numeric(unlist(lapply(strsplit(piRNAmappings$piRNAfastaID,split = "-"),'[[',1))) #Only take first part of piRNAfastaID collumn (everything afte the '-')
  
  
  return(piRNAmappings)
}

##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
##---------------------------(Visualization of) piRNAs mapping to CFAV--------------------
##----------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------
outputdir = "/home/zwhitfield/Desktop/ForMarkGenomePaper/FrozenData/Aag2_assembly/piRNAdata/piRNAmapping/"
genomeSize =15000

# xlimBoundries = c(9750,10250) #Genome positions to graph
xlimBoundries = c(1,genomeSize) #Genome positions to graph
ylimBoundries = c(-2000,2000) #Min and max y axis values (antisense reads are counted in negative direction)

# sampleNames = c("CFAV_PITpaper",
#                 "SINV",
#                 "DENV",
#                 "PCLV_GP",
#                 "PCLV_RDRP",
#                 "PCLV_NC",
#                 "Kamiti",
#                 "Xincheng",
#                 "CoastalPV")
sampleNames = c("Kamiti")

regionsOfInterest = list(c(400,600))

# conditions = c("dsFluc-B_S2","dsAgo3-B_S6","dsPiwi4-B_S8")
conditions = c("dsFluc-B_S2", "I251_dsFluc-B")

# pdf(paste(outputdir,"piRNAmapTo", "all_1to15000_withSelectTy3.pdf",sep=""),height = 24, width = 11)

#*****Note: genomeSeq in following function means the genome sequence the piRNA mapped to. THIS INCLUDES MISMATCHES, so this is the rev. complement of the piRNA in the antisense case. It is the exact sequence in the sense case.
#piRNAstart is the 5' position IN THE REFERENCE the piRNA maps to. In the case of antisense, this is the position mapping to the 3' end of the piRNA
#piRNAfastaID is a unique ID asigned within the FASTA file (ie before mapping). This allows for comparison of piRNAs found in two different mapping/bowtie outputs

senseCounts = c()
antiSenseCounts = c()

for (currentSample in sampleNames)
  {
  print (currentSample)
  
  # pdf(paste(outputdir, currentSample, "/piRNAmapTo", currentSample, "_", xlimBoundries[1], "to", xlimBoundries[2],"_SINV_DENVdatasetComparison.pdf",sep=""),height = 8, width = 11)
  # pdf(paste(outputdir, currentSample, "/piRNAmapTo", currentSample, "_", regionsOfInterest[[match(currentSample,sampleNames)]][1], "to", regionsOfInterest[[match(currentSample,sampleNames)]][2],".pdf",sep=""),height = 8, width = 11)
  
  # par(mfrow=c(length(conditions),1)) #create one row for each graph to be made in 'conditions'
  par(mar=c(8,8,6,2))
  
  for (currentCondition in conditions)
    {
    print (currentCondition)
    piRNAsMappingToEVEs = getData(paste(outputdir, "Aag_EVEs","/",currentCondition,"_c_", "Aag_EVEs", "_v1.map",sep = ""))
    piRNAsMappingToEVEs_IDsOnly = piRNAsMappingToEVEs$piRNAfastaIDonly
    
    piRNAsMappingToEVEs_specificEVE = piRNAsMappingToEVEs[piRNAsMappingToEVEs$Contig == "gi|752455949|gb|AJG39332.1|_nucleopasid_protein_[Wutai_Mosquito_Virus]_000212F_396963_396283",]
    
    piRNAmappingToVirus1Genome = getData(paste(outputdir, currentSample,"/",currentCondition,"_c_", currentSample, "_v3.map",sep = ""))
    
    piRNAmappingToVirus1Genome$inEVEcolor = "red"
    piRNAmappingToVirus1Genome[piRNAmappingToVirus1Genome$piRNAfastaIDonly %in% piRNAsMappingToEVEs_IDsOnly, "inEVEcolor"] = "blue"
    
    
    piRNAmappingToVirus1Genome_SENSE = piRNAmappingToVirus1Genome[piRNAmappingToVirus1Genome$piRNAstrand == "+",]
    piRNAmappingToVirus1Genome_ANTISENSE = piRNAmappingToVirus1Genome[piRNAmappingToVirus1Genome$piRNAstrand == "-",]
    piRNAmappingToVirus1Genome_ANTISENSE$fivePrimeEnd = piRNAmappingToVirus1Genome_ANTISENSE$piRNAstart + nchar(piRNAmappingToVirus1Genome_ANTISENSE$genomeSeq) -1 #Get antisense piRNAs 5' positions in genome. 3'most positon + length of piRNA -1 = 5'position
    
    senseCounts = append(senseCounts, sum(piRNAmappingToVirus1Genome_SENSE$Counts))
    antiSenseCounts = append(antiSenseCounts, sum(piRNAmappingToVirus1Genome_ANTISENSE$Counts))
    
    
    coverageVec_SENSE = numeric(genomeSize) #will hold coverage of sense piRNAs 5' ends at each position in genome. This create vector of 0's that is length of genomeSize
    coverageVec_ANTISENSE = numeric(genomeSize) #will hold coverage of antisense piRNAs 5' ends at each position in genome. This create vector of 0's that is length of genomeSize

    coverageVec_SENSE_mapToEVEs = numeric(genomeSize) #will hold coverage of sense piRNAs 5' ends at each position in genome. This create vector of 0's that is length of genomeSize
    coverageVec_ANTISENSE_mapToEVEs = numeric(genomeSize) #will hold coverage of antisense piRNAs 5' ends at each position in genome. This create vector of 0's that is length of genomeSize
    
    for (currentLine in 1:length(piRNAmappingToVirus1Genome_SENSE$Counts)) #Loop through each sense piRNA
      {
      if (length(piRNAmappingToVirus1Genome_SENSE$piRNAstrand)>0 && piRNAmappingToVirus1Genome_SENSE$piRNAstrand[currentLine] == "+")
        {
        piRNA5primePosition = piRNAmappingToVirus1Genome_SENSE$piRNAstart[currentLine] #get 5' positon of current piRNA
        currentCount = piRNAmappingToVirus1Genome_SENSE$Counts[currentLine] #get count of current piRNA
        
        coverageVec_SENSE[piRNA5primePosition] = coverageVec_SENSE[piRNA5primePosition] + currentCount #add count of current piRNA to appropriate coverageVec position.
        
        if (piRNAmappingToVirus1Genome_SENSE$piRNAfastaIDonly[currentLine] %in% piRNAsMappingToEVEs_IDsOnly)
          {
          coverageVec_SENSE_mapToEVEs[piRNA5primePosition] = coverageVec_SENSE_mapToEVEs[piRNA5primePosition] + currentCount #add count of current piRNA to appropriate coverageVec position.
          }
        
        }
      }
    
    for (currentLine in 1:length(piRNAmappingToVirus1Genome_ANTISENSE$Counts)) #Loop through each antisense piRNA
      {
      if (length(piRNAmappingToVirus1Genome_ANTISENSE$piRNAstrand)>0 && piRNAmappingToVirus1Genome_ANTISENSE$piRNAstrand[currentLine] == "-")
        {
        piRNA5primePosition = piRNAmappingToVirus1Genome_ANTISENSE$fivePrimeEnd[currentLine] #get 5' positon of current piRNA
        currentCount = piRNAmappingToVirus1Genome_ANTISENSE$Counts[currentLine] #get count of current piRNA
        
        coverageVec_ANTISENSE[piRNA5primePosition] = coverageVec_ANTISENSE[piRNA5primePosition] - currentCount #subtract count of current piRNA to appropriate coverageVec position.
        
        if (piRNAmappingToVirus1Genome_ANTISENSE$piRNAfastaIDonly[currentLine] %in% piRNAsMappingToEVEs_IDsOnly)
          {
          coverageVec_ANTISENSE_mapToEVEs[piRNA5primePosition] = coverageVec_ANTISENSE_mapToEVEs[piRNA5primePosition] - currentCount #add count of current piRNA to appropriate coverageVec position.
          }
        
        
        }
      }
    

    plot(NA,NA, 
         xlim = c(xlimBoundries[1],xlimBoundries[2]),
         # xlim = c(regionsOfInterest[[match(currentSample,sampleNames)]][1],regionsOfInterest[[match(currentSample,sampleNames)]][2]),
         ylim = c(min(coverageVec_ANTISENSE),max(coverageVec_SENSE)),
         # ylim = c(ylimBoundries[1],ylimBoundries[2]),
         xlab = "",
         ylab = "",
         main = "",
         type = "l",
         pch = 21, col="black", bg = "black", cex = 0.5, 
         xaxt = "n",
         yaxt = "n")
    
    title(main= paste(currentCondition, " piRNAs mapping to ", currentSample, sep = ""), line = 2, cex.main=2)
    title(xlab= paste(currentSample," genome position", sep = ""), line = 6, cex.lab=2)
    title(ylab="piRNA counts (by 5' position)", line = 6, cex.lab=2)
    xInterval = seq(from = xlimBoundries[1], to = xlimBoundries[2], by = 1000) #For tick marks
    axis(side = 1, at = xInterval, cex.axis = 2, las = 2, labels = FALSE)#las=2 makes labels perpendicular to axis
    xInterval = seq(from = xlimBoundries[1], to = xlimBoundries[2], by = 1000)#For tick mark labels
    axis(side = 1, at = xInterval, cex.axis = 2, las = 2)#las=2 makes labels perpendicular to axis
    axis(side = 2, cex.axis = 2, las = 2)
    
    lines(1:genomeSize,coverageVec_SENSE, col = "red") #Plot coverage of sense piRNAs
    lines(1:genomeSize,coverageVec_ANTISENSE, col = "red") #Plot coverage of anti-sense piRNAs

    lines(1:genomeSize,coverageVec_SENSE_mapToEVEs, col = "blue") #Plot coverage of sense piRNAs
    lines(1:genomeSize,coverageVec_ANTISENSE_mapToEVEs, col = "blue") #Plot coverage of anti-sense piRNAs
    
    legend("topright", legend = c("piRNAs NOT mapping to EVEs", "piRNAs mapping to EVEs"), fill = c("red","blue"))
    
    #Identify sense/anti-sense piRNA 'pairs' that have 10bp offset at 5' ends
    coverageCutoff = 250 #Use to filter piRNA 'pair' calls
    counter = 1
    for (currentPos in 1:length(coverageVec_ANTISENSE))
      {
      if (coverageVec_SENSE[currentPos]>coverageCutoff && -1*(coverageVec_ANTISENSE[currentPos+9])>coverageCutoff) #Filter using coverage of 5' ends. Looking for 10bp-offset so add 9 to sense 5' position to check antisense coverage
        {
        #Add rectangle highlight each pair
        rect(xleft = currentPos, xright = currentPos+9,
             ytop = coverageVec_SENSE[currentPos], ybottom = coverageVec_ANTISENSE[currentPos+9],
             border = NA,
             col= rgb(0,0,1.0,alpha=0.1))
        
        piRNAmappingToVirus1Genome[piRNAmappingToVirus1Genome$piRNAstart == currentPos & piRNAmappingToVirus1Genome$piRNAstrand == "+", "PhasedPair"] = paste("Pair",counter,sep="") #Assign sense piRNAs with the current 5'bp position to a sense/anti-sense 'pair'
        piRNAmappingToVirus1Genome[piRNAmappingToVirus1Genome$piRNAstart+nchar(piRNAmappingToVirus1Genome$genomeSeq)-1 == currentPos+9 & piRNAmappingToVirus1Genome$piRNAstrand == "-", "PhasedPair"] = paste("Pair",counter,sep="") #Assign antisense piRNAs with a 10bp offset of current 5'bp position to same sense/anti-sense 'pair'
        counter = counter + 1 #increase so can assign next hit to new 'pair' assignment
        }
      }
    }
  # dev.off()
}

warnings()

# 
# plot(antiSenseCounts,senseCounts, 
#      xlim = c(1,max(antiSenseCounts,senseCounts)), 
#      ylim = c(1,max(antiSenseCounts,senseCounts)),
#      log = "xy")
# abline(a = 0, b = 1)
# text(x = antiSenseCounts, y = senseCounts, labels = sampleNames)
# 
# plot(antiSenseCounts,senseCounts, 
#      xlim = c(1,max(antiSenseCounts,senseCounts)), 
#      ylim = c(1,max(antiSenseCounts,senseCounts)))
# abline(a = 0, b = 1)
# text(x = antiSenseCounts, y = senseCounts, labels = sampleNames)
# 
# ratio = senseCounts/antiSenseCounts
# hist(ratio, breaks = 20)
