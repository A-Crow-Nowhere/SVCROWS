#---
#title: "SVCROWS1.2.0"
#author: "Noah Brown"
#date: "2024-11-25"
#---

#       ███████
#   ██████  ████
# ██████████████
#      ██████████
#       ████████████████
#         ██████████████████████
#         ████████████████████████████
#           ██████████SV████████████████
#           ██████████████████████████████████
#             ████████████████████████████████████
#             ████████████████████████████████████████████
#               ████████████████████████████████████████
#                     ████████████████
#                       ██    ██
#                      ██    ██
#                     ██    ██
#               ██████ ███████





#Generate the featurelist in the internal running format.
CreateFL <- function(inputlist, xs, xl)
{
  #Column Names for the Consensus list
  print("Generating Feature List...")
  FL_columns <- c(
    #General
    "Chr",
    "Start",
    "End",
    "Length",
    "Type",
    "ID",
    "Var1",
    "Var2",
    "Var3",
    "IsKnown", #whole column is removed in post
    #For calculations
    "Size",
    "BPBoundSize",
    "ROPercentPass",

    #Record keeping for RO
    "ROCount",

    "BPStartLeft",
    "BPStartRight",
    "BPStartCount",

    "BPEndLeft",
    "BPEndRight",
    "BPEndCount",

    "NumReads",#whole column is removed in post
    "QScore",#whole column is removed in post

    "Frequency",
    "Rarity",
    "Matches"
  )
  #Make Tibble with correct column names (i.e. the consensus list)
  FLlist = data.frame(matrix(nrow = 0, ncol = length(FL_columns)))
  colnames(FLlist) = FL_columns

  for(i in seq_len(nrow(inputlist)))
  {
    empty_row <- as.data.frame(matrix(NA, ncol = ncol(FLlist), nrow = 1))
    colnames(empty_row) <- colnames(FLlist)
    FLlist <- rbind(FLlist, empty_row)

    FLlist$Chr[i] <- inputlist$Chr[i]
    FLlist$Start[i] <- inputlist$Start[i]
    FLlist$End[i] <- inputlist$End[i]
    FLlist$Length[i] <- (FLlist$End[i] - FLlist$Start[i])
    FLlist$Type[i] <- inputlist$Type[i]
    FLlist$ID[i] <- i
    FLlist$Var1[i] <- inputlist$Var1[i]
    FLlist$Var2[i] <- inputlist$Var2[i]
    FLlist$Var3[i] <- inputlist$Var3[i]
    FLlist$IsKnown[i] <- FALSE

    if (FLlist$Length[i] <= xs)
    {
      FLlist$Size[i] <- "Small"
      Sizes <- SizeDetermination(xs)
    }
    else if (FLlist$Length[i] <= xl)
    {
      FLlist$Size[i] <- "Medium"
      Sizes <- SizeDetermination(FLlist$Length[i])
    }
    else
    {
      FLlist$Size[i] <- "Large"
      Sizes <- SizeDetermination(xl)
    }

    FLlist$BPBoundSize[i] <- Sizes[1]
    FLlist$ROPercentPass[i] <- Sizes[2]

    FLlist$ROCount[i] <- 0

    FLlist$BPStartLeft[i] <- round(FLlist$Start[i] - (Sizes[1]/2))
    FLlist$BPStartRight[i] <- round(FLlist$Start[i] + (Sizes[1]/2))
    FLlist$BPStartCount[i] <- 0


    FLlist$BPEndLeft[i] <- round(FLlist$End[i] - (Sizes[1]/2))
    FLlist$BPEndRight[i] <- round(FLlist$End[i] + (Sizes[1]/2))
    FLlist$BPEndCount[i] <- 0

    FLlist$NumReads[i] <- 0
    FLlist$QScore[i] <- 0

    FLlist$Frequency[i] <- 0
    FLlist$Rarity[i] <- "Unassigned"
    FLlist$Matches[i] <- ""

  }

  print("...Done")
  return(FLlist)
}



#Main loop for Hunt mode
SVCrowsHunt <- function(QueryListX, FeatureListX, PerSampleX, BPfactor)
{

  #Type safe inputs for tables so nothing gets confused
  WorkingCon <<- FeatureListX
  QueryList <<- QueryListX
  PerSampleList <<- PerSampleX

  PerSampleList <- PerSampleList %>% mutate(MatchFeatureIDs = "")
  WorkingCon <- WorkingCon %>% mutate(IsKnown = FALSE)
  QueryList <- QueryList %>% mutate(IsKnown = FALSE)

  #Begin Query loop. Each query is assessed only once, Either a match is found
  #and is iterated, or it finds no matches and is added to the end of the
  #Consensus lsit
  for(i in seq_len(nrow(QueryList)))
  {
    QueryPosition <<- i
    print(paste("hQP: ",QueryPosition))

    #Set variables for convenience
    qChr <- QueryList$Chr[i]
    qStart <- QueryList$Start[i]
    qEnd <- QueryList$End[i]
    qLength <- QueryList$Length[i]
    qType <- QueryList$Type[i]
    qID <- QueryList$ID[i]
    qV1 <- QueryList$Var1[i]
    qV2 <- QueryList$Var2[i]
    qV3 <- QueryList$Var3[i]
    qIsKnown <- QueryList$IsKnown[i]
    qScore <- QueryList$QScore[i]
    qReads <- QueryList$NumReads[i]

    #Breakpoints (BP) can only iterate once. So if there are multiple
    #consensus list entries with a matching breakpoint, but different enough
    #sizes to be called unique, it will only iterate the largest of those
    #entries. When these get called as true, the query can no longer iterate BP.
    BPStartHasPassed <<- FALSE
    BPEndHasPassed <<- FALSE

    #Stops the Consensus List (CL) while loop. When a match is found, or a query
    #reaches the end of this CL, it will be added, and move on to the next
    #Query
    QueryAdded <<- FALSE

    #split the working consensus into entries that match the Chr

    twc <- WorkingCon
    MatchList <<- twc %>%
      filter(
        Chr == qChr,
        BPEndRight >= (qStart - y1LargeBound),
        BPStartLeft <= (qEnd + y1LargeBound)
      )

    # Filter conditions combined into a single filter function call for MinusMatchList
    MinusMatchList <<- twc %>%
      filter(
        Chr != qChr |
          BPEndRight < (qStart - y1LargeBound) |
          BPStartLeft > (qEnd + y1LargeBound)
      )
    suppressWarnings(remove(WorkingCon))
    for(j in seq_len(nrow(MatchList)))
    {
      #Set variables for convenience, all of these get passed to subsequent
      #functions. rX stands for 'reference'. Reference is the current CL entry
      #the query is being compared to.

      rChr <<- MatchList$Chr[j]
      rStart <- MatchList$Start[j]
      rEnd <- MatchList$End[j]
      rLength <- MatchList$Length[j]
      rType <- MatchList$Type[j]
      rID <- MatchList$ID[j]
      rV1 <- MatchList$Var1[j]
      rV2 <- MatchList$Var2[j]
      rV3 <- MatchList$Var3[j]
      rIsKnown <- MatchList$IsKnown[j]

      rReads <- MatchList$NumReads[j]
      rScore <- MatchList$QScore[j]

      rROPcntPass <- MatchList$ROPercentPass[j]

      rStartLeft <- MatchList$BPStartLeft[j]
      rStartRight <- MatchList$BPStartRight[j]
      rEndLeft <- MatchList$BPEndLeft[j]
      rEndRight <- MatchList$BPEndRight[j]

      #Reset the variables storing status of a query finding a match to FALSE
      #by default
      #In RO
      ROPass <- FALSE

      #In BP
      BPStartPass <- FALSE
      BPEndPass <- FALSE



      #Chomosome and Type must match to begin comparison
      #Breakpoint passes ('start' is the 5' BP, end is the 3' BP.)
      BPStartPass <- BreakPointStartCalc(qStart, rStartLeft, rStartRight)
      BPEndPass <- BreakPointEndCalc(qEnd, rEndLeft, rEndRight)

      #If query has same breakpoint as ref, drop overlap threshold to minimum
      #I.e. Breakpoint consideration
      TempROPercentPass <<- rROPcntPass
      if ((BPStartPass || BPEndPass) && BPfactor)
      {
        rROPcntPass = y2SmallOV
      }
      #calculate ROPass
      ROPass <- RecpricolOverlapCalc(qStart, rStart, qEnd, rEnd, rIsKnown, qIsKnown, qLength, rLength, rROPcntPass, TRUE, j)
      #Reset Passing threshold
      rROPcntPass <-  TempROPercentPass

      #print(ROPass)

      #If there is no match on anything, and we aren't at the end of the list,
      #just iterate the Consensus Position
      if (BPStartPass == FALSE && BPEndPass == FALSE && ROPass == FALSE)
      {

      }

      #If there is no match, and we are at the end of the list, add the new
      #entry
      else if (BPStartPass == FALSE && BPEndPass == FALSE && ROPass == FALSE && (j == nrow(MatchList)))
      {

      }
      else
      {
        #Based on Passing Criteria, enter the descion tree
        rowIn <- Decision(MatchList[j,], QueryList[QueryPosition,],
                          qIsKnown, rIsKnown, BPStartPass, BPEndPass,
                          BPStartHasPassed, BPEndHasPassed, ROPass,
                          PerSampleList[QueryPosition,],
                          qLength, rID, rReads, rScore, qReads, qScore, j, TRUE)

        #print(rowIn$Crow)

        #If a match is found in the query, update the CL and PSL
        if ((BPStartPass || BPEndPass || ROPass))
        {
          #Update rows based on results of decision
          tempWorkC <- MatchList
          tempWorkC[j,] <- rowIn$Crow
          MatchList <<- tempWorkC

          catmatch <- paste0(rowIn$Prow$MatchFeatureIDs[1], ":", rowIn$Crow$Var3[1])
          #print(catmatch)
          rowIn$Prow$MatchFeatureIDs[1] <- catmatch


          tempPerSamp <- PerSampleList
          suppressWarnings(remove(PerSampleList))
          tempPerSamp[QueryPosition,] <- rowIn$Prow
          PerSampleList <<- tempPerSamp

          remove(rowIn, tempWorkC, catmatch, tempPerSamp)
        }
      }
    }
    #print(MatchList)
    temp <- ReMergeLists(MatchList, MinusMatchList)
    WorkingCon <<- temp
  }

  FCLout1 <- FCLprocessing(WorkingCon)
  FCLout2 <- FCLout1[-c(10,21,22,24)]

  FPSLout1 <- PSLprocessing(PerSampleList, FCLout1)
  FPSLout2 <- FPSLout1[-c(10,15,16)]


  #print("3")
  outs <- list(FCL = FCLout2, FPSL = FPSLout2, AFPSL = FPSLout2, QL = QueryList)
  #print("4")
  return(outs)
}


#' Run SVCROWS in "Hunt" mode
#'
#' @description
#' This mode uses the framework of overlap comparison detailed for "scavange" mode, but applies it to overlapping features in the reference genome. For example, if you want to know what set of genes are captured in your Query List, you can input a 'Feature list, which contains information about all of the genes in the genome. Consider using for intergenic regions, Transcriptional start sites, enhancers, etc.
#'
#'
#' @param InputQueryList Query List in the designated format (see user manual)
#' @param FeatureList Similar to the Query List, with slight differences. Contains a list of all features the user wants to query the Query List against. (see user manual)
#' @param OutputDirectory Directory to write output files (see user manual)
#' @param BPfactor Boolean. By Default TRUE, uses breakpoints as a secondary piece of information to call RO overlapping regions. Regardless, matching regions will be counted in the FCL.
#' @param DefaultSizes Boolean. By default False, which uses the other 6 inputs given by the user. If set to True, SVCROWS will use the 2nd and 4th quartiles to generate values for the six variables (see user manual for details).
#' @param xs Int. Small SV size cutoff
#' @param xl Int. Large SV size cutoff
#' @param y1s Int. Small SV BP region size
#' @param y1l Int. Large SV BP region size
#' @param y2s Int. Small SV RO Threshold
#' @param y2l Int. Small SV RO Threshold
#'
#' @export
#'
#' @examples Hunt("~/user/R/SVCROWSin", "~/user/R/FeatureList", "~/user/R/SVCROWSout", TRUE, FALSE, 5000, 25000, 500, 2500, 50, 80)
#' @examples Hunt("~/user/R/SVCROWSin", "~/user/R/FeatureList", "~/user/R/SVCROWSout", TRUE, TRUE)
Hunt <- function(InputQueryList, FeatureList, OutputDirectory, BPfactor = TRUE, DefaultSizes = FALSE,  xs = NA, xl = NA, y1s = NA, y1l = NA, y2s = NA, y2l = NA)
{

  library(dplyr)
  library(tidyverse)
  library(purrr)

  #size cutoffs for CNVs (basepairs)
  xSmallSVL <<- xs
  xLargeSVL <<- xl

  #Boundry sizes for overlapping SV start sites.
  #Size of the whole boundry (basepairs)
  y1SmallBound <<- y1s
  y1LargeBound <<- y1l

  #Percent overlap threshold for small and large size SVs (%)
  y2SmallOV <<- y2s
  y2LargeOV <<- y2l


  FL <- read_tsv(FeatureList, show_col_types = FALSE)
  if (!DefaultSizes)
  {
    FeatureListIn <- CreateFL(FL, xSmallSVL, xLargeSVL)
  }


  inputlist <- InputQueryList
  files <- list.files(inputlist, full.names = TRUE)

  for (file in files)
  {
    if (DefaultSizes)
    {
      RawInputTable <- read_tsv(file, show_col_types = FALSE)
      quants <- quantile(abs(RawInputTable$Length))
      xSmallSVL <<- round(as.numeric(quants[2]))
      xLargeSVL <<- round(as.numeric(quants[4]))

      y1SmallBound <<- round(xSmallSVL / 10)
      y1LargeBound <<- round(xLargeSVL / 10)

      y2SmallOV <<- 40
      y2LargeOV <<- 70

      FeatureListIn <- CreateFL(FL, xSmallSVL, xLargeSVL)
    }
    filename <<- basename(file)
    first4 <<- substring(filename, 1, 20)

    OutPuts <- RunHunt(file, FeatureListIn, BPfactor)

    WriteHunt(OutPuts, first4, OutputDirectory)
  }
}

WriteHunt <- function(OutPuts, name, OutputDirectory)
{
  FinalConsensusList <<- OutPuts$FCL
  FinalPerSampleList <<- OutPuts$FPSL
  AdjustedFinalPerSampleList <<- OutPuts$AFPSL
  QueryList <<- OutPuts$QL


  QLName <- paste0(OutputDirectory,"/",name,".hunted.QL.tsv")
  FCLName <- paste0(OutputDirectory,"/",name,".hunted.FCL.tsv")
  FPSLName <- paste0(OutputDirectory,"/",name,".hunted.FPSL.tsv")
  AFPSLName <- paste0(OutputDirectory,"/",name,".hunted.AFPSL.tsv")


  write_tsv(QueryList, file = QLName)
  write_tsv(FinalConsensusList, file = FCLName)
  write_tsv(FinalPerSampleList, file = FPSLName)
  write_tsv(AdjustedFinalPerSampleList, file = AFPSLName)

  suppressWarnings(rm(list = names(which(!unlist(eapply(.GlobalEnv,
                                       \(x) inherits(x, what = "function")))))))

  print(paste0("Finished processing ",name))
}


RunHunt <- function(file, FeatureListIn, BPfactor)
{

  QueryListIn <- CreateQL(file)
  PerSampleIn <- CreateOQL(QueryListIn)
  OP <- SVCrowsHunt(QueryListIn, FeatureListIn, PerSampleIn, BPfactor)
  return(OP)
}



