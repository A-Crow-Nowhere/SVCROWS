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



###User Inputs and Input Query list
# ordered by size largest to smallest prefers larger SVs // vice versa.
# Also ordered with "known" CNVs first, to be added first to consensus.
# Read user inputs, create necessary input and output files. User is required to
# Input 6 variables to categorize SVs based on their size. User can run the
# program in without "weighted sizes" by entering the same number into each pair
# of variables below.
CreateQL <- function(fileIn)
{
  #File path to table containing all CNVs from all samples
  #See SV crows Style guide for format
  RawInputTable <- read_tsv(fileIn, show_col_types = FALSE) #add path

  #Order SVs by size, put CNVs designated as "known" at the top
  #Split Tables by Known status
  Known <<- RawInputTable %>% filter(IsKnown == TRUE)
  NotKnown <- RawInputTable %>% filter(IsKnown == FALSE)

  #Arrange Tables
  Known <- Known %>% arrange(desc(abs(Length)))
  NotKnown <- NotKnown %>% arrange(desc(abs(Length)))

  knownfix <<- RemoveInnerKnowns(Known)


  if (is.null(knownfix))
  {
    return(NotKnown)
  }
  else
  {
    #Reconcatonate the tables
    knownfix <- knownfix %>% arrange(desc(abs(Length)))
    QList <- rbind(Known, NotKnown)
    return(QList)
  }
}


##Make Knowns more relevant
# If a known cnv is completley overlapped by a larger one, then it will never be
# counted in this analysis. i.e. they will never be iterated. Here we 'merge'
# those CNVs into the larger one by removing the smaller
RemoveInnerKnowns <- function(knowns)
{
  #double loop scans through the list for overlaps
  tempKnowns <- knowns
  for(i in seq_len(nrow(tempKnowns)))
  {
    if (i == nrow(tempKnowns))
    {
      return(tempKnowns)
    }
    else
    {
      j = i + 1
      while (j <= (nrow(tempKnowns)))
      {
        if (tempKnowns$Chr[i] == tempKnowns$Chr[j] && tempKnowns$Type[i] == tempKnowns$Type[j] && tempKnowns$Start[i] < tempKnowns$Start[j] && tempKnowns$End[i] > tempKnowns$End[j])
        {
          tempKnowns <- tempKnowns[-j,]
        }

        j = j + 1
      }
    }
  }
}


###Create Consensus List
# Consensus List is the final output list that merges all of the individual CNVs
# together if a match is not found. If a match is found, The matching variables
# will iterate by 1, indicating that something else matched to it.
CreateCL <- function()
{
  #Column Names for the Consensus list
  CL_columns <- c(
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
    "IsKnown",

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

    "NumReads",
    "QScore",

    "Frequency",
    "Rarity",
    "Matches"
  )

  #Make Tibble with correct column names (i.e. the consensus list)
  CList = data.frame(matrix(nrow = 0, ncol = length(CL_columns)))
  colnames(CList) = CL_columns

  return(CList)
}


###Create Output Query Table
# Same as the Query list with column added for which Consensus ID was matched to
# by the query. Reffered to as the 'PerSampleList'.
CreateOQL <- function(QL)
{

  OQL <- QL
  #0s in the output indicate a match was not found
  OQL <- OQL %>% mutate(BPStartParentID = 0)
  OQL <- OQL %>% mutate(BPEndParentID = 0)
  OQL <- OQL %>% mutate(ROParentID = 0)
  OQL <- OQL %>% mutate(Rarity = "Unassigned")

  return(OQL)
}



##Part 2: Main Loop
# Handles List iteration, calls other functions as relevant. Returns the final
# consensus list and the per-sample list.
SVCrowsScavenge <- function(QueryListX, ConsensusListX, PerSampleX, BPfactor, ExpandRORegion)
{

  WorkingCon <<- ConsensusListX
  QueryList <<- QueryListX
  PerSampleList <<- PerSampleX

  #first entry in the list is tricky, because its empty. Set this variable TRUE
  #when the list is first empty
  EmptyList <<- TRUE

  #Begin Query loop. Each query is assessed only once, Either a match is found
  #and is iterated, or it finds no matches and is added to the end of the
  #Consensus lsit
  for(i in seq_len(nrow(QueryList)))
  {
    QueryPosition <<- i
    print(paste("sQP: ",QueryPosition))

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

    #add the first query to the CL by default
    if (EmptyList)
    {
      EmptyList <<- FALSE
      rowIn <- AddQuery(qIsKnown, qLength, QueryList[QueryPosition,])
      WorkingCon <<- rbind(WorkingCon, rowIn)
      suppressWarnings(remove(rowIn))

      #temp files are used to not interfere with original files as much as
      #possible
      tempProw <<- PerSampleList[QueryPosition,]

      tempProw$ROParentID[1] <- qID
      tempProw$BPStartParentID[1] <- qID
      tempProw$BPEndParentID[1] <- qID

      tempPerSamp <- PerSampleList
      tempPerSamp[QueryPosition,] <- tempProw
      PerSampleList <<- tempPerSamp

      suppressWarnings(remove(tempProw))
    }


    else
    {
      #split the working consensus into entries that match the Chr and Type
      #for speed
      twc <- WorkingCon
      MatchList <<- twc %>%
        filter(
          Chr == qChr,
          Type == qType,
         # !Var3 == qV3, # Optional conditions for list matching cna be applied
          #Var2 != qV2, #
          #nothing ouside of these bounds are relevant to search
          BPEndRight >= (qStart - y1LargeBound),
          BPStartLeft <= (qEnd + y1LargeBound)
        )

      # Filter conditions combined into a single filter function call for MinusMatchList
      MinusMatchList <<- twc %>%
        filter(
          Chr != qChr |
            Type != qType |
            #Var3 == qV3 | # Optional conditions for list matching cna be applied
            #Var2 == qV2 | # Still must be inverses of the previous block^
            BPEndRight < (qStart - y1LargeBound) |
            BPStartLeft > (qEnd + y1LargeBound)
        )
      suppressWarnings(remove(WorkingCon))

      #If it is the  first of it's kind, add the object to the match list
      if (as.numeric(nrow(MatchList)) == 0)
      {
        rowIn <- AddQuery(qIsKnown, qLength, QueryList[QueryPosition,])

        tempML <- MatchList
        tempML <- rbind(tempML, rowIn)
        MatchList <<- tempML

        tempProw <<- PerSampleList[QueryPosition,]
        tempProw$ROParentID[1] <- qID
        tempProw$BPStartParentID[1] <- qID
        tempProw$BPEndParentID[1] <- qID

        tempPerSamp <- PerSampleList
        tempPerSamp[QueryPosition,] <- tempProw
        PerSampleList <<- tempPerSamp
        suppressWarnings(remove(tempProw))

        TwC <- ReMergeLists(MatchList, MinusMatchList)
        WorkingCon <<- TwC

        ConsensusPosition <<- 1
      }

      #There are matches to assess
      else
      {
        #Index for the CL, not a for loop because it needs to end prematurely.
        ConsensusPosition <<- 1
        while(!QueryAdded)
        {
          #Set variables for convenience, all of these get passed to subsequent
          #functions. rX stands for 'reference'. Reference is the current CL entry
          #the query is being compared to.
          rChr <<- MatchList$Chr[ConsensusPosition]
          rStart <- MatchList$Start[ConsensusPosition]
          rEnd <- MatchList$End[ConsensusPosition]
          rLength <- MatchList$Length[ConsensusPosition]
          rType <- MatchList$Type[ConsensusPosition]
          rID <- MatchList$ID[ConsensusPosition]
          rV1 <- MatchList$Var1[ConsensusPosition]
          rV2 <- MatchList$Var2[ConsensusPosition]
          rV3 <- MatchList$Var3[ConsensusPosition]
          rIsKnown <- MatchList$IsKnown[ConsensusPosition]

          rReads <- MatchList$NumReads[ConsensusPosition]
          rScore <- MatchList$QScore[ConsensusPosition]

          rROPcntPass <- MatchList$ROPercentPass[ConsensusPosition]

          rStartLeft <- MatchList$BPStartLeft[ConsensusPosition]
          rStartRight <- MatchList$BPStartRight[ConsensusPosition]
          rEndLeft <- MatchList$BPEndLeft[ConsensusPosition]
          rEndRight <- MatchList$BPEndRight[ConsensusPosition]

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
          ROPass <- RecpricolOverlapCalc(qStart, rStart, qEnd, rEnd, rIsKnown, qIsKnown, qLength, rLength, rROPcntPass, FALSE, ConsensusPosition)
          #Reset Passing threshold
          rROPcntPass <-  TempROPercentPass

          #If there is no match on anything, and we aren't at the end of the list,
          #just iterate the Consensus Position
          if (BPStartPass == FALSE && BPEndPass == FALSE && ROPass == FALSE && (ConsensusPosition != nrow(MatchList)))
          {
            tempPos <- ConsensusPosition
            ConsensusPosition <<- (tempPos + 1)
          }

          #If there is no match, and we are at the end of the list, add the new
          #entry
          else if (BPStartPass == FALSE && BPEndPass == FALSE && ROPass == FALSE && (ConsensusPosition == nrow(MatchList)))
          {

            rowIn <- AddQuery(qIsKnown, qLength, QueryList[QueryPosition,])
            MatchList <<- rbind(MatchList, rowIn)

            tempProw <<- PerSampleList[QueryPosition,]

            tempProw$ROParentID[1] <- qID
            tempProw$BPStartParentID[1] <- qID
            tempProw$BPEndParentID[1] <- qID

            tempPerSamp <- PerSampleList
            tempPerSamp[QueryPosition,] <- tempProw
            PerSampleList <<- tempPerSamp
            suppressWarnings(remove(tempProw))

            ReMergeLists(MatchList, MinusMatchList)

            ConsensusPosition <<- 1
            QueryAdded <<- TRUE
          }
          else
          {

            #Based on Passing Criteria, enter the descion tree
            rowIn <- Decision(MatchList[ConsensusPosition,], QueryList[QueryPosition,],
                              qIsKnown, rIsKnown, BPStartPass, BPEndPass,
                              BPStartHasPassed, BPEndHasPassed, ROPass,
                              PerSampleList[QueryPosition,],
                              qLength, rID, rReads, rScore, qReads, qScore, ConsensusPosition, FALSE)



            if(ExpandRORegion && ROPass)
            {
              if((abs(rowIn$Crow$Start[1]) > abs(qStart)) || (abs(rowIn$Crow$End[1]) < abs(qEnd)))
              {
                if (rowIn$Crow$Length[1] < 0)
                {
                  makeitnegative <- TRUE
                }

                if(abs(rowIn$Crow$Start[1]) > abs(qStart))
                {
                  rowIn$Crow$Start[1] <- abs(qStart)
                }
                if(abs(rowIn$Crow$End[1]) < abs(qEnd))
                {
                  rowIn$Crow$End[1] <- abs(qEnd)
                }
                stringIn <- rowIn$Crow$Matches[1]
                StringOut <- sub("^(.*):", "\\1+", stringIn)
                rowIn$Crow$Matches[1] <- StringOut


                rowIn$Crow$Length[1] <- (rowIn$Crow$End[1] - rowIn$Crow$Start[1])
                if(makeitnegative)
                {
                  rowIn$Crow$Length[1] <- (rowIn$Crow$Length[1] * -1)
                }


                newLength <- rowIn$Crow$Length[1]
                if (abs(newLength) <= xSmallSVL)
                {
                  rowIn$Crow$Size[1] <- "Small"
                  Sizes <- SizeDetermination(xSmallSVL)
                }
                else if (abs(newLength) <= xLargeSVL)
                {
                  rowIn$Crow$Size[1] <- "Medium"
                  Sizes <- SizeDetermination(newLength)
                }
                else
                {
                  rowIn$Crow$Size[1] <- "Large"
                  Sizes <- SizeDetermination(xLargeSVL)
                }

                rowIn$Crow$BPStartLeft[1] <- round(rowIn$Crow$Start[1] - (Sizes[1]/2))
                rowIn$Crow$BPStartRight[1] <- round(rowIn$Crow$Start[1] + (Sizes[1]/2))

                rowIn$Crow$BPEndLeft[1] <- round(rowIn$Crow$End[1] - (Sizes[1]/2))
                rowIn$Crow$BPEndRight[1] <- round(rowIn$Crow$End[1] + (Sizes[1]/2))

                rowIn$Crow$BPBoundSize[1] <- Sizes[1]
              }
            }


            #If a match is found in the query, update the CL and PSL
            if ((BPStartPass || BPEndPass || ROPass))
            {
              #Update rows based on results of decision
              tempWorkC <- MatchList
              tempWorkC[ConsensusPosition,] <- rowIn$Crow
              MatchList <<- tempWorkC

              tempPerSamp <- PerSampleList
              tempPerSamp[QueryPosition,] <- rowIn$Prow
              PerSampleList <<- tempPerSamp


              #if RO is the thing that passed, this Query is done.
              #Also accepts both BPs matching (RO is iterated in Decision())
              if (ROPass || (BPStartPass && BPEndPass))
              {
                QueryAdded <<- TRUE
              }

              #only the breakpoints passed and not at end, we can continue searching for a match
              else if ((BPStartPass || BPEndPass) && !ROPass && (ConsensusPosition != nrow(MatchList)))
              {
                tempPos <- ConsensusPosition
                ConsensusPosition <<- (tempPos + 1)
              }

              #a breakpoint passed and it is the end of the CL, update and then add
              else if (((BPStartPass || BPEndPass)) && (ConsensusPosition == nrow(MatchList)))
              {
                rowIn <- AddQuery(qIsKnown, qLength, QueryList[QueryPosition,])

                tempML <- MatchList
                tempML <- rbind(tempML, rowIn)
                MatchList <<- tempML

                tempProw <<- PerSampleList[QueryPosition,]

                if(BPStartPass)
                {
                  tempProw$BPStartParentID[1] <- rID
                }
                else
                {
                  tempProw$BPStartParentID[1] <- qID
                }

                if(BPEndPass)
                {
                  tempProw$BPEndParentID[1] <- rID
                }
                else
                {
                  tempProw$BPEndParentID[1] <- qID
                }

                tempProw$ROParentID[1] <- qID

                tempPerSamp <- PerSampleList
                tempPerSamp[QueryPosition,] <- tempProw
                PerSampleList <<- tempPerSamp

                TwC <- ReMergeLists(MatchList, MinusMatchList)
                WorkingCon <<- TwC

                QueryAdded <<- TRUE
              }
              suppressWarnings(remove(rowIn, TwC))
            }
          }
          TwC <- ReMergeLists(MatchList, MinusMatchList)
          WorkingCon <<- TwC
          suppressWarnings(remove(TwC))
        }
      }
    }
  }

  print("Done with comparisons, generating outputs")
  FCLout <<- FCLprocessing(WorkingCon)
  FPSLout <<- PSLprocessing(PerSampleList, FCLout)
  AFPSLout <<- AdjustPSL(FPSLout)
  outs <- list(FCL = FCLout, FPSL = FPSLout, AFPSL = AFPSLout, QL = QueryList)
  return(outs)
}

#Lists are split for speed, must be remerged in between resets
ReMergeLists <- function(Mlist, MMlist)
{

  #The While loop concludes, remerge the match and nonmatch lists
  tempWC <- rbind(Mlist, MMlist)

  Known <- tempWC %>% filter(IsKnown == TRUE)
  NotKnown <- tempWC %>% filter(IsKnown == FALSE)
  Known <- Known %>% arrange(desc(abs(Length)))
  NotKnown <- NotKnown %>% arrange(desc(abs(Length)))

  tempWC <- rbind(Known, NotKnown)

  return(tempWC)
}

# If query start is within the bounds of the acceptable reference start, return
# true. A query start cannot match with a reference end.
BreakPointStartCalc <- function(qStart, rStartLeft, rStartRight)
{
  if((qStart >= rStartLeft && qStart <=  rStartRight))
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

# If query end is within the bounds of the acceptable reference end, return
# true. A query end cannot match with a reference start.
BreakPointEndCalc <- function(qEnd, rEndLeft, rEndRight)
{
  if((qEnd >= rEndLeft && qEnd <= rEndRight))
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

#function to calculate RO bounds/ check if overlap passes a threshold.
RecpricolOverlapCalc <- function(qStart, rStart, qEnd, rEnd, rIsKnown, qIsKnown, qLength, rLength, rROPcntPass, isHunt, MLposition)
{
  if((qStart < rStart && qEnd < rStart) || (qStart > rEnd && qEnd > rEnd) && !rIsKnown)
  {
    return(FALSE)
  }
  if((qStart < rStart && qEnd > rStart) && rIsKnown)
  {
    return(TRUE)
  }

  #if Query is known and doesnt overlap with a known, automatic fail.
  if(!rIsKnown && qIsKnown)
  {
    return(FALSE)
  }

  # calc left and right bounds of the overlap
  ROLeftBound <- max(qStart, rStart)
  RORightBound <- min(qEnd, rEnd)
  ROBoundLength <- RORightBound - ROLeftBound
  ROPercentPass <- WorkingCon$ROPercentPass[MLposition]

  # Calc overlap ratio for the query and reference
  query_overlap_ratio <- ROBoundLength / abs(qLength)
  ref_overlap_ratio <- ROBoundLength / abs(rLength)

  # if the overlap length is 0 or negative, return FALSE. This should never
  # with valid input
  if (ROBoundLength <= 0)
  {
    return(FALSE)
  }

  #Hunt mode only considers the overlap of the Feature list, and has no Knowns
  if(isHunt)
  {
    pass <- (100*ref_overlap_ratio >= rROPcntPass)
    return(pass)
  }
  else
  {
    #if both are known, query must overlap reference by 50% or more to call a
    #pass. 10% is arbitrary, it does not have to match recipricoly.
    if (rIsKnown && qIsKnown && ref_overlap_ratio <= .1)
    {
      return(FALSE)
    }
    else if (rIsKnown && qIsKnown && ref_overlap_ratio > .1)
    {
      return(TRUE)
    }

    #if neither is known, must calculate ROpass
    else
    {
      # determine if both overlap ratios exceed the ROPercentPass threshold
      pass <- (100*query_overlap_ratio >= rROPcntPass) && (100*ref_overlap_ratio >= rROPcntPass)
      return(pass)
    }

  }
}

##Part 3: Merge Information

###Decision
# This controls the fate of the query based on the results of BPStart/EndPass and
# ROPass. ROPass = TRUE is the only way a CL entry will be iterated. An entry will
# only be added if it is at the end of the consensus list.
Decision <- function(rRow, qRow, qIsKnown, rIsKnown, BPStartPass, BPEndPass, BPStartHasPassed, BPEndHasPassed, ROPass, pRow, qLength, rID, rReads, rScore, qReads, qScore, MLposition, isHunt)
{
  #If Something passes, something must be updated
  if (BPStartPass == TRUE || BPEndPass == TRUE || ROPass == TRUE)
  {
    #only one pass
    if (BPStartPass && !BPEndPass && !ROPass)
    {
      if (!BPStartHasPassed)
      {
        updatedrows <- BPStartIterate(rRow, pRow, rID)
        BPStartHasPassed <<- TRUE
      }
      else
      {
        updatedrows <- list(Crow = rRow, Prow = pRow)
      }
      return(updatedrows)
    }
    else if (BPEndPass && !BPStartPass && !ROPass)
    {
      if (!BPStartHasPassed)
      {
        updatedrows <- BPEndIterate(rRow, pRow, rID)
        BPEndHasPassed <<- TRUE
      }
      else
      {
        updatedrows <- list(Crow = rRow, Prow = pRow)
      }
      return(updatedrows)
    }
    else if (!BPEndPass && !BPStartPass && ROPass)
    {
      updatedrows <- ROIterate(rRow, pRow, rID, rReads, rScore, qReads, qScore)
      return(updatedrows)
    }

    #Two Pass
    else if(BPStartPass && !BPEndPass && ROPass)
    {
      if (!BPStartHasPassed)
      {
        itStart <- BPStartIterate(rRow, pRow, rID)
        BPStartHasPassed <<- TRUE
      }
      else
      {
        itStart <- list(Crow = rRow, Prow = pRow)
      }
      itStartRo <- ROIterate(itStart$Crow, itStart$Prow, rID, rReads, rScore, qReads, qScore)
      return(itStartRo)
    }

    else if(BPEndPass && !BPStartPass && ROPass)
    {
      if (!BPEndHasPassed)
      {
        itEnd <- BPEndIterate(rRow, pRow, rID)
        BPEndHasPassed <<- TRUE
      }
      else
      {
        itEnd <- list(Crow = rRow, Prow = pRow)
      }
      itEndRo <- ROIterate(itEnd$Crow, itEnd$Prow, rID, rReads, rScore, qReads, qScore)
      return(itEndRo)
    }

    #All Three Pass
    else if((BPEndPass && BPStartPass && ROPass) || (BPStartPass && BPEndPass && !ROPass))
    {
      if (!BPStartHasPassed)
      {
        itStart <- BPStartIterate(rRow, pRow, rID)
        BPStartHasPassed <<- TRUE
      }
      else
      {
        itStart <- list(Crow = rRow, Prow = pRow)
      }
      #######################
      if (!BPEndHasPassed)
      {
        itStartEnd <- BPEndIterate(itStart$Crow, itStart$Prow, rID)
        BPEndHasPassed <<- TRUE
      }
      else
      {
        itStartEnd <- itStart
      }
      #######################

      itStartEndRO <- ROIterate(itStartEnd$Crow, itStartEnd$Prow, rID, rReads, rScore, qReads, qScore)
      return(itStartEndRO)
    }

    #One of the BP passes, but have already been iterated (i.e. nothing changes)
    # else
    # {
    #   print("aha")
    #   tempPos <- ConsensusPosition
    #   ConsensusPosition <<- (tempPos + 1)
    #   old <- list(Crow = rRow, Prow = pRow)
    #   return(old)
    # }
  }

  #End of the consensus list reached without a match, add new query
  else if (MLposition == nrow(WorkingCon))
  {
    print("it is here")
    if(isHunt)
    {
      return(rRow)
    }
    else
    {
      NewRow <- AddQuery(qIsKnown, qLength, qRow)
      QueryAdded <<- TRUE
      MLposition <<- 1
      new <- list(Crow = NewRow, Prow = pRow)
      return(new)
    }
  }

  else
  {
    print("uh oh")
    return(rRow)
  }

}


###Break Point Start Iterate
# Add 1 to the BreakpointStartCount in the CL, Update PSL with the matched parent
# ID
BPStartIterate <- function(rRow, pRow, rID)
{
  tempRow <- rRow
  tempRow$BPStartCount[1] <- (tempRow$BPStartCount[1] + 1)
  UpdateRrow <- tempRow

  tempPerSamp <- pRow
  tempPerSamp$BPStartParentID[1] <- rID
  UpdateProw <- tempPerSamp

  outlist <- list(Crow = UpdateRrow, Prow = UpdateProw)
  return(outlist)
}


###Break Point End Iterate
# Add 1 to the BreakpointEndCount in the CL, Update PSL with the matched parent
# ID
BPEndIterate <- function(rRow, pRow, rID)
{
  tempRow <- rRow
  tempRow$BPEndCount[1] <- (tempRow$BPEndCount[1] + 1)
  UpdateRrow <- tempRow

  tempPerSamp <- pRow
  tempPerSamp$BPEndParentID[1] <- rID
  UpdateProw <- tempPerSamp

  outlist <- list(Crow = UpdateRrow, Prow = UpdateProw)
  return(outlist)
}


###Recpricol Overlap Iterate
# Add 1 to the ReciprocolOverlapCount in the CL, Update PSL with the matched parent
# ID
ROIterate <- function(rRow, pRow, rID, rNR, rQS, qNR, qQS)
{

  addNR <- (rNR + qNR)
  avgQS <- ((rRow$QScore[1] * rRow$ROCount[1]) + qQS) / (rRow$ROCount[1] + 1)
  catmatch <- paste0(rRow$Matches[1], ":", pRow$ID[1])

  tempRow <- rRow
  tempRow$ROCount[1] <- (tempRow$ROCount[1] + 1)
  tempRow$NumReads[1] <- addNR
  tempRow$QScore[1] <- avgQS
  tempRow$Matches[1] <- catmatch
  UpdateRrow <- tempRow

  tempPerSamp <- pRow
  tempPerSamp$ROParentID[1] <- rID
  UpdateProw <- tempPerSamp

  outlist <- list(Crow = UpdateRrow, Prow = UpdateProw)
  return(outlist)
}


###Add Query to end of consesnsus list
# The current query row is used to construct a new CL entry. First we must
# calculate the size of the new CL, and calculate the resulting boundries.
AddQuery <- function(isKnown, qLength, qRow)
{
  if (abs(qLength) <= xSmallSVL)
  {
    sizeCategory <- "Small"
    Sizes <- SizeDetermination(xSmallSVL)
  }
  else if (abs(qLength) <= xLargeSVL)
  {
    sizeCategory <- "Medium"
    Sizes <- SizeDetermination(qLength)
  }
  else
  {
    sizeCategory <- "Large"
    Sizes <- SizeDetermination(xLargeSVL)
  }
  NewEntry <- NewConsensusEntry(Sizes, sizeCategory, isKnown, qRow)
  return(NewEntry)
}


###Determine the size of the boundries for RO and BP
# Returns the sizes of the boundries for RO and BP in a list.
SizeDetermination <- function(size)
{
  #define BP variables in y = mx + B
  BPslope <- (y1LargeBound - y1SmallBound)/(xLargeSVL - xSmallSVL)
  BPintercept <- (y1LargeBound - (xLargeSVL * BPslope))

  BPSize <- (abs(size)*BPslope) + BPintercept
  BPSize <- round(BPSize)

  #define RO variables in y = mx + B
  ROslope <- (y2LargeOV - y2SmallOV)/(xLargeSVL - xSmallSVL)
  ROintercept <- (y2LargeOV - (xLargeSVL * ROslope))

  ROSize <- (abs(size)*ROslope) + ROintercept
  ROSize <- round(ROSize)

  return(c(BPSize, ROSize))
}


###Build the new Concensus entry as new row
NewConsensusEntry <- function(sizes, sizeCategory, isKnown, qRow)
{
  #Create a new consensus entry with calculated values and determined bounds
  in.Chr <- qRow$Chr[1]
  in.Start <- qRow$Start[1]
  in.End <- qRow$End[1]
  in.Length <- qRow$Length[1]
  in.Type <- qRow$Type[1]
  in.ID <- qRow$ID[1]
  in.Var1 <- qRow$Var1[1]
  in.Var2 <- qRow$Var2[1]
  in.Var3 <- qRow$Var3[1]
  in.IsKnown <- isKnown
  in.NumReads <- qRow$NumReads[1]
  in.QScore <- qRow$QScore[1]
  in.Size <- sizeCategory
  in.BPBoundSize <- sizes[1]
  in.ROPercentPass <- sizes[2]

  in.ROCount <- 1

  in.BPStartLeft <- round(in.Start - (sizes[1]/2))
  in.BPStartRight <- round(in.Start + (sizes[1]/2))
  in.BPStartCount <- 1


  in.BPEndLeft <- round(in.End - (sizes[1]/2))
  in.BPEndRight <- round(in.End + (sizes[1]/2))
  in.BPEndCount <- 1

  in.Frequency <- 0
  in.Rarity <- "Unassigned"
  in.Matches <- in.ID

  NewRow <- data.frame(Chr = in.Chr, Start = in.Start, End = in.End, Length = in.Length, Type = in.Type, ID = in.ID, Var1 = in.Var1, Var2 = in.Var2, Var3 = in.Var3, IsKnown = in.IsKnown, QScore = in.QScore, NumReads = in.NumReads, Size = in.Size, BPBoundSize = in.BPBoundSize, ROPercentPass = in.ROPercentPass, ROCount = in.ROCount, BPStartLeft = in.BPStartLeft, BPStartRight = in.BPStartRight, BPStartCount = in.BPStartCount, BPEndLeft = in.BPEndLeft, BPEndRight = in.BPEndRight, BPEndCount = in.BPEndCount, Frequency = in.Frequency, Rarity = in.Rarity, Matches = in.Matches)

  return(NewRow)
}


###Part4 Post Processing
# Calculate frequency and rarity assignments. The breakoffs for these categories
# can be changed by the user
FCLprocessing <- function(rFCLin)
{
  UniqueSamples <- rFCLin %>% summarize(unique_count = n_distinct(Var3))

  #calculate frequency. Frequency > 1 indicates more than one CNV in a sample
  #aligns reference. Especially noticable with known CNVs.
  for(i in seq_len(nrow(rFCLin)))
  {
    rFCLin$Frequency[i] <- as.double((rFCLin$ROCount[i]/UniqueSamples))
  }

  for(i in seq_len(nrow(rFCLin)))
  {
    if(rFCLin$IsKnown[i])
    {
      rFCLin$Rarity[i] <- "Known"
    }
    else if(rFCLin$ROCount[i] == 1)
    {
      rFCLin$Rarity[i] <- "Unique"
    }
    else if(rFCLin$Frequency[i] < 0.10)
    {
      rFCLin$Rarity[i] <- "Rare"
    }
    else if(rFCLin$Frequency[i] >= 0.10)
    {
      rFCLin$Rarity[i] <- "Common"
    }
    else
    {
      rFCLin$Rarity[i] <- "NA"
    }
  }
  return(rFCLin)
}


PSLprocessing <- function(rPLSin, FCLin)
{

  for(i in seq_len(nrow(rPLSin)))
  {
    rowindex <- which(FCLin$ID == rPLSin$ROParentID[i])
    if(suppressWarnings(is.null(rowindex) || length(rowindex) == 0 || is.na(FCLin$ID[i]) || is.na(rPLSin$ROParentID[i])))
    {
      rPLSin$Rarity[i] <- "Unassigned"
    }
    else
    {
      rPLSin$Rarity[i] <- FCLin$Rarity[rowindex]
    }
  }
  return(rPLSin)
}


# If multiple SVs from the same sample match to the same consensus entry, this
#artificially inflates the count. Ajusted PSL marks repeat matches as "clone"s,
# will keep the largest SV as the true one.
AdjustPSL <- function(FPSLout)
{
  header <- colnames(FPSLout)
  AFPSL = data.frame(matrix(nrow = 0, ncol = length(header)))
  colnames(AFPSL) = header

  splitdata <- split(FPSLout,FPSLout$Var3)

  result <- map(splitdata, function(df)
  {

    df <- df %>% arrange(desc(abs(Length)))

    i <- 1
    for(i in seq_len(nrow(df)))
    {
      if (i == nrow(df))
      {
        AFPSL <<- rbind(AFPSL, df)
      }
      else
      {
        j = i + 1
        while (j <= (nrow(df)))
        {
          if (df$ROParentID[i] == df$ROParentID[j] && df$Rarity[j] != "Clone")
          {
            df$ROParentID[j] <- (df$ROParentID[j] * -1)
            df$Rarity[j] <- "Clone"
          }

          j = j + 1
        }
      }
    }
    return(AFPSL)
  })
  return(AFPSL)
}




#' Run SVCROWS in "Scavenge" mode
#'
#' @description This is the main function of SVCROWS. This is a program uses reciprocal overlap (RO) to determine if two regions have a significant level of overlap. This program allows for user input to determine the stringency of the comparisons, while also giving options to use other pieces of evidence to support RO calls. In principle, the program uses the concept of "Weighted Sizes" to weight RO stringency.
#'
#'
#' @param InputQueryList Query List in the designated format (see user manual)
#' @param OutputDirectory Directory to write output files (see user manual)
#' @param ExpandRORegion If TRUE, When entries match, use the minimum and maximum breakpoints of either to define a new region (which will be used for all future comparisons). When FALSE, original bounds are kept
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
#' @examples Scavenge("~/user/R/SVCROWSin", "~/user/R/SVCROWSout", FALSE, TRUE, FALSE, 5000, 25000, 500, 2500, 50, 80)
#' @examples Scavenge("~/user/R/SVCROWSin", "~/user/R/SVCROWSout", TRUE, TRUE, TRUE)
Scavenge <- function(InputQueryList, OutputDirectory, ExpandRORegion = FALSE, BPfactor = TRUE, DefaultSizes = FALSE,  xs = NA, xl = NA, y1s = NA, y1l = NA, y2s = NA, y2l = NA)
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
    }
    filename <<- basename(file)
    first4 <<- substring(filename, 1, 20)

    OutPuts <- RunScavenge(file, FeatureListIn, BPfactor, ExpandRORegion)

    WriteScavenge(OutPuts, first4, OutputDirectory)
  }
}

WriteScavenge <- function(OutPuts, name, OutputDirectory)
{
  FinalConsensusList <<- OutPuts$FCL
  FinalPerSampleList <<- OutPuts$FPSL
  AdjustedFinalPerSampleList <<- OutPuts$AFPSL
  QueryList <<- OutPuts$QL

  FinalConsensusList <- FinalConsensusList %>% rename(Abundance = Frequency)


  QLName <- paste0(OutputDirectory,"/",name,".scavenged.QL.tsv")
  FCLName <- paste0(OutputDirectory,"/",name,".scavenged.FCL.tsv")
  FPSLName <- paste0(OutputDirectory,"/",name,".scavenged.FPSL.tsv")
  AFPSLName <- paste0(OutputDirectory,"/",name,".scavenged.AFPSL.tsv")


  write_tsv(QueryList, file = QLName)
  write_tsv(FinalConsensusList, file = FCLName)
  write_tsv(FinalPerSampleList, file = FPSLName)
  write_tsv(AdjustedFinalPerSampleList, file = AFPSLName)

  suppressWarnings(rm(list = names(which(!unlist(eapply(.GlobalEnv,
                                       \(x) inherits(x, what = "function")))))))

  print(paste0("Finished processing ",name))
}

##Run Mode: Scavenge. This is for comparing SVs in the same list together, and
##devloping the Consensus list. Main funtion of SVCROWS.

RunScavenge <- function(file, FeatureListIn, BPfactor, ExpandRORegion)
{

  QueryListIn <- CreateQL(file)
  ConsensusListIn <- CreateCL()
  PerSampleIn <- CreateOQL(QueryListIn)
  OP <- SVCrowsScavenge(QueryListIn, ConsensusListIn, PerSampleIn, BPfactor, ExpandRORegion)
  return(OP)
}










