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
#' Convert a consensus style list to a query style list
#' @description
#' ! re-IDs the input !
#'
#' @param fileIn Consensus List (with path) , in the CL format. Typically the output of a previous run
#' @param fileOut Query List (with path), in the QL format.
#'
#' @export
#'
#' @examples ConsensusToQuery("~/user/R/ConsensusList.tsv", "~/user/R/QueryListOut.tsv")
ConsensusToQuery <- function(fileIn, fileOut)
{
  CList <- read_tsv(fileIn)

  CList <- CList %>% arrange(desc(abs(Length)))
  QList <- CList[, c(1,2,3,4,5,6,7,8,9,10,11,12)]
  colnames(QList)[12] <- "NumReads"
  colnames(QList)[11] <- "QScore"

  for(i in seq_len(nrow(QList)))
  {
    QList$NumReads[i] <- 1
    QList$QScore[i] <- 1
    #QList$ID[i] <- CList$ID[i] #keeps IDs consistant
    QList$ID[i] <- i #RE-IDs the entire list
    #QList$NumReads[i] <- CList$NumReads[i]
    #QList$QScore[i] <- CList$QScore[i]
  }
  write_tsv(Biglist, file = fileOut)
}


#' Take directory of CL files and concatenate, sort, and re-ID them. This preserves all values from the original CL
#'
#' @description
#' If a user has many samples they may want to run each sample by itself first, to remove redundancy within a sample, and then combine all sample CLs to do an dataset-wide analysis.
#' I.e. single sample QL -> Single sample CL -> ConcatSamples (multisample) -> ConsensusToQuery.
#'
#' @param Indirectory Specify directory of CL files and only CL files There is no Type checking, it will just fail if not in the correct format.
#' @param Outdirectory Specify an directory including the desired name of the output file to write the concatenated list.
#'
#' @return A Query List (QL) style output using the origional CL information
#' @export
#'
#' @examples ConcatSamples(~user/R/BunchofCLs, ~User/R/concatOut/QL.out.tsv)
ConcatSamples <- function(Indirectory, Outdirectory)
{
  path <- directory

  files <- list.files(path, full.names = TRUE)

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

  Biglist <- CList

  for (file in files)
  {
    toadd <- read_tsv(file)
    Biglist <- rbind(Biglist, toadd)
  }

  Biglist <- Biglist[, c(1,2,3,4,5,6,7,8,9,10,11,12)]

  Biglist <- Biglist %>% arrange(desc(abs(Length)))
  for(i in seq_len(nrow(Biglist)))
  {
    Biglist$ID[i] <- i
  }

  name <- Outdirectory
  write_tsv(Biglist, file = name)
}





#' SVCSummarize
#'
#' @param OutputDirectory Directory where the output files are being written to from either Hunt() or Scavenge()
#' @param OutputFile What the summary file will be written as
#'
#' @export
#'
#' @examples SVCSummarize("~/R/SVCROWS2/SVCROWSout", "~/R/SVCROWS2/SVCROWSout/summary.tsv")
SVCSummarize <- function(OutputDirectory, OutputFile)
{
  # This code will loop through a folder containing all of your output files from multiple SVCROWS Scavenge runs and extract specific summary statistics from the Final Consensus Lists (FCLs).

  #Creates a list of files from the folder that end in ".FCL.tsv"
  file_list <- list.files(path = OutputDirectory, pattern = "\\.FCL\\.tsv$", full.names = TRUE)
  # Create a data frame to store the data extracted from the FCLs.
  results <- data.frame(FileName = character(),
                        TotalRegions = integer(),
                        MeanCNVSize = numeric(),
                        MedianCNVSize = numeric(),
                        CommonCNVFreq = numeric(),
                        RareCNVFreq = numeric(),
                        SmallCNVProportion = numeric(),
                        MediumCNVProportion = numeric(),
                        LargeCNVProportion = numeric(),
                        DuplicationCount = integer(),
                        DeletionCount = integer(),
                        DecomplexityFactor = numeric(),
                        stringsAsFactors = FALSE)
  for (file in file_list)
    {
    data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    #Mean and median of the (absolute value) CNV lengths.
    mean_cnv_size <- mean(abs(data$Length), na.rm = TRUE)
    median_cnv_size <- median(abs(data$Length), na.rm = TRUE)
    #Optional summary statistics.
    total_rows <- nrow(data)
    common_cnv_freq <- sum(data$Rarity == "Common", na.rm = TRUE) / total_rows
    rare_cnv_freq <- sum(data$Rarity %in% c("Unique", "Rare"), na.rm = TRUE) / total_rows
    small_cnv_prop <- sum(data$Size == "Small", na.rm = TRUE) / total_rows
    medium_cnv_prop <- sum(data$Size == "Medium", na.rm = TRUE) / total_rows
    large_cnv_prop <- sum(data$Size == "Large", na.rm = TRUE) / total_rows
    duplication_count <- sum(data$Type == "DUP", na.rm = TRUE)
    deletion_count <- sum(data$Type == "DEL", na.rm = TRUE)
    decomplexity_factor <- sum(data$ROCount, na.rm = TRUE) / total_rows

     #Copies the results from the functions above into the dataframe
    results <- rbind(results, data.frame(FileName = basename(file),
                                         TotalRegions = total_rows,
                                         MeanCNVSize = mean_cnv_size,
                                         MedianCNVSize = median_cnv_size,
                                         CommonCNVFreq = common_cnv_freq,
                                         RareCNVFreq = rare_cnv_freq,
                                         SmallCNVProportion = small_cnv_prop,
                                         MediumCNVProportion = medium_cnv_prop,
                                         LargeCNVProportion = large_cnv_prop,
                                         DuplicationCount = duplication_count,
                                         DeletionCount = deletion_count,
                                         DecomplexityFactor = decomplexity_factor,
                                         stringsAsFactors = FALSE))
  }
  #Creates a TSV with the resultant summary statistics.
  write.table(results, OutputFile, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Statistical analyses completed. Results saved to ", OutputDirectory, "\n")
}
