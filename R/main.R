##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '
## Aug 2011 nanomed project (Mike Milones human cd4/cd8 data)
## Bonneau lab - "Aviv Madar" <am2654@nyu.edu>, 
## NYU - Center for Genomics and Systems Biology
##  .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.
## /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ /|/ \|\ / / \ \ / / \ \
##`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   ' '

#----------------#
# 0-
#----------------#
cat("\nNanomed run:\n")
# set global variables 
cat("\n- setting global variables\n")
rm(list=ls())
GLOBAL <- list()
GLOBAL[["run.these.steps"]] <- 1:2
# stamp the date on this run
x <- unlist(strsplit(date()," +",perl=TRUE))
GLOBAL[["date.is"]] <- paste(x[2],x[3],x[5],sep="_")

#----------------#
# 1-
#----------------#
# Normalize the dustin data together with the immgen data.
#    -- script:
#       - r_scripts/normalize_immgen_with_dustin_data.R
#    -- input: 
#       - input/immgen_Mar_22_2011
#    -- output:
#       - input/nanomed_and_dustin_data_este.RData
#       - input/nanomed_and_dustin_data_matrix.RData
if(any(GLOBAL[["run.these.steps"]]==1)) {
  cat("\n- Normalizing milone data\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  date.is <- GLOBAL[["date.is"]]
  source("r_scripts/normalize_data.R")
}


#----------------#
# 2-
#----------------#
# Calculate differential expression zscores (heuristic as we don't have enough repeats)
#    -- script:
#       - r_scripts/normalize_immgen_with_dustin_data.R
#    -- input: 
#       - input/immgen_Mar_22_2011
#    -- output:
#       - input/nanomed_and_dustin_data_este.RData
#       - input/nanomed_and_dustin_data_matrix.RData
if(any(GLOBAL[["run.these.steps"]]==2)) {
  cat("\n- Calculate differential expression zscores\n")
  rm(list=ls()[-which(ls()=="GLOBAL")])
  ## date.is <- GLOBAL[["date.is"]]
  date.is <- "Aug_22_2011"
  source("r_scripts/diff_exp.R")
}
