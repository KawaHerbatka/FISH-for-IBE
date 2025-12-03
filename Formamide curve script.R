#Formamide curve creator for oligoNdesign results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#needed installations

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
install.packages("readr")

  #to install tidyverse on ubuntu it might be required to run in terminal:
    # sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev libharfbuzz-dev libfribidi-dev libfontconfig1-dev libtiff5-dev
  # and then, in R
    # install.packages("ragg")

install.packages("tidyverse")

#loading libraries

library(readr)
library(DECIPHER)
library(tidyverse)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#specifying parameters 
  # filepath is the path to the file with SORTED and FILTERED oligos from OligoNdesign
  # testfilePath should include a path to fasta file that contains sequence(s) with a gene from excluding/non-target set.
  # Ran - end range of formamide concentration for drawing the formamide curve. Should be within range (0;100>
  # n,m - the range from n to m oligos for which to create the formamide curve. It might be advantageous to add the starting 
  # efficiency to the dataframe and sort by it, in order to choose best candidates
  # oligoList - if you want to specify given oligos instead of their range, then you can add them to the list. F. eg. c(1,5,8)
  # IMPORTANT!! If you want to use oligoList instead of range n,m then you need to set isList to TRUE
  # If you want to test one probe against many nontargets you need to provide several nontarget sequences in the file under
  # testfilePath and set the parameter manyOffTargets to TRUE. Only the first probe will be tested in this manner. You can
  # change the tested probe by using oligoList. If manyOffTargets is set to FALSE, all probes will be tested against the first
  # nontarget sequence from the testfilePath file.

filepath = "/path/to/probe/data"
testfilePath = "/path/to/test/sequences.fasta"

Ran = 90

n = 1
m = 4

oligoList = c(6,8,13)

isList = TRUE
manyOffTarg = TRUE

#reading data

df<-readr::read_tsv(filepath)
test_df <- as.data.frame(seqinr::read.fasta(file = testfilePath, as.string = TRUE, forceDNAtolower = FALSE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ADDING START EFFICIENCY TO THE TABLE

probe <- df$sequence[]
targets <- reverseComplement(DNAStringSet(probe))

f <- function(FA)
  CalculateEfficiencyFISH(probe, targets, temp=46, P=250e-9, ions=1,FA)[,"HybEff"]
efficiency_all <- c(unlist(lapply(0, f)), ncol=nrow(df))

df$start_eff <- head(efficiency_all,-1)

#sorting the table by the start efficiency, descending

df = arrange(df, desc(df$start_eff))

#Alternatively, sort by other parameters, whatever you prefer

#df = arrange(df, df$class, desc(df$start_eff))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##GRAPH - several probes, one (first) sequence as non-target

FA_range <- 0:Ran # [FA] (%, v/v)

if(manyOffTarg){
  col = 1
  if (!isList){ 
    probe <- df$sequence[n]
  } else {
    probe <- df$sequence[oligoList[1]]
  }
} else {
  if (!isList){ 
    probe <- df$sequence[n:m]
    col = m-n+1
  } else {
    probe <- df$sequence[oligoList]
    col = length(oligoList)
  }
}

targets <- reverseComplement(DNAStringSet(probe))

f <- function(FA)
  CalculateEfficiencyFISH(probe, targets,
                          temp=46, P=250e-9, ions=1, FA)[, "HybEff"]
efficiency <- matrix(unlist(lapply(FA_range, f)), ncol=col, byrow=TRUE)
matplot(FA_range, efficiency, ylim=c(0,1), ylab="Hybridization Efficiency",
        xlab=expression(paste("[Formamide] (%, v/v)", sep="")),
        type="l", lwd=2, col="Blue", main="Formamide Curve", lty=1)

make_imputs <- function(a,b){
  vec <- c(as.character(b[1]))
  for(x in 1:length(a))
  {
    if (x != 1)
    {
      vec <- append(vec, c(as.character(b[1])))
    }
  }
  return(vec)
}

if(!manyOffTarg){
  nontargets <- DNAStringSet(as.character(make_imputs(probe,test_df)))
} else
{
  probe <- DNAStringSet(as.character(make_imputs(test_df,probe)))
  nontargets <- DNAStringSet(as.character(test_df))
  col = length(test_df)
}


f <- function(FA)
  CalculateEfficiencyFISH(probe, nontargets,
                          temp=46, P=250e-9, ions=1, FA)[, "HybEff"]
efficiency <- matrix(unlist(lapply(FA_range, f)), ncol=col, byrow=TRUE)
matlines(FA_range, efficiency, col="Red", lwd=2, lty=3)
abline(h=0.5, lty=2, lwd=2, col="Orange")
abline(v=35, lty=2, lwd=2, col="Green")
legend("topright", legend=c("Targets", "Non-Targets", "50% Efficiency",
                            "Experimental [FA]"), col=c("Blue", "Red", "Orange", "Green"),
       lwd=c(2, 2, 2, 2), lty=c(1, 3, 2, 2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# to check for hairpins and self-dimers - you can either use oligoNdesign command testThorough or refer to this website:
# http://oligocalc.eu







