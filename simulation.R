library(ChIPsim)

arguments <- commandArgs(trailingOnly=TRUE)
name <- arguments[1]
if (arguments[2] == "-") {
    simulate_background <- TRUE
} else {
    simulate_background <- FALSE
    bed_file <- read.table(arguments[2])
}
genome <- readDNAStringSet(arguments[3])
number_of_reads <- as.numeric(arguments[4])
read_length <- as.numeric(arguments[5])

###### PARAMETERS ######

#background_shape <- 0
background_shape <- 10
#background_scale <- 0
background_scale <- 20

binding_shape <- 7
binding_scale <- 4
binding_enrichment <- 2
binding_r <- 1.5

feature_args_shape <- 1
feature_args_scale <- 2
#feature_args_scale <- 20
control_length <- 50

readDens_args_bind <- 30
#readDens_args_bind <- 50
readDens_args_minLength <- 150
readDens_args_maxLength <- 250
readDens_args_meanLength <- 200

###### GENOME ######

chromosome_lengths = list()
for (chromosome in names(genome)) {
    chromosome_lengths[[chromosome]] <- length(genome[[chromosome]])
}

#######################

# The following part is largely based on ChIPsim vignette.
# Most important difference is using your own coodinates
#  (provided as argument in bed file) instead of simulating them.
# Slight modifications include:
# - fix a typo in for loop in dfReads() function
# - add N nucleotide

# Some parts of it are not important,
#  because we don't want to generate features randomly.
# But apparently we can't omit that, or at least I don't know how.

randomQuality <- function(read, ...){
  paste(sample(unlist(strsplit(rawToChar(as.raw(84:104)),"")),
               nchar(read), replace = TRUE), collapse="")
}

dfReads <- function(readPos, readNames, sequence, readLen, ...){
  ## create vector to hold read sequences and qualities
  readSeq <- character(sum(sapply(readPos, sapply, length)))
  readQual <- character(sum(sapply(readPos, sapply, length)))
  idx <- 1
  ## process read positions for each chromosome and strand
  for(k in 1:length(readPos)){ ## chromosome
    for(i in 1:2){ ## strand
      for(j in 1:length(readPos[[k]][[i]])){
        ## get (true) sequence
        readSeq[idx] <- as.character(readSequence(readPos[[k]][[i]][j], sequence[[k]], 
                                                  strand=ifelse(i==1, 1, -1), readLen=readLen))
        ## get quality
        readQual[idx] <- randomQuality(readSeq[idx])
        ## introduce sequencing errors
        readSeq[idx] <- readError(readSeq[idx], decodeQuality(readQual[idx]))
        idx <- idx + 1
      }
    }
  }
  writeFASTQ(readSeq, readQual, unlist(readNames), paste(name, ".fastq", sep = ""))
  data.frame(name=unlist(readNames), sequence=readSeq, quality=readQual,
             stringsAsFactors = FALSE)
}

backgroundFeature <- function(start, length=1000, shape=background_shape, scale=background_scale){
  weight <- rgamma(1, shape, scale)
  params <- list(start = start, length = length, weight = weight, overlap = 0)
  class(params) <- c("Background", "SimulatedFeature")
  return(params)
}

bindingFeature <- function(start, length,
                           shape=binding_shape, scale=binding_scale,
                           enrichment=binding_enrichment, r=binding_r) {
  avgWeight <- shape * scale * enrichment
  lowerBound <- ((r-1) * avgWeight)
  weight <- actuar::rpareto1(1, r, lowerBound)
  params <- list(start = start, length = length, weight = weight, overlap = 0)
  class(params) <- c("Binding", "SimulatedFeature")
  return(params)
}

generator <- list(Binding = bindingFeature, Background = backgroundFeature)
transition <- list(Binding=c(Background=1), Background=c(Binding=0.05, Background=0.95))
transition <- lapply(transition, "class<-", "StateDistribution")
init <- c(Binding=0, Background=1)
class(init) <- "StateDistribution"

constRegion <- function(weight, length) {
  rep(weight, length)
}
featureDensity.Binding <- function(feature, ...) constRegion(feature$weight, feature$length)
featureDensity.Background <- function(feature, ...) constRegion(feature$weight, feature$length)

fragLength <- function(x, minLength, maxLength, meanLength, ...){
  sd <- (maxLength - minLength)/4
  prob <- dnorm(minLength:maxLength, mean = meanLength, sd = sd)
  prob <- prob/sum(prob)
  prob[x - minLength + 1]
}

parameters <- defaultControl()
parameters$readSequence$readLen <- read_length

defaultErrorProb <- function(){
        prob <- list(A=c(0, 0.14, 0.05, 0.05, 0.01),
                     C=c(0.13, 0, 0.02, 0.04, 0.01),
                     G=c(0.04, 0.08, 0, 0.12, 0.01),
                     T=c(0.08, 0.15, 0.09, 0, 0.01),
                     N=c(0.01, 0.01, 0.01, 0.01, 0))
        prob <- lapply(prob, "names<-", c("A", "C", "G", "T", "N"))
        prob
}

readError <- function(read, qual, alphabet=c("A", "C", "G", "T", "N"), prob=defaultErrorProb(), ...){
        read <- gsub(paste("[^", paste(alphabet, collapse="", sep=""),"]", sep=""), "N", read)
        ## determine location of error
        errorPos <- runif(length(qual)) < qual
        if(any(errorPos)){
                for(i in which(errorPos)){
                        transProb <- prob[[substr(read, i, i)]]
                        substr(read, i, i) <- sample(names(transProb), 1, prob=transProb)
                }
        }
        read
}

parameters$readSequence$errorFun <- readError

myFunctions <- ChIPsim::defaultFunctions()
myFunctions$readSequence <- dfReads

# are those two below needed?
#featureArgs <- list(generator = generator,
#                    transition=transition, init=init, start = 0,
#                    length = 1e6, globals = list(shape = feature_args_shape, scale = feature_args_scale),
#                    experimentType = "TFExperiment", 
#                    lastFeat = c(Binding = FALSE, Background = TRUE),
#                    control = list(Binding=list(length=control_length)))

#readDensArgs <- list(fragment = fragLength, 
#                    bind = readDens_args_bind,
#                    minLength = readDens_args_minLength,
#                    maxLength = readDens_args_maxLength,
#                    meanLength = readDens_args_maxLength)

my_features <- list()

if (simulate_background) {
    for (chromosome in names(genome)) {
        my_features[[chromosome]] <- list()
        background_start <- 1
        background_length <- chromosome_lengths[[chromosome]] - background_start
        background_feature <- backgroundFeature(background_start, background_length)
        my_features[[chromosome]][[1]] <- background_feature
    }
} else {
    colnames(bed_file) <- c("chromosome", "start", "end") #, "name", "score")
    chromosome <- FALSE
    chromosome_nr <- 0
    for (peak_number in 1:nrow(bed_file)) {
      if (bed_file[peak_number, "chromosome"] != chromosome) {
        if (chromosome) {
          # we're at the end of a chromosome
          background_start <- previous_end
          background_length <- chromosome_lengths[[chromosome]] - background_start
          background_feature <- backgroundFeature(background_start, background_length)
          my_features[[chromosome_nr]][[feature_nr]] <- background_feature
        }
        # new chromosome
        chromosome <- bed_file[peak_number, "chromosome"]
        chromosome_nr <- chromosome_nr + 1
        previous_end <- 2
        feature_nr <- 1
        my_features[[chromosome_nr]] <- list()
      }
      peak_start <- bed_file[peak_number, "start"]
      peak_length <- bed_file[peak_number, "end"] - bed_file[peak_number, "start"]
      background_start <- previous_end - 1
      background_length <- peak_start - previous_end
      background_feature <- backgroundFeature(background_start, background_length)
      binding_feature <- bindingFeature(peak_start, peak_length)
      my_features[[chromosome_nr]][[feature_nr]] <- background_feature
      my_features[[chromosome_nr]][[feature_nr+1]] <- binding_feature
      feature_nr <- feature_nr + 2
      previous_end <- bed_file[peak_number, "end"]
    }
    background_start <- previous_end
    background_length <- chromosome_lengths[[chromosome]] - previous_end
    background_feature <- backgroundFeature(background_start, background_length)
    my_features[[chromosome_nr]][[feature_nr]] <- background_feature
}

class(my_features) <- c("ChIPseq", "SimulatedExperiment")
features <- my_features
save(features, file = paste(name, "features.rdata", sep = "_"))

simulated <- ChIPsim::simChIP(number_of_reads, genome, file = name,
                functions = myFunctions, control = parameters, load = TRUE)


