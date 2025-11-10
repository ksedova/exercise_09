rm(list=ls())

setwd('D:/university/Amasters/firstsemester/PRG/exercise_09')

# *De Novo* Genome Assembly

## The Greedy Shortest Common Superstring
### Task
#* In R, implement a function `GreedySuperstring()` according to the pseudocode.

#* Input:
#  * `S` A `DNAStringSet` object of strings (reads).

#* Output:
#  * `S` A `DNAStringSet` object of the shortest common superstring (contig).

#> **Hint:** 
#  > Create also functions:
#  > * `Overlap()` to calculate overlap between two sequences.
#> * `OverlapMatrix()` to create a matrix of overlaps among all sequences in `S`.

#```
#GreedySuperstring(S)
#1   while length of S > 1
#2     overlapMat <- OverlapMatrix(S)
#3     if max(overlapMat) = 0
#4       return S
#5     else
#  6       seq1, seq2 ← Two sequences from S with the longest overlap
#7       Merge seq1 and seq2 and add the new sequence to S
#8       Remove seq1 and seq2 from S
#9   return S
#```
library(Biostrings)

Overlap <- function(a, b){
  max_overlap <- min(nchar(a), nchar(b))
  for (overlap in seq(max_overlap, 1)) {
    if (substr(a, nchar(a) - overlap + 1, nchar(a)) == substr(b, 1, overlap)) {
      return(overlap)
    }
  }
  return(0)
}

Overlap("GATTACA", "TACAGA")

OverlapMatrix <- function(S){
  n <- length(S)
  mat <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j) {
        mat[i, j] <- Overlap(as.character(S[i]), as.character(S[j]))
      }
    }
  }
  return(mat)
}

GreedySuperstring <- function(s){
  while (length(s) > 1){
    overlapMat <- OverlapMatrix(s)
    if (max(overlapMat) == 0){
      return(s)
    } else {
      idx <- which(overlapMat == max(overlapMat), arr.ind = TRUE)[1, ]
      i <- idx[1]
      j <- idx[2]
      
      seq1 <- as.character(s[i])
      seq2 <- as.character(s[j])
      overlap <- overlapMat[i, j]
      
      #merge sequences
      merged_seq <- paste0(substr(seq1, 1, nchar(seq1) - overlap), seq2)
      
      # update S
      s <- s[-c(i, j)]
      s <- append(s, DNAStringSet(merged_seq))
    }
  }
  return(s)
}



reads <- DNAStringSet(c("ATTAC", "TACAG", "AGATT"))
reads1 <- DNAStringSet(c('CATGC', 'CTAAGT', 'GCTA', 'TTCA', 'ATGCATC'))
GreedySuperstring(reads1)

  