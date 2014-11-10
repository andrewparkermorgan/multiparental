## utilities for converting R/qtl genotypes to haplotype objects

library(plyr)
library(qtl)
library(GenomicRanges)

## find markers with spurious evidence for double-recombination event in >X samples
find.suspicious.markers <- function(x, min.samples = 1, ...) {
	
	lapply(x$geno, function(c) {
		mk <- apply(c$data, 1, function(z) {
			rl <- Rle(z)
			start(rl)[ runLength(rl) == 1 & !is.na(runValue(rl)) ]
		})
		mnames <- names(c$map)[ unlist(mk) ]
		tbl <- table(mnames)
		names(tbl)[ tbl >= min.samples ]
	})
	
}

## convert observed genotype data (NOT genoprobs) to haplotype blocks
## currently only works for backcross; else phasing is an issue
rqtl.to.haplotypes <- function(x, seqlengths = multiparental:::seqlengths, phase = "maternal", ...) {
	
	if (length(attr(x, "alleles")) > 2)
		warning("This function only designed to work with a backcross; this might not be a backcross.")
	
	## loop on chromosomes
	rez <- ldply(names(x$geno), function(c) {
		
		## loop on individuals
		haps <- adply(x$geno[[c]]$data, 1, function(z) {
			mnames <- names(z)
			rl <- Rle(z)
			istart <- start(rl)
			iend <- end(rl)
			starts <- mk[ mnames[istart],"start" ]
			ends <- mk[ mnames[iend],"start" ]
			alleles <- attr(x, "alleles")[ runValue(rl) ]
			data.frame(seqnames = c, start = starts, end = ends, origin = phase, strain = alleles)
		})
		colnames(haps)[1] <- "id"
		return(haps)
		
	}, .progress = "text")
	
	roll.haplotypes(rez, seqlengths = seqlengths, ...)
	
}

## convert observed genotype data (NOT genoprobs) to recombination events
## currently only works for backcross; else phasing is an issue
rqtl.to.recombinations <- function(x, seqlengths = multiparental:::seqlengths, phase = "maternal", ...) {
	
	if (length(attr(x, "alleles")) > 2)
		warning("This function only designed to work with a backcross; this might not be a backcross.")
	
	## loop on chromosomes
	rez <- ldply(names(x$geno), function(c) {
		
		## loop on individuals
		recombs <- adply(x$geno[[c]]$data, 1, function(z) {
			mnames <- names(z)
			rl <- Rle(z)
			nblocks <- length(runValue(rl))
			if (nblocks == 1) {
				return(NULL)
			}
			else {
				iright <- start(rl)[-1]
				ileft <- end(rl)[ -nblocks ]
				lefts <- mk[ mnames[ileft],"start" ]
				rights <- mk[ mnames[iright],"start" ]
				froms <- attr(x, "alleles")[ runValue(rl)[ -nblocks ] ]
				tos <- attr(x, "alleles")[ runValue(rl)[-1] ]
				return( data.frame(seqnames = c, start = lefts, end = rights, origin = phase,
													 from = froms, to = tos) )
			}
		})
		colnames(recombs)[1] <- "id"
		return(recombs)
		
	}, .progress = "text")
	
	roll.recombinations(rez, seqlengths = seqlengths, ...)
	
}