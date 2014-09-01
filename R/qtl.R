library(GenomicRanges)
library(plyr)

mapping.intervals <- function(haps, ...) {
	
	## get disjoint intervals covering whole genome, within each of which no individual is recombinant
	
	rez <- ldply(haps, function(h) ldply(h, as.data.frame),
							 .progress = "text")
	gr <- makeGRangesFromDataFrame(rez[ ,c("seqnames","start","end") ], ignore.strand = TRUE)
	bins <- disjoin(gr)
	return(bins)
	
}

calc.hap.dosage <- function(haps, bins, alleles = NULL, mode = "additive", ...) {
	
	## returns 3d array: samples x bins x strains
	
	if (!(mode == "additive") | is.null(alleles))
		stop()
	
	.haps <- haps
	
	haps <- flatten.haplotypes(haps, "parent")
	# print(haps)
	rez <- laply(haps, function(h) {
		olap <- as.data.frame(findOverlaps(bins, h))
		olap$strain <- h$strain[ olap$subjectHits ]
		olap$strain <- factor(olap$strain, levels = alleles)
		olap$bin <- factor(olap$queryHits, levels = 1:length(bins))
		xtabs(~ bin + strain, olap)
	}, .progress = "text")
	
	attr(rez, "intervals") <- bins
	return(rez)
	
}

hk.regress <- function(phe, gty, ...) {
	
	## given phenotype vector and haplotype-dosage matrix, do HK-style regressino
	
	d <- data.frame(y = phe, as.data.frame(gty))
	m0 <- lm(y ~ 1, d)
	m1 <- lm(y ~ ., d)
	# print(summary(m1))
	n <- length(resid(m0))
	rss0 <- sum(resid(m0)^2)
	rss <- sum(resid(m1)^2)
	lod <- as.numeric(n/2)*log10(rss0/rss)
	rez <- c(coef(m1), lod)
	names(rez)[ length(rez) ] <- "LOD"
	
	return(list(lod = lod, coef = coef(m1),
							n = n, df.resid = m1$df.resid))
	
}

map.single.locus <- function(phe, gty, ...) {
	
	## run QTL scan on haplotype dosage object returned by calc.hap.dosage()
	
	lod <- aaply(gty, 2, function(x) {
		rez <- hk.regress(phe, x)
		return(rez$lod)
	}, .progress = "text")
	
	bins <- attr(gty, "intervals")
	rez <- data.frame(chr = seqnames(bins), pos = mid(ranges(bins)), lod = lod)
	return(rez)
	
}