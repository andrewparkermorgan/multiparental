library(GenomicRanges)
library(plyr)

forbidden <- c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
							 "isCircular", "start", "end", "width", "element")

roll.haplotypes <- function(haps, seqlengths = seqlengths, seal = FALSE, ...) {
	
	dlply(haps, .(id), function(d) {
		rez <- dlply(d, .(origin), function(p) {
			seg <- with(p, GRanges(seqnames = Rle(seqnames), seqlengths = seqlengths,
														 ranges = IRanges(start = start, end = end),
														 strain = strain, id = id))
			seg <- sort(seg)
			if (seal)
				return( seal.segments(seg) )
			else
				return(seg)
		})
		class(rez) <- c("haplotypes", class(rez))
		return(rez)
	}, .progress = "text")
	
}

roll.recombinations <- function(haps, seqlengths = seqlengths, ...) {
	
	dlply(haps, .(id), function(d) {
		rez <- dlply(d, .(origin), function(p) {
			gr <- with(p, GRanges(seqnames = Rle(seqnames), seqlengths = seqlengths,
														ranges = IRanges(start = start, end = end),
														from = from, to = to, id = id))
			values(gr) <- cbind(values(gr), p[ ,setdiff(colnames(p), c(forbidden,colnames(values(gr)))) ])
			return(gr)
		})
		# class(rez) <- c("haplotypes", class(rez))
		return(rez)
	})
	
}

invert.recombinations <- function(recombs, seal = FALSE, ...) {
	
	lapply(recombs, function(r.all) {
		segments <- GRanges()
		for (c in unique(runValue(seqnames(r.all)))) {
			r <- r.all[ seqnames(r.all) == c ]
			if (seal) {
				# 'seal' haplotype segments together -- eg. get rid of uncertainty around recomb events
				g <- GRanges(seqnames = c, seqinfo = seqinfo(r),
										 ranges = IRanges(start = c(1, pmin( mid(ranges(r))+1, seqlengths(r)[c] )),
										 								 end = c(mid(ranges(r)), seqlengths(r)[c])),
										 strain = c( as.character(values(r)$from), rev(as.character(values(r)$to))[1] )
				)
				segments <- c(segments, g)
			}
			else {
				if (length(r)) {
					if (length(r) > 1) {
						i <- 1:(length(r)-1)
						# un-overlap possibly overlapping recomb events
						end(r)[i] <- pmin( end(r)[i], start(r)[i+1]-1 )
						start(r)[i+1] <- pmax( end(r)[i]+2, start(r)[i+1] )
					}
					r <- sort(r)
					g <- gaps(r)
					g <- g[ strand(g) == "*" & (as(seqnames(g), "character") %in% as(seqnames(r), "character")) ]
					values(g)$strain <- c( as.character(values(r)$from), rev(as.character(values(r)$to))[1] )
					segments <- c(segments, g)
				}
			}
		}
		# print(segments)
		return(segments)
	})
	
}

# overlap.with.meta <- function(query, subject, metavars = colnames(values(query)), ...) {
# 	
# 	flag <- (query %over% subject)
# 	meta.flag <- matrix(ncol = 0, nrow = length(query))
# 	
# 	if (length(subject) > 1 & !(length(subject) == length(query))) {
# 		warning("Result will probably not be well-defined due to shape of query vs subject ranges.")
# 	}
# 	
# 	for (i in colnames(values(query))) {
# 		if (i %in% colnames(values(subject))) {
# 			# hits.this.subject <- sapply(1:length(subject), function(j) flag & (as.character(values(query)[,i]) == as.character(values(subject)[j,i])) )
# 			# print(hits.this.subject)
# 			hits.this.subject <- flag & (values(query)[,i] == values(subject)[,i])
# 			meta.flag <- cbind(meta.flag, hits.this.subject)
# 		}
# 	}
# 	
# 	return(meta.flag)
# 	
# }
# 
# grep.haplotypes <- function(haps, target, strain.col = "strain", ...) {
# 	
# 	## how many haplotype segments in <haps> overlap <target>? 
# 	
# 	stopifnot(length(target) > 0)
# 	
# 	llply(haps, function(h) {
# 		rez <- lapply(h, function(hh) {
# 			sapply(1:length(target), function(t) overlap.with.meta(hh, target[t], metavars = strain.col))
# 		})
# 		return( sum(sapply(rez, any)) )
# 	}, .progress = "text")
# 	
# }
# 
# grep.recombinations <- function(recombs, target, from = NULL, to = NULL, ...) {
# 	
# 	stopifnot(length(target) > 0)
# 	
# 	llply(recombs, function(h) {
# 		sapply(h, function(hh) {
# 			i <- rep(TRUE, length(hh))
# 			if (!is.null(from))
# 				i <- i & (values(hh)$from == from)
# 			if (!is.null(to))
# 				i <- i | (values(hh)$to == to)
# 			countOverlaps(target, hh[i])
# 		})
# 	}, .progress = "text")
# 	
# }

overlap.with.meta <- function(haps, target, haps.flat = NULL, meta.cols = NULL, by = .(id), ...) {
	
	## given a list-of-lists of rolled haplotypes, find those which overlap <target> AND have same metadata
	## (if a flattened haplotypes object is already available, allow it to be passed in to avoid re-flattening)
	
	all.gr <- GRanges()
	if (!is.null(haps.flat))
		all.gr <- haps.flat
	else
		all.gr <- flatten.haplotypes(haps, "all")
	
	olap <- all.gr %over% target
	if (is.null(meta.cols))
		meta.cols <- intersect( colnames(values(target)), colnames(values(all.gr)) )
	for (col in meta.cols) {
		olap <- olap & (values(all.gr)[ ,col ] == values(target)[ ,col ])
	}
	all.gr$olap <- olap
	all.df <- as.data.frame(all.gr)
	rez <- daply(all.df, by, summarize, sum(as.numeric(olap)))
	
	return(rez)
	
}

grep.haplotypes <- function(haps, target, meta.cols = "strain", ...) {
	
	## how many haplotype segments in <haps> overlap <target>? 
	
	stopifnot(length(target) > 0)
	overlap.with.meta(haps, target, meta.cols = meta.cols, ...)
	
}

grep.recombinations <- function(recombs, target, meta.cols = c("to","from"), ...) {
	
	## how many recombination events in <recombs> overlap <target>? 
	
	stopifnot(length(target) > 0)
	overlap.with.meta(recombs, target, meta.cols = meta.cols, ...)
	
}

find.recombinants <- function(...) {
	
	rez <- grep.recombinations(...)
	sapply(rez, sum) > 0
	
}

seal.segments <- function(gr, ...) {
	
	## "seal" haplotype segments by joining them at midpoint of uncertainty interval

	df <- as.data.frame(sort(gr))
	rez <- ddply(df, .(seqnames), function(d) {
		mids <- floor((d$start[-1] + d$end[ -nrow(d) ])/2)
		d$start[1] <- 1
		d$end[-nrow(d)] <- mids
		d$start[-1] <- mids + 1
		d$end[ nrow(d) ] <- seqlengths(gr)[ d$seqnames[1] ]
		d$width <- NULL # not an allowed column name to makeGRangesFromDataFrame()
		return(d)
	})
	
	makeGRangesFromDataFrame(rez, keep.extra.columns = TRUE, ignore.strand = TRUE)
	
}

flatten.haplotypes <- function(haps, by = c("parent","all"), ...) {
	
	## flatten the nested-list haplotype structure either into a list of GRanges() (by = "parent") or a single GRanges(by = "all")
	
	if (by == "all") {
		rez <- ldply(haps, function(h) ldply(h, as.data.frame),
								 .progress = "text")
		if ("width" %in% colnames(rez))
			rez$width <- NULL
		makeGRangesFromDataFrame(rez, keep.extra.columns = TRUE, ignore.strand = TRUE)
	}
	else if (by == "parent") {
		rez <- llply(haps, function(h) {
			hh <- ldply(h, as.data.frame)
			if ("width" %in% colnames(hh))
				hh$width <- NULL
			makeGRangesFromDataFrame(hh, keep.extra.columns = TRUE, ignore.strand = TRUE)
		}, .progress = "text")
		return(rez)
	}
	
}