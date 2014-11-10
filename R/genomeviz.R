library(GenomicRanges)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)
library(plyr)

CC.STRAINS <- setNames( toupper(letters[1:8]), c("AJ","B6","S129","NOD","NZO","CAST","PWK","WSB") )
CC.COLORS <- c(A="#DAA520",D="#1010F0",B="#404040",C="#F08080",H="#9000E0",G="#F00000",E="#00A0F0",F="#00A000")[ CC.STRAINS ]

cc.strains <- setNames( names(CC.STRAINS), CC.STRAINS )
cc.colors <- setNames( CC.COLORS[ CC.STRAINS ], names(CC.STRAINS) )

trans_genome <- function(scale = 1e6, ...) {
	
	trans_new("genome coordinates",
						transform = function(x) x/scale,
						inverse = function(X) X*scale )
	
}

scale_x_genome <- function(scale = 1e6, unit = "Mbp", ...) {
	
	gtrans <- trans_genome(scale, unit)
	scale_x_continuous(trans = gtrans, labels = function(x) gtrans$transform(x), ...)
	
}

.scale_CC <- function(fn, colors = cc.colors, ...) {
	fn(values = colors, ...)
}

scale_fill_CC <- function(...) {
	.scale_CC(scale_fill_manual, ...)
}

scale_colour_CC <- function(...) {
	.scale_CC(scale_colour_manual, ...)
}

ggtitle.n <- function(title = NULL, n = 0, units = "", ...) {
	ggtitle(paste0(title, "\n", "(n = ", n, " ", units, ")\n"))
}

theme_genome <- function(...) {
	theme_bw() %+replace%
		theme(plot.title = element_text(face = "bold"))
}

gggenome <- function(...) {
	ggplot(...) +
		scale_x_genome() +
		xlab("\nposition (Mbp)") +
		theme_genome()
}

plot.haplotypes <- function(haps, flatten = TRUE, fill = "strain", omit.chrs = c("chrY","chrM"), ...) {
	
	## flatten list-of-lists, if necessary
	haps.df <- data.frame()
	if (flatten | !is.data.frame(haps))
		haps.df <- ldply(haps, as.data.frame)
	else
		haps.df <- haps
	
	## keep only desired chromosomes (1-19,X by default)
	haps.df <- subset(haps.df, !(seqnames %in% omit.chrs))
	chroms <- setdiff( levels(haps.df$seqnames), omit.chrs )
	haps.df$seqnames <- factor( haps.df$seqnames, levels = chroms )
	
	haps.df$seqnames <- with(haps.df, factor(as.character(seqnames), levels = rev(levels(seqnames))) )
	haps.df$updown <- ifelse(haps.df$origin == levels(haps.df$origin)[1], 1, -1)
	haps.df$ypos <- as.numeric(haps.df$seqnames)
	haps.df <- transform(haps.df, ymin = 0.4*updown+ypos, ymax = 0.05*updown+ypos)
	haps.df$fill <- haps.df[ ,fill ]
	
	p <- gggenome(haps.df) +
		geom_rect(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = fill)) +
		scale_y_continuous(breaks = 1:length(chroms), labels = rev(chroms)) +
		theme(axis.title.y = element_blank(), plot.title = element_text(face = "bold"))
	
	if ("id" %in% colnames(haps.df))
		p <- p + ggtitle(paste0(haps.df$id[1], "\n"))
	
	return(p)
	
}

plot.stacked.haps <- function(.haps, chroms = "chr2", which = "both",
															sort.order = setNames(1:length(.haps), names(.haps)),
															strain.col = "strain", label = TRUE, xlim = NULL,
															spacing = 0.05, gap = 0.1, group.by = NULL, ...) {
	
	## copy original input
	haps <- .haps
	haps.df <- data.frame()
	
	## extract desired haplotype for each sample
	if (length(which) == length(haps)) {
		## given a vector <which> parallel to <haps>, use its i-th element to subset the i-th segment set
		haps <- lapply(1:length(which), function(i) {
			rez <- haps[[i]][ which[i] ]
			return(rez)
			# if (!is.null(rez))
				# rez[ rez %over% x.target ]
			# else
				# return(x.target)
		})
		names(haps) <- names(.haps)
		## roll up haplotypes into dataframe for plotting, adding sort order
		haps.df <- ldply(haps, ldply, as.data.frame)
		colnames(haps.df)[1] <- "origin"
	}
	else if (which[1] != "both") {
		## subset using first element of <which>
		haps <- lapply(haps, "[[", which[1])
		## roll up haplotypes into dataframe for plotting, adding sort order
		haps.df <- ldply(haps, ldply, as.data.frame)
		colnames(haps.df)[1] <- "origin"
	}
	else {
		## roll up haplotypes into dataframe for plotting, adding sort order
		haps.df <- ldply(haps, ldply, as.data.frame)
	}
	
	shrink <- ifelse(which[1] == "both", 0.5, 1)
	haps.df$ypos <- sort.order[ as.character(haps.df$id) ]
	if (which[1] == "both") {
		haps.df$updown <- ifelse(haps.df$origin == "maternal", 1, -1)
	}
	else {
		haps.df$ypos <- haps.df$ypos - 0.5
		haps.df$updown <- 1
	}
	haps.df$ymin <- with(haps.df, ypos + pmin(0, updown*(shrink-gap/2)) + ifelse(updown == 1, spacing, 0))
	haps.df$ymax <- with(haps.df, ypos + pmax(0, updown*(shrink-gap/2)) - ifelse(updown == -1, spacing, 0))
	
	## construct stacked-haplotype plot
	if (!is.null(chroms))
		haps.df <- subset(haps.df, seqnames %in% chroms)
	rez <- ggplot(haps.df) +
		geom_rect(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = strain)) +
		xlab("\nposition (Mbp)") +
		scale_x_genome() +
		scale_y_continuous(breaks = sort.order, limits = c(0,max(sort.order)+1)) +
		theme_bw() +
		theme(axis.title.y = element_blank())
	if (!label) {
		rez <- rez + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
	}
	
	## tack on the sorting order, for use outside this function
	attr(rez, "sort.order") <- sort.order
	
	return(rez)
	
}

ggplot.GRanges <- function(x, spacing = 0.1, ...) {
	
	df <- as.data.frame(sort(x))
	df$.idx <- 1:nrow(df)-1
	p <- gggenome(df, aes(x = start, xend = end,
												xmin = start, xmax = end,
												y = .idx, yend = .idx,
												ymin = (.idx - 0.5),
												ymax = (.idx + 0.5))
								)
	p <- p + theme(axis.ticks.y = element_blank(),
								 axis.text.y = element_blank(),
								 axis.title.y = element_blank())
	
	return(p)
	
}