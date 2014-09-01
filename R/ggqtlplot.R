library(ggplot2)
library(plyr)

ggqtlplot <- function(x, chr = NULL, space = 5, band.cols = c(NA, "grey90"), alphas = c(0.1, 0.05, 0.01), thresh.lty = c("dotdash", "dashed","dotted"), perms = NULL,
											lod.col = "lod", title = NULL, marker.rug = TRUE, scale = setNames(1e6, "Mbp"), superpose.on = NULL, ...) {
	
	# require(qtl)
	# stopifnot("scanone" %in% class(x))
	
	if (is.null(chr)) {
		chr.bounds <- transform(ddply(x, .(chr), summarize, max = max(pos), min = min(pos)),
														offset = cumsum(c(0, max[ -length(max) ] + space)))
		chr.bounds <- ddply(chr.bounds, .(chr), transform, mid = (offset + (offset+max))/2,
												right = max(offset + max), fill = factor((as.integer(chr)+1) %% 2))
		x <- merge(chr.bounds, x)
		x <- transform(x, xpos = pos + offset)
		# print(head(x))
		# print(chr.bounds)
	}
	else {
		if (is(x, "cross")) {
			x <- subset(x, chr = chr) # careful: this is qtl::subset.cross(), not base::subset()
		}
		else {
			this.chr <- chr
			x <- subset.data.frame(x, chr == this.chr)
		}
		x$xpos <- x$pos/scale
	}
	if (!is.null(superpose.on)) {
		x$superpose <- x[,superpose.on]
	}
	# colnames(x)[ which(colnames(x) == lod.col) ] <- "lod"
	x$lod <- x[ ,lod.col ]
	
	## base plot
	p <- ggplot(x)
	
	## layer 1: banded colors and axes
	if (is.null(chr))
		## if specific chromosome not specified, do banded colors and set x-axis ticks to middle of chomosomes
		p <- p +	geom_rect(data = chr.bounds, aes(xmin = offset, xmax = right, ymin = -Inf, ymax = Inf, fill = fill)) +
							scale_fill_manual(values = band.cols) + scale_x_continuous(breaks = chr.bounds$mid, labels = gsub("^chr","", chr.bounds$chr)) +
							xlab("\nchromosome")
	else
		## no specific chromosome not specified; x-axis is position in Mbp
		p <- p + xlab(paste0("\nposition (", names(scale)[1], ")"))
	
	## layer 2: LOD curves
	if (!is.null(superpose.on))
		p <- p +	geom_line(aes(x = xpos, y = lod, group = factor(chr):factor(superpose), colour = superpose))
	else
		p <- p +	geom_line(aes(x = xpos, y = lod, group = chr))
						
	## layer 3: guides, y-axis label, theme...
	p <- p +	guides(fill = FALSE) + ylab("LOD score\n") +
						theme_bw()
	
	## option to add title
	if (!is.null(title))
		p <- p + ggtitle(title) + theme(plot.title = element_text(face = "bold"))
	
	## option to add markers as hashes ('rug') in bottom of plot
	if (marker.rug)
		p <- p + geom_point(aes(x = xpos, y = -0.1), pch = "|")
		
	lod.threshold.line <- function(z, thresh, lty, type = c("A","X"), ...) {
		d <- NULL
		if (type == "A") {
			A <- subset(z, !(chr %in% c("X","Y","M")))
			d <- data.frame( left = min(A$offset), right = max(A$right) )
		} else if (type == "X") {
			d <- data.frame(subset(z, chr == "X")[ 1,c("offset","right") ])
		}
		colnames(d) <- c("left","right")
		d$thresh <- thresh
		geom_segment(data = d, aes(x = left, xend = right,
															 y = thresh, yend = thresh), lty = lty)
	}
	
	if (!is.null(perms)) {
		if ("scanoneperm" %in% class(perms)) {
			thresh <- summary(perms, alpha = alphas)
			for (i in 1:nrow(thresh$A)) {
				p <- p + lod.threshold.line(chr.bounds, thresh$A[i], thresh.lty[i], type = "A") +
					lod.threshold.line(chr.bounds, thresh$X[i], thresh.lty[i], type = "X")
			}
		}
	}
	
	return(p)
	
}