### utility functions for genotype matrices

## given a vector of values, return the <n> most frequent ones
## NB: for integer or character data only; floating-point will cause problems
topn <- function(x, n = 1, ...) {
	
	x <- as.vector(x)
	tbl <- sort(table(x, useNA = "ifany"), decreasing = TRUE)
	top <- names(tbl)[1]
	
	if (is.numeric(x))
		top <- as.numeric(top)
	
	return(top)
	
}

## recode a matrix (markers x samples) of genotypes as 0/1/2/NA
recode.genotypes <- function(gty, calls = c("A","C","G","T","H","N"), ...) {
	
	if (!is.matrix(gty))
		gty <- as.matrix(gty)
	
	recoded <- matrix(0, nrow = nrow(gty), ncol = ncol(gty))
	colnames(recoded) <- colnames(gty)
	
	for (i in 1:nrow(gty)) {
		
		row <- gty[ i, ]
		tbl <- sort( table( factor(row[ (row != "H" & row != "N") ], levels = calls) ), decreasing = TRUE )
		maj <- names(tbl)[1]
		new.row <- rep(2, length(row))
		new.row[ row == "N" ] <- NA
		new.row[ row == "H" ] <- 1
		new.row[ row == maj ] <- 0
		recoded[ i, ] <- new.row
		
	}
	
	return(recoded)
	
}

## dump a genotype matrix (markers x samples) to a file suitable for R/qtl's "csvs" format, returning the filename
## <gty> = genotype matrix; <map> = dataframe(-like) object with at least the following columns: marker name, chromosome
## NB: map and genotype matrix are assumed to have same sort order
as.rqtl.geno <- function(gty, map, file = tempfile(), ...) {
	
	if (nrow(gty) != nrow(map)) {
		stop("Map and genotype matrix don't match.")
	}
	if (is.null(colnames(gty))) {
		stop("Genotype matrix must have column names.")
	}
	if (mode(gty) != "character") {
		stop("Genotype matrix must be of mode character.")
	}
	
	ids <- colnames(gty)
	gty <- cbind( as.character(map[,2]), gty )
	gty <- cbind( as.character(map[,1]), gty )
	tgty <- t(gty)
	tgty <- cbind( c("id","",ids), tgty )
	
	write.table(tgty, file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
	return(file)
	
}

## convert a genotype matrix (markers x samples) to a format compatible with PLINK's ped format
## ped format: family, ID, dad, mom, sex, phenotype, {genotypes: allele1 allele2}_1, ... {genotypes}_k
as.plink.geno <- function(gty, recode = FALSE, ...) {
	
	## check if genotype matrix appears to have been recoded
	.gty <- gty
	if (mode(.gty) != "numeric" | recode) {
		gty <- recode.genotypes(.gty, ...)
	}
	gty <- gty + 1 # since 0 = missing data in PLINK
	gty[ which(is.na(gty)) ] <- 0
	
	new.gty <- matrix(0, nrow = ncol(gty), ncol = 2*nrow(gty))
	rownames(new.gty) <- colnames(gty)
	for (i in 0:(nrow(gty)-1)) {
		this.mk <- gty[ (i+1), ]
		new.mk <- vapply(this.mk, function(x) if(x == 3) c(1,2) else rep(x,2), numeric(2))
		new.gty[ ,(2*i+1):(2*i+2) ] <- t(new.mk)
	}
	
	return(new.gty)
	
}