#
# oncoPrint : plot a genotype
#
# This function sorts the matrix for better visualization of mutual exclusivity across genes
exclusivity.sort <- function(M) {
	geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
	scoreCol <- function(x) {
		score <- 0;
		for(i in 1:length(x)) {
			if(x[i]) {
				score <- score + 2^(length(x)-i);
			}
		}
		return(score);
	}
	scores <- apply(M[geneOrder, ], 2, scoreCol);
	sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
	
	res = list()
	res$geneOrder = geneOrder
	res$sampleOrder = sampleOrder
	res$M = M[geneOrder, sampleOrder]
	
	return(res);
}
 
#
# This is the plotting function
#
oncoprint <- function(x, excl.sort=TRUE, col.cluster=FALSE, row.cluster=FALSE, device.new=FALSE, file=NA, ann.stage=TRUE, ann.score=TRUE, stage.color='YlOrRd', score.color = 'Purples',  null.color='darkgray', border.color='white', font.size=7, font.column = 3, title= paste('Genotypes'), ...) 
{
	if (!require('pheatmap')) {
    	install.packages('pheatmap', dependencies = TRUE)
    	library(pheatmap)
  	}

	if (!require('RColorBrewer')) {
    	install.packages('RColorBrewer', dependencies = TRUE)
    	library(RColorBrewer)
  	}

	is.compliant(x, 'oncoprint', stage=ann.stage)
	
	# We reverse the heatmap under the assumption that ncol(data) << nrow(data)
  	data = t(x$genotypes)
  	nc = ncol(data) 
  	nr = nrow(data)
  		
	# Sort data, if required. 
	if(excl.sort) {
		sorted.data = exclusivity.sort(data)
		data = sorted.data$M	
	}	
	cn = colnames(data)
  	rn = rownames(data)
	
	# Heatmap score annotation: total 1s per sample
	nmut = colSums(data)
	if(ann.score == TRUE && ann.stage == FALSE) annotation = data.frame(score=nmut)
	if(ann.score == FALSE && ann.stage == TRUE) annotation = data.frame(stage=x$stages[cn, ])
	if(ann.score == TRUE && ann.stage == TRUE)  annotation = data.frame(stage=x$stages[cn, ], score=nmut)
	
	if(ann.score == TRUE || ann.stage == TRUE) rownames(annotation) = cn
		
	# annotation colors
	score.gradient = (colorRampPalette(brewer.pal(6, score.color))) (max(nmut))
	annotation_colors = list(score = score.gradient)
	if(ann.stage == TRUE){ 
		num.stages = length(unlist(unique(x$stage)))
		stage.color.attr = brewer.pal(n=num.stages, name=stage.color)		
		annotation_colors = list(stage=stage.color.attr, score=score.gradient)
	}	
	
	# Display also event frequency, which gets computed now
	genes.freq = rowSums(data)/nc
	
	# Augment information to make type-dependent colored plots
	for(i in 1:nr)
	{
		a.type = x$annotations[rownames(data)[i],]
		idx.type = which(rownames(x$types) == a.type[1])
		idx.samples = names(which(data[i,]==1))
	    data[i,idx.samples] = idx.type
	}

	# Create new device output, if required
	if(device.new == TRUE) dev.new(width=ncol(data) * .3, height=nrow(data) * .7)

	# Map gradient
	map.gradient = c(null.color, x$types[,1])
	 
   	# Augment gene names with frequencies and prepare legend labels
	gene.names = x$annotations[rownames(data),2]
	rownames(data) = paste(round(100 * genes.freq, 2) ,'% ', gene.names, sep='')
	legend.labels = c('No alt.', unique(x$annotations[rn,1]))

	# Augment title
	title = paste(title, ' (', nc,' samples,  ', nr, ' events)', sep='')
	
	# Pheatmap
	if(ann.score == TRUE || ann.stage == TRUE)  
	 	pheatmap(data, 
		 	scale = "none", 
			col = map.gradient, 
			cluster_cols = col.cluster,
			cluster_rows = row.cluster,
			main= title,
			fontsize= font.size,
			fontsize_col= font.column,
			annotation = annotation,
			annotation_colors = annotation_colors,	
			border_color = border.color,
			border=T,
			margins=c(10,10),
			cellwidth = 6, 
			cellheigth = 20,
			legend=T,
			legend_breaks = c(0:max(data)),
			legend_labels = legend.labels,
			drop_levels=T,
			...
		)
	else
		pheatmap(data, 
			 	scale = "none", 
				col = map.gradient, 
				cluster_cols = col.cluster,
				cluster_rows = row.cluster,
				main= title,
				fontsize= font.size,
				fontsize_col= font.column,
				border_color= border.color,
				border=T,
				margins=c(10,10),
				cellwidth = 6, 
				legend=T,
				legend_breaks = c(0:max(data)),
				legend_labels = legend.labels,
				...
			)
			
	# save to file
	if(!is.na(file)) 
	{
		current.device <- dev.cur()
        dev.off(dev.copy(device = pdf, file = file))
        dev.set(current.device)
        cat(paste('Oncoprint written to file:', file,'\n'))
	}	
}