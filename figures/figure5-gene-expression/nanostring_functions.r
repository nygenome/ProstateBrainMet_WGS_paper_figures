#!/nfs/sw/R/R-3.6.1/bin/Rscript
## Common functions

extractGenes = function( X , geneNames) {
  ## this function is a private function of mergeCodesets. 
  if( !is.data.frame(X)) stop("No data frame provided")
  if( !is.character(geneNames) | !is.null(ncol(geneNames)) ) stop("geneNames must be a vector of character")
  geneData = X[ X[,1] %in% geneNames, ]
  if( nrow(geneData) > 0 ) return( X[ X[,1] %in% geneNames,  ]) else return(NULL)
}
mergeCodesets = function ( nanoData, geneNames ) {
  if( !is.list( nanoData )) stop( "nanoData must be a list of codesets")
  work.data = lapply( nanoData, extractGenes, geneNames )
  listOfValidData = which(!unlist(lapply( work.data, is.null )))
  work.xprl = work.data[[listOfValidData[1]]]
  # if( nrow(work.xprl) != length(geneNames) ) stop("List of genes do not match the nanoData")
  if( length( work.data ) > 1 ) {
    for( i in listOfValidData[2]:length(work.data)) {
      if( length(work.xprl$Gene)>0 & length(work.data[[i]][,1])>0 ) { ## only if there is a list of genes 
        if( length( work.data[[i]][,1] != nrow(work.xprl) ) ) { 
          warning("Different gene list... perform partial matching")
          sel.genes = intersect( work.data[[i]][,1], work.xprl$Gene)
          print(sel.genes)
        } else if( !(all( work.data[[i]][,1]==work.xprl$Gene) ) ) {
          warning("Different gene list... perform partial matching")
          sel.genes = work.xprl$Gene
        }
        if( ncol(work.data[[i]]) > 2 ) work.xprl = cbind( work.xprl[ work.xprl$Gene %in% sel.genes, ], work.data[[i]][ work.data[[i]]$Gene %in% sel.genes,-(1:2)])
      } else { stop ( "No genes available") }
    }
  } 
  return(work.xprl)
}

generateBoxPlots = function( nanoData, sel.genes, smpl.ann.red=NULL, pdfFileName = NULL, clrs = NULL, beeswarm = FALSE, pvalue=FALSE, comparisons=NULL, thresholdLine=NULL, smpl.names=NULL, x.angle=0 ) {
  if( is.null( smpl.ann.red ) ) stop("smpl.ann.red must be defined")
  if( is.null( smpl.ann.red$Class)) stop("smpl.ann.red$Class must be defined")
  if( pvalue & is.null(comparisons)) stop("if pvalue is true, comparisons should not be null. Comparison is a list of pairwise groups, eg. list( c('group1', 'group2'), c('group2','group3')) ")
  require(ggbeeswarm); require(ggpubr)
  if( !is.data.frame(nanoData)) {
    genes.selected = unique(unlist(lapply( nanoData, function( X ) { sel.genes[ sel.genes %in% X$Gene ] } )))
    work.xprl = mergeCodesets( nanoData, genes.selected)
  } else work.xprl = nanoData
  row.names( work.xprl ) = work.xprl$Gene
  work.xprl = work.xprl[,-(1:2)]
  
  plot.df = as.data.frame(t(log2(work.xprl)))
  plot.df$Class = smpl.ann.red$Class[ match(row.names(plot.df), smpl.ann.red$NanostringID )]
  plot.df$ID = rownames(plot.df)
  plot.df = plot.df[ !is.na(plot.df$Class), ]
  plot.m = melt(data = plot.df, measure.vars = sel.genes);head(plot.m)
  plot.m$highlight = FALSE
  if( !is.null(smpl.names) ) {
    plot.m$highlight[ plot.m$ID %in% smpl.names ]=TRUE
  }
  if(beeswarm) outliers = NA else outliers=NULL
  gv = ggplot( data = plot.m, mapping=aes(y=value, x=Class, fill=Class) ) + theme_bw(base_size = 20) + theme(panel.grid = element_blank()) + geom_boxplot(varwidth = TRUE, show.legend = F, outlier.shape = outliers )   # + geom_text( aes(1.5, 1, label="test"))
  if( !is.null( clrs )) gv = gv + scale_fill_manual(values=clrs)
  if( beeswarm ) gv = gv + geom_boxplot(varwidth = TRUE, show.legend = F ) + geom_beeswarm(priority = "ascending", cex=1.7, groupOnX = T, show.legend = F,  color="#980043")
  if( pvalue )  gv = gv  +  stat_compare_means(method = "wilcox.test", comparisons = comparisons)
  gv = gv +  scale_shape_identity() + geom_point(data = plot.m[plot.m$highlight, ], aes(x=Class, y=value, shape=23), fill="yellow", color="black", shape=23, size =2 )
  gv = gv +   facet_wrap(~variable) + xlab(label="") + theme(axis.text.y=element_text(size=10), axis.title = element_text(size=15), axis.text.x=element_text(size=12, angle=x.angle), strip.text=element_text(size=30)) + ylab(label="log2(normalized counts) ")
  if( !is.null(thresholdLine) ) {
    gv = gv + geom_hline(aes(yintercept=thresholdLine), colour="dark gray", linetype="dashed") 
  }
  print(gv)
  if( !is.null(pdfFileName) ) {
    if( file.exists(dirname(pdfFileName))) dev.print(pdf, file = pdfFileName ) else warning("The directory does not exists.", pdfFileName )
  }
  return(plot.m)
}

generateWaterfallPlots = function( nanoData, geneName,  sortByTarget = T, scaledByMean = T, smpl.names = NULL, whereTo=0, clrs = NULL, pdfFileName = NULL, ...) {
  if( is.null( smpl.ann.red ) ) stop("smpl.ann.red must be defined")
  if( is.null( smpl.ann.red$Class)) stop("smpl.ann.red$Class must be defined")
  if( is.null( geneName ) ) stop("geneName cannot be null")
  if( length( geneName ) > 1 ) stop("only 1 gene at the time")
  sel.genes=geneName
  if( !is.data.frame(nanoData)) {
    genes.selected = unique(unlist(lapply( nanoData, function( X ) { sel.genes[ sel.genes %in% X$Gene ] } )))
    work.xprl = mergeCodesets( nanoData, genes.selected)
  } else work.xprl = nanoData
  if( nrow( work.xprl ) == 0 ) { print(sprintf("Gene not found: (%s)", geneName)); return(-1) }
  row.names( work.xprl ) = work.xprl$Gene
  work.xprl = work.xprl[,-(1:2)]
  
  if( is.null(clrs) ) clrs = RColorBrewer::brewer.pal(4, "Dark2")
  plot.df = as.data.frame(t(log2(work.xprl)))
  g.mean = apply( plot.df, 2, mean )
  plot.df$Class = smpl.ann.red$Class[ match(row.names(plot.df), smpl.ann.red$NanostringID )]
  plot.df = plot.df[ !is.na(plot.df$Class), ]
  # plot.m = melt(data = plot.df,measure.vars = sel.genes);head(plot.m)
  plot.wf = NULL
  if( sortByTarget ) {
    for (l in 1:nlevels(plot.df$Class)) {
      tmp = droplevels( plot.df[ plot.df$Class==levels(plot.df$Class)[l], ])
      tmp = tmp[ sort.list( tmp[,1] ), ]
      if( is.null( plot.wf ) ) plot.wf = tmp else plot.wf = rbind( plot.wf, tmp )
    }
    if( scaledByMean )  plot.wf[,1] = plot.wf[,1] - g.mean
     
  } else {
    plot.wf = plot.df[ sort.list( plot.df[,1]), ]
  }
  if( scaledByMean ) {
    y.label = "log2(normalized counts) - mean"
  } else {
    y.label = "log2( normalized counts )"
  }
  
  gv = ggplot( data = plot.wf, mapping=aes_string(y=paste("`",geneName,"`",sep=""), x=1:nrow(plot.wf), fill="Class") ) + theme_bw(base_size = 20)+ geom_bar( stat="identity", show.legend = T)   
  
  if( !is.null( clrs )) gv = gv + scale_fill_manual(values=clrs)
  if( !is.null( smpl.names ) ) {
    idx.samples = NULL
    for( s.name in smpl.names) {
      ##idx.smpl <- grep( smpl.name, names(gene.xpr)[idx.s] )
      idx.smpl <- grep( s.name, rownames(plot.wf) )
      if(length(idx.smpl)>0) {
        idx.samples = c(idx.samples, idx.smpl)
      } else {
        warning( sprintf("Can't find %s", s.name))
      }
    }
    if( length(idx.samples ) > 0) {
      gv = gv + scale_shape_identity() + geom_point( aes(x=idx.samples, y=whereTo, shape=23), color="red", fill="red") 
    }
  }
  
  gv = gv +    xlab(label="") + theme(axis.text.y=element_text(size=10), axis.title = element_text(size=15), axis.text.x=element_text(size=12), strip.text=element_text(size=30)) + ylab(label= y.label)
  if( hasArg(title) ) gv = gv + ggtitle( geneName )
  print(gv)
  if( !is.null(pdfFileName) ) {
    if( file.exists(dirname(pdfFileName))) dev.print(pdf, file = pdfFileName ) else warning("The directory does not exists.", pdfFileName )
  } 
  return( plot.wf )
}

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

runCommand = function ( cmd , local=FALSE, silent=T) {
  if( !silent ) print(cmd)
  return( system(cmd, intern=T) )
}


