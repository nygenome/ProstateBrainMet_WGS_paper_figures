#!/nfs/sw/R/R-3.6.1/bin/Rscript
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper

################################################################# /COPYRIGHT ###
################################################################################
## Plot jabba track and/or coverage track and/or allele track by a bed file
libs = c('optparse', 'GenomicRanges', 'gGnome', 'gTrack', 'gUtils')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))



alpha = function(col, alpha) {
  col.rgb = col2rgb(col)
  out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
  names(out) = names(col)
  return(out)
}



events.long = c('Simple (TRA/INV)', 'Deletion', 'Duplication', 'Chromoplexy', 'Rigma', 'Pyrgo', 'Chromothripsis', 'TIC', 'QRP', 'BFB', 'Double Minute', 'Complex DM', 'Tyfonas')
events.short = c('simple', 'del', 'dup', 'chromoplexy', 'rigma', 'pyrgo', 'chromothripsis', 'tic', 'qrp', 'bfb', 'dm', 'cpxdm','tyfonas')
event.colors  = c(`Simple (TRA/INV)`='#B2DF8A', `Deletion`='#A6CEE3', `Duplication`='#FB9A99',
                  `Chromoplexy`='#33A02C', `Rigma`='#1F78B4', `Pyrgo`='#E31A1C', `Chromothripsis`='darkgreen',
                  `TIC`='orchid4', `QRP`='gold2', `BFB`='firebrick', `Tyfonas`='gray25', `Double Minute`='darkorange',
                  `Complex DM`='darkorange3')



GENCODE = '/gpfs/commons/projects/TCGA/gdc-awg/TGCT/WGS/jabba/annotations/gencode.composite.collapsed.rds'
gencode = track.gencode(gencode=GENCODE, cached.dir=dirname(GENCODE), cached.path=GENCODE, cex.label=0.5, xaxis.cex.label=1, xaxis.unit=1e6, xaxis.suffix='MB')



## Get arguments
option_list = list(
  make_option(c("-s", "--samples"),         type='character', help="Sample IDs (comma delimited)"),
  make_option(c("-j", "--jabba_gg"),        type='character', help="Path to jabba.events.rds file"),
  make_option(c("-i", "--id_map"),          type='character', help="Map between old samples IDs in --samples and manuscript IDs"),
  make_option(c("-c", "--coverage"),        type='character', help="Coverage GRanges object (RDS)"),
  make_option(c("-f", "--field"),           type='character', help="Field in --coverage file to use when plotting coverage", default='foreground'),
  make_option(c("-b", "--bed"),             type='character', help="BED file with regions for plotting"),
  make_option(c("-p", "--padding"),         type='numeric',   help="Padding for taking subgraph associated with BED intervals"),
  make_option(c("-o", "--out_file"),        type='character', help="Figure output file (PNG)"))
opt = parse_args(OptionParser(option_list=option_list))



## Expand arguments
opt$samples = unlist(strsplit(opt$samples, ',',fixed=T))
opt$jabba_gg = unlist(strsplit(opt$jabba_gg, ',',fixed=T))
opt$coverage = unlist(strsplit(opt$coverage, ',',fixed=T))

if (!is.null(opt$jabba_ag)) {
  opt$jabba_ag = unlist(strsplit(opt$jabba_ag, ',',fixed=T))
}

## Update sample names
id.map = read.csv(opt$id_map, h=F, stringsAsFactors=F, sep='\t', col.names=c('old','new'))
opt$samples = id.map$new[match(opt$samples, id.map$old)]


## If this is PM12, force a specific ordering
if ('WCM12_PR_1' %in% opt$samples) {
  
  new.order = match(c('WCM12_PR_1', 'WCM12_BR_3', 'WCM12_BR_2'), opt$samples)
  opt$samples = opt$samples[new.order]
  opt$jabba_gg = opt$jabba_gg[new.order]
  opt$coverage = opt$coverage[new.order]

  if (!is.null(opt$jabba_ag)) {
    opt$jabba_ag = opt$jabba_ag[new.order]
  }

}


## Read BED, format intervals
## This is the first mandatory track
bed = read.csv(opt$bed, h=F, stringsAsFactors=F, sep='\t', col.names=c('chr','start','end','label'))
bed = makeGRangesFromDataFrame(bed, keep.extra.columns=T)

names(bed) = bed$label

bed.gtrack = gTrack(bed, 
                    col=alpha('black', 0.5), name=' ', height=3,  labels.suppress=F, gr.labelfield='label',
                    grl.labelfield='label',  vadj.label=-2, cex.label=3, xaxis.cex.tick=1,
                    xaxis.cex.label=1, xaxis.unit=1e6, xaxis.suffix='MB')


## Read tumor coverage
cov.gt = vector(mode='list', length=length(opt$samples))
for (i in 1:length(opt$samples)) {
  cov.gt[[i]] = gTrack(readRDS(opt$coverage[i]), y.field=opt$field, col=alpha('black', 0.2), name=' ')  
}



## Add genome graph colored by event type
jabba.gt = vector(mode='list', length=length(opt$samples))
jabba.fp = vector(mode='list', length=length(opt$samples))

shown.events = c()

for (i in 1:length(opt$samples)) {
  
  ## Read genome graph
  jabba.gg = readRDS(opt$jabba_gg[i])
  
  ## Reset node/edge coloring
  jabba.gg$nodes$mark(col='gray')
  jabba.gg$edges[type=='ALT']$mark(col='gray50')
  
  ## Mark edges by event type 
  for (e in 1:length(events.short)) {
    jabba.gg$edges[!is.na(get(events.short[e]))]$mark(col=event.colors[events.long[e]])
  }
  
  ## Build gtrack
  jabba.fp[[i]] = jabba.gg$copy$subgraph(bed, opt$padding)$footprint + 1E3
  shown.events = union(shown.events, event.colors[event.colors %in% na.omit(jabba.gg$copy$subgraph(bed, opt$padding)$edges$dt$col)])

  if (i == length(opt$samples)) {
    shown.events = event.colors[event.colors %in% shown.events]
    print(shown.events)
    jabba.gt[[i]] = jabba.gg$gtrack(name=opt$samples[i], colormaps=list(Events=shown.events))
  } else {
    jabba.gt[[i]] = jabba.gg$gtrack(name=opt$samples[i])
  }
  
  
}


## Plot union of footprints
jabba.fp = do.call(gr.reduce, jabba.fp)



## Combine gTracks
for (i in length(opt$samples):1) {

  if (i == length(opt$samples)) {
    plt.gtrack = c(cov.gt[[i]], jabba.gt[[i]])
  } else {
    plt.gtrack = c(plt.gtrack, cov.gt[[i]], jabba.gt[[i]])
  }
  
}



##########
## Plot ## 
##########

if (substr(opt$out_file, nchar(opt$out_file)-2, nchar(opt$out_file)) == 'png') {
  png(opt$out_file, width=2000, height=1000)  
} else {
  svg(opt$out_file, width=13, height=15)
}


plot(c(gencode, plt.gtrack), jabba.fp)
  
dev.off()

message(opt$out_file)
