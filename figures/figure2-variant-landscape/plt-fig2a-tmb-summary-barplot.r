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
## Plot coding mutation burden for all samples

libs = c('optparse', 'reshape2', 'ggplot2')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200, scipen=999)



## Reformat a coding consequence string, replacing 
## underscores with spaces and capitalizing the first letter
reformat.conseqeunce = function(x) { 

  x = gsub('_', ' ', x, fixed=T)
  x = unlist(strsplit(x, ''))
  x[1] = toupper(x[1])
  x = paste(x, collapse='')

  return(x)

}



## Read arguments
option_list = list(
  make_option(c("-i", "--in_file"),      type='character', help="Output of init-tmb-summary-table.r"),
  make_option(c("-m", "--id_map"),       type='character', help="Map between IDs"),
  make_option(c("-O", "--sample_order"), type='character', help="Optional sample ordering"),
  make_option(c("-o", "--out_file"),     type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read input file, reformat 
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
dta = melt(dta, id.vars=c('sample','patient'), variable.name='consequence', value.name='count')
dta$sample = gsub('--.*', '', dta$sample)



## Read ID mapping 
id.map = read.csv(opt$id_map, h=F, stringsAsFactors=F, sep='\t', col.names=c('old','new'))
dta$sample = id.map$new[match(dta$sample, id.map$old)]


## Optionally update ordering
if (!is.null(opt$sample_order)) {

  opt$sample_order = scan(opt$sample_order, what=character(), quiet=T)
  dta$sample = factor(dta$sample, levels=opt$sample_order)

  dta = dta[order(dta$sample), ]
  dta$patient = factor(dta$patient, levels=unique(dta$patient))

}



## Subset to coding mutations
consequence.keep = c('transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant',
                     'stop_lost', 'start_lost', 'inframe_insertion', 'inframe_deletion', 'missense_variant')

dta = dta[dta$consequence %in% consequence.keep, ]
dta$Consequence = factor(dta$consequence, levels=rev(consequence.keep))
levels(dta$Consequence) = sapply(levels(dta$Consequence), reformat.conseqeunce)



##########
## Plot ##
##########

pdf(opt$out_file, width=15, height=10)

ymax = 1.1 * max(tapply(dta$count, dta$sample, sum))

ggplot(dta, aes(x=sample, y=count, fill=Consequence)) +
  geom_bar(position='stack', stat='identity') + 
  facet_grid(. ~ patient, space='free', scale='free') + 
  scale_fill_brewer(palette='Paired') +
  scale_y_continuous(limits=c(0, ymax), expand=c(0,0)) +
  xlab('Sample') +
  ylab('Coding SNV/INDEL Count') +
  theme_classic() + 
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

dev.off()
message(opt$out_file)



## Output source data
out.file.txt = gsub('\\.pdf$', '.txt', opt$out_file)

col.order = c('sample','patient','Consequence','count')
dta = dta[, col.order]
colnames(dta) = tolower(colnames(dta))
dta$patient = gsub('^PM','WCM',dta$patient)

write.table(dta, out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)
