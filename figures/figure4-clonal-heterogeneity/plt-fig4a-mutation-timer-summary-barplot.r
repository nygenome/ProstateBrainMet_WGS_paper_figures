#!/nfs/sw/R/R-3.6.1/bin/Rscript
################################################################################
### COPYRIGHT ##################################################################

# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2023) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Author: William F. Hooper & Timothy R. Chu 

################################################################# /COPYRIGHT ###
################################################################################
## Plot mutation timing summary
libs = c('optparse', 'ggplot2', 'patchwork','stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200)



col = c(`subclonal`="#CC6677", `clonal [late]`="#DDCC77", `clonal [NA]`="#6699CC", `clonal [early]`="#117733")



## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),  type='character', help="Tab delimited mutationtimer summary with columns sample,timing,count"),
  make_option(c("-m", "--id_map"),   type='character', help="ID Map with old/new sample names"),
  make_option(c("-p", "--pp"),       type='character', help="CSV with columns sample,purity,ploidy,patient,tumor,normal"),
  make_option(c("-o", "--out_file"), type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data, add purity/ploidy/date info
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F, sep='\t')
pp = read.csv(opt$pp, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','patient','site','purity','ploidy'))

dta[,c('purity','ploidy','patient','site')] = pp[match(dta$sample, pp$tumor),c('purity','ploidy','patient','site')]

dta$site[!dta$site %in% c('Prostate','Brain')] = 'Other'
dta$site = factor(dta$site, levels=c('Prostate','Brain','Other'))
dta$timing = factor(dta$timing, levels=c('clonal [early]', 'clonal [NA]', 'clonal [late]', 'subclonal'))

## Convert counts to proportions
dta$total = tapply(dta$count, dta$sample, sum)[dta$sample]
dta$proportion = dta$count / dta$total
dta$patient = gsub('[A-Z]+', '', gsub('(-|_).*','',dta$sample))


## Update IDs
id.map = read.csv(opt$id_map, h=F, stringsAsFactors=F, sep='\t', col.names=c('old','new'))
dta$sample = id.map$new[match(dta$sample, id.map$old)]

## Sort by proportion subclonal
subc = dta[dta$timing == 'subclonal', ]
dta$sample = factor(dta$sample, levels=subc$sample[order(subc$proportion, decreasing=T)])

## Deduplicate for purity/ploidy tracks
dta.dedup = dta[!duplicated(dta$sample),]

colnames(dta)[colnames(dta) == 'timing'] = 'Timing'



##########
## Plot ##
##########

pdf(opt$out_file, width=14, height=8)

## Barplots
purity = ggplot(dta.dedup, aes(x=sample, y=purity)) + 
  geom_bar(stat='identity') + 
  facet_grid(. ~ site, scales='free', space='free') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ylab('Purity') +
  xlab('Sample') +
  theme(axis.text.x=element_blank(),
        strip.background=element_blank(),
        strip.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.1, "lines"))
        # plot.margin=margin(0, 0, 0, 0, "cm"))


ploidy = ggplot(dta.dedup, aes(x=sample, y=ploidy)) + 
  geom_bar(stat='identity') + 
  facet_grid(. ~ site, scales='free', space='free') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ylab('Ploidy') +
  xlab('Sample') +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background=element_blank(),
        strip.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        plot.margin=margin(0, 0, 0, 0, "cm"))


mut.props = ggplot(dta, aes(x=sample, y=count, fill=Timing)) + 
  geom_bar(stat='identity', position='fill') + 
  facet_grid(. ~ site, scales='free', space='free') +
  scale_fill_manual(values=col) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ylab('Proportion') +
  xlab('Sample') +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.x = element_blank())


mut.counts = ggplot(dta, aes(x=sample, y=count, fill=Timing)) + 
  geom_bar(stat='identity', position='stack') + 
  facet_grid(. ~ site, scales='free', space='free') +
  scale_fill_manual(values=col) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  ylab('Count') +
  xlab('Sample') +
  ggtitle('Mutations assigned to each MutationTimer category') +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.x = element_blank(),
        legend.position='none')


mut.counts / mut.props / purity / ploidy + plot_layout(heights=c(6,3,1,1))

dev.off()
message(opt$out_file)


## Write out source data
out.file.txt = gsub('\\.pdf$', '.txt', opt$out_file)

colnames(dta) = tolower(colnames(dta))
dta$patient = paste0('WCM', dta$patient)

col.order = c('sample','patient','site','purity','ploidy','timing','count','proportion')
dta = dta[, col.order]

write.table(dta, out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)
