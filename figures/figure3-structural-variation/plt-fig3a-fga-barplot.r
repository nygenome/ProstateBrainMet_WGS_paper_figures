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
## Generate barplot of FGA scores for each sample, as well as an in-line summary boxplot

libs = c('optparse', 'ggplot2', 'patchwork')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



allDup = function (x) duplicated(x) | duplicated(x, fromLast = TRUE)



## Adds an annotation track using geom_bar
annotationTrack = function(x, ano_field, x_field, facet_on=NULL, col_map=NULL) {
  
  ## Set height of all bars to 1, handle stacked barplots
  x$y = 1
  x = x[!duplicated(x[,x_field]), ]
  x[,ano_field] = as.factor(x[,ano_field])
  
  plt = ggplot(x, aes_string(fill=ano_field, y='y', x=x_field)) + 
    geom_bar(position="stack", stat="identity") +
    theme_void() + 
    theme(strip.background=element_blank(), 
          strip.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.x=element_blank(), 
          plot.margin=margin(0, 0, 0, 0, "cm"))
  
  if (!is.null(facet_on)) {
    facet.frm = as.formula(paste0('~',facet_on))
    plt = plt + facet_grid(facet.frm, scale='free', space='free') 
  }
  
  if (!is.null(col_map)) {
    plt = plt + scale_fill_manual(values=col_map)
  }
  
  plt
  
}



META_SITE_COLS = c(`Brain`='#FF7F00', `Prostate`='#90EE90', `Non-brain`='#1E90FF')



## Get arguments
option_list = list(
  make_option(c("-f", "--fga"),           type='character', help="FGA data"),
  make_option(c("-t", "--tn_file"),       type='character', help="Tab delimited tumor-normal pairing file"),
  make_option(c("-p", "--purity_ploidy"), type='character', help="Tab-delimited purity ploidy file"),
  make_option(c("-i", "--id_map"),        type='character', help="Map between old and new IDs"),
  make_option(c("-s", "--site_colors"),   type='character', help="Site colors"),
  make_option(c("-o", "--out_file"),      type='character', help="Figure output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



## Read data
fga = read.csv(opt$fga, h=F, stringsAsFactors=F, sep='\t', col.names=c('sample','fga'))
tn = read.csv(opt$tn_file, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','patient','site'))
pp = read.csv(opt$purity_ploidy, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','patient','site', 'purity','ploidy'))
id = read.csv(opt$id_map, h=F, stringsAsFactors=F, sep='\t', col.names=c('old', 'new'))
site.cols = read.csv(opt$site_colors, h=F, stringsAsFactors=F, sep='\t', col.names=c('site', 'color'))
site.cols = structure(site.cols$color, names=site.cols$site)


tn[, c('purity','ploidy')] = pp[match(tn$tumor, pp$tumor), c('purity','ploidy')]
tn$sample = id$new[match(tn$tumor, id$old)]
tn$pair_name = paste0(tn$tumor,'--',tn$normal)


## Align metadata
fga[, c('sample', 'patient','site','purity','ploidy')] = tn[match(fga$sample, tn$pair_name), c('sample', 'patient','site','purity','ploidy')]
fga$site = gsub('_', ' ', fga$site, fixed=T)

fga$site_simplified = sapply(fga$site, function(x) ifelse(x %in% c('Brain', 'Prostate'), x, 'Non-brain'))
fga$site_simplified = factor(fga$site_simplified, levels=c('Prostate', 'Non-brain', 'Brain'))

fga$Localization = ifelse(fga$site_simplified == 'Prostate', 'Prostate', 'Metastatic')
fga$Localization = factor(fga$Localization, levels=c('Prostate', 'Metastatic'))

fga$patient = gsub('PM', 'WCM', fga$patient, fixed=T)
fga$patient = factor(fga$patient, levels=names(sort(tapply(fga$fga, fga$patient, median), decreasing=T)))
fga$sample = factor(fga$sample, levels=fga$sample[order(fga$fga, decreasing=T)])



##########
## Plot ##
##########

pdf(opt$out_file, height=5, width=14)

## Boxplot
fga.bxp = fga
fga.bxp$site_simplified = factor(fga.bxp$site_simplified, levels=c('Non-brain','Brain','Prostate'))
bxp.site.simple = ggplot(fga.bxp, aes(x=site_simplified, y=fga, fill=site_simplified)) + 
                    geom_boxplot(outlier.shape=NA) + 
                    geom_jitter(width=0.2, height=0) +
                    scale_fill_manual(values=META_SITE_COLS) +
                    scale_y_continuous(limits=0:1, expand=c(0,0), position='right') +
                    xlab(' ') +
                    ylab(' ') +
                    labs(fill='Localization') +
                    theme_bw() +
                    theme(legend.position='right',
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          plot.margin=margin(0, 0, 0, 0, "cm"))

bxp.site.local = ggplot(fga, aes(x=Localization, y=fga, fill=Localization)) + 
                    geom_boxplot(outlier.shape=NA) + 
                    geom_jitter(width=0.2, height=0) +
                    scale_y_continuous(limits=0:1, expand=c(0,0), position='right') +
                    xlab(' ') +
                    ylab(' ') +
                    theme_bw() +
                    theme(legend.position='none', 
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank())


## Barplots
plt.site.main = ggplot(fga, aes(x=sample, y=fga, fill=site_simplified)) + 
  geom_bar(stat='identity') + 
  facet_grid(. ~ site_simplified, space='free', scale='free') +
  scale_y_continuous(limits=0:1, expand=c(0,0)) +
  xlab('Site') +
  ylab('FGA') +
  theme_bw() + 
  theme(axis.text.x=element_blank(), 
        axis.title.x=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.ticks.x=element_blank())

plt.site.ano = annotationTrack(x=fga, x_field='sample', ano_field='site', facet_on='site_simplified', col_map=site.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        legend.position='bottom')



plt.patient.main = ggplot(fga, aes(x=sample, y=fga, fill=site_simplified)) + 
                    geom_bar(stat='identity') + 
                    scale_fill_manual(values=META_SITE_COLS) +
                    facet_grid(.  ~ patient, space='free', scale='free') +
                    scale_y_continuous(limits=0:1, expand=c(0,0)) +
                    xlab('Site') +
                    ylab('FGA') +
                    theme_bw() + 
                    theme(axis.text.x=element_blank(), 
                          axis.title.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          strip.background = element_blank(),
                          panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),
                          strip.text.x = element_blank(),
                          plot.margin=margin(0, 0, 0, 0, "cm"),
                          panel.spacing=unit(2,'mm'),
                          legend.position='none')

plt.patient.ano = annotationTrack(x=fga, x_field='sample', ano_field='site', facet_on='patient', col_map=site.cols) +
                    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
                          legend.position='bottom') + 
                    labs(fill='Site')



## Draw plots
((plt.patient.main / plt.patient.ano) + plot_layout(heights=c(9,0.85)) | 
(bxp.site.simple / plot_spacer() + plot_layout(heights=c(9,0)))) + 
plot_layout(widths=c(8,1))


dev.off()
message(opt$out_file)




## Write list with sample ordering (for use in other figures)
out.file.txt = gsub('\\.pdf$', '.txt', opt$out_file)
out.txt = fga$sample[order(fga$patient, fga$sample)]
out.txt = paste(out.txt, collapse='\n')

sink(out.file.txt)
cat(out.txt,'\n')
sink()
message(out.file.txt)


## Write source data
out.file.txt = gsub('\\.pdf$', '.sourcedata.txt', opt$out_file)

col.order = c('sample','patient','site_simplified','fga')
fga = fga[, col.order]

write.table(fga, out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)
