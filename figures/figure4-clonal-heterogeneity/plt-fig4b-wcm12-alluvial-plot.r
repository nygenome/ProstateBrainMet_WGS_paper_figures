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
## Compare the clonality/timing of shared mutations in multiple VCFs 
libs = c('optparse', 'VariantAnnotation', 'gUtils', 'ggplot2', 'ggalluvial')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(width=200)



## Read MutationTimer VCF
readMtVcf = function(f, name) {
  
  x = VariantAnnotation::readVcf(f)
  
  res = rowRanges(x)
  mcols(res) = NULL
  mcols(res)[, paste0(name,'_timing')] = info(x)$CLS   
  
  return(res)
  
}



## Get arguments
option_list = list(
  make_option(c("-v", "--vcf"),          type='character', help="Comma-separated list of VCFs (plotted in the order in which they appear)"),
  make_option(c("-n", "--vcf_name"),     type='character', help="Comma-separated list of VCF sample names for plotting"),
  make_option(c("-o", "--out_file_pdf"), type='character', help="Output PDF"))
opt = parse_args(OptionParser(option_list=option_list))



opt$vcf = unlist(strsplit(opt$vcf, ','))
opt$vcf_name = unlist(strsplit(opt$vcf_name, ','))


col = c(uncalled='#000000', `subclonal`="#CC6677", `clonal [late]`="#DDCC77", `clonal [NA]`="#6699CC", `clonal [early]`="#117733")


## Read data, select variants with timing data
vcf = mapply(FUN=readMtVcf, f=opt$vcf, name=paste0('vcf',1:length(opt$vcf)))

for (i in 1:length(vcf)) {
  vcf[[i]] = vcf[[i]][!is.na(mcols(vcf[[i]])[,paste0('vcf',i,'_timing')])]
}

## Get union
vcf.union = Reduce(f=function(x,y) gr.merge(x, y, all=T), x=vcf)



## Format for alluvial plot
vcf.union.alluvial = data.frame(timing=unlist(mcols(vcf.union)[, paste0('vcf',1:length(opt$vcf),'_timing')]),
                                id=rep(1:length(vcf.union), length(opt$vcf)),
                                vcf=rep(paste0('vcf',1:length(opt$vcf)), each=length(vcf.union)))

vcf.union.alluvial = vcf.union.alluvial[!is.na(vcf.union.alluvial$timing), ]
vcf.union.alluvial$timing = factor(vcf.union.alluvial$timing, levels=rev(names(col)))

## Rename VCFs for plot
for (i in 1:length(opt$vcf)) {
  vcf.union.alluvial$vcf[vcf.union.alluvial$vcf == paste0('vcf',i)] = opt$vcf_name[i]    
}

vcf.union.alluvial$vcf = factor(vcf.union.alluvial$vcf, levels=opt$vcf_name)




##############
## Plotting ##
##############

pdf(opt$out_file_pdf, width=10, height=8)

ggplot(vcf.union.alluvial, aes(x=vcf, stratum=timing, alluvium=id, fill=timing, label=timing)) +
       scale_x_discrete(expand = c(.1, .1)) +
       scale_fill_manual(values=col) +
       geom_flow() +
       geom_stratum(alpha = 0.8) +
       geom_text(stat="stratum", size=3) +
       geom_text(stat='stratum', aes(label = after_stat(n)), nudge_y=-350, size=3) +
       ggtitle('Timing/clonality changes across all mutations') +
       xlab('Sample') + ylab('Mutations') +
       theme_bw() + 
       theme(legend.position='false') +
       theme(panel.grid.major=element_blank(), 
             panel.grid.minor=element_blank())  

dev.off()
message(opt$out_file_pdf)


## Write out source data
out.file.txt = gsub('\\.pdf$', '.txt', opt$out_file_pdf)

by.timing = tapply(vcf.union.alluvial$timing, vcf.union.alluvial$id, paste, collapse='->')
by.sample = tapply(vcf.union.alluvial$vcf, vcf.union.alluvial$id, paste, collapse='->')

if(!identical(names(by.timing), names(by.sample))) stop()

res = as.data.frame.matrix(table(by.timing, by.sample))
write.table(res, out.file.txt, row.names=T, col.names=T, quote=F, sep='\t')
message(out.file.txt)
