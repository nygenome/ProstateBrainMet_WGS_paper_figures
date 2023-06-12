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
## Generate junction burden heatmap

libs = c('optparse', 'ComplexHeatmap', 'circlize', 'dplyr', 'stringr')
invisible(suppressPackageStartupMessages(sapply(libs, require, character.only=T)))
options(scipen=999, width=200)



read.tp53.status = function(f) {

  x = read.csv(f, h=T, stringsAsFactors=F, sep='\t')
  x = x[x$gene == 'TP53', ]
  
  x$full_deletion = sapply(x$cn, function(y) any(as.numeric(unlist(strsplit(y,',',fixed=T))) == 0))

  x$status = x$snv_indel

  ## Naming in file doesn't match convention in manuscript oncoprint
  ## loss in file refers to any copy loss
  x$status[x$snv_indel == '' & x$loss] = 'Deletion'
  x$status[x$snv_indel == '' & x$full_deletion] = 'Loss'
  x$status[x$status == ''] = 'WT'
  return(x)

}



read.colormap = function(f) {
  
  x = read.csv(f,  h=F, stringsAsFactors=F, sep='\t', col.names=c('name','color'))
  x = structure(x$color, names=x$name)
  return(x)
  
}



TP53_COL = c(`WT`='grey87', 
             `Moderate impact SNV/INDEL`='#F7B851', 
             `High impact SNV/INDEL`='#000000', 
             `Deletion`='#3980FA', 
             `Loss`='#053389')

TMPRSS2_COL = c(`No fusion`='grey87', 
                `Fusion`='grey30')

WGD_COL = c(`no WGD`='grey87', 
          `WGD`='grey30')




## Get arguments
option_list = list(
  make_option(c("-i", "--in_file"),        type='character', help="Input junction summary by type"),
  make_option(c("-p", "--pp"),             type='character', help="Tab-delimited headerless file with columns tumor,normal,patient,site,purity,ploidy"),
  make_option(c("-m", "--id_map"),         type='character', help="Optional file to map sample IDs. Tab-delimited headerless with columns old,new"),
  make_option(c("-s", "--site_map"),       type='character', help="Tab-delimited headerless with columns site,color"),
  make_option(c("-a", "--patient_map"),    type='character', help="Tab-delimited headerless with columns patient,color"),
  make_option(c("-t", "--tp53_status"),    type='character', help="Tab-delimited headerless file for TP53 mutated/WT status"),
  make_option(c("-e", "--tmprss2_status"), type='character', help="Tab-delimited headerless file for Tmprss2:erg fusion status pos/neg"),
  make_option(c("-O", "--sample_order"),   type='character', help="Optional ordering of samples"),
  make_option(c("-o", "--out_file"),       type='character', help="Output directory"))
opt = parse_args(OptionParser(option_list=option_list))



## Read tumor-normal pairs
pp = read.csv(opt$pp, h=F, stringsAsFactors=F, sep='\t', col.names=c('tumor','normal','patient','site','purity','ploidy'))
pp$pair_name = paste0(pp$tumor,'--',pp$normal)
pp$site = gsub('_', ' ', pp$site, fixed=T)


## Read metadata and align
tp53 = read.tp53.status(opt$tp53_status)
tmprss2 = read.csv(opt$tmprss2_status, h=F, stringsAsFactors=F, sep='\t', col.names=c('pair_name', 'status'))

pp$tp53 = tp53$status[match(pp$pair_name, tp53$sample)]
pp$tmprss2 = tmprss2$status[match(pp$pair_name, tmprss2$pair_name)]


## Update names and enforce some orderings
pp$tmprss2 = factor(ifelse(pp$tmprss2 == 'Positive', 'Fusion', 'No fusion'), levels=c('No fusion', 'Fusion'))
pp$tp53 = gsub('MODERATE', 'Moderate impact SNV/INDEL', pp$tp53)
pp$tp53 = gsub('HIGH', 'High impact SNV/INDEL', pp$tp53)
pp$tp53 = factor(pp$tp53, levels=c('WT', 'Moderate impact SNV/INDEL', 'High impact SNV/INDEL', 'Deletion', 'Loss'))



## Read data, select tumor-normal pairs, convert to matrix
dta = read.csv(opt$in_file, h=T, stringsAsFactors=F)
junctions = select(dta, sample, njunc)
junctions[c('tumor', 'normal')] = str_split_fixed(junctions$sample, '--', 2)

rownames(dta) = dta$sample
dta = dta[pp$pair_name, setdiff(colnames(dta), 'sample')]
dta.orig = dta
dta = as.matrix(dta[, setdiff(colnames(dta), c('unclassified', 'njunc'))])



## Optionally update sample names
if (!is.null(opt$id_map)) {
  
  id.map = read.csv(opt$id_map, h=F, stringsAsFactors=F, sep='\t', col.names=c('old', 'new'))
  pp$name = id.map$new[match(pp$tumor, id.map$old)]
  pp$patient = gsub('PM','WCM',pp$patient)

} else {
  
  tn$name = tn$tumor
  
}
rownames(dta) = pp$name
rownames(dta.orig) = pp$name


## Transpose
dta = t(dta)

## Update event names
rownames(dta) = gsub('invdup','Inverted Duplication',rownames(dta))
rownames(dta) = gsub('tra','Translocation',rownames(dta))
rownames(dta) = gsub('inv','Inversion',rownames(dta))
rownames(dta) = gsub('del','Deletion',rownames(dta))
rownames(dta) = gsub('dup','Duplication',rownames(dta))
rownames(dta) = gsub('chromoplexy','Chromoplexy',rownames(dta))
rownames(dta) = gsub('rigma','Rigma',rownames(dta))
rownames(dta) = gsub('pyrgo','Pyrgo',rownames(dta))
rownames(dta) = gsub('chromothripsis','Chromothripsis',rownames(dta))
rownames(dta) = gsub('tic','TIC',rownames(dta))
rownames(dta) = gsub('qrp','QRP',rownames(dta))
rownames(dta) = gsub('cpxdm','Complex DM',rownames(dta))
rownames(dta) = gsub('bfb','BFB',rownames(dta))
rownames(dta) = gsub('dm','Double Minute',rownames(dta))
rownames(dta) = gsub('tyfonas','Tyfonas',rownames(dta))

##Combine complex DM and DM into just DM
dta["Double Minute", ] <- dta["Double Minute", ] + dta["Complex DM", ]
dta = dta[rownames(dta) != "Complex DM", ]

## Take ln(x+1)
dta = log1p(dta)



#################
## Annotations ##
#################

## Read color maps
site.colors = read.colormap(opt$site_map)
patient.colors = read.colormap(opt$patient_map)

if (!is.null(opt$id_map)) {

  names(patient.colors) = gsub('PM','WCM', names(patient.colors))

}


## Ordering the rows and columns
row_order = c('Duplication', 'Deletion', 'Inversion', 'Translocation',
              'BFB', 'Chromoplexy', 'Chromothripsis','Double Minute',
              'Inverted Duplication', 'Pyrgo', 'QRP', 'Rigma', 'TIC', 'Tyfonas')
row.subsections <- c(5,9)

#row_sums = rowSums(dta)
row_order = c('Deletion', 'Chromoplexy', 'Chromothripsis', 'Duplication',
              'TIC', 'Inversion', 'QRP','Double Minute',
              'Rigma', 'BFB', 'Tyfonas', 'Inverted Duplication', 'Pyrgo', 'Translocation')


## Optionally 
if (!is.null(opt$sample_order)) {

  sample.order = scan(opt$sample_order, what=character(), quiet=T)
  
  pp = pp[match(sample.order, pp$name), ]
  dta = dta[, match(sample.order, colnames(dta))]

  dta.orig = dta.orig[match(sample.order, rownames(dta.orig)), ]

  pp$patient = factor(pp$patient, levels=unique(pp$patient))

}


top.ano = ComplexHeatmap::HeatmapAnnotation(`Junction Count`=ComplexHeatmap::anno_barplot(junctions$njunc, baseline=0, gp=gpar(fill='black')),
                                            `TP53 Somatic Status`=pp$tp53,
                                            `TMPRSS2:ERG Fusion Status`=pp$tmprss2,
                                            Site=pp$site,
                                            gap=unit(c(1, 0, 0, 0), "mm"),
                                            gp = gpar(col="black"),
                                            col=list(`TP53 Somatic Status`=TP53_COL,
                                                     `TMPRSS2:ERG Fusion Status`=TMPRSS2_COL,
                                                     `WGD Status`=WGD_COL,
                                                     Site=site.colors))


burden.col = circlize::colorRamp2(breaks = c(0, quantile(dta, 0.95)), colors = c("grey95", "#d94801"))



##########
## Plot ##
##########

pdf(opt$out_file, height=10, width=18)

ComplexHeatmap::Heatmap(dta,
        col=burden.col,
        top_annotation=top.ano,
        border=T,
        row_names_side='left',
        column_title=NULL,
        heatmap_legend_param = list(title = "ln(burden)"),
        show_column_names=T,
        show_row_names=T,
        cluster_rows=F,
        cluster_columns=F,
        column_split=pp$patient,
        column_gap = unit(2.6, "mm"),
        show_heatmap_legend=T,
        row_order=row_order,
        column_order=1:ncol(dta),
        row_split=data.frame(rep(c("Simple", "Complex"), row.subsections)))

dev.off()
message(opt$out_file)


## Write out source data
out.file.txt = gsub('\\.pdf$', '.txt', opt$out_file)

dta.orig$dm = dta.orig$dm + dta.orig$cpxdm
dta.orig$total = dta.orig$njunc

dta.orig = dta.orig[, !colnames(dta.orig) %in% c('cpxdm', 'njunc', 'unclassified')]
colnames(dta.orig) = paste0(colnames(dta.orig),'_junctions')

res = cbind(dta.orig, pp)
col.order = c('name','patient','tp53','tmprss2','site', colnames(dta.orig))
res = res[, col.order]

write.table(res, out.file.txt, row.names=F, col.names=T, quote=F, sep='\t')
message(out.file.txt)
